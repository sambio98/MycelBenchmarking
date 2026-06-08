"""Core benchmark loop: load Q&A pairs, run each through Biomni, write results.

Design notes (verified against snap-stanford/Biomni ``main``):

* A single ``A1`` agent leaks conversation state across ``go()`` calls (a
  hardcoded ``thread_id`` plus a lifetime ``MemorySaver``) and exposes no
  ``reset()``. To keep every question independent we build a **fresh agent per
  question**. This is cheap because ``A1`` builds no embeddings at init.
* ``A1.go(prompt)`` returns ``(log: list[str], answer: str)``.
* Timeouts do **not** raise -- they come back as ``"TIMEOUT: ..."`` text inside
  the answer. Hard errors (API/rate-limit, ``GraphRecursionError``) *can* raise,
  so each call is wrapped in ``try/except``.
* The ~11 GB data lake downloads once; ``expected_data_lake_files=[]`` skips it
  entirely (used both for ``--skip-data-lake`` and for every agent after the
  first, which already triggered the one-time download).
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone

from tqdm import tqdm

# Prefix -> Biomni ``source`` string, mirroring Biomni's own inference.
_PREFIX_SOURCE = {
    "claude-": "Anthropic",
    "anthropic.claude-": "Bedrock",
    "azure-": "AzureOpenAI",
    "gpt-oss": "Ollama",
    "gpt-": "OpenAI",
    "o1": "OpenAI",
    "o3": "OpenAI",
    "o4": "OpenAI",
    "gemini-": "Gemini",
}

# Provider -> environment variable that must hold its API key.
_SOURCE_API_KEY_ENV = {
    "OpenAI": "OPENAI_API_KEY",
    "AzureOpenAI": "OPENAI_API_KEY",
    "Anthropic": "ANTHROPIC_API_KEY",
    "Gemini": "GEMINI_API_KEY",
    "Groq": "GROQ_API_KEY",
}


def infer_source(model: str) -> str | None:
    """Infer the Biomni provider ``source`` from a model id, or ``None``."""
    lowered = model.lower()
    for prefix, source in _PREFIX_SOURCE.items():
        if lowered.startswith(prefix):
            return source
    return None


def required_api_key_env(source: str | None) -> str | None:
    """Return the env var that must be set for ``source`` (or ``None``)."""
    if source is None:
        return None
    return _SOURCE_API_KEY_ENV.get(source)


def iter_qa_pairs(path, question_field, answer_field, id_field):
    """Yield ``(id, question, reference_answer)`` from a JSONL file.

    Blank lines are skipped. Missing ids are auto-assigned by position. A missing
    question field is a fatal error (we cannot benchmark without a question).
    """
    with open(path, encoding="utf-8") as fh:
        index = 0
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_no}: invalid JSON: {exc}") from exc
            if question_field not in obj:
                raise ValueError(
                    f"{path}:{line_no}: missing question field {question_field!r}. "
                    f"Use --question-field to point at the right key."
                )
            item_id = obj.get(id_field, index)
            yield str(item_id), obj[question_field], obj.get(answer_field)
            index += 1


def load_done_ids(output_path) -> set[str]:
    """Return the set of ids already present in an existing output file."""
    done: set[str] = set()
    try:
        with open(output_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    done.add(str(json.loads(line)["id"]))
                except (json.JSONDecodeError, KeyError):
                    continue
    except FileNotFoundError:
        pass
    return done


def _build_agent(A1, *, data_path, model, source, timeout, use_tool_retriever, skip_data_lake):
    """Instantiate a fresh ``A1`` agent for a single question."""
    kwargs = dict(
        path=data_path,
        llm=model,
        timeout_seconds=timeout,
        use_tool_retriever=use_tool_retriever,
    )
    if source is not None:
        kwargs["source"] = source
    # ``[]`` skips the data-lake download/check; ``None`` triggers the one-time
    # download. After the first agent has downloaded, callers pass skip=True.
    kwargs["expected_data_lake_files"] = [] if skip_data_lake else None
    return A1(**kwargs)


def run(
    *,
    input_path,
    output_path,
    model,
    source,
    data_path,
    temperature,
    timeout,
    question_field,
    answer_field,
    id_field,
    limit,
    resume,
    save_log,
    use_tool_retriever,
    skip_data_lake,
):
    """Run the benchmark and stream results to ``output_path`` as JSONL."""
    # Import here so ``--help`` and arg validation work without Biomni installed.
    from biomni.agent import A1

    if temperature is not None:
        # Biomni reads temperature from its runtime config rather than A1 kwargs.
        from biomni.config import default_config

        default_config.temperature = temperature

    done_ids = load_done_ids(output_path) if resume else set()
    if done_ids:
        print(f"[resume] {len(done_ids)} item(s) already in {output_path}; skipping them.")

    pairs = list(iter_qa_pairs(input_path, question_field, answer_field, id_field))
    if limit is not None:
        pairs = pairs[:limit]

    pending = [p for p in pairs if p[0] not in done_ids]
    print(
        f"Loaded {len(pairs)} question(s) from {input_path}; "
        f"{len(pending)} to run with model={model!r} source={source!r}."
    )

    # The first agent does the one-time data-lake check/download (unless skipped);
    # every later agent reuses the already-present files.
    download_pending = not skip_data_lake

    n_done = n_error = n_timeout = 0
    with open(output_path, "a", encoding="utf-8") as out:
        for item_id, question, reference in tqdm(pending, desc="Biomni", unit="q"):
            agent_skip_lake = skip_data_lake or not download_pending
            started = time.time()
            error = None
            answer = None
            log = None
            try:
                agent = _build_agent(
                    A1,
                    data_path=data_path,
                    model=model,
                    source=source,
                    timeout=timeout,
                    use_tool_retriever=use_tool_retriever,
                    skip_data_lake=agent_skip_lake,
                )
                download_pending = False  # first successful build settled the lake
                log, answer = agent.go(question)
            except Exception as exc:  # noqa: BLE001 - keep the loop alive
                error = f"{type(exc).__name__}: {exc}"

            duration = round(time.time() - started, 3)
            timed_out = bool(answer) and "TIMEOUT" in answer[:200].upper()

            record = {
                "id": item_id,
                "question": question,
                "reference_answer": reference,
                "biomni_answer": answer,
                "model": model,
                "source": source,
                "duration_seconds": duration,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "error": error,
                "timed_out": timed_out,
            }
            if save_log:
                record["log"] = log

            out.write(json.dumps(record, ensure_ascii=False) + "\n")
            out.flush()

            if error is not None:
                n_error += 1
            elif timed_out:
                n_timeout += 1
            else:
                n_done += 1

    print(
        f"\nDone. ok={n_done} errored={n_error} timed_out={n_timeout} "
        f"-> {output_path}",
        file=sys.stderr,
    )
    return {"ok": n_done, "errored": n_error, "timed_out": n_timeout}
