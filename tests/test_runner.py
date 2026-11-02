"""Tests for the benchmark runner that do not require a live LLM.

A fake ``A1`` is patched into ``biomni.agent`` so the loop logic (record shape,
error/timeout handling, resume, data-lake-once, --no-log) can be exercised
offline.
"""

import json

import pytest

# Skip rather than error when biomni is absent: a collection-time ImportError
# here aborts the whole suite, including the NPBench-C grading tests, which have
# no biomni dependency.
ba = pytest.importorskip("biomni.agent")

from biomni_benchmark.runner import (
    infer_source,
    iter_qa_pairs,
    load_done_ids,
    required_api_key_env,
    run,
)


class FakeA1:
    instances = []

    def __init__(self, **kwargs):
        FakeA1.instances.append(kwargs)

    def go(self, prompt):
        if "boom" in prompt:
            raise RuntimeError("simulated API failure")
        if "slow" in prompt:
            return ["step1"], "TIMEOUT: Code execution timed out after 5 seconds"
        return ["gen", "exec", "final"], f"ANSWER to: {prompt}"


@pytest.fixture(autouse=True)
def patch_agent(monkeypatch):
    FakeA1.instances = []
    monkeypatch.setattr(ba, "A1", FakeA1)


@pytest.mark.parametrize(
    "model,source",
    [("gpt-4o", "OpenAI"), ("claude-sonnet-4-5", "Anthropic"), ("gemini-2.0-flash", "Gemini")],
)
def test_infer_source(model, source):
    assert infer_source(model) == source


def test_infer_source_unknown():
    assert infer_source("llama-3") is None
    assert required_api_key_env(None) is None


def test_iter_qa_pairs_and_done_ids(tmp_path):
    inp = tmp_path / "qa.jsonl"
    inp.write_text(
        json.dumps({"id": "x", "question": "Q1", "answer": "A1"})
        + "\n\n"  # blank line tolerated
        + json.dumps({"question": "Q2"})  # missing id -> positional
        + "\n"
    )
    pairs = list(iter_qa_pairs(inp, "question", "answer", "id"))
    assert pairs == [("x", "Q1", "A1"), ("1", "Q2", None)]

    out = tmp_path / "out.jsonl"
    out.write_text(json.dumps({"id": "x"}) + "\n" + json.dumps({"id": "y"}) + "\n")
    assert load_done_ids(out) == {"x", "y"}
    assert load_done_ids(tmp_path / "missing.jsonl") == set()


def _make_input(tmp_path):
    inp = tmp_path / "qa.jsonl"
    inp.write_text(
        json.dumps({"id": "a", "question": "normal q"})
        + "\n"
        + json.dumps({"id": "b", "question": "boom q"})
        + "\n"
        + json.dumps({"id": "c", "question": "slow q", "answer": "ref"})
        + "\n"
    )
    return inp


def _run(tmp_path, out, **overrides):
    kwargs = dict(
        input_path=str(_make_input(tmp_path)),
        output_path=str(out),
        model="gpt-4o",
        source="OpenAI",
        data_path="./biomni_data",
        temperature=None,
        timeout=600,
        question_field="question",
        answer_field="answer",
        id_field="id",
        limit=None,
        resume=False,
        save_log=True,
        use_tool_retriever=False,
        skip_data_lake=True,
    )
    kwargs.update(overrides)
    return run(**kwargs)


def test_run_records_and_classification(tmp_path):
    out = tmp_path / "out.jsonl"
    stats = _run(tmp_path, out)
    assert stats == {"ok": 1, "errored": 1, "timed_out": 1}

    recs = [json.loads(line) for line in out.read_text().splitlines()]
    by_id = {r["id"]: r for r in recs}
    assert by_id["a"]["biomni_answer"] == "ANSWER to: normal q"
    assert by_id["a"]["error"] is None and by_id["a"]["log"] == ["gen", "exec", "final"]
    assert by_id["b"]["error"].startswith("RuntimeError")
    assert by_id["c"]["timed_out"] is True
    assert by_id["c"]["reference_answer"] == "ref"
    # all required metadata keys present
    for r in recs:
        assert set(r) >= {"id", "question", "model", "source", "duration_seconds", "timestamp"}


def test_resume_skips_done(tmp_path):
    out = tmp_path / "out.jsonl"
    _run(tmp_path, out)
    before = len(out.read_text().splitlines())
    stats = _run(tmp_path, out, resume=True)
    assert stats == {"ok": 0, "errored": 0, "timed_out": 0}
    assert len(out.read_text().splitlines()) == before  # nothing appended


def test_data_lake_downloaded_once(tmp_path):
    out = tmp_path / "out.jsonl"
    _run(tmp_path, out, skip_data_lake=False)
    lakes = [kw["expected_data_lake_files"] for kw in FakeA1.instances]
    assert lakes[0] is None  # first agent triggers the one-time download
    assert all(x == [] for x in lakes[1:])  # later agents reuse files


def test_no_log_omits_trace(tmp_path):
    out = tmp_path / "out.jsonl"
    _run(tmp_path, out, save_log=False, limit=1)
    rec = json.loads(out.read_text().splitlines()[0])
    assert "log" not in rec
