# MycelBenchmarking

A single CLI command — **`biomni-benchmark`** — for benchmarking
[**Biomni**](https://github.com/snap-stanford/Biomni) (snap-stanford's biomedical
AI agent) over your own Q&A pairs, with the **underlying LLM fully configurable**
(OpenAI, Anthropic, and the other providers Biomni supports).

It runs every question in a JSONL file through Biomni and **collects the raw
outputs** — Biomni's answer, its full execution trace, timing, and metadata —
into a JSONL file so you can score them however you like afterward.

---

## Table of contents

- [Why this exists](#why-this-exists)
- [How it works](#how-it-works)
- [Install](#install)
- [Configure API keys](#configure-api-keys)
- [Input format](#input-format)
- [Quick start](#quick-start)
- [Choosing the model](#choosing-the-model)
- [All CLI flags](#all-cli-flags)
- [Output format](#output-format)
- [Resuming and crash safety](#resuming-and-crash-safety)
- [Cost and performance notes](#cost-and-performance-notes)
- [Scoring the results](#scoring-the-results-later)
- [Development & tests](#development--tests)
- [Troubleshooting](#troubleshooting)
- [Project layout](#project-layout)

---

## Why this exists

Biomni is a general-purpose biomedical agent: given a natural-language task it
plans, writes and executes code, calls bioinformatics tools, and returns an
answer. To evaluate it on a question set you normally have to write a bespoke
driver script, wire up the agent, handle failures, and decide how to swap the
LLM. This repo packages all of that into **one reproducible command** so you can
point Biomni at a file of Q&A pairs and a model of your choice and walk away.

The design goal is **clean data collection**: each question is answered
independently, every result is written the moment it finishes, and failures
never abort the run — so a benchmark over hundreds of questions produces a
complete, analyzable JSONL even if some items error or time out.

## How it works

For each question the tool:

1. Builds a **fresh Biomni `A1` agent** configured with your chosen model.
2. Calls `agent.go(question)` and captures `(log, answer)`.
3. Writes one JSON record (answer + trace + metadata) to the output file and
   flushes immediately.

A few Biomni-specific details shape the implementation (verified against the
installed Biomni source):

- **Fresh agent per question.** A single `A1` instance leaks conversation state
  across `go()` calls (it uses a hardcoded LangGraph thread id backed by an
  in-memory checkpointer) and exposes no `reset()`. Re-instantiating per question
  guarantees each answer is independent. This is cheap — `A1` builds no
  embeddings at construction time.
- **Data lake handled once.** On the first real run Biomni downloads a ~11 GB
  data lake into `--data-path`. The first agent triggers that one-time download;
  every subsequent agent is constructed with `expected_data_lake_files=[]` so it
  reuses the already-present files instead of re-checking S3. `--skip-data-lake`
  skips the download entirely for quick text-only runs.
- **Timeouts vs. errors.** Biomni returns timeouts as `"TIMEOUT: ..."` *text*
  (not an exception), so those are detected from the answer and flagged with
  `timed_out: true`. Hard failures (API/rate-limit errors,
  `GraphRecursionError`, etc.) are caught per question, recorded in the `error`
  field, and the loop continues.

## Install

Requires **Python ≥ 3.11**.

```bash
git clone <this-repo>
cd MycelBenchmarking
pip install -e .
```

This installs Biomni plus this package and registers the `biomni-benchmark`
command on your `PATH`.

## Configure API keys

Copy the template and fill in the key for the provider you want to test (the
file is loaded automatically via `python-dotenv`):

```bash
cp .env.example .env
# edit .env:  OPENAI_API_KEY=...   and/or   ANTHROPIC_API_KEY=...
```

Provider → required environment variable:

| Provider (`--source`) | Env var |
| --- | --- |
| `OpenAI` | `OPENAI_API_KEY` |
| `AzureOpenAI` | `OPENAI_API_KEY` (+ `OPENAI_ENDPOINT`) |
| `Anthropic` | `ANTHROPIC_API_KEY` |
| `Gemini` | `GEMINI_API_KEY` |
| `Groq` | `GROQ_API_KEY` |
| `Ollama` / `Custom` | none (local server) |

The CLI checks the relevant key up front and exits with a clear message if it's
missing, so you never start a long run that fails on the first API call.

## Input format

A **JSONL** file — one JSON object per line. By default the keys are
`question`, `answer` (the reference answer, kept for later scoring) and `id`:

```json
{"id": "q1", "question": "Which enzyme unwinds DNA during replication?", "answer": "Helicase."}
{"id": "q2", "question": "Name the four DNA nucleobases.", "answer": "Adenine, thymine, guanine, cytosine."}
```

- `id` is optional — if absent, items are numbered by position.
- `answer` (the reference) is optional and only carried through to the output for
  your later scoring; Biomni never sees it.
- If your keys differ, map them with `--question-field`, `--answer-field`,
  `--id-field`.

A ready-to-run sample lives at [`data/sample_qa.jsonl`](data/sample_qa.jsonl).

## Quick start

Smoke test — first question only, no data-lake download:

```bash
biomni-benchmark \
  --input data/sample_qa.jsonl \
  --output /tmp/out.jsonl \
  --model gpt-4o \
  --limit 1 \
  --skip-data-lake
```

Full benchmark run:

```bash
biomni-benchmark \
  --input data/qa.jsonl \
  --output results.jsonl \
  --model gpt-4o
```

## Choosing the model

The `--model` flag selects the LLM Biomni runs on. The provider is **inferred
from the model name** (`gpt-*`/`o3` → OpenAI, `claude-*` → Anthropic,
`gemini-*` → Gemini); override it with `--source` when inference can't tell
(e.g. a self-hosted model):

```bash
# OpenAI
biomni-benchmark -i data/qa.jsonl -o results_gpt4o.jsonl -m gpt-4o

# Anthropic
biomni-benchmark -i data/qa.jsonl -o results_claude.jsonl -m claude-sonnet-4-5

# Gemini
biomni-benchmark -i data/qa.jsonl -o results_gemini.jsonl -m gemini-2.0-flash

# Self-hosted / OpenAI-compatible endpoint
biomni-benchmark -i data/qa.jsonl -o results_local.jsonl \
  -m my-local-model --source Custom
```

To compare models, run the command once per model into separate output files
(the `model` and `source` are recorded in every record too).

## All CLI flags

| Flag | Default | Purpose |
| --- | --- | --- |
| `--input`, `-i` | *(required)* | Input JSONL file of Q&A pairs. |
| `--output`, `-o` | *(required)* | Output JSONL file (appended to). |
| `--model`, `-m` | *(required)* | LLM Biomni runs on, e.g. `gpt-4o`, `claude-sonnet-4-5`. |
| `--source` | inferred | Provider: `OpenAI`, `AzureOpenAI`, `Anthropic`, `Gemini`, `Bedrock`, `Ollama`, `Groq`, `Custom`. |
| `--data-path` | `./biomni_data` | Where Biomni stores its data lake. |
| `--temperature` | Biomni default | Sampling temperature. |
| `--timeout` | `600` | Per-question execution timeout (seconds). |
| `--question-field` | `question` | JSONL key holding the question. |
| `--answer-field` | `answer` | JSONL key holding the reference answer. |
| `--id-field` | `id` | JSONL key holding a stable item id. |
| `--limit N` | all | Only run the first N questions (great for testing). |
| `--resume` | off | Skip ids already present in the output file. |
| `--no-log` | off | Don't store Biomni's full execution trace. |
| `--no-tool-retriever` | off | Disable Biomni's per-question tool-retrieval LLM call. |
| `--skip-data-lake` | off | Skip the ~11 GB download (text-only; some tools unavailable). |

Run `biomni-benchmark --help` for the same reference at the command line.

## Output format

Results stream to `--output`, one JSON object per line, appended as each
question finishes:

```json
{
  "id": "q1",
  "question": "Which enzyme unwinds DNA during replication?",
  "reference_answer": "Helicase.",
  "biomni_answer": "The enzyme is helicase, which ...",
  "model": "gpt-4o",
  "source": "OpenAI",
  "duration_seconds": 42.1,
  "timestamp": "2026-06-08T12:00:00+00:00",
  "error": null,
  "timed_out": false,
  "log": ["...step-by-step execution trace..."]
}
```

| Field | Meaning |
| --- | --- |
| `id` | Item id (from your file or positional). |
| `question` | The prompt sent to Biomni. |
| `reference_answer` | Your gold answer, carried through (null if absent). |
| `biomni_answer` | Biomni's final answer string (null if the call raised). |
| `model` / `source` | The LLM and provider used. |
| `duration_seconds` | Wall-clock time for this question. |
| `timestamp` | UTC ISO-8601 time the record was written. |
| `error` | `null`, or `"<ExceptionType>: <message>"` on a hard failure. |
| `timed_out` | `true` if Biomni reported a timeout. |
| `log` | Full execution trace (omitted with `--no-log`). |

Because it's JSONL, you can load it directly:

```python
import json
records = [json.loads(l) for l in open("results.jsonl")]
```

## Resuming and crash safety

Each result is written and flushed immediately, so an interrupted run keeps
everything completed so far. Re-run the **same command with `--resume`** to skip
ids already in the output file and continue where you left off:

```bash
biomni-benchmark -i data/qa.jsonl -o results.jsonl -m gpt-4o --resume
```

## Cost and performance notes

- **Data lake**: ~11 GB downloaded once into `--data-path` on the first real
  run. Use `--skip-data-lake` for quick checks where you don't need
  file-backed tools.
- **Tool retriever**: with the default (on), Biomni makes one extra LLM call per
  question to select relevant tools. `--no-tool-retriever` removes that overhead
  and cost if you don't need dynamic tool selection.
- **Timeouts**: `--timeout` bounds each question; long agentic runs may need it
  raised above the 600 s default.

## Scoring the results (later)

This tool deliberately **does not score** — it collects raw outputs so you can
apply whatever rubric fits your data (exact match, multiple-choice extraction,
or an LLM-as-judge). A judge step is easy to add on top: read `results.jsonl`,
compare `biomni_answer` against `reference_answer`, and write a score per record.
Open an issue/PR if you'd like a built-in `--judge` mode.

## Development & tests

```bash
pip install -e ".[dev]"
pytest -q
```

The test suite (`tests/test_runner.py`) stubs the Biomni agent, so it runs fully
**offline** (no API key, no data lake) and covers: source inference, JSONL
parsing, record shape, error/timeout classification, `--resume`, the
data-lake-once behavior, and `--no-log`.

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| `could not infer a provider from model '...'` | Pass `--source` explicitly. |
| `OPENAI_API_KEY is not set ...` | Export the key or add it to `.env`. |
| First run is slow / downloading gigabytes | That's the one-time data lake; use `--skip-data-lake` to avoid it. |
| Many records have `timed_out: true` | Raise `--timeout`. |
| A few records have `error` set | Inspect the message; the run still completed the rest. Re-run with `--resume` to retry just those after removing them, or investigate the specific questions. |

## Project layout

```
MycelBenchmarking/
├── pyproject.toml              # package metadata, deps, console_scripts entry
├── README.md
├── .env.example                # OPENAI_API_KEY / ANTHROPIC_API_KEY templates
├── .gitignore
├── data/
│   └── sample_qa.jsonl         # 3-item smoke-test dataset
├── src/
│   └── biomni_benchmark/
│       ├── __init__.py
│       ├── cli.py              # argparse front end, .env loading, validation
│       └── runner.py           # load JSONL → run Biomni per item → write JSONL
└── tests/
    └── test_runner.py          # offline tests (stubbed agent)
```
