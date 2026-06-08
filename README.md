# MycelBenchmarking

A single CLI command — `biomni-benchmark` — for benchmarking
[**Biomni**](https://github.com/snap-stanford/Biomni) (snap-stanford's
biomedical AI agent) over your own Q&A pairs, with the **underlying LLM fully
configurable** (OpenAI, Anthropic, and the other providers Biomni supports).

It runs every question through Biomni and **collects the raw outputs** (Biomni's
answer, its execution trace, timing, and metadata) into a JSONL file so you can
score them however you like afterward.

## Install

Requires Python ≥ 3.11.

```bash
pip install -e .
```

This installs Biomni and registers the `biomni-benchmark` command.

## Configure API keys

Copy `.env.example` to `.env` and fill in the key for the provider you want to
test (the file is loaded automatically):

```bash
cp .env.example .env
# then edit .env: OPENAI_API_KEY=... and/or ANTHROPIC_API_KEY=...
```

## Input format

A JSONL file with one object per line. By default the keys are `question`,
`answer` (the reference answer, kept for later scoring), and `id`:

```json
{"id": "q1", "question": "Which enzyme unwinds DNA during replication?", "answer": "Helicase."}
```

If your keys differ, point at them with `--question-field`, `--answer-field`,
and `--id-field`. A sample file lives at `data/sample_qa.jsonl`.

## Usage

Run the whole benchmark with one command. The `--model` flag chooses the LLM
Biomni runs on — the provider is inferred from the model name (or set it
explicitly with `--source`):

```bash
# Benchmark Biomni on an OpenAI model
biomni-benchmark --input data/qa.jsonl --output results.jsonl --model gpt-4o

# ...or on an Anthropic model
biomni-benchmark --input data/qa.jsonl --output results.jsonl --model claude-sonnet-4-5
```

Cheap smoke test (no ~11 GB data lake, first item only):

```bash
biomni-benchmark -i data/sample_qa.jsonl -o /tmp/out.jsonl -m gpt-4o --limit 1 --skip-data-lake
```

### Useful flags

| Flag | Purpose |
| --- | --- |
| `--model`, `-m` | LLM Biomni runs on, e.g. `gpt-4o`, `claude-sonnet-4-5`. **Required.** |
| `--source` | Provider (`OpenAI`, `Anthropic`, `Gemini`, `Bedrock`, `Ollama`, `Custom`, …). Default: inferred from the model. |
| `--data-path` | Where Biomni stores its data lake (default `./biomni_data`). |
| `--temperature` | Sampling temperature. |
| `--timeout` | Per-question timeout in seconds (default 600). |
| `--limit N` | Only run the first N questions. |
| `--resume` | Skip ids already present in the output file. |
| `--no-log` | Don't store Biomni's full execution trace. |
| `--no-tool-retriever` | Disable Biomni's per-question tool-retrieval LLM call. |
| `--skip-data-lake` | Skip the ~11 GB download (text-only; some tools unavailable). |
| `--question-field` / `--answer-field` / `--id-field` | Override the JSONL keys. |

## Output format

Results stream to the `--output` file (one JSON object per line, appended as each
question finishes — so a crash mid-run keeps everything completed so far, and
`--resume` picks up where you left off):

```json
{
  "id": "q1",
  "question": "...",
  "reference_answer": "...",
  "biomni_answer": "...",
  "model": "gpt-4o",
  "source": "OpenAI",
  "duration_seconds": 42.1,
  "timestamp": "2026-06-08T12:00:00+00:00",
  "error": null,
  "timed_out": false,
  "log": ["...step traces..."]
}
```

## Notes

- Each question runs in a **fresh Biomni agent** so there's no context bleed
  between questions.
- On the first real run Biomni downloads a **~11 GB data lake** into
  `--data-path` (one time). Use `--skip-data-lake` for quick text-only checks.
- Timeouts are recorded with `timed_out: true` (Biomni returns them as text, not
  an exception); hard failures are captured in `error` without aborting the run.
