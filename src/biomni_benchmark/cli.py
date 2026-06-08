"""Command-line entry point for ``biomni-benchmark``.

Run an entire JSONL file of biomedical Q&A pairs through Biomni on a configurable
LLM and collect the raw outputs for later scoring::

    biomni-benchmark --input data/qa.jsonl --output results.jsonl --model gpt-4o
    biomni-benchmark --input data/qa.jsonl --output results.jsonl --model claude-sonnet-4-5
"""

from __future__ import annotations

import argparse
import os
import sys

from dotenv import load_dotenv

from . import __version__
from .runner import infer_source, required_api_key_env, run


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="biomni-benchmark",
        description="Benchmark Biomni over a JSONL file of Q&A pairs with a configurable LLM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", "-i", required=True, help="Input JSONL file of Q&A pairs.")
    p.add_argument("--output", "-o", required=True, help="Output JSONL file (appended to).")
    p.add_argument(
        "--model",
        "-m",
        required=True,
        help="LLM Biomni should run on, e.g. 'gpt-4o' or 'claude-sonnet-4-5'.",
    )
    p.add_argument(
        "--source",
        default=None,
        help="Biomni provider (OpenAI, Anthropic, Gemini, Bedrock, Ollama, Custom, ...). "
        "Default: inferred from the model name.",
    )
    p.add_argument(
        "--data-path",
        default="./biomni_data",
        help="Directory where Biomni stores its data lake.",
    )
    p.add_argument("--temperature", type=float, default=None, help="Sampling temperature.")
    p.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Per-question execution timeout in seconds.",
    )
    p.add_argument("--question-field", default="question", help="JSON key holding the question.")
    p.add_argument(
        "--answer-field",
        default="answer",
        help="JSON key holding the reference answer (kept for later scoring).",
    )
    p.add_argument("--id-field", default="id", help="JSON key holding a stable item id.")
    p.add_argument("--limit", type=int, default=None, help="Only run the first N questions.")
    p.add_argument(
        "--resume",
        action="store_true",
        help="Skip ids already present in the output file.",
    )
    p.add_argument(
        "--no-log",
        action="store_true",
        help="Do not store Biomni's full execution trace in the output.",
    )
    p.add_argument(
        "--no-tool-retriever",
        action="store_true",
        help="Disable Biomni's per-question tool-retrieval LLM call.",
    )
    p.add_argument(
        "--skip-data-lake",
        action="store_true",
        help="Skip the ~11GB data-lake download (text-only; some tools may be unavailable).",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv: list[str] | None = None) -> int:
    load_dotenv()
    args = build_parser().parse_args(argv)

    if not os.path.exists(args.input):
        print(f"error: input file not found: {args.input}", file=sys.stderr)
        return 2

    source = args.source or infer_source(args.model)
    if source is None:
        print(
            f"error: could not infer a provider from model {args.model!r}; "
            f"pass --source explicitly (e.g. --source OpenAI).",
            file=sys.stderr,
        )
        return 2

    key_env = required_api_key_env(source)
    if key_env and not os.environ.get(key_env):
        print(
            f"error: {key_env} is not set, but it is required for source {source!r}. "
            f"Export it or put it in a .env file.",
            file=sys.stderr,
        )
        return 2

    run(
        input_path=args.input,
        output_path=args.output,
        model=args.model,
        source=source,
        data_path=args.data_path,
        temperature=args.temperature,
        timeout=args.timeout,
        question_field=args.question_field,
        answer_field=args.answer_field,
        id_field=args.id_field,
        limit=args.limit,
        resume=args.resume,
        save_log=not args.no_log,
        use_tool_retriever=not args.no_tool_retriever,
        skip_data_lake=args.skip_data_lake,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
