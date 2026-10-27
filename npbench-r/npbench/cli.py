"""Command line entry points for the parts of NPBench-R that exist.

`validate`, `fetch`, `truth`, `selfgrade` and `grade` all work. `eval` belongs
to build-order step 4, the Inspect AI harness, and is deliberately absent rather
than stubbed.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from .build.manifest import fetch, load_manifest
from .build.validate import report, solution_steps, validate_task
from .grading import load_grading_config, score_task


def _validate(args: argparse.Namespace) -> int:
    findings = validate_task(args.task_dir)
    errors = [f for f in findings if f.level == "error"]
    print(report(findings) if findings else f"{args.task_dir.name}: bundle contract OK")
    return 1 if errors else 0


def _grade(args: argparse.Namespace) -> int:
    config = load_grading_config(args.task_dir / "grading.yaml")
    result = score_task(config, args.task_dir / "truth", args.outputs)
    print(result.summary())
    if args.expect_perfect and result.task_score < 1.0:
        print(
            "\nExpected 1.00. A reference solution that does not score 1.00 against "
            "its own truth means the grader or the truth is wrong; resolve that "
            "before building anything else.",
            file=sys.stderr,
        )
        return 1
    return 0


def _fetch(args: argparse.Namespace) -> int:
    entries = load_manifest(args.task_dir / "data.manifest.json")
    for line in fetch(args.task_dir, entries, force=args.force):
        print(line)
    print(f"{len(entries)} entries staged and verified")
    return 0


def _run_solution(task_dir: Path, out: Path, jobs: int) -> int:
    """Run the task's declared solution steps, writing artifacts to `out`."""
    for step in solution_steps(task_dir):
        command = step.format(task=task_dir, out=out, jobs=jobs)
        print(f"$ {command}", flush=True)
        result = subprocess.run(shlex.split(command))
        if result.returncode != 0:
            print(f"step failed: {command}", file=sys.stderr)
            return result.returncode
    return 0


def _truth(args: argparse.Namespace) -> int:
    return _run_solution(args.task_dir, args.task_dir / "truth", args.jobs)


def _selfgrade(args: argparse.Namespace) -> int:
    """Run the reference solution into a scratch directory and grade it.

    CLAUDE.md build order, step 3: this must come out at 1.00. It checks the
    grader and grading.yaml against the truth they were built from, not the
    science -- but if it fails, nothing downstream can be trusted.
    """
    with tempfile.TemporaryDirectory() as scratch:
        out = Path(scratch) / "outputs"
        failed = _run_solution(args.task_dir, out, args.jobs)
        if failed:
            return failed
        args.outputs = out
        args.expect_perfect = True
        return _grade(args)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="npbench", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="check a task bundle against the contract")
    validate.add_argument("task_dir", type=Path)
    validate.set_defaults(handler=_validate)

    grade = sub.add_parser("grade", help="grade an outputs directory against a task's truth")
    grade.add_argument("task_dir", type=Path)
    grade.add_argument("outputs", type=Path)
    grade.add_argument(
        "--expect-perfect",
        action="store_true",
        help="exit non-zero unless the task scores 1.00 (used by `make selfgrade`)",
    )
    grade.set_defaults(handler=_grade)

    fetch_cmd = sub.add_parser("fetch", help="stage and verify a task's data from its manifest")
    fetch_cmd.add_argument("task_dir", type=Path)
    fetch_cmd.add_argument("--force", action="store_true", help="refetch files already staged")
    fetch_cmd.set_defaults(handler=_fetch)

    truth = sub.add_parser("truth", help="run the reference solution and write truth/")
    truth.add_argument("task_dir", type=Path)
    truth.add_argument("--jobs", type=int, default=4)
    truth.set_defaults(handler=_truth)

    selfgrade = sub.add_parser(
        "selfgrade", help="run the reference solution into a scratch dir and expect 1.00"
    )
    selfgrade.add_argument("task_dir", type=Path)
    selfgrade.add_argument("--jobs", type=int, default=4)
    selfgrade.set_defaults(handler=_selfgrade)

    args = parser.parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
