"""Command line entry points for the parts of NPBench-R that exist.

`validate` and `grade` work today. The fetch, truth and eval commands belong to
build-order steps 3 and 4 and are deliberately absent rather than stubbed: a
command that pretends to fetch data is worse than one that is not there.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .build.validate import report, validate_task
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

    args = parser.parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
