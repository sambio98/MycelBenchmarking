"""Command line entry points for the parts of NPBench-R that exist.

`validate`, `fetch`, `truth`, `selfgrade`, `grade` and `calibrate` are all
implemented here. Running an eval is `make eval`, which shells out to
`inspect eval` against npbench/harness/task.py rather than wrapping it.

`calibrate` reads the scored results of runs that already happened; it does not
run models.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from .audit import dossier, render
from .build.manifest import fetch, load_manifest
from .build.truthdigest import write as write_truth_digest
from .calibrate import CalibrationError, calibrate, load_runs, report as calibration_report
from .build.validate import report, solution_steps, validate_task
from .grading import TruthError, load_grading_config, score_task


def _validate(args: argparse.Namespace) -> int:
    findings = validate_task(args.task_dir)
    errors = [f for f in findings if f.level == "error"]
    print(report(findings) if findings else f"{args.task_dir.name}: bundle contract OK")
    return 1 if errors else 0


def _grade(args: argparse.Namespace) -> int:
    config = load_grading_config(args.task_dir / "grading.yaml")
    truth = getattr(args, "truth", None) or args.task_dir / "truth"
    if not truth.is_dir():
        print(
            f"{truth} does not exist. Run `make truth TASK={args.task_dir.name}` first.",
            file=sys.stderr,
        )
        return 2
    try:
        result = score_task(config, truth, args.outputs)
    except TruthError as exc:
        # A truth-side fault is a bundle bug, not an agent failure, and must not
        # be reported as a score of zero.
        print(f"the task bundle is broken, not the submission:\n  {exc}", file=sys.stderr)
        return 3
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


def _calibrate(args: argparse.Namespace) -> int:
    try:
        runs = load_runs(args.runs_dir)
    except CalibrationError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(calibration_report(calibrate(runs)))
    return 0


def _digest(args: argparse.Namespace) -> int:
    """Rewrite a task's committed fingerprint of truth/.

    Run this after deliberately regenerating truth -- a new reference release,
    a changed method -- and commit the result, so the change is visible in a
    diff instead of being discovered later as an unexplained score shift.
    """
    for task_dir in args.task_dir:
        entries = write_truth_digest(task_dir)
        print(f"{task_dir.name}: {len(entries)} artifact(s) fingerprinted")
        for entry in entries:
            print(f"  {entry.path:34} {entry.rows:>8} rows  {entry.sha256[:12]}")
    return 0


def _audit(args: argparse.Namespace) -> int:
    """Per-task dossiers for the human curation pass.

    Exits 1 when any task has an unjudged finding, so it can gate a review the
    way `validate` gates the contract. A finding a reviewer has already judged
    is recorded in that task's task.yaml and does not fail the command.
    """
    reports = [dossier(task_dir) for task_dir in sorted(args.task_dir)]
    print("\n\n".join(render(report) for report in reports))
    open_findings = sum(len(report.findings) for report in reports)
    if open_findings:
        print(f"\n{open_findings} finding(s) across {len(reports)} task(s) need review",
              file=sys.stderr)
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
        "--truth", type=Path, default=None,
        help="grade against this directory instead of <task_dir>/truth",
    )
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

    audit = sub.add_parser(
        "audit",
        help="per-task dossier and the checks for the human curation pass",
    )
    audit.add_argument("task_dir", type=Path, nargs="+")
    audit.set_defaults(handler=_audit)

    digest = sub.add_parser(
        "digest",
        help="rewrite the committed fingerprint of a task's truth/ (truth/ is gitignored)",
    )
    digest.add_argument("task_dir", type=Path, nargs="+")
    digest.set_defaults(handler=_digest)

    calibration = sub.add_parser(
        "calibrate",
        help="per-artifact pass rates and drop recommendations over scored runs",
    )
    calibration.add_argument("runs_dir", type=Path)
    calibration.set_defaults(handler=_calibrate)

    args = parser.parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
