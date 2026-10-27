"""The scorer: read the agent's outputs off the VM and grade them.

Split in two so the part that matters can be tested without a sandbox.
`collect` takes a reader callable; the Inspect scorer passes the sandbox's.

A missing file is not an error here. It is simply not collected, and
score_artifact turns that into 0.00 with a reason, which is what CLAUDE.md
step 4 requires: a failing run scores zero without crashing.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Awaitable, Callable

from inspect_ai.scorer import Score, Target, mean, scorer
from inspect_ai.solver import TaskState
from inspect_ai.util import sandbox

from ..grading import GradingConfig, TaskScore, load_grading_config, score_task

Reader = Callable[[str], Awaitable[str]]

# Not graded, but collected: it feeds the failure-mode judge (PROJECT.md §6).
METHODS = "METHODS.md"


async def collect(config: GradingConfig, destination: Path, read: Reader) -> list[str]:
    """Copy every declared artifact out of the agent's workspace. Returns what was found."""
    found = []
    for name in [a.path for a in config.artifacts] + [METHODS]:
        try:
            content = await read(f"outputs/{name}")
        except Exception:
            continue  # absent or unreadable: scored as a zero, not a crash
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content)
        found.append(name)
    return found


def summarise(result: TaskScore, collected: list[str]) -> dict:
    """Per-artifact detail for the eval log, including the format axis."""
    return {
        "task_score": result.task_score,
        "artifacts_passed": result.n_hard_passed,
        "artifacts_total": result.n_hard,
        "format_compliance": result.format_compliance,
        "collected": collected,
        "methods_written": METHODS in collected,
        "artifacts": [
            {
                "path": a.path,
                "score": a.score,
                "passed": a.passed,
                "soft": a.soft,
                "format_compliant": a.format_compliant,
                "renamed_columns": a.renamed_columns,
                "reason": a.reason,
                "metrics": {m.name: m.score for m in a.metrics},
            }
            for a in result.artifacts + result.soft_artifacts
        ],
    }


@scorer(metrics=[mean()])
def artifact_scorer(task_dir: Path):
    """Grade the fraction of hard artifacts that pass, per PROJECT.md §7."""
    config = load_grading_config(task_dir / "grading.yaml")

    async def score(state: TaskState, target: Target) -> Score:
        async def read(path: str) -> str:
            return await sandbox().read_file(path, text=True)

        with tempfile.TemporaryDirectory() as scratch:
            outputs = Path(scratch)
            collected = await collect(config, outputs, read)
            result = score_task(config, task_dir / "truth", outputs)
        return Score(
            value=result.task_score,
            answer=f"{result.n_hard_passed}/{result.n_hard} artifacts",
            explanation=result.summary(),
            metadata=summarise(result, collected),
        )

    return score
