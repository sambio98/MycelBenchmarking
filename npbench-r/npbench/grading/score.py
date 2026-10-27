"""Task-level aggregation.

The task score is the fraction of hard artifacts that pass (PROJECT.md §1, §7).
Soft artifacts — anything depending on a live service, per §9 — are scored and
reported but kept out of the headline number. Format compliance is a separate
axis, so a model that produces correct science under wrong column names is
visible as exactly that rather than as a failure (§3, §7).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .artifact import ArtifactScore, score_artifact
from .config import GradingConfig


@dataclass
class TaskScore:
    task_score: float                       # fraction of hard artifacts passing
    n_hard: int
    n_hard_passed: int
    format_compliance: float                # fraction of artifacts needing no column rename
    artifacts: list[ArtifactScore] = field(default_factory=list)
    soft_artifacts: list[ArtifactScore] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"task score      {self.task_score:.2f}  "
            f"({self.n_hard_passed}/{self.n_hard} artifacts passed)",
            f"format compliance {self.format_compliance:.2f}",
            "",
        ]
        for a in self.artifacts + self.soft_artifacts:
            mark = "PASS" if a.passed else "FAIL"
            tags = " [soft]" if a.soft else ""
            if not a.format_compliant:
                tags += " [renamed: " + ", ".join(
                    f"{k}->{v}" for k, v in sorted(a.renamed_columns.items())
                ) + "]"
            detail = f"  {a.reason}" if a.reason else ""
            lines.append(f"{mark} {a.score:.2f}  {a.path}{tags}{detail}")
        return "\n".join(lines)


def score_task(config: GradingConfig, truth_dir: Path, outputs_dir: Path) -> TaskScore:
    """Grade every declared artifact of one task."""
    scored = [
        score_artifact(spec, truth_dir / spec.truth, outputs_dir / spec.path)
        for spec in config.artifacts
    ]
    hard = [a for a in scored if not a.soft]
    soft = [a for a in scored if a.soft]
    passed = sum(1 for a in hard if a.passed)

    # Only artifacts that got far enough to have their columns read can be judged
    # on format. One that was missing outright is neither compliant nor not.
    judged = [a for a in scored if a.metrics]
    return TaskScore(
        task_score=passed / len(hard) if hard else 0.0,
        n_hard=len(hard),
        n_hard_passed=passed,
        format_compliance=(
            sum(1 for a in judged if a.format_compliant) / len(judged) if judged else 0.0
        ),
        artifacts=hard,
        soft_artifacts=soft,
    )
