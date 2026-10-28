"""Difficulty calibration over scored runs.

PROJECT.md §8: difficulty is measured, not designed. Build more candidate
artifacts than needed, run two or three strong models on each at least twice,
and keep the artifacts whose pass rate falls between roughly 0.1 and 0.8.
Artifacts every model passes and artifacts no model passes both carry zero
information about model differences.

This is that arithmetic. It does not run models; it reads the scored results of
runs that already happened, one JSON file per run, in the shape
`npbench.harness.scorer.summarise` produces plus the identifiers of the run:

    {"task": "02_salinispora_bgc", "model": "anthropic/claude-opus-5",
     "epoch": 1, "task_score": 0.43, "artifacts": [{"path": ..., "passed": ...}]}

It refuses to call a pass rate meaningful when there are too few runs behind
it. Two runs give a rate that can only be 0.0, 0.5 or 1.0, and reading that as
"outside the band" would drop a good artifact on no evidence.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean

BAND_LOW = 0.10
BAND_HIGH = 0.80

# Below this many runs a pass rate is reported but never acted on. PROJECT.md
# §8 asks for two or three models run at least twice, so four is the floor at
# which a rate can even distinguish "most" from "some".
MIN_RUNS_FOR_A_VERDICT = 4

INFORMATIVE = "informative"
TOO_EASY = "too_easy"
TOO_HARD = "too_hard"
UNDECIDED = "undecided"


class CalibrationError(ValueError):
    """A run record is not usable."""


@dataclass(frozen=True)
class Run:
    task: str
    model: str
    epoch: int
    task_score: float
    passed: dict[str, bool]        # artifact path -> passed
    soft: frozenset[str]           # artifact paths excluded from the headline


@dataclass
class ArtifactVerdict:
    task: str
    path: str
    n_runs: int
    n_passed: int
    verdict: str

    @property
    def pass_rate(self) -> float:
        return self.n_passed / self.n_runs if self.n_runs else 0.0


@dataclass
class Calibration:
    artifacts: list[ArtifactVerdict] = field(default_factory=list)
    model_scores: dict[tuple[str, str], float] = field(default_factory=dict)
    n_runs: int = 0

    def drops(self) -> list[ArtifactVerdict]:
        """Artifacts to drop, worst first: the ones carrying no information."""
        candidates = [a for a in self.artifacts if a.verdict in {TOO_EASY, TOO_HARD}]
        return sorted(candidates, key=lambda a: (abs(a.pass_rate - 0.5), a.task, a.path),
                      reverse=True)

    def task_scores_after_drops(self) -> dict[tuple[str, str], float]:
        """Recompute each model's task score over the artifacts that survive.

        The headline number PROJECT.md §8 targets at 0.3 to 0.5 is the score on
        the kept artifacts, not on the candidate set, so dropping changes it and
        the new value is what should be reported.
        """
        dropped = {(a.task, a.path) for a in self.drops()}
        kept: dict[tuple[str, str], list[float]] = {}
        for run in self._runs:
            scores = [
                1.0 if passed else 0.0
                for path, passed in run.passed.items()
                if path not in run.soft and (run.task, path) not in dropped
            ]
            if scores:
                kept.setdefault((run.model, run.task), []).append(mean(scores))
        return {key: mean(values) for key, values in sorted(kept.items())}

    _runs: list[Run] = field(default_factory=list, repr=False)


def load_run(path: Path) -> Run:
    try:
        raw = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise CalibrationError(f"{path}: not valid JSON: {exc}") from exc
    for field_name in ("task", "model", "artifacts"):
        if field_name not in raw:
            raise CalibrationError(f"{path}: missing required field {field_name!r}")
    passed: dict[str, bool] = {}
    soft: set[str] = set()
    for entry in raw["artifacts"]:
        if "path" not in entry or "passed" not in entry:
            raise CalibrationError(f"{path}: an artifact entry lacks 'path' or 'passed'")
        if entry["path"] in passed:
            raise CalibrationError(f"{path}: artifact {entry['path']!r} appears twice")
        passed[entry["path"]] = bool(entry["passed"])
        if entry.get("soft"):
            soft.add(entry["path"])
    return Run(
        task=str(raw["task"]),
        model=str(raw["model"]),
        epoch=int(raw.get("epoch", 1)),
        task_score=float(raw.get("task_score", 0.0)),
        passed=passed,
        soft=frozenset(soft),
    )


def load_runs(directory: Path) -> list[Run]:
    files = sorted(directory.glob("*.json"))
    if not files:
        raise CalibrationError(f"{directory}: holds no run JSON files")
    return [load_run(f) for f in files]


def classify(n_runs: int, n_passed: int) -> str:
    if n_runs < MIN_RUNS_FOR_A_VERDICT:
        return UNDECIDED
    rate = n_passed / n_runs
    if rate > BAND_HIGH:
        return TOO_EASY
    if rate < BAND_LOW:
        return TOO_HARD
    return INFORMATIVE


def calibrate(runs: list[Run]) -> Calibration:
    """Per-artifact pass rates and per-model task scores over a set of runs."""
    tally: dict[tuple[str, str], list[bool]] = {}
    for run in runs:
        for path, passed in run.passed.items():
            if path in run.soft:
                continue  # soft artifacts are reported, never in the headline
            tally.setdefault((run.task, path), []).append(passed)

    verdicts = [
        ArtifactVerdict(task=task, path=path, n_runs=len(results),
                        n_passed=sum(results), verdict=classify(len(results), sum(results)))
        for (task, path), results in sorted(tally.items())
    ]
    scores: dict[tuple[str, str], list[float]] = {}
    for run in runs:
        scores.setdefault((run.model, run.task), []).append(run.task_score)

    result = Calibration(
        artifacts=verdicts,
        model_scores={k: mean(v) for k, v in sorted(scores.items())},
        n_runs=len(runs),
    )
    result._runs = list(runs)
    return result


def report(result: Calibration) -> str:
    lines = [f"{result.n_runs} run(s)", ""]
    counts = {v: sum(1 for a in result.artifacts if a.verdict == v)
              for v in (INFORMATIVE, TOO_EASY, TOO_HARD, UNDECIDED)}
    lines.append(f"artifacts: {len(result.artifacts)}  "
                 + "  ".join(f"{k}={v}" for k, v in counts.items() if v))
    lines.append("")
    lines.append(f"{'task':28}{'artifact':34}{'runs':>5}{'rate':>7}  verdict")
    for a in result.artifacts:
        lines.append(f"{a.task:28}{a.path:34}{a.n_runs:5}{a.pass_rate:7.2f}  {a.verdict}")

    if counts[UNDECIDED]:
        lines += ["", f"{counts[UNDECIDED]} artifact(s) have fewer than "
                      f"{MIN_RUNS_FOR_A_VERDICT} runs. A pass rate from two runs can only "
                      f"be 0.00, 0.50 or 1.00, so these are reported and not acted on."]

    drops = result.drops()
    lines += ["", "recommended drops (PROJECT.md section 8), least informative first:"]
    lines += [f"  {a.task} {a.path} (rate {a.pass_rate:.2f}, {a.verdict})" for a in drops] \
        or ["  none; every decided artifact is in the 0.10 to 0.80 band"]

    lines += ["", f"{'model':38}{'task':28}{'as built':>9}{'after drops':>13}"]
    after = result.task_scores_after_drops()
    for (model, task), score in result.model_scores.items():
        revised = after.get((model, task))
        lines.append(f"{model:38}{task:28}{score:9.2f}"
                     + (f"{revised:13.2f}" if revised is not None else f"{'n/a':>13}"))
    return "\n".join(lines)
