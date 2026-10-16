"""The grader: a pure function of (submission, gold, spec).

Purity contract (NPBench-C spec 3.3). Same inputs produce a byte-identical
result on any machine at any time:

  * no network access
  * no LLM call anywhere in the scoring path
  * no wall-clock, no RNG, no hash-ordering or filesystem-ordering dependence
  * every float reduction in a fixed order, every score rounded to 6 decimals
  * inputs are never mutated

``grade()`` takes three in-memory mappings so that the pure core is trivially
testable 100x over. Disk loading lives in npbench_c.io, outside this module.
"""

from __future__ import annotations

import json
import re
from typing import Any, Mapping

from . import primitives as P
from .ladder import (
    LadderResult,
    chance_floor,
    discriminating_chance_floor,
    score_ladder,
    score_rung,
)
from .version import grader_version

__all__ = ["grade", "GradingSpecError", "resolve"]


class GradingSpecError(ValueError):
    """Raised when a campaign's grading spec is malformed. Always a campaign
    bug, never an agent error."""


_INDEX_RE = re.compile(r"^(.*?)\[(\d+)\]$")


def resolve(root: Mapping[str, Any], path: str) -> Any:
    """Resolve a dotted path with optional list indices: ``a.b[2].c``.

    Raises GradingSpecError on a missing gold path (campaign bug) but returns
    the sentinel ``_MISSING`` for a missing submission path, which scorers treat
    as a wrong answer rather than an error.
    """
    cur: Any = root
    for part in path.split("."):
        if not part:
            raise GradingSpecError(f"empty path segment in {path!r}")
        m = _INDEX_RE.match(part)
        idx = None
        if m:
            part, idx = m.group(1), int(m.group(2))
        if part:
            if not isinstance(cur, Mapping) or part not in cur:
                return _MISSING
            cur = cur[part]
        if idx is not None:
            if not isinstance(cur, (list, tuple)) or idx >= len(cur):
                return _MISSING
            cur = cur[idx]
    return cur


class _Missing:
    __slots__ = ()

    def __repr__(self) -> str:  # pragma: no cover - debug aid only
        return "<MISSING>"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _Missing)

    def __hash__(self) -> int:
        return hash("<MISSING>")


_MISSING = _Missing()


# Scorers that return a float score directly.
_SCORERS = {
    "exact": P.exact,
    "enum_exact": P.enum_exact,
    "set_exact": P.set_exact,
    "set_f1": P.set_f1,
    "multilabel_f1": P.multilabel_f1,
    "interval_score": P.interval_score,
    "allowlist_score": P.allowlist_score,
    "abs_tol_match": P.abs_tol_match,
    "floor_match": P.floor_match,
    "lins_ccc": P.lins_ccc,
    "overlap_f1": P.overlap_f1,
}

# Scorers returning bool, adapted to a 0/1 score.
_BOOL_SCORERS = {
    "inchikey_match": P.inchikey_match,
}


def _bind(arg: Any, submission: Mapping[str, Any], gold: Mapping[str, Any]) -> Any:
    """Bind one scorer argument.

    A string of the form ``pred:<path>`` or ``gold:<path>`` is resolved against
    the submission or the gold respectively. Anything else is a literal taken
    from the spec (a measured tolerance, a vocabulary, an allowlist).
    """
    if isinstance(arg, str):
        if arg.startswith("pred:"):
            return resolve(submission, arg[5:])
        if arg.startswith("gold:"):
            value = resolve(gold, arg[5:])
            if value is _MISSING:
                raise GradingSpecError(f"gold path not found: {arg[5:]!r}")
            if value is None:
                # Gold is never null. A null here is a campaign bug, and
                # silently scoring it would hide a broken gold artifact.
                raise GradingSpecError(f"gold path is null: {arg[5:]!r}")
            return value
    return arg


def _score_component(
    comp: Mapping[str, Any],
    submission: Mapping[str, Any],
    gold: Mapping[str, Any],
) -> float:
    name = comp.get("name")
    scorer_name = comp.get("scorer")
    if not name or not scorer_name:
        raise GradingSpecError(f"component needs 'name' and 'scorer': {comp!r}")

    args = {k: _bind(v, submission, gold) for k, v in sorted(comp.get("args", {}).items())}

    # A missing or null prediction is a wrong answer, scored zero, never an
    # exception -- an agent that omits a field has failed the rung, and the
    # grader must not crash on it. None must be caught alongside _MISSING: a
    # measurer records an omitted numeric field as JSON null, which *resolves*
    # successfully to None, and passing that to a numeric scorer raises and
    # takes down the whole grade. Gold nulls already raised in _bind, so any
    # None reaching here came from the submission.
    if any(v is _MISSING or v is None for v in args.values()):
        return P.rnd(0.0)

    if scorer_name in _BOOL_SCORERS:
        try:
            return P.rnd(1.0 if _BOOL_SCORERS[scorer_name](**args) else 0.0)
        except ValueError as exc:
            raise GradingSpecError(f"component {name!r}: {exc}") from exc

    if scorer_name not in _SCORERS:
        raise GradingSpecError(f"unknown scorer: {scorer_name!r}")
    try:
        return _SCORERS[scorer_name](**args)
    except TypeError as exc:
        raise GradingSpecError(f"component {name!r} bad arguments: {exc}") from exc


def grade(
    submission: Mapping[str, Any],
    gold: Mapping[str, Any],
    spec: Mapping[str, Any],
) -> dict:
    """Grade one campaign submission. Pure.

    ``spec`` is the campaign's grading.json: a rung list, each rung carrying
    components with their scorers, bound arguments and pass threshold.
    """
    rung_specs = spec.get("rungs")
    if not rung_specs:
        raise GradingSpecError("grading spec has no rungs")

    scored = []
    chances = []
    for rs in sorted(rung_specs, key=lambda r: r["ordinal"]):
        comps = rs.get("components")
        if not comps:
            raise GradingSpecError(f"rung {rs.get('rung_id')!r} has no components")
        component_scores = {
            c["name"]: _score_component(c, submission, gold) for c in comps
        }
        scored.append(
            score_rung(
                rung_id=rs["rung_id"],
                ordinal=int(rs["ordinal"]),
                components=component_scores,
                composition=rs.get("composition", "product"),
                pass_threshold=float(rs.get("pass_threshold", 1.0)),
                weights=rs.get("weights"),
            )
        )
        chances.append(float(rs.get("chance_level", 0.0)))

    ladder: LadderResult = score_ladder(scored)
    floor = chance_floor(chances) if any(chances) else 0.0
    disc_floor = (discriminating_chance_floor(chances)
                  if any(chances) and len(chances) >= 2 else 0.0)

    return {
        "campaign_id": spec.get("campaign_id"),
        "grader_version": grader_version(),
        "score": ladder.score,
        "depth": ladder.depth,
        "rungs_cleared": ladder.rungs_cleared,
        "monotonic": ladder.monotonic,
        "chance_floor": floor,
        "discriminating_chance_floor": disc_floor,
        "rungs": [
            {
                "rung_id": r.rung_id,
                "ordinal": r.ordinal,
                "passed": r.passed,
                "composite": r.composite,
                "pass_threshold": r.pass_threshold,
                "composition": r.composition,
                "components": dict(sorted(r.components.items())),
            }
            for r in ladder.rungs
        ],
    }


def grade_json(submission: Mapping, gold: Mapping, spec: Mapping) -> str:
    """Canonical JSON serialisation of a grade, for byte-identity testing."""
    return json.dumps(grade(submission, gold, spec), sort_keys=True, indent=2)
