"""Rung ladders (NPBench-C spec 8.3a).

A campaign is scored as an ordered ladder of rungs rather than as one scalar.
The canonical shape generalises across every NP/BGC family:

  R1 Execute         toolchain runs; output contract honoured
  R2 Correct         primary artifact matches gold within measured tolerance
  R3 Commit          derived claim or decision that depends on R2
  R4 Counterfactual  answer under a perturbation that makes the published
                     answer wrong

The headline campaign score is the highest *contiguous* rung cleared, divided by
the rung count. Contiguity matters: without it an agent can be credited for
clearing R4 by luck while failing R2, which is the signature of a guessed
answer rather than a computed one. ``rungs_cleared`` is reported alongside as a
diagnostic, so monotonicity violations surface instead of hiding.

Each rung is independently runnable in two modes:
  cold -- the agent receives the campaign objective and attempts every rung
  warm -- the agent starts from gold artifacts of rung k-1 and attempts rung k

Warm mode is what makes a rung a standalone (prompt, gold) pair, and what gives
a per-rung difficulty profile; in cold mode a failure at R2 censors everything
above it, so cold runs alone cannot calibrate a ladder.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping, Sequence

from .compose import compose
from .primitives import rnd

MAX_CHANCE_FLOOR = 0.10
MIN_R1_CLEAR_RATE = 0.80

__all__ = [
    "RungResult",
    "LadderResult",
    "score_rung",
    "score_ladder",
    "chance_floor",
    "discriminating_chance_floor",
    "check_monotonicity",
    "MAX_CHANCE_FLOOR",
    "MIN_R1_CLEAR_RATE",
]


@dataclass(frozen=True)
class RungResult:
    rung_id: str
    ordinal: int
    components: Mapping[str, float]
    composite: float
    passed: bool
    pass_threshold: float
    composition: str


@dataclass(frozen=True)
class LadderResult:
    rungs: Sequence[RungResult]
    depth: int                 # highest contiguous rung cleared
    rungs_cleared: int         # total cleared, contiguous or not
    score: float               # depth / n_rungs
    monotonic: bool            # no rung cleared above a failed one
    chance_floor: float
    detail: Mapping[str, object] = field(default_factory=dict)


def score_rung(
    rung_id: str,
    ordinal: int,
    components: Mapping[str, float],
    composition: str = "product",
    pass_threshold: float = 1.0,
    weights: Mapping[str, float] | None = None,
) -> RungResult:
    """Score one rung. Components compose conjunctively by default."""
    composite = compose(components, composition, weights)
    return RungResult(
        rung_id=rung_id,
        ordinal=ordinal,
        components=dict(sorted(components.items())),
        composite=composite,
        passed=composite >= pass_threshold,
        pass_threshold=pass_threshold,
        composition=composition,
    )


def score_ladder(rungs: Sequence[RungResult]) -> LadderResult:
    """Reduce scored rungs to a campaign score.

    Rungs are sorted by ordinal first, so the result does not depend on the
    order in which the caller happened to score them.
    """
    if not rungs:
        raise ValueError("a ladder must have at least one rung")
    ordered = sorted(rungs, key=lambda r: r.ordinal)
    ordinals = [r.ordinal for r in ordered]
    if ordinals != list(range(1, len(ordered) + 1)):
        raise ValueError(f"rung ordinals must be 1..n with no gaps, got {ordinals}")

    depth = 0
    for r in ordered:
        if r.passed:
            depth += 1
        else:
            break

    cleared = sum(1 for r in ordered if r.passed)
    return LadderResult(
        rungs=tuple(ordered),
        depth=depth,
        rungs_cleared=cleared,
        score=rnd(depth / len(ordered)),
        monotonic=(cleared == depth),
        chance_floor=0.0,
        detail={"n_rungs": len(ordered)},
    )


def chance_floor(rung_chance: Sequence[float]) -> float:
    """Expected ladder score for an agent that guesses every rung.

    With independent per-rung guess probabilities p_k, the contiguous depth
    satisfies P(depth >= k) = prod_{j<=k} p_j, so

        E[score] = (1/n) * sum_k prod_{j<=k} p_j

    A campaign whose chance floor exceeds MAX_CHANCE_FLOOR hands out free
    credit and loses control of the headline; the readiness gate fails it. This
    is the rule that most protects the benchmark's calibration: a four-rung
    ladder whose R3 is a three-way enum is not a hard campaign no matter how
    hard R1 and R2 are.
    """
    if not rung_chance:
        raise ValueError("need at least one rung chance level")
    for p in rung_chance:
        if not (0.0 <= p <= 1.0):
            raise ValueError(f"chance level out of range [0,1]: {p}")
    n = len(rung_chance)
    cumulative = []
    acc = 1.0
    for p in rung_chance:
        acc *= p
        cumulative.append(acc)
    return rnd(math.fsum(cumulative) / n)


def discriminating_chance_floor(rung_chance: Sequence[float]) -> float:
    """Chance floor over the *discriminating* rungs only -- R2 upward.

    R1 is an execution precondition, not a measurement: its design target is an
    85-95% clear rate, and clearing it without tools is often legitimate (for a
    construct-design campaign, reverse-translating a protein needs the genetic
    code and nothing else). So a four-rung ladder hands any no-tool agent 0.25
    for free, and a cap of 0.10 on the *full* floor is unsatisfiable by
    construction -- it would condemn every ladder built to the documented shape.

    The cap therefore applies here, to the rungs that are supposed to separate
    systems. The full floor still matters, as the baseline the no-tool leakage
    ablation is compared against.
    """
    if len(rung_chance) < 2:
        raise ValueError("need at least two rungs to have a discriminating subset")
    return chance_floor(list(rung_chance)[1:])


def check_monotonicity(clear_rates: Sequence[float], tolerance: float = 0.0) -> bool:
    """Verify that measured clear rates are non-increasing in rung ordinal.

    Called by the readiness gate against internal agent runs, never at scoring
    time. A violation means the ladder is mis-ordered -- a campaign bug caught
    mechanically rather than by a reviewer counting.
    """
    for a, b in zip(clear_rates, clear_rates[1:]):
        if b > a + tolerance:
            return False
    return True
