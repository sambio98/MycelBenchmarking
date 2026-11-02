"""Statistics for the internal sweep.

n=3 repeats per system is fixed by project decision, which means the unit of
analysis cannot be the campaign: three observations cannot support a claim.
With a rung ladder the unit becomes the rung-run -- 32 campaigns x 4 rungs x 3
repeats is 384 observations -- and confidence intervals come from a cluster
bootstrap that resamples *campaigns*, which absorbs the correlation between
rungs of one campaign and between repeats of one run.

The resulting resolution is honest and should be pre-registered: NPBench-C
separates systems differing by roughly 8 points or more at n=3, and makes no
claim about finer ordering.

Deterministic: the bootstrap is seeded, so a reported interval is reproducible.
"""

from __future__ import annotations

import math
import random
from typing import Mapping, Sequence

BOOTSTRAP_RESAMPLES = 2000
BOOTSTRAP_SEED = 20261001

# Repeated identical scores do not give an SE of exactly zero: three runs of
# 0.05 leave a residue around 1e-18 from binary floating point. Scores carry six
# decimals, so anything below this is no spread at all, and treating it as real
# would turn a degenerate case into a z of order 1e16.
SE_EPSILON = 1e-9

__all__ = [
    "mean",
    "standard_error",
    "pass_at_1",
    "pass_at_k",
    "cluster_bootstrap_ci",
    "leakage_z",
    "SE_EPSILON",
    "BOOTSTRAP_RESAMPLES",
    "BOOTSTRAP_SEED",
]


def mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("mean of an empty sequence")
    return math.fsum(values) / len(values)


def standard_error(values: Sequence[float]) -> float:
    """Standard error of the mean. Zero for a single observation."""
    n = len(values)
    if n < 2:
        return 0.0
    m = mean(values)
    variance = math.fsum((v - m) ** 2 for v in values) / (n - 1)
    return math.sqrt(variance / n)


def pass_at_1(outcomes: Sequence[bool]) -> float:
    """Fraction of individual runs that cleared."""
    if not outcomes:
        raise ValueError("no outcomes")
    return sum(1 for o in outcomes if o) / len(outcomes)


def pass_at_k(outcomes_by_cell: Sequence[Sequence[bool]]) -> float:
    """Fraction of cells where *every* repeat cleared (Pass^k).

    Reported alongside Pass@1: the gap between them is the reliability
    measurement, which at n=3 is the only variance story available without
    paying for more repeats.
    """
    if not outcomes_by_cell:
        raise ValueError("no cells")
    return sum(1 for cell in outcomes_by_cell if cell and all(cell)) / len(outcomes_by_cell)


def cluster_bootstrap_ci(
    clusters: Mapping[str, Sequence[float]],
    alpha: float = 0.05,
    resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[float, float, float]:
    """Percentile CI for a mean, resampling whole clusters with replacement.

    ``clusters`` maps a cluster key (campaign id, or template id once campaigns
    are instantiated from templates) to that cluster's observations. Resampling
    clusters rather than observations is what keeps the interval honest: four
    instantiations of one template share an oracle and a failure mode, so they
    are not four independent campaigns.

    Returns (point_estimate, ci_low, ci_high).
    """
    keys = sorted(clusters)
    if not keys:
        raise ValueError("no clusters")
    for k in keys:
        if not clusters[k]:
            raise ValueError(f"cluster {k!r} has no observations")

    def grand_mean(selected: Sequence[str]) -> float:
        # Mean of cluster means: each campaign weighs equally regardless of how
        # many runs it happened to contribute.
        return mean([mean(list(clusters[k])) for k in selected])

    point = grand_mean(keys)
    if len(keys) == 1:
        return point, point, point

    rng = random.Random(seed)
    draws = []
    for _ in range(resamples):
        sample = [keys[rng.randrange(len(keys))] for _ in keys]
        draws.append(grand_mean(sample))
    draws.sort()
    lo_idx = int(math.floor((alpha / 2) * (resamples - 1)))
    hi_idx = int(math.ceil((1 - alpha / 2) * (resamples - 1)))
    return point, draws[lo_idx], draws[hi_idx]


def leakage_z(leakage_scores: Sequence[float], chance: float) -> tuple[float, float, float]:
    """How far a no-tool ablation system sits above the chance floor.

    Spec 14 kills a campaign only when leakage exceeds chance by more than two
    standard errors -- a binary kill on any above-chance performance was too
    strict and removed good tasks that merely have legitimate priors.

    Returns (leakage_mean, standard_error, z). z is +inf when the estimate is
    above chance with zero observed spread, and -inf when at or below chance
    with zero spread, so the caller never has to special-case a degenerate SE.
    """
    m = mean(leakage_scores)
    se = standard_error(leakage_scores)
    if se < SE_EPSILON:
        return m, 0.0, math.inf if m > chance else -math.inf
    return m, se, (m - chance) / se
