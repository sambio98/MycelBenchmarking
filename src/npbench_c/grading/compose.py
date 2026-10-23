"""Score composition (NPBench-C spec 8.3).

BixBench3 composes an artifact score as the product of its component metrics and
passes at >= 0.80. Product is correct where components are genuinely conjunctive
-- wrong genes make the right numbers meaningless -- and wrong where they are
independent quality dimensions. Over a deep ladder a bare product collapses
toward zero and destroys discrimination: eighteen terms at 0.9 each score 0.15,
which no agent quality can lift to 0.50.

NPBench-C therefore composes conjunctively *within* a rung and additively
*across* rungs. A rung is the unit inside which all-or-nothing is the honest
semantics.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence

from .primitives import rnd

GEOMETRIC_FLOOR = 0.05

__all__ = ["compose", "GEOMETRIC_FLOOR"]


def _ordered(components: Mapping[str, float]) -> list[tuple[str, float]]:
    """Fix reduction order by component name, so composition is order-stable
    regardless of how the caller's mapping was built."""
    return sorted(components.items())


def compose(
    components: Mapping[str, float],
    mode: str,
    weights: Mapping[str, float] | None = None,
) -> float:
    """Compose component scores into one artifact score.

    mode="product"   conjunctive components; any zero zeroes the artifact
    mode="geometric" independent quality dimensions; weighted geometric mean,
                     floored at GEOMETRIC_FLOOR so one bad dimension cannot
                     annihilate an otherwise-sound artifact
    mode="min"       gate semantics: the artifact is as good as its worst part

    Component-wise scores are always reported alongside the composite
    (spec 8.3); the composite exists for ranking, the components for diagnosis.
    """
    if not components:
        raise ValueError("cannot compose an empty component set")
    items = _ordered(components)
    for name, value in items:
        if not (0.0 <= value <= 1.0):
            raise ValueError(f"component {name!r} out of range [0,1]: {value}")

    if mode == "product":
        out = 1.0
        for _, value in items:
            out *= value
        return rnd(out)

    if mode == "min":
        return rnd(min(value for _, value in items))

    if mode == "geometric":
        w = dict(weights or {})
        ws = [float(w.get(name, 1.0)) for name, _ in items]
        if any(x < 0 for x in ws):
            raise ValueError("weights must be non-negative")
        total = math.fsum(ws)
        if total <= 0:
            raise ValueError("weights must not sum to zero")
        # Floor each component before taking logs so a single zero yields the
        # floor rather than negative infinity.
        acc = math.fsum(
            ws[i] * math.log(max(value, GEOMETRIC_FLOOR))
            for i, (_, value) in enumerate(items)
        )
        return rnd(math.exp(acc / total))

    raise ValueError(f"unknown composition mode: {mode!r}")
