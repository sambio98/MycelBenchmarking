"""Deterministic scoring primitives for NPBench-C.

Purity contract (NPBench-C spec 3.3):
  * no network, no LLM, no wall-clock, no RNG
  * no dependence on set/dict iteration order in any numeric path
  * all floating-point reductions happen in a fixed (sorted) order
  * every returned score is rounded to PRECISION decimals

Only the Python standard library is importable from this module. Third-party
numeric libraries are excluded from the scoring path on purpose: a numpy or
pandas version bump must never be able to move a published score.
"""

from __future__ import annotations

import math
import re
from typing import Hashable, Iterable, Sequence

PRECISION = 6

__all__ = [
    "PRECISION",
    "rnd",
    "exact",
    "enum_exact",
    "set_exact",
    "set_f1",
    "multilabel_f1",
    "interval_score",
    "allowlist_score",
    "inchikey_match",
    "abs_tol_match",
    "floor_match",
    "lins_ccc",
    "overlap_f1",
]


def rnd(x: float) -> float:
    """Round to the reporting precision. Applied to every primitive's output."""
    return round(float(x), PRECISION)


# --------------------------------------------------------------------------
# exact / categorical
# --------------------------------------------------------------------------

def exact(pred: object, gold: object) -> float:
    """1.0 iff pred equals gold under Python equality."""
    return rnd(1.0 if pred == gold else 0.0)


def enum_exact(pred: object, gold: object, vocabulary: Iterable[Hashable]) -> float:
    """1.0 iff pred == gold and pred is a member of the declared vocabulary.

    A prediction outside the vocabulary scores 0.0 rather than raising: an agent
    inventing an enum member is wrong, not a grader error. A *gold* value outside
    the vocabulary is a campaign bug and does raise.
    """
    vocab = set(vocabulary)
    if gold not in vocab:
        raise ValueError(f"gold value {gold!r} is not in the declared vocabulary")
    if pred not in vocab:
        return rnd(0.0)
    return rnd(1.0 if pred == gold else 0.0)


def set_exact(pred: Iterable[Hashable], gold: Iterable[Hashable]) -> float:
    """1.0 iff the two sets are equal. Used where partial credit is meaningless."""
    return rnd(1.0 if set(pred) == set(gold) else 0.0)


def _f1(tp: int, n_pred: int, n_gold: int) -> float:
    if n_pred == 0 and n_gold == 0:
        # Both empty is a correct answer, not an undefined one. This is the
        # zero-confound negative-control case (spec 7.3) and must score 1.0.
        return 1.0
    if tp == 0:
        return 0.0
    precision = tp / n_pred
    recall = tp / n_gold
    return 2.0 * precision * recall / (precision + recall)


def set_f1(pred: Iterable[Hashable], gold: Iterable[Hashable]) -> float:
    """Unweighted F1 over two sets of identifiers."""
    p, g = set(pred), set(gold)
    return rnd(_f1(len(p & g), len(p), len(g)))


def multilabel_f1(
    pred: Iterable[Hashable],
    gold: Iterable[Hashable],
    vocabulary: Iterable[Hashable],
) -> float:
    """F1 over a closed vocabulary (spec 7.3).

    Predictions outside the vocabulary are dropped rather than counted as false
    positives, so that a schema violation is caught by contract validation and
    not silently folded into this score. Gold outside the vocabulary raises.
    """
    vocab = set(vocabulary)
    g = set(gold)
    if not g <= vocab:
        raise ValueError(f"gold labels outside vocabulary: {sorted(g - vocab)!r}")
    p = set(pred) & vocab
    return rnd(_f1(len(p & g), len(p), len(g)))


# --------------------------------------------------------------------------
# numeric
# --------------------------------------------------------------------------

def abs_tol_match(pred: float, gold: float, abs_tol: float) -> float:
    """1.0 iff |pred - gold| <= abs_tol. abs_tol must come from a measured
    determinism report, never from judgement (spec 8.5)."""
    if abs_tol < 0:
        raise ValueError("abs_tol must be non-negative")
    # The difference is rounded to the declared precision before comparison.
    # Without this, |1.0 - 1.1| evaluates to 0.10000000000000009 in binary
    # floating point and a value sitting exactly on the stated tolerance
    # boundary is marked wrong -- a brittle match that contradicts the
    # documented semantics and would read to an auditor as grader noise.
    diff = rnd(abs(float(pred) - float(gold)))
    return rnd(1.0 if diff <= rnd(abs_tol) else 0.0)


def floor_match(pred: float, floor: float) -> float:
    """1.0 iff pred >= floor.

    The honest criterion where many answers are valid and only a quality level
    is required -- a designed construct's codon adaptation, for instance. The
    floor is set below the worst legitimate method measured during construction,
    so an agent is never failed for using a different valid optimiser. Grading
    such a quantity against the oracle's own value would be trajectory lock-in.
    """
    try:
        return rnd(1.0 if float(pred) >= float(floor) else 0.0)
    except (TypeError, ValueError):
        return rnd(0.0)


def interval_score(
    pred_lo: float,
    pred_hi: float,
    gold: float,
    max_width: float,
) -> float:
    """Score a predicted interval against a scalar gold value (spec 7.2).

    1.0  gold inside [pred_lo, pred_hi] and width <= max_width
    0.5  gold inside and max_width < width <= 2 * max_width
    0.0  otherwise

    Penalises wrongness *and* vagueness. ``max_width`` is 3x the observed oracle
    spread from the determinism audit, so both thresholds are measured rather
    than chosen. An inverted interval (lo > hi) scores 0.0 -- it is a malformed
    answer, not a grader error.
    """
    if max_width <= 0:
        raise ValueError("max_width must be positive")
    lo, hi, g = float(pred_lo), float(pred_hi), float(gold)
    if math.isnan(lo) or math.isnan(hi) or math.isinf(lo) or math.isinf(hi):
        return rnd(0.0)
    if lo > hi:
        return rnd(0.0)
    if not (lo <= g <= hi):
        return rnd(0.0)
    width = hi - lo
    if width <= max_width:
        return rnd(1.0)
    if width <= 2.0 * max_width:
        return rnd(0.5)
    return rnd(0.0)


def lins_ccc(pred: Sequence[float], gold: Sequence[float]) -> float:
    """Lin's concordance correlation coefficient, fixed-order reduction.

    Reductions use math.fsum over the input order (which the caller has already
    fixed by sorting on the artifact's declared sort key), so the result is
    bit-stable across platforms.
    """
    if len(pred) != len(gold):
        raise ValueError("pred and gold must have equal length")
    n = len(pred)
    if n == 0:
        raise ValueError("cannot compute CCC on empty vectors")
    if n == 1:
        return rnd(1.0 if float(pred[0]) == float(gold[0]) else 0.0)

    x = [float(v) for v in pred]
    y = [float(v) for v in gold]
    mx = math.fsum(x) / n
    my = math.fsum(y) / n
    vx = math.fsum((xi - mx) ** 2 for xi in x) / n
    vy = math.fsum((yi - my) ** 2 for yi in y) / n
    cov = math.fsum((x[i] - mx) * (y[i] - my) for i in range(n)) / n

    denom = vx + vy + (mx - my) ** 2
    if denom == 0.0:
        # Both vectors constant and identical.
        return rnd(1.0)
    return rnd(2.0 * cov / denom)


# --------------------------------------------------------------------------
# intervals on a coordinate axis (genomic ranges)
# --------------------------------------------------------------------------

def _total_overlap(a: Sequence[tuple], b: Sequence[tuple]) -> int:
    """Total overlapping length between two interval sets. Half-open [start, end).

    Both sides are merged first so that overlapping input intervals cannot
    double-count. Iteration is over sorted tuples, so the result is
    order-independent.
    """
    def merge(ivs):
        out = []
        for s, e in sorted((int(s), int(e)) for s, e in ivs):
            if e <= s:
                continue
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return out

    ma, mb = merge(a), merge(b)
    total = 0
    i = j = 0
    while i < len(ma) and j < len(mb):
        lo = max(ma[i][0], mb[j][0])
        hi = min(ma[i][1], mb[j][1])
        if hi > lo:
            total += hi - lo
        if ma[i][1] < mb[j][1]:
            i += 1
        else:
            j += 1
    return total


def _total_length(ivs: Sequence[tuple]) -> int:
    out = 0
    prev_end = None
    for s, e in sorted((int(s), int(e)) for s, e in ivs):
        if e <= s:
            continue
        s = max(s, prev_end) if prev_end is not None and s < prev_end else s
        if e > s:
            out += e - s
            prev_end = e
    return out


def overlap_f1(pred: Sequence[tuple], gold: Sequence[tuple]) -> float:
    """Base-level F1 between two sets of half-open coordinate intervals.

    Used for BGC boundary calling, where identifier F1 is too coarse (a cluster
    called 2 kb short is not simply 'wrong') and exact match is too brittle.
    """
    inter = _total_overlap(pred, gold)
    lp, lg = _total_length(pred), _total_length(gold)
    return rnd(_f1(inter, lp, lg))


# --------------------------------------------------------------------------
# chemistry
# --------------------------------------------------------------------------

_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def inchikey_match(pred: str, gold: str, block: int = 1) -> bool:
    """Compare InChIKeys at the declared block level (spec 8.1).

    block=1  connectivity skeleton only -- the right granularity when
             stereochemistry is not determinable from the data
    block=2  connectivity + stereochemistry
    block=3  full key including the protonation/version flag

    Chemistry is always keyed on InChIKey, never on canonical SMILES: RDKit's
    canonical SMILES changed at 2022.09, 2023.03, 2023.09 and 2026.03.1, and
    canonicalisation is not guaranteed idempotent.
    """
    if block not in (1, 2, 3):
        raise ValueError("block must be 1, 2 or 3")
    if not _INCHIKEY_RE.match(gold or ""):
        raise ValueError(f"gold is not a well-formed InChIKey: {gold!r}")
    if not _INCHIKEY_RE.match(pred or ""):
        return False
    return pred.split("-")[:block] == gold.split("-")[:block]


# --------------------------------------------------------------------------
# method appropriateness
# --------------------------------------------------------------------------

def allowlist_score(
    pred: Hashable,
    acceptable: Iterable[Hashable],
    unacceptable: Iterable[Hashable],
) -> float:
    """Score a declared method against an oracle-computed allowlist.

    1.0  the method is in ``acceptable``  -- demonstrably recovers gold on this data
    0.0  the method is in ``unacceptable`` -- demonstrably does not
    0.5  the method is in neither         -- untested on this campaign

    This replaces enum-exact matching on ``test_used`` / ``correction``. Grading
    conformity to the oracle's own choice of test would be trajectory lock-in
    (spec 3.4, 20.7): DESeq2 Wald, DESeq2 LRT and limma-voom can all be correct
    for one question. The two lists are *computed*, not curated -- every
    candidate method is executed against gold during campaign construction and
    sorted by whether it reproduces the gold verdict within tolerance. The 0.5
    band is deliberately visible: a campaign accumulating many 0.5s has an
    under-populated candidate set, which the readiness gate reports.
    """
    acc, unacc = set(acceptable), set(unacceptable)
    overlap = acc & unacc
    if overlap:
        raise ValueError(f"method appears in both lists: {sorted(overlap)!r}")
    if not acc:
        raise ValueError("acceptable list is empty; allowlist was never computed")
    if pred in acc:
        return rnd(1.0)
    if pred in unacc:
        return rnd(0.0)
    return rnd(0.5)
