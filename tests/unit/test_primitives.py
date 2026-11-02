"""Property tests for the scoring primitives.

Case generation is deterministic (seeded or enumerated) -- a test suite with an
unseeded RNG cannot support a claim about a pure grader.
"""

from __future__ import annotations

import itertools
import math
import random

import pytest

from npbench_c.grading import primitives as P


# ---------------------------------------------------------------- interval_score

@pytest.mark.parametrize("lo,hi,gold,width,expected", [
    (1.0, 2.0, 1.5, 1.0, 1.0),      # contained, width at limit
    (1.0, 2.0, 1.5, 0.5, 0.5),      # contained, width 2x limit
    (1.0, 2.0, 1.5, 0.4, 0.0),      # contained, too vague
    (1.0, 2.0, 2.5, 1.0, 0.0),      # not contained
    (1.5, 1.5, 1.5, 0.1, 1.0),      # point interval, exactly right
    (2.0, 1.0, 1.5, 1.0, 0.0),      # inverted interval is malformed
    (1.0, 2.0, 1.0, 1.0, 1.0),      # gold on the lower boundary
    (1.0, 2.0, 2.0, 1.0, 1.0),      # gold on the upper boundary
])
def test_interval_score_cases(lo, hi, gold, width, expected):
    assert P.interval_score(lo, hi, gold, width) == expected


def test_interval_score_rejects_unbounded_hedge():
    """The whole point of the primitive: [-inf, inf] must not score."""
    assert P.interval_score(-math.inf, math.inf, 0.0, 1.0) == 0.0
    assert P.interval_score(-1e308, 1e308, 0.0, 1.0) == 0.0


def test_interval_score_monotone_in_width():
    """Widening an interval can never raise the score."""
    rng = random.Random(20260101)
    for _ in range(300):
        gold = rng.uniform(-100, 100)
        centre = gold + rng.uniform(-1, 1)
        scores = []
        for half in (0.05, 0.2, 0.5, 1.0, 2.0, 5.0):
            scores.append(P.interval_score(centre - half, centre + half, gold, 1.0))
        contained = [s for s in scores if s > 0]
        assert contained == sorted(contained, reverse=True)


def test_interval_score_rejects_nonpositive_width():
    with pytest.raises(ValueError):
        P.interval_score(0.0, 1.0, 0.5, 0.0)


# ---------------------------------------------------------------- f1 family

def test_empty_sets_score_one():
    """A correct empty answer is correct -- the zero-confound control case."""
    assert P.set_f1([], []) == 1.0
    assert P.multilabel_f1([], [], {"a", "b"}) == 1.0


def test_multilabel_f1_rejects_hallucinated_labels():
    vocab = {"batch_effect", "pseudoreplication", "underpowered"}
    # Finding the real one and inventing another costs precision.
    score = P.multilabel_f1({"batch_effect", "underpowered"}, {"batch_effect"}, vocab)
    # Compared at the grader's declared precision, not tighter: every score is
    # rounded to P.PRECISION before it is ever reported or compared.
    assert score == round(2 / 3, P.PRECISION)


def test_multilabel_f1_drops_out_of_vocab_predictions():
    vocab = {"a", "b"}
    assert P.multilabel_f1({"a", "zzz"}, {"a"}, vocab) == 1.0


def test_multilabel_f1_raises_on_bad_gold():
    with pytest.raises(ValueError):
        P.multilabel_f1({"a"}, {"not_in_vocab"}, {"a", "b"})


def test_set_f1_exhaustive_small_universe():
    """All 2^4 x 2^4 subset pairs: F1 is symmetric and bounded."""
    universe = ["a", "b", "c", "d"]
    subsets = [set(c) for r in range(5) for c in itertools.combinations(universe, r)]
    for p in subsets:
        for g in subsets:
            s = P.set_f1(p, g)
            assert 0.0 <= s <= 1.0
            assert s == P.set_f1(g, p)
            assert (s == 1.0) == (p == g)


# ---------------------------------------------------------------- floor / tolerance

@pytest.mark.parametrize("pred,floor,expected", [
    (0.95, 0.90, 1.0), (0.90, 0.90, 1.0), (0.89, 0.90, 0.0), (None, 0.9, 0.0),
])
def test_floor_match(pred, floor, expected):
    assert P.floor_match(pred, floor) == expected


def test_abs_tol_match_boundary():
    """A value exactly on the tolerance boundary must pass.

    Regression test: |1.0 - 1.1| is 0.10000000000000009 in binary floating
    point, so an unrounded comparison fails a boundary case that the contract
    says should pass.
    """
    assert P.abs_tol_match(1.0, 1.1, 0.1) == 1.0
    assert P.abs_tol_match(1.0, 1.2, 0.1) == 0.0
    assert P.abs_tol_match(0.3, 0.1 + 0.2, 0.0) == 1.0


# ---------------------------------------------------------------- allowlist

def test_allowlist_score_bands():
    acc, unacc = {"deseq2_wald", "deseq2_lrt"}, {"ttest_on_raw_counts"}
    assert P.allowlist_score("deseq2_lrt", acc, unacc) == 1.0
    assert P.allowlist_score("ttest_on_raw_counts", acc, unacc) == 0.0
    assert P.allowlist_score("limma_voom", acc, unacc) == 0.5


def test_allowlist_rejects_contradictory_or_empty_lists():
    with pytest.raises(ValueError):
        P.allowlist_score("x", {"a"}, {"a"})
    with pytest.raises(ValueError):
        P.allowlist_score("x", set(), {"a"})


# ---------------------------------------------------------------- CCC

def test_ccc_identity_and_shift():
    v = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert P.lins_ccc(v, v) == 1.0
    assert P.lins_ccc([x + 10 for x in v], v) < 1.0


def test_ccc_reduction_order_is_fixed():
    """fsum over a fixed order: the same vector pair always gives the same bits."""
    rng = random.Random(7)
    a = [rng.uniform(-1e6, 1e6) for _ in range(500)]
    b = [rng.uniform(-1e6, 1e6) for _ in range(500)]
    assert len({P.lins_ccc(a, b) for _ in range(50)}) == 1


def test_ccc_rejects_mismatched_length():
    with pytest.raises(ValueError):
        P.lins_ccc([1.0], [1.0, 2.0])


# ---------------------------------------------------------------- intervals on an axis

def test_overlap_f1_cases():
    assert P.overlap_f1([(0, 100)], [(0, 100)]) == 1.0
    assert P.overlap_f1([(0, 100)], [(200, 300)]) == 0.0
    assert P.overlap_f1([(0, 100)], [(0, 200)]) == round(2 / 3, P.PRECISION)


def test_overlap_f1_merges_overlapping_inputs():
    """Double-reporting a region must not inflate the score."""
    assert P.overlap_f1([(0, 60), (40, 100)], [(0, 100)]) == 1.0


def test_overlap_f1_order_invariant():
    pred = [(500, 900), (0, 100), (200, 300)]
    gold = [(250, 600), (50, 120)]
    base = P.overlap_f1(pred, gold)
    for _ in range(20):
        random.Random(1).shuffle(pred)
        assert P.overlap_f1(pred, gold) == base


# ---------------------------------------------------------------- chemistry

def test_inchikey_block_levels():
    a = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    b = "RYYVLZVUVIJVGH-ABCDEFGHIJ-N"
    assert P.inchikey_match(a, b, block=1) is True
    assert P.inchikey_match(a, b, block=2) is False
    assert P.inchikey_match(a, a, block=3) is True


def test_inchikey_malformed_prediction_is_wrong_not_an_error():
    assert P.inchikey_match("CCO", "RYYVLZVUVIJVGH-UHFFFAOYSA-N") is False


def test_inchikey_malformed_gold_raises():
    with pytest.raises(ValueError):
        P.inchikey_match("RYYVLZVUVIJVGH-UHFFFAOYSA-N", "not-a-key")


# ---------------------------------------------------------------- precision

def test_every_primitive_rounds_to_declared_precision():
    rng = random.Random(99)
    for _ in range(200):
        a = [rng.uniform(0, 1) for _ in range(9)]
        b = [rng.uniform(0, 1) for _ in range(9)]
        for value in (P.lins_ccc(a, b), P.set_f1(set("abc"), set("bcd"))):
            assert value == round(value, P.PRECISION)
