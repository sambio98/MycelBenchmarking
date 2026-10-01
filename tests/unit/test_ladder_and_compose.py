"""Ladder and composition behaviour."""

from __future__ import annotations

import itertools

import pytest

from npbench_c.grading.compose import GEOMETRIC_FLOOR, compose
from npbench_c.grading.ladder import (
    MAX_CHANCE_FLOOR,
    chance_floor,
    discriminating_chance_floor,
    check_monotonicity,
    score_ladder,
    score_rung,
)


def _rungs(passes):
    return [
        score_rung(f"r{i+1}", i + 1, {"c": 1.0 if ok else 0.0})
        for i, ok in enumerate(passes)
    ]


def test_depth_is_contiguous_not_a_count():
    """Clearing R4 while failing R2 must not earn credit for R4.

    Non-contiguous success is the signature of a guessed answer rather than a
    computed one, so depth stops at the first failure and the discrepancy is
    surfaced via monotonic=False.
    """
    result = score_ladder(_rungs([True, False, True, True]))
    assert result.depth == 1
    assert result.rungs_cleared == 3
    assert result.monotonic is False
    assert result.score == 0.25


def test_all_ladder_outcomes_enumerated():
    for passes in itertools.product([True, False], repeat=4):
        r = score_ladder(_rungs(passes))
        expected_depth = 0
        for ok in passes:
            if ok:
                expected_depth += 1
            else:
                break
        assert r.depth == expected_depth
        assert r.score == round(expected_depth / 4, 6)
        assert r.monotonic == (r.rungs_cleared == r.depth)


def test_ladder_is_order_independent():
    rungs = _rungs([True, True, False, True])
    base = score_ladder(rungs)
    for perm in itertools.permutations(rungs):
        assert score_ladder(list(perm)).score == base.score


def test_ladder_rejects_gapped_ordinals():
    bad = [score_rung("r1", 1, {"c": 1.0}), score_rung("r3", 3, {"c": 1.0})]
    with pytest.raises(ValueError):
        score_ladder(bad)


def test_chance_floor_catches_a_three_way_enum_rung():
    """A four-rung ladder whose R3 is a 3-way enum hands out free credit."""
    loose = chance_floor([0.5, 0.1, 1 / 3, 0.1])
    tight = chance_floor([0.3, 0.02, 0.05, 0.02])
    assert loose > MAX_CHANCE_FLOOR
    assert tight <= MAX_CHANCE_FLOOR


def test_chance_floor_is_monotone_in_guessability():
    assert chance_floor([0.9, 0.9, 0.9, 0.9]) > chance_floor([0.9, 0.9, 0.9, 0.1])


def test_discriminating_floor_excludes_an_easy_r1():
    """An execution rung any no-tool system clears must not consume the cap."""
    levels = [1.0, 0.02, 0.001, 0.001]
    assert chance_floor(levels) > MAX_CHANCE_FLOOR          # full floor: ~0.26
    assert discriminating_chance_floor(levels) <= MAX_CHANCE_FLOOR


def test_discriminating_floor_still_catches_a_small_enum():
    """Excluding R1 must not let a guessable R3 through."""
    assert discriminating_chance_floor([1.0, 0.5, 1 / 3, 0.5]) > MAX_CHANCE_FLOOR


def test_discriminating_floor_needs_two_rungs():
    with pytest.raises(ValueError):
        discriminating_chance_floor([0.5])


def test_monotonicity_check():
    assert check_monotonicity([0.9, 0.7, 0.4, 0.15]) is True
    assert check_monotonicity([0.9, 0.4, 0.7, 0.15]) is False


def test_product_collapses_over_depth():
    """The quantitative reason a bare product cannot reach a 0.50 headline.

    Eighteen conjunctive terms at 0.9 each -- a genuinely competent agent --
    compose to 0.15. No agent quality lifts that to half marks, which is why
    composition is conjunctive within a rung and additive across rungs.
    """
    assert compose({f"c{i:02d}": 0.9 for i in range(18)}, "product") < 0.16
    assert compose({f"c{i:02d}": 0.9 for i in range(6)}, "product") > 0.53


def test_geometric_floors_a_single_zero():
    assert compose({"a": 1.0, "b": 1.0, "c": 0.0}, "geometric") > GEOMETRIC_FLOOR
    assert compose({"a": 1.0, "b": 1.0, "c": 0.0}, "product") == 0.0


def test_compose_is_mapping_order_independent():
    forward = {"a": 0.5, "b": 0.7, "c": 0.9}
    reverse = {"c": 0.9, "b": 0.7, "a": 0.5}
    for mode in ("product", "geometric", "min"):
        assert compose(forward, mode) == compose(reverse, mode)


def test_compose_rejects_out_of_range_and_empty():
    with pytest.raises(ValueError):
        compose({"a": 1.5}, "product")
    with pytest.raises(ValueError):
        compose({}, "product")
    with pytest.raises(ValueError):
        compose({"a": 1.0}, "harmonic")
