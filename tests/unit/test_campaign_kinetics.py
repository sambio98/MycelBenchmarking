"""BRENDA kinetics: a corpus checked against itself, with the units implicit.

Km, kcat and kcat/Km are stored independently and the third is determined by the
first two, so no external standard and no curator is needed -- the arithmetic is
the gold. What makes it a task rather than a division is that BRENDA records no
units: they belong to the field, and a system that assumes the wrong one is wrong
by exactly the factor the campaign measures.

These tests pin the parsing rules that a careless pass would get wrong, the key
that makes two measurements comparable, and the band order that would otherwise
silently empty the tighter band.
"""

from __future__ import annotations

import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.kinetics_consistency import analyse as kin
from npbench_c.templates.kinetics_consistency.build import _assert_has_content
from npbench_c.templates.kinetics_consistency.core import (
    FIELDS,
    Band,
    Kinetics,
    KineticsError,
    parse_value,
)

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "kinetics-brenda-ec1_1-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def kinetics():
    params = kin.load_params(CAMPAIGN)
    return Kinetics.load(CAMPAIGN, params["kinetics_input"])


# ---------------------------------------------------------------- the gates

def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_all_rungs_are_g1():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


# --------------------------------------------------------------- the parser

@pytest.mark.parametrize("raw,expected", [
    ("0.05 {benzyl alcohol}", (0.05, "benzyl alcohol", "parsed")),
    ("1e-3 {X}", (0.001, "X", "parsed")),
    ("12", (12.0, "", "parsed")),
    # The sentinel, in every form it is written. Numeric comparison, not a string
    # prefix: -9990 is a measurement, -9.99e2 is the sentinel.
    ("-999", (None, "", "sentinel")),
    ("-999.0", (None, "", "sentinel")),
    ("-9.99e2 {X}", (None, "X", "sentinel")),
    ("-9990 {X}", (-9990.0, "X", "parsed")),
    # Not coerced into a number nobody wrote down.
    ("0.05 - 0.1 {X}", (None, "", "unparsed")),
    ("", (None, "", "unparsed")),
])
def test_value_parsing(raw, expected):
    assert parse_value(raw, "-999") == expected


def test_the_sentinel_is_present_and_excluded(kinetics):
    """BRENDA writes -999 for a missing measurement: large, negative, and in
    thousands of records. Counting it as a number is the single most damaging
    thing a naive pass can do here."""
    gold = _gold()
    assert sum(gold["sentinel_by_field"].values()) > 100
    for field in FIELDS:
        assert (gold["records_by_field"][field]
                == gold["parsed_by_field"][field]
                + gold["sentinel_by_field"][field]
                + gold["unparsed_by_field"][field])
    # Nothing negative survives into the triples, which is what the exclusion is
    # for -- a kinetic constant is a positive quantity.
    for values in kinetics.triples().values():
        assert all(v > -999 for v in values.values())


def test_ranges_are_counted_not_averaged():
    assert sum(_gold()["unparsed_by_field"].values()) > 0


# ------------------------------------------------------------------ the key

def test_the_triple_key_is_tight(kinetics):
    """Dropping the reference would pair values from different papers; dropping
    the protein would pair different enzymes. Either way the disagreement rate
    stops meaning anything, so the loose key must give a different -- larger --
    set of complete triples."""
    tight = {k: v for k, v in kinetics.triples().items() if len(v) == 3}
    loose: dict = {}
    for key, values in kinetics.triples().items():
        loose.setdefault((key[0], key[2]), {}).update(values)
    loose_complete = {k: v for k, v in loose.items() if len(v) == 3}
    assert len(loose_complete) > 0
    assert len(tight) != len(loose_complete), (
        "the loose key gives the same answer, so the campaign's counting rule is "
        "not doing the work it claims")


# ---------------------------------------------------------------- the bands

def test_band_order_matters_and_is_declared(kinetics):
    """'near' contains 'agree'. Testing them in the other order would empty the
    tighter band, so the first-match-wins order is part of the rule."""
    names = [b.name for b in kinetics.rules.bands]
    assert names.index("agree") < names.index("near")
    assert kinetics.rules.band_of(1.0) == "agree"
    assert kinetics.rules.band_of(1.2) == "near"


def test_a_reciprocal_band_matches_both_directions():
    band = Band("factor_1000", 900.0, 1100.0, reciprocal=True)
    assert band.matches(1000.0) and band.matches(1 / 1000.0)
    plain = Band("x", 900.0, 1100.0, reciprocal=False)
    assert plain.matches(1000.0) and not plain.matches(1 / 1000.0)


def test_the_unit_mixup_bands_have_real_instances():
    """The campaign's claim is that the implicit convention is load-bearing. If
    no triple landed at a factor of 60 or 1000 the claim would have no evidence
    in this corpus, and the build refuses in that case."""
    bands = _gold()["band_counts"]
    assert bands["factor_60"] > 0 and bands["factor_1000"] > 0
    assert _gold()["unit_mixup_total"] == bands["factor_60"] + bands["factor_1000"]
    # Most of the corpus agrees, which is what makes the exceptions identifiable.
    assert _gold()["agreement_rate"] > 0.5


def test_the_examples_show_both_directions():
    examples = _gold()["band_examples"]["factor_1000"]
    assert examples
    ratios = [e["ratio"] for e in examples]
    assert all(r > 500 or r < 0.002 for r in ratios)


# -------------------------------------------------------- the counterfactual

def test_the_identity_variant_reproduces_the_headline():
    """A control: a system that computes R4 by a different route than R3 will
    disagree with itself here before it disagrees with gold anywhere else."""
    gold = _gold()
    assert (gold["unit_error_variants"]["v0_identity"]["agreement_rate"]
            == gold["agreement_rate"])


def test_a_unit_mistake_collapses_the_agreement_rate():
    gold = _gold()
    baseline = gold["agreement_rate"]
    for name, rate in gold["unit_error_rates"].items():
        if name == "v0_identity":
            continue
        assert rate < baseline / 10, (name, rate, baseline)


# -------------------------------------------------------------- build refusal

@pytest.mark.parametrize("mutate,message", [
    (lambda g: g["sentinel_by_field"].update({f: 0 for f in FIELDS}), "sentinel"),
    (lambda g: g["completeness_counts"].__setitem__("3", 3), "complete triples"),
    (lambda g: g.__setitem__("unit_mixup_total", 0), "unit-mix-up band"),
    (lambda g: g["unit_error_variants"]["v0_identity"].__setitem__(
        "agreement_rate", 0.1), "identity variant"),
])
def test_build_refuses_gold_without_content(mutate, message):
    doctored = json.loads(json.dumps(_gold()))
    mutate(doctored)
    with pytest.raises(KineticsError, match=message):
        _assert_has_content(doctored)


# --------------------------------------------------- sandbox self-containment

def test_sandbox_carries_no_gold(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "cold")
    assert sorted(p.name for p in s.iterdir()) == ["inputs", "reference",
                                                   "submission", "task.yaml"]


def test_campaign_is_solvable_from_its_sandbox(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    recomputed = kin.analyse(s)
    gold = _gold()
    for key in ("records_by_field", "sentinel_by_field", "completeness_counts",
                "band_counts", "agreement_rate", "unit_error_variants"):
        assert recomputed[key] == gold[key], key


def test_the_catalog_deviation_and_the_credentials_correction_are_recorded():
    """Two things this campaign overturned: the catalogued task, and the
    catalogued reason it was blocked."""
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    deviation = task["catalog_deviation"]
    assert "harmonisation across unit conventions" in deviation
    assert "credentials" in deviation
    provenance = json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())
    assert provenance["license"] == "CC BY 4.0"
    assert "not a registration wall" in provenance["access"]["gate"].lower() or \
        "NOT a registration wall" in provenance["access"]["gate"]
    assert provenance["access"]["acceptance_recorded"]
