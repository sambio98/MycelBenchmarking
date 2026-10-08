"""End-to-end behaviour of the NRPS mass-balance campaign."""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.mass_balance_nrps.core import (
    AssemblyError,
    FormulaError,
    MassBalance,
)

CAMPAIGNS = pathlib.Path(__file__).resolve().parents[2] / "campaigns"
CAMPAIGN = CAMPAIGNS / "massbalance-nrps-malleobactin-01"
SEVADICIN = CAMPAIGNS / "massbalance-nrps-sevadicin-01"
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


@pytest.fixture(scope="module")
def mb():
    """The shared template oracle, bound to the malleobactin instantiation.

    Campaigns of this template contain no code: one oracle serves every
    instantiation, parameterised by campaign directory.
    """
    return MassBalance.load(CAMPAIGN)


def _gold(campaign: pathlib.Path = CAMPAIGN):
    return json.loads((campaign / "gold" / "gold.json").read_text())


def _grade_oracle():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    return grade(m, _gold(), s)


# ------------------------------------------------------------------ gates

def test_oracle_scores_exactly_one():
    r = _grade_oracle()
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True


def test_discriminating_chance_floor_within_limit():
    assert _grade_oracle()["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_all_rungs_are_g1():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


# ------------------------------------------------------------------ formula arithmetic

def test_formula_round_trip_is_hill_order(mb):
    F = mb
    # Written out of order, it must come back in Hill order: C, H, then
    # alphabetical. Presentation is normalised so chemistry is what is graded.
    assert F.render(F.parse("O11N6H32C18")) == "C18H32N6O11"
    assert F.render(F.parse("H2O")) == "H2O"


@pytest.mark.parametrize("bad", ["", "   ", "c6h12o6", "C6H12O6!", "Xx4", "C0H2", "6CH"])
def test_formula_parser_rejects_malformed_input(mb, bad):
    with pytest.raises(FormulaError):
        mb.parse(bad)


def test_residual_keeps_its_sign(mb):
    F = mb
    # A deficit must stay negative: dropping the sign would lose the whole
    # point of R3, where a -H2O residual means cyclisation and a positive one
    # means the product carries atoms the modules do not supply.
    assert F.subtract(F.parse("H2O"), F.parse("H4O2")) == {"H": -2, "O": -1}


def test_render_refuses_negative_counts(mb):
    with pytest.raises(FormulaError):
        mb.render({"H": -2})


# ------------------------------------------------------------------ assembly logic

def test_naive_assembly_matches_hand_calculation(mb):
    """Four monomers, three peptide bonds, verified independently."""
    gold = _gold()
    assert gold["n_modules"] == 4
    assert gold["peptide_bonds"] == 3
    total = mb.add(*(mb.monomer_formula(m) for m in gold["monomers"]))
    expected = mb.subtract(total, mb.scale(mb.water, 3))
    assert mb.render(expected) == gold["naive_assembly_formula"] == "C18H32N6O11"


def test_every_annotated_congener_is_analysed():
    """Taking compounds[0] would bake an undeclared choice into gold.

    This entry annotates four congeners, including a dimer, so 'the' product is
    not well defined and a correct answer about another congener must not be
    marked wrong.
    """
    gold = _gold()
    assert len(gold["compounds"]) == 4
    assert gold["compounds"] == sorted(gold["compounds"])
    assert "malleobactin D" in gold["per_compound"]
    # The dimer's residual is far larger than the monomeric congeners'.
    dimer = gold["per_compound"]["malleobactin D"]["residual"]
    mono = gold["per_compound"]["malleobactin A"]["residual"]
    assert dimer["C"] > 2 * mono["C"]


def test_verdict_classification_rules(mb):
    assert mb.classify({}) == "balanced"
    assert mb.classify({"H": -2, "O": -1}) == "cyclisation_only"
    assert mb.classify({"C": 5, "H": 8}) == "additional_chemistry_required"
    assert mb.classify({"H": -4}) == "deficit_unexplained"
    # Mixed signs resolve by the declared precedence, not by argument.
    assert mb.classify({"C": 1, "H": -3}) == "additional_chemistry_required"


def test_deletion_removes_one_monomer_less_one_water(mb):
    """Deleting a module removes that monomer minus the water it no longer
    condenses, so the mass shift is -(monomer - H2O)."""
    gold = _gold()
    water = mb.monoisotopic_mass(mb.water)
    for index, name in enumerate(gold["monomers"], start=1):
        block = gold["single_module_deletions"][str(index)]
        expected = -(mb.monoisotopic_mass(mb.monomer_formula(name)) - water)
        assert block["mass_shift"] == pytest.approx(expected, abs=1e-3)


def test_deletion_is_by_index_not_by_name(mb):
    """A repeated monomer must not collapse two deletions into one."""
    out = mb.single_module_deletions(["serine", "serine", "aspartic acid"])
    assert sorted(out) == ["1", "2", "3"]
    assert out["3"]["deleted_monomer"] == "aspartic acid"


# ------------------------------------------------------------------ isolation & solvability

def test_warm_r4_bundle_withholds_the_deletion_answers(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r4", "gold/warm/r4")
    blob = "\n".join(p.read_text(errors="replace") for p in s.rglob("*") if p.is_file())
    gold = _gold()
    # task.yaml legitimately describes single_module_deletions in the output
    # contract -- the agent must know what to produce -- so the test is on the
    # answer *values*, not on the field name.
    for block in gold["single_module_deletions"].values():
        if block.get("feasible"):
            assert str(block["mass_shift"]) not in blob
            assert block["formula"] not in blob
    # ...while R3's answer, which R4 is entitled to, is present.
    handed = json.loads((s / "analysis.json").read_text())
    assert set(handed["per_compound"]) == set(gold["compounds"])


def test_campaign_is_solvable_from_its_sandbox(tmp_path, mb):
    """Recompute the assembly using only files the agent can see."""
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    mono = json.loads((s / "reference" / "monomer_formulas.json").read_text())["monomers"]
    entry = json.loads((s / "inputs" / "entry.json").read_text())

    names = [m["a_domain"]["substrates"][0]["name"]
             for m in entry["biosynthesis"]["modules"]]
    total = mb.add(*(mb.parse(mono[n]["formula"]) for n in names))
    assembly = mb.subtract(total, mb.scale(mb.water, len(names) - 1))
    assert mb.render(assembly) == _gold()["naive_assembly_formula"]
    assert not (s / "oracle").exists()


# ------------------------------------------------------------ second instantiation

SEVADICIN_READY = (SEVADICIN / "gold" / "gold.json").is_file()
sevadicin_only = pytest.mark.skipif(not SEVADICIN_READY,
                                    reason="sevadicin instantiation not built")


@sevadicin_only
def test_sevadicin_exercises_the_balanced_branch():
    """The template's two instantiations must not test the same verdict.

    Every malleobactin congener is additional_chemistry_required, so that
    instantiation alone never exercises the balanced branch. Sevadicin's
    residual is exactly zero: Phe + Ala + Trp minus two waters is C23H26N4O4,
    which is the product.
    """
    gold = _gold(SEVADICIN)
    assert gold["naive_assembly_formula"] == "C23H26N4O4"
    assert gold["per_compound"]["sevadicin"]["residual"] == {}
    assert gold["per_compound"]["sevadicin"]["verdict"] == "balanced"
    assert {b["verdict"] for b in _gold()["per_compound"].values()} == {
        "additional_chemistry_required"}


@sevadicin_only
def test_both_instantiations_share_one_oracle_and_carry_no_code():
    """The authoring unit is the template: a campaign is data plus a build."""
    for campaign in (CAMPAIGN, SEVADICIN):
        assert not (campaign / "oracle").exists()
        assert (campaign / "reference" / "monomer_formulas.json").is_file()
        assert (campaign / "grading.json").is_file()
        assert not list(campaign.rglob("*.py"))


@sevadicin_only
def test_sevadicin_oracle_scores_one():
    m = json.loads((SEVADICIN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((SEVADICIN / "grading.json").read_text())
    r = grade(m, _gold(SEVADICIN), s)
    assert r["score"] == 1.0 and r["depth"] == 4
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


@sevadicin_only
def test_monomer_tables_differ_between_instantiations():
    """Reference data is per-campaign; only the oracle is shared."""
    a = json.loads((CAMPAIGN / "reference" / "monomer_formulas.json").read_text())
    b = json.loads((SEVADICIN / "reference" / "monomer_formulas.json").read_text())
    assert set(a["monomers"]) != set(b["monomers"])
    assert set(b["monomers"]) == {"alanine", "phenylalanine", "tryptophan"}


def test_assembly_refuses_a_single_monomer(mb):
    with pytest.raises(AssemblyError):
        mb.naive_assembly(["alanine"])
