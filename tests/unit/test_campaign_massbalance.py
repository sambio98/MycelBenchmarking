"""End-to-end behaviour of the NRPS mass-balance campaign."""

from __future__ import annotations

import json
import pathlib
import sys

import pytest

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "massbalance-nrps-malleobactin-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


@pytest.fixture(scope="module")
def oracle_modules():
    sys.path.insert(0, str(CAMPAIGN / "oracle"))
    import assembly
    import formula
    return formula, assembly


def _gold():
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


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

def test_formula_round_trip_is_hill_order(oracle_modules):
    F, _ = oracle_modules
    # Written out of order, it must come back in Hill order: C, H, then
    # alphabetical. Presentation is normalised so chemistry is what is graded.
    assert F.render(F.parse("O11N6H32C18")) == "C18H32N6O11"
    assert F.render(F.parse("H2O")) == "H2O"


@pytest.mark.parametrize("bad", ["", "   ", "c6h12o6", "C6H12O6!", "Xx4", "C0H2", "6CH"])
def test_formula_parser_rejects_malformed_input(oracle_modules, bad):
    F, _ = oracle_modules
    with pytest.raises(F.FormulaError):
        F.parse(bad)


def test_residual_keeps_its_sign(oracle_modules):
    F, _ = oracle_modules
    # A deficit must stay negative: dropping the sign would lose the whole
    # point of R3, where a -H2O residual means cyclisation and a positive one
    # means the product carries atoms the modules do not supply.
    assert F.subtract(F.parse("H2O"), F.parse("H4O2")) == {"H": -2, "O": -1}


def test_render_refuses_negative_counts(oracle_modules):
    F, _ = oracle_modules
    with pytest.raises(F.FormulaError):
        F.render({"H": -2})


# ------------------------------------------------------------------ assembly logic

def test_naive_assembly_matches_hand_calculation(oracle_modules):
    """Four monomers, three peptide bonds, verified independently."""
    F, A = oracle_modules
    gold = _gold()
    assert gold["n_modules"] == 4
    assert gold["peptide_bonds"] == 3
    total = F.add(*(A.monomer_formula(m) for m in gold["monomers"]))
    expected = F.subtract(total, F.scale(F.WATER, 3))
    assert F.render(expected) == gold["naive_assembly_formula"] == "C18H32N6O11"


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


def test_verdict_classification_rules(oracle_modules):
    _, A = oracle_modules
    assert A.classify({}) == "balanced"
    assert A.classify({"H": -2, "O": -1}) == "cyclisation_only"
    assert A.classify({"C": 5, "H": 8}) == "additional_chemistry_required"
    assert A.classify({"H": -4}) == "deficit_unexplained"
    # Mixed signs resolve by the declared precedence, not by argument.
    assert A.classify({"C": 1, "H": -3}) == "additional_chemistry_required"


def test_deletion_removes_one_monomer_less_one_water(oracle_modules):
    """Deleting a module removes that monomer minus the water it no longer
    condenses, so the mass shift is -(monomer - H2O)."""
    F, A = oracle_modules
    gold = _gold()
    water = F.monoisotopic_mass(F.WATER)
    for index, name in enumerate(gold["monomers"], start=1):
        block = gold["single_module_deletions"][str(index)]
        expected = -(F.monoisotopic_mass(A.monomer_formula(name)) - water)
        assert block["mass_shift"] == pytest.approx(expected, abs=1e-3)


def test_deletion_is_by_index_not_by_name(oracle_modules):
    """A repeated monomer must not collapse two deletions into one."""
    _, A = oracle_modules
    out = A.single_module_deletions(["serine", "serine", "aspartic acid"])
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


def test_campaign_is_solvable_from_its_sandbox(tmp_path, oracle_modules):
    """Recompute the assembly using only files the agent can see."""
    F, _ = oracle_modules
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    mono = json.loads((s / "reference" / "monomer_formulas.json").read_text())["monomers"]
    entry = json.loads((s / "inputs" / "entry.json").read_text())

    names = [m["a_domain"]["substrates"][0]["name"]
             for m in entry["biosynthesis"]["modules"]]
    total = F.add(*(F.parse(mono[n]["formula"]) for n in names))
    assembly = F.subtract(total, F.scale(F.WATER, len(names) - 1))
    assert F.render(assembly) == _gold()["naive_assembly_formula"]
    assert not (s / "oracle").exists()
