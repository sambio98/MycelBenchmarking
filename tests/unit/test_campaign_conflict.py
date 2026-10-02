"""The constraint-conflict campaign: abstention gold that is provable.

The property under test is not "our optimiser could not find a design" -- that
proves nothing -- but that NO sequence encoding the protein can satisfy the
constraint set, by exact arithmetic over per-residue codon choices.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.construct_design.core import AA_CODONS, ConstructEngine, gc_count
from npbench_c.templates.construct_design.params import (
    CONFLICT_CONSTRAINTS,
    VERDICTS,
    ConflictSpec,
    read_protein,
)

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "constraint-conflict-ecoli-rebh-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _gold():
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def engine():
    return ConstructEngine.load(CAMPAIGN)


# ---------------------------------------------------------------- gates

def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_all_rungs_are_g1():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


# ------------------------------------------------- the infeasibility proof

def test_the_conflict_is_genuinely_pairwise(engine):
    """Each constraint alone must be satisfiable, or it is not a pair conflict."""
    gold = _gold()
    spec = ConflictSpec.load(CAMPAIGN)
    protein = read_protein(CAMPAIGN, spec.protein_input)

    # forbidden_codons alone: every residue must retain an admissible codon.
    restricted = engine.achievable_gc_bounds(protein, spec.forbidden_codons)
    assert restricted["feasible_codon_assignment"] is True
    assert restricted["residues_without_codon"] == []

    # gc_global_min alone: reachable with every codon available.
    unrestricted = engine.achievable_gc_bounds(protein)
    assert spec.gc_global_min <= unrestricted["gc_max"]

    # Together: impossible.
    assert spec.gc_global_min > restricted["gc_max"]
    assert gold["each_constraint_satisfiable_alone"] == {
        "forbidden_codons": True, "gc_global_min": True}


def test_bounds_are_attained_not_merely_valid(engine):
    """The GC bounds must be exact, since the proof rests on them.

    Codon choices are independent, so the per-residue extremes are jointly
    realisable: a sequence attaining gc_max exists. Construct it and check.
    """
    spec = ConflictSpec.load(CAMPAIGN)
    protein = read_protein(CAMPAIGN, spec.protein_input)
    forbidden = set(spec.forbidden_codons)

    best = "".join(
        max((c for c in AA_CODONS[aa] if c not in forbidden),
            key=lambda c: (gc_count(c), c))
        for aa in protein
    )
    attained = engine.gc_fraction(best)
    assert round(attained, 6) == _gold()["restricted_gc_max"]
    assert engine.translate(best) == protein
    # ...and it still falls short of the requirement, which is the whole point.
    assert attained < spec.gc_global_min


def test_margins_are_not_a_knife_edge():
    """An auditor recomputing the proof must not land on a rounding boundary."""
    gold = _gold()
    assert gold["witness"]["margin"] > 0.01
    assert gold["unrestricted_gc_max"] - gold["gc_global_min"] > 0.01


def test_verdict_and_pair_are_the_constructed_ones():
    gold = _gold()
    assert gold["verdict"] == "infeasible_gc_unreachable"
    assert gold["verdict"] in VERDICTS
    assert gold["conflicting_constraints"] == ["forbidden_codons", "gc_global_min"]


def test_only_gc_relevant_relaxations_resolve_the_conflict():
    """R4's discriminating content: the five composition constraints restrict
    which sequences are admissible without changing what GC any sequence can
    reach, so relaxing them cannot help."""
    relax = _gold()["single_constraint_relaxations"]
    assert set(relax) == set(CONFLICT_CONSTRAINTS)
    resolving = {c for c, v in relax.items() if v["gc_requirement_reachable"]}
    assert resolving == {"forbidden_codons", "gc_global_min"}
    for c in set(CONFLICT_CONSTRAINTS) - resolving:
        assert relax[c]["max_achievable_gc"] == _gold()["restricted_gc_max"]


# ---------------------------------------------------- fabrication and build guards

def test_fabricating_a_design_is_graded():
    """Supplying a design when none can exist is the behaviour this campaign
    exists to catch, so it must be a graded component and not just advice."""
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    r3 = next(r for r in spec["rungs"] if r["rung_id"] == "r3")
    comp = next(c for c in r3["components"] if c["name"] == "no_design_claimed")
    assert comp["args"]["gold"] == 0


def test_build_refuses_a_non_pairwise_or_satisfiable_conflict(tmp_path):
    """The build must not emit abstention gold it cannot prove."""
    from npbench_c.templates.construct_design.conflict import build as conflict_build
    from npbench_c.templates.construct_design.core import OracleFailure

    campaign = tmp_path / "fake"
    (campaign / "inputs").mkdir(parents=True)
    (campaign / "reference").mkdir()
    for name in ("codon_tables.json", "cloning_strategies.json",
                 "composition_limits.json", "conflict_constraints.json"):
        src = CAMPAIGN / "reference" / name
        (campaign / "reference" / name).write_text(src.read_text())
    (campaign / "inputs" / "rebh.faa").write_text(
        (CAMPAIGN / "inputs" / "rebh.faa").read_text())

    task = (CAMPAIGN / "task.yaml").read_text()
    # A requirement already reachable under the forbidden set is not infeasible.
    (campaign / "task.yaml").write_text(task.replace("gc_global_min: 0.62",
                                                     "gc_global_min: 0.40"))
    with pytest.raises(OracleFailure, match="still reachable"):
        conflict_build.analyse(campaign)

    # A requirement above even the unrestricted maximum is unsatisfiable alone,
    # so the conflict would not be pairwise.
    (campaign / "task.yaml").write_text(task.replace("gc_global_min: 0.62",
                                                     "gc_global_min: 0.95"))
    with pytest.raises(OracleFailure, match="unsatisfiable on"):
        conflict_build.analyse(campaign)


def test_required_sibling_is_recorded():
    """Abstention measured without a sufficient-evidence control rewards
    reflexive abstention: a system always answering "infeasible" would score
    full marks here. The requirement must be in the artifact, not just intent.
    """
    import yaml
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    assert "required_sibling" in task
    assert "feasible" in task["required_sibling"]


# ---------------------------------------------------------- isolation

def test_warm_r4_bundle_withholds_the_relaxation_answers(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r4", "gold/warm/r4")
    handed = json.loads((s / "assessment.json").read_text())
    assert "single_constraint_relaxations" not in handed
    assert handed["verdict"] == _gold()["verdict"]


def test_campaign_is_solvable_from_its_sandbox(tmp_path, engine):
    """Recompute the restricted bound from sandbox files alone."""
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    spec = json.loads((s / "reference" / "conflict_constraints.json").read_text())
    protein = "".join(
        l.strip() for l in (s / "inputs" / "rebh.faa").read_text().splitlines()[1:]) + "*"
    forbidden = set(spec["forbidden_codons"])
    total = sum(max(gc_count(c) for c in AA_CODONS[aa] if c not in forbidden)
                for aa in protein)
    assert round(total / (3 * len(protein)), 6) == _gold()["restricted_gc_max"]
    assert spec["gc_global_min"] == _gold()["gc_global_min"]
    assert not (s / "oracle").exists()
