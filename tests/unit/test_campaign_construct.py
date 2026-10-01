"""End-to-end behaviour of the construct-design campaign."""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR

CAMPAIGN = pathlib.Path(__file__).resolve().parents[2] / "campaigns" / "construct-ecoli-rebh-01"
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _grade_oracle():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    g = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    return grade(m, g, s)


def test_oracle_scores_exactly_one():
    """Spec 2.7: the reference solution passes 1.0 through the real grader, or
    the campaign does not ship."""
    r = _grade_oracle()
    assert r["score"] == 1.0
    assert r["depth"] == 4
    assert r["monotonic"] is True


def test_chance_floor_within_limit():
    """The cap applies to the discriminating rungs (R2 upward).

    The full floor is ~0.26 for this ladder and that is expected, not a defect:
    R1 asks for a schema-valid ORF, which a no-tool system produces by
    reverse-translating the protein. The internal sweep measured exactly that
    (no-tool clears R1 on 3/3 runs), which is why R1's chance_level is 1.0.
    Capping the full floor at 0.10 would be unsatisfiable for any ladder built
    to the documented shape.
    """
    r = _grade_oracle()
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR
    assert r["chance_floor"] > r["discriminating_chance_floor"]


def test_all_rungs_are_g1():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


def test_host_swap_flips_the_binding_constraint():
    """R4 is parametric-proof only because the answer genuinely changes with the
    host. If the binding constraint were host-invariant, an agent could carry
    the R3 answer across and clear R4 without recomputing anything."""
    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())
    assert gold["binding_constraint_flips"] == 1
    assert gold["hosts"]["ecoli"]["binding_constraint"] != \
           gold["hosts"]["streptomyces"]["binding_constraint"]


def test_graded_quantities_are_stable_across_optimiser_policies():
    """Anything graded must not depend on which legitimate optimiser ran.

    cai_full and the per-constraint costs are measured to swing far more than
    the between-constraint signal, which is why they are excluded from grading
    and recorded in task.yaml with the measurement that excluded them.
    """
    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())
    for host in ("ecoli", "streptomyces"):
        h = gold["hosts"][host]
        assert h["cai_spread_across_policies"] > 0.01, "spread should be real, not assumed"
        assert h["cai_floor"] < h["cai_min_across_policies"], \
            "the CAI floor must sit below the worst legitimate policy"


def test_unconstrained_violations_are_optimiser_independent():
    """R3's highest-information field is a pure function of protein and tables:
    the unconstrained optimum involves no search at all."""
    import sys
    sys.path.insert(0, str(CAMPAIGN / "oracle"))
    import construct as C
    from measure import read_fasta

    protein = read_fasta(CAMPAIGN / "inputs" / "rebh.faa") + "*"
    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())
    for host, strategy, key in (
        ("ecoli_bl21", "pet28a_ndei_xhoi", "ecoli"),
        ("streptomyces_coelicolor", "pset152_ndei_xhoi", "streptomyces"),
    ):
        opt = C.unconstrained_optimum(protein, host)
        counts = C.violation_counts(opt, host, strategy)
        assert counts == gold["hosts"][key]["unconstrained_violations"]
