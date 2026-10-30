"""MIBiG release reconciliation: separating representation from substance.

The campaign's content is that a raw diff of two releases across a schema
migration is true and useless. These tests pin the distinction and the traps
the corpus itself supplies.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.mibig_diff.core import DiffError, ReleaseDiff, load_release

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "mibig-diff-3_1-to-4_0-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _gold():
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def diff():
    return ReleaseDiff.load(CAMPAIGN)


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


# ------------------------------------------- the representation/substance gap

def test_the_taxon_id_trap_is_in_the_gold():
    """3.1 stores the taxon id as a string, 4.0 as an integer. Raw comparison
    calls every shared entry re-annotated; normalised comparison finds 8."""
    gold = _gold()
    assert gold["raw_difference_counts"]["ncbi_tax_id"] == gold["shared_entries"]
    verdicts = gold["classified_verdict_counts"]["ncbi_tax_id"]
    assert verdicts["substantive_change"] == 8
    assert verdicts["representation_only"] == gold["shared_entries"] - 8
    assert "unchanged" not in verdicts, "every entry's raw type differs"


def test_class_vocabularies_never_agree_raw_but_mostly_agree_mapped():
    gold = _gold()
    assert gold["raw_difference_counts"]["biosynthetic_class"] == gold["shared_entries"]
    v = gold["classified_verdict_counts"]["biosynthetic_class"]
    assert v["representation_only"] > 2000
    assert v["substantive_change"] < 50
    assert v["no_counterpart_in_target"] == 73


def test_r2_and_r3_are_not_the_same_question():
    """If no field's raw count differed from its classified substantive count,
    R3 would restate R2 and the campaign would have no normalisation content.
    The build refuses to emit gold in that case; this pins that it holds."""
    gold = _gold()
    gaps = gold["raw_versus_classified_gap"]
    assert any(v > 0 for v in gaps.values())
    assert gaps["ncbi_tax_id"] == gold["shared_entries"] - 8
    assert gaps["locus_accession"] == 0, "accessions are genuinely unchanged"


def test_a_field_with_no_change_is_reported_as_such(diff):
    gold = _gold()
    assert gold["raw_difference_counts"]["locus_accession"] == 0
    assert gold["classified_verdict_counts"]["locus_accession"] == {
        "unchanged": gold["shared_entries"]}


# ---------------------------------------------------------- verdict semantics

def test_precedence_puts_a_retired_term_ahead_of_substantive_change(diff):
    """A retired vocabulary term is also a substantive change, so an unstated
    tie-break would be a grader defect rather than a hard task."""
    order = list(diff.verdict_precedence)
    assert order.index("no_counterpart_in_target") < order.index("substantive_change")
    assert order.index("absent_in_source") < order.index("substantive_change")
    assert order.index("unchanged") < order.index("representation_only")


def test_absence_is_distinguished_from_change():
    """1015 shared entries carry no compound formula in 3.1 at all, which is
    not the same as a changed formula."""
    v = _gold()["classified_verdict_counts"]["compound_formulas"]
    assert v["absent_in_source"] > 1000
    assert v["substantive_change"] < v["absent_in_source"]


# -------------------------------------------------------------- perturbations

def test_perturbations_exercise_the_normalisation_boundary():
    """p1 is the campaign in miniature: making the source taxon id an integer
    changes nothing substantive, so the verdict moves from representation_only
    to unchanged. A system that does not model the distinction gets it wrong."""
    p = _gold()["perturbations"]
    assert p["p1_taxid_type_only"]["baseline_verdict"] == "representation_only"
    assert p["p1_taxid_type_only"]["verdict"] == "unchanged"
    assert p["p2_taxid_value"]["verdict"] == "substantive_change"
    assert p["p3_formula_altered"]["verdict"] == "substantive_change"
    assert p["p4_class_retired"]["verdict"] == "no_counterpart_in_target"
    assert p["p5_organism_removed"]["verdict"] == "absent_in_source"


def test_perturbation_verdicts_span_several_classes():
    """A counterfactual whose cases all land in one class tests little."""
    verdicts = {v["verdict"] for v in _gold()["perturbations"].values()}
    assert len(verdicts) >= 4


def test_build_refuses_an_unknown_perturbation_operation():
    from npbench_c.templates.mibig_diff.build import _apply
    record = {"cluster": {"ncbi_tax_id": "1"}}
    with pytest.raises(DiffError, match="unknown perturbation"):
        _apply(record, "not_an_operation")


# ------------------------------------------------- isolation and solvability

def test_warm_r4_bundle_withholds_the_perturbation_verdicts(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r4", "gold/warm/r4")
    handed = json.loads((s / "classified.json").read_text())
    assert "perturbations" not in handed
    assert handed["classified_verdict_counts"] == _gold()["classified_verdict_counts"]


def test_warm_r3_bundle_hands_over_raw_counts_but_not_classified(tmp_path):
    """R3's work is the normalisation, so its warm start is the naive diff."""
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r3", "gold/warm/r3")
    handed = json.loads((s / "raw_diff.json").read_text())
    assert handed["raw_difference_counts"] == _gold()["raw_difference_counts"]
    assert "classified_verdict_counts" not in handed


def test_campaign_is_solvable_from_its_sandbox(tmp_path):
    """Recompute the headline counts from sandbox files alone."""
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    diff = ReleaseDiff.load(s)
    import yaml
    params = yaml.safe_load((s / "task.yaml").read_text())["template_params"]
    source = load_release(s / params["source_archive"])
    target = load_release(s / params["target_archive"])
    corpus = diff.corpus_diff(source, target)
    gold = _gold()
    assert corpus["shared_entries"] == gold["shared_entries"]
    assert corpus["classified_verdict_counts"] == gold["classified_verdict_counts"]
    assert not (s / "oracle").exists()


def test_the_leakage_ablation_and_the_competence_probe_are_distinct():
    """Conflating them makes the leakage gate ask the wrong question.

    stub-nonormalise does the work and fails one rung, so it is a tooled
    competence probe. The ablation must withhold computation itself, which for
    an S0 campaign is the only thing there is to withhold.
    """
    import yaml
    harness = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())["harness"]
    by_name = {s["name"]: s for s in harness["stub_levels"]}
    assert by_name["stub-nonormalise"].get("mode") in (None, "tooled")
    assert by_name["stub-noncompute"]["mode"] == "no_tool"
