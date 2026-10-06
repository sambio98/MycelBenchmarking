"""ChEMBL selectivity: measured-or-computed only, and the licence travels with it.

Two things make this campaign different from the rest. Its gold is derived from a
ShareAlike source, so it is the benchmark's first `share_alike` campaign and the
obligation has to travel with the files. And its admission rules are the whole
task: every one reads a structured ChEMBL field, and the rejection reasons have to
partition the rejected set exactly or the breakdown means nothing.
"""

from __future__ import annotations

import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.readiness.gate import LICENSE_CLASSES
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.chembl_selectivity import analyse as sel
from npbench_c.templates.chembl_selectivity.build import _assert_has_content
from npbench_c.templates.chembl_selectivity.core import (
    Admission,
    Selectivity,
    SelectivityError,
    median,
    minimum,
)

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "selectivity-saureus-topoisomerase-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def selectivity():
    params = sel.load_params(CAMPAIGN)
    return Selectivity.load(CAMPAIGN, params["activities_input"],
                            params["focal_target"], params["reference_target"])


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


# ------------------------------------------------------------- the licence

def test_this_is_the_share_alike_campaign_and_carries_its_notice():
    """The ChEMBL decision, made concrete. The obligation attaches to these data
    files and travels with them, which is why the notice is in the artifact rather
    than in a design document."""
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    assert task["license_class"] == "share_alike"
    assert task["license_class"] in LICENSE_CLASSES
    notice = task["license_notice"]
    assert {"source", "license", "attribution", "redistribution"} <= set(notice)
    assert "CC BY-SA 3.0" in notice["license"]
    assert "ChEMBL" in notice["attribution"]
    # The reason the obligation is local has to be stated where a redistributor
    # reads it, not only where a reviewer does.
    assert "collection" in notice["redistribution"].lower()


# ----------------------------------------------------------- the admission

def test_the_rejection_reasons_partition_the_rejected_set(selectivity):
    """The reasons are tested in a declared order and each activity takes the
    first that applies. If they did not partition, two rules would be
    double-counting and every number downstream would be suspect."""
    gold = _gold()
    for target, rejected in gold["rejected_by_target"].items():
        of_type = sum(1 for a in selectivity.of_type()
                      if a["target_chembl_id"] == target)
        assert gold["admitted_by_target"][target] + sum(rejected.values()) == of_type


def test_every_rejection_rule_that_matters_here_has_instances():
    """A rule with no instances is untested by this campaign. Three have many;
    the confidence rule has none, which is recorded rather than hidden."""
    gold = _gold()
    totals: dict[str, int] = {}
    for rejected in gold["rejected_by_target"].values():
        for reason, n in rejected.items():
            totals[reason] = totals.get(reason, 0) + n
    for reason in ("censored_or_missing_relation", "unconvertible_units",
                   "flagged_by_chembl", "potential_duplicate"):
        assert totals[reason] > 0, reason
    assert totals["low_confidence_assay"] == 0


def test_the_confidence_rule_needs_a_join_and_is_a_no_op_here(selectivity):
    """confidence_score lives on the ASSAY, not the activity -- a system looking
    for it on the activity finds nothing and admits everything. On this pair every
    assay scores 7, so the rule rejects nothing, and the campaign says so."""
    assert all("confidence_score" not in a for a in selectivity.activities)
    pair = (selectivity.focal_target, selectivity.reference_target)
    # The no-op claim is about the assays the campaign actually reads: the ones
    # behind the declared-type activities on the pair. The input ships the assay
    # table for every activity type on both targets, so it is wider than that,
    # and a blanket claim over it would be a different (and false) statement.
    on_pair = {a["assay_chembl_id"] for a in selectivity.of_type()
               if a["target_chembl_id"] in pair and a.get("assay_chembl_id")}
    assert len(on_pair) < len(selectivity.assays)
    assert all(selectivity.assays[i]["confidence_score"] == 7 for i in on_pair)
    assert _gold()["assay_confidence_scores"] == {"7": len(on_pair)}
    assert _gold()["assay_confidence_scores"] == {"7": 105}
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    assert "confidence" in task["catalog_deviation"]


def test_a_censored_relation_is_not_a_measurement(selectivity):
    """A '>' is a bound. Admitting one as a value is how a potency table acquires
    numbers nobody measured, and the v1 variant shows the size of that error."""
    censored = [a for a in selectivity.of_type()
                if a.get("standard_relation") in (">", "<")]
    assert censored
    for activity in censored[:20]:
        assert selectivity.reject_reason(activity, Admission()) == \
            "censored_or_missing_relation"
    gold = _gold()
    assert (gold["variant_eligible"]["v1_admit_censored"]
            > gold["variant_eligible"]["v0_baseline"])


# --------------------------------------------------------- the aggregation

def test_the_declared_aggregators():
    assert median([3.0, 1.0, 2.0]) == 2.0
    # Even count: the mean of the two middle values, declared and rounded.
    assert median([1.0, 2.0, 3.0, 4.0]) == 2.5
    assert minimum([3.0, 1.0, 2.0]) == 1.0
    with pytest.raises(SelectivityError):
        median([])


def test_the_aggregator_choice_changes_the_answer():
    """Which is why it is declared rather than assumed: the minimum keeps the most
    potent replicate and shifts the verdict for sixteen molecules."""
    gold = _gold()
    assert (gold["variant_focal_counts"]["v4_minimum_not_median"]
            != gold["variant_focal_counts"]["v0_baseline"])


# --------------------------------------------------------- the selectivity

def test_both_directions_are_populated():
    """If one side were empty the verdict would be a foregone conclusion, and the
    build refuses in that case."""
    gold = _gold()
    assert gold["focal_selective_counts"]["100"] > 0
    assert gold["reference_selective_counts"]["2"] > 0


def test_the_committed_set_matches_its_count():
    gold = _gold()
    highest = str(max(int(k) for k in gold["focal_selective_counts"]))
    assert (len(gold["focal_selective_at_highest_threshold"])
            == gold["focal_selective_counts"][highest])


def test_the_ratio_direction_is_declared_and_inverting_it_is_catchable(selectivity):
    """Inverting the direction gets every magnitude right and every verdict
    backwards, which is what the stub-inverted probe does and what R3 catches."""
    ratios = selectivity.ratios(Admission())
    flipped = {m: 1.0 / r for m, r in ratios.items()}
    original = selectivity.verdicts(ratios)
    inverted = selectivity.verdicts(flipped)
    assert original["focal_selective_counts"] != inverted["focal_selective_counts"]


def test_the_extreme_molecule_has_a_declared_tie_break():
    gold = _gold()
    assert gold["most_selective_molecule"] in \
        gold["focal_selective_at_highest_threshold"]
    assert gold["most_selective_ratio"] > 100


# -------------------------------------------------------------- build refusal

@pytest.mark.parametrize("mutate,message", [
    (lambda g: g["focal_selective_counts"].__setitem__("100", 0),
     "highest\nthreshold|highest threshold"),
    (lambda g: g["reference_selective_counts"].__setitem__("2", 0),
     "foregone conclusion"),
    (lambda g: g["admission_variants"]["v0_baseline"].__setitem__(
        "focal_selective_counts", {"2": 1, "10": 1, "100": 1}), "baseline variant"),
])
def test_build_refuses_gold_without_content(mutate, message, selectivity):
    doctored = json.loads(json.dumps(_gold()))
    doctored["_of_type_by_target"] = {
        t: sum(1 for a in selectivity.of_type() if a["target_chembl_id"] == t)
        for t in sorted(doctored["admitted_by_target"])}
    mutate(doctored)
    with pytest.raises(SelectivityError, match=message):
        _assert_has_content(doctored)


# --------------------------------------------------- sandbox self-containment

def test_sandbox_carries_no_gold(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "cold")
    assert sorted(p.name for p in s.iterdir()) == ["inputs", "reference",
                                                   "submission", "task.yaml"]


def test_campaign_is_solvable_from_its_sandbox(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    recomputed = sel.analyse(s)
    gold = _gold()
    for key in ("admitted_by_target", "rejected_by_target", "eligible_molecules",
                "focal_selective_counts", "focal_selective_at_highest_threshold",
                "admission_variants"):
        assert recomputed[key] == gold[key], key
