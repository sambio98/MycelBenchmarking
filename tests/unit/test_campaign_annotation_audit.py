"""Gene-function annotation audit: the benchmark's G2 doctrine turned on MIBiG.

MIBiG records three independent evidence axes, and conflating them is the error
the campaign targets. These tests pin the axes apart and pin the build's refusal
to emit gold where the distinction would not bite.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.mibig_annotation_audit.core import (
    AnnotationAudit,
    AuditError,
    Policy,
)
from npbench_c.templates.mibig_diff.core import load_release

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "mibig-gene-function-evidence-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _gold():
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def audit():
    return AnnotationAudit.load(CAMPAIGN)


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


# ------------------------------------------------- the three evidence axes

def test_admissible_is_strictly_below_backed_for_most_categories():
    """If it were not, R3 would restate R2 and the campaign would have no
    evidence-policy content. The build refuses to emit gold in that case."""
    gold = _gold()
    gaps = gold["backed_minus_admissible"]
    assert all(v >= 0 for v in gaps.values())
    assert sum(gaps.values()) > 0
    assert sum(1 for v in gaps.values() if v > 0) >= 5


def test_the_focal_category_loses_most_of_its_backed_annotations():
    """Resistance/immunity is the focal category and the starkest case: a
    function-level evidence record is not the same as an experimentally
    supported cluster-product link."""
    gold = _gold()
    focal = gold["focal_function"]
    assert gold["backed_by_function"][focal] > gold["admissible_by_function"][focal]
    assert gold["admissible_by_function"][focal] == len(gold["focal_admissible_entries"])


def test_counting_unit_is_the_annotation_function_pair(audit):
    """A gene carrying two functions contributes two. The contract says so, and
    the census must match it."""
    gold = _gold()
    total = sum(gold["backed_by_function"].values()) + \
        sum(gold["unbacked_by_function"].values())
    assert total == gold["total_function_annotations"]


def test_locus_and_function_evidence_vocabularies_are_distinct(audit):
    """The two axes use different method names, which is why conflating them is
    a mistake a system can actually make."""
    policy = json.loads((CAMPAIGN / "reference" / "evidence_policy.json").read_text())
    function_methods = set(policy["function_evidence_methods"])
    locus_methods = set(policy["locus_evidence_methods"])
    assert "Knock-out" in function_methods
    assert "Knock-out studies" in locus_methods
    assert function_methods != locus_methods


def test_all_three_conditions_are_required(audit):
    """Dropping any one condition changes the admissible total, so none of them
    is decorative."""
    corpus = load_release(CAMPAIGN / "inputs" / "mibig_json_4.0.tar.gz")
    baseline = audit.admissible(corpus, audit.baseline)
    total = sum(baseline["admissible_by_function"].values())

    loosened_locus = Policy(
        name="any_locus", function_evidence=audit.baseline.function_evidence,
        locus_evidence=frozenset(json.loads(
            (CAMPAIGN / "reference" / "evidence_policy.json").read_text()
        )["locus_evidence_methods"]),
        require_active=True)
    loosened_status = Policy(
        name="any_status", function_evidence=audit.baseline.function_evidence,
        locus_evidence=audit.baseline.locus_evidence, require_active=False)

    assert sum(audit.admissible(corpus, loosened_locus)[
        "admissible_by_function"].values()) > total
    assert sum(audit.admissible(corpus, loosened_status)[
        "admissible_by_function"].values()) > total


# ------------------------------------------------------------- counterfactual

def test_every_policy_variant_moves_the_total():
    """A counterfactual whose variants all give the baseline answer tests
    nothing, and the build refuses to emit such gold."""
    gold = _gold()
    baseline = gold["admissible_total"]
    totals = {k: v["total"] for k, v in gold["policy_variants"].items()}
    assert len(totals) >= 4
    assert all(t != baseline for t in totals.values())
    # Tightening to knockout-only must shrink it; adding a method must grow it.
    assert totals["v1_knockout_only"] < baseline
    assert totals["v2_plus_activity_assay"] > baseline


def test_build_refuses_gold_with_no_policy_content(tmp_path):
    from npbench_c.templates.mibig_annotation_audit import build as audit_build

    campaign = tmp_path / "flat"
    (campaign / "inputs").mkdir(parents=True)
    (campaign / "reference").mkdir()
    (campaign / "reference" / "evidence_policy.json").write_text(
        (CAMPAIGN / "reference" / "evidence_policy.json").read_text())
    (campaign / "inputs" / "mibig_json_4.0.tar.gz").write_bytes(
        (CAMPAIGN / "inputs" / "mibig_json_4.0.tar.gz").read_bytes())

    task = (CAMPAIGN / "task.yaml").read_text()
    # Every variant identical to the baseline leaves R4 with no content.
    flat = task.replace('''    v1_knockout_only:
      function_evidence_accepted: ["Knock-out"]''',
                        '''    v1_knockout_only:
      function_evidence_accepted: ["Knock-out", "Heterologous expression"]''')
    flat = flat.replace('''      function_evidence_accepted: ["Knock-out", "Heterologous expression", "Activity assay"]''',
                        '''      function_evidence_accepted: ["Knock-out", "Heterologous expression"]''')
    flat = flat.replace('''      locus_evidence_accepted: ["Knock-out studies", "Heterologous expression",
                                "Enzymatic assays",
                                "Gene expression correlated with compound production",
                                "Correlation of genomic and metabolomic data",
                                "Homology-based prediction", "In vitro expression"]''',
                        '''      locus_evidence_accepted: ["Knock-out studies", "Heterologous expression"]''')
    flat = flat.replace('      require_active_status: false', '      require_active_status: true')
    (campaign / "task.yaml").write_text(flat)
    with pytest.raises(AuditError, match="no counterfactual content"):
        audit_build.analyse(campaign)


def test_census_rejects_a_function_outside_the_vocabulary(audit):
    """An unknown function name is a schema surprise, not something to bucket
    silently: the census raises so a MIBiG release that adds a term cannot
    quietly shrink the counts."""
    corpus = {"BGC9999999": {"status": "active", "genes": {"annotations": [
        {"id": "x", "functions": [{"function": {"name": "Teleportation"},
                                   "evidence": []}]}]}}}
    with pytest.raises(AuditError, match="outside the declared vocabulary"):
        audit.census(corpus)


# ---------------------------------------------------- isolation and provenance

def test_warm_r4_bundle_withholds_the_variant_answers(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r4", "gold/warm/r4")
    handed = json.loads((s / "baseline.json").read_text())
    assert "policy_variants" not in handed
    assert handed["admissible_by_function"] == _gold()["admissible_by_function"]


def test_warm_r3_bundle_hands_over_the_census_but_not_admissibility(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r3", "gold/warm/r3")
    handed = json.loads((s / "census.json").read_text())
    assert handed["backed_by_function"] == _gold()["backed_by_function"]
    assert "admissible_by_function" not in handed


def test_campaign_is_solvable_from_its_sandbox(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    audit = AnnotationAudit.load(s)
    import yaml
    params = yaml.safe_load((s / "task.yaml").read_text())["template_params"]
    corpus = load_release(s / params["corpus_archive"])
    census = audit.census(corpus)
    gold = _gold()
    assert census["backed_by_function"] == gold["backed_by_function"]
    assert audit.admissible(corpus, audit.baseline)["admissible_by_function"] == \
        gold["admissible_by_function"]
    assert not (s / "oracle").exists()


def test_catalog_deviation_is_recorded():
    """The catalogued T-L3-5 task is not what this campaign does, and the
    artifact must say so rather than the difference living only in a commit."""
    import yaml
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    assert "catalog_deviation" in task
    assert "build:image" in task["catalog_deviation"].replace(" ", "")
    excluded = {e["field"] for e in task["excluded_from_grading"]}
    assert any("target" in f for f in excluded)
