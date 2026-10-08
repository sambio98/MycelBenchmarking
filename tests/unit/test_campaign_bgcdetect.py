"""BGC detection: the first S2 campaign, and three artefacts it had to refuse.

Everything hard about this campaign is upstream of the tool. Slicing the genome to
keep the run cheap breaks it twice, the two resources count coordinates
differently, and half the called regions have no curated counterpart and are not
thereby wrong. Each has a test, because each would otherwise produce a
plausible-looking answer.

The tool itself is run once per module through a cached fixture: a minimal pass
over 8.67 Mb takes about 75 seconds.
"""

from __future__ import annotations

import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.bgc_detection import analyse as bgc
from npbench_c.templates.bgc_detection.build import _assert_has_content
from npbench_c.templates.bgc_detection.core import (
    Detection,
    DetectionError,
    Interval,
    Rules,
    delete_genes,
    read_genbank,
    tool_prefix,
)
from npbench_c.templates.bgc_detection.measure import measure

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "bgcdetect-scoelicolor-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")

#: antiSMASH over a whole chromosome is slow enough that the tests run it only
#: when asked. The gold, the oracle and the sweep all exercise it; what this
#: controls is whether the unit suite re-runs it.
RUN_SLOW = __import__("os").environ.get("NPBENCH_SLOW_TOOLS") == "1"
slow = pytest.mark.skipif(
    not RUN_SLOW,
    reason="a minimal antiSMASH pass over 8.67 Mb takes ~75s; set "
           "NPBENCH_SLOW_TOOLS=1 to re-measure")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


@pytest.fixture(scope="module")
def detection():
    params = bgc.load_params(CAMPAIGN)
    return Detection.load(CAMPAIGN, params["genbank_input"],
                          params["curated_loci_input"])


# ---------------------------------------------------------------- the gates


def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_it_is_declared_s2_and_says_why():
    task = _task()
    assert task["min_tool_surface"] == "S2"
    assert "first S2 campaign" in task["tool_surface_note"]
    assert task["snapshots_required"]["image"][0] == "antismash@8.0.4"
    # The budget is deliberately larger than the benchmark's usual 25 minutes.
    assert task["budgets"]["wall_clock_minutes"] >= 45
    assert "budget_note" in task["budgets"]


# -------------------------------------------------- the slice that cannot be


def test_no_region_touches_a_record_edge():
    """The reason this campaign ships a complete 6.8 MB chromosome.

    A region reaching a record edge has its boundary set by the record rather
    than by the detection rules, and the boundary is what R2 and R3 grade. On an
    81 kb window around the actinorhodin cluster antiSMASH returned one region
    running to the record end with contig_edge=True; on the whole chromosome
    every region comes back False, and the build refuses gold otherwise.
    """
    gold = _gold()
    assert gold["contig_edge_regions"] == 0
    assert all(row["contig_edge"] is False for row in gold["region_table"])
    assert gold["record_length"] == 8667507

    with pytest.raises(DetectionError, match="contig_edge=True"):
        _assert_has_content(dict(gold, contig_edge_regions=3))


def test_the_input_is_the_whole_record_and_the_reason_is_recorded():
    provenance = json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())
    assert "unsliced" in provenance["note"]
    assert provenance["accession"] == "AL645882.2"
    rules = json.loads((CAMPAIGN / "reference" / "detection_rules.json").read_text())
    note = rules["detection"]["whole_record_note"]
    assert "translation extends out of record" in note
    assert "the slice bought nothing" in note


# ----------------------------------------------------- the two conventions


def test_the_coordinate_conversion_is_verified_not_assumed():
    """antiSMASH writes 0-based half-open; MIBiG writes 1-based inclusive.

    The second half is checked against the source annotation rather than
    inferred: BGC0000194's `to` is 5,535,091 and SCO5092 ends at 5,535,091 in
    AL645882.2. A MIBiG locus therefore spans to-from+1 bases, one more than its
    own numbers suggest at a glance.
    """
    curated = json.loads(
        (CAMPAIGN / "inputs" / "mibig_loci_AL645882.2.json").read_text())
    act = next(l for l in curated["loci"] if l["bgc"] == "BGC0000194")
    assert act["to"] == 5535091
    interval = Interval.from_mibig(act["from"], act["to"])
    assert interval.start == act["from"] - 1
    assert interval.end == act["to"]
    assert interval.length == act["to"] - act["from"] + 1
    assert interval.as_one_based() == (act["from"], act["to"])

    rules = json.loads((CAMPAIGN / "reference" / "detection_rules.json").read_text())
    assert "SCO5092" in rules["detection"]["coordinate_note"]


def test_the_antismash_convention_is_read_as_written():
    interval = Interval.from_antismash("[10992:81283](+)")
    assert (interval.start, interval.end) == (10992, 81283)
    assert interval.length == 70291
    assert interval.as_one_based() == (10993, 81283)
    with pytest.raises(DetectionError, match="unparseable"):
        Interval.from_antismash("(+)")


def test_the_region_table_reports_both_conventions():
    """Where an off-by-one hides. A submission that gets one form right and the
    other shifted has not understood the convention, and grading one would let
    that pass -- which is what the `offbyone` stub level demonstrates."""
    for row in _gold()["region_table"]:
        assert row["start_one_based"] == row["start"] + 1
        assert row["end_one_based"] == row["end"]
        assert row["length"] == row["end"] - row["start"]


# ------------------------------------------------------ what was reconciled


def test_every_curated_locus_is_detected_so_extent_is_the_question():
    """A detected/missed verdict would be a constant on this accession. What
    varies is by how much the called region exceeds the curated one."""
    gold = _gold()
    census = gold["verdict_census"]
    assert census["not_detected"] == 0
    assert census["marginal"] == 0
    assert sum(census.values()) == gold["loci_total"] == 15
    assert census["contained_much_larger"] > census["contained_generous"] > 0
    assert 0.0 < gold["jaccard_median"] < 1.0
    assert gold["jaccard_max"] - gold["jaccard_min"] > 0.5


def test_half_the_called_regions_have_no_curated_counterpart():
    """Not false positives. MIBiG records what has been characterised, which is
    not the same as what exists, and the campaign says so in its exclusions."""
    gold = _gold()
    assert gold["regions_without_a_locus"] > gold["regions"] // 3
    assert (gold["regions_matching_a_locus"] + gold["regions_without_a_locus"]
            == gold["regions"])
    excluded = [e["field"] for e in _task()["excluded_from_grading"]]
    assert any("really" in f for f in excluded)


def test_one_region_holds_two_curated_loci():
    """The coelimycin cluster contains the SCB1 butyrolactone genes, so a
    one-to-one join would silently drop one of the two."""
    gold = _gold()
    assert gold["regions_holding_several_loci"] == [21]
    shared = [row for row in gold["per_locus"] if row["region"] == 21]
    assert len(shared) == 2
    assert {row["bgc"] for row in shared} == {"BGC0000038", "BGC0000849"}
    # And the pair spans the extremes of the Jaccard range, because one is a
    # 1.7 kb locus inside a 75.6 kb region.
    assert min(row["jaccard"] for row in shared) < 0.05


def test_the_per_locus_rows_are_in_a_declared_order():
    """Ordered by MIBiG accession -- a key the grader reproduces from the shipped
    locus table alone, without needing any coordinate."""
    rows = _gold()["per_locus"]
    assert [r["bgc"] for r in rows] == sorted(r["bgc"] for r in rows)
    rules = json.loads((CAMPAIGN / "reference" / "detection_rules.json").read_text())
    assert rules["detection"]["per_locus_order"] == "MIBiG accession ascending"


# ------------------------------------------------------- the counterfactual


def test_the_deletion_variants_probe_different_dependencies():
    """Four kinds of outcome, not four counts of one.

    v1 and v2 agreeing is the measured result rather than a redundancy: deleting
    one of the two core type II PKS genes is as fatal as deleting both, so the
    rule is conjunctive. v5 changes nothing and v6 replaces a region instead of
    removing it.
    """
    summary = _gold()["deletion_summary"]
    assert summary["v0_identity"]["regions_lost"] == 0
    assert summary["v0_identity"]["regions_gained"] == 0

    # Conjunctive rule: one gene is enough to abolish the region.
    assert (summary["v1_t2pks_core_both"]["lost_products"]
            == summary["v2_t2pks_core_one"]["lost_products"] == ["T2PKS"])

    # Negative control: a biosynthetic-additional gene is not a core gene.
    assert summary["v5_additional_not_core"]["regions_lost"] == 0
    assert summary["v5_additional_not_core"]["regions_gained"] == 0

    # A shared core gene in a multi-product region: the region survives smaller.
    multi = summary["v6_multiproduct_one_product"]
    assert multi["regions_lost"] == 1 and multi["regions_gained"] == 1
    assert "butyrolactone" in multi["lost_products"]
    assert "butyrolactone" not in multi["gained_products"]

    signatures = {(o["regions_lost"], o["regions_gained"],
                   tuple(o["lost_products"]))
                  for name, o in summary.items() if name != "v0_identity"}
    assert len(signatures) >= 3


def test_the_build_refuses_gold_whose_rungs_say_nothing():
    gold = _gold()
    _assert_has_content(dict(gold))

    quiet = dict(gold, regions=3)
    with pytest.raises(DetectionError, match="anecdote"):
        _assert_has_content(quiet)

    drifted = dict(gold, deletion_outcomes={
        **gold["deletion_outcomes"],
        "v0_identity": {**gold["deletion_outcomes"]["v0_identity"],
                        "regions_lost": 2}})
    with pytest.raises(DetectionError, match="identity control"):
        _assert_has_content(drifted)

    # Every variant doing the same thing: the defect that fired on my first
    # attempt, when all four deletions lost exactly one region.
    uniform = dict(gold, deletion_outcomes={
        name: (o if name == "v0_identity"
               else {**o, "regions_lost": 1, "regions_gained": 0,
                     "lost_products": ["T2PKS"], "gained_products": []})
        for name, o in gold["deletion_outcomes"].items()})
    with pytest.raises(DetectionError, match="distinct outcomes"):
        _assert_has_content(uniform)

    flat = dict(gold, jaccard_median=1.0)
    with pytest.raises(DetectionError, match="measures nothing"):
        _assert_has_content(flat)


# ------------------------------------------------------- the gene deletion


def test_deleting_a_gene_removes_its_features_and_nothing_else(detection):
    """The sequence is untouched: antiSMASH scans the CDS translations in an
    annotated record, so removing the feature removes the protein from detection
    while every other coordinate stays where it was."""
    text = read_genbank(detection.genbank_path)
    trimmed, removed = delete_genes(text, ["SCO5087", "SCO5088"])
    assert sorted(removed) == [("CDS", "SCO5087"), ("CDS", "SCO5088"),
                               ("gene", "SCO5087"), ("gene", "SCO5088")]
    assert len(trimmed) < len(text)
    # The sequence block is byte-identical: only features were touched.
    assert text[text.index("\nORIGIN"):] == trimmed[trimmed.index("\nORIGIN"):]
    assert trimmed.count("\n     CDS ") == text.count("\n     CDS ") - 2


def test_deleting_nothing_is_a_no_op(detection):
    text = read_genbank(detection.genbank_path)
    same, removed = delete_genes(text, [])
    assert same == text and removed == []


def test_a_deletion_naming_an_absent_gene_is_refused(detection):
    """A variant that deletes nothing is not a counterfactual, and silently
    deleting nothing would read as 'the region survived'."""
    text = read_genbank(detection.genbank_path)
    with pytest.raises(DetectionError, match="no such gene"):
        delete_genes(text, ["SCO5087", "SCO99999"])


# -------------------------------------------------- the external resource


def test_the_detection_rules_are_pinned_and_not_shipped():
    """The rule set IS the gold for a region-detection campaign -- 58 rules at v5,
    88 at v7.1, 103 at 8.0.4 -- and it is AGPL package data, so the campaign
    records a fingerprint instead of the files."""
    rules = Rules.load(CAMPAIGN)
    assert len(rules.rule_paths) == 3
    assert _gold()["detection_rule_fingerprint"].startswith(rules.rule_fingerprint)
    shipped = {p.name for p in (CAMPAIGN / "inputs").iterdir()}
    assert not any(n.endswith(".txt") for n in shipped)
    declared = _task()["constraints"]["external_resources"]
    assert [e["name"] for e in declared] == ["antismash_detection_rules"]
    assert "AGPL" in declared[0]["reason"]


@image_present
def test_the_fingerprint_is_computed_from_the_image(detection):
    assert detection.rule_fingerprint().startswith(
        detection.rules.rule_fingerprint)


@slow
@image_present
def test_the_projection_is_what_gets_graded(detection):
    """Not the HTML, not the 20 MB JSON. Per region the coordinates, the sorted
    products and the contig-edge flag."""
    regions = detection.regions()
    assert len(regions) == _gold()["regions"]
    assert detection.region_table(regions) == _gold()["region_table"]
    assert all(not r.contig_edge for r in regions)


# ------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    import shutil

    stripped = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, stripped,
                    ignore=shutil.ignore_patterns("gold", "oracle",
                                                  "oracle_submission"))
    assert not (stripped / "gold").exists()
    submission = CAMPAIGN / "oracle_submission"
    assert measure(stripped, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_survives_a_missing_or_broken_report(tmp_path):
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}
    (tmp_path / "detection_report.json").write_text("{not json")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}
    (tmp_path / "detection_report.json").write_text("[]")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}


def test_an_off_by_one_in_either_convention_is_caught(tmp_path):
    """The failure the `offbyone` stub level is built from."""
    report = json.loads(
        (CAMPAIGN / "oracle_submission" / "detection_report.json").read_text())
    report["region_table"] = [{**row, "start_one_based": row["start"]}
                              for row in report["region_table"]]
    (tmp_path / "detection_report.json").write_text(json.dumps(report))
    out = measure(CAMPAIGN, tmp_path)
    assert out["region_table"] != _gold()["region_table"]
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    r2 = next(r for r in spec["rungs"] if r["rung_id"] == "r2")
    assert "region_table" in {c["name"] for c in r2["components"]}


# ------------------------------------------------------------ gold isolation


def test_the_sandbox_hands_over_no_answers(tmp_path):
    sandbox = build_sandbox(CAMPAIGN, tmp_path / "s")
    assert {p.name for p in sandbox.iterdir()} == {
        "inputs", "reference", "task.yaml", "submission"}
    blob = json.dumps({p.name: p.read_text()
                       for p in (sandbox / "reference").iterdir()})
    gold = _gold()
    assert str(gold["regions"]) not in blob
    assert str(gold["region_bases"]) not in blob


def test_the_warm_bundles_withhold_their_own_rung():
    gold = _gold()
    for rung, forbidden in (("r2", "region_table"),
                            ("r3", "verdict_census"),
                            ("r4", "deletion_summary")):
        bundle = CAMPAIGN / "gold" / "warm" / rung
        blob = json.dumps({p.name: json.loads(p.read_text())
                           for p in bundle.iterdir()})
        assert json.dumps(gold[forbidden]) not in blob
