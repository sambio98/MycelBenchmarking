"""GCF cutoff sensitivity: the partition is a function of a number you choose.

Two things here produce a wrong answer without producing an error. The clustering
table lists only clusters that joined a family, so reading it as the partition
loses a third of the corpus and makes the structure look cleaner. And a family
label is a name, so comparing labels calls two identical partitions different.
Both have tests, and both have a stub level.

The tool is slow enough that the suite runs it only on request; the gold, the
oracle and the sweep all exercise it.
"""

from __future__ import annotations

import itertools
import json
import os
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.gcf_cutoff import analyse as gcf
from npbench_c.templates.gcf_cutoff.build import _assert_has_content
from npbench_c.templates.gcf_cutoff.core import (
    Clustering,
    CutoffError,
    Rules,
    tool_prefix,
)
from npbench_c.templates.gcf_cutoff.measure import measure

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "gcfcutoff-mibig-32-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")

slow = pytest.mark.skipif(
    os.environ.get("NPBENCH_SLOW_TOOLS") != "1",
    reason="a BiG-SCAPE run over the grid takes about a minute; set "
           "NPBENCH_SLOW_TOOLS=1 to re-measure")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


@pytest.fixture(scope="module")
def clustering():
    params = gcf.load_params(CAMPAIGN)
    return Clustering.load(CAMPAIGN, params["cluster_dir"],
                           params["designed_families_input"])


# ---------------------------------------------------------------- the gates


def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_it_is_declared_s2_with_one_tool_run():
    task = _task()
    assert task["min_tool_surface"] == "S2"
    assert task["snapshots_required"]["image"][0] == "bigscape@2.0.3"
    assert "one invocation" in task["tool_surface_note"].lower()


# ------------------------------------------------- the reconstruction trap


def test_the_partition_covers_the_corpus_and_the_table_does_not(clustering):
    """The campaign's first trap, measured on the shipped corpus.

    The clustering table lists only clusters that joined a family. At the baseline
    cutoff a third of the corpus is absent from it, so the partition has to be
    reconstructed: listed records form their groups and every omitted record is a
    singleton.
    """
    gold = _gold()
    partition = gold["baseline_partition"]
    flat = sorted(m for group in partition for m in group)
    assert flat == sorted(clustering.corpus)
    assert len(flat) == gold["corpus_size"] == len(clustering.corpus)
    # The listed records are strictly fewer than the corpus, which is what makes
    # the reconstruction load-bearing rather than a formality.
    listed = sum(len(g) for g in partition if len(g) > 1)
    assert listed == gold["clustered_records"] < gold["corpus_size"]
    assert gold["singletons"] > 0


def test_the_build_refuses_a_corpus_that_never_exercises_it():
    gold = _gold()
    _assert_has_content(dict(gold))
    with pytest.raises(CutoffError, match="singleton reconstruction"):
        _assert_has_content(dict(gold, singletons=0))
    with pytest.raises(CutoffError, match="all singletons"):
        _assert_has_content(dict(gold, groups=gold["corpus_size"]))


def test_the_reconstruction_rejects_a_record_outside_the_corpus(clustering):
    """A clustering table naming something the corpus does not hold means the
    join key is wrong, and every count below it would be silently off."""
    with pytest.raises(CutoffError, match="does not cover the corpus"):
        clustering.partition({"NOT_IN_CORPUS": ("other", "FAM_00001")})


# --------------------------------------------- a label is not the answer


def test_the_graded_artifact_is_the_partition_not_the_labels():
    """Two runs can agree on which clusters belong together and disagree on what
    the groups are called. The partition is graded; the labels are excluded."""
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    assert "baseline_partition" in graded
    blob = json.dumps(_gold()["baseline_partition"])
    assert "FAM_" not in blob
    excluded = [e["field"] for e in _task()["excluded_from_grading"]]
    assert any("label" in f for f in excluded)


def test_the_partition_is_sorted_canonically():
    """Groups as sorted member lists, sorted among themselves, so a reordered
    feature list is not mistaken for a different partition."""
    partition = _gold()["baseline_partition"]
    assert partition == sorted(sorted(g) for g in partition)


def test_agreement_is_counted_over_pairs(clustering):
    """A pair-counting measure rather than a label comparison. Recomputed here
    from the shipped families and the graded partition."""
    gold = _gold()
    designed = clustering.designed_pairs()
    called = set()
    for group in gold["baseline_partition"]:
        called |= {tuple(sorted(p)) for p in itertools.combinations(sorted(group), 2)}
    assert gold["designed_pairs"] == len(designed)
    assert gold["recovered_pairs"] == len(called & designed)
    assert gold["false_join_pairs"] == len(called - designed)
    assert gold["missed_pairs"] == len(designed - called)
    union = gold["recovered_pairs"] + gold["false_join_pairs"] + gold["missed_pairs"]
    assert gold["pair_jaccard"] == round(gold["recovered_pairs"] / union, 6)


# ----------------------------------------------- the families are independent


def test_the_families_come_from_chemistry_not_from_clustering():
    """What makes the reconciliation a comparison rather than a restatement: the
    groups are entries sharing an InChIKey connectivity block, selected by a rule
    recorded with the input."""
    payload = json.loads(
        (CAMPAIGN / "inputs" / "designed_families.json").read_text())
    assert "InChIKey connectivity block" in payload["note"]
    assert "connectivity" in payload["selection_rule"]
    assert payload["excluded_placeholder_names"]
    assert "not necessarily homologous" in payload["excluded_placeholder_reason"]


def test_the_families_and_controls_partition_the_corpus(clustering):
    members = sorted({m for v in clustering.designed.values() for m in v}
                     | set(clustering.controls))
    assert members == sorted(clustering.corpus)
    for members_of in clustering.designed.values():
        assert not set(members_of) & set(clustering.controls)


def test_the_tools_default_cutoff_recovers_nothing_exactly():
    """The campaign's headline. BiG-SCAPE's documented default is the baseline, and
    at it no chemistry-derived family comes back as exactly one group -- every one
    is split, which is why the verdict vocabulary is about the shape of the
    disagreement."""
    gold = _gold()
    assert gold["baseline_cutoff"] == "0.3"
    assert gold["families_exact"] == 0
    assert gold["verdict_census"]["exact"] == 0
    assert gold["verdict_census"]["fragmented"] == len(gold["designed_families"])
    assert all(block["groups_spanned"] > 1 for block in gold["per_family"].values())
    assert not gold["controls_absorbed"]


# ------------------------------------------------------------- the sweep


def test_agreement_is_not_monotone_in_the_cutoff():
    """Too low fragments the families, too high absorbs a negative control, and
    the best cutoff is strictly inside the grid. The build refuses gold if either
    end is invisible, because then the optimum might lie outside it."""
    gold = _gold()
    sweep = gold["cutoff_sweep_summary"]
    ordered = sorted(sweep, key=float)
    jaccards = [sweep[c]["pair_jaccard"] for c in ordered]

    assert gold["best_cutoff"] not in (ordered[0], ordered[-1])
    assert max(jaccards) == gold["best_pair_jaccard"] < 1.0
    # Rising below the optimum, and strictly lower above it.
    best = ordered.index(gold["best_cutoff"])
    assert jaccards[:best] == sorted(jaccards[:best])
    assert max(jaccards[best + 1:]) < jaccards[best]

    assert gold["first_cutoff_with_a_false_join"] is not None
    assert all(sweep[c]["false_join_pairs"] == 0 for c in ordered
               if float(c) < float(gold["first_cutoff_with_a_false_join"]))


def test_the_upper_bound_is_a_control_being_absorbed():
    """What stops the sweep from just preferring the largest cutoff."""
    gold = _gold()
    first = gold["first_cutoff_absorbing_a_control"]
    assert first is not None
    assert float(first) > float(gold["best_cutoff"])
    assert gold["cutoff_sweep"][first]["controls_absorbed"]
    assert gold["first_cutoff_with_a_false_join"] == first


def test_a_family_no_cutoff_recovers_reports_null():
    """A real answer, not a gap: the tool does not reproduce that chemistry group
    at any declared cutoff."""
    gold = _gold()
    recovered = gold["first_cutoff_recovering_each_family"]
    assert set(recovered) == set(gold["designed_families"])
    assert any(v is None for v in recovered.values())
    assert any(v is not None for v in recovered.values())
    for name, cutoff in recovered.items():
        if cutoff is None:
            assert all(row["family_verdicts"][name] != "exact"
                       for row in gold["cutoff_sweep"].values())
        else:
            assert gold["cutoff_sweep"][cutoff]["family_verdicts"][name] == "exact"


def test_interval_score_is_not_used_and_the_reason_is_recorded():
    """The catalogued R4 named it. It penalises vagueness against a max_width
    taken from a measured oracle spread, and a deterministic tool on a declared
    grid has no spread to take one from -- so any max_width would be chosen rather
    than measured, which is what the primitive exists to avoid."""
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    scorers = {c["scorer"] for r in spec["rungs"] for c in r["components"]}
    assert "interval_score" not in scorers
    assert "interval_score" in _task()["catalog_deviation"]
    assert "max_width" in _task()["catalog_deviation"]


def test_the_build_refuses_a_sweep_that_shows_nothing():
    gold = _gold()
    sweep = gold["cutoff_sweep"]
    flat = dict(gold, cutoff_sweep={
        c: {**row, "pair_jaccard": 0.5} for c, row in sweep.items()})
    with pytest.raises(CutoffError, match="does not depend on the cutoff"):
        _assert_has_content(flat)

    no_overmerge = dict(gold, first_cutoff_with_a_false_join=None)
    with pytest.raises(CutoffError, match="never reaches over-merging"):
        _assert_has_content(no_overmerge)

    at_edge = dict(gold, best_cutoff=max(sweep, key=float))
    with pytest.raises(CutoffError, match="highest declared one"):
        _assert_has_content(at_edge)

    perfect = dict(gold, best_pair_jaccard=1.0)
    with pytest.raises(CutoffError, match="no disagreement to grade"):
        _assert_has_content(perfect)

    none_recovered = dict(gold, first_cutoff_recovering_each_family={
        k: None for k in gold["first_cutoff_recovering_each_family"]})
    with pytest.raises(CutoffError, match="empty everywhere"):
        _assert_has_content(none_recovered)


# ----------------------------------------------------- the declared grid


def test_the_baseline_must_be_in_the_grid():
    rules = Rules.load(CAMPAIGN)
    assert rules.baseline_cutoff in rules.cutoffs
    assert len(rules.cutoffs) >= 4
    assert _gold()["cutoffs"] == [rules.label(c) for c in rules.cutoffs]


def test_every_cutoff_comes_from_one_run():
    """Declared, because it is what makes two cutoffs comparable: the distances
    behind them are identical by construction."""
    rules = json.loads((CAMPAIGN / "reference" / "cutoff_rules.json").read_text())
    assert "One BiG-SCAPE run covers the whole grid" in \
        rules["clustering"]["cutoffs_note"]
    assert rules["clustering"]["include_gbk"] == "BGC"
    assert "no valid input GBKs" in rules["clustering"]["include_gbk_note"]


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
    (tmp_path / "clustering_report.json").write_text("{not json")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}
    (tmp_path / "clustering_report.json").write_text("[]")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}


def test_the_measurer_canonicalises_group_order(tmp_path):
    """A submission listing the same partition in another order must compare
    equal: the order is a rendering, not a result."""
    report = json.loads(
        (CAMPAIGN / "oracle_submission" / "clustering_report.json").read_text())
    report["baseline_partition"] = [list(reversed(g)) for g in
                                    reversed(report["baseline_partition"])]
    (tmp_path / "clustering_report.json").write_text(json.dumps(report))
    out = measure(CAMPAIGN, tmp_path)
    assert out["baseline_partition"] == _gold()["baseline_partition"]


def test_the_measurer_keeps_a_null_recovery_cutoff_null(tmp_path):
    """Normalising null to a string would turn 'no cutoff recovers this family'
    into a wrong answer about which one does."""
    report = json.loads(
        (CAMPAIGN / "oracle_submission" / "clustering_report.json").read_text())
    (tmp_path / "clustering_report.json").write_text(json.dumps(report))
    out = measure(CAMPAIGN, tmp_path)
    recovered = out["sweep_boundaries"]["first_cutoff_recovering_each_family"]
    assert any(v is None for v in recovered.values())
    assert recovered == _gold()["sweep_boundaries"][
        "first_cutoff_recovering_each_family"]


# ------------------------------------------------------------ gold isolation


def test_the_sandbox_hands_over_no_answers(tmp_path):
    sandbox = build_sandbox(CAMPAIGN, tmp_path / "s")
    assert {p.name for p in sandbox.iterdir()} == {
        "inputs", "reference", "task.yaml", "submission"}
    blob = json.dumps({p.name: p.read_text()
                       for p in (sandbox / "reference").iterdir()})
    gold = _gold()
    assert json.dumps(gold["baseline_partition"]) not in blob
    assert str(gold["pair_jaccard"]) not in blob


def test_the_warm_bundles_withhold_their_own_rung():
    gold = _gold()
    for rung, forbidden in (("r2", "baseline_partition"),
                            ("r3", "verdict_census"),
                            ("r4", "cutoff_sweep_summary")):
        bundle = CAMPAIGN / "gold" / "warm" / rung
        blob = json.dumps({p.name: json.loads(p.read_text())
                           for p in bundle.iterdir()})
        assert json.dumps(gold[forbidden]) not in blob


@slow
@image_present
def test_the_tool_run_reproduces_the_graded_partition(clustering):
    assignments = clustering.assignments()
    label = clustering.rules.label(clustering.rules.baseline_cutoff)
    assert clustering.partition(assignments[label]) == _gold()["baseline_partition"]
