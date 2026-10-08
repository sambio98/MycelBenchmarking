"""Annotation transfer: the closest hit is often the wrong donor, and how often
depends on three declared rules.

What this campaign can get silently wrong is not the search -- the search is the
visible work and it is easy to check. It is the bookkeeping around it: a tied top
bitscore resolved by whatever the binary emitted first, a dash in an EC number
treated as a value, and a self hit counted as support. Each has a test, because
each produces a report that looks more complete than gold rather than obviously
broken.
"""

from __future__ import annotations

import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.annotation_transfer import analyse as atr
from npbench_c.templates.annotation_transfer.build import _assert_has_content
from npbench_c.templates.annotation_transfer.core import (
    Hit,
    Rules,
    Transfer,
    TransferError,
    ec_agreement,
    ec_level,
    tool_prefix,
)
from npbench_c.templates.annotation_transfer.measure import measure

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "transfer-ec1_14-uniprot-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


@pytest.fixture(scope="module")
def transfer():
    params = atr.load_params(CAMPAIGN)
    return Transfer.load(CAMPAIGN, params["sequence_input"],
                         params["annotation_input"])


@pytest.fixture(scope="module")
def hits(transfer):
    """One all-against-all search for the whole module."""
    return transfer.hits()


# ---------------------------------------------------------------- the gates


def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_all_rungs_are_g1_and_the_surface_is_declared():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4
    task = _task()
    assert task["min_tool_surface"] == "S1"
    assert task["snapshots_required"]["image"] == ["diamond@2.2.8"]
    # Unlike the EC domain audit, nothing is read from the image but the binary.
    assert "external_resources" not in task["constraints"]


# ------------------------------------------------------- a dash is not a value


def test_a_partial_ec_number_makes_no_claim():
    assert ec_level("1.14.14.1", 4) == "1.14.14.1"
    assert ec_level("1.14.14.1", 3) == "1.14.14"
    assert ec_level("1.14.-.-", 2) == "1.14"
    assert ec_level("1.14.-.-", 3) is None
    assert ec_level("1.14.11.-", 3) == "1.14.11"
    assert ec_level("1.14.11.-", 4) is None
    # Not four fields: refuse to guess what it meant rather than pad it.
    assert ec_level("1.14.14", 3) is None
    assert ec_level("", 1) is None


def test_undecidable_is_not_the_same_as_disagreeing():
    """The distinction the campaign turns on. A query annotated to the subclass
    cannot disagree at the sub-subclass, and counting it as an error would inflate
    the headline by hundreds of queries."""
    assert ec_agreement(["1.14.14.1"], ["1.14.14.1"], 4) is True
    assert ec_agreement(["1.14.14.1"], ["1.14.14.5"], 4) is False
    assert ec_agreement(["1.14.14.1"], ["1.14.14.5"], 3) is True
    assert ec_agreement(["1.14.-.-"], ["1.14.14.1"], 3) is None
    assert ec_agreement(["1.14.11.-"], ["1.14.14.1"], 3) is False
    # Multi-EC is any-vs-any: a multifunctional enzyme's annotation is a set.
    assert ec_agreement(["1.14.11.2", "1.14.14.1"], ["1.14.14.1"], 4) is True


def test_the_set_really_contains_the_partial_annotations():
    """The rule above is only load-bearing because the data needs it."""
    gold = _gold()
    assert gold["partial_ec_mentions"] > 1000
    assert gold["entries_with_multiple_ec"] > 100
    assert gold["verdict_census"]["undecidable"] > 500


# ------------------------------------------------------------- the tie-break


def test_the_tie_break_is_declared_and_decides_real_queries():
    """602 of 4,552 queries share their top bitscore, so about one in eight has
    its closest hit chosen by the rule rather than by the score. Without a
    published rule gold would rest on the binary's emission order."""
    gold = _gold()
    assert gold["tied_top_queries"] > 300
    assert gold["tie_break"] == "bitscore descending, then subject accession ascending"
    rules = json.loads((CAMPAIGN / "reference" / "transfer_rules.json").read_text())
    assert "Load-bearing" in rules["transfer"]["tie_break_note"]


def test_the_tie_break_orders_by_score_then_accession(transfer):
    tied = (Hit("Q00002", 90.0, 500.0), Hit("Q00001", 55.0, 500.0),
            Hit("Q00003", 99.0, 900.0))
    assert [h.subject for h in transfer.ranked(tied)] == ["Q00003", "Q00001", "Q00002"]
    # Identity does not enter the ordering: the closest hit is the best-scoring
    # one, and percent identity over a short alignment is not that.
    assert transfer.ranked(tied)[1].identity == 55.0


@image_present
def test_the_declared_invocation_excludes_the_self_hit(transfer, hits):
    """Dropped while the table is parsed, not filtered later, so no downstream
    rule can transfer an annotation from a protein to itself and call it
    support."""
    assert hits
    for query, found in hits.items():
        assert all(h.subject != query for h in found)
    assert set(hits) <= set(transfer.sequences)
    assert len(hits) == _gold()["queries_with_a_hit"]


# ------------------------------------------------------------- the findings


@image_present
def test_agreement_falls_off_with_identity(transfer, hits):
    """The campaign's central claim, and the thing a broken search would destroy.

    Measured at EC level 4: 0.41 in the lowest populated band against 1.00 above
    90% identity. The build refuses to emit gold if this ordering does not hold.
    """
    curve = transfer.agreement_by_identity(hits, 4)
    populated = [b for b in curve.values() if b["decidable"] >= 10]
    assert len(populated) >= 4
    assert populated[-1]["rate"] > populated[0]["rate"]
    assert curve == _gold()["agreement_by_identity"]


def test_agreement_falls_off_with_ec_depth():
    """Level 3 agrees 0.974 and level 4 agrees 0.923, so the deeper the claim the
    more often transfer gets it wrong."""
    depth = _gold()["agreement_by_depth"]
    assert set(depth) == {"3", "4"}
    assert depth["4"]["rate"] < depth["3"]["rate"]
    assert 0.0 < depth["4"]["rate"] < 1.0
    assert depth["4"]["decidable"] < depth["3"]["decidable"]


def test_the_verdict_census_partitions_the_queries():
    """The verdicts are tested in a declared order and each query takes the first
    that applies, so the cells must sum to the queries with a hit. If they do not,
    two cells are double-counting and every rate below is suspect."""
    gold = _gold()
    assert sum(gold["verdict_census"].values()) == gold["queries_with_a_hit"]
    assert gold["verdict_census"]["no_hit"] == 0
    assert len(gold["transfer_errors_at_subsubclass"]) == \
        gold["verdict_census"]["wrong_at_subsubclass"]


def test_levels_one_and_two_are_excluded_because_they_are_constant():
    """1.0 by construction: the set is defined by its EC prefix, so a hit inside
    it shares the first two levels with every query. A component whose value is
    fixed adds weight to a conjunctive product and no information."""
    gold = _gold()
    constant = gold["constant_depth_agreement"]
    assert {k: v["rate"] for k, v in constant.items()} == {"1": 1.0, "2": 1.0}
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    assert "constant_depth_agreement" not in graded
    excluded = [e["field"] for e in _task()["excluded_from_grading"]]
    assert any("level 1" in f for f in excluded)


# ------------------------------------------------------- the counterfactual


def test_removing_the_closest_hit_moves_the_verdicts():
    """The catalogued R4, and it needs no second search: removing a hit means
    reading further down a ranking already computed."""
    gold = _gold()
    changed = gold["removal_changed"]
    assert changed["v0_identity"] == 0
    assert changed["v1_remove_top"] > 300
    assert changed["v2_remove_top_two"] > changed["v1_remove_top"]
    assert gold["removal_census"]["v0_identity"] == gold["verdict_census"]
    # The move that matters: transfers that looked supported stop being so.
    assert gold["removal_transitions"]["v1_remove_top"]["supported->wrong_at_serial"] > 50


@image_present
def test_removal_reads_the_declared_ranking(transfer, hits):
    """A removal that read an unordered hit list would drop an arbitrary hit
    rather than the closest one, and the identity control would still pass."""
    query = next(q for q in sorted(hits) if len(hits[q]) >= 3)
    ranked = transfer.ranked(hits[query])
    top = ranked[0].subject
    assert transfer.verdict(query, hits[query])[1] == top
    assert transfer.verdict(query, hits[query], frozenset({top}))[1] == ranked[1].subject


def test_the_build_refuses_gold_whose_rungs_say_nothing():
    gold = _gold()
    _assert_has_content(dict(gold))

    flat = dict(gold, verdict_census={**gold["verdict_census"],
                                      "wrong_at_serial": 0,
                                      "wrong_at_subsubclass": 0})
    with pytest.raises(TransferError, match="transfer errors"):
        _assert_has_content(flat)

    perfect = dict(gold, agreement_by_depth={
        "3": gold["agreement_by_depth"]["3"],
        "4": {"decidable": 100, "agree": 100, "rate": 1.0}})
    with pytest.raises(TransferError, match="measures nothing"):
        _assert_has_content(perfect)

    no_curve = dict(gold, agreement_by_identity={
        k: {**v, "rate": 0.9} for k, v in gold["agreement_by_identity"].items()})
    with pytest.raises(TransferError, match="rise with identity"):
        _assert_has_content(no_curve)

    drifted = dict(gold, removal_census={
        **gold["removal_census"],
        "v0_identity": {**gold["verdict_census"], "supported": 1}})
    with pytest.raises(TransferError, match="identity control"):
        _assert_has_content(drifted)

    inert = dict(gold, removal_changed={**gold["removal_changed"],
                                        "v1_remove_top": 3,
                                        "v2_remove_top_two": 4})
    with pytest.raises(TransferError, match="counterfactual shows nothing"):
        _assert_has_content(inert)

    no_ties = dict(gold, tied_top_queries=0)
    with pytest.raises(TransferError, match="tie-break is a "):
        _assert_has_content(no_ties)


# ------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    """Run against a campaign copy with no gold/ at all, and it must agree."""
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
    (tmp_path / "transfer_report.json").write_text("{not json")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}
    (tmp_path / "transfer_report.json").write_text("[]")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}


def test_the_measurer_reads_a_missing_rate_as_missing(tmp_path):
    """An empty identity band has no rate. A submission that invents one has made
    a claim about nothing, and that must not compare equal to gold's None."""
    report = json.loads(
        (CAMPAIGN / "oracle_submission" / "transfer_report.json").read_text())
    bands = report["agreement_by_identity"]
    empty = [k for k, v in bands.items() if v["decidable"] == 0]
    if not empty:
        pytest.skip("every declared identity band is populated on this set")
    bands[empty[0]] = {"decidable": 0, "agree": 0, "rate": 1.0}
    (tmp_path / "transfer_report.json").write_text(json.dumps(report))
    out = measure(CAMPAIGN, tmp_path)
    assert out["agreement_by_identity"][empty[0]]["rate"] == 1.0
    assert _gold()["agreement_by_identity"][empty[0]]["rate"] is None


# ------------------------------------------------------------ gold isolation


def test_the_sandbox_hands_over_no_answers(tmp_path):
    sandbox = build_sandbox(CAMPAIGN, tmp_path / "s")
    assert {p.name for p in sandbox.iterdir()} == {
        "inputs", "reference", "task.yaml", "submission"}
    blob = json.dumps({p.name: p.read_text()
                       for p in (sandbox / "reference").iterdir()})
    gold = _gold()
    for key in ("queries_with_a_hit", "entries_total"):
        assert str(gold[key]) not in blob


def test_the_warm_bundles_withhold_their_own_rung():
    gold = _gold()
    for rung, forbidden in (("r2", "tied_top_queries"),
                            ("r3", "verdict_census"),
                            ("r4", "removal_census")):
        bundle = CAMPAIGN / "gold" / "warm" / rung
        blob = json.dumps({p.name: json.loads(p.read_text())
                           for p in bundle.iterdir()})
        assert json.dumps(gold[forbidden]) not in blob
