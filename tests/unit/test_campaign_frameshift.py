"""A planted frameshift: what a caller reports when a gene is demonstrably broken.

"Which locus is a pseudogene" is a judgement with no oracle, so this campaign
breaks the gene itself and grades what the caller says about a gene that is known
to be broken by construction.

Three things here produce a wrong answer without producing an error.

The coordinate lift past the deletion is one addition per call. Skip it and every
call downstream of the deleted base moves by exactly one, so the baseline matches
almost nothing: the target reads as abolished whatever actually happened, and the
collateral columns run into the thousands instead of zero. That is the competence
probe.

The mode is a flag. Single mode trains on its input, so one deleted base out of
8.7 megabases moves calls across the whole chromosome and the planted
frameshift's own signature -- one call lost, two gained -- is not separable
inside it. The companion gene-calling campaign pins the opposite mode on this
same chromosome, each for a stated measurement, and the noise arm here is what
makes that a measurement rather than a preference.

And the length ladder is the campaign's own measurement, not a convention. The
caller's recovery of a broken gene turns over somewhere between 170 and 280 bases
on this chromosome, and inside that band it recovers some genes and not others. A
ladder of round numbers steps over the band and makes the per-length verdict list
one threshold restated once per target; the build refuses gold in that case.

The caller is fast -- twenty seconds over the chromosome -- so the suite runs it
where it needs it, and leans on gold where the answer is already recorded.
"""

from __future__ import annotations

import copy
import dataclasses
import json
import pathlib
import shutil

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.templates.frameshift_response import analyse as fsa
from npbench_c.templates.frameshift_response.build import (
    COUNTERFACTUAL_KEYS,
    INVENTORY_KEYS,
    RESPONSE_KEYS,
    VERDICT_KEYS,
    _assert_has_content,
)
from npbench_c.templates.frameshift_response.core import (
    Call,
    Frameshift,
    FrameshiftError,
    Rules,
)
from npbench_c.templates.frameshift_response.measure import measure
from npbench_c.templates.frameshift_response.solve import render
from npbench_c.tools.resources import tool_prefix

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "frameshift-scoelicolor-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


def _rules() -> Rules:
    return Rules.load(CAMPAIGN)


def _spec() -> dict:
    return json.loads((CAMPAIGN / "grading.json").read_text())


@pytest.fixture(scope="module")
def shift() -> Frameshift:
    params = fsa.load_params(CAMPAIGN)
    return Frameshift.load(CAMPAIGN, params["contig_input"])


# ------------------------------------------------------------ what is shipped


def test_the_mutated_contigs_are_not_shipped():
    """They follow from the declared rule. Shipping them would hand over both the
    target positions and the coordinate lift, which is the whole experiment."""
    shipped = sorted(p.name for p in (CAMPAIGN / "inputs").iterdir())
    assert shipped == ["contig.fna.gz", "provenance.json"]
    assert "not shipped" in _task()["license_notice"]["redistribution"]


def test_the_provenance_does_not_carry_the_graded_values():
    """`contig_length` and `contig_sha256` are graded at R1."""
    provenance = json.loads(
        (CAMPAIGN / "inputs" / "provenance.json").read_text())
    gold = _gold()
    flat = json.dumps(provenance)
    assert str(gold["contig_length"]) not in flat
    assert gold["contig_sha256"] not in flat


def test_the_input_is_the_same_contig_as_the_gene_calling_campaign():
    """Two campaigns on one chromosome with opposite modes, so the comparison
    between them is about the mode and nothing else."""
    other = CAMPAIGN.parent / "genecall-scoelicolor-01"
    if not (other / "gold" / "gold.json").is_file():
        pytest.skip("gene-calling campaign gold not generated")
    theirs = json.loads((other / "gold" / "gold.json").read_text())
    gold = _gold()
    assert gold["contig_sha256"] == theirs["contig_sha256"]
    assert gold["baseline_mode"] != theirs["baseline_mode"]


# -------------------------------------------------------------- the rules file


def test_the_baseline_and_the_noise_mode_are_different_modes(tmp_path):
    """The noise arm exists to compare two modes. Pointing both at one mode
    compares a run with itself and reports a ratio of zero."""
    copy = tmp_path / "campaign"
    (copy / "reference").mkdir(parents=True)
    spec = json.loads((CAMPAIGN / "reference" / "frameshift_rules.json").read_text())
    spec["frameshift"]["noise_mode"] = spec["frameshift"]["baseline_mode"]
    (copy / "reference" / "frameshift_rules.json").write_text(json.dumps(spec))
    with pytest.raises(FrameshiftError, match="the same mode"):
        Rules.load(copy)


def test_an_undeclared_mode_is_refused(tmp_path):
    copy = tmp_path / "campaign"
    (copy / "reference").mkdir(parents=True)
    spec = json.loads((CAMPAIGN / "reference" / "frameshift_rules.json").read_text())
    spec["frameshift"]["baseline_mode"] = "normal"
    (copy / "reference" / "frameshift_rules.json").write_text(json.dumps(spec))
    with pytest.raises(FrameshiftError, match="is not a declared mode"):
        Rules.load(copy)


def test_the_target_ladder_spans_the_measured_turnover_band():
    """The dense part of the ladder is a measurement: every ten bases from 170 to
    280, which is where this caller stops recovering broken genes reliably on
    this chromosome. Declared at a uniform step, so no length inside the band is
    chosen over another."""
    lengths = list(_rules().target_lengths)
    band = [n for n in lengths if 170 <= n <= 280]
    assert band == list(range(170, 290, 10))
    assert min(lengths) <= 100 and max(lengths) >= 20000
    rules = json.loads(
        (CAMPAIGN / "reference" / "frameshift_rules.json").read_text())["frameshift"]
    assert "turnover_band_provenance" in rules


# --------------------------------------------------- the planted frameshift


def test_every_deleted_base_is_inside_the_gene_it_breaks():
    for label, row in _gold()["target_table"].items():
        assert row["start"] <= row["deleted_at"] <= row["end"], label


def test_the_indel_position_is_a_fraction_and_never_lands_outside(shift):
    """The rule reads the same for a 100 bp gene and a 22 kb one, which means the
    extremes have to be clamped: at 0.0 and 1.0 it must still be inside."""
    target = Call(start=1000, end=1099, strand="+", partial="00")
    for fraction in (0.0, 0.05, 0.5, 0.95, 1.0):
        position = shift.indel_position(target, fraction)
        assert target.start <= position <= target.end, fraction
    assert shift.indel_position(target, 0.0) == 1000
    assert shift.indel_position(target, 0.5) == 1000 + 49


def test_no_declared_length_audits_a_gene_another_length_already_took():
    gold = _gold()
    spans = [(r["start"], r["end"], r["strand"])
             for r in gold["target_table"].values()]
    assert len(set(spans)) == len(spans)


def test_a_target_that_collides_is_refused(shift):
    """With one call on the contig, every declared length wants it. Auditing one
    gene twice makes a length stand for nothing, so the rule refuses."""
    only = (Call(start=1, end=300, strand="+", partial="00"),)
    tiny = dataclasses.replace(shift.rules, target_lengths=(300, 301))
    with pytest.raises(FrameshiftError, match="no call is left to stand for"):
        dataclasses.replace(shift, rules=tiny).targets(only)


@pytest.mark.parametrize("want", [100, 100000])
def test_a_declared_length_nothing_on_the_contig_is_near_is_refused(shift, want):
    """Both directions, because the obvious one-sided test cannot catch the
    second case: the gap between a 300 bp call and a declared 100 kb is bounded
    by 100 kb, so a difference test passes and the campaign would audit a 300 bp
    gene as its 100 kb target."""
    only = (Call(start=1, end=300, strand="+", partial="00"),)
    far = dataclasses.replace(shift.rules, target_lengths=(want,))
    with pytest.raises(FrameshiftError, match="more than a factor of two away"):
        dataclasses.replace(shift, rules=far).targets(only)


def test_the_declared_ladder_itself_clears_the_proximity_rule(shift):
    """The rule has to be satisfiable on the shipped contig, and it is the
    outermost rungs of the ladder that test it: 100 bp against a 99 bp call and
    20000 against a 22392 bp one."""
    gold = _gold()
    for label, row in gold["target_table"].items():
        want = int(label)
        assert want / 2 <= row["length"] <= 2 * want, (label, row["length"])


# -------------------------------------------------------------- the response


def test_no_target_survives_its_own_frameshift():
    """A planted single-base deletion cannot leave the caller reporting the gene's
    exact coordinates. A target reported unchanged means the deletion did not
    land where the campaign thinks it did, or the lift is wrong -- and both
    failures produce plausible-looking tables."""
    gold = _gold()
    assert gold["verdict_census"]["unchanged"] == 0
    assert all(row["verdict"] != "unchanged"
               for row in gold["response_table"].values())


def test_the_verdict_is_classified_by_overlap_and_not_by_a_window(shift):
    """A frameshifted gene can produce an ORF running PAST where the original
    ended. Any fixed window around the target misclassifies it, and the ladder
    contains such a case: a new call longer than the gene it came from."""
    gold = _gold()
    overshoot = [(label, row["target_length"], row["new_call_lengths"])
                 for label, row in gold["response_table"].items()
                 if row["new_call_lengths"]
                 and max(row["new_call_lengths"]) > row["target_length"]]
    assert overshoot, gold["response_table"]

    # And the classification itself: one new call reaching far past the target's
    # end is still an overlap, so it is truncated_to_one and not abolished.
    target = Call(start=1000, end=1300, strand="+", partial="00")
    mutated = (Call(start=1100, end=2000, strand="+", partial="00"),)
    block = shift.response((target,), target, mutated)
    assert block["verdict"] == "truncated_to_one"
    assert block["bases_of_the_target_still_called"] == 201


def test_a_new_call_that_misses_the_target_is_not_the_target_recovered(shift):
    target = Call(start=1000, end=1300, strand="+", partial="00")
    elsewhere = (Call(start=5000, end=5300, strand="+", partial="00"),)
    block = shift.response((target,), target, elsewhere)
    assert block["verdict"] == "abolished"
    assert block["other_calls_gained"] == 1
    assert block["other_calls_lost"] == 0


def test_the_collateral_does_not_count_the_target_itself(shift):
    """`other_calls_lost` is what a single deleted base cost ELSEWHERE. Counting
    the broken gene in it would make every row report at least one."""
    target = Call(start=1000, end=1300, strand="+", partial="00")
    neighbour = Call(start=2000, end=2300, strand="+", partial="00")
    block = shift.response((target, neighbour), target, (neighbour,))
    assert block["verdict"] == "abolished"
    assert block["other_calls_lost"] == 0


def test_the_bases_still_called_never_exceed_the_overlap(shift):
    """The coverage column counts the target's OWN bases that a new call covers,
    so it is bounded by the gene's length however long the new ORF is."""
    gold = _gold()
    for label, row in gold["response_table"].items():
        assert row["bases_of_the_target_still_called"] >= 0, label
        if row["overlapping_new_calls"] <= 1:
            assert (row["bases_of_the_target_still_called"]
                    <= row["target_length"]), label


# ---------------------------------------------------------------- the verdicts


def test_the_census_counts_every_declared_verdict_always():
    """The vocabulary is closed, so a verdict that does not occur is reported as
    zero rather than missing. A census whose keys depend on the data would make
    the measurer's shape depend on the answer."""
    gold, rules = _gold(), _rules()
    assert set(gold["verdict_census"]) == set(rules.verdicts)
    assert sum(gold["verdict_census"].values()) == len(gold["response_table"])


def test_three_verdicts_occur_on_the_ladder():
    """Below the turnover nothing is recovered, above it everything splits, and
    inside it one half is sometimes callable. If fewer than three terms occur the
    ladder does not cross the band."""
    census = _gold()["verdict_census"]
    assert len([v for v, n in census.items() if n]) >= 3
    assert census["split_into_two"] and census["abolished"]


def test_the_verdict_is_not_a_threshold_in_length():
    """The campaign's finding, and the reason R3 grades a per-length list rather
    than a cut point. Whether each half of a broken gene is callable depends on a
    start codon and a ribosome-binding site falling in the right place, not only
    on how much room there is -- so inside the turnover band the verdict
    alternates, and the list does not collapse to two cut points."""
    gold = _gold()
    assert gold["verdict_is_monotone_in_length"] is False
    ordered = [gold["verdict_by_length"][label]
               for label in sorted(gold["verdict_by_length"], key=int)]
    changes = sum(1 for a, b in zip(ordered, ordered[1:]) if a != b)
    assert changes > len(set(ordered)) - 1, ordered


# ------------------------------------------------- the build's own refusals
#
# One case per refusal in `_assert_has_content`, which the refusal contract in
# test_build_refusals.py requires and checks. Each mutation is the smallest one
# that reaches its own refusal: a refusal sitting behind an earlier one would
# never fire, so the ORDER of these is part of the build's contract and the cases
# below are written against it.


def _mutated(**overrides) -> dict:
    gold = copy.deepcopy(_gold())
    gold.update(overrides)
    return gold


def _quiet_contig() -> dict:
    return _mutated(baseline_calls=10)


def _too_few_targets() -> dict:
    gold = copy.deepcopy(_gold())
    gold["target_length_summary"]["targets"] = 3
    return gold


def _narrow_ladder() -> dict:
    gold = copy.deepcopy(_gold())
    summary = gold["target_length_summary"]
    summary["longest"] = summary["shortest"] * 5
    return gold


def _one_gene_audited_twice() -> dict:
    gold = copy.deepcopy(_gold())
    labels = sorted(gold["target_table"], key=int)
    first, second = gold["target_table"][labels[0]], gold["target_table"][labels[1]]
    for field in ("start", "end", "strand"):
        second[field] = first[field]
    return gold


def _deletion_outside_the_gene() -> dict:
    gold = copy.deepcopy(_gold())
    row = gold["target_table"][sorted(gold["target_table"], key=int)[0]]
    row["deleted_at"] = row["start"] - 1
    return gold


def _mismatched_tables() -> dict:
    gold = copy.deepcopy(_gold())
    label = sorted(gold["response_table"], key=int)[0]
    gold["response_table"]["not_a_declared_length"] = \
        gold["response_table"].pop(label)
    return gold


def _every_target_survives() -> dict:
    gold = copy.deepcopy(_gold())
    for row in gold["response_table"].values():
        row["verdict"] = "unchanged"
    gold["verdict_census"] = {"unchanged": len(gold["response_table"]),
                              "split_into_two": 0, "truncated_to_one": 0,
                              "abolished": 0}
    return gold


def _one_verdict_only() -> dict:
    gold = copy.deepcopy(_gold())
    for row in gold["response_table"].values():
        row["verdict"] = "abolished"
    gold["verdict_census"] = {"unchanged": 0, "split_into_two": 0,
                              "truncated_to_one": 0,
                              "abolished": len(gold["response_table"])}
    return gold


def _census(**counts) -> dict:
    """A gold whose verdicts are exactly the given counts, rows and census
    agreeing, so the refusal under test is the one that fires."""
    gold = copy.deepcopy(_gold())
    verdicts = [v for v, n in counts.items() for _ in range(n)]
    rows = sorted(gold["response_table"], key=int)
    assert len(verdicts) == len(rows)
    for label, verdict in zip(rows, verdicts):
        gold["response_table"][label]["verdict"] = verdict
    gold["verdict_census"] = {v: counts.get(v, 0) for v in _rules().verdicts}
    gold["verdict_by_length"] = {label: gold["response_table"][label]["verdict"]
                                 for label in rows}
    gold["verdict_is_monotone_in_length"] = False
    return gold


def _nothing_splits() -> dict:
    n = len(_gold()["response_table"])
    return _census(unchanged=1, truncated_to_one=1, abolished=n - 2)


def _nothing_is_abolished() -> dict:
    n = len(_gold()["response_table"])
    return _census(unchanged=1, truncated_to_one=1, split_into_two=n - 2)


def _miscounted_census() -> dict:
    gold = copy.deepcopy(_gold())
    gold["verdict_census"]["split_into_two"] += 1
    return gold


def _short_sweep() -> dict:
    gold = copy.deepcopy(_gold())
    keep = sorted(gold["offset_sweep"])[:3]
    gold["offset_sweep"] = {k: gold["offset_sweep"][k] for k in keep}
    return gold


def _sweep_with_a_dead_offset() -> dict:
    gold = copy.deepcopy(_gold())
    row = gold["offset_sweep"][sorted(gold["offset_sweep"])[0]]
    row["new_call_lengths"] = []
    return gold


def _sweep_that_does_not_move() -> dict:
    gold = copy.deepcopy(_gold())
    for row in gold["offset_sweep"].values():
        row["new_call_lengths"] = [500, 700]
    return gold


def _sweep_in_the_wrong_order() -> dict:
    gold = copy.deepcopy(_gold())
    for n, key in enumerate(sorted(gold["offset_sweep"], key=float)):
        gold["offset_sweep"][key]["new_call_lengths"] = [100 - n, 700]
    return gold


def _no_noise_arm() -> dict:
    return _mutated(noise_arm={})


def _noise_arm_that_explains_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    for row in gold["noise_arm"].values():
        row["noise_to_signal"] = None
    return gold


@pytest.mark.parametrize("mutate,pattern", [
    (_quiet_contig, "a record this quiet is an anecdote"),
    (_too_few_targets, "cannot span a length range"),
    (_narrow_ladder, "less than a factor of twenty"),
    (_one_gene_audited_twice, "two declared lengths selected the same call"),
    (_deletion_outside_the_gene, "is outside the gene"),
    (_mismatched_tables, "describe different gene sets"),
    (_every_target_survives, "every target survives its own frameshift"),
    (_one_verdict_only, "is not separating"),
    (_nothing_splits, "no target splits in two"),
    (_nothing_is_abolished, "no target is abolished"),
    (_miscounted_census, "the verdict census counts"),
    (_short_sweep, "cannot show where the split tracks the indel"),
    (_sweep_with_a_dead_offset, "some offset produces no overlapping call"),
    (_sweep_that_does_not_move, "lengths repeat"),
    (_sweep_in_the_wrong_order, "does not grow as the indel moves"),
    (_no_noise_arm, "the noise arm is empty"),
    (_noise_arm_that_explains_nothing, "explains any call"),
])
def test_the_build_refuses_gold_that_says_nothing(mutate, pattern):
    with pytest.raises(FrameshiftError, match=pattern):
        _assert_has_content(mutate())


def test_the_real_gold_passes_every_refusal():
    """The other half of the contract: seventeen mutations above raise, and the
    shipped gold raises nothing."""
    _assert_has_content(_gold())


def test_a_monotone_ladder_is_refused():
    """The guard that keeps R3 from being a threshold wearing a list's clothing.
    A verdict list in contiguous blocks is recoverable from its cut points, so
    R3 would grade one number restated once per target."""
    gold = dict(_gold())
    gold["verdict_is_monotone_in_length"] = True
    with pytest.raises(FrameshiftError, match="collapses to a threshold"):
        _assert_has_content(gold)


def test_the_monotone_test_is_contiguity_and_not_sortedness(shift):
    """Monotone here means each verdict occupies one contiguous block. A list
    whose verdicts happen to be in alphabetical order is not the claim."""
    blocks = {"1": {"verdict": "abolished"}, "2": {"verdict": "abolished"},
              "3": {"verdict": "truncated_to_one"},
              "4": {"verdict": "split_into_two"}}
    assert shift.verdict_is_monotone_in_length(blocks) is True
    alternating = {"1": {"verdict": "abolished"},
                   "2": {"verdict": "truncated_to_one"},
                   "3": {"verdict": "abolished"},
                   "4": {"verdict": "split_into_two"}}
    assert shift.verdict_is_monotone_in_length(alternating) is False


def test_the_ladder_is_read_as_lengths_and_not_as_labels(shift):
    """The labels are the declared lengths as strings, and sorting them as text
    puts 10000 between 1200 and 150. The monotone test is the one computation
    that depends on the order, so it sorts by integer -- and this is a case where
    the two orders give different answers.

    Not a claim about gold.json, which is canonically key-sorted like every other
    artifact here; a claim about the function that produces the ladder.
    """
    rows = {"2": {"verdict": "abolished"},
            "3": {"verdict": "split_into_two"},
            "10": {"verdict": "abolished"}}
    assert shift.verdict_is_monotone_in_length(rows) is False
    by_label = {k: rows[k] for k in sorted(rows)}       # "10", "2", "3"
    assert [r["verdict"] for r in by_label.values()] != \
        [r["verdict"] for r in rows.values()]
    assert list(shift.verdict_by_length(rows)) == ["2", "3", "10"]


# --------------------------------------------------------- the counterfactual


def test_the_offset_sweep_grades_the_split_point_and_not_the_verdict():
    """The catalogued R4 asked for the new call as the indel moves. Measured, the
    verdict does not move: every declared offset on a long gene gives two
    overlapping ORFs. What moves is where the gene splits."""
    gold = _gold()
    verdicts = {row["verdict"] for row in gold["offset_sweep"].values()}
    assert verdicts == {"split_into_two"}
    excluded = " ".join(str(b["field"]) for b in _task()["excluded_from_grading"])
    assert "verdict as a function of the indel" in excluded


def test_the_split_tracks_the_indel():
    """The first fragment grows as the deletion walks along the gene, and the
    lengths are distinct. If they repeated, the axis would measure nothing."""
    gold = _gold()
    rows = [gold["offset_sweep"][f"{f:g}"] for f in _rules().offset_fractions]
    firsts = [row["new_call_lengths"][0] for row in rows]
    assert firsts == sorted(firsts)
    assert len(set(firsts)) == len(firsts)
    into = [row["bases_into_the_gene"] for row in rows]
    assert into == sorted(into)


def test_the_sweep_shares_its_midpoint_with_the_rest_of_the_campaign():
    """The baseline fraction is on the offset grid, so the sweep and the response
    table agree at one point and a disagreement there is a bug and not an axis."""
    gold, rules = _gold(), _rules()
    assert rules.indel_fraction in rules.offset_fractions
    label = str(rules.offset_target_length)
    shared = gold["offset_sweep"][f"{rules.indel_fraction:g}"]
    row = gold["response_table"][label]
    assert shared["deleted_at"] == row["deleted_at"]
    assert shared["verdict"] == row["verdict"]
    assert shared["new_call_lengths"] == row["new_call_lengths"]


def test_the_noise_arm_measures_the_mode_choice_rather_than_asserting_it():
    """The campaign pins meta mode, and this is why: in the retraining mode a
    single deleted base moves an order of magnitude more calls elsewhere than the
    broken gene itself explains."""
    gold = _gold()
    ratios = [row["noise_to_signal"] for row in gold["noise_arm"].values()
              if row["noise_to_signal"] is not None]
    assert ratios
    assert max(ratios) > 10.0, gold["noise_arm"]


def test_the_baseline_mode_does_not_move_calls_the_deletion_cannot_explain():
    """The other half of the comparison, and the reason the arm is a ratio: in
    the pinned mode the collateral is a handful of calls across the whole
    chromosome, not fifty."""
    gold = _gold()
    collateral = max(row["other_calls_lost"] + row["other_calls_gained"]
                     for row in gold["response_table"].values())
    noise = min(row["other_calls_moved"] for row in gold["noise_arm"].values())
    assert collateral < noise, (collateral, noise)


def test_a_noise_arm_that_is_quieter_than_the_signal_is_refused():
    gold = dict(_gold())
    gold["noise_arm"] = {label: {**row, "noise_to_signal": 0.5}
                         for label, row in gold["noise_arm"].items()}
    with pytest.raises(FrameshiftError, match="not the noisy mode"):
        _assert_has_content(gold)


def test_the_noise_arm_spans_the_ladder():
    """Four targets rather than the whole ladder, one from each part of it, so a
    ratio that depended on gene length would show."""
    gold, rules = _gold(), _rules()
    assert set(gold["noise_arm"]) == {str(n) for n in rules.noise_target_lengths}
    assert set(rules.noise_target_lengths) <= set(rules.target_lengths)
    assert min(rules.noise_target_lengths) < 170
    assert max(rules.noise_target_lengths) >= 20000


# ---------------------------------------------------------------- the ladder


def test_the_oracle_clears_every_rung(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "frameshift_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["score"] == 1.0
    assert result["depth"] == 4
    assert result["monotonic"] is True
    assert result["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_an_empty_submission_clears_nothing(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["depth"] == 0


def test_the_warm_bundles_hand_over_exactly_the_preceding_rung():
    gold = _gold()
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", RESPONSE_KEYS, "responses.json"),
                             ("r4", VERDICT_KEYS, "verdicts.json")):
        payload = json.loads(
            (CAMPAIGN / "gold" / "warm" / rung / name).read_text())
        assert set(payload) == set(keys)
        assert all(payload[k] == gold[k] for k in keys)


def test_the_rung_key_sets_do_not_overlap():
    sets = [set(INVENTORY_KEYS), set(RESPONSE_KEYS), set(VERDICT_KEYS),
            set(COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right


def test_the_echoes_are_reported_and_not_graded():
    graded = {c["name"] for r in _spec()["rungs"] for c in r["components"]}
    gold = _gold()
    for echo in ("accession", "baseline_mode", "noise_mode", "indel_fraction",
                 "fixed_flags"):
        assert echo not in graded
        assert echo in gold


def test_the_tool_identity_is_graded_as_one_claim():
    """The version, the mode, the mode held back and the flags are one claim
    about what made the calls, so they are graded together rather than four
    times."""
    gold = _gold()
    assert set(gold["environment"]) == {"prodigal_version", "baseline_mode",
                                        "noise_mode", "fixed_flags"}
    graded = {c["name"] for r in _spec()["rungs"] for c in r["components"]}
    assert "environment" in graded
    assert "prodigal_version" not in graded


def test_this_campaign_pins_no_external_resource():
    """An ab initio caller reads only the sequence it is given, so the version
    and the mode are the whole of the pin."""
    assert _task()["constraints"]["external_resources"] == []


def test_the_judgement_the_catalog_asked_for_is_excluded():
    """The catalogued R3 was "which locus is a pseudogene and why". That has no
    oracle, and scoring a system against a curator's reading is the human
    judgement this benchmark keeps out of grading."""
    excluded = [str(b["field"]) for b in _task()["excluded_from_grading"]]
    assert any("pseudogene" in field for field in excluded)
    assert "planted" in _task()["structural_elements"]


# --------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    copy = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, copy)
    shutil.rmtree(copy / "gold")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "frameshift_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    assert measure(copy, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_keys_every_table_by_the_declared_vocabulary(tmp_path):
    """A submission that invented its own labels scores zero rather than being
    read on its own terms."""
    report = render(_gold())
    report["response_table"] = {f"gene_{k}": v
                                for k, v in report["response_table"].items()}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "frameshift_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["depth"] == 1


def test_a_missing_length_list_is_not_equal_to_an_empty_one(tmp_path):
    """`new_call_lengths` is legitimately empty for an abolished gene, so a
    submission that omitted the field must not read as that answer."""
    report = render(_gold())
    for row in report["response_table"].values():
        row.pop("new_call_lengths")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "frameshift_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["depth"] == 1


def test_a_null_ratio_is_preserved_and_not_read_as_zero(tmp_path):
    """The noise arm's ratio is null where a target explains no call, and null is
    not the same answer as a ratio of zero."""
    report = render(_gold())
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "frameshift_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    measured = measure(CAMPAIGN, submission)
    for label, row in measured["noise_arm"].items():
        assert row["noise_to_signal"] == _gold()["noise_arm"][label]["noise_to_signal"]

    report["noise_arm"] = {label: {**row, "noise_to_signal": 0}
                           for label, row in report["noise_arm"].items()}
    (submission / "frameshift_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    assert grade(measure(CAMPAIGN, submission), _gold(), _spec())["depth"] == 3


# --------------------------------------------------- the unlifted probe


def test_the_lifted_run_moves_nothing_the_deletion_cannot_explain():
    """What the lift buys, read off gold: in the pinned mode the mutated run
    matches the baseline everywhere except the gene that was broken. That is the
    baseline the unlifted probe is measured against, and it is why the collateral
    columns are a diagnostic rather than noise."""
    gold = _gold()
    noisy = sorted(label for label, row in gold["response_table"].items()
                   if row["other_calls_lost"] or row["other_calls_gained"])
    assert len(noisy) == 2, noisy
    assert all(gold["response_table"][label]["other_calls_lost"] <= 1
               for label in noisy)
    assert all(row["calls_after"] > 500
               for row in gold["response_table"].values())


#: The two targets whose verdict the missing lift flips, measured across the whole
#: ladder: 19 of the 21 survive it. Both are in the dense band, and both flip for
#: the same reason -- a call already inside the target's span is displaced too, so
#: it counts as new. Two rows rather than twenty-one are run here, because each
#: row costs a caller pass over 8.67 Mb; the ladder-wide figure is in
#: docs/npbench-c.md, measured by this same construction.
UNLIFTED_FLIPS = {"240": "truncated_to_one", "250": "split_into_two"}


@image_present
@pytest.mark.parametrize("label,wrong_verdict", sorted(UNLIFTED_FLIPS.items()))
def test_the_missing_lift_flips_these_two_verdicts(label, wrong_verdict):
    """The part of the probe that is NOT robust, pinned per row.

    Skipping the lift leaves 19 of 21 verdicts intact, because overlap with a
    gene's span survives a one-base displacement. These two do not survive it,
    and a change in either direction is a change in what the probe measures.
    """
    params = fsa.load_params(CAMPAIGN)
    shift = Frameshift.load(CAMPAIGN, params["contig_input"])
    rules = shift.rules
    baseline = shift.calls(rules.baseline_mode)
    target = shift.targets(baseline)[label]
    position = shift.indel_position(target)
    unlifted = tuple(sorted(
        (dataclasses.replace(
            call,
            start=call.start - (1 if call.start > position else 0),
            end=call.end - (1 if call.end > position else 0))
         for call in shift.calls(rules.baseline_mode, position)),
        key=lambda c: c.key))
    block = shift.response(baseline, target, unlifted)

    assert block["verdict"] == wrong_verdict
    assert block["verdict"] != _gold()["verdict_by_length"][label]
    assert block["other_calls_lost"] > 1000


@image_present
def test_the_unlifted_probe_clears_the_first_rung_and_fails_the_second(tmp_path):
    """The probe's whole point: it is wrong at R2 and indistinguishable from a
    correct run at R1, because R1 is about the contig and the untouched baseline
    and needs no lift."""
    params = fsa.load_params(CAMPAIGN)
    shift = Frameshift.load(CAMPAIGN, params["contig_input"])
    rules = shift.rules
    baseline = shift.calls(rules.baseline_mode)
    targets = shift.targets(baseline)

    label = str(rules.offset_target_length)
    target = targets[label]
    position = shift.indel_position(target)
    unlifted = tuple(sorted(
        (dataclasses.replace(
            call,
            start=call.start - (1 if call.start > position else 0),
            end=call.end - (1 if call.end > position else 0))
         for call in shift.calls(rules.baseline_mode, position)),
        key=lambda c: c.key))
    block = shift.response(baseline, target, unlifted)

    truth = _gold()["response_table"][label]
    assert truth["verdict"] == "split_into_two"

    # The VERDICT survives, and that is the probe's point rather than a weakness:
    # overlap with a 22 kb span is robust to a one-base displacement, so the
    # verdict column reads correct while everything that compares two runs is
    # wrong. R2 is where the two coordinate systems have to meet.
    assert block["verdict"] == truth["verdict"]
    assert block["other_calls_lost"] > 1000 and truth["other_calls_lost"] == 0
    assert block["other_calls_gained"] > 1000 and truth["other_calls_gained"] == 0
    # And the lift's absence shows in the target's own row too: a spurious third
    # overlapping call, lengths short by one, coverage past the gene's end.
    assert block["overlapping_new_calls"] > truth["overlapping_new_calls"]
    assert block["new_call_lengths"] != truth["new_call_lengths"]
    assert (block["bases_of_the_target_still_called"]
            != truth["bases_of_the_target_still_called"])
