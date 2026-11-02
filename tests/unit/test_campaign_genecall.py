"""Gene calling under fragmentation: the flag matters more than the assembly.

Two things here produce a wrong answer without producing an error.

The caller reports positions on the FRAGMENT, and comparing those with
whole-contig positions requires adding each fragment's offset. Skip the lift and
almost every call looks lost and almost every one looks new, with counts that
stay in the right order of magnitude and tables the right shape. That is the
competence probe.

And the mode is a flag, not a detail: single mode trains on the input it is
given, so fragmenting the contig changes the training set in one mode and not
the other. Reading the campaign as being about assembly quality without checking
the mode column gets the finding backwards.

The caller is fast -- fifteen seconds over the chromosome -- so the suite runs it
where it needs it.
"""

from __future__ import annotations

import dataclasses
import json
import pathlib
import shutil

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.templates.gene_calling import analyse as gca
from npbench_c.templates.gene_calling.build import (
    COUNTERFACTUAL_KEYS,
    FRAGMENTATION_KEYS,
    INVENTORY_KEYS,
    TABLE_KEYS,
    _assert_has_content,
)
from npbench_c.templates.gene_calling.core import (
    Call,
    Calling,
    CallingError,
    Rules,
    read_contig,
)
from npbench_c.templates.gene_calling.measure import measure
from npbench_c.templates.gene_calling.solve import render
from npbench_c.tools.resources import tool_prefix

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "genecall-scoelicolor-01")
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


@pytest.fixture(scope="module")
def calling() -> Calling:
    params = gca.load_params(CAMPAIGN)
    return Calling.load(CAMPAIGN, params["contig_input"])


# ------------------------------------------------------------ what is shipped


def test_no_annotation_is_shipped():
    """A gene-calling campaign cannot hand over a record carrying CDS features,
    because that is the answer. The detection campaign on this same accession
    ships the annotated GenBank; this one ships the bases."""
    shipped = sorted(p.name for p in (CAMPAIGN / "inputs").iterdir())
    assert shipped == ["contig.fna.gz", "provenance.json"]
    provenance = json.loads(
        (CAMPAIGN / "inputs" / "provenance.json").read_text())
    assert "CDS features and none of them is shipped" in \
        provenance["annotation_discarded"]
    excluded = " ".join(str(b["field"]) for b in _task()["excluded_from_grading"])
    assert "curated CDS annotation" in excluded


def test_the_provenance_does_not_carry_the_graded_values():
    """`contig_length` and `contig_sha256` are graded at R1, so a provenance that
    stated them would make both free components."""
    provenance = json.loads(
        (CAMPAIGN / "inputs" / "provenance.json").read_text())
    assert "length" not in provenance
    assert "sequence_sha256" not in provenance
    gold = _gold()
    flat = json.dumps(provenance)
    assert str(gold["contig_length"]) not in flat
    assert gold["contig_sha256"] not in flat


def test_the_input_is_one_unfragmented_record(calling):
    assert calling.length == _gold()["contig_length"]
    assert set(calling.sequence) <= set("ACGTN")


def test_a_multi_record_input_is_refused(tmp_path):
    """The campaign's subject is one contig becoming several, so an input that is
    already fragmented has no unfragmented arm to compare against."""
    path = tmp_path / "two.fna"
    path.write_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(CallingError, match="more than one record"):
        read_contig(path)


def test_a_sequence_outside_the_alphabet_is_refused(tmp_path):
    path = tmp_path / "odd.fna"
    path.write_text(">a\nACGTX\n")
    with pytest.raises(CallingError, match="declared alphabet"):
        read_contig(path)


def test_the_digest_is_over_the_bases_and_not_the_file(tmp_path, calling):
    """A FASTA rewrapped at another width is the same contig and has to
    fingerprint the same, or the digest is a claim about formatting."""
    rewrapped = tmp_path / "rewrapped.fna"
    seq = calling.sequence
    rewrapped.write_text(f">{calling.accession}\n"
                         + "\n".join(seq[i:i + 61] for i in range(0, len(seq), 61))
                         + "\n")
    identifier, bases = read_contig(rewrapped)
    assert identifier == calling.accession and bases == seq


# ------------------------------------------------------ the declared rule


def test_the_breakpoint_rule_partitions_the_contig(calling):
    rules = calling.rules
    for count in rules.fragment_counts:
        cuts = rules.breakpoints(calling.length, count)
        assert len(cuts) == count - 1
        assert list(cuts) == sorted(cuts)
        bounds = list(zip([0, *cuts], [*cuts, calling.length]))
        assert sum(hi - lo for lo, hi in bounds) == calling.length
        assert all(hi > lo for lo, hi in bounds)


def test_the_breakpoint_census_matches_the_rule():
    gold = _gold()
    for count, cuts in gold["breakpoint_census"].items():
        assert cuts == int(count) - 1


def test_the_fragment_fasta_reconstructs_the_contig(calling):
    """And carries no offset in its headers: the offsets follow from the declared
    rule, and putting them in the file would hand over the coordinate
    contract."""
    text = calling.fasta(calling.rules.baseline_fragments)
    headers = [l for l in text.splitlines() if l.startswith(">")]
    assert len(headers) == calling.rules.baseline_fragments
    assert all("offset" not in h for h in headers)
    rebuilt = "".join(l.strip() for l in text.splitlines() if not l.startswith(">"))
    assert rebuilt == calling.sequence


def test_an_undeclared_mode_or_fragment_count_is_refused(calling):
    with pytest.raises(CallingError, match="not a declared mode"):
        calling.calls("anonymous", 1)
    with pytest.raises(CallingError, match="not on the declared fragment grid"):
        calling.calls(None, 7)


def test_a_fragment_count_the_contig_cannot_carry_is_refused(calling):
    with pytest.raises(CallingError, match="not an assembly"):
        calling.rules.breakpoints(calling.length, 0)
    with pytest.raises(CallingError, match="shorter than a codon"):
        calling.rules.breakpoints(10, 20)


# ---------------------------------------------------- the coordinate contract


def test_the_loss_verdict_reads_the_two_conventions_together(calling):
    """A cut is a 0-based offset and a call is 1-based inclusive, so the cut
    falls strictly inside the call when `start <= cut < end`. Off by one either
    way and a gene spanning a breakpoint is classified as merely near it."""
    call = Call(start=100, end=200, strand="+", partial="00",
                start_type="ATG", fragment="frag_0001")
    assert calling.loss_verdict(call, [150]) == "crosses_a_cut"
    assert calling.loss_verdict(call, [100]) == "crosses_a_cut"
    # a cut at the call's last base is past it: the call ends there intact
    assert calling.loss_verdict(call, [200]) == "near_a_cut"
    assert calling.loss_verdict(call, [99]) == "near_a_cut"
    window = calling.rules.locality_window
    assert calling.loss_verdict(call, [200 + window]) == "near_a_cut"
    assert calling.loss_verdict(
        call, [200 + window + 1]) == "distant_from_every_cut"


def test_a_loss_cannot_be_classified_without_a_cut(calling):
    call = Call(start=100, end=200, strand="+", partial="00",
                start_type="ATG", fragment="frag_0001")
    with pytest.raises(CallingError, match="no cut to classify"):
        calling.loss_verdict(call, [])


def test_the_call_table_is_in_parent_coordinates():
    """Every call on the whole contig lies inside it, and the table is sorted by
    the declared key."""
    gold = _gold()
    rows = gold["call_table"]
    assert len(rows) == gold["call_census"]["calls"]
    assert all(1 <= r["start"] <= r["end"] <= gold["contig_length"] for r in rows)
    keys = [(r["start"], r["end"], r["strand"]) for r in rows]
    assert keys == sorted(keys)
    assert len(set(keys)) == len(keys)


def test_the_lost_table_rows_are_whole_contig_calls():
    """Lost calls come from the baseline side, so each must be in the call
    table -- a lost row carrying a fragment coordinate would be the unlifted
    defect showing up at R3."""
    gold = _gold()
    table = {(r["start"], r["end"], r["strand"]) for r in gold["call_table"]}
    lost = [(r["start"], r["end"], r["strand"])
            for r in gold["fragmented_lost_table"]]
    assert lost and all(k in table for k in lost)
    assert len(lost) == gold["fragmented_lost"]


# -------------------------------------------------------------- the inventory


def test_the_unfragmented_contig_has_no_partial_call():
    """A complete record has no edge for a call to run off. If any call were
    partial here the shipped input is not the whole chromosome, and R3 would be
    comparing two cut assemblies."""
    gold = _gold()
    assert gold["partial_census"]["00"] == gold["call_census"]["calls"]
    assert not any(n for flag, n in gold["partial_census"].items() if flag != "00")


def test_the_census_adds_up():
    gold = _gold()
    census = gold["call_census"]
    assert (census["on_the_forward_strand"] + census["on_the_reverse_strand"]
            == census["calls"])
    assert census["shortest"] <= census["median_length"] <= census["longest"]
    assert 0 < census["coding_fraction"] < 1
    assert census["coding_bases"] == sum(
        r["end"] - r["start"] + 1 for r in gold["call_table"])


def test_the_start_type_census_covers_every_call():
    gold = _gold()
    assert sum(gold["start_type_census"].values()) == gold["call_census"]["calls"]
    assert len(gold["start_type_census"]) >= 2


# ------------------------------------------------------- the fragmentation


def test_fragmentation_loses_and_gains_calls():
    gold = _gold()
    assert gold["fragmented_lost"] > 0
    assert gold["fragmented_gained"] > 0
    assert gold["fragmented_retained"] > 0
    assert (gold["fragmented_retained"] + gold["fragmented_lost"]
            == gold["fragmented_baseline_calls"])
    assert (gold["fragmented_retained"] + gold["fragmented_gained"]
            == gold["fragmented_fragmented_calls"])


def test_most_losses_are_geometry_and_some_are_not():
    """The campaign's point. A call spanning a breakpoint cannot be called whole
    -- that is geometry and it dominates. A call no breakpoint comes near is the
    caller reacting to a changed training set, and it is a different finding."""
    verdicts = _gold()["fragmented_loss_verdicts"]
    assert verdicts["crosses_a_cut"] > 0
    assert verdicts["crosses_a_cut"] > verdicts["near_a_cut"]
    assert verdicts["distant_from_every_cut"] > 0
    assert sum(verdicts.values()) == _gold()["fragmented_lost"]


def test_fragmentation_creates_partial_calls():
    """The flags the whole contig cannot produce: a call running off the left or
    the right end of its fragment."""
    gold = _gold()
    after = gold["fragmented_partial_census"]
    assert after["01"] > 0 and after["10"] > 0
    assert sum(after.values()) == gold["fragmented_fragmented_calls"]


# ------------------------------------------------------- the counterfactual


def test_the_grid_crosses_both_modes_with_the_series():
    gold, rules = _gold(), _rules()
    expected = {f"{m}|{n}" for m in rules.modes for n in rules.fragment_counts}
    assert set(gold["fragmentation_grid"]) == expected


def test_the_unfragmented_cells_are_the_controls():
    """Each compares a mode's whole contig against itself, so it must lose and
    gain nothing. Anything else means the comparison is reading two runs."""
    gold = _gold()
    for label, cell in gold["fragmentation_grid"].items():
        if label.endswith("|1"):
            assert cell["lost"] == 0 and cell["gained"] == 0
            assert cell["retained"] == cell["baseline_calls"]


def test_losses_rise_with_fragmentation_in_both_modes():
    gold, rules = _gold(), _rules()
    for mode in rules.modes:
        losses = [gold["fragmentation_grid"][f"{mode}|{n}"]["lost"]
                  for n in sorted(rules.fragment_counts)]
        assert losses == sorted(losses), (mode, losses)
        assert losses[0] == 0 and losses[-1] > 0


def test_only_the_retraining_mode_loses_calls_a_cut_cannot_explain():
    """The claim the mode axis exists for, measured per cell: single mode trains
    on its input, so fragmenting it changes calls nowhere near a cut. Meta mode
    does not train, and its distant column is empty at every fragmentation."""
    gold, rules = _gold(), _rules()
    distant = {
        mode: [gold["fragmentation_grid"][f"{mode}|{n}"]["loss_verdicts"]
               ["distant_from_every_cut"] for n in sorted(rules.fragment_counts)]
        for mode in rules.modes}
    assert sum(distant["single"]) > 0, distant
    assert sum(distant["meta"]) == 0, distant


def test_the_flag_moves_more_calls_than_the_assembly():
    """The headline, and it is a comparison rather than an assertion: the two
    modes on identical input bytes disagree on several times as many calls as
    cutting the contig into the baseline number of pieces costs."""
    gold = _gold()
    difference = gold["mode_difference"]["single_vs_meta"]
    moved = difference["only_in_single"] + difference["only_in_meta"]
    fragmentation = gold["fragmented_lost"] + gold["fragmented_gained"]
    assert moved > fragmentation
    assert difference["shared"] > 0


def test_the_grid_compares_each_mode_against_its_own_whole_contig():
    """Comparing meta-on-fragments against single-on-whole would add the mode
    difference to the fragmentation difference and report the sum as
    fragmentation."""
    gold = _gold()
    for label, cell in gold["fragmentation_grid"].items():
        mode = label.split("|")[0]
        whole = gold["fragmentation_grid"][f"{mode}|1"]
        assert cell["baseline_calls"] == whole["fragmented_calls"], label


# ---------------------------------------------------------------- the ladder


def test_the_oracle_clears_every_rung(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "calling_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(),
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["score"] == 1.0
    assert result["depth"] == 4
    assert result["monotonic"] is True
    assert result["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_an_empty_submission_clears_nothing(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    result = grade(measure(CAMPAIGN, submission), _gold(),
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["depth"] == 0


def test_the_warm_bundles_hand_over_exactly_the_preceding_rung():
    gold = _gold()
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", TABLE_KEYS, "call_table.json"),
                             ("r4", FRAGMENTATION_KEYS, "fragmentation.json")):
        payload = json.loads(
            (CAMPAIGN / "gold" / "warm" / rung / name).read_text())
        assert set(payload) == set(keys)
        assert all(payload[k] == gold[k] for k in keys)


def test_the_rung_key_sets_do_not_overlap():
    sets = [set(INVENTORY_KEYS), set(TABLE_KEYS), set(FRAGMENTATION_KEYS),
            set(COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right


def test_the_echoes_are_reported_and_not_graded():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    for echo in ("accession", "baseline_mode", "baseline_fragments",
                 "locality_window", "fixed_flags"):
        assert echo not in graded
        assert echo in _gold()


def test_this_campaign_pins_no_external_resource():
    """Every other tool-running campaign here reads a reference database from the
    image and pins it by path and digest. An ab initio caller reads only the
    sequence, so the version and the mode are the whole of the pin."""
    assert _task()["constraints"]["external_resources"] == []
    rules = json.loads(
        (CAMPAIGN / "reference" / "calling_rules.json").read_text())
    assert "external_resources" not in rules


# --------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    copy = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, copy)
    shutil.rmtree(copy / "gold")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "calling_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    assert measure(copy, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_projects_rows_without_reordering_them(tmp_path):
    gold = _gold()
    report = render(gold)
    report["call_table"] = list(reversed(report["call_table"]))
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "calling_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), gold,
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["depth"] == 1


# --------------------------------------------------------- the unlifted probe


def test_skipping_the_coordinate_lift_looks_plausible_and_is_wrong(calling):
    """The competence probe, computed here rather than asserted. Undo the lift on
    the fragmented run and the comparison reports almost everything as lost and
    almost everything as new -- counts in the right order of magnitude, tables
    the right shape, and nothing saying a convention was missed."""
    rules = calling.rules
    whole = calling.calls(rules.baseline_mode, 1)
    fragmented = calling.calls(rules.baseline_mode, rules.baseline_fragments)
    cuts = rules.breakpoints(calling.length, rules.baseline_fragments)
    offsets = calling.offsets(rules.baseline_fragments)
    unlifted = tuple(
        dataclasses.replace(call,
                            start=call.start - offsets[call.fragment],
                            end=call.end - offsets[call.fragment])
        for call in fragmented)

    honest = calling.reconcile(whole, fragmented, cuts)
    broken = calling.reconcile(whole, unlifted, cuts)
    assert broken["lost"] > 20 * honest["lost"]
    assert broken["gained"] > 20 * honest["gained"]
    # The shape is identical: same keys, same column names, same vocabulary.
    assert set(broken) == set(honest)
    assert set(broken["loss_verdicts"]) == set(honest["loss_verdicts"])


def test_the_task_declares_the_stub_levels_the_module_implements():
    from npbench_c.sweep.stubs.stub_genecall import LEVELS
    declared = {level["level"] for level in _task()["harness"]["stub_levels"]}
    assert declared == set(LEVELS)


# ------------------------------------------------------------ the gold refusals


def test_the_build_refuses_a_contig_too_quiet_to_fragment():
    gold = _gold()
    with pytest.raises(CallingError, match="an anecdote"):
        _assert_has_content({**gold, "call_census": {
            **gold["call_census"], "calls": 100}})


def test_the_build_refuses_a_single_stranded_or_fully_coding_contig():
    gold = _gold()
    census = gold["call_census"]
    with pytest.raises(CallingError, match="a free column"):
        _assert_has_content({**gold, "call_census": {
            **census, "on_the_reverse_strand": 0}})
    with pytest.raises(CallingError, match="not an assembly"):
        _assert_has_content({**gold, "call_census": {
            **census, "coding_fraction": 1.0}})


def test_the_build_refuses_an_input_that_is_already_fragmented():
    gold = _gold()
    census = gold["call_census"]
    partial = {**gold["partial_census"], "00": census["calls"] - 5, "10": 5}
    with pytest.raises(CallingError, match="not the whole record"):
        _assert_has_content({**gold, "partial_census": partial})


def test_the_build_refuses_a_constant_start_type_census():
    gold = _gold()
    with pytest.raises(CallingError, match="constant dressed up"):
        _assert_has_content({**gold, "start_type_census": {
            "ATG": gold["call_census"]["calls"]}})


def test_the_build_refuses_a_fragmentation_that_changes_nothing():
    gold = _gold()
    with pytest.raises(CallingError, match="subject does not exist"):
        _assert_has_content({**gold, "fragmented_lost": 0})
    with pytest.raises(CallingError, match="comparison is one-sided"):
        _assert_has_content({**gold, "fragmented_gained": 0})
    with pytest.raises(CallingError, match="broken run rather than"):
        _assert_has_content({**gold, "fragmented_retained": 0})


def test_the_build_refuses_losses_the_geometry_cannot_produce():
    """No lost call spanning a breakpoint means the coordinate lift or the
    breakpoint rule is wrong -- and both failures produce plausible counts."""
    gold = _gold()
    verdicts = gold["fragmented_loss_verdicts"]
    moved = {**verdicts, "crosses_a_cut": 0,
             "near_a_cut": verdicts["near_a_cut"] + verdicts["crosses_a_cut"]}
    with pytest.raises(CallingError, match="the geometry guarantees"):
        _assert_has_content({**gold, "fragmented_loss_verdicts": moved})


def test_the_build_refuses_a_collapsed_loss_vocabulary():
    gold = _gold()
    total = gold["fragmented_lost"]
    sole = {k: (total if k == "crosses_a_cut" else 0)
            for k in gold["fragmented_loss_verdicts"]}
    with pytest.raises(CallingError, match="not separating anything"):
        _assert_has_content({**gold, "fragmented_loss_verdicts": sole})


def test_the_build_refuses_a_lost_table_that_disagrees_with_its_count():
    gold = _gold()
    with pytest.raises(CallingError, match="one of them is wrong"):
        _assert_has_content({**gold,
                             "fragmented_lost_table": gold["fragmented_lost_table"][:-1]})


def test_the_build_refuses_a_grid_too_small_or_without_a_control():
    gold = _gold()
    label, cell = next(iter(gold["fragmentation_grid"].items()))
    with pytest.raises(CallingError, match="cannot cross two modes"):
        _assert_has_content({**gold, "fragmentation_grid": {label: cell}})
    cut_only = {l: c for l, c in gold["fragmentation_grid"].items()
                if not l.endswith("|1")}
    with pytest.raises(CallingError, match="nothing in it is a control"):
        _assert_has_content({**gold, "fragmentation_grid": cut_only})


def test_the_build_refuses_a_control_cell_that_moved():
    gold = _gold()
    grid = {l: ({**c, "lost": 2} if l.endswith("|1") else c)
            for l, c in gold["fragmentation_grid"].items()}
    with pytest.raises(CallingError, match="impossible unless"):
        _assert_has_content({**gold, "fragmentation_grid": grid})


def test_the_build_refuses_a_fragmentation_series_that_is_not_a_series():
    gold = _gold()
    grid = {l: ({**c, "lost": 9} if not l.endswith("|1") else c)
            for l, c in gold["fragmentation_grid"].items()}
    with pytest.raises(CallingError, match="is not a series"):
        _assert_has_content({**gold, "fragmentation_grid": grid})


def test_the_build_refuses_modes_that_agree_on_everything():
    gold = _gold()
    with pytest.raises(CallingError, match="not a parameter on this contig"):
        _assert_has_content({**gold, "mode_difference": {
            "single_vs_meta": {"shared": 7000, "only_in_single": 0,
                               "only_in_meta": 0}}})
    with pytest.raises(CallingError, match="cannot separate the flag"):
        _assert_has_content({**gold, "mode_difference": {}})
