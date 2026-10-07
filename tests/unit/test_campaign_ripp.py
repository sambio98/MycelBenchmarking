"""RiPP precursor annotation: three records of one peptide, checked against each
other.

MIBiG gives a precursor a gene, a core sequence and sometimes a cleavage
coordinate; the reference record gives the translation. Four things here produce
a wrong answer without producing an error, and each has a test.

`core_sequence` arrives in four JSON shapes, and `core_sequence[0]` -- the
natural thing to write -- silently takes one residue where the field is a bare
string. Some cores are in a different case from the translations, so a reader
that does not fold case reports them absent. One shape is the printed form of a
list, brackets and quotes included, which is not a peptide and must not be read
as one. And the cleavage coordinate has no documented convention, so a reader
that assumes one restates its own assumption.

The tool is fast here -- two scans of an eight-model panel over sequences of
fifty residues -- so the suite runs it where it needs it.
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
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.ripp_precursor import analyse as rpa
from npbench_c.templates.ripp_precursor.build import (
    COUNTERFACTUAL_KEYS,
    EVIDENCE_KEYS,
    INVENTORY_KEYS,
    LOCALISATION_KEYS,
    _assert_has_content,
)
from npbench_c.templates.ripp_precursor.core import (
    Corpus,
    PrecursorError,
    Rules,
    cds_translations,
)
from npbench_c.templates.ripp_precursor.measure import measure
from npbench_c.templates.ripp_precursor.solve import render
from npbench_c.tools.resources import tool_prefix

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "ripp-precursor-mibig-4_0-01")
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
def corpus() -> Corpus:
    params = rpa.load_params(CAMPAIGN)
    return Corpus.load(CAMPAIGN, params["annotations_input"],
                       params["reference_records_input"])


# ----------------------------------------------------------- the GenBank read


def test_a_cds_translation_is_reachable_by_every_declared_name(corpus):
    """A precursor is named by protein accession in some entries and by gene name
    in others, and the audit cannot work unless both resolve."""
    by_shape = {"accession": 0, "name": 0}
    for record in corpus.records:
        key = "accession" if record.gene[:1].isupper() and "." in record.gene \
            else "name"
        by_shape[key] += 1
        assert record.translation, record.gene
    # Both naming styles have to be present or the qualifier fallback chain is
    # untested here and inherited from another instantiation.
    assert by_shape["accession"] > 0 and by_shape["name"] > 0, by_shape


def test_the_feature_table_ends_at_column_zero():
    """The parser's stop condition. A feature line starts at column 5 and a
    qualifier at 21, so a `startswith(" " * 6)` test ends the table at the first
    feature and finds nothing -- which is exactly how this was first written."""
    text = "\n".join([
        "LOCUS       X                 100 bp    DNA     linear",
        "FEATURES             Location/Qualifiers",
        "     CDS             1..30",
        '                     /protein_id="AAA00000.1"',
        '                     /translation="MGPVVVFDC"',
        "ORIGIN",
        "        1 atg",
        "//",
    ])
    assert cds_translations(text) == {"AAA00000.1": "MGPVVVFDC"}


def test_a_wrapped_translation_is_rejoined_without_spaces():
    text = "\n".join([
        "FEATURES             Location/Qualifiers",
        "     CDS             1..30",
        '                     /gene="bmbA"',
        '                     /translation="MGPVVVFDCMTADFLNDDPNNAELSALEMEELES',
        '                     WGAWDGEATS"',
        "ORIGIN",
    ])
    assert cds_translations(text)["bmbA"] == (
        "MGPVVVFDCMTADFLNDDPNNAELSALEMEELESWGAWDGEATS")


def test_a_cds_without_a_translation_is_not_indexed():
    text = "\n".join([
        "FEATURES             Location/Qualifiers",
        "     gene            1..30",
        '                     /gene="orfA"',
        "     CDS             1..30",
        '                     /gene="orfA"',
        '                     /pseudo',
        "ORIGIN",
    ])
    assert cds_translations(text) == {}


# ---------------------------------------------------- the declared transforms


def test_every_declared_shape_occurs_in_the_corpus():
    """The ingestion policy is only worth declaring where the field is
    inconsistent, and the build refuses gold otherwise. This states which shapes
    it is inconsistent in."""
    shapes = _gold()["core_sequence_shape_census"]
    assert set(shapes) >= {"bare_string", "list_of_one", "list_of_several",
                           "printed_list_in_a_string"}
    assert all(n > 0 for n in shapes.values())


def test_ingestion_reads_the_two_legitimate_shapes_the_same_way():
    rules = Rules.load(CAMPAIGN)
    assert rules.ingest("GPVVVFDC") == ["GPVVVFDC"]
    assert rules.ingest(["GPVVVFDC"]) == ["GPVVVFDC"]
    assert rules.ingest(["TLPVPTLC", "TVPTLC"]) == ["TLPVPTLC", "TVPTLC"]


def test_case_folding_is_the_policy_and_not_the_comparison():
    """The baseline folds case and the strict policy does not, so a lower-case
    core is absent under one and located under the other. Folding inside the
    comparison instead would make the baseline its own counterfactual."""
    rules = Rules.load(CAMPAIGN)
    assert rules.ingest("ggplageeiggfnvpgisee", "p0_strict") == [
        "ggplageeiggfnvpgisee"]
    assert rules.ingest("ggplageeiggfnvpgisee", "p1_case_folded") == [
        "GGPLAGEEIGGFNVPGISEE"]

    sweep = _gold()["policy_sweep"]
    strict = sweep["p0_strict"]["verdict_census"]
    folded = sweep["p1_case_folded"]["verdict_census"]
    assert strict["core_absent"] > 0
    assert folded["core_absent"] == 0
    assert sweep["p1_case_folded"]["located"] > sweep["p0_strict"]["located"]


def test_a_printed_list_is_malformed_rather_than_parsed():
    """The field sometimes holds `['GGAGH...']` as a STRING. The baseline reads
    it as one core, whose brackets and quotes are not amino acids, so the record
    is malformed -- which is the honest verdict. Repairing it is R4's business."""
    rules = Rules.load(CAMPAIGN)
    printed = "['GGAGHVPEYFVGIGTPISFYG']"
    assert rules.ingest(printed, "p1_case_folded") == [printed.upper()]
    assert rules.well_formed(printed) is False
    assert rules.ingest(printed, "p2_repair_printed_list") == [
        "GGAGHVPEYFVGIGTPISFYG"]

    sweep = _gold()["policy_sweep"]
    assert sweep["p1_case_folded"]["verdict_census"]["malformed_core"] > 0
    assert sweep["p2_repair_printed_list"]["verdict_census"]["malformed_core"] == 0


def test_a_bracket_inside_a_core_is_not_reinterpreted():
    """The repair matches the WHOLE field against a printed list. A core that
    merely contains a bracket stays malformed, because reinterpreting it would be
    guessing at a different accident."""
    rules = Rules.load(CAMPAIGN)
    assert rules.ingest("GPVV['FDC", "p2_repair_printed_list") == ["GPVV['FDC"]
    assert rules.well_formed("GPVV['FDC") is False


def test_lower_case_is_well_formed_and_a_bracket_is_not():
    """Case is a normalisation; an alphabet violation is a defect. Conflating
    them reports a case difference as malformed data."""
    rules = Rules.load(CAMPAIGN)
    assert rules.well_formed("ggplageeiggfnvpgisee") is True
    assert rules.well_formed("GGPLAGEEIGGFNVPGISEE") is True
    assert rules.well_formed("['GG']") is False
    assert rules.well_formed("") is False


def test_a_non_string_core_sequence_is_refused_rather_than_coerced():
    rules = Rules.load(CAMPAIGN)
    with pytest.raises(PrecursorError, match="only a string or a list"):
        rules.ingest({"core": "GPVVVFDC"})
    with pytest.raises(PrecursorError, match="only a string or a list"):
        rules.ingest([1, 2, 3])


# -------------------------------------------------------- the localisation


def test_the_verdict_vocabulary_is_closed_and_populated(corpus):
    gold = _gold()
    rules = corpus.rules
    # gold.json is serialised with sorted keys, so the declared ORDER is not
    # observable there -- and does not need to be: a census is compared as a
    # mapping. What has to hold is that the vocabulary is exactly the declared
    # one, with no term invented and none dropped.
    assert set(gold["localisation_census"]) == set(rules.verdict_order)
    populated = [v for v, n in gold["localisation_census"].items() if n]
    assert len(populated) >= 4, populated


def test_an_unlocatable_core_reports_nulls_rather_than_a_guess():
    """Where the core has no single position the three derived coordinates are
    null. Filling them with the nearest match would make R3 reconcile a
    coordinate against a number nobody measured."""
    rows = _gold()["localisation_table"]
    unlocated = [r for r in rows if r["core_start"] is None]
    assert unlocated, "no record exercises the null path"
    for row in unlocated:
        assert row["leader_length"] is None and row["follower_length"] is None
        assert row["verdict"] in ("malformed_core", "core_absent",
                                 "core_ambiguous", "core_is_whole_precursor")
    for row in rows:
        if row["core_start"] is not None:
            assert row["leader_length"] == row["core_start"]
            assert row["occurrences"] == 1


def test_a_core_that_is_the_whole_precursor_is_its_own_verdict():
    """Some entries paste the entire precursor into core_sequence. It localises
    perfectly at offset zero, so a reader that only checks 'is it findable'
    reports a leader of length zero and a clean record."""
    gold = _gold()
    assert gold["localisation_census"]["core_is_whole_precursor"] > 0
    whole = [r for r in gold["localisation_table"]
             if r["verdict"] == "core_is_whole_precursor"]
    assert all(r["core_start"] is None for r in whole)


def test_the_corpus_holds_a_core_with_more_than_one_occurrence():
    """A tandem-repeat precursor carries its core twice, so 'where is the core'
    has no answer. Taking the first match would be a declared choice the campaign
    has not declared."""
    assert _gold()["localisation_census"]["core_ambiguous"] > 0


# ---------------------------------------------------- the coordinate convention


def test_the_convention_is_elected_decisively():
    gold = _gold()
    tally = gold["convention_tally"]
    ranked = sorted((b["agrees"] for b in tally.values()), reverse=True)
    assert ranked[0] > 0
    assert ranked[0] > ranked[1], "the election is a tie-break, not a finding"
    assert tally[gold["elected_convention"]]["agrees"] == ranked[0]
    assert any(b["disagrees"] for b in tally.values())


def test_the_elected_convention_agrees_with_every_derivable_record():
    """What the measurement actually says: on every record where both a
    coordinate and a position exist, they agree. The three that do not agree are
    the three with no position -- which is a statement about those annotations
    and not about the convention."""
    gold = _gold()
    elected = gold["convention_tally"][gold["elected_convention"]]
    assert elected["disagrees"] == 0
    assert elected["underivable"] > 0
    assert (elected["agrees"] + elected["underivable"]
            == gold["records_declaring_a_cleavage_location"])


def test_every_candidate_convention_is_scored(corpus):
    gold = _gold()
    assert set(gold["convention_tally"]) == set(corpus.rules.conventions)
    for block in gold["convention_tally"].values():
        assert (sum(block.values())
                == gold["records_declaring_a_cleavage_location"])


def test_the_convention_arithmetic_is_what_it_claims(corpus):
    rules = corpus.rules
    location = {"from": 39, "to": 40}
    assert rules.convention_start("to_as_0based_core_start", location) == 40
    assert rules.convention_start("to_as_1based_core_start", location) == 39
    assert rules.convention_start("from_as_0based_core_start", location) == 39
    assert rules.convention_start("from_as_1based_core_start", location) == 38
    with pytest.raises(PrecursorError, match="not a declared convention"):
        rules.convention_start("to_as_2based_core_start", location)


# --------------------------------------------------------- the domain evidence


def test_no_declared_family_marks_the_boundary():
    """The measurement R3 exists for. If some Pfam family's envelope stopped
    exactly where the leader does, 'use the domain call as the boundary' would be
    sound advice; measured here, none does."""
    gold = _gold()
    assert gold["records_with_a_panel_hit"] > 0
    assert gold["boundary_offset_summary"]["exactly_at_the_boundary"] == 0
    relations = gold["boundary_relation_census"]
    assert relations["core_only"] == 0
    assert relations["spans_boundary"] > 0


def test_every_panel_hit_is_placed_against_the_measured_boundary():
    gold = _gold()
    starts = {(r["bgc"], r["gene"]): r["core_start"]
              for r in gold["localisation_table"]}
    for row in gold["domain_table"]:
        start = starts[(row["bgc"], row["gene"])]
        if start is None:
            assert row["relation"] == "core_not_located"
            assert row["boundary_offset"] is None
            continue
        assert row["boundary_offset"] == row["env_to"] - start
        if row["env_to"] <= start:
            assert row["relation"] == "leader_only"
        elif row["env_from"] >= start + 1:
            assert row["relation"] == "core_only"
        else:
            assert row["relation"] == "spans_boundary"


def test_the_panel_census_covers_the_declared_panel(corpus):
    gold = _gold()
    assert set(gold["panel_census"]) == set(corpus.rules.panel)
    assert sum(gold["panel_census"].values()) == len(gold["domain_table"])


# ------------------------------------------------------- the counterfactuals


def test_the_policy_sweep_moves_and_every_policy_is_distinct(corpus):
    gold = _gold()
    assert set(gold["policy_sweep"]) == set(corpus.rules.policies)
    censuses = {json.dumps(b["verdict_census"], sort_keys=True)
                for b in gold["policy_sweep"].values()}
    assert len(censuses) == len(gold["policy_sweep"])
    located = [b["located"] for _, b in sorted(gold["policy_sweep"].items())]
    assert located == sorted(located), "the policies are ordered by strictness"


def test_the_leader_carries_the_domain_signal_for_most_records_and_not_all():
    """R4's finding, and its own control. If every call went with its leader the
    rung would have one answer for every record; if none did, the deletion would
    show nothing."""
    census = _gold()["deletion_census"]
    assert census["lost"] > 0 and census["retained"] > 0
    assert census["lost"] > census["retained"]
    assert census["unchanged_absent"] > 0, "no negative control"
    assert census["not_applicable"] > 0


def test_the_deletion_table_carries_only_the_records_that_moved_or_held():
    gold = _gold()
    assert {r["outcome"] for r in gold["deletion_table"]} <= {
        "retained", "changed", "lost", "gained"}
    assert len(gold["deletion_table"]) == (
        gold["deletion_census"]["retained"] + gold["deletion_census"]["changed"]
        + gold["deletion_census"]["lost"] + gold["deletion_census"]["gained"])
    for row in gold["deletion_table"]:
        if row["outcome"] == "lost":
            assert row["before"] and not row["after"]
        if row["outcome"] == "retained":
            assert row["before"] == row["after"] and row["before"]


def test_a_record_without_a_located_core_is_left_out_of_the_second_scan(corpus):
    """Passing it through whole would compare a full precursor against the
    others' cores, which is not the counterfactual the campaign declares."""
    table = corpus.localisation_table()
    fasta = corpus.fasta(table)
    located = {(r["bgc"], r["gene"]) for r in table
               if r["core_start"] is not None}
    headers = {line[1:] for line in fasta.splitlines() if line.startswith(">")}
    assert headers == {f"{bgc}|{gene}" for bgc, gene in located}
    assert len(headers) < len(corpus.records)


# ------------------------------------------------------------ the selection


def test_each_exclusion_clause_removes_something():
    census = _gold()["selection_census"]
    assert census["entries_excluded_by_status"] > 0
    assert census["entries_absent_from_reference_set"] > 0
    assert (census["corpus_entries"] + census["entries_excluded_by_status"]
            + census["entries_absent_from_reference_set"]
            - census["entries_excluded_by_both"]
            == census["entries_with_a_named_precursor"])


def test_the_two_exclusion_clauses_coincide_exactly_on_this_data():
    """A measured coupling worth stating: every named entry MIBiG has retired is
    also absent from the antiSMASH reference set, and vice versa. The counts are
    taken independently for exactly this reason -- in sequence, one of the two
    clauses would have reported zero and the coincidence would have been
    invisible."""
    census = _gold()["selection_census"]
    assert (census["entries_excluded_by_both"]
            == census["entries_excluded_by_status"]
            == census["entries_absent_from_reference_set"])


def test_the_reference_archive_ships_entries_the_audit_never_uses():
    """Shipping only the entries that pass the reference clause would make the
    count of those that do not readable off a listing instead of computed from
    the two resources -- a free component."""
    census = _gold()["selection_census"]
    assert census["reference_records_shipped"] > census["corpus_entries"]
    assert census["reference_records_shipped"] < census["entries_in_class"]


def test_the_library_is_not_redistributed():
    shipped = {p.name for p in (CAMPAIGN / "inputs").iterdir()}
    assert not any(name.endswith(".hmm") for name in shipped), shipped
    declared = json.loads(
        (CAMPAIGN / "reference" / "precursor_rules.json").read_text())
    pfam = declared["external_resources"]["pfam_a_hmm"]
    assert len(pfam["sha256"]) == 64
    assert "size, not licence" in pfam["note"]


# -------------------------------------------------------------- the ladder


def test_the_oracle_clears_every_rung(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "precursor_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(),
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["score"] == 1.0
    assert result["depth"] == 4
    assert result["monotonic"] is True
    # R1's chance level is nominal at 1.0 -- any system that reads the inputs
    # and runs the tool clears it -- so the floor that has to be low is the one
    # computed over R2 upward.
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
                             ("r3", LOCALISATION_KEYS, "localisation.json"),
                             ("r4", EVIDENCE_KEYS, "evidence.json")):
        payload = json.loads(
            (CAMPAIGN / "gold" / "warm" / rung / name).read_text())
        assert set(payload) == set(keys)
        assert all(payload[k] == gold[k] for k in keys)


def test_the_rung_key_sets_do_not_overlap():
    """A key handed over by one rung's bundle and graded at that same rung would
    be a free answer, and the readiness gate's own check is the backstop. This is
    the structural statement."""
    sets = [set(INVENTORY_KEYS), set(LOCALISATION_KEYS), set(EVIDENCE_KEYS),
            set(COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right


def test_the_echoes_of_the_rules_file_are_reported_and_not_graded():
    """Their value is fixed by the campaign's own declaration, so grading them
    adds weight to a conjunctive product and no information."""
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    for echo in ("mibig_release", "biosynthetic_class", "baseline_policy",
                 "threshold", "panel_size"):
        assert echo not in graded
        assert echo in _gold()
    excluded = " ".join(
        str(block["field"]) for block in _task()["excluded_from_grading"])
    assert "baseline_policy" in excluded


# --------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    """Run it against a campaign copy with gold removed. The measurer's key lists
    all come from the campaign's declared vocabularies, so this has to hold."""
    copy = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, copy)
    shutil.rmtree(copy / "gold")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "precursor_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    assert measure(copy, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_projects_rows_without_reordering_them(tmp_path):
    """Gold's row order is declared, so a submission that sorts differently has
    the order wrong and must not be silently repaired."""
    gold = _gold()
    report = render(gold)
    report["localisation_table"] = list(reversed(report["localisation_table"]))
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "precursor_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), gold,
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["depth"] == 1


def test_a_null_coordinate_survives_the_projection(tmp_path):
    """core_start is legitimately null in gold, so coercing it to a number would
    make a submission that answered null disagree with gold."""
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "precursor_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    measured = measure(CAMPAIGN, submission)
    assert any(row["core_start"] is None
               for row in measured["localisation_table"])
    assert any(row["boundary_offset"] is None
               for row in measured["domain_table"])


def test_a_censored_census_still_equals_gold(tmp_path):
    """Gold omits a shape with no members, so a submission that spells it out as
    zero is the same answer."""
    report = render(_gold())
    report["core_sequence_shape_census"] = {
        **report["core_sequence_shape_census"], "other": 0}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "precursor_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(),
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["depth"] == 4


# ------------------------------------------------------------ the gold refusals


def test_the_build_refuses_a_corpus_too_small_to_audit():
    gold = _gold()
    thin = {**gold, "selection_census": {**gold["selection_census"],
                                         "corpus_records": 9}}
    with pytest.raises(PrecursorError, match="anecdote"):
        _assert_has_content(thin)


def test_the_build_refuses_an_unexercised_selection_clause():
    gold = _gold()
    for field in ("entries_excluded_by_status",
                  "entries_absent_from_reference_set"):
        broken = {**gold, "selection_census": {**gold["selection_census"],
                                               field: 0}}
        with pytest.raises(PrecursorError, match="not exercised|different sets"):
            _assert_has_content(broken)


def test_the_build_refuses_a_consistent_field():
    gold = _gold()
    with pytest.raises(PrecursorError, match="field is consistent"):
        _assert_has_content({**gold, "core_sequence_shape_census":
                             {"bare_string": 67}})
    with pytest.raises(PrecursorError, match="same case"):
        _assert_has_content({**gold, "core_sequence_case_census":
                             {"upper": 67}})


def test_the_build_refuses_an_undecidable_or_tied_convention():
    gold = _gold()
    flat = {c: {"agrees": 0, "disagrees": 0, "underivable": 30}
            for c in gold["convention_tally"]}
    with pytest.raises(PrecursorError, match="cannot be reconciled"):
        _assert_has_content({**gold, "convention_tally": flat})

    tied = {c: {"agrees": 27, "disagrees": 3, "underivable": 0}
            for c in gold["convention_tally"]}
    with pytest.raises(PrecursorError, match="declared tie-break"):
        _assert_has_content({**gold, "convention_tally": tied})

    agreeing = {c: {"agrees": 30 - i, "disagrees": 0, "underivable": i}
                for i, c in enumerate(gold["convention_tally"])}
    with pytest.raises(PrecursorError, match="cannot tell the readings apart"):
        _assert_has_content({**gold, "convention_tally": agreeing})


def test_the_build_refuses_a_fully_located_corpus():
    gold = _gold()
    located = [{**row, "core_start": 1, "leader_length": 1,
                "follower_length": 1} for row in gold["localisation_table"]]
    with pytest.raises(PrecursorError, match="coordinate it cannot derive"):
        _assert_has_content({**gold, "localisation_table": located})


def test_the_build_refuses_an_empty_panel_or_one_sided_boundary():
    gold = _gold()
    with pytest.raises(PrecursorError, match="matches nothing"):
        _assert_has_content({**gold, "records_with_a_panel_hit": 0})
    one_sided = {r: (len(gold["domain_table"]) if r == "spans_boundary" else 0)
                 for r in gold["boundary_relation_census"]}
    with pytest.raises(PrecursorError, match="not distinguishing anything"):
        _assert_has_content({**gold, "boundary_relation_census": one_sided})


def test_the_build_refuses_a_policy_sweep_that_does_not_move():
    gold = _gold()
    same = {name: {**block,
                   "verdict_census": gold["localisation_census"]}
            for name, block in gold["policy_sweep"].items()}
    with pytest.raises(PrecursorError, match="two names"):
        _assert_has_content({**gold, "policy_sweep": same})


def test_the_build_refuses_a_deletion_counterfactual_with_one_answer():
    gold = _gold()
    census = gold["deletion_census"]
    with pytest.raises(PrecursorError, match="no negative control"):
        _assert_has_content({**gold, "deletion_census":
                             {**census, "unchanged_absent": 0}})
    with pytest.raises(PrecursorError, match="changes no domain call"):
        _assert_has_content({**gold, "deletion_census":
                             {**census, "lost": 0, "changed": 0, "gained": 0}})
    with pytest.raises(PrecursorError, match="free answer"):
        _assert_has_content({**gold, "deletion_census":
                             {**census, "retained": 0}})


# --------------------------------------------------------------- the stub probe


def test_the_first_element_shortcut_breaks_the_localisation(corpus):
    """`core_sequence[0]` is right for the list shapes and takes one residue
    where the field is a bare string. It raises no error on any record, and the
    broken rows still localise -- just ambiguously. This is the competence probe
    `stub-firstshape` runs."""
    broken = dataclasses.replace(
        corpus, records=tuple(
            dataclasses.replace(r, raw_core_sequence=r.raw_core_sequence[0])
            for r in corpus.records))
    table = broken.localisation_table()
    census = broken.localisation_census(table)
    assert census != _gold()["localisation_census"]
    # The symptom: a third of the corpus becomes a single residue that occurs
    # many times over, so it is reported as ambiguous rather than as a parse
    # failure.
    assert census["core_ambiguous"] > _gold()[
        "localisation_census"]["core_ambiguous"]
    assert len(table) == len(corpus.records)


def test_the_task_declares_the_stub_levels_the_module_implements():
    from npbench_c.sweep.stubs.stub_ripp import LEVELS
    declared = {level["level"] for level in _task()["harness"]["stub_levels"]}
    assert declared == set(LEVELS)
