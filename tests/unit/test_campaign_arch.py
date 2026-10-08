"""Domain architecture: the tool's output is not a domain list yet.

Four things here produce a wrong answer without producing an error, and each has
a test.

Pfam's nested families hit the same residues, so listing every hit as a domain
gives an architecture longer than the protein has -- the mistake the catalogued
design made by asking for a "non-overlapping" table. Ordering by the alignment
span or by the envelope can reorder the architecture. Pfam's own intra-clan rule
looks like the right resolution and is not available, because the pinned library
carries no clan annotation. And `AMP-binding` and `AMP-binding_C` are two halves
of one domain, so counting both double-counts every module.

The tool is fast enough here that the suite runs it where it needs it; the gold,
the oracle and the sweep all exercise it.
"""

from __future__ import annotations

import json
import pathlib
import shutil

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.templates.domain_architecture import analyse as dar
from npbench_c.templates.domain_architecture.build import (
    ARCHITECTURE_KEYS,
    COUNTERFACTUAL_KEYS,
    INVENTORY_KEYS,
    RECONCILIATION_KEYS,
    _assert_has_content,
)
from npbench_c.templates.domain_architecture.core import (
    ArchitectureError,
    Corpus,
    Hit,
    Rules,
    read_proteins,
)
from npbench_c.templates.domain_architecture.measure import measure
from npbench_c.templates.domain_architecture.solve import render
from npbench_c.tools.resources import tool_prefix

CAMPAIGNS = pathlib.Path(__file__).resolve().parents[2] / "campaigns"
CAMPAIGN = CAMPAIGNS / "arch-nrps-mibig-4_0-01"
PKS = CAMPAIGNS / "arch-pks-mibig-4_0-01"
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")


def _gold(campaign: pathlib.Path = CAMPAIGN) -> dict:
    return json.loads((campaign / "gold" / "gold.json").read_text())


def _task(campaign: pathlib.Path = CAMPAIGN) -> dict:
    return yaml.safe_load((campaign / "task.yaml").read_text())


def _rules(campaign: pathlib.Path = CAMPAIGN) -> Rules:
    return Rules.load(campaign)


@pytest.fixture(scope="module")
def corpus() -> Corpus:
    params = dar.load_params(CAMPAIGN)
    return Corpus.load(CAMPAIGN, params["annotations_input"],
                       params["proteins_input"])


# ------------------------------------------------------------- the input read


def test_the_shipped_fasta_parses_to_cluster_and_gene():
    params = dar.load_params(CAMPAIGN)
    proteins = read_proteins(CAMPAIGN / params["proteins_input"])
    assert proteins
    for bgc, gene in proteins:
        assert bgc.startswith("BGC") and gene
    assert all(seq and seq.isupper() for seq in proteins.values())


def test_a_header_that_is_not_cluster_pipe_gene_is_refused(tmp_path):
    """The FASTA is written by this template's own prepare_input, so a header
    that does not parse is a broken input rather than something to tolerate."""
    path = tmp_path / "bad.faa"
    path.write_text(">sp|P12345|SOMETHING\nMGPV\n")
    with pytest.raises(ArchitectureError, match="unparseable FASTA header"):
        read_proteins(path)


def test_the_shipped_fasta_covers_more_than_the_corpus(corpus):
    """Shipping only the selected genes would make the corpus size and every drop
    readable off the input instead of computed from the two resources, and would
    let the NRPS and PKS instantiations reveal each other's corpus by size."""
    census = _gold()["selection_census"]
    assert census["shipped_proteins"] > census["corpus_genes"]
    assert census["drop_other_module_type"] > 0


# ---------------------------------------------------- the curated record's shape


def test_the_curated_coordinates_are_unusable_and_the_campaign_says_so():
    """Nearly every typed domain block carries a placeholder location. A
    coordinate-level reconciliation would be computed over a hundredth of the
    record and reported as though it covered all of it."""
    usability = _gold()["curated_coordinate_usability"]
    assert usability["placeholder"] > 10 * usability.get("usable", 0)
    excluded = " ".join(str(b["field"]) for b in _task()["excluded_from_grading"])
    assert "curated domain coordinates" in excluded


def test_the_active_flag_is_unstated_on_many_modules():
    """So "active modules" is not a quantity this database supports. The census
    of the flag's real values is graded; nothing downstream filters on it."""
    flags = _gold()["curated_active_flag"]
    assert "None" in flags and flags["None"] > 0
    excluded = " ".join(str(b["field"]) for b in _task()["excluded_from_grading"])
    assert "active" in excluded


def test_a_module_naming_several_genes_is_excluded_and_counted():
    """Its per-gene counts cannot be apportioned without a rule nothing in the
    record supplies, so it is dropped -- and the drop is reported."""
    census = _gold()["selection_census"]
    assert census["modules_spanning_several_genes"] > 0
    assert census["genes_named_by_a_single_gene_module"] > census["corpus_genes"]


# ----------------------------------------------------------- the overlap problem


def test_the_raw_hits_overlap():
    """The number the catalogued "non-overlapping domain table" assumed was zero.
    If it were, the resolution policy would be ceremony and one of R4's axes a
    no-op; the build refuses gold in that case."""
    overlaps = _gold()["overlap_census"]
    assert overlaps["overlapping_pairs"] > 0
    assert overlaps["genes_with_an_overlap"] > 0
    assert overlaps["hits_total"] > _gold()["architecture_census"]["domains_total"]


def test_resolution_drops_hits_on_many_genes():
    table = _gold()["architecture_table"]
    dropped = [row for row in table if row["resolved_away"]]
    assert dropped
    for row in table:
        assert row["resolved_away"] == row["hits"] - row["domains"]
        # The string is empty exactly when there is no domain. Its FAMILY COUNT
        # is deliberately not asserted by splitting it: a Pfam name can contain
        # the separator, so the split over-counts -- which is the trap the test
        # below states, and which this assertion fell into on its first draft.
        assert bool(row["architecture"]) == (row["domains"] > 0)


def test_resolution_keeps_the_better_scoring_hit_of_a_nested_pair(corpus):
    """The rule, on a constructed pair: a short domain wholly inside a long one,
    the shorter scoring higher. Measuring the overlap against the LONGER span
    would never trigger, which is why the fraction is of the shorter."""
    inner = Hit(model="Inner", ali=(100, 140), env=(98, 142), score=200.0)
    outer = Hit(model="Outer", ali=(50, 400), env=(48, 402), score=100.0)
    kept = corpus.resolve([inner, outer])
    assert [h.model for h in kept] == ["Inner"]

    kept_all = corpus.resolve([inner, outer], policy="keep_all")
    assert [h.model for h in kept_all] == ["Outer", "Inner"]


def test_two_domains_that_merely_abut_are_both_kept(corpus):
    left = Hit(model="Left", ali=(10, 100), env=(8, 102), score=50.0)
    right = Hit(model="Right", ali=(101, 200), env=(99, 202), score=40.0)
    assert [h.model for h in corpus.resolve([left, right])] == ["Left", "Right"]


def test_an_undeclared_span_or_policy_is_refused(corpus):
    with pytest.raises(ArchitectureError, match="not a declared coordinate span"):
        corpus.resolve([], span="hmm")
    with pytest.raises(ArchitectureError, match="not a declared resolution policy"):
        corpus.resolve([], policy="keep_the_longest")


def test_pfam_clan_annotation_is_absent_and_the_rules_say_why():
    """Pfam's own intra-clan rule is the obvious resolution and it is not
    available: the pinned library carries no clan lines at all. Recording that is
    what keeps the campaign's rule from being misattributed to Pfam."""
    declared = json.loads(
        (CAMPAIGN / "reference" / "architecture_rules.json").read_text())
    note = declared["audit"]["resolution_policies_note"]
    assert "clan" in note and "not" in note.lower()
    assert "clan" in declared["external_resources"]["pfam_a_hmm"]["note"]


@image_present
def test_the_pinned_library_really_has_no_clan_lines():
    """The claim above, against the file. Asserted on a bounded prefix because
    the library is 1.5 GB and a clan line, if present, appears in every model's
    header."""
    library = _rules().resources["pfam_a_hmm"].resolve()
    with open(library, "r", errors="replace") as handle:
        head = handle.read(4 << 20)
    assert "\nNAME " in head, "sanity: this is an HMMER library"
    assert "\nCL   " not in head and "\nCL  " not in head


# -------------------------------------------------------------- the architecture


def test_the_architecture_is_rich_rather_than_a_domain_name():
    census = _gold()["architecture_census"]
    assert census["distinct_architectures"] >= 20
    assert census["domains_per_gene"]["max"] >= 10
    assert census["genes"] == len(_gold()["architecture_table"])


def test_the_corpus_is_multidomain():
    table = _gold()["architecture_table"]
    multidomain = sum(1 for row in table if row["domains"] >= 3)
    assert multidomain > len(table) // 2


def test_a_gene_with_no_domain_reports_nulls_rather_than_zeros():
    table = _gold()["architecture_table"]
    for row in table:
        if row["domains"]:
            assert row["first_domain_start"] is not None
            assert row["last_domain_end"] is not None
        else:
            assert row["first_domain_start"] is None
            assert row["last_domain_end"] is None
            assert row["architecture"] == ""


# ------------------------------------------------------------ the reconciliation


def test_the_reconciliation_disagrees_in_both_directions():
    """Perfect agreement would mean the curated record is recomputable from
    sequence and R3 grades nothing. The build refuses that, and refuses a single
    offset too -- a convention error rather than a distribution."""
    verdicts = _gold()["primary_verdicts"]
    assert verdicts["agrees"] > 0
    assert verdicts["computed_more"] > 0 or verdicts["computed_fewer"] > 0
    assert len(_gold()["primary_delta_census"]) >= 3


def test_the_delta_census_matches_the_table():
    gold = _gold()
    from collections import Counter
    counted = Counter(str(row["primary_delta"]) for row in gold["reconciliation_table"])
    assert {k: v for k, v in counted.items() if v} == gold["primary_delta_census"]
    for row in gold["reconciliation_table"]:
        assert row["primary_delta"] == row["computed_primary"] - row["curated_primary"]


def test_an_uncurated_domain_type_is_its_own_verdict():
    """A gene whose record names no condensation domain reports `not_curated`
    rather than a disagreement with zero -- there is nothing to disagree with."""
    gold = _gold()
    assert gold["secondary_verdicts"]["not_curated"] > 0
    for row in gold["reconciliation_table"]:
        if row["curated_secondary"] == 0:
            assert row["secondary_verdict"] == "not_curated"
        else:
            assert row["secondary_verdict"] != "not_curated"


def test_the_two_reconciled_domain_types_are_separate_claims():
    """The primary type is the one a module carries exactly once, so the curated
    column is reliable; the secondary one is recorded less consistently and its
    agreement is its own measurement. The two must not be the same number."""
    gold = _gold()
    a = gold["primary_verdicts"]
    c = gold["secondary_verdicts"]
    decided = c["agrees"] + c["computed_more"] + c["computed_fewer"]
    assert decided > 0
    assert a != c


def test_the_module_decomposition_is_a_different_claim_from_the_domain_count():
    """A gene can carry the right number of adenylation domains and not assemble
    into complete modules at all, so the two verdicts must not be the same
    column under two names."""
    gold = _gold()
    assert gold["module_verdicts"] != gold["primary_verdicts"]
    assert gold["complete_modules_total"] == sum(
        row["complete_modules"] for row in gold["module_table"])
    assert gold["module_verdicts"]["agrees"] > 0


def test_the_module_pattern_is_counted_without_overlapping_itself(corpus):
    """Three modules back to back are three, not five: the walk advances past a
    match rather than one family at a time.

    The pattern is fed as hits rather than as an architecture string on purpose
    -- see the next test.
    """
    pattern = list(corpus.rules.module_pattern)
    key = corpus.keys[0]
    hits = {key: [Hit(model=model, ali=(100 * i + 1, 100 * i + 80),
                      env=(100 * i, 100 * i + 81), score=100.0 - i)
                  for i, model in enumerate(pattern * 3)]}
    result = corpus.module_decomposition(hits)
    row = next(r for r in result["module_table"]
               if (r["bgc"], r["gene"]) == key)
    assert row["complete_modules"] == 3


def test_a_pfam_family_name_can_contain_the_architecture_separator(corpus):
    """`AMP-binding` holds a hyphen, so splitting the hyphen-joined architecture
    shreds every family name and the declared pattern matches nothing. The
    decomposition reads the resolved hits instead; the string is for reading.

    This is not hypothetical -- it is how the module decomposition was first
    written, and the build's own refusal ("the declared module pattern is found
    on no gene") is what caught it.
    """
    assert any("-" in family for family in corpus.rules.module_pattern)
    architecture = "-".join(corpus.rules.module_pattern)
    assert architecture.split("-") != list(corpus.rules.module_pattern)
    assert _gold()["complete_modules_total"] > 0
    # And on real gold: splitting over-counts, so a reader that recovers the
    # architecture that way reports more domains than the gene has.
    row = next(r for r in _gold()["architecture_table"] if r["domains"] > 3)
    assert len(row["architecture"].split("-")) > row["domains"]


def test_the_strata_partition_the_corpus_and_behave_differently():
    """A stratification every stratum answers the same way is decoration, and the
    build refuses it. Where a family declares one stratum the census has one
    populated name and the check does not apply."""
    gold = _gold()
    strata = gold["strata_census"]
    assert sum(strata.values()) == gold["architecture_census"]["genes"]
    by_stratum = gold["secondary_verdicts_by_stratum"]
    assert set(strata) <= set(by_stratum)
    if len(strata) > 1:
        shapes = {json.dumps(by_stratum[name], sort_keys=True) for name in strata}
        assert len(shapes) > 1, by_stratum


def test_a_stratum_whose_type_means_absence_grades_the_zero(corpus):
    """Where the module type means the secondary domain is absent -- a trans-AT
    module has no acyltransferase -- a curated zero is a CLAIM, so the row gets a
    real verdict rather than `not_curated`. Which strata those are is read from
    the per-module-type census, not asserted."""
    gold = _gold()
    per_type = gold["domains_by_module_type"]
    secondary = gold["secondary_domain"]
    asserting = [name for name, types in corpus.rules.strata.items()
                 if sum(per_type.get(t, {}).get(secondary, 0) for t in types) == 0]
    if not asserting:
        pytest.skip("this family declares no stratum that asserts absence")
    for name in asserting:
        block = gold["secondary_verdicts_by_stratum"][name]
        assert block["not_curated"] == 0, name
        assert block["agrees"] > 0, name
    rows = {r["stratum"] for r in gold["reconciliation_table"]
            if r["secondary_verdict"] == "not_curated"}
    assert not (rows & set(asserting))


# ------------------------------------------------------------ the counterfactual


def test_each_declared_reading_is_a_cell_and_some_move():
    gold = _gold()
    rules = _rules()
    expected = {f"{c}|{p}|{s}" for c in rules.cutoffs
                for p in rules.resolution_policies for s in rules.coordinate_spans}
    assert set(gold["policy_sweep"]) == expected
    moved = [cell["genes_differing_from_baseline"]
             for cell in gold["policy_sweep"].values()]
    assert max(moved) > 0
    baseline = f"{rules.baseline_cutoff}|{rules.baseline_resolution}|{rules.baseline_span}"
    assert gold["policy_sweep"][baseline]["genes_differing_from_baseline"] == 0


def test_not_resolving_overlaps_is_the_axis_that_moves_most():
    """The measured shape of R4's first axis: the resolution policy changes
    hundreds of architectures, while the coordinate span changes almost none."""
    gold = _gold()
    rules = _rules()
    base = rules.baseline_cutoff
    unresolved = gold["policy_sweep"][f"{base}|keep_all|{rules.baseline_span}"]
    other_span = gold["policy_sweep"][
        f"{base}|{rules.baseline_resolution}|"
        f"{[s for s in rules.coordinate_spans if s != rules.baseline_span][0]}"]
    assert unresolved["genes_differing_from_baseline"] > 100
    assert (unresolved["genes_differing_from_baseline"]
            > other_span["genes_differing_from_baseline"])
    assert unresolved["domains_total"] > gold["architecture_census"]["domains_total"]


def test_the_curated_cutoffs_are_interchangeable_and_the_evalue_is_not():
    """The contrast the E-value is in the grid for: a per-model curated threshold
    and a global E-value are not substitutes, and the familiar flag is the one
    that moves the answer."""
    gold = _gold()
    rules = _rules()
    tail = f"|{rules.baseline_resolution}|{rules.baseline_span}"
    curated = [gold["policy_sweep"][c + tail]["genes_differing_from_baseline"]
               for c in rules.cutoffs if c.startswith("--cut")]
    evalue = [gold["policy_sweep"][c + tail]["genes_differing_from_baseline"]
              for c in rules.cutoffs if not c.startswith("--cut")]
    assert evalue, "the grid declares no E-value cutoff"
    assert max(curated) < 50
    assert max(evalue) > 100


def test_the_family_to_domain_mapping_is_load_bearing():
    """Pfam splits the adenylation domain into two families, so counting both
    double-counts every module. The build refuses gold unless the best and worst
    declared mappings differ by at least a factor of two."""
    gold = _gold()
    agreements = {name: block["verdicts"]["agrees"]
                  for name, block in gold["mapping_sweep"].items()}
    assert len(agreements) >= 3
    assert len(set(agreements.values())) == len(agreements)
    assert max(agreements.values()) >= 2 * min(agreements.values())
    assert agreements[_rules().baseline_mapping] == max(agreements.values())


def test_the_mapping_sweep_leaves_the_secondary_mapping_out():
    """It is the mapping for the OTHER reconciled domain type, not an alternative
    reading of the primary one, so sweeping it would mix two axes. Its key is the
    secondary domain's own name, and the rules loader checks that -- when the two
    drifted apart during this template's generalisation, the secondary
    reconciliation silently computed zero for every gene."""
    rules = _rules()
    assert rules.secondary_domain in rules.domain_mappings
    assert rules.secondary_domain not in _gold()["mapping_sweep"]
    assert rules.baseline_mapping != rules.secondary_domain


def test_the_overlap_fraction_is_declared_and_not_swept():
    rules = _rules()
    assert 0.0 < rules.overlap_fraction < 1.0
    excluded = " ".join(str(b["field"]) for b in _task()["excluded_from_grading"])
    assert "overlap fraction" in excluded


def test_the_deleted_domain_counterfactual_is_refused_with_its_measurement():
    """It was built and measured before being dropped: the outcome is one answer
    repeated. The exclusion records that rather than omitting the rung silently."""
    blocks = {str(b["field"]): b["reason"] for b in _task()["excluded_from_grading"]}
    key = next(k for k in blocks if "deleted" in k)
    assert "measured" in blocks[key] and "1,112" in blocks[key]


# ---------------------------------------------------------------- the ladder


def test_the_oracle_clears_every_rung(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "architecture_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(),
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["score"] == 1.0
    assert result["depth"] == 4
    assert result["monotonic"] is True
    # R1's chance level is nominal at 1.0, so the floor that has to be low is the
    # one computed over R2 upward.
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
                             ("r3", ARCHITECTURE_KEYS, "architecture.json"),
                             ("r4", RECONCILIATION_KEYS, "reconciliation.json")):
        payload = json.loads(
            (CAMPAIGN / "gold" / "warm" / rung / name).read_text())
        assert set(payload) == set(keys)
        assert all(payload[k] == gold[k] for k in keys)


def test_the_rung_key_sets_do_not_overlap():
    sets = [set(INVENTORY_KEYS), set(ARCHITECTURE_KEYS),
            set(RECONCILIATION_KEYS), set(COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right


def test_the_echoes_of_the_rules_file_are_reported_and_not_graded():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    for echo in ("mibig_release", "module_family", "baseline_cutoff",
                 "baseline_resolution", "baseline_span", "baseline_mapping",
                 "panel_size"):
        assert echo not in graded
        assert echo in _gold()


# --------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    copy = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, copy)
    shutil.rmtree(copy / "gold")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "architecture_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    assert measure(copy, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_projects_rows_without_reordering_them(tmp_path):
    gold = _gold()
    report = render(gold)
    report["architecture_table"] = list(reversed(report["architecture_table"]))
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "architecture_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), gold,
                   json.loads((CAMPAIGN / "grading.json").read_text()))
    assert result["depth"] == 1


def test_a_null_coordinate_survives_the_projection_when_gold_has_one(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    report = render(_gold())
    report["architecture_table"] = [
        {**report["architecture_table"][0], "domains": 0, "architecture": "",
         "first_domain_start": None, "last_domain_end": None}]
    (submission / "architecture_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    measured = measure(CAMPAIGN, submission)
    assert measured["architecture_table"][0]["first_domain_start"] is None


# ------------------------------------------------------------ the gold refusals


def test_the_build_refuses_a_corpus_too_small_to_audit():
    gold = _gold()
    with pytest.raises(ArchitectureError, match="anecdote"):
        _assert_has_content({**gold, "selection_census": {
            **gold["selection_census"], "corpus_genes": 10}})


def test_the_build_refuses_an_unexercised_selection_clause():
    gold = _gold()
    for field in ("drop_excluded_by_status", "drop_other_module_type"):
        with pytest.raises(ArchitectureError, match="not exercised"):
            _assert_has_content({**gold, "selection_census": {
                **gold["selection_census"], field: 0}})


def test_the_build_refuses_an_architecture_with_nothing_in_it():
    gold = _gold()
    with pytest.raises(ArchitectureError, match="one protein repeated"):
        _assert_has_content({**gold, "architecture_census": {
            **gold["architecture_census"], "distinct_architectures": 3}})
    with pytest.raises(ArchitectureError, match="multidomain"):
        _assert_has_content({**gold, "architecture_census": {
            **gold["architecture_census"],
            "domains_per_gene": {"min": 1, "median": 1.0, "max": 2}}})


def test_the_build_refuses_a_corpus_with_no_overlaps():
    gold = _gold()
    with pytest.raises(ArchitectureError, match="central complication is absent"):
        _assert_has_content({**gold, "overlap_census": {
            **gold["overlap_census"], "overlapping_pairs": 0}})
    with pytest.raises(ArchitectureError, match="same answer"):
        _assert_has_content({**gold, "architecture_table": [
            {**row, "resolved_away": 0} for row in gold["architecture_table"]]})


def test_the_build_refuses_a_reconciliation_with_nothing_to_grade():
    gold = _gold()
    verdicts = gold["primary_verdicts"]
    with pytest.raises(ArchitectureError, match="broken mapping"):
        _assert_has_content({**gold, "primary_verdicts":
                             {**verdicts, "agrees": 0}})
    with pytest.raises(ArchitectureError, match="recomputable from sequence"):
        _assert_has_content({**gold, "primary_verdicts": {
            **verdicts, "computed_more": 0, "computed_fewer": 0}})
    with pytest.raises(ArchitectureError, match="single offset"):
        _assert_has_content({**gold, "primary_delta_census": {"0": 1000, "1": 50}})


def test_the_build_refuses_a_secondary_reconciliation_of_one_kind_of_row():
    """A corpus in a single stratum where no gene leaves the secondary domain
    uncurated has only one kind of secondary row. The check counts POPULATED
    strata, not the by-stratum census's keys: that census always carries a
    "mixed" row, so counting its keys made this pass vacuously."""
    gold = _gold()
    single = {**gold,
              "secondary_verdicts": {**gold["secondary_verdicts"],
                                     "not_curated": 0},
              "strata_census": {"only": gold["architecture_census"]["genes"]}}
    with pytest.raises(ArchitectureError, match="single stratum"):
        _assert_has_content(single)


def test_the_build_refuses_a_cosmetic_counterfactual():
    gold = _gold()
    flat = {cell: {**block, "genes_differing_from_baseline": 0}
            for cell, block in gold["policy_sweep"].items()}
    with pytest.raises(ArchitectureError, match="no axis"):
        _assert_has_content({**gold, "policy_sweep": flat})
    same = {cell: {"domains_total": 1, "distinct_architectures": 1,
                   "a_domain_agreement": 1, "genes_differing_from_baseline": 5}
            for cell in gold["policy_sweep"]}
    with pytest.raises(ArchitectureError, match="not independent"):
        _assert_has_content({**gold, "policy_sweep": same})


def test_the_build_refuses_a_mapping_choice_that_does_not_matter():
    gold = _gold()
    sweep = gold["mapping_sweep"]
    tied = {name: {**block, "verdicts": {**block["verdicts"], "agrees": 100}}
            for name, block in sweep.items()}
    with pytest.raises(ArchitectureError, match="two names"):
        _assert_has_content({**gold, "mapping_sweep": tied})
    close = {}
    for i, (name, block) in enumerate(sorted(sweep.items())):
        close[name] = {**block, "verdicts": {**block["verdicts"],
                                             "agrees": 100 + i}}
    with pytest.raises(ArchitectureError, match="factor of two"):
        _assert_has_content({**gold, "mapping_sweep": close})


# --------------------------------------------------------------- the stub probe


def test_the_task_declares_the_stub_levels_the_module_implements():
    from npbench_c.sweep.stubs.stub_arch import LEVELS
    declared = {level["level"] for level in _task()["harness"]["stub_levels"]}
    assert declared == set(LEVELS)


# ------------------------------------------------- the other instantiation (PKS)

pks_built = pytest.mark.skipif(not (PKS / "gold" / "gold.json").is_file(),
                               reason="the PKS instantiation is not built")


@pytest.fixture(scope="module")
def pks_corpus() -> Corpus:
    params = dar.load_params(PKS)
    return Corpus.load(PKS, params["annotations_input"],
                       params["proteins_input"])


@pks_built
def test_the_two_instantiations_share_a_template_and_differ_only_by_rules():
    """The instantiation axis is a declared rules file, not a forked template.
    Both campaigns name the same commands and the same stub module; what differs
    is which constants file the build copies in."""
    nrps, pks = _task(CAMPAIGN), _task(PKS)
    assert nrps["template_id"] == pks["template_id"] == "domain-architecture"
    assert nrps["harness"] == pks["harness"]
    assert (nrps["template_params"]["rules_source"]
            != pks["template_params"]["rules_source"])
    # and the pinned table the agent reads is at the same path in both
    assert nrps["constraints"]["pinned_tables"] == pks["constraints"]["pinned_tables"]


@pks_built
def test_the_pks_family_accepts_several_module_types():
    """A modular PKS gene normally carries a loading module of its own type
    alongside its elongation modules, so the family is a set rather than one
    type -- and the strata, not the selection rule, keep the architectures
    apart."""
    rules = _rules(PKS)
    assert len(rules.accepted_module_types) > 1
    assert len(_rules(CAMPAIGN).accepted_module_types) == 1
    covered = {t for types in rules.strata.values() for t in types}
    assert covered == set(rules.accepted_module_types)


@pks_built
def test_the_trans_at_stratum_finds_no_acyltransferase():
    """This instantiation's own control, and the reason it is not a relabelling.
    A trans-AT module has no acyltransferase -- the enzyme acts from a separate
    protein -- so the curated count is zero by definition and the computed
    architecture should find none. Measured, not assumed: if the tool found some,
    the finding would be about the tool."""
    gold = _gold(PKS)
    by_stratum = gold["secondary_verdicts_by_stratum"]
    assert "trans_at" in by_stratum
    trans = by_stratum["trans_at"]
    assert trans["agrees"] > 0
    assert trans["computed_more"] == 0 and trans["computed_fewer"] == 0
    # and the zero is graded as a claim, not reported as a silence
    assert trans["not_curated"] == 0
    rows = [r for r in gold["reconciliation_table"] if r["stratum"] == "trans_at"]
    assert rows and all(r["computed_secondary"] == 0 for r in rows)
    assert all(r["curated_secondary"] == 0 for r in rows)


@pks_built
def test_the_cis_at_stratum_does_find_one_and_disagrees_sometimes():
    gold = _gold(PKS)
    cis = gold["secondary_verdicts_by_stratum"]["cis_at"]
    assert cis["agrees"] > 0
    assert cis["computed_more"] + cis["computed_fewer"] > 0
    assert cis != gold["secondary_verdicts_by_stratum"]["trans_at"]


@pks_built
def test_a_gene_spanning_both_strata_is_mixed():
    gold = _gold(PKS)
    assert gold["strata_census"].get("mixed", 0) > 0
    assert "mixed" not in _rules(PKS).strata


@pks_built
def test_the_ketosynthase_is_split_across_three_pfam_families():
    """The mapping axis, harder here than for the adenylation domain: counting
    all three families triples every module, so agreement collapses."""
    sweep = _gold(PKS)["mapping_sweep"]
    agreements = {name: b["verdicts"]["agrees"] for name, b in sweep.items()}
    assert len(agreements) >= 4
    best = max(agreements.values())
    worst = min(agreements.values())
    assert best >= 10 * worst, agreements
    widest = max(sweep, key=lambda n: len(sweep[n]["families"]))
    assert len(sweep[widest]["families"]) == 3
    assert agreements[widest] == worst


@pks_built
def test_the_pks_module_pattern_cannot_match_a_trans_at_gene():
    """The cis-AT module core includes the acyltransferase, which a trans-AT gene
    does not encode -- so the decomposition separates the strata too, through a
    second and independent measure."""
    gold = _gold(PKS)
    pattern = _rules(PKS).module_pattern
    assert "Acyl_transf_1" in pattern
    strata = {(r["bgc"], r["gene"]): r["stratum"]
              for r in gold["reconciliation_table"]}
    trans = [r for r in gold["module_table"]
             if strata[(r["bgc"], r["gene"])] == "trans_at"]
    assert trans and all(r["complete_modules"] == 0 for r in trans)
    cis = [r for r in gold["module_table"]
           if strata[(r["bgc"], r["gene"])] == "cis_at"]
    assert any(r["complete_modules"] > 0 for r in cis)


@pks_built
def test_the_pks_oracle_clears_every_rung(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "architecture_report.json").write_text(
        json.dumps(render(_gold(PKS)), indent=2, sort_keys=True) + "\n")
    result = grade(measure(PKS, submission), _gold(PKS),
                   json.loads((PKS / "grading.json").read_text()))
    assert result["score"] == 1.0
    assert result["depth"] == 4
    assert result["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


@pks_built
def test_the_pks_measurer_never_reads_gold(tmp_path):
    copy = tmp_path / "campaign"
    shutil.copytree(PKS, copy)
    shutil.rmtree(copy / "gold")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "architecture_report.json").write_text(
        json.dumps(render(_gold(PKS)), indent=2, sort_keys=True) + "\n")
    assert measure(copy, submission) == measure(PKS, submission)


@pks_built
def test_the_pks_overlaps_are_heavier_than_the_nrps_ones():
    """Worth stating rather than assuming: the reductase families this corpus
    carries overlap far more than the NRPS panel's do, so the resolution policy
    matters on most of the corpus instead of a quarter of it."""
    pks = _gold(PKS)["overlap_census"]
    nrps = _gold(CAMPAIGN)["overlap_census"]
    pks_rate = pks["genes_with_an_overlap"] / _gold(PKS)["architecture_census"]["genes"]
    nrps_rate = (nrps["genes_with_an_overlap"]
                 / _gold(CAMPAIGN)["architecture_census"]["genes"])
    assert pks_rate > nrps_rate
    assert pks_rate > 0.5
