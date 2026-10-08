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

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "arch-nrps-mibig-4_0-01")
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
    verdicts = _gold()["a_domain_verdicts"]
    assert verdicts["agrees"] > 0
    assert verdicts["computed_more"] > 0 or verdicts["computed_fewer"] > 0
    assert len(_gold()["a_delta_census"]) >= 3


def test_the_delta_census_matches_the_table():
    gold = _gold()
    from collections import Counter
    counted = Counter(str(row["a_delta"]) for row in gold["reconciliation_table"])
    assert {k: v for k, v in counted.items() if v} == gold["a_delta_census"]
    for row in gold["reconciliation_table"]:
        assert row["a_delta"] == row["computed_a_domains"] - row["curated_a_domains"]


def test_an_uncurated_domain_type_is_its_own_verdict():
    """A gene whose record names no condensation domain reports `not_curated`
    rather than a disagreement with zero -- there is nothing to disagree with."""
    gold = _gold()
    assert gold["c_domain_verdicts"]["not_curated"] > 0
    for row in gold["reconciliation_table"]:
        if row["curated_c_domains"] == 0:
            assert row["c_verdict"] == "not_curated"
        else:
            assert row["c_verdict"] != "not_curated"


def test_the_condensation_reconciliation_is_worse_than_the_adenylation_one():
    """The measured contrast: a curator records one adenylation domain per module
    by construction, and records a condensation domain only sometimes -- and when
    they do, agreement with the tool is far lower."""
    gold = _gold()
    a = gold["a_domain_verdicts"]
    c = gold["c_domain_verdicts"]
    a_rate = a["agrees"] / (a["agrees"] + a["computed_more"] + a["computed_fewer"])
    curated_c = c["agrees"] + c["computed_more"] + c["computed_fewer"]
    c_rate = c["agrees"] / curated_c
    assert curated_c > 0
    assert c_rate < a_rate


def test_the_module_decomposition_is_a_different_claim_from_the_domain_count():
    """A gene can carry the right number of adenylation domains and not assemble
    into complete modules at all, so the two verdicts must not be the same
    column under two names."""
    gold = _gold()
    assert gold["module_verdicts"] != gold["a_domain_verdicts"]
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


def test_the_mapping_sweep_leaves_the_condensation_mapping_out():
    """It is the mapping for the OTHER reconciled domain type, not an alternative
    reading of the adenylation domain, so sweeping it would mix two axes."""
    assert "condensation_domain" in _rules().domain_mappings
    assert "condensation_domain" not in _gold()["mapping_sweep"]


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
    for echo in ("mibig_release", "module_type", "baseline_cutoff",
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
    verdicts = gold["a_domain_verdicts"]
    with pytest.raises(ArchitectureError, match="broken mapping"):
        _assert_has_content({**gold, "a_domain_verdicts":
                             {**verdicts, "agrees": 0}})
    with pytest.raises(ArchitectureError, match="recomputable from sequence"):
        _assert_has_content({**gold, "a_domain_verdicts": {
            **verdicts, "computed_more": 0, "computed_fewer": 0}})
    with pytest.raises(ArchitectureError, match="single offset"):
        _assert_has_content({**gold, "a_delta_census": {"0": 1000, "1": 50}})
    with pytest.raises(ArchitectureError, match="not_curated"):
        _assert_has_content({**gold, "c_domain_verdicts": {
            **gold["c_domain_verdicts"], "not_curated": 0}})


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
