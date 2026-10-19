"""Remote homology: grading the search, never the homology.

The question the catalogue asked -- is this pair homologous -- has no oracle at
the identities this campaign is about. So the campaign grades what a pinned
pairwise search REPORTS against what InterPro itself asserts, band by band, and
excludes the catalogued claim by name.

Three things here produce a wrong answer without producing an error.

The query is a member of the corpus it is searched against, so every query
returns itself at 100 per cent identity. That self-hit is the query, not a
finding, and keeping it inflates every count by exactly one in a way that stays
internally consistent everywhere. That is the competence probe.

The label has three values and the third is not "non-homologous". A pair sharing
no InterPro entry is unclassified: the resource records what has been classified,
which is not the same as what is related. Reading the third value as a negative
turns a census into a claim the data does not support.

And a loss rate without its denominator is not a measurement. The identity bands
differ in size by two orders of magnitude on some queries, so the ladder reports
the rate and the baseline band sizes together, and the absolute counts and the
rates tell different stories on purpose.

The search is slow enough -- one index build per run -- that the suite leans on
gold where the answer is already recorded and runs the tool only where it must.
"""

from __future__ import annotations

import copy
import json
import pathlib
import shutil

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.templates.remote_homology import analyse as rha
from npbench_c.templates.remote_homology.build import (
    COUNTERFACTUAL_KEYS,
    INVENTORY_KEYS,
    RECOVERY_KEYS,
    VERDICT_KEYS,
    _assert_has_content,
)
from npbench_c.templates.remote_homology.core import (
    Corpus,
    Hit,
    HomologyError,
    Rules,
    _ordered,
    read_corpus,
)
from npbench_c.templates.remote_homology.measure import measure
from npbench_c.templates.remote_homology.solve import render
from npbench_c.tools.resources import tool_prefix

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "homology-twilight-ipr036188-01")
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
def corpus() -> Corpus:
    params = rha.load_params(CAMPAIGN)
    return Corpus.load(CAMPAIGN, params["corpus_input"], params["labels_input"])


# ------------------------------------------------------------ what is shipped


def test_the_search_database_is_not_shipped():
    """It follows from the declared invocation over the shipped corpus. Shipping
    it would hand over the search's own setup, and an index is also where a
    different --max-seqs would hide."""
    shipped = sorted(p.name for p in (CAMPAIGN / "inputs").iterdir())
    assert shipped == ["corpus_ipr036188_reviewed_2026_03.fasta.gz",
                       "interpro_labels_110_0.json.gz", "provenance.json"]
    assert "not shipped" in _task()["license_notice"]["redistribution"]


def test_the_label_table_is_shipped_and_says_why():
    """The campaign cannot be posed without it, and it is not the answer: the
    answer is what the search reports, which the table does not contain."""
    provenance = json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())
    assert "cannot be posed without it" in provenance["labels"]["why_it_is_shipped"]
    assert provenance["labels"]["license"].startswith("CC0")
    assert provenance["corpus"]["license"].startswith("CC BY")


def test_the_provenance_does_not_carry_the_graded_values():
    """`corpus_sequences` and `corpus_sha256` are graded at R1. The provenance
    carries a digest of the FILE, which is a different number by construction."""
    provenance = json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())
    gold = _gold()
    flat = json.dumps(provenance)
    assert gold["corpus_sha256"] not in flat
    assert str(gold["corpus_sequences"]) not in flat


def test_no_pfam_release_is_involved():
    """The image carries Pfam 35.0 for three other campaigns. This campaign's
    label is InterPro 110.0, a different release of a different resource, and
    conflating them would make the label unreproducible."""
    assert _task()["constraints"]["external_resources"] == []
    assert "Pfam" in _task()["constraints"]["external_resources_note"]
    assert _gold()["label_resource"] == "InterPro 110.0"


# -------------------------------------------------------------- the rules file


def _rules_copy(tmp_path, **overrides):
    copy_dir = tmp_path / "campaign"
    (copy_dir / "reference").mkdir(parents=True, exist_ok=True)
    spec = json.loads((CAMPAIGN / "reference" / "homology_rules.json").read_text())
    spec["homology"].update(overrides)
    (copy_dir / "reference" / "homology_rules.json").write_text(json.dumps(spec))
    return copy_dir


def test_a_band_partition_that_does_not_match_its_labels_is_refused(tmp_path):
    with pytest.raises(HomologyError, match="the same partition"):
        Rules.load(_rules_copy(tmp_path, band_labels=["a", "b"]))


def test_a_ladder_without_a_zero_window_is_refused(tmp_path):
    """Without the identity control every loss is uncalibrated."""
    with pytest.raises(HomologyError, match="no zero-length window"):
        Rules.load(_rules_copy(tmp_path, shuffle_windows=[60, 120, 240]))


def test_a_ladder_out_of_order_is_refused(tmp_path):
    with pytest.raises(HomologyError, match="not in increasing order"):
        Rules.load(_rules_copy(tmp_path, shuffle_windows=[0, 240, 120]))


def test_repeated_windows_are_refused(tmp_path):
    with pytest.raises(HomologyError, match="windows repeat"):
        Rules.load(_rules_copy(tmp_path, shuffle_windows=[0, 60, 60]))


def test_an_entry_type_that_is_both_family_and_superfamily_is_refused(tmp_path):
    """One entry would make a pair same-family and superfamily-only at once, so
    the two labels would stop being exclusive."""
    with pytest.raises(HomologyError, match="declared as both"):
        Rules.load(_rules_copy(tmp_path, superfamily_types=["family"]))


def test_too_few_queries_is_refused(tmp_path):
    with pytest.raises(HomologyError, match="cannot carry a per-query census"):
        Rules.load(_rules_copy(tmp_path, queries=2))


def test_an_empty_label_type_list_is_refused(tmp_path):
    with pytest.raises(HomologyError, match="neither may be empty"):
        Rules.load(_rules_copy(tmp_path, family_types=[]))


# ------------------------------------------------------------- the corpus read


def test_the_corpus_is_the_declared_superfamily(corpus):
    gold = _gold()
    assert corpus.size == gold["corpus_sequences"]
    assert corpus.sequences_sha256() == gold["corpus_sha256"]
    assert all(corpus.rules.corpus_superfamily in s
               for s in corpus.superfamilies.values())


def test_the_digest_is_over_the_sequences_and_not_the_file(corpus, tmp_path):
    """A corpus rewrapped at another line width must fingerprint the same, which
    is what makes the digest a parse check rather than a checksum."""
    path = tmp_path / "rewrapped.faa"
    path.write_text("".join(
        f">{a}\n" + "\n".join(s[i:i + 23] for i in range(0, len(s), 23)) + "\n"
        for a, s in sorted(corpus.sequences.items())))
    assert read_corpus(path) == dict(corpus.sequences)


def test_a_duplicate_accession_is_refused(tmp_path):
    path = tmp_path / "dup.faa"
    path.write_text(">A\nMKV\n>A\nMKL\n")
    with pytest.raises(HomologyError, match="appears twice"):
        read_corpus(path)


def test_a_sequence_outside_the_declared_alphabet_is_refused(tmp_path):
    path = tmp_path / "odd.faa"
    path.write_text(">A\nMKV*\n")
    with pytest.raises(HomologyError, match="outside the\n?.*alphabet"):
        read_corpus(path)


def test_an_empty_corpus_is_refused(tmp_path):
    path = tmp_path / "empty.faa"
    path.write_text("\n")
    with pytest.raises(HomologyError, match="holds no sequence"):
        read_corpus(path)


# --------------------------------------------------------------- the query rule


def test_no_two_families_share_a_query():
    """One accession can be the first member of two families, and auditing it
    twice would make a family stand for nothing."""
    accessions = [row["accession"] for row in _gold()["query_table"].values()]
    assert len(set(accessions)) == len(accessions)


def test_the_query_rule_is_reproducible_from_the_shipped_files(corpus):
    """The measurer recomputes the query set rather than reading it out of a
    submission, so the rule has to give the same answer every time."""
    first = corpus.queries()
    assert first == corpus.queries()
    gold = _gold()
    assert {entry: row["accession"] for entry, row in gold["query_table"].items()} \
        == first


def test_the_declared_queries_are_the_largest_families(corpus):
    sizes = {entry: n for entry, n in corpus.family_census().items()}
    chosen = set(corpus.queries())
    smallest_chosen = min(sizes[entry] for entry in chosen)
    bigger_unchosen = [e for e, n in sizes.items()
                       if e not in chosen and n > smallest_chosen]
    # A bigger family may be skipped only because every member it could offer was
    # already taken by a bigger one.
    for entry in bigger_unchosen:
        members = {a for a, es in corpus.families.items() if entry in es}
        assert members <= set(corpus.queries().values()), entry


# ------------------------------------------------------------- the label join


def test_the_third_label_is_not_a_claim_of_non_homology(corpus):
    """InterPro records what has been classified, which is not the same as what
    is related, so the vocabulary has `unclassified_pair` and no rung treats it
    as a negative."""
    labels = {corpus.pair_label(a, b) for a, b in
              list(zip(sorted(corpus.sequences)[:40], sorted(corpus.sequences)[1:41]))}
    assert labels <= {"same_family", "superfamily_only", "unclassified_pair"}
    excluded = " ".join(str(b["field"]) for b in _task()["excluded_from_grading"])
    assert "non-homologous" in excluded


def test_every_corpus_pair_shares_the_declared_superfamily(corpus):
    """The corpus was selected by that superfamily, so `unclassified_pair` cannot
    occur between two corpus members -- which is why the campaign's question is
    never whether two members are related at all."""
    accessions = sorted(corpus.sequences)[:60]
    labels = {corpus.pair_label(a, b)
              for i, a in enumerate(accessions) for b in accessions[i + 1:]}
    assert "unclassified_pair" not in labels


def test_the_band_boundaries_are_half_open_with_a_sentinel(corpus):
    rules = corpus.rules
    assert rules.band(0.0) == rules.band_labels[0]
    assert rules.band(0.2499) == rules.band_labels[0]
    assert rules.band(0.25) == rules.band_labels[1]
    assert rules.band(0.3499) == rules.band_labels[1]
    assert rules.band(0.35) == rules.band_labels[-1]
    assert rules.band(1.0) == rules.band_labels[-1]


# ----------------------------------------------------------------- the recovery


def test_the_search_misses_relationships_the_resource_asserts():
    """The campaign's subject. If the search recovered everything, the gap it
    exists to measure would not occur and the build refuses."""
    missed = {a: row["same_family_missed"] for a, row in _gold()["recovery"].items()}
    assert any(missed.values()), missed
    rates = [row["recovery_rate"] for row in _gold()["recovery"].values()]
    assert all(r is not None and 0.0 < r <= 1.0 for r in rates), rates


def test_the_declared_family_size_is_not_the_asserted_relative_count():
    """Two numbers that look like the same number and are not.

    A query stands for one declared family, and `family_members` is that
    family's size. But a protein commonly carries several InterPro families, and
    `same_family_asserted` counts every corpus member sharing ANY family with the
    query. Measured here they differ on every declared query and in both
    directions, so reading one for the other would make the recovery rate look
    wrong -- and "fixing" either to match the other would be the bug.
    """
    gold = _gold()
    pairs = [(row["family_members"], gold["recovery"][row["accession"]]["same_family_asserted"])
             for row in gold["query_table"].values()]
    assert all(declared != asserted for declared, asserted in pairs), pairs
    assert any(declared > asserted for declared, asserted in pairs), pairs
    assert any(declared < asserted for declared, asserted in pairs), pairs


def test_the_asserted_sets_may_coincide_across_queries(corpus):
    """InterPro families are not a partition, so two declared queries can carry a
    family in common and then assert the same relative set. The recovery rows are
    still per query -- different proteins, different hit lists -- but they are not
    independent samples, and a reading that treats them as such would be wrong.
    """
    gold = _gold()
    asserted = {row["accession"]: gold["recovery"][row["accession"]]["same_family_asserted"]
                for row in gold["query_table"].values()}
    hits = {accession: gold["baseline_hits"][accession] for accession in asserted}
    # Measured on this corpus: at least two queries share an asserted count,
    # while no two share a hit count.
    assert len(set(asserted.values())) < len(asserted), asserted
    assert len(set(hits.values())) == len(hits), hits


def test_the_recovery_arithmetic_closes():
    for accession, row in _gold()["recovery"].items():
        assert row["same_family_detected"] + row["same_family_missed"] \
            == row["same_family_asserted"], accession
        assert row["same_family_detected"] + row["superfamily_only_detected"] \
            == row["hits"], accession


def test_the_band_census_counts_every_band_and_label_always():
    """A census whose keys depend on the data would make the measurer's shape
    depend on the answer.

    A claim about the key SET, not their order: gold.json is canonically
    key-sorted like every other artifact here, so the file's order is
    lexicographic and says nothing about the declared partition.
    """
    gold, rules = _gold(), _rules()
    for accession, cells in gold["band_census"].items():
        assert set(cells) == set(rules.band_labels), accession
        for band, cell in cells.items():
            assert set(cell) == {"same_family", "superfamily_only",
                                 "unclassified_pair"}, (accession, band)


def test_the_bands_are_read_in_declared_order_and_not_alphabetically(corpus):
    """Where order matters it is the declared order, and this is a case where the
    two differ: sorted alphabetically the bands read 25_to_35, over_35, under_25,
    which would make the loss-rate ordering test compare the wrong pairs.

    A claim about the producing functions, which is where order is load-bearing.
    """
    rules = corpus.rules
    assert list(rules.band_labels) != sorted(rules.band_labels)
    hits = (Hit("Q", "a", 0.10, 100, 1e-9, 90.0),
            Hit("Q", "b", 0.30, 100, 1e-9, 90.0),
            Hit("Q", "c", 0.90, 100, 1e-9, 90.0))
    census = corpus.band_census(hits, "Q")
    assert list(census) == list(rules.band_labels)
    # And the ordering test reads the rates distant-first, so a run that fell
    # only in the weakest band must not read as ordered against a stronger one.
    assert _ordered([0.9, 0.5, 0.1]) is True
    assert _ordered([0.1, 0.5, 0.9]) is False


def test_the_census_totals_match_the_hit_counts():
    gold = _gold()
    for accession, cells in gold["band_census"].items():
        total = sum(v for cell in cells.values() for v in cell.values())
        assert total == gold["baseline_hits"][accession], accession
        banded = {band: sum(cell.values()) for band, cell in cells.items()}
        assert banded == gold["baseline_band_sizes"][accession], accession


def test_the_twilight_band_holds_both_labels():
    """The campaign is named for the band where identity alone cannot separate a
    family relationship from a superfamily one. If no band held both, the
    campaign would be measuring a clean separation and the build refuses."""
    gold, rules = _gold(), _rules()
    overlapping = [(a, band) for a, cells in gold["band_census"].items()
                   for band, cell in cells.items()
                   if cell["same_family"] and cell["superfamily_only"]]
    assert overlapping
    middle = rules.band_labels[1]
    assert any(band == middle for _, band in overlapping), overlapping


# ----------------------------------------------------------------- the verdicts


@pytest.mark.parametrize("cell,expected", [
    ({"same_family": 3, "superfamily_only": 2, "unclassified_pair": 0}, "overlapping"),
    ({"same_family": 3, "superfamily_only": 0, "unclassified_pair": 0}, "separable"),
    ({"same_family": 0, "superfamily_only": 4, "unclassified_pair": 0}, "search_only"),
    ({"same_family": 0, "superfamily_only": 0, "unclassified_pair": 0}, "unresolved"),
])
def test_the_verdict_is_computed_from_the_census(corpus, cell, expected):
    assert corpus.verdict(cell) == expected


def test_an_unresolved_band_is_a_fact_about_the_search(corpus):
    """`unresolved` means the band holds no hit at all. That is a statement about
    what the search reported, not a verdict about homology, and the vocabulary
    note says so."""
    rules = json.loads((CAMPAIGN / "reference" / "homology_rules.json").read_text())
    assert "not a verdict about homology" in rules["homology"]["verdicts_note"]
    assert corpus.verdict({"same_family": 0, "superfamily_only": 0,
                           "unclassified_pair": 7}) == "unresolved"


def test_the_verdict_census_counts_every_cell():
    gold, rules = _gold(), _rules()
    assert set(gold["verdict_census"]) == set(rules.verdicts)
    cells = sum(len(c) for c in gold["band_census"].values())
    assert sum(gold["verdict_census"].values()) == cells


def test_the_verdicts_agree_with_the_census_they_came_from(corpus):
    gold = _gold()
    for accession, cells in gold["band_census"].items():
        for band, cell in cells.items():
            assert gold["band_verdicts"][accession][band] == corpus.verdict(cell)


# ----------------------------------------------------------- the counterfactual


def test_the_zero_window_control_reproduces_the_baseline():
    """The unperturbed query re-run must return its own baseline hit set. Anything
    else means the search is not deterministic under this invocation, and every
    loss above it is uncalibrated."""
    gold = _gold()
    control = gold["shuffle_ladder"]["0"]
    for accession, row in control.items():
        assert row["lost"] == 0, accession
        assert row["hits"] == gold["baseline_hits"][accession], accession


def test_the_widest_window_costs_something():
    gold = _gold()
    top = max(gold["shuffle_ladder"], key=int)
    assert any(row["lost"] for row in gold["shuffle_ladder"][top].values())


def test_the_ladder_covers_every_declared_window_and_query():
    gold, rules = _gold(), _rules()
    assert set(gold["shuffle_ladder"]) == {str(w) for w in rules.shuffle_windows}
    accessions = {row["accession"] for row in gold["query_table"].values()}
    for window, rows in gold["shuffle_ladder"].items():
        assert set(rows) == accessions, window


def test_the_distant_hits_die_first_where_anything_dies():
    """The R4 claim, measured per cell rather than asserted: a local shuffle costs
    the weakest hits first, so the loss rates fall from the distant band to the
    close one. Reported per cell as a boolean, because a cell with no loss carries
    no ordering."""
    gold = _gold()
    cells = [(w, a, row) for w, rows in gold["shuffle_ladder"].items()
             for a, row in rows.items() if row["lost"]]
    assert cells, "no cell loses anything, so the ladder measures nothing"
    assert all(row["band_order_held"] for _, _, row in cells), \
        [(w, a) for w, a, row in cells if not row["band_order_held"]]


def test_a_rate_without_a_denominator_is_null_rather_than_zero():
    """A band with no baseline hit has no rate. Reporting zero would say the
    shuffle cost nothing there, which is a different claim from having nothing to
    lose."""
    gold = _gold()
    for window, rows in gold["shuffle_ladder"].items():
        for accession, row in rows.items():
            sizes = gold["baseline_band_sizes"][accession]
            for band, rate in row["loss_rate_by_band"].items():
                if sizes[band] == 0:
                    assert rate is None, (window, accession, band)
                else:
                    assert rate is not None, (window, accession, band)


def test_the_loss_rates_are_fractions_of_the_reported_denominators():
    gold = _gold()
    for window, rows in gold["shuffle_ladder"].items():
        for accession, row in rows.items():
            sizes = gold["baseline_band_sizes"][accession]
            for band, lost in row["lost_by_band"].items():
                if sizes[band]:
                    assert abs(row["loss_rate_by_band"][band]
                               - lost / sizes[band]) < 1e-6, (window, accession, band)


def test_the_band_order_test_skips_bands_with_no_rate():
    """Counting an absent band as zero would invent an ordering."""
    assert _ordered([0.9, None, 0.1]) is True
    assert _ordered([None, 0.2, 0.5]) is False
    assert _ordered([None, None, None]) is True


def test_the_shuffle_touches_only_the_declared_window(corpus):
    accession = sorted(corpus.queries().values())[0]
    original = corpus.sequences[accession]
    for window in (60, 120):
        mutated = corpus.shuffled(accession, window)
        assert len(mutated) == len(original)
        assert sorted(mutated) == sorted(original)
        middle = len(original) // 2
        low = max(0, middle - window // 2)
        high = min(len(original), low + window)
        assert mutated[:low] == original[:low]
        assert mutated[high:] == original[high:]
        assert mutated != original


def test_the_shuffle_is_reproducible_and_query_specific(corpus):
    """A shuffle has to be reproducible for the counterfactual to be gold, and a
    seed derived from the query means the same window on two queries is not the
    same permutation."""
    accessions = sorted(corpus.queries().values())[:2]
    assert corpus.shuffled(accessions[0], 60) == corpus.shuffled(accessions[0], 60)
    assert corpus.shuffled(accessions[0], 60) != corpus.shuffled(accessions[0], 120)
    first = corpus.shuffled(accessions[0], 60)
    second = corpus.shuffled(accessions[1], 60)
    middle = lambda s: s[len(s) // 2 - 30:len(s) // 2 + 30]
    assert middle(first) != middle(second)


def test_a_window_longer_than_the_query_is_refused(corpus):
    """At that point the ladder's top rung is a whole-sequence shuffle, which is
    the degenerate perturbation the catalogued R4 asked for."""
    accession = sorted(corpus.queries().values())[0]
    with pytest.raises(HomologyError, match="longer than"):
        corpus.shuffled(accession, len(corpus.sequences[accession]) + 1)


# ------------------------------------------------- the build's own refusals


def _mutated(**overrides) -> dict:
    gold = copy.deepcopy(_gold())
    gold.update(overrides)
    return gold


def _tiny_corpus() -> dict:
    return _mutated(corpus_sequences=12)


def _two_queries() -> dict:
    gold = copy.deepcopy(_gold())
    keep = sorted(gold["query_table"])[:2]
    gold["query_table"] = {k: gold["query_table"][k] for k in keep}
    return gold


def _one_query_twice() -> dict:
    gold = copy.deepcopy(_gold())
    entries = sorted(gold["query_table"])
    gold["query_table"][entries[1]]["accession"] = \
        gold["query_table"][entries[0]]["accession"]
    return gold


def _thin_family() -> dict:
    gold = copy.deepcopy(_gold())
    gold["query_table"][sorted(gold["query_table"])[0]]["family_members"] = 4
    return gold


def _mismatched_hit_counts() -> dict:
    gold = copy.deepcopy(_gold())
    accession = sorted(gold["baseline_hits"])[0]
    gold["baseline_hits"]["not_a_declared_query"] = gold["baseline_hits"].pop(accession)
    return gold


def _a_silent_query() -> dict:
    gold = copy.deepcopy(_gold())
    gold["baseline_hits"][sorted(gold["baseline_hits"])[0]] = 0
    return gold


def _nothing_missed() -> dict:
    gold = copy.deepcopy(_gold())
    for row in gold["recovery"].values():
        row["same_family_missed"] = 0
    return gold


def _a_query_with_no_relatives() -> dict:
    gold = copy.deepcopy(_gold())
    gold["recovery"][sorted(gold["recovery"])[0]]["same_family_asserted"] = 0
    return gold


def _no_overlap_anywhere() -> dict:
    gold = copy.deepcopy(_gold())
    for cells in gold["band_census"].values():
        for cell in cells.values():
            cell["superfamily_only"] = 0
    return gold


def _one_verdict_only() -> dict:
    gold = copy.deepcopy(_gold())
    cells = sum(len(c) for c in gold["band_census"].values())
    gold["verdict_census"] = {"separable": cells, "overlapping": 0,
                              "search_only": 0, "unresolved": 0}
    return gold


def _no_overlapping_verdict() -> dict:
    gold = copy.deepcopy(_gold())
    cells = sum(len(c) for c in gold["band_census"].values())
    gold["verdict_census"] = {"separable": cells - 1, "overlapping": 0,
                              "search_only": 1, "unresolved": 0}
    return gold


def _miscounted_verdicts() -> dict:
    gold = copy.deepcopy(_gold())
    gold["verdict_census"]["overlapping"] += 1
    return gold


def _no_control_rung() -> dict:
    gold = copy.deepcopy(_gold())
    gold["shuffle_ladder"] = {k: v for k, v in gold["shuffle_ladder"].items()
                              if k != "0"}
    return gold


def _a_noisy_control() -> dict:
    gold = copy.deepcopy(_gold())
    accession = sorted(gold["shuffle_ladder"]["0"])[0]
    gold["shuffle_ladder"]["0"][accession]["lost"] = 3
    return gold


def _a_ladder_that_costs_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    top = max(gold["shuffle_ladder"], key=int)
    for row in gold["shuffle_ladder"][top].values():
        row["lost"] = 0
    return gold


def _a_two_rung_ladder() -> dict:
    gold = copy.deepcopy(_gold())
    keep = sorted(gold["shuffle_ladder"], key=int)[:2]
    gold["shuffle_ladder"] = {k: gold["shuffle_ladder"][k] for k in keep}
    for row in gold["shuffle_ladder"][keep[-1]].values():
        row["lost"] = 1
    return gold


@pytest.mark.parametrize("mutate,pattern", [
    (_tiny_corpus, "cannot carry a per-band census"),
    (_two_queries, "cannot carry a per-query census"),
    (_one_query_twice, "two declared families selected the same query"),
    (_thin_family, "fewer than ten corpus members"),
    (_mismatched_hit_counts, "describe different\n?.*query sets"),
    (_a_silent_query, "which is a broken run"),
    (_nothing_missed, "recovers every relationship"),
    (_a_query_with_no_relatives, "no asserted family relative"),
    (_no_overlap_anywhere, "does not occur"),
    (_one_verdict_only, "not separating anything"),
    (_no_overlapping_verdict, "no band is classified as overlapping"),
    (_miscounted_verdicts, "the verdict census counts"),
    (_no_control_rung, "no zero-window rung"),
    (_a_noisy_control, "zero-window control loses hits"),
    (_a_ladder_that_costs_nothing, "costs no query a single hit"),
    (_a_two_rung_ladder, "cannot show a response to window"),
])
def test_the_build_refuses_gold_that_says_nothing(mutate, pattern):
    with pytest.raises(HomologyError, match=pattern):
        _assert_has_content(mutate())


def test_the_real_gold_passes_every_refusal():
    _assert_has_content(_gold())


# ---------------------------------------------------------------- the ladder


def test_the_oracle_clears_every_rung(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "homology_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["score"] == 1.0
    assert result["depth"] == 4
    assert result["monotonic"] is True
    assert result["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_an_empty_submission_clears_nothing(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    assert grade(measure(CAMPAIGN, submission), _gold(), _spec())["depth"] == 0


def test_the_warm_bundles_hand_over_exactly_the_preceding_rung():
    gold = _gold()
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", RECOVERY_KEYS, "recovery.json"),
                             ("r4", VERDICT_KEYS, "verdicts.json")):
        payload = json.loads(
            (CAMPAIGN / "gold" / "warm" / rung / name).read_text())
        assert set(payload) == set(keys)
        assert all(payload[k] == gold[k] for k in keys)


def test_the_rung_key_sets_do_not_overlap():
    sets = [set(INVENTORY_KEYS), set(RECOVERY_KEYS), set(VERDICT_KEYS),
            set(COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right


def test_the_echoes_are_reported_and_not_graded():
    graded = {c["name"] for r in _spec()["rungs"] for c in r["components"]}
    gold = _gold()
    for echo in ("label_resource", "corpus_superfamily", "sensitivity", "evalue",
                 "max_seqs", "threads", "format_fields"):
        assert echo not in graded
        assert echo in gold


def test_the_search_identity_is_graded_as_one_claim():
    """The binary's version, the sensitivity, the reporting threshold, the cap,
    the thread count and the output projection are one claim about what made the
    hits. The flags are declared in the rules file; the version is not, which is
    what keeps the component from being a copy."""
    gold = _gold()
    assert set(gold["environment"]) == {"mmseqs_version", "sensitivity", "evalue",
                                        "max_seqs", "threads", "format_fields"}
    graded = {c["name"] for r in _spec()["rungs"] for c in r["components"]}
    assert "environment" in graded
    assert "mmseqs_version" not in graded
    rules = json.loads((CAMPAIGN / "reference" / "homology_rules.json").read_text())
    assert gold["mmseqs_version"] not in json.dumps(rules)


def test_the_catalogued_judgement_is_excluded_by_name():
    excluded = [str(b["field"]) for b in _task()["excluded_from_grading"]]
    assert any("really homologous" in field for field in excluded)


def test_nothing_is_sampled_filtered_or_planted():
    """This campaign declared `structural_elements: [planted]` and plants
    nothing: its provenance says the whole reviewed set is the corpus. The
    perturbation it does carry is a shuffle it performs at run time, under a
    declared window, which is a rule rather than an artefact. The label is gone;
    the two claims it stood for are asserted here instead.
    """
    provenance = json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())
    rule = provenance["corpus"]["selection_rule"]
    assert "Nothing is sampled or filtered" in rule
    rules = json.loads(
        (CAMPAIGN / "reference" / "homology_rules.json").read_text())["homology"]
    assert rules["shuffle_windows"][0] == 0, (
        "the zero-width shuffle is the control that makes the perturbation a "
        "declared rule rather than a planted corpus")


# --------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    copy_dir = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, copy_dir)
    shutil.rmtree(copy_dir / "gold")
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "homology_report.json").write_text(
        json.dumps(render(_gold()), indent=2, sort_keys=True) + "\n")
    assert measure(copy_dir, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_recomputes_the_query_set_rather_than_trusting_it(tmp_path):
    """A submission that invented its own query accessions scores zero rather
    than being read on its own terms."""
    report = render(_gold())
    report["recovery"] = {f"X{k}": v for k, v in report["recovery"].items()}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "homology_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    assert grade(measure(CAMPAIGN, submission), _gold(), _spec())["depth"] == 1


def test_a_missing_rate_is_not_equal_to_a_rate_of_zero(tmp_path):
    report = render(_gold())
    for rows in report["shuffle_ladder"].values():
        for row in rows.values():
            row["loss_rate_by_band"] = {
                band: (0.0 if rate is None else rate)
                for band, rate in row["loss_rate_by_band"].items()}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "homology_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    measured = measure(CAMPAIGN, submission)
    gold = _gold()
    nulls = [(w, a, b) for w, rows in gold["shuffle_ladder"].items()
             for a, row in rows.items()
             for b, rate in row["loss_rate_by_band"].items() if rate is None]
    if not nulls:
        pytest.skip("no null rate in gold on this corpus")
    assert grade(measured, gold, _spec())["depth"] == 3


def test_the_measurer_keys_the_ladder_by_the_declared_windows(tmp_path):
    report = render(_gold())
    report["shuffle_ladder"] = {f"w{k}": v
                                for k, v in report["shuffle_ladder"].items()}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / "homology_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    assert grade(measure(CAMPAIGN, submission), _gold(), _spec())["depth"] == 3


# ---------------------------------------------------------- the self-hit probe


def test_the_self_hit_is_excluded_from_every_count():
    """The query is a member of the corpus it is searched against, so it returns
    itself at 100 per cent identity. That hit is the query, not a finding."""
    gold = _gold()
    for accession, cells in gold["band_census"].items():
        # A self-hit would be a 100 per cent same-family hit in the top band.
        top = gold["baseline_band_sizes"][accession][_rules().band_labels[-1]]
        assert sum(cells[_rules().band_labels[-1]].values()) == top
    for accession, row in gold["recovery"].items():
        assert row["same_family_detected"] <= row["same_family_asserted"]


@image_present
def test_keeping_the_self_hit_shifts_every_count_by_exactly_one(corpus):
    """The competence probe, measured rather than asserted: the symptom is an
    off-by-one that stays internally consistent, which is why it is caught at R1
    where the counts are rather than at R2."""
    accessions = sorted(corpus.queries().values())
    baseline = corpus.search({a: corpus.sequences[a] for a in accessions})
    gold = _gold()
    for accession in accessions:
        with_self = sum(1 for h in baseline if h.query == accession)
        assert with_self == gold["baseline_hits"][accession] + 1, accession
        self_hits = [h for h in baseline
                     if h.query == accession and h.target == accession]
        assert len(self_hits) == 1
        assert self_hits[0].identity == pytest.approx(1.0)
