"""MS2 dereplication: grading the retrieval, and the split that makes it one.

What this campaign claims is simple enough to state in a sentence and was false
the first time it was stated: a correct answer is another lab's measurement of
the same compound. The split first shipped as "the query is one contributor's
spectrum and every other spectrum of that compound joins the library", which
leaves the query's own lab's replicates in the library -- nineteen of the hundred
graded compounds had one, and ten of the ninety correct answers came from one at
the same instrument type too. So the exclusion is a rule now, enforced in the
input builder and again in the grading core's loader, and the tests below hold
both ends of it.

Three other things here produce a wrong answer without producing an error.

The declared ranking is score, then matched peaks, then accession. Ordering by
matched peaks first is a defensible reading of a spectral comparison and it is
not this one: it moves most of the top matches. That is one competence probe.

The declared confidence rule is a PAIR of floors and the first row meeting both
wins. Reading it as a score threshold alone moves a sixth of the annotations
upward into the most confident class. That is the other.

And the two silent states are not one state. A query whose precursor window is
empty has no candidate; a query whose candidates all score zero has no match,
because a sparse score matrix omits a zero-scoring pair rather than reporting it.
Collapsing them loses the distinction R4's census is built on.

The scoring is slow enough that the suite scores the corpus ONCE, in a
module-scoped fixture, and derives every ranking variant from that one matrix.
Everything else reads gold, which already records the answer.
"""

from __future__ import annotations

import copy
import gzip
import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.sweep.stubs.stub_derep import _peaks_first
from npbench_c.templates.ms2_dereplication import analyse as derep
from npbench_c.templates.ms2_dereplication.build import (
    COMMITMENT_KEYS,
    COUNTERFACTUAL_KEYS,
    EXECUTE_KEYS,
    RETRIEVAL_KEYS,
    _assert_has_content,
)
from npbench_c.templates.ms2_dereplication.core import (
    Corpus,
    DereplicationError,
    Rules,
    declared_ranking,
)
from npbench_c.templates.ms2_dereplication.measure import measure
from npbench_c.templates.ms2_dereplication.solve import REPORT, render
from npbench_c.tools.resources import tool_prefix

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "derep-massbank-2026_03-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


def _spec() -> dict:
    return json.loads((CAMPAIGN / "grading.json").read_text())


def _provenance() -> dict:
    return json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())


def _rules() -> Rules:
    return Rules.load(CAMPAIGN)


@pytest.fixture(scope="module")
def corpus() -> Corpus:
    params = derep.load_params(CAMPAIGN)
    return Corpus.load(CAMPAIGN, params["queries_input"], params["library_input"])


@pytest.fixture(scope="module")
def scored(corpus) -> dict:
    """The candidate scores, once. Every ranking variant comes from this."""
    tasks = [{"query": query, "candidates": corpus.candidates(query)}
             for query in corpus.queries]
    raw = corpus._score(tasks, corpus.rules.peak_tolerance)
    return {accession: [[str(name), float(score), int(matches)]
                        for name, score, matches in rows]
            for accession, rows in raw.items()}


# ------------------------------------------------------------ what is shipped


def test_only_the_two_corpora_and_their_provenance_are_shipped():
    """No annotation of the queries, and no perturbed search. Both follow from
    the declared rules, and either one would hand over a rung."""
    assert sorted(p.name for p in (CAMPAIGN / "inputs").iterdir()) == [
        "library.json.gz", "provenance.json", "queries.json.gz"]


def test_the_queries_carry_no_annotation_of_their_own_identity(corpus):
    """`compound_block` is on both sides because the library needs identities to
    be a reference library -- but a query's is gold, so it must not be readable
    as an annotation. It is not: it is the compound the query IS, which is the
    answer, and the test that matters is that the query is not in the library."""
    library = {row["accession"] for row in corpus.library}
    assert not [query for query in corpus.queries
                if query["accession"] in library]


def test_every_query_compound_is_in_the_library(corpus):
    blocks = {row["compound_block"] for row in corpus.library}
    assert all(query["compound_block"] in blocks for query in corpus.queries)


def test_no_query_is_answerable_from_its_own_contributor(corpus):
    """The campaign's whole construct, checked on the shipped files."""
    for query in corpus.queries:
        own = [row for row in corpus.library
               if row["compound_block"] == query["compound_block"]]
        assert own, query["accession"]
        assert not [row for row in own
                    if row["contributor"] == query["contributor"]]


def test_the_loader_refuses_a_corpus_answerable_from_its_own_contributor(tmp_path):
    """The rule the docstring used to assert on its own. A replicate of the
    query's own lab in the library makes a correct answer a lookup."""
    params = derep.load_params(CAMPAIGN)
    queries = json.loads(gzip.open(CAMPAIGN / params["queries_input"], "rt").read())
    library = json.loads(gzip.open(CAMPAIGN / params["library_input"], "rt").read())
    victim = queries[0]
    replicate = dict(victim)
    replicate["accession"] = victim["accession"] + "-REPLICATE"

    campaign = tmp_path / "campaign"
    (campaign / "inputs").mkdir(parents=True)
    (campaign / "reference").mkdir()
    (campaign / "reference" / "dereplication_rules.json").write_text(
        (CAMPAIGN / "reference" / "dereplication_rules.json").read_text())
    with gzip.open(campaign / "inputs" / "queries.json.gz", "wt") as handle:
        json.dump(queries, handle)
    with gzip.open(campaign / "inputs" / "library.json.gz", "wt") as handle:
        json.dump(library[:200] + [replicate], handle)

    with pytest.raises(DereplicationError, match="same lab's replicate"):
        Corpus.load(campaign, "inputs/queries.json.gz", "inputs/library.json.gz")


def test_the_provenance_reports_the_instrument_support_rather_than_claiming_it():
    """The contributor exclusion guarantees another lab and not another
    instrument. Measured: every graded query does have cross-instrument support
    on this corpus, which is a fact about the corpus and not a rule."""
    support = _provenance()["instrument_support"]
    assert support == {"cross_instrument": 100, "same_instrument_only": 0}
    assert sum(support.values()) == _provenance()["graded_queries"]


def test_the_provenance_states_the_census_and_not_the_filter():
    """`corpus_census` is excluded from grading because provenance states it.
    `candidate_census` is graded because provenance cannot: the filter is work
    over the shipped files rather than a property of them."""
    blob = (CAMPAIGN / "inputs" / "provenance.json").read_text()
    gold = _gold()
    assert str(gold["corpus_census"]["library_spectra"]) in blob
    assert str(gold["corpus_census"]["graded_queries"]) in blob
    for value in gold["candidate_census"].values():
        assert str(value) not in blob


# --------------------------------------------------------------- the pipeline


def test_the_precursor_filter_is_the_declared_window(corpus):
    """Recomputed independently of `candidates()`: a candidate is a library
    spectrum within the declared window, and nothing else is."""
    window = corpus.rules.precursor_window
    query = corpus.queries[0]
    expected = sorted(row["accession"] for row in corpus.library
                      if abs(row["precursor_mz"] - query["precursor_mz"]) <= window)
    assert sorted(row["accession"]
                  for row in corpus.candidates(query)) == expected


def test_the_candidate_census_is_the_filter_and_not_the_library(corpus):
    gold = _gold()
    sizes = [len(corpus.candidates(query)) for query in corpus.queries]
    assert derep.candidate_census(sizes) == gold["candidate_census"]
    assert gold["candidate_census"]["total"] < gold["corpus_census"]["library_spectra"]


def test_the_declared_ranking_is_total_and_deterministic():
    rows = [["b", 0.5, 3], ["a", 0.5, 3], ["c", 0.9, 1], ["d", 0.5, 9]]
    assert [row[0] for row in sorted(rows, key=declared_ranking)] == \
        ["c", "d", "a", "b"]


@image_present
def test_the_ranking_tie_breaks_change_no_answer_but_the_ties_are_real(scored):
    """Measured, and the reason the rule stays anyway. Eight queries have tied
    candidate scores, so the order matchms returned would decide an answer on a
    corpus where a tie reached the top; here none does, which is why dropping
    the tie-break is NOT the stub ladder's probe."""
    tied = [accession for accession, rows in scored.items()
            if len(rows) > 1
            and len({round(row[1], 6) for row in rows}) < len(rows)]
    assert len(tied) == 8

    top = _gold()["top_matches"]
    for accession, rows in scored.items():
        if not rows:
            continue
        assert sorted(rows, key=lambda row: -row[1])[0][0] == top[accession]
        assert sorted(rows, key=lambda row: (-row[1], -row[2]))[0][0] == \
            top[accession]


@image_present
def test_ranking_by_matched_peaks_first_moves_most_of_the_answers(scored):
    """The first competence probe, and the measurement that makes it one."""
    top = _gold()["top_matches"]
    moved = [accession for accession, rows in scored.items()
             if rows and sorted(rows, key=_peaks_first)[0][0] != top[accession]]
    assert len(moved) == 65


def test_the_confidence_rule_needs_both_floors(corpus):
    """Both floors, and the first row meeting both wins."""
    rules = corpus.rules
    assert rules.classify(0.95, 6) == "confirmed"
    assert rules.classify(0.95, 5) == "tentative"   # the score alone is not enough
    assert rules.classify(0.95, 2) == "unsupported"


@image_present
def test_reading_the_confidence_rule_as_a_score_threshold_moves_a_sixth_upward(
        corpus, scored):
    """The second competence probe, and the measurement that makes it one. Every
    one of the moves is upward into the most confident class: the direction that
    makes a system look more certain than its own evidence supports."""
    rules = corpus.rules
    gold = _gold()
    moved = []
    for accession, rows in scored.items():
        if not rows:
            continue
        _, score, matches = sorted(rows, key=declared_ranking)[0]
        declared = rules.classify(score, matches)
        loose = next((label for floor, _, label in rules.confidence_rule
                      if score >= floor), rules.confidence_fallback)
        assert declared == gold["confidence_by_query"][accession]
        if loose != declared:
            moved.append((declared, loose))
    assert len(moved) == 17
    assert {pair[1] for pair in moved} == {"confirmed"}


def test_a_window_with_candidates_that_all_score_zero_is_not_an_empty_window(
        corpus, monkeypatch):
    """`no_match` and `no_candidate` are different states. matchms omits a
    zero-scoring pair from its sparse matrix, so a query with candidates and no
    score rows has been measured rather than skipped."""
    assert {"no_candidate", "no_match"} < set(derep.RETRIEVAL_VERDICTS)
    assert set(_gold()["retrieval_census"]) == set(derep.RETRIEVAL_VERDICTS)

    monkeypatch.setattr(Corpus, "_score", lambda self, tasks, tolerance: {})
    rows = corpus.annotations()
    census = derep.retrieval(rows)["retrieval_census"]
    assert census["annotated"] == 0
    assert census["no_match"] == sum(1 for row in rows if row["candidates"])
    assert census["no_candidate"] == sum(1 for row in rows
                                         if not row["candidates"])
    assert census["no_candidate"] == _gold()["retrieval_census"]["no_candidate"]


# ------------------------------------------------------------- the arithmetic


def test_the_retrieval_arithmetic_closes():
    gold = _gold()
    graded = gold["corpus_census"]["graded_queries"]
    assert sum(gold["retrieval_census"].values()) == graded
    assert len(gold["top_matches"]) == graded
    assert len(gold["correct_queries"]) + len(gold["misannotated_queries"]) == \
        gold["retrieval_census"]["annotated"]
    assert not set(gold["correct_queries"]) & set(gold["misannotated_queries"])


def test_a_query_with_no_candidate_has_no_top_match():
    gold = _gold()
    silent = [accession for accession, match in gold["top_matches"].items()
              if match is None]
    assert len(silent) == (gold["retrieval_census"]["no_candidate"]
                           + gold["retrieval_census"]["no_match"])


def test_the_confidence_census_is_the_tally_of_the_per_query_classes():
    gold = _gold()
    tally: dict[str, int] = {label: 0 for label in gold["confidence_census"]}
    for label in gold["confidence_by_query"].values():
        tally[label] += 1
    assert tally == gold["confidence_census"]
    assert sum(tally.values()) == gold["corpus_census"]["graded_queries"]


def test_a_silent_query_is_classified_as_silent_and_not_as_unsupported():
    gold = _gold()
    for accession, match in gold["top_matches"].items():
        if match is None:
            assert gold["confidence_by_query"][accession] in (
                "no_candidate", "no_match")


def test_the_ambiguous_set_is_the_margin_summary_s_own_count():
    gold = _gold()
    margin = gold["margin_summary"]
    assert margin["within_the_declared_margin"] == len(gold["ambiguous_queries"])
    assert margin["queries_with_a_competing_compound"] > \
        len(gold["ambiguous_queries"])


def test_the_ambiguous_set_is_a_subset_of_the_annotated_queries():
    gold = _gold()
    annotated = {accession for accession, match in gold["top_matches"].items()
                 if match is not None}
    assert set(gold["ambiguous_queries"]) <= annotated


def test_the_margin_median_is_over_the_competitor_population_not_the_corpus():
    """The denominator mistake this campaign is most likely to be audited for.
    The median is over the queries that have a competing compound at all, which
    is fewer than the graded queries and more than the ambiguous set."""
    gold = _gold()
    margin = gold["margin_summary"]
    assert margin["queries_with_a_competing_compound"] < \
        gold["corpus_census"]["graded_queries"]
    assert 0.0 < margin["median"] <= 1.0


def test_the_pinned_tolerance_is_the_identity_control():
    gold, rules = _gold(), _rules()
    assert gold["tolerance_movement"][f"{rules.peak_tolerance}"] == []
    others = [queries for key, queries in gold["tolerance_movement"].items()
              if key != f"{rules.peak_tolerance}"]
    assert all(others), "every other declared tolerance moves at least one answer"


def test_the_movement_sets_name_graded_queries_only():
    gold = _gold()
    graded = set(gold["top_matches"])
    for queries in gold["tolerance_movement"].values():
        assert set(queries) <= graded
        assert queries == sorted(queries)


def test_the_tolerance_grid_is_the_declared_one():
    gold, rules = _gold(), _rules()
    assert sorted(gold["tolerance_movement"]) == \
        sorted(f"{tolerance}" for tolerance in rules.tolerance_grid)


def test_the_adduct_census_covers_every_query_and_nothing_survives():
    """The measured shape, and the reason the rung grades a census: a misassigned
    adduct makes the search silent rather than wrong, so the per-query survival
    bit the catalogue asked for has one answer."""
    gold = _gold()
    census = gold["adduct_misassignment_census"]
    assert sum(census.values()) == gold["corpus_census"]["graded_queries"]
    assert census["survives"] == 0
    assert census["no_candidate"] > census["lost"] > 0


def test_the_adduct_shift_is_the_declared_mass_difference(corpus):
    rules = corpus.rules
    shift = rules.sodium_mass - rules.proton_mass
    assert 21.9 < shift < 22.0
    assert rules.assumed_adduct == "[M+H]+"
    assert rules.misread_adduct == "[M+Na]+"


# ------------------------------------------------------------- gold and grading


def test_the_rung_key_sets_partition_the_report():
    """A warm bundle hands over the rung below's answers and nothing else, so the
    key sets must be disjoint and must cover the whole report."""
    sets = [set(keys) for keys in (EXECUTE_KEYS, RETRIEVAL_KEYS, COMMITMENT_KEYS,
                                   COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right
    assert set().union(*sets) == set(render(_gold()))


def test_every_graded_component_resolves_in_gold():
    gold = _gold()
    for rung in _spec()["rungs"]:
        for component in rung["components"]:
            reference = component["args"]["gold"]
            if isinstance(reference, str) and reference.startswith("gold:"):
                assert reference.split(":", 1)[1] in gold


def test_the_echoes_are_reported_and_not_graded():
    graded = {component["args"]["gold"] for rung in _spec()["rungs"]
              for component in rung["components"]
              if isinstance(component["args"]["gold"], str)}
    assert "gold:corpus_census" not in graded
    assert "gold:misannotated_queries" not in graded
    assert {"corpus_census", "misannotated_queries"} <= set(render(_gold()))
    excluded = {entry["field"] for entry in _task()["excluded_from_grading"]}
    assert "corpus_census" in excluded
    assert "misannotated_queries" in excluded


def test_the_oracle_submission_scores_one(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / REPORT).write_text(json.dumps(render(_gold())))
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["score"] == 1.0
    assert result["depth"] == 4


def test_a_missing_report_scores_nothing(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["depth"] == 0


def test_the_measurer_reads_the_query_set_from_the_inputs(tmp_path):
    """A submission that answered for a query set of its own scores zero rather
    than being read on its own terms."""
    report = render(_gold())
    report["top_matches"] = {f"INVENTED-{i}": "MSBNK-SOMETHING"
                             for i in range(100)}
    report["confidence_by_query"] = {f"INVENTED-{i}": "confirmed"
                                     for i in range(100)}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / REPORT).write_text(json.dumps(report))
    measurement = measure(CAMPAIGN, submission)
    assert set(measurement["top_matches"]) == set(_gold()["top_matches"])
    assert all(value is None for value in measurement["top_matches"].values())
    result = grade(measurement, _gold(), _spec())
    assert result["depth"] == 1


def test_a_set_answer_is_normalised_by_order_only(tmp_path):
    """Sorted, because the claim is a set; not deduplicated, because naming a
    query twice is a different claim from the one gold makes."""
    gold = _gold()
    report = render(gold)
    report["correct_queries"] = list(reversed(gold["correct_queries"]))
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / REPORT).write_text(json.dumps(report))
    assert grade(measure(CAMPAIGN, submission), gold, _spec())["score"] == 1.0

    report["correct_queries"] = gold["correct_queries"] + gold["correct_queries"][:1]
    (submission / REPORT).write_text(json.dumps(report))
    assert grade(measure(CAMPAIGN, submission), gold, _spec())["depth"] == 1


# ----------------------------------------------------------- build refusals


def _mutated(**overrides) -> dict:
    gold = copy.deepcopy(_gold())
    gold.update(overrides)
    return gold


def _a_different_library_version() -> dict:
    gold = copy.deepcopy(_gold())
    gold["environment"]["matchms_version"] = "0.0.1"
    return gold


def _a_tiny_library() -> dict:
    gold = copy.deepcopy(_gold())
    gold["corpus_census"]["library_compounds"] = 12
    return gold


def _a_miscounted_retrieval() -> dict:
    gold = copy.deepcopy(_gold())
    gold["retrieval_census"]["annotated"] += 1
    return gold


def _a_filter_with_no_competition() -> dict:
    gold = copy.deepcopy(_gold())
    gold["candidate_census"]["max"] = 1
    return gold


def _a_broken_split() -> dict:
    gold = copy.deepcopy(_gold())
    keep = gold["correct_queries"][:40]
    gold["correct_queries"] = keep
    gold["misannotated_queries"] = sorted(
        set(gold["top_matches"]) - set(keep)
        - {a for a, m in gold["top_matches"].items() if m is None})
    gold["retrieval_census"] = {"annotated": 98, "no_candidate": 2, "no_match": 0}
    return gold


def _nothing_misannotated() -> dict:
    gold = copy.deepcopy(_gold())
    gold["correct_queries"] = sorted(set(gold["correct_queries"])
                                     | set(gold["misannotated_queries"]))
    gold["misannotated_queries"] = []
    return gold


def _a_query_counted_twice() -> dict:
    gold = copy.deepcopy(_gold())
    gold["misannotated_queries"] = sorted(
        set(gold["misannotated_queries"]) | {gold["correct_queries"][0]})
    return gold


def _a_missing_top_match() -> dict:
    gold = copy.deepcopy(_gold())
    gold["top_matches"].pop(sorted(gold["top_matches"])[0])
    return gold


def _an_empty_confidence_class() -> dict:
    gold = copy.deepcopy(_gold())
    census = gold["confidence_census"]
    census["confirmed"], census["tentative"] = 0, census["confirmed"] + census["tentative"]
    return gold


def _a_miscounted_confidence_census() -> dict:
    gold = copy.deepcopy(_gold())
    gold["confidence_census"]["confirmed"] += 1
    return gold


def _a_short_confidence_table() -> dict:
    gold = copy.deepcopy(_gold())
    gold["confidence_by_query"].pop(sorted(gold["confidence_by_query"])[0])
    return gold


def _too_few_competitors() -> dict:
    gold = copy.deepcopy(_gold())
    gold["margin_summary"]["queries_with_a_competing_compound"] = 4
    return gold


def _no_ambiguity() -> dict:
    gold = copy.deepcopy(_gold())
    gold["ambiguous_queries"] = []
    gold["margin_summary"]["within_the_declared_margin"] = 0
    return gold


def _every_competitor_inside_the_margin() -> dict:
    gold = copy.deepcopy(_gold())
    competitors = gold["margin_summary"]["queries_with_a_competing_compound"]
    gold["ambiguous_queries"] = sorted(gold["correct_queries"][:competitors])
    gold["margin_summary"]["within_the_declared_margin"] = competitors
    return gold


def _a_two_point_grid() -> dict:
    gold = copy.deepcopy(_gold())
    keep = sorted(gold["tolerance_movement"])[:2]
    gold["tolerance_movement"] = {k: gold["tolerance_movement"][k] for k in keep}
    return gold


def _a_moving_control() -> dict:
    gold = copy.deepcopy(_gold())
    gold["tolerance_movement"][f"{_rules().peak_tolerance}"] = \
        gold["correct_queries"][:3]
    return gold


def _a_grid_that_moves_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    gold["tolerance_movement"] = {k: [] for k in gold["tolerance_movement"]}
    return gold


def _a_miscounted_adduct_census() -> dict:
    gold = copy.deepcopy(_gold())
    gold["adduct_misassignment_census"]["lost"] += 1
    return gold


def _an_adduct_shift_that_silences_nobody() -> dict:
    gold = copy.deepcopy(_gold())
    census = gold["adduct_misassignment_census"]
    census["lost"] += census["no_candidate"]
    census["no_candidate"] = 0
    return gold


def _an_adduct_shift_that_costs_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    graded = gold["corpus_census"]["graded_queries"]
    survives = len(gold["correct_queries"])
    gold["adduct_misassignment_census"] = {
        "no_candidate": 5, "lost": graded - survives - 5, "survives": survives}
    return gold


@pytest.mark.parametrize("mutate,pattern", [
    (_a_different_library_version, "a library the campaign does not declare"),
    (_a_tiny_library, "cannot carry a dereplication question"),
    (_a_miscounted_retrieval, "the retrieval census covers"),
    (_a_filter_with_no_competition, "two candidates"),
    (_a_broken_split, "broken split rather than a hard search"),
    (_nothing_misannotated, "every annotated query retrieves its own compound"),
    (_a_query_counted_twice, "both correct and misannotated"),
    (_a_missing_top_match, "queries carry a top match"),
    (_an_empty_confidence_class, "are empty on this corpus"),
    (_a_miscounted_confidence_census, "the confidence census counts"),
    (_a_short_confidence_table, "queries carry a confidence"),
    (_too_few_competitors, "a competing compound exists for only"),
    (_no_ambiguity, "ambiguous set is empty"),
    (_every_competitor_inside_the_margin, "the margin separates nothing"),
    (_a_two_point_grid, "cannot show a response to the"),
    (_a_moving_control, "not deterministic under this invocation"),
    (_a_grid_that_moves_nothing, "could have left undeclared"),
    (_a_miscounted_adduct_census, "the adduct census counts"),
    (_an_adduct_shift_that_silences_nobody, "leaves every query a candidate"),
    (_an_adduct_shift_that_costs_nothing, "costs the search nothing"),
])
def test_the_build_refuses_gold_that_says_nothing(mutate, pattern):
    with pytest.raises(DereplicationError, match=pattern):
        _assert_has_content(CAMPAIGN, mutate())


def test_the_real_gold_passes_every_refusal():
    _assert_has_content(CAMPAIGN, _gold())
