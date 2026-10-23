"""Everything the dereplication campaign grades, computed from the pinned rules.

The oracle's analysis and the measurer's reading share no code: gold is this
file, and `measure.py` reads a submission. Nothing here reads gold.

Each rung's projection is a module-level function rather than a block inside
`analyse`, because the stub ladder needs exactly the first N of them and a stub
that reimplemented a census would drift from the one gold is built with.
"""

from __future__ import annotations

import json
import pathlib
import statistics

import yaml

from .core import Corpus, DereplicationError, PRECISION

#: What one query's search can come to. `no_candidate` is an empty precursor
#: window; `no_match` is a window with candidates that all scored zero, which
#: matchms omits from its sparse score matrix rather than reporting as a zero.
RETRIEVAL_VERDICTS = ("annotated", "no_candidate", "no_match")

#: What the adduct misassignment does to one annotation. Measured on this corpus
#: it mostly empties the window, which is why the rung grades the census and not
#: a per-query survival bit.
ADDUCT_OUTCOMES = ("no_candidate", "lost", "survives")


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


# ------------------------------------------------------------------ R1: execute

def environment(corpus: Corpus) -> dict:
    """What ran the comparison. The version is the half the sandbox cannot state."""
    return {
        "matchms_version": corpus.matchms_version(),
        "similarity": corpus.rules.similarity,
        "peak_tolerance": corpus.rules.peak_tolerance,
        "precursor_window": corpus.rules.precursor_window,
    }


def corpus_census(corpus: Corpus) -> dict:
    """The two files as parsed. Reported for the record, not graded: all three
    numbers are in `inputs/provenance.json`, which the agent is given."""
    return {
        "graded_queries": len(corpus.queries),
        "library_spectra": len(corpus.library),
        "library_compounds": len({row["compound_block"] for row in corpus.library}),
    }


def candidate_census(sizes: list[int]) -> dict:
    """The precursor filter applied to the whole library, in three numbers.

    This is R1's work. The filter is part of the declared pipeline, so a system
    that only counted the shipped spectra has not run it, and the window's shape
    over this library is stated nowhere the agent can read.

    Takes the candidate-set sizes rather than annotation rows, because the filter
    is upstream of the search: R1 is answerable without scoring anything, and the
    stub ladder's first level says so by never calling the tool for it.
    """
    return {"median": statistics.median(sizes), "max": max(sizes),
            "total": sum(sizes)}


# ------------------------------------------------------------------ R2: correct

def retrieval(rows: list[dict]) -> dict:
    """What the search answered, per query and in aggregate."""
    census = {verdict: 0 for verdict in RETRIEVAL_VERDICTS}
    correct, wrong = [], []
    for row in rows:
        if row["verdict"] not in census:
            raise DereplicationError(
                f"{row['query']} came back with the verdict {row['verdict']!r}, "
                "which is not one this template declares")
        census[row["verdict"]] += 1
        if row["verdict"] == "annotated":
            (correct if row["correct"] else wrong).append(row["query"])
    return {
        "retrieval_census": census,
        "correct_queries": sorted(correct),
        "misannotated_queries": sorted(wrong),
        "top_matches": {row["query"]: row["top_match"] for row in rows},
    }


# ------------------------------------------------------------------- R3: commit

def commitments(rows: list[dict], rules) -> dict:
    """The declared confidence class per query, and the ambiguity it admits to."""
    tally = {label: 0 for label in rules.confidence_vocabulary}
    for row in rows:
        if row["confidence"] not in tally:
            raise DereplicationError(
                f"{row['query']} was classified {row['confidence']!r}, which is "
                "not in the declared confidence vocabulary")
        tally[row["confidence"]] += 1

    margins = [row["margin_to_other_compound"] for row in rows
               if row["margin_to_other_compound"] is not None]
    ambiguous = sorted(row["query"] for row in rows
                       if row["margin_to_other_compound"] is not None
                       and row["margin_to_other_compound"] <= rules.ambiguity_margin)
    return {
        "confidence_by_query": {row["query"]: row["confidence"] for row in rows},
        "confidence_census": {label: tally[label] for label in sorted(tally)},
        "margin_summary": {
            "queries_with_a_competing_compound": len(margins),
            "median": round(statistics.median(margins), PRECISION) if margins else None,
            "within_the_declared_margin": len(ambiguous),
        },
        "ambiguous_queries": ambiguous,
    }


# ----------------------------------------------------------- R4: counterfactual

def tolerance_movement(corpus: Corpus, rows: list[dict]) -> dict:
    """Per declared tolerance, the queries whose answer is not the pinned one.

    The pinned tolerance is the identity control and its set is empty by
    construction; the build refuses gold if it is not.
    """
    by_query = {row["query"]: row for row in rows}
    movement = {}
    for tolerance in corpus.rules.tolerance_grid:
        key = f"{tolerance}"
        if tolerance == corpus.rules.peak_tolerance:
            movement[key] = []
            continue
        moved = []
        for row in corpus.annotations(tolerance=tolerance):
            base = by_query[row["query"]]
            if (row["top_match"], row["verdict"]) != (base["top_match"],
                                                      base["verdict"]):
                moved.append(row["query"])
        movement[key] = sorted(moved)
    return movement


def adduct_census(corpus: Corpus) -> dict:
    """What the search reports when the precursor is read as the wrong adduct.

    The spectra do not change. The precursor window moves, which is exactly what
    a misassigned adduct does to a real dereplication run.
    """
    rules = corpus.rules
    centres = {query["accession"]:
               query["precursor_mz"] - rules.sodium_mass + rules.proton_mass
               for query in corpus.queries}
    census = {label: 0 for label in ADDUCT_OUTCOMES}
    for row in corpus.annotations(centres=centres):
        if row["verdict"] in ("no_candidate", "no_match"):
            census["no_candidate"] += 1
        elif row["correct"]:
            census["survives"] += 1
        else:
            census["lost"] += 1
    return census


# ------------------------------------------------------------------- the whole

def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["queries_input"], params["library_input"])
    rows = corpus.annotations()

    gold = {
        "environment": environment(corpus),
        "corpus_census": corpus_census(corpus),
        "candidate_census": candidate_census(
            [row["candidates"] for row in rows]),
        **retrieval(rows),
        **commitments(rows, corpus.rules),
        "tolerance_movement": tolerance_movement(corpus, rows),
        "adduct_misassignment_census": adduct_census(corpus),
    }
    if not gold["correct_queries"]:
        raise DereplicationError(
            "no query retrieves its own compound; the split or the library is "
            "wrong and nothing above R1 would be gradeable")
    return gold


def main(argv: list[str] | None = None) -> int:
    import sys
    campaign = pathlib.Path((argv or sys.argv[1:])[0])
    gold = analyse(campaign)
    bulky = ("top_matches", "correct_queries", "misannotated_queries",
             "ambiguous_queries", "tolerance_movement", "confidence_by_query")
    print(json.dumps({k: v for k, v in gold.items() if k not in bulky}, indent=1))
    print("correct", len(gold["correct_queries"]),
          "wrong", len(gold["misannotated_queries"]),
          "ambiguous", len(gold["ambiguous_queries"]))
    print("tolerance movement:",
          {k: len(v) for k, v in gold["tolerance_movement"].items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
