"""Everything the dereplication campaign grades, computed from the pinned rules.

The oracle's analysis and the measurer's reading are deliberately the same code
path only as far as the SEARCH: gold is this file, the measurer reads a
submission. Nothing here reads gold.
"""

from __future__ import annotations

import json
import pathlib
import statistics

import yaml

from .core import Corpus, DereplicationError, PRECISION


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["queries_input"], params["library_input"])
    rules = corpus.rules

    baseline = corpus.annotations()
    by_query = {row["query"]: row for row in baseline}
    candidate_sizes = [row["candidates"] for row in baseline]

    census = {"annotated": 0, "no_candidate": 0, "no_match": 0}
    confidence: dict[str, int] = {label: 0 for label in rules.confidence_vocabulary}
    correct, wrong = [], []
    for row in baseline:
        census[row["verdict"] if row["verdict"] != "annotated" else "annotated"] += 1
        confidence[row["confidence"]] = confidence.get(row["confidence"], 0) + 1
        if row["verdict"] == "annotated":
            (correct if row["correct"] else wrong).append(row["query"])

    margins = [row["margin_to_other_compound"] for row in baseline
               if row["margin_to_other_compound"] is not None]
    ambiguous = sorted(row["query"] for row in baseline
                       if row["margin_to_other_compound"] is not None
                       and row["margin_to_other_compound"] <= rules.ambiguity_margin)

    # ---- R4a: the tolerance grid, as the set of answers that move
    movement = {}
    for tolerance in rules.tolerance_grid:
        if tolerance == rules.peak_tolerance:
            movement[f"{tolerance}"] = []
            continue
        moved = []
        for row in corpus.annotations(tolerance=tolerance):
            base = by_query[row["query"]]
            if (row["top_match"], row["verdict"]) != (base["top_match"],
                                                      base["verdict"]):
                moved.append(row["query"])
        movement[f"{tolerance}"] = sorted(moved)

    # ---- R4b: the adduct misassignment, as a census
    centres = {query["accession"]:
               query["precursor_mz"] - rules.sodium_mass + rules.proton_mass
               for query in corpus.queries}
    shifted = corpus.annotations(centres=centres)
    adduct = {"no_candidate": 0, "lost": 0, "survives": 0}
    for row in shifted:
        if row["verdict"] in ("no_candidate", "no_match"):
            adduct["no_candidate"] += 1
        elif row["correct"]:
            adduct["survives"] += 1
        else:
            adduct["lost"] += 1

    gold = {
        "environment": {
            "matchms_version": corpus.matchms_version(),
            "similarity": "CosineGreedy",
            "peak_tolerance": rules.peak_tolerance,
            "precursor_window": rules.precursor_window,
        },
        "corpus_census": {
            "graded_queries": len(corpus.queries),
            "library_spectra": len(corpus.library),
            "library_compounds": len({row["compound_block"]
                                      for row in corpus.library}),
            "candidate_set_median": statistics.median(candidate_sizes),
            "candidate_set_max": max(candidate_sizes),
        },
        "retrieval_census": census,
        "correct_queries": sorted(correct),
        "misannotated_queries": sorted(wrong),
        "top_matches": {row["query"]: row["top_match"] for row in baseline},
        "confidence_census": {k: v for k, v in sorted(confidence.items()) if v},
        "margin_summary": {
            "queries_with_a_competing_compound": len(margins),
            "median": round(statistics.median(margins), PRECISION) if margins else None,
            "within_the_declared_margin": len(ambiguous),
        },
        "ambiguous_queries": ambiguous,
        "tolerance_movement": movement,
        "adduct_misassignment_census": adduct,
    }
    if not correct:
        raise DereplicationError(
            "no query retrieves its own compound; the split or the library is "
            "wrong and nothing above R1 would be gradeable")
    return gold


def main(argv: list[str] | None = None) -> int:
    import sys
    campaign = pathlib.Path((argv or sys.argv[1:])[0])
    gold = analyse(campaign)
    print(json.dumps({k: v for k, v in gold.items()
                      if k not in ("top_matches", "correct_queries",
                                   "misannotated_queries", "ambiguous_queries",
                                   "tolerance_movement")}, indent=1))
    print("correct", len(gold["correct_queries"]),
          "wrong", len(gold["misannotated_queries"]),
          "ambiguous", len(gold["ambiguous_queries"]))
    print("tolerance movement:",
          {k: len(v) for k, v in gold["tolerance_movement"].items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
