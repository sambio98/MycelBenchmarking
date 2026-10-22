"""Measure a dereplication submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies in `reference/dereplication_rules.json` -- the confidence classes,
the tolerance grid -- from this template's own closed column lists, or from the
shipped queries. The query accessions are the one thing the measurer cannot
declare in advance, and it reads them from `inputs/queries.json.gz` rather than
out of the submission, so a submission that answered for a query set of its own
scores zero instead of being read on its own terms.

Usage:
  python -m npbench_c.templates.ms2_dereplication.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.ms2_dereplication.analyse import (
    ADDUCT_OUTCOMES,
    RETRIEVAL_VERDICTS,
    load_params,
)
from npbench_c.templates.ms2_dereplication.core import Corpus

REPORT = "dereplication_report.json"

ENVIRONMENT_KEYS = ("matchms_version", "similarity", "peak_tolerance",
                    "precursor_window")
CANDIDATE_KEYS = ("median", "max", "total")
MARGIN_KEYS = ("queries_with_a_competing_compound", "median",
               "within_the_declared_margin")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _opt_float(value: object) -> float | None:
    """A median that is legitimately null in gold, where no competitor exists."""
    return None if value is None else _float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _query_set(value: object) -> list[str]:
    """A set answer, normalised to gold's order and nothing more.

    Sorted because the claim is a set and gold is sorted; NOT deduplicated,
    because a submission that names the same query twice has made a different
    claim from the one gold makes.
    """
    return sorted(str(v) for v in value) if isinstance(value, list) else []


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["queries_input"], params["library_input"])
    rules = corpus.rules
    accessions = sorted(query["accession"] for query in corpus.queries)
    classes = sorted(rules.confidence_vocabulary)
    tolerances = [f"{tolerance}" for tolerance in rules.tolerance_grid]

    path = pathlib.Path(submission_dir) / REPORT
    out: dict = {"present": {"report": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["report"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["report"] = 0
        return out

    # ---- R1: what ran the comparison, and the precursor filter it applied
    environment = sub.get("environment")
    environment = environment if isinstance(environment, dict) else {}
    out["environment"] = {
        "matchms_version": _str(environment.get("matchms_version")),
        "similarity": _str(environment.get("similarity")),
        "peak_tolerance": _float(environment.get("peak_tolerance")),
        "precursor_window": _float(environment.get("precursor_window")),
    }
    candidates = sub.get("candidate_census")
    candidates = candidates if isinstance(candidates, dict) else {}
    out["candidate_census"] = {
        "median": _float(candidates.get("median")),
        "max": _int(candidates.get("max")),
        "total": _int(candidates.get("total")),
    }

    # ---- R2: what the search retrieved
    out["retrieval_census"] = _int_map(sub.get("retrieval_census"),
                                       sorted(RETRIEVAL_VERDICTS))
    matches = sub.get("top_matches")
    matches = matches if isinstance(matches, dict) else {}
    out["top_matches"] = {accession: _str(matches.get(accession))
                          for accession in accessions}
    out["correct_queries"] = _query_set(sub.get("correct_queries"))

    # ---- R3: the declared class per query, and the ambiguity admitted to
    classified = sub.get("confidence_by_query")
    classified = classified if isinstance(classified, dict) else {}
    out["confidence_by_query"] = {accession: _str(classified.get(accession))
                                  for accession in accessions}
    out["confidence_census"] = _int_map(sub.get("confidence_census"), classes)
    margin = sub.get("margin_summary")
    margin = margin if isinstance(margin, dict) else {}
    out["margin_summary"] = {
        "queries_with_a_competing_compound":
            _int(margin.get("queries_with_a_competing_compound")),
        "median": _opt_float(margin.get("median")),
        "within_the_declared_margin": _int(margin.get("within_the_declared_margin")),
    }
    out["ambiguous_queries"] = _query_set(sub.get("ambiguous_queries"))

    # ---- R4: the tolerance grid and the misassigned adduct
    movement = sub.get("tolerance_movement")
    movement = movement if isinstance(movement, dict) else {}
    out["tolerance_movement"] = {tolerance: _query_set(movement.get(tolerance))
                                 for tolerance in tolerances}
    out["adduct_misassignment_census"] = _int_map(
        sub.get("adduct_misassignment_census"), sorted(ADDUCT_OUTCOMES))
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>",
              file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
