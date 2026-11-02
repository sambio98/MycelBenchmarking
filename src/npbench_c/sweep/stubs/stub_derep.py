"""Deterministic stub systems for the MS2 dereplication ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each one calls matchms through the tool prefix.

  reader      clears R1: what ran, and the precursor filter over the library
  retrieve    clears R2: what the search retrieved, per query
  commit      clears R3: the declared confidence class and the ambiguity
  complete    clears R4
  peaksfirst  COMPETENCE probe, tooled: ranks candidates by MATCHED PEAKS before
              score. Clears R1 and fails R2.
  scorefloor  COMPETENCE probe, tooled: reads the declared confidence rule as a
              score threshold and drops its matched-peak floor. Clears R1 and R2,
              fails R3.
  noncompute  LEAKAGE ablation: never scores a spectrum and never opens an input.

Two probes, because this campaign has two independent places to be plausibly
wrong and the measurements put them on different rungs.

**`peaksfirst`.** "More matched peaks is a better match" is a defensible reading
of a spectral comparison and it is not the declared rule, which orders by cosine
score first. Measured on the shipped hundred: 65 of the 98 annotated queries get
a different top match, and only 2 of them change from correct to incorrect. So
the retrieval census is untouched, R1's candidate census is untouched -- the
precursor filter is upstream of any ranking -- and what fails is the per-query
answer at R2, which is where the search is.

**`scorefloor`.** The declared confidence rule is a pair of floors, a score AND a
matched-peak count, and the first row whose BOTH are met wins. Reading it as a
score threshold alone is the single-condition slip. Measured: 17 of the 98
annotated queries move class, 15 from `tentative` and 2 from `unsupported`, all
of them upward to `confirmed` -- the direction that makes a system look more
certain than its own evidence supports. R1 and R2 are untouched because the class
is downstream of the retrieval, and two of R3's four components fail: the margin
summary and the ambiguous set do not depend on the confidence rule.

Neither probe is the campaign's third plausible slip. Dropping the declared
ranking's tie-break keys -- ordering by score alone, which leaves equal scores to
whatever order matchms returned -- moves nothing at all on this corpus: 8 queries
have tied candidate scores and none of them at the top. That measurement is
recorded on `core.declared_ranking`, which is where the rule lives, and the rule
stays: it is what makes the answer independent of the tool's iteration order on a
corpus where the tie does reach the top. It is not a fixture, because a probe
that fails nothing clears every rung and measures the harness rather than the
system.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import pathlib

from npbench_c.templates.ms2_dereplication import analyse as derep
from npbench_c.templates.ms2_dereplication.core import Corpus

LEVELS = ("reader", "retrieve", "commit", "complete", "peaksfirst", "scorefloor",
          "noncompute")
REPORT = "dereplication_report.json"

#: The probe's wrong order: matched peaks down, then score down, then accession.
#: Internally consistent -- it is a total order and it is deterministic -- and
#: not the declared one.
def _peaks_first(row) -> tuple:
    return (-row[2], -row[1], row[0])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    out = pathlib.Path(args.submission)
    out.mkdir(parents=True, exist_ok=True)
    target = out / REPORT

    def write(payload: dict) -> int:
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return 0

    if args.level == "noncompute":
        return write({
            "environment": {}, "corpus_census": {}, "candidate_census": {},
            "retrieval_census": {}, "top_matches": {}, "confidence_by_query": {},
            "confidence_census": {}, "margin_summary": {},
            "ambiguous_queries": [], "tolerance_movement": {},
            "adduct_misassignment_census": {},
        })

    params = derep.load_params(sandbox)
    corpus = Corpus.load(sandbox, params["queries_input"], params["library_input"])

    # R1. The candidate census is the precursor filter alone, so this level
    # never scores a spectrum -- only the version query touches the tool.
    report = {
        "environment": derep.environment(corpus),
        "corpus_census": derep.corpus_census(corpus),
        "candidate_census": derep.candidate_census(
            [len(corpus.candidates(query)) for query in corpus.queries]),
    }
    if args.level == "reader":
        return write(report)

    # R2.
    ranking = _peaks_first if args.level == "peaksfirst" else None
    rows = corpus.annotations(ranking=ranking)
    report.update(derep.retrieval(rows))
    if args.level in ("retrieve", "peaksfirst"):
        return write(report)

    # R3. The probe keeps the declared score floors and drops the matched-peak
    # floor, which is the whole of the slip: the rule is otherwise untouched and
    # every other number in the report still closes.
    rules = corpus.rules
    if args.level == "scorefloor":
        loose = dataclasses.replace(rules, confidence_rule=tuple(
            (floor, 0, label) for floor, _, label in rules.confidence_rule))
        classified = dataclasses.replace(corpus, rules=loose).annotations()
        report.update(derep.commitments(classified, loose))
        return write(report)
    report.update(derep.commitments(rows, rules))
    if args.level == "commit":
        return write(report)

    # R4.
    report["tolerance_movement"] = derep.tolerance_movement(corpus, rows)
    report["adduct_misassignment_census"] = derep.adduct_census(corpus, rows)
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
