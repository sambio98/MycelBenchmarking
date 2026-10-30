"""Deterministic stub systems for the remote-homology ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each runs the pinned searcher over the shipped corpus.

  reader      clears R1: the corpus, the declared queries, and one search
  recover     clears R2: what the search found of what the resource asserts
  verdict     clears R3: the per-band verdicts
  complete    clears R4
  selfhit     COMPETENCE probe, tooled: runs the search correctly, reports the
              hit counts correctly, and keeps the SELF-HIT in the label join.
              Clears R1 and fails R2.
  noncompute  LEAKAGE ablation: never searches and never opens the corpus.

`selfhit` is this campaign's probe because the query is a member of the corpus it
is searched against. Every query therefore returns itself at 100 per cent
identity, and that self-hit is not a finding about remote homology: it is the
query.

There are two distinguishable places to get that wrong, and the difference
matters enough to state. Keeping the self-hit in the HIT COUNTS inflates every
per-query total by one and shifts the top identity band, which R1 grades.
Keeping it in the LABEL JOIN leaves the counts right and adds one
`same_family` hit and one detected relative -- because a protein shares every one
of its own families with itself -- which R2 grades.

This stub does the second. The first is measured and recorded in
docs/npbench-c.md, with a test pinning the off-by-one, but it is not the fixture:
the harness's own R1 clear-rate gate exists to catch a campaign whose first rung
is unclearable, and a fixture built to fail R1 drags that rate below its floor
and makes the gate report a fixture artefact as a campaign defect. Measured here
before the change: R1 pass@1 of 0.75 against a 0.8 minimum.

So the symptom is an off-by-one in the reasoning layer, internally consistent
everywhere: the recovery rate moves by a fraction of a per cent, the census keeps
its shape, every arithmetic identity in the report still closes, and no number
looks wrong on its own. It clears R1, where the search is, and fails R2, where
the search meets the label.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.remote_homology import analyse as rha
from npbench_c.templates.remote_homology.core import Corpus

LEVELS = ("reader", "recover", "verdict", "complete", "selfhit", "noncompute")
REPORT = "homology_report.json"


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
            "label_resource": "", "corpus_superfamily": "", "sensitivity": 0.0,
            "evalue": 0.0, "max_seqs": 0, "threads": 0, "format_fields": [],
            "mmseqs_version": "", "corpus_sequences": 0,
            "corpus_sha256": "0" * 64, "query_table": {}, "baseline_hits": {},
            "baseline_band_sizes": {}, "recovery": {}, "band_census": {},
            "band_verdicts": {}, "verdict_census": {}, "shuffle_ladder": {},
        })

    params = rha.load_params(sandbox)
    corpus = Corpus.load(sandbox, params["corpus_input"], params["labels_input"])
    rules = corpus.rules
    queries = corpus.queries()
    accessions = sorted(queries.values())
    baseline = corpus.search({a: corpus.sequences[a] for a in accessions})

    report = {
        "label_resource": "InterPro 110.0",
        "corpus_superfamily": rules.corpus_superfamily,
        "sensitivity": rules.sensitivity,
        "evalue": rules.evalue,
        "max_seqs": rules.max_seqs,
        "threads": rules.threads,
        "format_fields": list(rules.format_fields),
        "mmseqs_version": corpus.mmseqs_version(),
        "corpus_sequences": corpus.size,
        "corpus_sha256": corpus.sequences_sha256(),
        "query_table": {
            entry: {"accession": accession,
                    "length": len(corpus.sequences[accession]),
                    "family_members": sum(
                        1 for entries in corpus.families.values() if entry in entries),
                    "family_name": corpus.entry_names.get(entry, "")}
            for entry, accession in sorted(queries.items())},
    }

    report["baseline_hits"] = {
        a: sum(1 for h in baseline if h.query == a and h.target != a)
        for a in accessions}
    report["baseline_band_sizes"] = corpus.baseline_band_sizes(baseline)
    if args.level == "reader":
        return write(report)

    if args.level == "selfhit":
        # The counts are right; the join is not. The query is counted as its own
        # same-family relative, which is true of any protein and tells you
        # nothing, so the detected count and the rate each move by one.
        report["recovery"] = {}
        report["band_census"] = {a: corpus.band_census(baseline, a)
                                 for a in accessions}
        for a in accessions:
            row = dict(corpus.recovery(baseline, a))
            asserted = row["same_family_asserted"]
            row["same_family_detected"] += 1
            row["recovery_rate"] = (round(row["same_family_detected"] / asserted, 6)
                                    if asserted else None)
            report["recovery"][a] = row
        return write(report)

    report["recovery"] = {a: corpus.recovery(baseline, a) for a in accessions}
    census = {a: corpus.band_census(baseline, a) for a in accessions}
    report["band_census"] = census
    if args.level == "recover":
        return write(report)

    report["band_verdicts"] = {
        a: {band: corpus.verdict(cells[band]) for band in rules.band_labels}
        for a, cells in sorted(census.items())}
    counts = {verdict: 0 for verdict in rules.verdicts}
    for cells in census.values():
        for band in rules.band_labels:
            counts[corpus.verdict(cells[band])] += 1
    report["verdict_census"] = counts
    if args.level == "verdict":
        return write(report)

    report["shuffle_ladder"] = corpus.ladder(baseline)
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
