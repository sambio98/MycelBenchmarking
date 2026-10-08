"""Deterministic stub systems for the BGC-detection ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each runs antiSMASH over the complete chromosome.

  reader       clears R1: the region inventory and which rules produced it
  boundaries   clears R2: the region table in both coordinate conventions
  reconcile    clears R3: the reconciliation against the curated boundaries
  complete     clears R4
  offbyone     COMPETENCE probe, tooled: runs the tool correctly and reports the
               1-based inclusive coordinates without the conversion, so every
               start is one too low. Clears R1 and fails R2 -- the search is
               right, every region is right, and one convention is off by one.
  noncompute   LEAKAGE ablation: never runs the tool and never opens the record.

`offbyone` is the probe this campaign needs. The visible work -- running antiSMASH
over 8.67 Mb and projecting 29 regions -- is correct, and what is wrong is a
single conversion the campaign publishes. The symptom is a report that looks
complete and is wrong by one base in half its columns, which no amount of staring
at the tool output would reveal.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.bgc_detection import analyse as bgc
from npbench_c.templates.bgc_detection.core import Detection

LEVELS = ("reader", "boundaries", "reconcile", "complete", "offbyone",
          "noncompute")
REPORT = "detection_report.json"


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
            "accession": "AL645882.2", "record_length": 8667507,
            "coordinate_convention": "0-based half-open",
            "search_flags": ["--minimal"],
            "antismash_version": "antiSMASH 8.0.4",
            "detection_rule_fingerprint": "0" * 16,
            "regions": 25, "distinct_products": 20, "product_census": {},
            "contig_edge_regions": 0, "region_bases": 1000000,
            "fraction_of_record": 0.1,
            "region_table": [], "region_span": {},
            "loci_total": 15, "verdict_census": {}, "per_locus": [],
            "regions_matching_a_locus": 14, "regions_without_a_locus": 11,
            "regions_holding_several_loci": [],
            "jaccard_median": 0.3, "jaccard_min": 0.0, "jaccard_max": 0.9,
            "deletion_outcomes": {},
        })

    params = bgc.load_params(sandbox)
    detection = Detection.load(sandbox, params["genbank_input"],
                               params["curated_loci_input"])
    rules = detection.rules
    baseline = detection.regions()

    report = {
        "accession": detection.accession,
        "record_length": detection.record_length,
        "antismash_version": detection.antismash_version(),
        "detection_rule_fingerprint": detection.rule_fingerprint(),
        "coordinate_convention": bgc.coordinate_convention(sandbox),
        "search_flags": list(rules.search_flags),
        **detection.region_census(baseline),
    }
    if args.level == "reader":
        return write(report)

    table = detection.region_table(baseline)
    if args.level == "offbyone":
        # The tool ran, the regions are right; the 1-based inclusive form is
        # reported without the +1 on the start.
        table = [{**row, "start_one_based": row["start"]} for row in table]
    report["region_table"] = table
    report["region_span"] = detection.span_summary(baseline)
    if args.level in ("boundaries", "offbyone"):
        return write(report)

    report.update(detection.reconciliation(baseline))
    if args.level == "reconcile":
        return write(report)

    report["deletion_outcomes"] = {
        name: detection.deletion_outcome(baseline, genes)
        for name, genes in sorted(rules.deletion_variants.items())}
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
