"""Deterministic stub systems for the gene-calling ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each runs the pinned caller over the contig.

  reader      clears R1: the contig's own census and one run of the caller over it
  caller      clears R2: the call table in parent coordinates
  reconcile   clears R3: the fragmented contig against the whole one
  complete    clears R4
  unlifted    COMPETENCE probe, tooled: runs the caller correctly on the declared
              fragments and reports FRAGMENT coordinates without adding each
              fragment's offset. Clears R1 and fails R3.
  noncompute  LEAKAGE ablation: never runs the caller and never opens the contig.

`unlifted` is the probe this campaign needs, because its symptom is a report that
is internally consistent and wrong. The coordinate lift is one addition per call;
skip it and the fragmented run's coordinates are all small numbers, so almost
every whole-contig call looks lost and almost every fragment call looks new. The
counts stay in the right order of magnitude, the tables are the right shape, and
nothing says a convention was missed. It clears R1 and R2 -- which are about the
WHOLE contig and need no lift -- and fails at R3, which is where the two
coordinate systems first have to meet.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import pathlib

from npbench_c.templates.gene_calling import analyse as gca
from npbench_c.templates.gene_calling.core import Calling

LEVELS = ("reader", "caller", "reconcile", "complete", "unlifted", "noncompute")
REPORT = "calling_report.json"


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
            "accession": "", "baseline_mode": "", "baseline_fragments": 0,
            "locality_window": 0, "fixed_flags": [], "prodigal_version": "",
            "contig_length": 0, "contig_sha256": "0" * 64, "gc_fraction": 0.5,
            "call_census": {}, "partial_census": {}, "start_type_census": {},
            "breakpoint_census": {},
            "call_table": [], "call_span": {},
            "fragmentation_summary": {}, "fragmented_loss_verdicts": {},
            "fragmented_lost_table": [], "fragmented_partial_census": {},
            "fragmented_gained_partial_census": {},
            "fragmented_call_census": {},
            "fragmentation_grid": {}, "mode_difference": {},
        })

    params = gca.load_params(sandbox)
    calling = Calling.load(sandbox, params["contig_input"])
    rules = calling.rules

    whole = calling.calls(rules.baseline_mode, 1)
    report = {
        "accession": calling.accession,
        "baseline_mode": rules.baseline_mode,
        "baseline_fragments": rules.baseline_fragments,
        "locality_window": rules.locality_window,
        "fixed_flags": list(rules.fixed_flags),
        "prodigal_version": calling.prodigal_version(),
        "contig_length": calling.length,
        "contig_sha256": calling.sequence_sha256(),
        "gc_fraction": calling.gc_fraction(),
        "call_census": calling.call_census(whole),
        "partial_census": calling.partial_census(whole),
        "start_type_census": calling.start_type_census(whole),
        "breakpoint_census": {
            str(n): len(rules.breakpoints(calling.length, n))
            for n in rules.fragment_counts},
    }
    if args.level == "reader":
        return write(report)

    report["call_table"] = calling.call_table(whole)
    report["call_span"] = calling.span(whole)
    if args.level == "caller":
        return write(report)

    fragmented = calling.calls(rules.baseline_mode, rules.baseline_fragments)
    cuts = rules.breakpoints(calling.length, rules.baseline_fragments)

    if args.level == "unlifted":
        # Every fragment call moved back to its own fragment's coordinates: the
        # lift undone, one subtraction per call, nothing else touched.
        offsets = calling.offsets(rules.baseline_fragments)
        unlifted = tuple(
            dataclasses.replace(call,
                                start=call.start - offsets[call.fragment],
                                end=call.end - offsets[call.fragment])
            for call in fragmented)
        block = calling.reconcile(whole, unlifted, cuts)
        report.update({f"fragmented_{k}": v for k, v in block.items()})
        report["fragmentation_summary"] = {
            k: block[k] for k in ("baseline_calls", "fragmented_calls", "retained",
                                  "lost", "gained", "retained_but_truncated")}
        report["fragmented_partial_census"] = calling.partial_census(unlifted)
        report["fragmented_call_census"] = calling.call_census(unlifted)
        return write(report)

    block = calling.reconcile(whole, fragmented, cuts)
    report.update({f"fragmented_{k}": v for k, v in block.items()})
    report["fragmentation_summary"] = {
        k: block[k] for k in ("baseline_calls", "fragmented_calls", "retained",
                              "lost", "gained", "retained_but_truncated")}
    report["fragmented_partial_census"] = calling.partial_census(fragmented)
    report["fragmented_call_census"] = calling.call_census(fragmented)
    if args.level == "reconcile":
        return write(report)

    report["fragmentation_grid"] = calling.grid()
    report["mode_difference"] = calling.mode_difference()
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
