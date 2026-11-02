"""Deterministic stub systems for the frameshift-response ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each runs the pinned caller over the contig and over mutated
copies of it.

  reader      clears R1: the contig, one baseline run, and the declared targets
  respond     clears R2: what the caller did to each broken gene
  verdict     clears R3: the verdict census and the per-length list
  complete    clears R4
  unlifted    COMPETENCE probe, tooled: deletes the right base, runs the caller
              correctly, and reports the mutated run's coordinates WITHOUT
              lifting them back past the deletion. Clears R1 and fails R2.
  noncompute  LEAKAGE ablation: never runs the caller and never opens the contig.

`unlifted` is this campaign's probe, and it is the cheapest way to get the whole
experiment wrong while looking right. The lift is one addition per call -- a
mutated coordinate at or after the deleted position stands for a parent
coordinate one higher -- and skipping it moves every call downstream of the
deletion by exactly one base. On this chromosome that is some seven thousand
calls displaced by a single base, so almost nothing matches the baseline.

What that does to the report is worth stating exactly, because it is not what one
would guess. The VERDICT survives: overlap with the target's own span is robust
to a one-base displacement, so a gene the caller split in two is still reported
as split in two. What breaks is everything that compares two runs.
`other_calls_lost` and `other_calls_gained` go from zero to about seven thousand
each, the new ORFs come out one base short, a call that was already inside the
target's span counts as a third new call because it moved too, and the coverage
column runs past the gene's end. Measured on the 22 kb target: the truth is two
new calls of 11233 and 11187 bases with no collateral at all, and the unlifted
run reports three new calls -- 11232, 11187 and 11013 -- with 4571 lost and 4570
gained.

Across the whole ladder the unlifted run gets 19 of the 21 verdicts right and
flips two -- the 240 bp gene reads truncated instead of abolished, the 249 bp one
split instead of truncated, both because the displacement turns a call already
inside the span into a spurious new one -- while its collateral columns run from
1,539 to 7,745 against gold's 0 on nineteen rows.

So the probe clears R1, which is about the contig and the untouched baseline and
needs no lift; fails R2 in every row; and fails R3 by two rows out of
twenty-one. That ratio is what makes it worth having: a system can be right about
what happened to nearly every gene and still not have run the experiment.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.frameshift_response import analyse as fsa
from npbench_c.templates.frameshift_response.core import Frameshift

LEVELS = ("reader", "respond", "verdict", "complete", "unlifted", "noncompute")
REPORT = "frameshift_report.json"


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
            "accession": "", "baseline_mode": "", "noise_mode": "",
            "indel_fraction": 0.0, "fixed_flags": [], "prodigal_version": "",
            "contig_length": 0, "contig_sha256": "0" * 64, "baseline_calls": 0,
            "target_table": {}, "target_length_summary": {},
            "response_table": {}, "verdict_census": {},
            "verdict_by_length": {}, "verdict_is_monotone_in_length": False,
            "offset_sweep": {}, "noise_arm": {},
        })

    params = fsa.load_params(sandbox)
    shift = Frameshift.load(sandbox, params["contig_input"])
    rules = shift.rules

    baseline = shift.calls(rules.baseline_mode)
    targets = shift.targets(baseline)
    report = {
        "accession": shift.accession,
        "baseline_mode": rules.baseline_mode,
        "noise_mode": rules.noise_mode,
        "indel_fraction": rules.indel_fraction,
        "fixed_flags": list(rules.fixed_flags),
        "prodigal_version": shift.prodigal_version(),
        "contig_length": shift.length,
        "contig_sha256": shift.sequence_sha256(),
        "baseline_calls": len(baseline),
        "target_table": {
            label: {"start": t.start, "end": t.end, "strand": t.strand,
                    "length": t.length,
                    "deleted_at": shift.indel_position(t)}
            for label, t in sorted(targets.items(), key=lambda kv: int(kv[0]))},
        "target_length_summary": {
            "targets": len(targets),
            "shortest": min(t.length for t in targets.values()),
            "longest": max(t.length for t in targets.values()),
        },
    }
    if args.level == "reader":
        return write(report)

    if args.level == "unlifted":
        # The lift undone: one subtraction per call, for the calls at or after
        # the deleted base. Nothing else is touched -- the right base is deleted,
        # the caller runs correctly, and the classification is the campaign's own.
        import dataclasses

        responses = {}
        for label, t in sorted(targets.items(), key=lambda kv: int(kv[0])):
            position = shift.indel_position(t)
            mutated = tuple(sorted(
                (dataclasses.replace(
                    call,
                    start=call.start - (1 if call.start > position else 0),
                    end=call.end - (1 if call.end > position else 0))
                 for call in shift.calls(rules.baseline_mode, position)),
                key=lambda c: c.key))
            responses[label] = {"deleted_at": position,
                                **shift.response(baseline, t, mutated)}
        report["response_table"] = responses
        return write(report)

    responses = shift.responses(baseline, targets)
    report["response_table"] = responses
    if args.level == "respond":
        return write(report)

    report["verdict_census"] = shift.verdict_census(responses)
    report["verdict_by_length"] = shift.verdict_by_length(responses)
    report["verdict_is_monotone_in_length"] = \
        shift.verdict_is_monotone_in_length(responses)
    if args.level == "verdict":
        return write(report)

    report["offset_sweep"] = shift.offset_sweep(baseline, targets)
    report["noise_arm"] = shift.noise_arm(targets)
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
