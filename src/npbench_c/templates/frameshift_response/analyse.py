"""Everything the frameshift-response ladder grades, from the sandbox plus the image.

The calls that are not pure computation are the caller runs: the two whole-contig
baselines, one per declared target under the baseline mode, four more for the
offset sweep, and one per target in the noise arm. About twenty runs at fifteen
seconds each on this chromosome, and no reference database at all.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.frameshift_response.core import Frameshift


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    shift = Frameshift.load(campaign, params["contig_input"])
    rules = shift.rules

    baseline = shift.calls(rules.baseline_mode)
    targets = shift.targets(baseline)
    responses = shift.responses(baseline, targets)

    return {
        "accession": shift.accession,
        "baseline_mode": rules.baseline_mode,
        "noise_mode": rules.noise_mode,
        "indel_fraction": rules.indel_fraction,
        "fixed_flags": list(rules.fixed_flags),
        "prodigal_version": shift.prodigal_version(),

        # ---- R1: the contig, the baseline run, and the declared targets
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

        # ---- R2: what the caller did to each broken gene
        "response_table": responses,

        # ---- R3: the verdicts
        "verdict_census": shift.verdict_census(responses),
        "verdict_by_length": shift.verdict_by_length(responses),
        "verdict_is_monotone_in_length":
            shift.verdict_is_monotone_in_length(responses),

        # ---- R4: the indel walked across a gene, and the mode that drowns it
        "offset_sweep": shift.offset_sweep(baseline, targets),
        "noise_arm": shift.noise_arm(targets),
    }
