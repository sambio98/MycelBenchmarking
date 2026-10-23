"""Everything the gene-calling ladder grades, from the sandbox plus the image.

The calls that are not pure computation are the gene-calling runs: one per
declared (mode, fragment count) cell, which is eight. Each is about fifteen
seconds on this chromosome, so the whole build is a couple of minutes and needs
no reference database at all -- the first tool-running campaign here with no
external resource to pin.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.gene_calling.core import Calling


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    calling = Calling.load(campaign, params["contig_input"])
    rules = calling.rules

    whole = calling.calls(rules.baseline_mode, 1)
    fragmented = calling.calls(rules.baseline_mode, rules.baseline_fragments)
    cuts = rules.breakpoints(calling.length, rules.baseline_fragments)

    return {
        "accession": calling.accession,
        "baseline_mode": rules.baseline_mode,
        "baseline_fragments": rules.baseline_fragments,
        "locality_window": rules.locality_window,
        "fixed_flags": list(rules.fixed_flags),
        "prodigal_version": calling.prodigal_version(),

        # ---- R1: the contig, and what one run of the tool over it says
        "contig_length": calling.length,
        "contig_sha256": calling.sequence_sha256(),
        "gc_fraction": calling.gc_fraction(),
        "call_census": calling.call_census(whole),
        "partial_census": calling.partial_census(whole),
        "start_type_census": calling.start_type_census(whole),
        "breakpoint_census": {
            str(n): len(rules.breakpoints(calling.length, n))
            for n in rules.fragment_counts},

        # ---- R2: the calls themselves, in parent coordinates
        "call_table": calling.call_table(whole),
        "call_span": calling.span(whole),

        # ---- R3: the same contig cut up, against the whole one
        **{f"fragmented_{k}": v for k, v in
           calling.reconcile(whole, fragmented, cuts).items()},
        "fragmented_partial_census": calling.partial_census(fragmented),
        "fragmented_call_census": calling.call_census(fragmented),

        # ---- R4: the declared grid, and the flag against the assembly
        "fragmentation_grid": calling.grid(),
        "mode_difference": calling.mode_difference(),
    }
