"""Everything the BGC-detection ladder grades, from the sandbox plus the image.

One antiSMASH run per declared deletion variant, including the identity control.
Each is about 75 seconds over the complete chromosome.
"""

from __future__ import annotations

import pathlib

import yaml

import json

from npbench_c.templates.bgc_detection.core import (
    RULES,
    Detection,
    DetectionError,
)


def coordinate_convention(campaign: pathlib.Path) -> str:
    """The one convention both sides are converted into before comparison."""
    spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
    return spec["detection"]["coordinate_convention"]


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    detection = Detection.load(campaign, params["genbank_input"],
                               params["curated_loci_input"])
    rules = detection.rules

    fingerprint = detection.rule_fingerprint()
    if not fingerprint.startswith(rules.rule_fingerprint):
        raise DetectionError(
            f"the detection rule files hash to {fingerprint} and the campaign "
            f"declares {rules.rule_fingerprint}. A region set is a function of "
            "the rules that produced it, so this is not the same campaign.")

    baseline = detection.regions()
    deletions = {}
    for name, genes in sorted(rules.deletion_variants.items()):
        deletions[name] = detection.deletion_outcome(baseline, genes)

    return {
        "accession": detection.accession,
        "record_length": detection.record_length,
        "antismash_version": detection.antismash_version(),
        "detection_rule_fingerprint": fingerprint,
        "coordinate_convention": coordinate_convention(campaign),
        "search_flags": list(rules.search_flags),
        **detection.region_census(baseline),
        "region_span": detection.span_summary(baseline),
        "region_table": detection.region_table(baseline),
        **detection.reconciliation(baseline),
        "deletion_outcomes": deletions,
    }
