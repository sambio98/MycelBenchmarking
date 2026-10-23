"""Everything the selectivity ladder grades, from the sandbox alone."""

from __future__ import annotations

import collections
import pathlib

import yaml

from npbench_c.templates.chembl_selectivity.core import (
    Admission,
    Selectivity,
    SelectivityError,
)


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    selectivity = Selectivity.load(campaign, params["activities_input"],
                                   params["focal_target"],
                                   params["reference_target"])
    baseline = Admission()
    of_type = selectivity.of_type()
    pair = (selectivity.focal_target, selectivity.reference_target)
    on_pair = [a for a in of_type if a["target_chembl_id"] in pair]

    assay_ids = sorted({a["assay_chembl_id"] for a in on_pair
                        if a.get("assay_chembl_id")})
    scores = collections.Counter(
        selectivity.assays.get(i, {}).get("confidence_score") for i in assay_ids)

    census = selectivity.admission_census(baseline)
    values = selectivity.aggregated(baseline)
    ratios = selectivity.ratios(baseline)
    verdicts = selectivity.verdicts(ratios)

    variants = {}
    for name, spec in sorted((params.get("admission_variants") or {}).items()):
        admission = Admission.from_spec(spec or {})
        block_census = selectivity.admission_census(admission)
        block_ratios = selectivity.ratios(admission)
        block_verdicts = selectivity.verdicts(block_ratios)
        variants[name] = {
            "aggregator": admission.aggregator,
            "admitted_total": sum(block_census["admitted_by_target"].values()),
            "eligible_molecules": block_verdicts["eligible_molecules"],
            "focal_selective_counts": block_verdicts["focal_selective_counts"],
            "reference_selective_counts":
                block_verdicts["reference_selective_counts"],
        }

    return {
        "release": selectivity.release,
        "focal_target": selectivity.focal_target,
        "reference_target": selectivity.reference_target,
        "activities_total": len(selectivity.activities),
        "activities_of_declared_type": len(on_pair),
        "distinct_assays": len(assay_ids),
        "assay_confidence_scores": {str(k): v for k, v in sorted(
            scores.items(), key=lambda kv: (kv[0] is None, kv[0]))},
        **census,
        "aggregated_pairs": len(values),
        "aggregated_by_target": {
            t: sum(1 for _, target in values if target == t) for t in sorted(pair)},
        **verdicts,
        "admission_variants": variants,
        # Audit material, not graded: see excluded_from_grading.
        "ratio_extremes": {
            "most_focal_selective": sorted(
                ((-r, m) for m, r in ratios.items()))[:5],
            "most_reference_selective": sorted(
                ((r, m) for m, r in ratios.items()))[:5],
        },
    }
