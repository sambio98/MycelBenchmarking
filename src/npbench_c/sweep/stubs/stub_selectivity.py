"""Deterministic stub systems for the selectivity ladder.

Harness fixtures. Sandbox-only.

  reader      clears R1: the activity and assay census, including the
              confidence_score distribution that lives on the assay
  admitted    clears R2: admission under the declared rules, with the rejection
              breakdown that must partition the rejected set
  ratios      clears R3: the selectivity verdicts
  complete    clears R4
  inverted    COMPETENCE probe, tooled: admits exactly the right activities and
              aggregates them correctly, then forms the ratio the other way
              round. Clears R2 and fails R3 -- every magnitude right and every
              verdict backwards, which is the error a declared ratio direction
              exists to catch.
  noncompute  LEAKAGE ablation: never opens the input.
"""

from __future__ import annotations

import argparse
import collections
import json
import pathlib

from npbench_c.templates.chembl_selectivity import analyse as sel
from npbench_c.templates.chembl_selectivity.core import Admission, Selectivity

LEVELS = ("reader", "admitted", "ratios", "complete", "inverted", "noncompute")
REPORT = "selectivity_report.json"


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
            "release": "ChEMBL_37", "focal_target": "?", "reference_target": "?",
            "activities_total": 1000, "activities_of_declared_type": 800,
            "distinct_assays": 50, "assay_confidence_scores": {},
            "admitted_by_target": {}, "rejected_by_target": {},
            "aggregated_pairs": 0, "aggregated_by_target": {},
            "eligible_molecules": 100, "focal_selective_counts": {},
            "reference_selective_counts": {},
            "focal_selective_at_highest_threshold": [],
            "most_selective_molecule": "CHEMBL1",
            "most_selective_ratio": 10.0,
            "admission_variants": {}, "ratio_extremes": {},
        })

    params = sel.load_params(sandbox)
    selectivity = Selectivity.load(sandbox, params["activities_input"],
                                   params["focal_target"],
                                   params["reference_target"])
    baseline = Admission()
    pair = (selectivity.focal_target, selectivity.reference_target)
    on_pair = [a for a in selectivity.of_type() if a["target_chembl_id"] in pair]
    assay_ids = sorted({a["assay_chembl_id"] for a in on_pair
                        if a.get("assay_chembl_id")})
    scores = collections.Counter(
        selectivity.assays.get(i, {}).get("confidence_score") for i in assay_ids)

    report = {
        "release": selectivity.release,
        "focal_target": selectivity.focal_target,
        "reference_target": selectivity.reference_target,
        "activities_total": len(selectivity.activities),
        "activities_of_declared_type": len(on_pair),
        "distinct_assays": len(assay_ids),
        "assay_confidence_scores": {str(k): v for k, v in sorted(
            scores.items(), key=lambda kv: (kv[0] is None, kv[0]))},
    }
    if args.level == "reader":
        return write(report)

    census = selectivity.admission_census(baseline)
    values = selectivity.aggregated(baseline)
    report.update({
        **census,
        "aggregated_pairs": len(values),
        "aggregated_by_target": {
            t: sum(1 for _, tg in values if tg == t) for t in sorted(pair)},
    })
    ratios = selectivity.ratios(baseline)
    report["eligible_molecules"] = len(ratios)
    if args.level == "admitted":
        return write(report)

    if args.level == "inverted":
        # Every magnitude right, every verdict backwards: the ratio is formed
        # focal-over-reference where the declared direction is the other way.
        flipped = {m: round(1.0 / r, 6) for m, r in ratios.items() if r > 0}
        report.update(selectivity.verdicts(flipped))
        return write(report)

    report.update(selectivity.verdicts(ratios))
    if args.level == "ratios":
        return write(report)

    report["admission_variants"] = {}
    for name, spec in sorted((params.get("admission_variants") or {}).items()):
        admission = Admission.from_spec(spec or {})
        block_census = selectivity.admission_census(admission)
        block = selectivity.verdicts(selectivity.ratios(admission))
        report["admission_variants"][name] = {
            "aggregator": admission.aggregator,
            "admitted_total": sum(block_census["admitted_by_target"].values()),
            "eligible_molecules": block["eligible_molecules"],
            "focal_selective_counts": block["focal_selective_counts"],
            "reference_selective_counts": block["reference_selective_counts"],
        }
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
