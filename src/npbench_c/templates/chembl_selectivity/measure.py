"""Measure a selectivity submission. Deterministic, agent-independent.

Nothing here reads gold: the target ids, the thresholds and the variant names all
come from the campaign's own declared vocabularies.

Usage:
  python -m npbench_c.templates.chembl_selectivity.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.chembl_selectivity.analyse import load_params
from npbench_c.templates.chembl_selectivity.core import Rules

REPORT = "selectivity_report.json"


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped: gold omits them."""
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    params = load_params(campaign)
    targets = tuple(sorted((params["focal_target"], params["reference_target"])))
    thresholds = tuple(str(t) for t in rules.thresholds)
    variants = tuple(sorted(params.get("admission_variants") or {}))

    path = pathlib.Path(submission_dir) / REPORT
    out: dict = {"present": {"report": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["report"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["report"] = 0
        return out

    # ---- R1: the activity and assay census
    for key in ("release", "focal_target", "reference_target"):
        out[key] = _str(sub.get(key))
    for key in ("activities_total", "activities_of_declared_type",
                "distinct_assays"):
        out[key] = _int(sub.get(key))
    out["assay_confidence_scores"] = _free_int_map(sub.get("assay_confidence_scores"))

    # ---- R2: admission and aggregation
    out["admitted_by_target"] = _int_map(sub.get("admitted_by_target"), targets)
    rejected = sub.get("rejected_by_target")
    rejected = rejected if isinstance(rejected, dict) else {}
    out["rejected_by_target"] = {
        t: _int_map(rejected.get(t), rules.rejection_order) for t in targets}
    out["aggregated_pairs"] = _int(sub.get("aggregated_pairs"))
    out["aggregated_by_target"] = _int_map(sub.get("aggregated_by_target"), targets)
    out["eligible_molecules"] = _int(sub.get("eligible_molecules"))

    # ---- R3: the selectivity verdict
    out["focal_selective_counts"] = _int_map(sub.get("focal_selective_counts"),
                                             thresholds)
    out["reference_selective_counts"] = _int_map(
        sub.get("reference_selective_counts"), thresholds)
    committed = sub.get("focal_selective_at_highest_threshold")
    out["focal_selective_at_highest_threshold"] = (
        sorted({str(m) for m in committed}) if isinstance(committed, list) else [])
    out["most_selective_molecule"] = _str(sub.get("most_selective_molecule"))
    out["most_selective_ratio"] = _float(sub.get("most_selective_ratio"))

    # ---- R4: the admission counterfactual
    declared = sub.get("admission_variants")
    declared = declared if isinstance(declared, dict) else {}
    admitted, eligible, focal, reference = {}, {}, {}, {}
    for name in variants:
        block = declared.get(name)
        block = block if isinstance(block, dict) else {}
        admitted[name] = _int(block.get("admitted_total"))
        eligible[name] = _int(block.get("eligible_molecules"))
        focal[name] = _int_map(block.get("focal_selective_counts"), thresholds)
        reference[name] = _int_map(block.get("reference_selective_counts"),
                                   thresholds)
    out["variant_admitted"] = admitted
    out["variant_eligible"] = eligible
    out["variant_focal_counts"] = focal
    out["variant_reference_counts"] = reference
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>",
              file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
