"""Deterministic stub systems for the kinetics-consistency ladder.

Harness fixtures. Sandbox-only.

  parser      clears R1: the per-field parse census, sentinels separated out
  triples     clears R2: the triple inventory under the declared key
  consistent  clears R3: the banded consistency verdict
  complete    clears R4
  micromolar  COMPETENCE probe, tooled: parses and keys everything correctly,
              then reads Km in micromolar where the convention is millimolar.
              Clears R1 and R2 and fails R3 -- which is the campaign's whole
              point, so the probe is the error the ladder is built to catch.
  noncompute  LEAKAGE ablation: never opens the input.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.kinetics_consistency import analyse as kin
from npbench_c.templates.kinetics_consistency.core import FIELDS, Kinetics

LEVELS = ("parser", "triples", "consistent", "complete", "micromolar", "noncompute")
REPORT = "kinetics_report.json"


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
        zero = {f: 0 for f in FIELDS}
        return write({
            "release": "2026.1", "ec_prefix": "1.1.", "distinct_ec_numbers": 400,
            "declared_units": {}, "records_by_field": zero,
            "parsed_by_field": zero, "sentinel_by_field": zero,
            "unparsed_by_field": zero, "triple_keys_total": 0,
            "completeness_counts": {"1": 0, "2": 0, "3": 0},
            "complete_triples": 0, "ec_with_complete_triples": 0,
            "distinct_substrates_in_complete_triples": 0, "checked": 0,
            "excluded_non_positive": 0, "band_counts": {},
            "agreement_rate": 0.5, "unit_mixup_total": 0,
            "unit_error_variants": {}, "band_examples": {},
        })

    params = kin.load_params(sandbox)
    kinetics = Kinetics.load(sandbox, params["kinetics_input"])
    table = kinetics.triples()

    report = {
        "release": kinetics.release, "ec_prefix": params["ec_prefix"],
        "distinct_ec_numbers": len(kinetics.corpus),
        "declared_units": dict(sorted(kinetics.rules.units.items())),
        **kinetics.census(),
    }
    if args.level == "parser":
        return write(report)

    complete = {k: v for k, v in table.items() if len(v) == 3}
    report.update({
        "triple_keys_total": len(table),
        "completeness_counts": kinetics.completeness(table),
        "complete_triples": len(complete),
        "ec_with_complete_triples": len({k[0] for k in complete}),
        "distinct_substrates_in_complete_triples": len({k[2] for k in complete}),
    })
    if args.level == "triples":
        return write(report)

    if args.level == "micromolar":
        # The declared convention is millimolar. Reading micromolar scales every
        # Km by a thousand and the agreement rate collapses -- which is exactly
        # what R4 measures deliberately and this probe does by mistake.
        wrong = kinetics.under_unit_error(table, 1000.0, 1.0)
        report.update({
            "checked": wrong["checked"], "excluded_non_positive": 0,
            "band_counts": {"agree": wrong["agree"], "near": 0, "factor_60": 0,
                            "factor_1000": 0, "other": wrong["checked"] - wrong["agree"]},
            "agreement_rate": wrong["agreement_rate"], "unit_mixup_total": 0,
        })
        return write(report)

    consistency = kinetics.consistency(table)
    bands = consistency["band_counts"]
    report.update({
        "checked": consistency["checked"],
        "excluded_non_positive": consistency["excluded_non_positive"],
        "band_counts": bands,
        "agreement_rate": consistency["agreement_rate"],
        "unit_mixup_total": bands.get("factor_60", 0) + bands.get("factor_1000", 0),
    })
    if args.level == "consistent":
        return write(report)

    report["unit_error_variants"] = {
        name: {"km_factor": float(spec["km_factor"]),
               "kcat_factor": float(spec["kcat_factor"]),
               **kinetics.under_unit_error(table, float(spec["km_factor"]),
                                           float(spec["kcat_factor"]))}
        for name, spec in sorted((params.get("unit_error_variants") or {}).items())}
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
