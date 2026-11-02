"""Everything the kinetics ladder grades, from the sandbox alone."""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.kinetics_consistency.core import (
    FIELDS,
    Kinetics,
    KineticsError,
)


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    kinetics = Kinetics.load(campaign, params["kinetics_input"])

    census = kinetics.census()
    table = kinetics.triples()
    completeness = kinetics.completeness(table)
    consistency = kinetics.consistency(table)

    complete = {k: v for k, v in table.items() if len(v) == 3}
    bands = consistency["band_counts"]
    variants = {}
    for name, spec in sorted((params.get("unit_error_variants") or {}).items()):
        variants[name] = {
            "km_factor": float(spec["km_factor"]),
            "kcat_factor": float(spec["kcat_factor"]),
            **kinetics.under_unit_error(table, float(spec["km_factor"]),
                                        float(spec["kcat_factor"])),
        }

    return {
        "release": kinetics.release,
        "ec_prefix": params["ec_prefix"],
        "distinct_ec_numbers": len(kinetics.corpus),
        "declared_units": dict(sorted(kinetics.rules.units.items())),
        **census,
        "triple_keys_total": len(table),
        "completeness_counts": completeness,
        "complete_triples": consistency["complete_triples"],
        "ec_with_complete_triples": len({k[0] for k in complete}),
        "distinct_substrates_in_complete_triples": len({k[2] for k in complete}),
        "checked": consistency["checked"],
        "excluded_non_positive": consistency["excluded_non_positive"],
        "band_counts": bands,
        "agreement_rate": consistency["agreement_rate"],
        "unit_mixup_total": bands.get("factor_60", 0) + bands.get("factor_1000", 0),
        "unit_error_variants": variants,
        # Audit material, not graded: the first few triples in each unit-mix-up
        # band, so a reader can see what the campaign is claiming.
        "band_examples": {
            name: kinetics.band_examples(table, name)
            for name in ("factor_60", "factor_1000")
        },
    }
