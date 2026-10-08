"""Build gold, warm bundles and the grading spec for the kinetics ladder.

Usage: python -m npbench_c.templates.kinetics_consistency.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.kinetics_consistency.analyse import analyse, load_params
from npbench_c.templates.kinetics_consistency.core import (
    CONSTANTS,
    RULES,
    KineticsError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.003, "r3": 0.001, "r4": 0.001}
#: One unit in the last reported decimal of a rate.
RATE_TOLERANCE = 0.000001

CENSUS_KEYS = ("release", "ec_prefix", "distinct_ec_numbers", "declared_units",
               "records_by_field", "parsed_by_field", "sentinel_by_field",
               "unparsed_by_field")
TRIPLE_KEYS = ("triple_keys_total", "completeness_counts", "complete_triples",
               "ec_with_complete_triples",
               "distinct_substrates_in_complete_triples")
CONSISTENCY_KEYS = ("checked", "excluded_non_positive", "band_counts",
                    "agreement_rate", "unit_mixup_total")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if not any(gold["sentinel_by_field"].values()):
        raise KineticsError(
            "no sentinel values in the corpus, so the -999 trap R1 is built "
            "around does not exist here")
    if gold["completeness_counts"]["3"] < 100:
        raise KineticsError(
            f"only {gold['completeness_counts']['3']} complete triples; the "
            "consistency rate would be an anecdote rather than a measurement")
    if not gold["unit_mixup_total"]:
        raise KineticsError(
            "no triple falls in a unit-mix-up band, so the campaign's claim that "
            "the implicit convention is load-bearing has no instances")
    variants = gold["unit_error_variants"]
    identity = variants.get("v0_identity")
    if identity is None or identity["agreement_rate"] != gold["agreement_rate"]:
        raise KineticsError(
            "the identity variant does not reproduce the headline agreement rate, "
            "so the R4 control is not controlling anything")
    moved = [v["agreement_rate"] for n, v in variants.items() if n != "v0_identity"]
    if not moved or max(moved) >= gold["agreement_rate"] / 10:
        raise KineticsError(
            f"a declared unit mistake barely moves the agreement rate "
            f"({moved} against {gold['agreement_rate']}), so R4 shows nothing")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", CENSUS_KEYS, "parse_census.json"),
                             ("r3", TRIPLE_KEYS, "triples.json"),
                             ("r4", CONSISTENCY_KEYS, "consistency.json")):
        (warm / rung).mkdir(parents=True)
        (warm / rung / name).write_text(
            json.dumps({k: gold[k] for k in keys}, indent=2, sort_keys=True) + "\n")


def _flatten_for_grading(gold: dict) -> dict:
    variants = gold["unit_error_variants"]
    gold["unit_error_rates"] = {k: v["agreement_rate"]
                                for k, v in sorted(variants.items())}
    gold["unit_error_agree"] = {k: v["agree"] for k, v in sorted(variants.items())}
    gold["unit_error_factors"] = {
        k: {"km_factor": v["km_factor"], "kcat_factor": v["kcat_factor"]}
        for k, v in sorted(variants.items())}
    gold["sentinel_and_unparsed"] = {
        "sentinel": gold["sentinel_by_field"], "unparsed": gold["unparsed_by_field"]}
    return gold


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    grading = task.get("grading") or {}
    chance = {**DEFAULT_CHANCE_LEVELS, **(grading.get("chance_levels") or {})}
    notes = grading.get("chance_level_notes") or {}

    def ex(name, path=None):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{path or name}"}}

    def tol(name, path=None, abs_tol=RATE_TOLERANCE):
        return {"name": name, "scorer": "abs_tol_match",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{path or name}",
                         "abs_tol": abs_tol}}

    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("release"), ex("ec_prefix"), ex("distinct_ec_numbers"),
          ex("declared_units"), ex("records_by_field"), ex("parsed_by_field"),
          ex("sentinel_and_unparsed")]

    r2 = [ex("triple_keys_total"), ex("completeness_counts"),
          ex("complete_triples"), ex("ec_with_complete_triples"),
          ex("distinct_substrates_in_complete_triples")]

    r3 = [ex("checked"), ex("excluded_non_positive"), ex("band_counts"),
          tol("agreement_rate"), ex("unit_mixup_total")]

    r4 = [ex("unit_error_factors"), ex("unit_error_agree"),
          ex("unit_error_rates")]

    rungs = []
    for rid, ordinal, label, components in (("r1", 1, "execute", r1),
                                            ("r2", 2, "correct", r2),
                                            ("r3", 3, "commit", r3),
                                            ("r4", 4, "counterfactual", r4)):
        rung = {"rung_id": rid, "ordinal": ordinal, "name": label, "gold_tier": "G1",
                "composition": "product", "pass_threshold": 1.0,
                "chance_level": chance[rid], "components": components}
        if rid in notes:
            rung["chance_level_note"] = notes[rid]
        rungs.append(rung)

    spec = {"campaign_id": task["campaign_id"], "template_id": task["template_id"],
            "grading_spec_version": "1.0.0", "rungs": rungs}
    (pathlib.Path(campaign) / "grading.json").write_text(
        json.dumps(spec, indent=2) + "\n")
    return spec


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)

    gold = _flatten_for_grading(analyse(campaign))
    _assert_has_content(gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: BRENDA {gold['release']}, EC {gold['ec_prefix']}- , "
          f"{gold['distinct_ec_numbers']} EC numbers")
    print("  field                records  parsed  sentinel  unparsed")
    for field in sorted(gold["records_by_field"]):
        print(f"    {field:<18} {gold['records_by_field'][field]:>7} "
              f"{gold['parsed_by_field'][field]:>7} "
              f"{gold['sentinel_by_field'][field]:>9} "
              f"{gold['unparsed_by_field'][field]:>9}")
    print(f"  {gold['triple_keys_total']} keys, completeness "
          f"{gold['completeness_counts']}, {gold['complete_triples']} complete "
          f"over {gold['ec_with_complete_triples']} EC numbers")
    print(f"  bands {gold['band_counts']}; agreement "
          f"{gold['agreement_rate']:.4f}; unit mix-ups {gold['unit_mixup_total']}")
    print("  unit-error counterfactual: "
          f"{ {k: v for k, v in gold['unit_error_rates'].items()} }")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
