"""Build gold, warm bundles and the grading spec for the selectivity ladder.

Usage: python -m npbench_c.templates.chembl_selectivity.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.chembl_selectivity.analyse import analyse, load_params
from npbench_c.templates.chembl_selectivity.core import (
    CONSTANTS,
    RULES,
    SelectivityError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.003, "r3": 0.001, "r4": 0.001}
RATIO_TOLERANCE = 0.000001

CENSUS_KEYS = ("release", "focal_target", "reference_target", "activities_total",
               "activities_of_declared_type", "distinct_assays",
               "assay_confidence_scores")
ADMISSION_KEYS = ("admitted_by_target", "rejected_by_target", "aggregated_pairs",
                  "aggregated_by_target", "eligible_molecules")
VERDICT_KEYS = ("focal_selective_counts", "reference_selective_counts",
                "focal_selective_at_highest_threshold",
                "most_selective_molecule", "most_selective_ratio")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    rejected = gold["rejected_by_target"]
    for target, counts in rejected.items():
        total = gold["admitted_by_target"][target] + sum(counts.values())
        # The rejection reasons are ordered and first-match-wins, so they must
        # partition the rejected set exactly. If they do not, two rules are
        # double-counting and every number downstream is suspect.
        if total != gold["_of_type_by_target"][target]:
            raise SelectivityError(
                f"{target}: admitted plus rejections is {total}, not "
                f"{gold['_of_type_by_target'][target]}; the rejection reasons do "
                "not partition the rejected set")
    if not any(sum(c.values()) for c in rejected.values()):
        raise SelectivityError(
            "no activity is rejected by any rule, so the admission contract R2 is "
            "built around does not bite on this input")
    highest = str(max(int(k) for k in gold["focal_selective_counts"]))
    if not gold["focal_selective_counts"][highest]:
        raise SelectivityError(
            "no molecule is selective for the focal target at the highest "
            "threshold, so R3's committed set is empty")
    if not gold["reference_selective_counts"]["2"]:
        raise SelectivityError(
            "no molecule favours the reference target even at twofold, so the "
            "verdict is a foregone conclusion in one direction")
    variants = gold["admission_variants"]
    baseline = variants.get("v0_baseline")
    if baseline is None or (baseline["focal_selective_counts"]
                            != gold["focal_selective_counts"]):
        raise SelectivityError(
            "the baseline variant does not reproduce R3's counts, so the R4 "
            "control is not controlling anything")
    moved = [n for n, v in variants.items() if n != "v0_baseline"
             and (v["eligible_molecules"] != baseline["eligible_molecules"]
                  or v["focal_selective_counts"] != baseline["focal_selective_counts"])]
    if len(moved) < len(variants) - 1:
        raise SelectivityError(
            f"only {moved} of the declared variants move anything; a variant that "
            "changes nothing is a component that cannot discriminate")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", CENSUS_KEYS, "activity_census.json"),
                             ("r3", ADMISSION_KEYS, "admitted.json"),
                             ("r4", VERDICT_KEYS, "selectivity.json")):
        (warm / rung).mkdir(parents=True)
        (warm / rung / name).write_text(
            json.dumps({k: gold[k] for k in keys}, indent=2, sort_keys=True) + "\n")


def _flatten_for_grading(gold: dict) -> dict:
    variants = gold["admission_variants"]
    gold["variant_admitted"] = {k: v["admitted_total"]
                                for k, v in sorted(variants.items())}
    gold["variant_eligible"] = {k: v["eligible_molecules"]
                                for k, v in sorted(variants.items())}
    gold["variant_focal_counts"] = {k: v["focal_selective_counts"]
                                    for k, v in sorted(variants.items())}
    gold["variant_reference_counts"] = {k: v["reference_selective_counts"]
                                        for k, v in sorted(variants.items())}
    return gold


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    grading = task.get("grading") or {}
    chance = {**DEFAULT_CHANCE_LEVELS, **(grading.get("chance_levels") or {})}
    notes = grading.get("chance_level_notes") or {}

    def ex(name):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{name}"}}

    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}}] + \
         [ex(k) for k in CENSUS_KEYS]
    r2 = [ex(k) for k in ADMISSION_KEYS]
    r3 = [ex("focal_selective_counts"), ex("reference_selective_counts"),
          {"name": "focal_selective_at_highest_threshold", "scorer": "set_exact",
           "args": {"pred": "pred:focal_selective_at_highest_threshold",
                    "gold": "gold:focal_selective_at_highest_threshold"}},
          ex("most_selective_molecule"),
          {"name": "most_selective_ratio", "scorer": "abs_tol_match",
           "args": {"pred": "pred:most_selective_ratio",
                    "gold": "gold:most_selective_ratio",
                    "abs_tol": RATIO_TOLERANCE}}]
    r4 = [ex("variant_admitted"), ex("variant_eligible"),
          ex("variant_focal_counts"), ex("variant_reference_counts")]

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
    # Per-target IC50 totals, needed only to check the rejection partition.
    params = load_params(campaign)
    from npbench_c.templates.chembl_selectivity.core import Selectivity

    selectivity = Selectivity.load(campaign, params["activities_input"],
                                   params["focal_target"],
                                   params["reference_target"])
    gold["_of_type_by_target"] = {
        t: sum(1 for a in selectivity.of_type() if a["target_chembl_id"] == t)
        for t in sorted((gold["focal_target"], gold["reference_target"]))}
    _assert_has_content(gold)
    gold.pop("_of_type_by_target")

    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: {gold['release']}, "
          f"{gold['activities_total']} activities, "
          f"{gold['activities_of_declared_type']} IC50 on the pair, "
          f"{gold['distinct_assays']} assays "
          f"(confidence {gold['assay_confidence_scores']})")
    for target in sorted(gold["admitted_by_target"]):
        print(f"  {target}: {gold['admitted_by_target'][target]} admitted, "
              f"rejected { {k: v for k, v in gold['rejected_by_target'][target].items() if v} }")
    print(f"  {gold['eligible_molecules']} molecules measured on both; "
          f"focal-selective {gold['focal_selective_counts']}, "
          f"reference-selective {gold['reference_selective_counts']}")
    print(f"  most selective: {gold['most_selective_molecule']} at "
          f"{gold['most_selective_ratio']:.1f}x")
    print("  variants (admitted / eligible / focal at 2x):")
    for name, v in sorted(gold["admission_variants"].items()):
        print(f"    {name:<24} {v['admitted_total']:>4} / "
              f"{v['eligible_molecules']:>4} / {v['focal_selective_counts']['2']:>4}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
