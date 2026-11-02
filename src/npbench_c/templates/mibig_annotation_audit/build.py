"""Build gold, warm bundles and the grading spec for the annotation audit.

Usage: python -m npbench_c.templates.mibig_annotation_audit.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.mibig_annotation_audit.core import (
    CONSTANTS,
    AnnotationAudit,
    AuditError,
    Policy,
)
from npbench_c.templates.mibig_diff.core import load_release

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.004, "r3": 0.002, "r4": 0.001}


def load_task(campaign: pathlib.Path) -> dict:
    return yaml.safe_load((pathlib.Path(campaign) / "task.yaml").read_text())


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / "evidence_policy.json",
                    reference / "evidence_policy.json")


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    audit = AnnotationAudit.load(campaign)
    params = load_task(campaign)["template_params"]
    corpus = load_release(campaign / params["corpus_archive"])

    census = audit.census(corpus)
    baseline = audit.admissible(corpus, audit.baseline)

    gap = {f: census["backed_by_function"][f] - baseline["admissible_by_function"][f]
           for f in audit.function_vocabulary}
    if not any(v > 0 for v in gap.values()):
        raise AuditError(
            "no category loses annotations between backed and admissible, so R3 "
            "would restate R2 and the campaign has no evidence-policy content")

    variants = {}
    for name, spec in sorted(params["policy_variants"].items()):
        policy = Policy(
            name=name,
            function_evidence=frozenset(spec["function_evidence_accepted"]),
            locus_evidence=frozenset(spec["locus_evidence_accepted"]),
            require_active=bool(spec.get("require_active_status", True)),
        )
        result = audit.admissible(corpus, policy)
        variants[name] = {
            "admissible_by_function": result["admissible_by_function"],
            "total": sum(result["admissible_by_function"].values()),
        }

    baseline_total = sum(baseline["admissible_by_function"].values())
    if all(v["total"] == baseline_total for v in variants.values()):
        raise AuditError(
            "every policy variant gives the same admissible total as the "
            "baseline, so R4 has no counterfactual content")

    focal = params["focal_function"]
    if focal not in audit.function_vocabulary:
        raise AuditError(f"focal_function {focal!r} is not in the vocabulary")

    return {
        "corpus_release": params["corpus_release"],
        "corpus_entries": len(corpus),
        "function_vocabulary": sorted(audit.function_vocabulary),
        **census,
        "baseline_policy": audit.baseline.name,
        "admissible_by_function": baseline["admissible_by_function"],
        "admissible_total": baseline_total,
        "backed_minus_admissible": dict(sorted(gap.items())),
        "focal_function": focal,
        "focal_admissible_entries": baseline["admissible_entries_by_function"][focal],
        "policy_variants": variants,
    }


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "corpus_loaded.json").write_text(json.dumps({
        "corpus_entries": gold["corpus_entries"],
        "entries_with_annotations": gold["entries_with_annotations"],
        "total_function_annotations": gold["total_function_annotations"],
        "function_vocabulary": gold["function_vocabulary"],
    }, indent=2, sort_keys=True) + "\n")

    # R3 starts from R2's deliverable: the raw backed/unbacked census. R3's work
    # is applying the three-condition admissibility policy on top of it.
    (warm / "r3").mkdir(parents=True)
    (warm / "r3" / "census.json").write_text(json.dumps({
        "backed_by_function": gold["backed_by_function"],
        "unbacked_by_function": gold["unbacked_by_function"],
        "function_evidence_method_totals": gold["function_evidence_method_totals"],
    }, indent=2, sort_keys=True) + "\n")

    # R4 adds the baseline admissibility -- and nothing about the variants.
    (warm / "r4").mkdir(parents=True)
    (warm / "r4" / "baseline.json").write_text(json.dumps({
        "baseline_policy": gold["baseline_policy"],
        "admissible_by_function": gold["admissible_by_function"],
        "backed_minus_admissible": gold["backed_minus_admissible"],
        "focal_function": gold["focal_function"],
    }, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = load_task(campaign)
    overrides = ((task.get("grading") or {}).get("chance_levels") or {})
    chance = {**DEFAULT_CHANCE_LEVELS, **overrides}
    notes = ((task.get("grading") or {}).get("chance_level_notes") or {})
    functions = gold["function_vocabulary"]

    def one(name, pred):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": 1}}

    def ex(name, pred, g):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}"}}

    def key(f):
        return f.replace(" ", "_").replace("/", "_")

    r1 = [one("report_present", "present.report"),
          ex("corpus_entries", "corpus_entries", "corpus_entries"),
          ex("entries_with_annotations", "entries_with_annotations",
             "entries_with_annotations"),
          ex("total_function_annotations", "total_function_annotations",
             "total_function_annotations"),
          {"name": "function_vocabulary", "scorer": "set_exact",
           "args": {"pred": "pred:function_vocabulary",
                    "gold": "gold:function_vocabulary"}}]

    r2 = [ex(f"backed_{key(f)}", f"backed_by_function.{f}", f"backed_by_function.{f}")
          for f in functions]
    r2 += [ex(f"unbacked_{key(f)}", f"unbacked_by_function.{f}",
              f"unbacked_by_function.{f}") for f in functions]
    r2.append(ex("method_totals", "function_evidence_method_totals",
                 "function_evidence_method_totals"))

    r3 = [ex(f"admissible_{key(f)}", f"admissible_by_function.{f}",
             f"admissible_by_function.{f}") for f in functions]
    r3.append(ex("backed_minus_admissible", "backed_minus_admissible",
                 "backed_minus_admissible"))
    r3.append({"name": "focal_admissible_entries", "scorer": "set_exact",
               "args": {"pred": "pred:focal_admissible_entries",
                        "gold": "gold:focal_admissible_entries"}})

    r4 = []
    for name in sorted(gold["policy_variants"]):
        r4.append(ex(f"variant_{name}", f"policy_variants.{name}.admissible_by_function",
                     f"policy_variants.{name}.admissible_by_function"))
        r4.append(ex(f"variant_{name}_total", f"policy_variants.{name}.total",
                     f"policy_variants.{name}.total"))

    rungs = []
    for rid, ordinal, label, comps in (("r1", 1, "execute", r1), ("r2", 2, "correct", r2),
                                       ("r3", 3, "commit", r3),
                                       ("r4", 4, "counterfactual", r4)):
        rung = {"rung_id": rid, "ordinal": ordinal, "name": label, "gold_tier": "G1",
                "composition": "product", "pass_threshold": 1.0,
                "chance_level": chance[rid], "components": comps}
        if rid in notes:
            rung["chance_level_note"] = notes[rid]
        rungs.append(rung)

    spec = {"campaign_id": task["campaign_id"], "template_id": task["template_id"],
            "grading_spec_version": "1.0.0", "rungs": rungs}
    (pathlib.Path(campaign) / "grading.json").write_text(json.dumps(spec, indent=2) + "\n")
    return spec


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)
    gold = analyse(campaign)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: MIBiG {gold['corpus_release']}, "
          f"{gold['corpus_entries']} entries, "
          f"{gold['entries_with_annotations']} with gene annotations, "
          f"{gold['total_function_annotations']} function annotations")
    print("  function                  backed  unbacked  admissible  gap")
    for f in gold["function_vocabulary"]:
        print(f"    {f:<24} {gold['backed_by_function'][f]:>6} "
              f"{gold['unbacked_by_function'][f]:>9} "
              f"{gold['admissible_by_function'][f]:>11} "
              f"{gold['backed_minus_admissible'][f]:>4}")
    print(f"  admissible total {gold['admissible_total']}; "
          f"focal {gold['focal_function']} -> "
          f"{len(gold['focal_admissible_entries'])} entries")
    print("  policy variants: "
          f"{ {k: v['total'] for k, v in sorted(gold['policy_variants'].items())} }")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
