"""Build gold, warm bundles and the grading spec for the release-diff ladder.

Usage: python -m npbench_c.templates.mibig_diff.build <campaign>
"""

from __future__ import annotations

import copy
import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.mibig_diff.core import CONSTANTS, DiffError, ReleaseDiff, load_release

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.005, "r3": 0.001, "r4": 0.002}


def load_task(campaign: pathlib.Path) -> dict:
    return yaml.safe_load((pathlib.Path(campaign) / "task.yaml").read_text())


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    for name in ("field_mapping.json", "class_mapping.json"):
        shutil.copyfile(CONSTANTS / name, reference / name)


def _apply(record: dict, operation: str) -> dict:
    """Apply one declared perturbation to a copy of a 3.1 record.

    Perturbations are declared operations rather than hand-edited files, so each
    case is regenerable and its expected verdict is computed by the oracle
    rather than asserted by an author.
    """
    out = copy.deepcopy(record)
    cluster = out["cluster"]
    if operation == "taxid_to_int":
        cluster["ncbi_tax_id"] = int(cluster["ncbi_tax_id"])
    elif operation == "taxid_to_other_value":
        cluster["ncbi_tax_id"] = "999999"
    elif operation == "first_formula_altered":
        compounds = cluster.get("compounds") or []
        target = next((c for c in compounds if c.get("molecular_formula")), None)
        if target is None:
            raise DiffError("base entry has no compound formula to perturb")
        target["molecular_formula"] = target["molecular_formula"] + "S"
    elif operation == "class_to_retired_term":
        cluster["biosyn_class"] = ["Alkaloid"]
    elif operation == "organism_removed":
        cluster.pop("organism_name", None)
    else:
        raise DiffError(f"unknown perturbation operation {operation!r}")
    return out


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    diff = ReleaseDiff.load(campaign)
    params = load_task(campaign)["template_params"]

    source = load_release(campaign / params["source_archive"])
    target = load_release(campaign / params["target_archive"])
    corpus = diff.corpus_diff(source, target)

    # The campaign only means anything if the normalisation actually changes the
    # answer somewhere; otherwise R3 is a restatement of R2.
    gaps = {f: corpus["raw_difference_counts"][f]
                - corpus["classified_verdict_counts"][f].get("substantive_change", 0)
            for f in sorted(diff.fields)}
    if not any(v > 0 for v in gaps.values()):
        raise DiffError(
            "no field differs between the raw and classified counts, so R3 would "
            "restate R2 and the campaign has no normalisation content")

    base = params["perturbation_base"]
    if base not in source or base not in target:
        raise DiffError(f"perturbation base {base} is not shared by both releases")

    perturbations = {}
    for name, operation in sorted(params["perturbations"].items()):
        mutated = _apply(source[base], operation["operation"])
        field = operation["field"]
        perturbations[name] = {
            "operation": operation["operation"],
            "field": field,
            "baseline_verdict": diff.classify_field(field, source[base], target[base]),
            "verdict": diff.classify_field(field, mutated, target[base]),
        }

    return {
        "source_release": params["source_release"],
        "target_release": params["target_release"],
        "fields": sorted(diff.fields),
        "retired_class_terms": sorted(diff.retired_terms),
        **corpus,
        "raw_versus_classified_gap": gaps,
        "perturbation_base": base,
        "perturbations": perturbations,
    }


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    # R2 starts from R1's deliverable: both corpora loaded and sized.
    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "corpora_loaded.json").write_text(json.dumps({
        "source_entries": gold["source_entries"],
        "target_entries": gold["target_entries"],
        "fields": gold["fields"],
    }, indent=2, sort_keys=True) + "\n")

    # R3 starts from R2's deliverable: the entry-set diff and the RAW counts,
    # which is exactly the naive comparison. R3's work is the normalisation.
    (warm / "r3").mkdir(parents=True)
    (warm / "r3" / "raw_diff.json").write_text(json.dumps({
        "only_in_source": gold["only_in_source"],
        "only_in_target": gold["only_in_target"],
        "shared_entries": gold["shared_entries"],
        "shared_status_in_target": gold["shared_status_in_target"],
        "raw_difference_counts": gold["raw_difference_counts"],
    }, indent=2, sort_keys=True) + "\n")

    # R4 adds R3's classification -- and nothing about the perturbations.
    (warm / "r4").mkdir(parents=True)
    (warm / "r4" / "classified.json").write_text(json.dumps({
        "classified_verdict_counts": gold["classified_verdict_counts"],
        "retired_class_terms": gold["retired_class_terms"],
        "perturbation_base": gold["perturbation_base"],
    }, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = load_task(campaign)
    overrides = ((task.get("grading") or {}).get("chance_levels") or {})
    chance = {**DEFAULT_CHANCE_LEVELS, **overrides}
    notes = ((task.get("grading") or {}).get("chance_level_notes") or {})
    diff = ReleaseDiff.load(pathlib.Path(campaign))

    def one(name, pred):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": 1}}

    def ex(name, pred, g):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}"}}

    r1 = [one("report_present", "present.report"),
          ex("source_entries", "source_entries", "source_entries"),
          ex("target_entries", "target_entries", "target_entries"),
          {"name": "fields_enumerated", "scorer": "set_exact",
           "args": {"pred": "pred:fields", "gold": "gold:fields"}}]

    r2 = [ex("only_in_source", "only_in_source", "only_in_source"),
          ex("only_in_target", "only_in_target", "only_in_target"),
          ex("shared_entries", "shared_entries", "shared_entries"),
          ex("shared_status", "shared_status_in_target", "shared_status_in_target")]
    r2 += [ex(f"raw_diff_{f}", f"raw_difference_counts.{f}",
              f"raw_difference_counts.{f}") for f in gold["fields"]]

    r3 = [ex(f"verdicts_{f}", f"classified_verdict_counts.{f}",
             f"classified_verdict_counts.{f}") for f in gold["fields"]]
    # retired_class_terms is NOT graded. The declared class mapping has to publish
    # which source terms have no counterpart in the target -- that is what makes
    # the campaign solvable -- so reporting them back is reading the shipped table,
    # not classifying anything. It stays in the report as audit material and is
    # recorded in excluded_from_grading.

    r4 = []
    for name in sorted(gold["perturbations"]):
        r4.append({"name": f"perturbation_{name}", "scorer": "enum_exact",
                   "args": {"pred": f"pred:perturbations.{name}.verdict",
                            "gold": f"gold:perturbations.{name}.verdict",
                            "vocabulary": list(diff.verdict_vocabulary)}})

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

    print(f"{campaign.name}: {gold['source_release']} -> {gold['target_release']}")
    print(f"  entries {gold['source_entries']} -> {gold['target_entries']}  "
          f"only_in_source={gold['only_in_source']} only_in_target={gold['only_in_target']} "
          f"shared={gold['shared_entries']}")
    print(f"  status of shared in target: {gold['shared_status_in_target']}")
    print("  field                 raw-diff  classified")
    for f in gold["fields"]:
        print(f"    {f:<20} {gold['raw_difference_counts'][f]:>8}  "
              f"{gold['classified_verdict_counts'][f]}")
    print(f"  perturbations: "
          f"{ {k: v['verdict'] for k, v in sorted(gold['perturbations'].items())} }")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
