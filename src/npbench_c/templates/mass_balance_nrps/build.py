"""Build gold, warm-start bundles and the grading spec. Pure stdlib.

Separate from build_reference so the readiness gate's determinism check can
re-run gold generation without a chemistry toolkit, and so regenerating gold
cannot silently re-derive the pinned monomer table.

Usage: python -m npbench_c.templates.mass_balance_nrps.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.mass_balance_nrps.core import MassBalance

# Prior estimates. The internal sweep's no-computation ablation measures them
# and a campaign records the measured value with its provenance; R1's estimate
# was wrong on both instantiations built so far, which is the point of measuring.
DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.01, "r3": 0.004, "r4": 0.0005}

MASS_TOLERANCE = 0.001


def build_gold(campaign: pathlib.Path) -> dict:
    mb = MassBalance.load(campaign)
    entry = json.loads((campaign / "inputs" / "entry.json").read_text())
    analysis = mb.analyse(entry)
    gold = {
        "accession": entry["accession"],
        **analysis,
        "single_module_deletions": mb.single_module_deletions(analysis["monomers"]),
    }
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    return gold


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    """Each bundle holds exactly the earlier rungs' deliverables.

    Never a whole gold file: handing a rung "the previous answer" by copying
    gold.json would also hand over the answer it is being asked to compute. The
    R4 bundle must not contain single_module_deletions.
    """
    warm = campaign / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "route_skeleton.json").write_text(json.dumps({
        "n_modules": gold["n_modules"],
        "monomers": gold["monomers"],
        "compounds": gold["compounds"],
    }, indent=2, sort_keys=True) + "\n")

    (warm / "r3").mkdir(parents=True)
    (warm / "r3" / "assembly.json").write_text(json.dumps({
        "n_modules": gold["n_modules"],
        "monomers": gold["monomers"],
        "monomer_formulas": gold["monomer_formulas"],
        "peptide_bonds": gold["peptide_bonds"],
        "naive_assembly_formula": gold["naive_assembly_formula"],
        "naive_assembly_monoisotopic_mass": gold["naive_assembly_monoisotopic_mass"],
        "compounds": gold["compounds"],
    }, indent=2, sort_keys=True) + "\n")

    (warm / "r4").mkdir(parents=True)
    (warm / "r4" / "analysis.json").write_text(json.dumps({
        "n_modules": gold["n_modules"],
        "monomers": gold["monomers"],
        "monomer_formulas": gold["monomer_formulas"],
        "peptide_bonds": gold["peptide_bonds"],
        "naive_assembly_formula": gold["naive_assembly_formula"],
        "per_compound": {
            name: {"product_formula": block["product_formula"],
                   "residual": block["residual"],
                   "verdict": block["verdict"]}
            for name, block in gold["per_compound"].items()
        },
    }, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    """Generate grading.json from gold, so the spec cannot drift from the data."""
    mb = MassBalance.load(campaign)
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    overrides = ((task.get("grading") or {}).get("chance_levels") or {})
    chance = {**DEFAULT_CHANCE_LEVELS, **overrides}
    notes = ((task.get("grading") or {}).get("chance_level_notes") or {})
    vocab = list(mb.rules["verdict_vocabulary"])

    def one(name, pred):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": 1}}

    def ex(name, pred, g):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}"}}

    def tol(name, pred, g):
        return {"name": name, "scorer": "abs_tol_match",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}",
                         "abs_tol": MASS_TOLERANCE}}

    r1 = [one("route_present", "present.route"),
          one("all_formulas_parse", "all_formulas_parse"),
          ex("n_modules", "n_modules", "n_modules"),
          ex("monomer_order", "monomers", "monomers")]

    r2 = [ex("peptide_bonds", "peptide_bonds", "peptide_bonds"),
          ex("naive_assembly_formula", "naive_assembly_formula", "naive_assembly_formula"),
          tol("naive_assembly_mass", "naive_assembly_monoisotopic_mass",
              "naive_assembly_monoisotopic_mass")]
    r2 += [ex(f"formula_{m.replace(' ', '_')}", f"monomer_formulas.{m}",
              f"monomer_formulas.{m}") for m in sorted(gold["monomer_formulas"])]

    r3 = [{"name": "compounds_enumerated", "scorer": "set_exact",
           "args": {"pred": "pred:compounds", "gold": "gold:compounds"}}]
    for c in gold["compounds"]:
        key = c.replace(" ", "_")
        r3.append(ex(f"residual_{key}", f"per_compound.{c}.residual_full",
                     f"per_compound.{c}.residual_full"))
        r3.append({"name": f"verdict_{key}", "scorer": "enum_exact",
                   "args": {"pred": f"pred:per_compound.{c}.verdict",
                            "gold": f"gold:per_compound.{c}.verdict",
                            "vocabulary": vocab}})

    r4 = []
    for i in range(1, gold["n_modules"] + 1):
        block = gold["single_module_deletions"][str(i)]
        if not block.get("feasible"):
            continue
        r4.append(ex(f"del{i}_formula", f"deletions.{i}.formula",
                     f"single_module_deletions.{i}.formula"))
        r4.append(tol(f"del{i}_mass_shift", f"deletions.{i}.mass_shift",
                      f"single_module_deletions.{i}.mass_shift"))
    if not r4:
        raise SystemExit("no feasible single-module deletion: R4 would be empty")

    rungs = []
    for rid, ordinal, name, comps in (("r1", 1, "execute", r1), ("r2", 2, "correct", r2),
                                      ("r3", 3, "commit", r3),
                                      ("r4", 4, "counterfactual", r4)):
        rung = {"rung_id": rid, "ordinal": ordinal, "name": name, "gold_tier": "G1",
                "composition": "product", "pass_threshold": 1.0,
                "chance_level": chance[rid], "components": comps}
        if rid in notes:
            rung["chance_level_note"] = notes[rid]
        rungs.append(rung)

    spec = {"campaign_id": task["campaign_id"],
            "template_id": task["template_id"],
            "grading_spec_version": "1.1.0",
            "rungs": rungs}
    (campaign / "grading.json").write_text(json.dumps(spec, indent=2) + "\n")
    return spec


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    gold = build_gold(campaign)
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)
    print(f"{campaign.name}: {gold['n_modules']} modules, "
          f"{len(gold['compounds'])} compound(s), "
          f"assembly {gold['naive_assembly_formula']}")
    for name, block in sorted(gold["per_compound"].items()):
        print(f"  {name:<18} {block['product_formula']:<16} "
              f"residual={block['residual']} -> {block['verdict']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
