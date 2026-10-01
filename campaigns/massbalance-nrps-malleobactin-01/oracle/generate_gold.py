"""Generate gold and the leak-free warm-start bundles. Tier G1 throughout."""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import assembly as A
import formula as F

CAMPAIGN = pathlib.Path(__file__).resolve().parent.parent


def main() -> int:
    entry = json.loads((CAMPAIGN / "inputs" / "entry.json").read_text())
    analysis = A.analyse(entry)
    deletions = A.single_module_deletions(analysis["monomers"])

    gold = {
        "accession": entry["accession"],
        **analysis,
        "single_module_deletions": deletions,
    }
    (CAMPAIGN / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")

    write_warm_bundles(gold)
    print(json.dumps({k: v for k, v in gold.items()
                      if k not in ("single_module_deletions",)},
                     indent=2, sort_keys=True))
    return 0


def write_warm_bundles(gold: dict) -> None:
    """Each bundle holds exactly the earlier rungs' deliverables.

    Never a whole gold file: handing a rung "the previous answer" by copying
    gold.json would also hand it the answer it is being asked to compute. The
    R4 bundle in particular must not contain single_module_deletions.
    """
    warm = CAMPAIGN / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    # R2 starts from R1's deliverable: the module order and monomer identities,
    # without any formula resolved.
    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "route_skeleton.json").write_text(json.dumps({
        "n_modules": gold["n_modules"],
        "monomers": gold["monomers"],
        "compounds": gold["compounds"],
    }, indent=2, sort_keys=True) + "\n")

    # R3 starts from R2's deliverable: formulas resolved and the assembly closed.
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

    # R4 adds R3's deliverable -- residual and verdict -- and nothing about the
    # deletion counterfactuals.
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


if __name__ == "__main__":
    raise SystemExit(main())
