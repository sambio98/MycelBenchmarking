"""Reference solution. Must score exactly 1.0 or the campaign does not ship."""

from __future__ import annotations

import json
import pathlib
import sys

CAMPAIGN = pathlib.Path(__file__).resolve().parent.parent


def main() -> int:
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else CAMPAIGN / "oracle_submission"
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())

    route = {
        "modules": [
            {"module": i + 1, "monomer": name,
             "formula": gold["monomer_formulas"][name]}
            for i, name in enumerate(gold["monomers"])
        ],
        "peptide_bonds": gold["peptide_bonds"],
        "naive_assembly_formula": gold["naive_assembly_formula"],
        "naive_assembly_monoisotopic_mass": gold["naive_assembly_monoisotopic_mass"],
        "per_compound": {
            name: {"residual": block["residual"], "verdict": block["verdict"]}
            for name, block in gold["per_compound"].items()
        },
        "single_module_deletions": {
            k: {"formula": v["formula"], "mass_shift": v["mass_shift"]}
            for k, v in gold["single_module_deletions"].items() if v.get("feasible")
        },
    }
    (out / "route.json").write_text(json.dumps(route, indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
