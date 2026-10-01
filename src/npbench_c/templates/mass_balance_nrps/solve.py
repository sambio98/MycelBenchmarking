"""Reference solution for a campaign of this template.

Must score exactly 1.0 through the real grader or the campaign does not ship.
Used by the readiness gate; never shipped to agents.

Usage: python -m npbench_c.templates.mass_balance_nrps.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: solve <campaign_dir> <submission_dir>", file=sys.stderr)
        return 2
    campaign, out = pathlib.Path(argv[0]), pathlib.Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((campaign / "gold" / "gold.json").read_text())

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
