"""Resolve monomer formulas once, at campaign construction, and pin them.

RDKit is used HERE ONLY, to turn the SMILES that MIBiG supplies for each
A-domain substrate into a molecular formula. The result is written to
reference/monomer_formulas.json with its SMILES source and the RDKit version
recorded, so:

  * the oracle and the agent read the same pinned table
  * neither needs a chemistry toolkit at run time
  * a toolkit version bump cannot move a published score

Formula derivation is atom counting plus valence-implied hydrogens, which is
stable across RDKit versions. Canonical SMILES is NOT stable across versions and
is never used as a key here or anywhere else.
"""

from __future__ import annotations

import json
import pathlib
import sys

CAMPAIGN = pathlib.Path(__file__).resolve().parent.parent


def main() -> int:
    from rdkit import Chem, __version__ as rdkit_version
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    entry = json.loads((CAMPAIGN / "inputs" / "entry.json").read_text())
    modules = entry["biosynthesis"]["modules"]

    monomers = {}
    for index, module in enumerate(modules, start=1):
        subs = (module.get("a_domain") or {}).get("substrates") or []
        if len(subs) != 1:
            print(f"module {index}: expected exactly one substrate, got {len(subs)}",
                  file=sys.stderr)
            return 1
        sub = subs[0]
        mol = Chem.MolFromSmiles(sub["structure"])
        if mol is None:
            print(f"module {index}: RDKit could not parse {sub['structure']!r}",
                  file=sys.stderr)
            return 1
        name = sub["name"]
        formula = CalcMolFormula(mol)
        if name in monomers and monomers[name]["formula"] != formula:
            print(f"monomer {name!r} resolves to two different formulas", file=sys.stderr)
            return 1
        monomers[name] = {
            "formula": formula,
            "smiles_source": sub["structure"],
            "proteinogenic": bool(sub.get("proteinogenic")),
        }

    out = {
        "description": "Monomer molecular formulas for the monomers this entry's "
                       "modules load. Derived from the SMILES MIBiG supplies for "
                       "each A-domain substrate.",
        "derived_with": f"rdkit {rdkit_version} CalcMolFormula",
        "derivation_note": "Atom counting plus valence-implied hydrogens, which is "
                           "stable across RDKit versions. Canonical SMILES is not "
                           "stable and is never used as a key.",
        "monomers": monomers,
    }
    path = CAMPAIGN / "reference" / "monomer_formulas.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(f"{len(monomers)} monomers -> {path.relative_to(CAMPAIGN)}")
    for name, spec in sorted(monomers.items()):
        print(f"  {spec['formula']:<12} {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
