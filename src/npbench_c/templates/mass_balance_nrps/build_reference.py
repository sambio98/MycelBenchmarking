"""Resolve monomer formulas once, at instantiation, and pin them.

RDKit is used HERE ONLY, to turn the SMILES MIBiG supplies for each A-domain
substrate into a molecular formula. The result is pinned into the campaign's
reference/monomer_formulas.json with its SMILES source and the RDKit version, so
the oracle and the agent read the same table, neither needs a chemistry toolkit
at run time, and a toolkit version bump cannot move a published score.

Formula derivation is atom counting plus valence-implied hydrogens, stable
across RDKit versions. Canonical SMILES is not stable and is never a key.

Usage: python -m npbench_c.templates.mass_balance_nrps.build_reference <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

from npbench_c.templates.mass_balance_nrps.core import CONSTANTS


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    reference = campaign / "reference"
    reference.mkdir(parents=True, exist_ok=True)

    # Template-level constants are copied in, so the campaign's sandbox is
    # self-contained while the template remains the single source of truth.
    for name in ("assembly_rules.json", "atomic_masses.json"):
        shutil.copyfile(CONSTANTS / name, reference / name)

    from rdkit import Chem, __version__ as rdkit_version
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    entry = json.loads((campaign / "inputs" / "entry.json").read_text())
    monomers: dict[str, dict] = {}
    for index, module in enumerate(entry["biosynthesis"]["modules"], start=1):
        subs = (module.get("a_domain") or {}).get("substrates") or []
        if len(subs) != 1:
            print(f"module {index}: expected one substrate, got {len(subs)}",
                  file=sys.stderr)
            return 1
        sub = subs[0]
        mol = Chem.MolFromSmiles(sub["structure"])
        if mol is None:
            print(f"module {index}: RDKit could not parse {sub['structure']!r}",
                  file=sys.stderr)
            return 1
        formula = CalcMolFormula(mol)
        name = sub["name"]
        if name in monomers and monomers[name]["formula"] != formula:
            print(f"monomer {name!r} resolves to two different formulas",
                  file=sys.stderr)
            return 1
        monomers[name] = {
            "formula": formula,
            "smiles_source": sub["structure"],
            "proteinogenic": bool(sub.get("proteinogenic")),
        }

    (reference / "monomer_formulas.json").write_text(json.dumps({
        "description": "Monomer molecular formulas for the monomers this entry's "
                       "modules load, derived from the SMILES MIBiG supplies for "
                       "each A-domain substrate.",
        "derived_with": f"rdkit {rdkit_version} CalcMolFormula",
        "derivation_note": "Atom counting plus valence-implied hydrogens, which is "
                           "stable across RDKit versions. Canonical SMILES is not "
                           "stable and is never used as a key.",
        "monomers": monomers,
    }, indent=2, sort_keys=True) + "\n")

    print(f"{len(monomers)} monomers -> {reference.name}/monomer_formulas.json")
    for name, spec in sorted(monomers.items()):
        print(f"  {spec['formula']:<12} {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
