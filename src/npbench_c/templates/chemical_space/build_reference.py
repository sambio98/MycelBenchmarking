"""Perceive every compound's chemistry once, at instantiation, and pin it.

RDKit is used HERE ONLY. It turns each structure MIBiG supplies into a molecular
formula, a standard InChIKey, a Bemis-Murcko scaffold and the two atom counts the
informativeness rule needs. The result is pinned into the campaign's
reference/chemistry_table.json with the RDKit version, so the oracle and the
agent read the same table, neither needs a chemistry toolkit at run time, and a
toolkit version bump cannot move a published score.

The table is the campaign's INPUT, not its answer. What it supplies is scaffold
perception, which needs a toolkit; what the campaign grades is the bookkeeping
over it -- which identity key, which filter, which join -- and none of that is
in the table.

Two things are deliberately not keys. Canonical SMILES is a toolkit output that
is not stable across versions, so scaffold_smiles ships as a human-readable
label and scaffold_inchikey is the key. And a compound's name is not unique
within an entry, so the record id is (accession, ordinal).

Usage: python -m npbench_c.templates.chemical_space.build_reference <campaign>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.chemical_space.core import CONSTANTS, TABLE
from npbench_c.templates.mibig_diff.core import load_release


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build_reference <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()

    import yaml

    params = yaml.safe_load((campaign / "task.yaml").read_text())["template_params"]
    corpus = load_release(campaign / params["corpus_archive"])

    from rdkit import Chem, RDLogger, __version__ as rdkit_version
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRings
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")

    records, failures = [], []
    for accession in sorted(corpus):
        entry = corpus[accession]
        if entry.get("status") != "active":
            continue
        for ordinal, compound in enumerate(entry.get("compounds") or []):
            smiles = compound.get("structure")
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                failures.append(f"{accession}:{ordinal}")
                continue
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            has_rings = scaffold.GetNumAtoms() > 0
            records.append({
                "record_id": f"{accession}:{ordinal}",
                "accession": accession,
                "ordinal": ordinal,
                "name": compound.get("name"),
                "formula": CalcMolFormula(mol),
                "inchikey": Chem.MolToInchiKey(mol),
                "compound_heavy_atoms": mol.GetNumHeavyAtoms(),
                "compound_ring_count": CalcNumRings(mol),
                "scaffold_smiles": Chem.MolToSmiles(scaffold) if has_rings else None,
                "scaffold_inchikey": Chem.MolToInchiKey(scaffold) if has_rings else None,
                "scaffold_heavy_atoms": (scaffold.GetNumHeavyAtoms() if has_rings else 0),
                "scaffold_ring_count": (CalcNumRings(scaffold) if has_rings else 0),
            })

    if failures:
        raise SystemExit(
            f"{len(failures)} structures did not parse: {failures[:10]}. The "
            "campaign assumes every shipped structure is perceivable; a parse "
            "failure is a corpus change that must be looked at, not skipped.")

    header = {
        "description": "Chemistry perceived once with RDKit and pinned. The "
                       "identity keys are inchikey, its first block, and "
                       "scaffold_inchikey; scaffold_smiles is a label and is "
                       "never a key. See space_rules.json.",
        "corpus_release": params["corpus_release"],
        "corpus_filter": "status == active and the compound carries a structure",
        "rdkit_version": rdkit_version,
        "scaffold_definition": "rdkit.Chem.Scaffolds.MurckoScaffold."
                               "GetScaffoldForMol, the Bemis-Murcko scaffold: "
                               "ring systems plus the linkers between them, "
                               "with all other substituents removed",
        "record_count": len(records),
    }
    # One record per line: readable in a terminal, diffable, and a third the
    # size of a fully indented table. The file is shipped to every agent, so its
    # size is a cost the campaign pays on every run.
    lines = [json.dumps(r, sort_keys=True, separators=(",", ":")) for r in records]
    body = ",\n  ".join(lines)
    head = json.dumps(header, indent=1, sort_keys=True)[1:-1].rstrip()
    reference = campaign / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    (reference / TABLE).write_text(
        "{" + head + ",\n \"records\": [\n  " + body + "\n ]\n}\n")

    distinct = len({r["inchikey"] for r in records})
    print(f"{campaign.name}: {len(records)} compound records from "
          f"{len({r['accession'] for r in records})} active entries")
    print(f"  RDKit {rdkit_version}; {distinct} distinct InChIKeys, "
          f"{len({r['inchikey'].split('-')[0] for r in records})} connectivity "
          f"blocks, "
          f"{len({r['scaffold_inchikey'] for r in records if r['scaffold_inchikey']})}"
          " scaffold keys, "
          f"{sum(1 for r in records if not r['scaffold_inchikey'])} acyclic")
    print(f"  wrote reference/{TABLE} "
          f"({(reference / TABLE).stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
