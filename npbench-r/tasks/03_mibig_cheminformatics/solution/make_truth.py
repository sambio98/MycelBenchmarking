"""Reference solution for tasks/03_mibig_cheminformatics.

Reads the frozen MIBiG 4.0 JSON release and writes the eight graded artifacts.

    python solution/make_truth.py --task tasks/03_mibig_cheminformatics --out <dir>

Every output is sorted on its row key, so two runs produce byte-identical files.

Three properties of the input shape the method, and all three are MIBiG's own
data rather than anything introduced here:

  * MIBiG states a `formula` and a `mass` per compound, and both disagree with
    the structure for a few hundred entries. The artifacts are computed from
    the SMILES, never copied from those fields.
  * A BGC can list the same compound name twice with different structures.
    Four such pairs exist, so the row key needs a tie-break rather than an
    assumption of uniqueness.
  * Hundreds of distinct structures collapse onto each other once
    stereochemistry is dropped, and hundreds more share a molecular formula.
"""

from __future__ import annotations

import argparse
import csv
import json
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Scaffolds import MurckoScaffold

# RDKit logs parse warnings to stderr; they are not failures and would bury the
# progress output.
RDLogger.DisableLog("rdApp.*")

ARCHIVE = "data/mibig_json_4.0.tar.gz"

# Morgan fingerprint and clustering parameters. Fixed here and in the prompt:
# Tanimoto similarity depends on radius and bit size, so a family assignment is
# only reproducible against these exact values.
MORGAN_RADIUS = 2
MORGAN_BITS = 2048
TANIMOTO_THRESHOLD = 0.50

MASS_DECIMALS = 4

# The pip pin is rdkit==2026.9.1; rdkit.__version__ prints it as "2026.09.1".
# Both forms appear in the wild, so check the one the library reports.
RDKIT_VERSION = "2026.09.1"


def check_version() -> None:
    import rdkit

    if rdkit.__version__ != RDKIT_VERSION:
        raise SystemExit(
            f"this task's truth was produced by RDKit {RDKIT_VERSION}, found "
            f"{rdkit.__version__}. Canonical SMILES, scaffold perception and "
            f"fingerprint bits are version-sensitive, so the artifacts would differ."
        )


def read_records(archive: Path) -> list[dict]:
    """Every MIBiG entry, in accession order, read straight from the tarball."""
    records = []
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith(".json"):
                handle = tar.extractfile(member)
                if handle is not None:
                    records.append(json.load(handle))
    records.sort(key=lambda r: r["accession"])
    return records


def canonical(smiles: str) -> tuple[Chem.Mol, str] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol, Chem.MolToSmiles(mol)


def flatten(mol: Chem.Mol) -> str:
    """Canonical SMILES with stereochemistry removed."""
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=False)))


def scaffold_of(mol: Chem.Mol) -> tuple[str, int, int]:
    """Bemis-Murcko scaffold SMILES, its ring count and its heavy-atom count.

    An acyclic molecule has an empty scaffold; that is RDKit's answer and it is
    kept rather than dropped, because those compounds are part of the universe.
    """
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    smiles = Chem.MolToSmiles(scaffold)
    if not smiles:
        return "", 0, 0
    return smiles, rdMolDescriptors.CalcNumRings(scaffold), scaffold.GetNumHeavyAtoms()


def single_linkage(labels: list[str], neighbours: dict[int, set[int]]) -> dict[str, str]:
    """Connected components, named for the lexicographically smallest member.

    Same convention as task 02: naming by rank would make every row depend on a
    tie-break nobody can reproduce, so the name comes from the component itself.
    """
    parent = list(range(len(labels)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, group in neighbours.items():
        for j in group:
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[max(ri, rj)] = min(ri, rj)

    members: dict[int, list[str]] = defaultdict(list)
    for i, label in enumerate(labels):
        members[find(i)].append(label)
    names = {root: min(group) for root, group in members.items()}
    return {labels[i]: names[find(i)] for i in range(len(labels))}


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    check_version()
    records = read_records(args.task / ARCHIVE)
    print(f"{len(records)} MIBiG entries")

    # One entry per (accession, compound name). Where a BGC lists a name twice,
    # the lexicographically smallest canonical SMILES wins, so the key is unique
    # without discarding a compound the agent should have seen.
    #
    # The key is case-folded. MIBiG has four exact duplicate names and two that
    # differ only in capitalisation ("spumigin E" and "Spumigin E" in
    # BGC0000430, "salinipostin G" in BGC0001458). Deduplicating on the raw name
    # leaves those two pairs in, and since the grader compares compound names
    # case-insensitively they then collide as duplicate row keys and fail the
    # artifact. The dedup key has to match the key the grader builds.
    best: dict[tuple[str, str], dict] = {}
    classes_of: dict[str, list[str]] = {}
    for record in records:
        accession = record["accession"]
        classes_of[accession] = sorted(
            {
                entry["class"]
                for entry in (record.get("biosynthesis") or {}).get("classes") or []
                if entry.get("class")
            }
        )
        for compound in record.get("compounds") or []:
            smiles = compound.get("structure")
            name = (compound.get("name") or "").strip()
            if not smiles or not name:
                continue
            parsed = canonical(smiles)
            if parsed is None:
                continue
            mol, canonical_smiles = parsed
            key = (accession, name.casefold())
            if key in best and best[key]["canonical_smiles"] <= canonical_smiles:
                continue
            best[key] = {
                "bgc_accession": accession,
                "compound_name": name,
                "canonical_smiles": canonical_smiles,
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "exact_mass": round(ExactMolWt(mol), MASS_DECIMALS),
                "n_heavy_atoms": mol.GetNumHeavyAtoms(),
                "n_rings": rdMolDescriptors.CalcNumRings(mol),
                "stated_formula": (compound.get("formula") or "").strip(),
                "_mol": mol,
            }

    compounds = [best[k] for k in sorted(best)]
    print(f"{len(compounds)} compounds with a usable structure")

    write_csv(
        args.out / "compounds.csv",
        ["bgc_accession", "compound_name", "canonical_smiles", "formula",
         "exact_mass", "n_heavy_atoms", "n_rings"],
        [{k: row[k] for k in ("bgc_accession", "compound_name", "canonical_smiles",
                              "formula", "exact_mass", "n_heavy_atoms", "n_rings")}
         for row in compounds],
    )

    write_csv(
        args.out / "formula_disagreements.csv",
        ["bgc_accession", "compound_name", "stated_formula", "computed_formula"],
        [
            {
                "bgc_accession": row["bgc_accession"],
                "compound_name": row["compound_name"],
                "stated_formula": row["stated_formula"],
                "computed_formula": row["formula"],
            }
            for row in compounds
            if row["stated_formula"] and row["stated_formula"] != row["formula"]
        ],
    )

    # From here on the unit is the distinct structure, not the MIBiG entry:
    # the same molecule is reported by several BGCs.
    structures: dict[str, Chem.Mol] = {}
    bgcs_of: dict[str, set[str]] = defaultdict(set)
    for row in compounds:
        structures.setdefault(row["canonical_smiles"], row["_mol"])
        bgcs_of[row["canonical_smiles"]].add(row["bgc_accession"])
    smiles_list = sorted(structures)
    print(f"{len(smiles_list)} distinct structures")

    scaffolds = {s: scaffold_of(structures[s]) for s in smiles_list}
    write_csv(
        args.out / "scaffolds.csv",
        ["canonical_smiles", "scaffold_smiles", "n_scaffold_rings", "n_scaffold_atoms"],
        [
            {"canonical_smiles": s, "scaffold_smiles": scaffolds[s][0],
             "n_scaffold_rings": scaffolds[s][1], "n_scaffold_atoms": scaffolds[s][2]}
            for s in smiles_list
        ],
    )

    per_scaffold: dict[str, set[str]] = defaultdict(set)
    scaffold_bgcs: dict[str, set[str]] = defaultdict(set)
    scaffold_classes: dict[str, set[str]] = defaultdict(set)
    for s in smiles_list:
        scaffold = scaffolds[s][0]
        per_scaffold[scaffold].add(s)
        scaffold_bgcs[scaffold] |= bgcs_of[s]
        for accession in bgcs_of[s]:
            scaffold_classes[scaffold] |= set(classes_of[accession])
    write_csv(
        args.out / "scaffold_counts.csv",
        ["scaffold_smiles", "n_structures", "n_bgcs", "n_classes"],
        [
            {"scaffold_smiles": k, "n_structures": len(v),
             "n_bgcs": len(scaffold_bgcs[k]), "n_classes": len(scaffold_classes[k])}
            for k, v in sorted(per_scaffold.items())
        ],
    )

    by_formula: dict[str, set[str]] = defaultdict(set)
    for s in smiles_list:
        by_formula[rdMolDescriptors.CalcMolFormula(structures[s])].add(s)
    write_csv(
        args.out / "formula_collisions.csv",
        ["formula", "n_structures", "n_bgcs"],
        [
            {"formula": f, "n_structures": len(group),
             "n_bgcs": len({b for s in group for b in bgcs_of[s]})}
            for f, group in sorted(by_formula.items())
            if len(group) > 1
        ],
    )

    by_flat: dict[str, set[str]] = defaultdict(set)
    for s in smiles_list:
        by_flat[flatten(structures[s])].add(s)
    write_csv(
        args.out / "stereo_collapse_groups.csv",
        ["flat_smiles", "canonical_smiles", "group_size"],
        [
            {"flat_smiles": flat, "canonical_smiles": s, "group_size": len(group)}
            for flat, group in sorted(by_flat.items()) if len(group) > 1
            for s in sorted(group)
        ],
    )

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_RADIUS, fpSize=MORGAN_BITS
    )
    fingerprints = [generator.GetFingerprint(structures[s]) for s in smiles_list]
    neighbours: dict[int, set[int]] = defaultdict(set)
    for i in range(len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[i], fingerprints[i + 1:]
        )
        for offset, value in enumerate(similarities):
            if value >= TANIMOTO_THRESHOLD:
                neighbours[i].add(i + 1 + offset)
    families = single_linkage(smiles_list, neighbours)
    size = Counter(families.values())
    print(f"{len(size)} similarity families")
    write_csv(
        args.out / "similarity_families.csv",
        ["canonical_smiles", "family_id", "family_size"],
        [{"canonical_smiles": s, "family_id": families[s], "family_size": size[families[s]]}
         for s in smiles_list],
    )

    # A BGC may carry several biosynthesis classes; a structure counts once for
    # each distinct class across the BGCs that report it.
    cell_structures: dict[tuple[str, int], set[str]] = defaultdict(set)
    cell_bgcs: dict[tuple[str, int], set[str]] = defaultdict(set)
    for s in smiles_list:
        rings = scaffolds[s][1]
        for accession in bgcs_of[s]:
            for klass in classes_of[accession]:
                cell_structures[(klass, rings)].add(s)
                cell_bgcs[(klass, rings)].add(accession)
    write_csv(
        args.out / "class_scaffold_summary.csv",
        ["biosynthesis_class", "n_scaffold_rings", "n_structures", "n_bgcs"],
        [
            {"biosynthesis_class": klass, "n_scaffold_rings": rings,
             "n_structures": len(group), "n_bgcs": len(cell_bgcs[(klass, rings)])}
            for (klass, rings), group in sorted(cell_structures.items())
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
