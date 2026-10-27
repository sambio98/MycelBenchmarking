"""Reference solution for tasks/04_mibig_nrps_pks_modules.

Parses the NRPS and PKS module annotations in the frozen MIBiG 4.0 release and
writes the eight graded artifacts.

    python solution/make_truth.py --task tasks/04_mibig_nrps_pks_modules --out <dir>

Every artifact derives from one parse of the module list, so the artifacts are
deliberately NOT independent: `modules.csv` feeds everything, and
`a_domain_substrates.csv` feeds five of the remaining six. That is the point.
Task 03's notes record why -- a task discriminates when its artifacts share
upstream dependencies, and task 03's did not.

Five properties of the input shape the method, all of them MIBiG's own:

  * Module names repeat within a BGC (12 pairs), so the row key is the module's
    1-based position in the `modules` list, which also forces the order to be
    preserved.
  * Over half the modules are marked `active: false`. They are annotated
    modules and are kept.
  * Almost every A-domain `location` is the sentinel `{from: -1, to: -1}`, so
    coordinates are not usable and are not asked for.
  * An A-domain carries between 0 and 9 substrates, and 295 carry no evidence.
  * Nine substrate names map to more than one distinct structure, `alanine`
    and `valine` among them, because MIBiG records D and L forms under the
    same name.
  * Twelve substrate names are spelled two ways, differing only in
    capitalisation ("Leucine" and "leucine"), so 281 verbatim names are 269
    substrates. Grouping on the raw name leaves duplicate row keys once the
    grader folds case -- the same trap task 03 hit.
"""

from __future__ import annotations

import argparse
import csv
import json
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Descriptors import ExactMolWt

RDLogger.DisableLog("rdApp.*")

ARCHIVE = "data/mibig_json_4.0.tar.gz"
RDKIT_VERSION = "2026.09.1"
MASS_DECIMALS = 4
# Strength ranking for adenylation-domain substrate evidence, strongest first.
# This is the task's rule and a design choice, not a universal truth; it is
# recorded in NOTES.md as such. In vitro assays that measure the A-domain
# directly rank above in vivo work, which ranks above observing the product,
# which ranks above inferring from the compound structure, which ranks above
# homology and bare sequence prediction.
#
# It covers all ten methods the 4.0 release uses. An earlier version listed
# only six, which let "Homology" outrank a direct enzyme assay -- the ranking
# has to be total over the data, or the fallback silently reorders it.
EVIDENCE_RANK = [
    "ATP-PPi exchange assay",
    "ACVS assay",
    "Enzyme-coupled assay",
    "Activity assay",
    "Heterologous expression",
    "Feeding study",
    "Mass spectrometry",
    "Structure-based inference",
    "Homology",
    "Sequence-based prediction",
]
NO_EVIDENCE = "none"

# The domain slots a module may carry, in a fixed order so the composition
# string is deterministic.
DOMAIN_SLOTS = ["a_domain", "c_domain", "ks_domain", "at_domain"]


def check_version() -> None:
    import rdkit

    if rdkit.__version__ != RDKIT_VERSION:
        raise SystemExit(
            f"this task's truth was produced by RDKit {RDKIT_VERSION}, found "
            f"{rdkit.__version__}. Substrate formulas and masses would differ."
        )


def read_records(archive: Path) -> list[dict]:
    records = []
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith(".json"):
                handle = tar.extractfile(member)
                if handle is not None:
                    records.append(json.load(handle))
    records.sort(key=lambda r: r["accession"])
    return records


def best_evidence(a_domain: dict) -> str:
    """The strongest evidence method on this A-domain, or "none"."""
    methods = {
        entry.get("method")
        for entry in (a_domain.get("evidence") or [])
        if entry.get("method")
    }
    for method in EVIDENCE_RANK:
        if method in methods:
            return method
    # An unranked method is reported verbatim rather than silently dropped,
    # so a new MIBiG release shows up as a new row instead of vanishing.
    return sorted(methods)[0] if methods else NO_EVIDENCE


def domain_types(module: dict) -> list[str]:
    """Every domain this module carries, as type strings, in a fixed order."""
    present = []
    for slot in DOMAIN_SLOTS:
        domain = module.get(slot)
        if domain:
            present.append(domain.get("type") or slot.removesuffix("_domain"))
    for domain in module.get("modification_domains") or []:
        present.append(domain.get("type") or "other")
    return present


def canonical_substrate_names(rows: list[dict]) -> dict[str, str]:
    """Map every verbatim substrate name onto one spelling per substrate.

    MIBiG spells twelve substrates two ways, differing only in capitalisation.
    The chosen spelling is the lexicographically smallest, which is
    deterministic and keeps MIBiG's own text rather than lowercasing names
    where capitalisation carries meaning, as in "N-methyl-L-phenylalanine".
    """
    variants: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        name = row["substrate_name"]
        variants[" ".join(name.split()).casefold()].add(name)
    chosen = {key: min(group) for key, group in variants.items()}
    return {
        name: chosen[" ".join(name.split()).casefold()]
        for group in variants.values()
        for name in group
    }


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

    # --- the one parse everything else depends on -------------------------
    modules: list[dict] = []
    substrates: list[dict] = []
    classes_of: dict[str, list[str]] = {}
    architecture: dict[str, list[str]] = defaultdict(list)
    composition: dict[str, list[str]] = defaultdict(list)

    for record in records:
        accession = record["accession"]
        biosynthesis = record.get("biosynthesis") or {}
        declared = biosynthesis.get("modules") or []
        if not declared:
            continue
        classes_of[accession] = sorted(
            {e["class"] for e in biosynthesis.get("classes") or [] if e.get("class")}
        )
        for position, module in enumerate(declared, start=1):
            present = domain_types(module)
            a_domain = module.get("a_domain")
            modules.append(
                {
                    "bgc_accession": accession,
                    "module_index": position,
                    "module_name": module.get("name") or "",
                    "module_type": module.get("type") or "",
                    "active": "true" if module.get("active") else "false",
                    "n_genes": len(module.get("genes") or []),
                    "has_a_domain": "true" if a_domain else "false",
                    "n_modification_domains": len(module.get("modification_domains") or []),
                    "n_domains": len(present),
                }
            )
            architecture[accession].append(module.get("type") or "")
            composition[accession].extend(present)

            if not a_domain:
                continue
            evidence = best_evidence(a_domain)
            for substrate in a_domain.get("substrates") or []:
                name = (substrate.get("name") or "").strip()
                if not name:
                    continue
                smiles = substrate.get("structure") or ""
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                substrates.append(
                    {
                        "bgc_accession": accession,
                        "module_index": position,
                        "substrate_name": name,
                        "proteinogenic": "true" if substrate.get("proteinogenic") else "false",
                        "substrate_formula": rdMolDescriptors.CalcMolFormula(mol) if mol else "",
                        "best_evidence": evidence,
                        "_canonical": Chem.MolToSmiles(mol) if mol else "",
                        "_mol": mol,
                    }
                )

    canonical_name = canonical_substrate_names(substrates)
    collapsed = sum(1 for k, v in canonical_name.items() if k != v)
    for row in substrates:
        row["substrate_name"] = canonical_name[row["substrate_name"]]

    print(f"{len(modules)} modules in {len(classes_of)} BGCs; "
          f"{len(substrates)} A-domain substrate assignments; "
          f"{collapsed} name spellings collapsed")

    modules.sort(key=lambda r: (r["bgc_accession"], r["module_index"]))
    write_csv(
        args.out / "modules.csv",
        ["bgc_accession", "module_index", "module_name", "module_type", "active",
         "n_genes", "has_a_domain", "n_modification_domains", "n_domains"],
        modules,
    )

    substrates.sort(key=lambda r: (r["bgc_accession"], r["module_index"], r["substrate_name"]))
    write_csv(
        args.out / "a_domain_substrates.csv",
        ["bgc_accession", "module_index", "substrate_name", "proteinogenic",
         "substrate_formula", "best_evidence"],
        [{k: r[k] for k in ("bgc_accession", "module_index", "substrate_name",
                            "proteinogenic", "substrate_formula", "best_evidence")}
         for r in substrates],
    )

    # --- derived from a_domain_substrates --------------------------------
    by_substrate: dict[str, list[dict]] = defaultdict(list)
    for row in substrates:
        by_substrate[row["substrate_name"]].append(row)
    write_csv(
        args.out / "substrate_counts.csv",
        ["substrate_name", "n_modules", "n_bgcs", "proteinogenic", "n_distinct_structures"],
        [
            {
                "substrate_name": name,
                "n_modules": len({(r["bgc_accession"], r["module_index"]) for r in group}),
                "n_bgcs": len({r["bgc_accession"] for r in group}),
                # MIBiG is not self-consistent on this flag for every name; the
                # value reported is the one the majority of entries carry, with
                # "true" winning a tie.
                "proteinogenic": "true"
                if sum(1 for r in group if r["proteinogenic"] == "true") * 2 >= len(group)
                else "false",
                "n_distinct_structures": len({r["_canonical"] for r in group if r["_canonical"]}),
            }
            for name, group in sorted(by_substrate.items())
        ],
    )

    by_evidence: dict[str, list[dict]] = defaultdict(list)
    for row in substrates:
        by_evidence[row["best_evidence"]].append(row)
    write_csv(
        args.out / "substrate_evidence_summary.csv",
        ["best_evidence", "n_assignments", "n_distinct_substrates", "n_bgcs"],
        [
            {
                "best_evidence": method,
                "n_assignments": len(group),
                "n_distinct_substrates": len({r["substrate_name"] for r in group}),
                "n_bgcs": len({r["bgc_accession"] for r in group}),
            }
            for method, group in sorted(by_evidence.items())
        ],
    )

    nonproteinogenic: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in substrates:
        if row["proteinogenic"] == "false" and row["_canonical"]:
            nonproteinogenic[(row["substrate_name"], row["_canonical"])].append(row)
    write_csv(
        args.out / "nonproteinogenic_substrates.csv",
        ["substrate_name", "canonical_smiles", "formula", "exact_mass",
         "n_heavy_atoms", "n_modules"],
        [
            {
                "substrate_name": name,
                "canonical_smiles": smiles,
                "formula": group[0]["substrate_formula"],
                "exact_mass": round(ExactMolWt(group[0]["_mol"]), MASS_DECIMALS),
                "n_heavy_atoms": group[0]["_mol"].GetNumHeavyAtoms(),
                "n_modules": len({(r["bgc_accession"], r["module_index"]) for r in group}),
            }
            for (name, smiles), group in sorted(nonproteinogenic.items())
        ],
    )

    cells: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in substrates:
        for klass in classes_of[row["bgc_accession"]]:
            cells[(klass, row["substrate_name"])].append(row)
    write_csv(
        args.out / "class_substrate_matrix.csv",
        ["biosynthesis_class", "substrate_name", "n_modules", "n_bgcs"],
        [
            {
                "biosynthesis_class": klass,
                "substrate_name": name,
                "n_modules": len({(r["bgc_accession"], r["module_index"]) for r in group}),
                "n_bgcs": len({r["bgc_accession"] for r in group}),
            }
            for (klass, name), group in sorted(cells.items())
        ],
    )

    # --- derived from the module order -----------------------------------
    accessions = sorted(architecture)
    a_domains_per_bgc = Counter(
        m["bgc_accession"] for m in modules if m["has_a_domain"] == "true"
    )
    arch_rows = [
        {
            "bgc_accession": accession,
            "n_modules": len(architecture[accession]),
            "n_a_domains": a_domains_per_bgc.get(accession, 0),
            "architecture": ">".join(architecture[accession]),
            "domain_composition": ";".join(sorted(composition[accession])),
        }
        for accession in accessions
    ]
    write_csv(
        args.out / "module_architecture.csv",
        ["bgc_accession", "n_modules", "n_a_domains", "architecture", "domain_composition"],
        arch_rows,
    )

    by_architecture: dict[str, list[dict]] = defaultdict(list)
    for row in arch_rows:
        by_architecture[row["architecture"]].append(row)
    print(f"{len(by_architecture)} distinct architectures over {len(accessions)} BGCs")
    write_csv(
        args.out / "architecture_counts.csv",
        ["architecture", "n_bgcs", "n_modules", "n_a_domains", "representative_bgc"],
        [
            {
                "architecture": arch,
                "n_bgcs": len(group),
                "n_modules": group[0]["n_modules"],
                "n_a_domains": sum(r["n_a_domains"] for r in group),
                "representative_bgc": min(r["bgc_accession"] for r in group),
            }
            for arch, group in sorted(by_architecture.items())
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
