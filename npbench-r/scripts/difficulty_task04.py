"""Measure task 04's difficulty by grading plausible wrong pipelines.

Same approach and caveat as the task 02 and 03 scripts: synthetic variants, a
cheap proxy for the model runs PROJECT.md §8 asks for, useful mainly for
finding artifacts nothing can fail.

Task 04 was built so that a parse mistake propagates rather than stopping at
one artifact. This is the script that says whether that worked.

    python scripts/difficulty_task04.py --task tasks/04_mibig_nrps_pks_modules \
        --truth tasks/04_mibig_nrps_pks_modules/truth --scratch <dir>
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import rdMolDescriptors  # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

RDLogger.DisableLog("rdApp.*")

MISTAKES = {
    "filters_inactive_modules": "drops the 2,089 modules marked active:false",
    "keys_on_module_name": "keys modules on name instead of list position",
    "sorts_modules_by_name": "sorts the module list instead of preserving its order",
    "raw_substrate_names": "groups substrates on the raw name, keeping case variants",
    "filters_on_location": "drops A-domains whose location is the -1 sentinel",
    "first_evidence_not_strongest": "takes the first evidence entry, not the strongest",
    "proteinogenic_from_first_entry": "takes the first entry's flag, not the majority",
    "counts_rows_not_modules": "counts substrate rows where distinct modules are asked for",
    "class_counted_once": "attributes a cluster to only its first biosynthesis class",
    "unsorted_domain_composition": "leaves domain_composition in slot order, unsorted",
}

_RECORDS: list[dict] = []


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref04", task_dir / "solution" / "make_truth.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build(task_dir: Path, out: Path, mistakes: set[str]) -> None:
    ref = load_solution(task_dir)
    if not _RECORDS:
        _RECORDS.extend(ref.read_records(task_dir / ref.ARCHIVE))

    modules: list[dict] = []
    substrates: list[dict] = []
    classes_of: dict[str, list[str]] = {}
    architecture: dict[str, list[str]] = defaultdict(list)
    composition: dict[str, list[str]] = defaultdict(list)

    for record in _RECORDS:
        accession = record["accession"]
        biosynthesis = record.get("biosynthesis") or {}
        declared = list(biosynthesis.get("modules") or [])
        if not declared:
            continue
        classes_of[accession] = sorted(
            {e["class"] for e in biosynthesis.get("classes") or [] if e.get("class")}
        )
        if "filters_inactive_modules" in mistakes:
            declared = [m for m in declared if m.get("active")]
        if "sorts_modules_by_name" in mistakes:
            declared = sorted(declared, key=lambda m: m.get("name") or "")

        for position, module in enumerate(declared, start=1):
            present = ref.domain_types(module)
            a_domain = module.get("a_domain")
            index = module.get("name") if "keys_on_module_name" in mistakes else position
            modules.append({
                "bgc_accession": accession, "module_index": index,
                "module_name": module.get("name") or "",
                "module_type": module.get("type") or "",
                "active": "true" if module.get("active") else "false",
                "n_genes": len(module.get("genes") or []),
                "has_a_domain": "true" if a_domain else "false",
                "n_modification_domains": len(module.get("modification_domains") or []),
                "n_domains": len(present),
            })
            architecture[accession].append(module.get("type") or "")
            composition[accession].extend(present)

            if not a_domain:
                continue
            location = a_domain.get("location") or {}
            if "filters_on_location" in mistakes and (
                location.get("from") == -1 or location.get("to") == -1
            ):
                continue
            if "first_evidence_not_strongest" in mistakes:
                entries = a_domain.get("evidence") or []
                evidence = (entries[0].get("method") or ref.NO_EVIDENCE) if entries \
                    else ref.NO_EVIDENCE
            else:
                evidence = ref.best_evidence(a_domain)
            for substrate in a_domain.get("substrates") or []:
                name = (substrate.get("name") or "").strip()
                if not name:
                    continue
                smiles = substrate.get("structure") or ""
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                substrates.append({
                    "bgc_accession": accession, "module_index": index,
                    "substrate_name": name,
                    "proteinogenic": "true" if substrate.get("proteinogenic") else "false",
                    "substrate_formula": rdMolDescriptors.CalcMolFormula(mol) if mol else "",
                    "best_evidence": evidence,
                    "_canonical": Chem.MolToSmiles(mol) if mol else "",
                    "_mol": mol,
                })

    if "raw_substrate_names" not in mistakes:
        canonical_name = ref.canonical_substrate_names(substrates)
        for row in substrates:
            row["substrate_name"] = canonical_name[row["substrate_name"]]

    modules.sort(key=lambda r: (r["bgc_accession"], str(r["module_index"])))
    substrates.sort(key=lambda r: (r["bgc_accession"], str(r["module_index"]),
                                   r["substrate_name"]))
    w = ref.write_csv
    w(out / "modules.csv",
      ["bgc_accession", "module_index", "module_name", "module_type", "active",
       "n_genes", "has_a_domain", "n_modification_domains", "n_domains"], modules)
    w(out / "a_domain_substrates.csv",
      ["bgc_accession", "module_index", "substrate_name", "proteinogenic",
       "substrate_formula", "best_evidence"],
      [{k: r[k] for k in ("bgc_accession", "module_index", "substrate_name",
                          "proteinogenic", "substrate_formula", "best_evidence")}
       for r in substrates])

    def module_count(group: list[dict]) -> int:
        if "counts_rows_not_modules" in mistakes:
            return len(group)
        return len({(r["bgc_accession"], r["module_index"]) for r in group})

    by_substrate: dict[str, list[dict]] = defaultdict(list)
    for row in substrates:
        by_substrate[row["substrate_name"]].append(row)
    w(out / "substrate_counts.csv",
      ["substrate_name", "n_modules", "n_bgcs", "proteinogenic", "n_distinct_structures"],
      [{"substrate_name": name, "n_modules": module_count(group),
        "n_bgcs": len({r["bgc_accession"] for r in group}),
        "proteinogenic": (group[0]["proteinogenic"]
                          if "proteinogenic_from_first_entry" in mistakes
                          else ("true" if sum(1 for r in group
                                              if r["proteinogenic"] == "true") * 2
                                >= len(group) else "false")),
        "n_distinct_structures": len({r["_canonical"] for r in group if r["_canonical"]})}
       for name, group in sorted(by_substrate.items())])

    by_evidence: dict[str, list[dict]] = defaultdict(list)
    for row in substrates:
        by_evidence[row["best_evidence"]].append(row)
    w(out / "substrate_evidence_summary.csv",
      ["best_evidence", "n_assignments", "n_distinct_substrates", "n_bgcs"],
      [{"best_evidence": m, "n_assignments": len(g),
        "n_distinct_substrates": len({r["substrate_name"] for r in g}),
        "n_bgcs": len({r["bgc_accession"] for r in g})}
       for m, g in sorted(by_evidence.items())])

    nonpro: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in substrates:
        if row["proteinogenic"] == "false" and row["_canonical"]:
            nonpro[(row["substrate_name"], row["_canonical"])].append(row)
    w(out / "nonproteinogenic_substrates.csv",
      ["substrate_name", "canonical_smiles", "formula", "exact_mass",
       "n_heavy_atoms", "n_modules"],
      [{"substrate_name": n, "canonical_smiles": s, "formula": g[0]["substrate_formula"],
        "exact_mass": round(ExactMolWt(g[0]["_mol"]), ref.MASS_DECIMALS),
        "n_heavy_atoms": g[0]["_mol"].GetNumHeavyAtoms(), "n_modules": module_count(g)}
       for (n, s), g in sorted(nonpro.items())])

    cells: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in substrates:
        klasses = classes_of[row["bgc_accession"]]
        if "class_counted_once" in mistakes:
            klasses = klasses[:1]
        for klass in klasses:
            cells[(klass, row["substrate_name"])].append(row)
    w(out / "class_substrate_matrix.csv",
      ["biosynthesis_class", "substrate_name", "n_modules", "n_bgcs"],
      [{"biosynthesis_class": k, "substrate_name": n, "n_modules": module_count(g),
        "n_bgcs": len({r["bgc_accession"] for r in g})}
       for (k, n), g in sorted(cells.items())])

    accessions = sorted(architecture)
    a_per_bgc = Counter(m["bgc_accession"] for m in modules if m["has_a_domain"] == "true")
    arch_rows = [
        {"bgc_accession": a, "n_modules": len(architecture[a]),
         "n_a_domains": a_per_bgc.get(a, 0),
         "architecture": ">".join(architecture[a]),
         "domain_composition": ";".join(
             composition[a] if "unsorted_domain_composition" in mistakes
             else sorted(composition[a]))}
        for a in accessions
    ]
    w(out / "module_architecture.csv",
      ["bgc_accession", "n_modules", "n_a_domains", "architecture", "domain_composition"],
      arch_rows)

    by_arch: dict[str, list[dict]] = defaultdict(list)
    for row in arch_rows:
        by_arch[row["architecture"]].append(row)
    w(out / "architecture_counts.csv",
      ["architecture", "n_bgcs", "n_modules", "n_a_domains", "representative_bgc"],
      [{"architecture": a, "n_bgcs": len(g), "n_modules": g[0]["n_modules"],
        "n_a_domains": sum(r["n_a_domains"] for r in g),
        "representative_bgc": min(r["bgc_accession"] for r in g)}
       for a, g in sorted(by_arch.items())])


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: drops inactive modules": {"filters_inactive_modules"},
    "one slip: keys modules on name": {"keys_on_module_name"},
    "one slip: sorts the module list": {"sorts_modules_by_name"},
    "one slip: raw substrate names": {"raw_substrate_names"},
    "one slip: filters on the -1 location": {"filters_on_location"},
    "one slip: first evidence, not strongest": {"first_evidence_not_strongest"},
    "one slip: proteinogenic from first entry": {"proteinogenic_from_first_entry"},
    "one slip: counts rows, not modules": {"counts_rows_not_modules"},
    "one slip: class counted once": {"class_counted_once"},
    "one slip: unsorted domain composition": {"unsorted_domain_composition"},
    "hasty agent (3 slips)": {"filters_inactive_modules", "first_evidence_not_strongest",
                              "raw_substrate_names"},
    "careless agent (5 slips)": {"filters_inactive_modules", "sorts_modules_by_name",
                                 "raw_substrate_names", "counts_rows_not_modules",
                                 "class_counted_once"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    config = load_grading_config(args.task / "grading.yaml")

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':42}{'task':>6}  " + "".join(f"{n[:12]:>14}" for n in names)
    print(header)
    print("-" * len(header))
    passes = Counter()
    for label, mistakes in VARIANTS.items():
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace(",", "")
        build(args.task, out, mistakes)
        result = score_task(config, args.truth, out)
        row = f"{label:42}{result.task_score:6.2f}  "
        for artifact in result.artifacts:
            row += f"{artifact.score:14.2f}"
            if artifact.passed:
                passes[artifact.path] += 1
        print(row, flush=True)

    print()
    print(f"artifact pass rate across the {len(VARIANTS)} variants "
          f"(PROJECT.md section 8 wants 0.1 to 0.8):")
    for artifact in config.artifacts:
        rate = passes[artifact.path] / len(VARIANTS)
        print(f"  {artifact.path:34} {rate:4.2f}  "
              f"{'ok' if 0.1 <= rate <= 0.8 else 'OUTSIDE BAND'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
