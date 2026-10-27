"""Measure task 03's difficulty by grading plausible wrong pipelines.

Same approach and same caveat as scripts/difficulty_task02.py: each variant is
a mistake a competent agent actually makes, applied to the same MIBiG release
the reference solution reads. It is a cheap proxy for the calibration
PROJECT.md §8 asks for, useful mainly for finding artifacts nothing can fail.

    python scripts/difficulty_task03.py --task tasks/03_mibig_cheminformatics \
        --truth tasks/03_mibig_cheminformatics/truth --scratch <dir>
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, DataStructs, RDLogger  # noqa: E402
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors  # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

RDLogger.DisableLog("rdApp.*")

MISTAKES = {
    "trusts_stated_formula": "copies MIBiG's formula field instead of computing it",
    "no_dedupe": "assumes (accession, name) is unique",
    "dedupe_on_raw_name": "deduplicates on the raw name, missing the case variants",
    "keyed_on_entries": "keys structure artifacts on MIBiG entries, not distinct structures",
    "strips_stereo_characters": "removes @ and slashes from the SMILES text to flatten it",
    "pairwise_not_transitive": "treats each similar pair as a family, no connected components",
    "families_numbered_by_rank": "names families FAM_0001.. by size rank instead of by member",
    "morgan_radius_3": "fingerprints at radius 3",
    "drops_acyclic": "drops compounds whose Murcko scaffold is empty",
    "class_counted_per_bgc": "counts a structure once per BGC rather than once per class",
}


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref03", task_dir / "solution" / "make_truth.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_RECORDS: list[dict] = []


def build(task_dir: Path, out: Path, mistakes: set[str]) -> None:
    ref = load_solution(task_dir)
    if not _RECORDS:
        _RECORDS.extend(ref.read_records(task_dir / ref.ARCHIVE))

    best: dict[tuple[str, str], dict] = {}
    classes_of: dict[str, list[str]] = {}
    rows: list[dict] = []
    for record in _RECORDS:
        accession = record["accession"]
        classes_of[accession] = sorted(
            {e["class"] for e in (record.get("biosynthesis") or {}).get("classes") or []
             if e.get("class")}
        )
        for compound in record.get("compounds") or []:
            smiles, name = compound.get("structure"), (compound.get("name") or "").strip()
            if not smiles or not name:
                continue
            parsed = ref.canonical(smiles)
            if parsed is None:
                continue
            mol, canonical_smiles = parsed
            stated = (compound.get("formula") or "").strip()
            computed = rdMolDescriptors.CalcMolFormula(mol)
            row = {
                "bgc_accession": accession,
                "compound_name": name,
                "canonical_smiles": canonical_smiles,
                "formula": stated if ("trusts_stated_formula" in mistakes and stated)
                else computed,
                "exact_mass": round(ExactMolWt(mol), ref.MASS_DECIMALS),
                "n_heavy_atoms": mol.GetNumHeavyAtoms(),
                "n_rings": rdMolDescriptors.CalcNumRings(mol),
                "stated_formula": stated,
                "computed_formula": stated if "trusts_stated_formula" in mistakes
                else computed,
                "_mol": mol,
            }
            if "no_dedupe" in mistakes:
                rows.append(row)
                continue
            key = (accession, name) if "dedupe_on_raw_name" in mistakes \
                else (accession, name.casefold())
            if key in best and best[key]["canonical_smiles"] <= canonical_smiles:
                continue
            best[key] = row
    if "no_dedupe" not in mistakes:
        rows = [best[k] for k in sorted(best)]
    rows.sort(key=lambda r: (r["bgc_accession"], r["compound_name"]))

    if "drops_acyclic" in mistakes:
        rows = [r for r in rows if ref.scaffold_of(r["_mol"])[0]]

    w = ref.write_csv
    w(out / "compounds.csv",
      ["bgc_accession", "compound_name", "canonical_smiles", "formula", "exact_mass",
       "n_heavy_atoms", "n_rings"],
      [{k: r[k] for k in ("bgc_accession", "compound_name", "canonical_smiles", "formula",
                          "exact_mass", "n_heavy_atoms", "n_rings")} for r in rows])
    w(out / "formula_disagreements.csv",
      ["bgc_accession", "compound_name", "stated_formula", "computed_formula"],
      [{"bgc_accession": r["bgc_accession"], "compound_name": r["compound_name"],
        "stated_formula": r["stated_formula"], "computed_formula": r["computed_formula"]}
       for r in rows if r["stated_formula"] and r["stated_formula"] != r["computed_formula"]])

    structures: dict[str, Chem.Mol] = {}
    bgcs_of: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        structures.setdefault(r["canonical_smiles"], r["_mol"])
        bgcs_of[r["canonical_smiles"]].add(r["bgc_accession"])

    # `units` is what the agent treats as one thing. Collapsing to distinct
    # structures is the method; leaving the MIBiG entries uncollapsed inflates
    # every count downstream, which is the realistic shape of this mistake.
    if "keyed_on_entries" in mistakes:
        units = [r["canonical_smiles"] for r in rows]
    else:
        units = sorted(structures)
    smiles_list = sorted(set(units))
    scaffolds = {s: ref.scaffold_of(structures[s]) for s in smiles_list}

    w(out / "scaffolds.csv",
      ["canonical_smiles", "scaffold_smiles", "n_scaffold_rings", "n_scaffold_atoms"],
      [{"canonical_smiles": s, "scaffold_smiles": scaffolds[s][0],
        "n_scaffold_rings": scaffolds[s][1], "n_scaffold_atoms": scaffolds[s][2]}
       for s in sorted(units)])

    per_scaffold: dict[str, list[str]] = defaultdict(list)
    scaffold_bgcs: dict[str, set[str]] = defaultdict(set)
    scaffold_classes: dict[str, set[str]] = defaultdict(set)
    for s in units:
        key = scaffolds[s][0]
        per_scaffold[key].append(s)
        scaffold_bgcs[key] |= bgcs_of[s]
        for accession in bgcs_of[s]:
            scaffold_classes[key] |= set(classes_of[accession])
    w(out / "scaffold_counts.csv",
      ["scaffold_smiles", "n_structures", "n_bgcs", "n_classes", "representative_smiles"],
      [{"scaffold_smiles": k, "n_structures": len(v), "n_bgcs": len(scaffold_bgcs[k]),
        "n_classes": len(scaffold_classes[k]), "representative_smiles": min(v)}
       for k, v in sorted(per_scaffold.items())])

    by_formula: dict[str, list[str]] = defaultdict(list)
    for s in units:
        by_formula[rdMolDescriptors.CalcMolFormula(structures[s])].append(s)
    w(out / "formula_collisions.csv",
      ["formula", "n_structures", "n_bgcs", "representative_smiles"],
      [{"formula": f, "n_structures": len(g),
        "n_bgcs": len({b for s in g for b in bgcs_of[s]}),
        "representative_smiles": min(g)}
       for f, g in sorted(by_formula.items()) if len(set(g)) > 1])

    def flat(smiles: str) -> str:
        if "strips_stereo_characters" in mistakes:
            # Textual stripping: plausible, and wrong. It does not recanonicalise,
            # so the result is often not even the canonical form of the flat molecule.
            return smiles.replace("@", "").replace("/", "").replace("\\", "")
        return ref.flatten(structures[smiles])

    by_flat: dict[str, set[str]] = defaultdict(set)
    for s in units:
        by_flat[flat(s)].add(s)
    w(out / "stereo_collapse_groups.csv", ["flat_smiles", "canonical_smiles", "group_size"],
      [{"flat_smiles": f, "canonical_smiles": s, "group_size": len(g)}
       for f, g in sorted(by_flat.items()) if len(g) > 1 for s in sorted(g)])

    radius = 3 if "morgan_radius_3" in mistakes else ref.MORGAN_RADIUS
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=ref.MORGAN_BITS)
    fps = [generator.GetFingerprint(structures[s]) for s in smiles_list]
    neighbours: dict[int, set[int]] = defaultdict(set)
    for i in range(len(fps)):
        for offset, value in enumerate(
            DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        ):
            if value >= ref.TANIMOTO_THRESHOLD:
                neighbours[i].add(i + 1 + offset)

    if "pairwise_not_transitive" in mistakes:
        families = {s: s for s in smiles_list}
        for i, group in neighbours.items():
            for j in group:
                families[smiles_list[j]] = min(families[smiles_list[j]], smiles_list[i])
    else:
        families = ref.single_linkage(smiles_list, neighbours)
    if "families_numbered_by_rank" in mistakes:
        counts = Counter(families.values())
        order = sorted(counts, key=lambda f: (-counts[f], f))
        rename = {f: f"FAM_{i + 1:04d}" for i, f in enumerate(order)}
        families = {s: rename[f] for s, f in families.items()}
    size = Counter(families[s] for s in units)
    w(out / "similarity_families.csv", ["canonical_smiles", "family_id", "family_size"],
      [{"canonical_smiles": s, "family_id": families[s], "family_size": size[families[s]]}
       for s in sorted(units)])

    # One append per (unit, distinct class). Appending once per reporting BGC
    # double-counts a structure that several BGCs report, which made the
    # no-mistake baseline score 0.99 instead of 1.00 and so invalidated this
    # column: a proxy measured against a wrong reference measures nothing.
    cells: dict[tuple[str, int], list[str]] = defaultdict(list)
    cell_bgcs: dict[tuple[str, int], set[str]] = defaultdict(set)
    for s in units:
        rings = scaffolds[s][1]
        accessions = sorted(bgcs_of[s])
        if "class_counted_per_bgc" in mistakes:
            accessions = accessions[:1]
        klasses = sorted({k for a in accessions for k in classes_of[a]})
        for klass in klasses:
            cells[(klass, rings)].append(s)
        for accession in accessions:
            for klass in classes_of[accession]:
                cell_bgcs[(klass, rings)].add(accession)
    w(out / "class_scaffold_summary.csv",
      ["biosynthesis_class", "n_scaffold_rings", "n_structures", "n_bgcs",
       "representative_smiles"],
      [{"biosynthesis_class": k, "n_scaffold_rings": r, "n_structures": len(g),
        "n_bgcs": len(cell_bgcs[(k, r)]), "representative_smiles": min(g)}
       for (k, r), g in sorted(cells.items())])


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: trusts MIBiG's formula field": {"trusts_stated_formula"},
    "one slip: no dedupe at all": {"no_dedupe"},
    "one slip: dedupes on the raw name": {"dedupe_on_raw_name"},
    "one slip: keyed on entries not structures": {"keyed_on_entries"},
    "one slip: strips stereo characters": {"strips_stereo_characters"},
    "one slip: no transitive closure": {"pairwise_not_transitive"},
    "one slip: families numbered by rank": {"families_numbered_by_rank"},
    "one slip: Morgan radius 3": {"morgan_radius_3"},
    "one slip: drops acyclic compounds": {"drops_acyclic"},
    "one slip: class counted per BGC": {"class_counted_per_bgc"},
    "hasty agent (3 slips)": {"trusts_stated_formula", "strips_stereo_characters",
                              "dedupe_on_raw_name"},
    "careless agent (5 slips)": {"trusts_stated_formula", "no_dedupe", "morgan_radius_3",
                                 "families_numbered_by_rank", "drops_acyclic"},
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
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace("'", "")
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
        print(f"  {artifact.path:32} {rate:4.2f}  "
              f"{'ok' if 0.1 <= rate <= 0.8 else 'OUTSIDE BAND'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
