"""Measured difficulty of tasks/07_npatlas_dereplication.

Synthetic-variant proxy, same caveat as the other tasks: it is not the model
calibration PROJECT.md section 8 asks for, only the cheap version of it.

The variant set deliberately spans three layers, because task 01 showed that a
set drawn only from a pipeline's parameters is blind to the parsing in front of
them: how the reference release is read, how the query set is chosen, and how
the search and the aggregates are computed.

Caching matters here. The 475 MB reference parse (~40 s) is cached per
reference configuration, the spectrum parse per peak-filter setting, and the
all-pairs network (~100 s) per spectrum set -- only one modelled mistake
changes the node universe, so the network is built at most twice.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import networkx as nx  # noqa: E402
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz  # noqa: E402
from matchms.importing import load_from_mgf  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

MISTAKES = {
    # Reading the reference.
    "dedupe_on_inchikey_block": "dedupes on the InChIKey's first block, losing stereoisomers",
    "no_dedupe": "never dedupes the release at all",
    "drops_unclassified": "drops compounds with no npclassifier block",
    "pathway_first_only": "takes the first pathway, not the sorted joined label",
    "uses_mol_weight": "searches against average mol_weight instead of exact_mass",
    "precomputed_adducts": "uses the release's own m_plus_h/m_plus_na instead of computing",
    "formula_from_smiles": "recomputes the formula from the SMILES instead of reading the field",
    # Choosing the queries.
    "includes_charge_2": "searches doubly-charged spectra as if singly charged",
    "charge_1_only": "drops the spectra whose charge the file does not state",
    "average_h_mass": "subtracts an average hydrogen, 1.00794, rather than a proton",
    "keeps_empty_spectra": "never discards the spectra with no peaks",
    # Searching.
    "tolerance_on_neutral": "sets the 10 ppm window from the neutral mass, not the m/z",
    "ppm_20": "searches at 20 ppm",
    "ppm_5": "searches at 5 ppm",
    "h_adduct_only": "forgets the sodium adduct",
    # Aggregating.
    "best_on_signed_ppm": "ranks the best match on signed, not absolute, ppm error",
    "pathway_split": "splits a multi-pathway label into one row per pathway",
    "components_numbered": "numbers families C0001.. instead of naming them by scan",
}

_ATLAS: dict[tuple, list[dict]] = {}
_SPECTRA: dict[bool, list] = {}
_EDGES: dict[int, list[dict]] = {}


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref_solution_07", task_dir / "solution" / "make_truth.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atlas_for(ref, task_dir: Path, mistakes: set[str]) -> list[dict]:
    """The reference table, cached per reference configuration."""
    key = ("dedupe_on_inchikey_block" in mistakes, "no_dedupe" in mistakes,
           "drops_unclassified" in mistakes, "pathway_first_only" in mistakes,
           "uses_mol_weight" in mistakes, "precomputed_adducts" in mistakes,
           "formula_from_smiles" in mistakes)
    if key in _ATLAS:
        return _ATLAS[key]
    print(f"    (reading the reference for {key})", flush=True)
    records = json.loads((task_dir / ref.ATLAS).read_text())

    seen: set[str] = set()
    compounds = []
    for record in sorted(records, key=lambda r: r["id"]):
        if "no_dedupe" not in mistakes:
            key_value = (record["inchikey"].split("-")[0]
                         if "dedupe_on_inchikey_block" in mistakes
                         else record["inchikey"])
            if key_value in seen:
                continue
            seen.add(key_value)
        organism = record["origin_organism"]
        classifier = record.get("npclassifier") or {}
        pathways = classifier.get("pathway_results") or []
        if not pathways and "drops_unclassified" in mistakes:
            continue
        if pathways:
            pathway = (pathways[0] if "pathway_first_only" in mistakes
                       else ";".join(sorted(set(pathways))))
        else:
            pathway = ref.UNCLASSIFIED
        mass_field = "mol_weight" if "uses_mol_weight" in mistakes else "exact_mass"
        if "formula_from_smiles" in mistakes:
            from rdkit import Chem, RDLogger
            from rdkit.Chem.rdMolDescriptors import CalcMolFormula
            RDLogger.DisableLog("rdApp.*")
            molecule = Chem.MolFromSmiles(record["smiles"])
            formula = CalcMolFormula(molecule) if molecule else record["mol_formula"]
        else:
            formula = record["mol_formula"]
        compounds.append({
            "npaid": record["npaid"],
            "mol_formula": formula,
            "exact_mass": round(float(record[mass_field]), ref.DECIMALS),
            "inchikey": record["inchikey"],
            "origin_type": organism["type"],
            "genus": organism["genus"],
            "pathway": pathway,
            # Only consulted by the precomputed_adducts variant.
            "m_plus_h": round(float(record["m_plus_h"]), ref.DECIMALS)
                        if record.get("m_plus_h") else None,
            "m_plus_na": round(float(record["m_plus_na"]), ref.DECIMALS)
                         if record.get("m_plus_na") else None,
        })
    _ATLAS[key] = compounds
    return compounds


def spectra_for(ref, task_dir: Path, mistakes: set[str]) -> list:
    """Filtered spectra, cached on whether the empty-spectrum filter ran."""
    filtered = "keeps_empty_spectra" not in mistakes
    if filtered in _SPECTRA:
        return _SPECTRA[filtered]
    print(f"    (parsing spectra; filtered={filtered})", flush=True)
    kept = []
    for spectrum in load_from_mgf(str(task_dir / ref.MGF)):
        if spectrum is None:
            continue
        if filtered and len(spectrum.peaks.mz) == 0:
            continue
        normalised = normalize_intensities(spectrum)
        if normalised is None:
            continue
        if filtered:
            trimmed = remove_peaks_around_precursor_mz(
                normalised, mz_tolerance=ref.PRECURSOR_WINDOW_DA)
            if trimmed is None or len(trimmed.peaks.mz) < ref.MIN_PEAKS:
                continue
            kept.append(trimmed)
        else:
            kept.append(normalised)
    kept.sort(key=lambda s: int(s.get("scans")))
    _SPECTRA[filtered] = kept
    return kept


def queries_for(ref, spectra: list, mistakes: set[str]) -> list[dict]:
    if "includes_charge_2" in mistakes:
        allowed = None
    elif "charge_1_only" in mistakes:
        allowed = ("1",)
    else:
        allowed = ref.QUERY_CHARGES
    rows = []
    for spectrum in spectra:
        charge = str(int(spectrum.get("charge") or 0))
        if allowed is not None and charge not in allowed:
            continue
        mz = float(spectrum.get("precursor_mz"))
        proton = 1.00794 if "average_h_mass" in mistakes else ref.PROTON_MASS
        rows.append({
            "scan": int(spectrum.get("scans")),
            "precursor_mz": round(mz, ref.DECIMALS),
            "charge": int(charge),
            "neutral_mass_h": round(mz - proton, ref.DECIMALS),
            "neutral_mass_na": round(mz - ref.SODIUM_SHIFT, ref.DECIMALS),
        })
    return rows


def annotate(ref, queries, compounds, mistakes: set[str]) -> list[dict]:
    import bisect

    ppm = 20.0 if "ppm_20" in mistakes else 5.0 if "ppm_5" in mistakes else ref.PPM_TOLERANCE
    adducts = [("[M+H]+", "neutral_mass_h", "m_plus_h")]
    if "h_adduct_only" not in mistakes:
        adducts.append(("[M+Na]+", "neutral_mass_na", "m_plus_na"))

    # The precomputed variant searches the release's own adduct m/z against the
    # measured m/z directly, rather than a neutral mass against exact_mass.
    precomputed = "precomputed_adducts" in mistakes
    rows = []
    for adduct, neutral_key, adduct_key in adducts:
        if precomputed:
            table = sorted((c[adduct_key], c["npaid"]) for c in compounds
                           if c[adduct_key] is not None)
        else:
            table = sorted((c["exact_mass"], c["npaid"]) for c in compounds)
        values = [m for m, _ in table]
        by_id = {c["npaid"]: c for c in compounds}
        for query in queries:
            observed = query["precursor_mz"] if precomputed else query[neutral_key]
            base = observed if "tolerance_on_neutral" in mistakes else query["precursor_mz"]
            window = base * ppm / 1e6
            lo = bisect.bisect_left(values, observed - window)
            hi = bisect.bisect_right(values, observed + window)
            for exact, npaid in table[lo:hi]:
                rows.append({
                    "scan": query["scan"],
                    "npaid": npaid,
                    "adduct": adduct,
                    "ppm_error": round((observed - exact) / exact * 1e6, 2),
                    "mol_formula": by_id[npaid]["mol_formula"],
                })
    rows.sort(key=lambda r: (r["scan"], r["adduct"], r["npaid"]))
    return rows


def best_for(annotations, compounds, mistakes: set[str]) -> list[dict]:
    by_id = {c["npaid"]: c for c in compounds}
    best: dict[int, dict] = {}
    for row in annotations:
        if "best_on_signed_ppm" in mistakes:
            key = (row["ppm_error"], row["npaid"])
            current = best.get(row["scan"])
            better = current is None or key < (current["ppm_error"], current["npaid"])
        else:
            key = (abs(row["ppm_error"]), row["npaid"])
            current = best.get(row["scan"])
            better = current is None or key < (abs(current["ppm_error"]), current["npaid"])
        if better:
            best[row["scan"]] = row
    rows = []
    for scan in sorted(best):
        row = best[scan]
        compound = by_id[row["npaid"]]
        rows.append({"scan": scan, "npaid": row["npaid"], "adduct": row["adduct"],
                     "ppm_error": row["ppm_error"], "mol_formula": compound["mol_formula"],
                     "genus": compound["genus"], "pathway": compound["pathway"]})
    return rows


def edges_for(ref, spectra: list) -> list[dict]:
    """The network, cached per distinct set of spectra that can carry an edge.

    An edge needs MIN_MATCHED_PEAKS matched peaks, so a spectrum with fewer
    peaks than that can never form one. Restricting the all-pairs comparison to
    spectra that clear the bar gives an identical edge list and is the
    difference between minutes and hours: the variant that keeps the 96,182
    empty spectra adds isolated nodes and no edges. The full list still goes to
    the component builder, so those nodes exist there and are dropped by the
    two-node rule.
    """
    eligible = [s for s in spectra if len(s.peaks.mz) >= ref.MIN_MATCHED_PEAKS]
    key = len(eligible)
    if key in _EDGES:
        return _EDGES[key]
    print(f"    (networking {key} of {len(spectra)} spectra)", flush=True)
    edges = ref.apply_top_k(ref.candidate_edges(eligible))
    _EDGES[key] = edges
    return edges


def build(task_dir: Path, out: Path, mistakes: set[str]) -> None:
    ref = load_solution(task_dir)
    compounds = atlas_for(ref, task_dir, mistakes)
    spectra = spectra_for(ref, task_dir, mistakes)
    queries = queries_for(ref, spectra, mistakes)
    annotations = annotate(ref, queries, compounds, mistakes)
    best = best_for(annotations, compounds, mistakes)
    w = ref.write_csv

    w(out / "atlas_compounds.csv",
      ["npaid", "mol_formula", "exact_mass", "origin_type", "genus", "pathway"], compounds)
    w(out / "queries.csv",
      ["scan", "precursor_mz", "charge", "neutral_mass_h", "neutral_mass_na"], queries)
    w(out / "annotations.csv",
      ["scan", "npaid", "adduct", "ppm_error", "mol_formula"], annotations)
    w(out / "best_annotation.csv",
      ["scan", "npaid", "adduct", "ppm_error", "mol_formula", "genus", "pathway"], best)

    # The variant that counted only the matched compounds per formula is gone,
    # along with the column it got wrong: a formula fixes the monoisotopic mass,
    # so every compound sharing a matched formula is matched too, and the column
    # was provably equal to its neighbour.
    w(out / "formula_candidates.csv", ["mol_formula", "n_compounds", "n_scans"],
      ref.formula_candidates(annotations, compounds))

    w(out / "genus_annotations.csv", ["genus", "n_compounds_matched", "n_scans"],
      ref.genus_annotations(annotations, compounds))

    if "pathway_split" in mistakes:
        by_id = {c["npaid"]: c for c in compounds}
        scans, hits = defaultdict(set), defaultdict(set)
        for row in annotations:
            for part in by_id[row["npaid"]]["pathway"].split(";"):
                scans[part].add(row["scan"])
                hits[part].add(row["npaid"])
        pathway_rows = [{"pathway": p, "n_compounds_matched": len(hits[p]),
                         "n_scans": len(scans[p])} for p in sorted(scans)]
    else:
        pathway_rows = ref.pathway_counts(annotations, compounds)
    w(out / "pathway_counts.csv", ["pathway", "n_compounds_matched", "n_scans"], pathway_rows)

    edges = edges_for(ref, spectra)
    rows = ref.component_annotations(spectra, edges, best, annotations)
    if "components_numbered" in mistakes:
        for index, row in enumerate(sorted(rows, key=lambda r: -r["n_nodes"]), start=1):
            row["component_id"] = f"C{index:04d}"
        rows.sort(key=lambda r: r["component_id"])
    w(out / "component_annotations.csv",
      ["component_id", "n_nodes", "n_annotated", "n_distinct_formulas", "top_pathway"], rows)


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: dedupe on InChIKey block": {"dedupe_on_inchikey_block"},
    "one slip: no dedupe at all": {"no_dedupe"},
    "one slip: drops the unclassified": {"drops_unclassified"},
    "one slip: first pathway only": {"pathway_first_only"},
    "one slip: average mass, not exact": {"uses_mol_weight"},
    "one slip: release's own adduct m/z": {"precomputed_adducts"},
    "one slip: searches charge 2 too": {"includes_charge_2"},
    "one slip: charge 1 only": {"charge_1_only"},
    "one slip: formula from the SMILES": {"formula_from_smiles"},
    "one slip: average H, not a proton": {"average_h_mass"},
    "one slip: keeps the empty spectra": {"keeps_empty_spectra"},
    "one slip: window from neutral mass": {"tolerance_on_neutral"},
    "one slip: 20 ppm": {"ppm_20"},
    "one slip: 5 ppm": {"ppm_5"},
    "one slip: forgets [M+Na]+": {"h_adduct_only"},
    "one slip: best match on signed ppm": {"best_on_signed_ppm"},
    "one slip: splits pathway labels": {"pathway_split"},
    "one slip: families numbered": {"components_numbered"},
    "hasty agent (3 slips)": {"charge_1_only", "h_adduct_only", "pathway_split"},
    "careless agent (5 slips)": {"dedupe_on_inchikey_block", "uses_mol_weight",
                                 "includes_charge_2", "ppm_20", "components_numbered"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    config = load_grading_config(args.task / "grading.yaml")

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':40}{'task':>6}  " + "".join(f"{n[:13]:>15}" for n in names)
    print(header)
    print("-" * len(header))
    passes = Counter()
    for label, mistakes in VARIANTS.items():
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace("/", "")
        build(args.task, out, mistakes)
        result = score_task(config, args.truth, out)
        scored = {a.path: a for a in result.artifacts + result.soft_artifacts}
        row = f"{label:40}{result.task_score:6.2f}  "
        for spec in config.artifacts:
            artifact = scored[spec.path]
            row += f"{artifact.score:15.2f}"
            if artifact.passed:
                passes[artifact.path] += 1
        print(row, flush=True)

    print()
    print(f"artifact pass rate across the {len(VARIANTS)} variants "
          f"(PROJECT.md section 8 wants 0.1 to 0.8):")
    for artifact in config.artifacts:
        rate = passes[artifact.path] / len(VARIANTS)
        verdict = ("soft, out of the headline" if artifact.soft
                   else "ok" if 0.1 <= rate <= 0.8 else "OUTSIDE BAND")
        print(f"  {artifact.path:34} {rate:4.2f}  {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
