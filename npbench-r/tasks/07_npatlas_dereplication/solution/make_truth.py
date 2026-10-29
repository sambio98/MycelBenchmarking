"""Reference solution for tasks/07_npatlas_dereplication.

Exact-mass dereplication of the Crusemann molecular network against a frozen
NP Atlas release, PROJECT.md §4 category 4.

    python solution/make_truth.py --task tasks/07_npatlas_dereplication --out <dir>

Why exact mass and not SIRIUS or CSI:FingerID. PROJECT.md §9 rules both out:
they call a web service whose model version determines the output, and SIRIUS
4, 5.8.x and 6.x disagree with each other. A GNPS library hit list is out for
the same reason -- the library is continuously recurated. What is left, and
what §9 names as safe, is matching against a database snapshot frozen inside
the bundle. So the reference is an NP Atlas release pinned by checksum, and
every number here is arithmetic over it.

The ambiguity is the point. Dereplication by mass alone cannot identify a
compound, and this task does not pretend otherwise: 2,925 of the 36,422
compounds share a stereochemistry-stripped skeleton with another, hundreds of
molecular formulae are shared by several compounds, and a 10 ppm window at
m/z 400 is 4 mDa wide. The artifacts are built to measure whether an agent
reports that ambiguity honestly or collapses it to one confident answer.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import matchms
import networkx as nx
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosineGreedy

MGF = "data/MSV000078836_specs_ms.mgf"
ATLAS = "data/NPAtlas_download.json"

MATCHMS_VERSION = "0.33.1"

# NP Atlas serves its current release from an unversioned URL, so the release
# cannot be read back out of the file. The manifest pins the bytes by sha256
# and this count is the corroborating fingerprint: v2024_09 has exactly this
# many records, and any other release will trip the check rather than silently
# producing different truth.
ATLAS_RELEASE = "2024_09"
ATLAS_RECORDS = 36454

# Networking parameters, identical to task 01. The network is rebuilt here
# because the question "which molecular families are unknown" is the point of
# dereplication, and that needs the families.
PRECURSOR_WINDOW_DA = 17.0
MIN_PEAKS = 6
FRAGMENT_TOLERANCE = 0.5
MIN_COSINE = 0.70
MIN_MATCHED_PEAKS = 6
MAX_PRECURSOR_DELTA = 2.0
TOP_K = 10

# Dereplication parameters.
PROTON_MASS = 1.007276      # [M+H]+
SODIUM_SHIFT = 22.989218    # [M+Na]+, sodium minus one electron
PPM_TOLERANCE = 10.0
QUERY_CHARGES = ("0", "1")  # 0 means the file states none; both are read as 1+
UNCLASSIFIED = "Unclassified"
DECIMALS = 4


def check_version() -> None:
    if matchms.__version__ != MATCHMS_VERSION:
        raise SystemExit(
            f"this task's truth was produced by matchms {MATCHMS_VERSION}, found "
            f"{matchms.__version__}. The similarity class was renamed between "
            f"releases and cosine scores are implementation-dependent."
        )


# ------------------------------------------------------------------ the reference


def load_atlas(path: Path) -> list[dict]:
    """The NP Atlas release, deduplicated on full InChIKey.

    32 records are a second name for a structure already in the release. The
    dedupe keeps the lowest `id`, which is the lowest NPAID, so the survivor
    does not depend on dict or file ordering.

    Deduplicating on the InChIKey's *first block* instead would be wrong here
    and is a modelled mistake: that block ignores stereochemistry, and 2,925
    compounds would lose a distinct stereoisomer to it.
    """
    records = json.loads(path.read_text())
    if len(records) != ATLAS_RECORDS:
        raise SystemExit(
            f"expected {ATLAS_RECORDS} records for NP Atlas {ATLAS_RELEASE}, found "
            f"{len(records)}. The download URL is unversioned, so a newer release "
            f"replaces it in place; refetch against the manifest checksum."
        )

    seen: set[str] = set()
    compounds = []
    for record in sorted(records, key=lambda r: r["id"]):
        key = record["inchikey"]
        if key in seen:
            continue
        seen.add(key)
        organism = record["origin_organism"]
        classifier = record.get("npclassifier") or {}
        pathways = classifier.get("pathway_results") or []
        compounds.append({
            "npaid": record["npaid"],
            "mol_formula": record["mol_formula"],
            "exact_mass": round(float(record["exact_mass"]), DECIMALS),
            "inchikey": key,
            "origin_type": organism["type"],
            "genus": organism["genus"],
            # A compound with no npclassifier block, or an empty pathway list,
            # is reported as Unclassified rather than dropped: 3,923 of the
            # release have none, and dropping them would silently shrink every
            # downstream count.
            "pathway": ";".join(sorted(set(pathways))) if pathways else UNCLASSIFIED,
        })
    return compounds


# ------------------------------------------------------------------- the queries


def load_and_filter(path: Path) -> tuple[list, int]:
    """Retained spectra sorted by scan, plus how many were loaded.

    The same filter as task 01: drop the empty spectra, normalise, remove peaks
    within 17 Da of the precursor, require six surviving peaks.
    """
    loaded = list(load_from_mgf(str(path)))
    retained = []
    for spectrum in loaded:
        if spectrum is None or len(spectrum.peaks.mz) == 0:
            continue
        normalised = normalize_intensities(spectrum)
        if normalised is None:
            continue
        trimmed = remove_peaks_around_precursor_mz(
            normalised, mz_tolerance=PRECURSOR_WINDOW_DA
        )
        if trimmed is None or len(trimmed.peaks.mz) < MIN_PEAKS:
            continue
        retained.append(trimmed)
    retained.sort(key=lambda s: int(s.get("scans")))
    return retained, len(loaded)


def query_rows(spectra: list) -> list[dict]:
    """One row per spectrum eligible for mass search.

    Charge 2 is excluded: its precursor m/z does not give the neutral mass
    under a singly-charged assumption, and guessing a charge the file states is
    worse than declining to search. Charge 0 means the file states none, and is
    read as singly charged -- which is the convention, not a certainty, and is
    why `queries.csv` keeps the charge column.
    """
    rows = []
    for spectrum in spectra:
        charge = str(int(spectrum.get("charge") or 0))
        if charge not in QUERY_CHARGES:
            continue
        mz = float(spectrum.get("precursor_mz"))
        rows.append({
            "scan": int(spectrum.get("scans")),
            "precursor_mz": round(mz, DECIMALS),
            "charge": int(charge),
            "neutral_mass_h": round(mz - PROTON_MASS, DECIMALS),
            "neutral_mass_na": round(mz - SODIUM_SHIFT, DECIMALS),
        })
    return rows


# --------------------------------------------------------------- the annotations


def annotate(queries: list[dict], compounds: list[dict]) -> list[dict]:
    """Every compound within 10 ppm of a query's neutral mass, for both adducts.

    The window is 10 ppm **of the measured precursor m/z**, so it is one width
    per query rather than one per candidate; the reported error is relative to
    the candidate's exact mass. Stating both is the point -- a tolerance
    defined on the neutral mass instead gives a slightly different window and
    is a modelled mistake.

    Candidates are found by bisection over a sorted mass list, which is exact.
    A linear scan over 36,422 compounds per query would give the same answer
    and is not worth the 139 million comparisons.
    """
    masses = sorted((c["exact_mass"], c["npaid"]) for c in compounds)
    mass_only = [m for m, _ in masses]
    by_id = {c["npaid"]: c for c in compounds}

    rows = []
    for query in queries:
        window = query["precursor_mz"] * PPM_TOLERANCE / 1e6
        for adduct, neutral in (("[M+H]+", query["neutral_mass_h"]),
                                ("[M+Na]+", query["neutral_mass_na"])):
            lo = bisect.bisect_left(mass_only, neutral - window)
            hi = bisect.bisect_right(mass_only, neutral + window)
            for exact, npaid in masses[lo:hi]:
                compound = by_id[npaid]
                rows.append({
                    "scan": query["scan"],
                    "npaid": npaid,
                    "adduct": adduct,
                    "ppm_error": round((neutral - exact) / exact * 1e6, 2),
                    "mol_formula": compound["mol_formula"],
                })
    rows.sort(key=lambda r: (r["scan"], r["adduct"], r["npaid"]))
    return rows


def best_annotations(annotations: list[dict], compounds: list[dict]) -> list[dict]:
    """One row per annotated scan: the nearest candidate.

    Ranked on absolute ppm error, then NPAID ascending. The second key is
    load-bearing: masses are rounded to four decimals, so exact ties on
    ppm_error are common and a sort without a tie-break is not reproducible.
    """
    by_id = {c["npaid"]: c for c in compounds}
    best: dict[int, dict] = {}
    for row in annotations:
        key = (abs(row["ppm_error"]), row["npaid"])
        current = best.get(row["scan"])
        if current is None or key < (abs(current["ppm_error"]), current["npaid"]):
            best[row["scan"]] = row
    rows = []
    for scan in sorted(best):
        row = best[scan]
        compound = by_id[row["npaid"]]
        rows.append({
            "scan": scan,
            "npaid": row["npaid"],
            "adduct": row["adduct"],
            "ppm_error": row["ppm_error"],
            "mol_formula": compound["mol_formula"],
            "genus": compound["genus"],
            "pathway": compound["pathway"],
        })
    return rows


# ----------------------------------------------------------------- the aggregates


def formula_candidates(annotations: list[dict], compounds: list[dict]) -> list[dict]:
    """Per molecular formula matched: how ambiguous that match is.

    `n_compounds` counts the distinct NP Atlas compounds carrying the formula
    in the whole release, which is the size of the ambiguity an annotation
    inherits. `n_scans` counts the queries it was offered to.

    An earlier version also reported how many of those compounds were matched,
    and a difficulty sweep showed that column could not be got wrong. The
    reason is arithmetic, not luck: a molecular formula fixes the monoisotopic
    mass, so every compound sharing a matched formula sits at the same mass and
    is inside the same window. Measured on the release, all 1,538 rows agreed.
    A column that is provably equal to its neighbour is not a test, so it is
    gone.
    """
    release_total = Counter(c["mol_formula"] for c in compounds)
    matched_scans = defaultdict(set)
    matched_formulas = set()
    for row in annotations:
        matched_formulas.add(row["mol_formula"])
        matched_scans[row["mol_formula"]].add(row["scan"])
    return [
        {
            "mol_formula": formula,
            "n_compounds": release_total[formula],
            "n_scans": len(matched_scans[formula]),
        }
        for formula in sorted(matched_formulas)
    ]


def genus_annotations(annotations: list[dict], compounds: list[dict]) -> list[dict]:
    by_id = {c["npaid"]: c for c in compounds}
    scans = defaultdict(set)
    hits = defaultdict(set)
    for row in annotations:
        genus = by_id[row["npaid"]]["genus"]
        scans[genus].add(row["scan"])
        hits[genus].add(row["npaid"])
    return [
        {"genus": genus, "n_compounds_matched": len(hits[genus]),
         "n_scans": len(scans[genus])}
        for genus in sorted(scans)
    ]


def pathway_counts(annotations: list[dict], compounds: list[dict]) -> list[dict]:
    """Counts per npclassifier pathway label, as the label appears.

    A compound carrying two pathways has one combined label, joined on ';' and
    sorted. Splitting the label and counting the compound once per pathway is a
    different question and a modelled mistake.
    """
    by_id = {c["npaid"]: c for c in compounds}
    scans = defaultdict(set)
    hits = defaultdict(set)
    for row in annotations:
        pathway = by_id[row["npaid"]]["pathway"]
        scans[pathway].add(row["scan"])
        hits[pathway].add(row["npaid"])
    return [
        {"pathway": pathway, "n_compounds_matched": len(hits[pathway]),
         "n_scans": len(scans[pathway])}
        for pathway in sorted(scans)
    ]


# -------------------------------------------------------------------- the network


def candidate_edges(spectra: list) -> list[dict]:
    similarity = ModifiedCosineGreedy(tolerance=FRAGMENT_TOLERANCE)
    scans = [int(s.get("scans")) for s in spectra]
    masses = [float(s.get("precursor_mz")) for s in spectra]
    edges = []
    for i in range(len(spectra)):
        for j in range(i + 1, len(spectra)):
            delta = masses[j] - masses[i]
            if abs(delta) > MAX_PRECURSOR_DELTA:
                continue
            result = similarity.pair(spectra[i], spectra[j])
            score = float(result["score"])
            matches = int(result["matches"])
            if score >= MIN_COSINE and matches >= MIN_MATCHED_PEAKS:
                a, b = (scans[i], scans[j]) if scans[i] < scans[j] else (scans[j], scans[i])
                edges.append({"scan_a": a, "scan_b": b,
                              "cosine": round(score, DECIMALS)})
    return edges


def apply_top_k(edges: list[dict]) -> list[dict]:
    neighbours: dict[int, list] = defaultdict(list)
    for edge in edges:
        neighbours[edge["scan_a"]].append((-edge["cosine"], edge["scan_b"], edge))
        neighbours[edge["scan_b"]].append((-edge["cosine"], edge["scan_a"], edge))
    kept_by = {}
    for scan, group in neighbours.items():
        group.sort(key=lambda item: (item[0], item[1]))
        kept_by[scan] = {(e["scan_a"], e["scan_b"]) for _, _, e in group[:TOP_K]}
    return [e for e in edges
            if (e["scan_a"], e["scan_b"]) in kept_by[e["scan_a"]]
            and (e["scan_a"], e["scan_b"]) in kept_by[e["scan_b"]]]


def component_annotations(
    spectra: list, edges: list[dict], best: list[dict], annotations: list[dict]
) -> list[dict]:
    """Per molecular family with two or more nodes: how much of it is known.

    Named by its smallest member scan, as in task 01, so the identifier does
    not depend on an ordering the agent chose.

    This is the artifact dereplication exists for. A family where no node has a
    candidate is the interesting one, and `n_annotated` of 0 is a real and
    common answer, not a failure.
    """
    graph = nx.Graph()
    graph.add_nodes_from(int(s.get("scans")) for s in spectra)
    graph.add_edges_from((e["scan_a"], e["scan_b"]) for e in edges)

    best_by_scan = {row["scan"]: row for row in best}
    formulas_by_scan = defaultdict(set)
    for row in annotations:
        formulas_by_scan[row["scan"]].add(row["mol_formula"])

    rows = []
    for component in nx.connected_components(graph):
        if len(component) < 2:
            continue
        members = sorted(component)
        annotated = [s for s in members if s in best_by_scan]
        formulas = set()
        for scan in members:
            formulas |= formulas_by_scan[scan]
        # The family's pathway label is the most common among its annotated
        # nodes, ties broken alphabetically so the answer is reproducible.
        labels = Counter(best_by_scan[s]["pathway"] for s in annotated)
        top = min(sorted(labels), key=lambda p: (-labels[p], p)) if labels else UNCLASSIFIED
        rows.append({
            "component_id": f"C{members[0]}",
            "n_nodes": len(members),
            "n_annotated": len(annotated),
            "n_distinct_formulas": len(formulas),
            "top_pathway": top,
        })
    rows.sort(key=lambda r: r["component_id"])
    return rows


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    check_version()

    compounds = load_atlas(args.task / ATLAS)
    print(f"NP Atlas {ATLAS_RELEASE}: {len(compounds)} compounds after dedupe", flush=True)
    write_csv(args.out / "atlas_compounds.csv",
              ["npaid", "mol_formula", "exact_mass", "origin_type", "genus", "pathway"],
              compounds)

    spectra, n_loaded = load_and_filter(args.task / MGF)
    print(f"{n_loaded} spectra loaded, {len(spectra)} retained", flush=True)
    queries = query_rows(spectra)
    print(f"{len(queries)} eligible for mass search", flush=True)
    write_csv(args.out / "queries.csv",
              ["scan", "precursor_mz", "charge", "neutral_mass_h", "neutral_mass_na"],
              queries)

    annotations = annotate(queries, compounds)
    print(f"{len(annotations)} candidate annotations", flush=True)
    write_csv(args.out / "annotations.csv",
              ["scan", "npaid", "adduct", "ppm_error", "mol_formula"], annotations)

    best = best_annotations(annotations, compounds)
    print(f"{len(best)} scans with at least one candidate", flush=True)
    write_csv(args.out / "best_annotation.csv",
              ["scan", "npaid", "adduct", "ppm_error", "mol_formula", "genus", "pathway"],
              best)

    write_csv(args.out / "formula_candidates.csv",
              ["mol_formula", "n_compounds", "n_scans"],
              formula_candidates(annotations, compounds))
    write_csv(args.out / "genus_annotations.csv",
              ["genus", "n_compounds_matched", "n_scans"],
              genus_annotations(annotations, compounds))
    write_csv(args.out / "pathway_counts.csv",
              ["pathway", "n_compounds_matched", "n_scans"],
              pathway_counts(annotations, compounds))

    edges = apply_top_k(candidate_edges(spectra))
    print(f"{len(edges)} edges after top-{TOP_K}", flush=True)
    write_csv(args.out / "component_annotations.csv",
              ["component_id", "n_nodes", "n_annotated", "n_distinct_formulas",
               "top_pathway"],
              component_annotations(spectra, edges, best, annotations))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
