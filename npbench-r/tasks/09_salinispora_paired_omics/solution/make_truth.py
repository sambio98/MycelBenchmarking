"""Stage 3 of the reference solution for tasks/09_salinispora_paired_omics.

Metabologenomic link ranking across 111 Salinispora strains that have both a
sequenced genome and an LC-MS/MS run, PROJECT.md §4 category 5 -- the paired
omics category, and the one §11 calls the differentiator worth building.

    python solution/make_truth.py --task tasks/09_salinispora_paired_omics --out <dir>

Run after stages 1 and 2 (run_gecco.py, extract_features.py).

How the pairing was possible at all. The clustered .mgf this deposit ships
records `FILENAME=specs_ms.pklbin` for all 100,263 of its spectra, so it cannot
say which strain a consensus spectrum came from, and MassIVE's own listing
endpoints are unreachable from a current certificate store -- its chain roots
in an emSign CA absent from both the system store and certifi. The Paired Omics
Data Platform curates exactly the missing mapping: 360 genome-to-file links for
this project, of which 111 strains also have an assembly staged for task 02.
strains.tsv is that join, built by scripts/build_task09_strains.py.

The three halves:

  * **Genomes.** GECCO over each assembly, then gene cluster families by
    single linkage over distinct-Pfam-domain Jaccard -- the same rule as task
    02, so a reader can compare the two inventories directly.
  * **Metabolomes.** Each strain's filtered MS2 spectra, from stage 2, matched
    against the deposit's own GNPS-clustered consensus spectra. A consensus
    spectrum is *present* in a strain when that strain has an MS2 scan within
    0.5 Da of its precursor whose modified cosine reaches 0.70. This is the
    reproducible way to get strain resolution: §9 names network topology from a
    fixed .mgf as the safe artifact, and the clustering is the depositors' own,
    frozen and checksummed.
  * **Links.** Standardised Metcalf scoring over strain co-occurrence, the
    scheme NPLinker implements.

What this does NOT claim. The paper's own prioritisation is not the truth here,
and neither are PoDP's nine curated BGC-to-molecular-family links: those came
from manual structural work and MIBiG comparison, which PROJECT.md §9 excludes
from grading. The task measures execution of a prescribed linking pipeline.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matchms
import networkx as nx
import numpy as np
from matchms import Spectrum
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosineGreedy

MGF = "data/MSV000078836_specs_ms.mgf"
MATCHMS_VERSION = "0.33.1"

# Consensus spectrum filtering and networking, identical to task 01.
PRECURSOR_WINDOW_DA = 17.0
MIN_PEAKS = 6
FRAGMENT_TOLERANCE = 0.5
MIN_COSINE = 0.70
MIN_MATCHED_PEAKS = 6
MAX_PRECURSOR_DELTA = 2.0
TOP_K = 10

# Gene cluster families, identical to task 02.
JACCARD_THRESHOLD = 0.50
MIN_DOMAINS = 3

# Presence of a consensus spectrum in a strain.
PRESENCE_WINDOW_DA = 0.5
PRESENCE_MIN_COSINE = 0.70

# Standardised Metcalf weights (Doroghazi et al. 2014, as NPLinker implements
# them). The asymmetry is the point: a gene cluster family present in a strain
# whose metabolome lacks the molecular family is evidence against the link,
# while the reverse is not -- a molecule can be made by a cluster the detector
# missed, so that cell scores zero rather than penalising.
METCALF_BOTH = 10
METCALF_GCF_ONLY = -10
METCALF_FAMILY_ONLY = 0
METCALF_NEITHER = 1

TOP_LINKS = 500
DECIMALS = 4


def check_version() -> None:
    if matchms.__version__ != MATCHMS_VERSION:
        raise SystemExit(
            f"this task's truth was produced by matchms {MATCHMS_VERSION}, found "
            f"{matchms.__version__}. The similarity class was renamed between "
            f"releases and cosine scores are implementation-dependent."
        )


def load_strains(path: Path) -> list[dict]:
    with path.open() as handle:
        return sorted(csv.DictReader(handle, delimiter="\t"),
                      key=lambda r: r["accession"])


def species_of(organism: str) -> str:
    return " ".join(organism.split()[:2])


# ------------------------------------------------------------------- the genomes


def domain_set(field: str) -> frozenset[str]:
    """Distinct Pfam accessions in a GECCO `domains` field, which is a multiset."""
    return frozenset(d for d in field.split(";") if d)


def read_clusters(work: Path, accession: str) -> list[dict]:
    path = work / "gecco" / accession / f"{accession}.clusters.tsv"
    if not path.is_file():
        raise SystemExit(f"missing GECCO output for {accession}: {path}; run stage 1")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assign_families(clusters: dict[str, frozenset[str]]) -> dict[str, str]:
    """Single-linkage families over Jaccard similarity. Returns cluster -> family.

    A cluster with fewer than MIN_DOMAINS distinct domains is too thin to place
    and is left unassigned rather than forced into a family. A family is named
    for the lexicographically smallest cluster it contains, so no identifier
    depends on a tie-break nobody can reproduce.
    """
    eligible = sorted(c for c, d in clusters.items() if len(d) >= MIN_DOMAINS)
    parent = {c: c for c in eligible}

    def find(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for i, a in enumerate(eligible):
        for b in eligible[i + 1:]:
            da, db = clusters[a], clusters[b]
            union_size = len(da | db)
            if union_size and len(da & db) / union_size >= JACCARD_THRESHOLD:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    return {c: find(c) for c in eligible}


# --------------------------------------------------------------- the metabolomes


def load_consensus(path: Path) -> list:
    """The deposit's clustered consensus spectra, filtered as in task 01."""
    kept = []
    for spectrum in load_from_mgf(str(path)):
        if spectrum is None or len(spectrum.peaks.mz) == 0:
            continue
        normalised = normalize_intensities(spectrum)
        if normalised is None:
            continue
        trimmed = remove_peaks_around_precursor_mz(
            normalised, mz_tolerance=PRECURSOR_WINDOW_DA)
        if trimmed is None or len(trimmed.peaks.mz) < MIN_PEAKS:
            continue
        kept.append(trimmed)
    kept.sort(key=lambda s: int(s.get("scans")))
    return kept


def strain_spectra(path: Path) -> tuple[list, int]:
    """One strain's filtered MS2 spectra, from stage 2's feature file."""
    data = np.load(path)
    precursors = data["precursors"]
    offsets = data["offsets"]
    peak_mz, peak_intensity = data["peak_mz"], data["peak_intensity"]
    spectra = []
    for index, precursor in enumerate(precursors):
        lo, hi = int(offsets[index]), int(offsets[index + 1])
        spectra.append(Spectrum(
            mz=peak_mz[lo:hi], intensities=peak_intensity[lo:hi],
            metadata={"precursor_mz": float(precursor)},
            metadata_harmonization=False,
        ))
    order = np.argsort(precursors, kind="stable")
    spectra = [spectra[i] for i in order]
    return spectra, int(data["n_ms2_scans"][0])


def presence(consensus: list, spectra: list) -> set[int]:
    """Which consensus spectra this strain contains.

    Only pairs inside PRESENCE_WINDOW_DA are compared, found by bisection over
    the strain's precursor-sorted scans; the first scan to reach the cosine
    threshold settles it, so the answer does not depend on how many others
    would also have matched.
    """
    if not spectra:
        return set()
    similarity = ModifiedCosineGreedy(tolerance=FRAGMENT_TOLERANCE)
    masses = np.array([float(s.get("precursor_mz")) for s in spectra])
    found = set()
    for spectrum in consensus:
        mass = float(spectrum.get("precursor_mz"))
        lo, hi = np.searchsorted(masses, [mass - PRESENCE_WINDOW_DA,
                                          mass + PRESENCE_WINDOW_DA])
        for index in range(int(lo), int(hi)):
            if float(similarity.pair(spectrum, spectra[index])["score"]) >= PRESENCE_MIN_COSINE:
                found.add(int(spectrum.get("scans")))
                break
    return found


def candidate_edges(spectra: list) -> list[tuple[int, int]]:
    """Molecular network over the consensus spectra, as in task 01."""
    similarity = ModifiedCosineGreedy(tolerance=FRAGMENT_TOLERANCE)
    scans = [int(s.get("scans")) for s in spectra]
    masses = [float(s.get("precursor_mz")) for s in spectra]
    edges = []
    for i in range(len(spectra)):
        for j in range(i + 1, len(spectra)):
            if abs(masses[j] - masses[i]) > MAX_PRECURSOR_DELTA:
                continue
            result = similarity.pair(spectra[i], spectra[j])
            score = float(result["score"])
            if score >= MIN_COSINE and int(result["matches"]) >= MIN_MATCHED_PEAKS:
                a, b = sorted((scans[i], scans[j]))
                edges.append((a, b, round(score, DECIMALS)))
    return edges


def apply_top_k(edges: list[tuple[int, int, float]]) -> list[tuple[int, int]]:
    neighbours: dict[int, list] = defaultdict(list)
    for a, b, score in edges:
        neighbours[a].append((-score, b, (a, b)))
        neighbours[b].append((-score, a, (a, b)))
    kept_by = {}
    for scan, group in neighbours.items():
        group.sort(key=lambda item: (item[0], item[1]))
        kept_by[scan] = {pair for _, _, pair in group[:TOP_K]}
    return [(a, b) for a, b, _ in edges
            if (a, b) in kept_by[a] and (a, b) in kept_by[b]]


# -------------------------------------------------------------------- the linking


def metcalf(gcf_strains: set[str], family_strains: set[str],
            all_strains: set[str]) -> int:
    both = len(gcf_strains & family_strains)
    gcf_only = len(gcf_strains - family_strains)
    family_only = len(family_strains - gcf_strains)
    neither = len(all_strains) - both - gcf_only - family_only
    return (both * METCALF_BOTH + gcf_only * METCALF_GCF_ONLY
            + family_only * METCALF_FAMILY_ONLY + neither * METCALF_NEITHER)


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
    parser.add_argument("--work", type=Path, default=None)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    check_version()

    work = args.work or args.task / "work"
    strains = load_strains(args.task / "strains.tsv")
    accessions = [s["accession"] for s in strains]
    all_strains = set(accessions)
    species = {s["accession"]: species_of(s["organism"]) for s in strains}

    # --- genomes -------------------------------------------------------------
    domains_of: dict[str, frozenset[str]] = {}
    cluster_rows, strain_of_cluster = [], {}
    for strain in strains:
        accession = strain["accession"]
        for record in read_clusters(work, accession):
            cluster_id = record["cluster_id"]
            domains = domain_set(record.get("domains", ""))
            domains_of[cluster_id] = domains
            strain_of_cluster[cluster_id] = accession
            cluster_rows.append({
                "cluster_id": cluster_id,
                "accession": accession,
                "product_type": record.get("type") or "Unknown",
                "n_domains": len(domains),
            })
    cluster_rows.sort(key=lambda r: r["cluster_id"])
    print(f"{len(cluster_rows)} clusters over {len(strains)} strains", flush=True)
    write_csv(args.out / "bgc_inventory.csv",
              ["cluster_id", "accession", "product_type", "n_domains"], cluster_rows)

    family_of = assign_families(domains_of)
    strains_of_gcf = defaultdict(set)
    clusters_of_gcf = defaultdict(list)
    for cluster_id, family_id in family_of.items():
        strains_of_gcf[family_id].add(strain_of_cluster[cluster_id])
        clusters_of_gcf[family_id].append(cluster_id)
    write_csv(args.out / "gene_cluster_families.csv",
              ["gcf_id", "n_clusters", "n_strains", "n_species", "representative_cluster"],
              [{"gcf_id": gcf,
                "n_clusters": len(clusters_of_gcf[gcf]),
                "n_strains": len(strains_of_gcf[gcf]),
                "n_species": len({species[a] for a in strains_of_gcf[gcf]}),
                "representative_cluster": min(clusters_of_gcf[gcf])}
               for gcf in sorted(strains_of_gcf)])

    # --- metabolomes ---------------------------------------------------------
    consensus = load_consensus(args.task / MGF)
    print(f"{len(consensus)} consensus spectra retained", flush=True)

    present_strains = defaultdict(set)
    strain_rows = []
    for index, strain in enumerate(strains, start=1):
        accession = strain["accession"]
        feature_file = work / "features" / f"{accession}.npz"
        if not feature_file.is_file():
            raise SystemExit(f"missing features for {accession}: {feature_file}; run stage 2")
        spectra, n_ms2 = strain_spectra(feature_file)
        found = presence(consensus, spectra)
        for scan in found:
            present_strains[scan].add(accession)
        strain_rows.append({
            "accession": accession,
            "species": species[accession],
            "strain": strain["strain"],
            "n_ms2_scans": n_ms2,
            "n_usable_ms2": len(spectra),
            "n_spectra_present": len(found),
        })
        print(f"[{index}/{len(strains)}] {accession}: {len(spectra)} usable MS2, "
              f"{len(found)} consensus spectra present", flush=True)
    write_csv(args.out / "strain_inventory.csv",
              ["accession", "species", "strain", "n_ms2_scans", "n_usable_ms2",
               "n_spectra_present"], strain_rows)

    write_csv(args.out / "spectrum_occurrence.csv",
              ["scan", "accession"],
              [{"scan": scan, "accession": accession}
               for scan in sorted(present_strains)
               for accession in sorted(present_strains[scan])])

    # --- molecular families --------------------------------------------------
    edges = apply_top_k(candidate_edges(consensus))
    graph = nx.Graph()
    graph.add_nodes_from(int(s.get("scans")) for s in consensus)
    graph.add_edges_from(edges)
    members_of_family = {}
    for component in nx.connected_components(graph):
        members = sorted(component)
        members_of_family[f"MF{members[0]}"] = members
    print(f"{len(edges)} edges, {len(members_of_family)} molecular families", flush=True)

    strains_of_family = {
        family: set().union(*(present_strains[scan] for scan in members)) if members else set()
        for family, members in members_of_family.items()
    }
    write_csv(args.out / "molecular_families.csv",
              ["mf_id", "n_spectra", "n_strains", "n_species"],
              [{"mf_id": family,
                "n_spectra": len(members_of_family[family]),
                "n_strains": len(strains_of_family[family]),
                "n_species": len({species[a] for a in strains_of_family[family]})}
               for family in sorted(members_of_family)])

    # --- links ---------------------------------------------------------------
    # Only families seen in at least one strain can be linked; a family nobody
    # detected would score identically against every cluster family and would
    # fill the table with ties that mean nothing.
    linkable = sorted(f for f in members_of_family if strains_of_family[f])
    scores = []
    for gcf in sorted(strains_of_gcf):
        for family in linkable:
            scores.append({
                "gcf_id": gcf,
                "mf_id": family,
                "metcalf_score": metcalf(strains_of_gcf[gcf],
                                         strains_of_family[family], all_strains),
                "n_shared_strains": len(strains_of_gcf[gcf] & strains_of_family[family]),
            })
    print(f"{len(scores)} candidate links", flush=True)
    write_csv(args.out / "metcalf_scores.csv",
              ["gcf_id", "mf_id", "metcalf_score", "n_shared_strains"], scores)

    # Ranked descending by score, ties broken on the identifiers so the cut is
    # reproducible; rank is dense, so tied links share a rank.
    ordered = sorted(scores, key=lambda r: (-r["metcalf_score"], r["gcf_id"], r["mf_id"]))
    top = []
    rank = 0
    previous = None
    for row in ordered[:TOP_LINKS]:
        if row["metcalf_score"] != previous:
            rank += 1
            previous = row["metcalf_score"]
        top.append({**row, "rank": rank})
    write_csv(args.out / "top_links.csv",
              ["rank", "gcf_id", "mf_id", "metcalf_score", "n_shared_strains"], top)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
