"""Reference solution for tasks/05_salinispora_concordance.

Asks whether *Salinispora* strains that are genomically close carry similar
biosynthetic repertoires -- taxonomy versus chemistry concordance, PROJECT.md
§4 category 8.

    python solution/make_truth.py --task tasks/05_salinispora_concordance --out <dir>

Two independent distance matrices are built over the same 120 assemblies and
then compared:

  * a **genomic** distance from shared canonical 21-mers, and
  * a **chemistry** distance from shared BGC families.

The artifacts share upstream dependencies on purpose. The k-mer parse feeds the
genomic matrix; the GECCO families feed the repertoires, the family breadth and
the chemistry matrix; and both matrices feed the pair table, which feeds the two
summaries. Task 04's notes record why that matters.

On the k-mer sketch: a MinHash would need its hash function pinned exactly, and
an agent that chose a different one would fail every distance for a reason that
is not science. Instead the sketch keeps every canonical 21-mer whose 2-bit
encoding is divisible by 1000. There is no hash function to agree on, only the
encoding and the modulus, and it separates the species slightly better than a
2,000-element MinHash did on the same data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

GECCO_VERSION = "0.11.0"

# k-mer sketch. A=0, C=1, G=2, T=3; a k-mer is the smaller of its forward and
# reverse-complement encodings, because contig orientation is arbitrary.
K = 21
MODULUS = 1000
_MASK = np.uint64((1 << (2 * K)) - 1)
_BASE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _BASE[_b] = _i

# BGC family rule, identical to task 02 so the two tasks agree on what a family
# is: distinct-Pfam-domain sets, Jaccard >= 0.50, single linkage, named for the
# lexicographically smallest member cluster.
JACCARD_THRESHOLD = 0.50
MIN_DOMAINS = 3

N_BINS = 10
DECIMALS = 6


def canonical_kmers(archive: Path) -> np.ndarray:
    """The kept canonical 21-mers of one gzipped assembly, sorted and unique."""
    kept = []
    with gzip.open(archive, "rb") as handle:
        for block in handle.read().split(b">")[1:]:
            sequence = block.split(b"\n", 1)[1].replace(b"\n", b"").upper()
            bases = _BASE[np.frombuffer(sequence, dtype=np.uint8)]
            if len(bases) < K:
                continue
            usable = bases != 255
            n = len(bases) - K + 1
            forward = np.zeros(n, dtype=np.uint64)
            reverse = np.zeros(n, dtype=np.uint64)
            valid = np.ones(n, dtype=bool)
            for offset in range(K):
                base = bases[offset:n + offset].astype(np.uint64)
                forward = ((forward << np.uint64(2)) | base) & _MASK
                reverse = (reverse >> np.uint64(2)) | (
                    (np.uint64(3) - base) << np.uint64(2 * (K - 1))
                )
                valid &= usable[offset:n + offset]
            canonical = np.minimum(forward, reverse)
            kept.append(canonical[valid & ((canonical % np.uint64(MODULUS)) == 0)])
    if not kept:
        return np.zeros(0, dtype=np.uint64)
    return np.unique(np.concatenate(kept))


def contig_stats(archive: Path) -> tuple[int, int]:
    """(number of contigs, total bases) for one gzipped assembly."""
    contigs = 0
    total = 0
    with gzip.open(archive, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                contigs += 1
            else:
                total += len(line.strip())
    return contigs, total


def domain_set(field: str) -> frozenset[str]:
    return frozenset(d for d in field.split(";") if d)


def assign_families(clusters: dict[str, frozenset[str]]) -> dict[str, str]:
    """Single-linkage families over Jaccard similarity of domain sets."""
    eligible = sorted(c for c, d in clusters.items() if len(d) >= MIN_DOMAINS)
    parent = {c: c for c in eligible}

    def find(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    index: dict[str, list[str]] = defaultdict(list)
    for cluster in eligible:
        for domain in clusters[cluster]:
            index[domain].append(cluster)

    checked: set[tuple[str, str]] = set()
    for cluster in eligible:
        for other in {o for d in clusters[cluster] for o in index[d]}:
            if other == cluster:
                continue
            pair = (cluster, other) if cluster < other else (other, cluster)
            if pair in checked:
                continue
            checked.add(pair)
            a, b = clusters[pair[0]], clusters[pair[1]]
            if len(a & b) / len(a | b) >= JACCARD_THRESHOLD:
                ra, rb = find(pair[0]), find(pair[1])
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    return {c: find(c) for c in eligible}


def jaccard_distance(a: np.ndarray, b: np.ndarray) -> float:
    shared = np.intersect1d(a, b, assume_unique=True).size
    union = a.size + b.size - shared
    return 1.0 if union == 0 else 1.0 - shared / union


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation of average ranks. Avoids a scipy dependency."""
    def ranks(v: np.ndarray) -> np.ndarray:
        order = np.argsort(v, kind="stable")
        out = np.empty(len(v), dtype=float)
        sorted_v = v[order]
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and sorted_v[j + 1] == sorted_v[i]:
                j += 1
            out[order[i:j + 1]] = (i + j) / 2 + 1
            i = j + 1
        return out

    rx, ry = ranks(x), ranks(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=None, help="default: <task>/work")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    work = args.work or args.task / "work"

    with (args.task / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted(s["accession"] for s in strains)
    species = {s["accession"]: " ".join(s["organism"].split()[:2]) for s in strains}

    # --- genomic side -----------------------------------------------------
    sketches: dict[str, np.ndarray] = {}
    profiles = []
    for accession in accessions:
        archive = args.task / "data" / "genomes" / f"{accession}.fna.gz"
        sketches[accession] = canonical_kmers(archive)
        contigs, total = contig_stats(archive)
        profiles.append({
            "accession": accession,
            "n_contigs": contigs,
            "assembly_bp": total,
            "n_kept_kmers": int(sketches[accession].size),
        })
    print(f"sketched {len(accessions)} assemblies", flush=True)
    write_csv(args.out / "genome_kmer_profiles.csv",
              ["accession", "n_contigs", "assembly_bp", "n_kept_kmers"], profiles)

    # --- chemistry side ---------------------------------------------------
    domains: dict[str, frozenset[str]] = {}
    accession_of: dict[str, str] = {}
    for accession in accessions:
        table = work / "gecco" / accession / f"{accession}.clusters.tsv"
        if not table.is_file():
            raise SystemExit(f"missing GECCO output for {accession}: {table}")
        with table.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                domains[row["cluster_id"]] = domain_set(row["domains"])
                accession_of[row["cluster_id"]] = accession

    families = assign_families(domains)
    size = Counter(families.values())
    print(f"{len(size)} families over {len(families)} assignable clusters", flush=True)
    write_csv(args.out / "bgc_families.csv",
              ["cluster_id", "accession", "family_id", "family_size"],
              [{"cluster_id": c, "accession": accession_of[c],
                "family_id": families[c], "family_size": size[families[c]]}
               for c in sorted(families)])

    repertoire: dict[str, set[str]] = {a: set() for a in accessions}
    clusters_per_genome = Counter()
    for cluster, family in families.items():
        repertoire[accession_of[cluster]].add(family)
    for cluster in domains:
        clusters_per_genome[accession_of[cluster]] += 1

    genomes_of_family: dict[str, set[str]] = defaultdict(set)
    for accession, found in repertoire.items():
        for family in found:
            genomes_of_family[family].add(accession)
    family_members = Counter(families.values())

    def largest_family(found: set[str]) -> str:
        """The genome's most populous family, ties to the smallest family_id."""
        if not found:
            return ""
        return min(found, key=lambda f: (-family_members[f], f))

    write_csv(args.out / "genome_repertoires.csv",
              ["accession", "n_clusters", "n_families", "n_unique_families",
               "largest_family_id"],
              [{"accession": a, "n_clusters": clusters_per_genome.get(a, 0),
                "n_families": len(repertoire[a]),
                "n_unique_families": sum(
                    1 for f in repertoire[a] if len(genomes_of_family[f]) == 1),
                "largest_family_id": largest_family(repertoire[a])}
               for a in accessions])

    clusters_of_family: dict[str, list[str]] = defaultdict(list)
    for cluster, family in families.items():
        clusters_of_family[family].append(cluster)
    write_csv(args.out / "family_taxonomic_breadth.csv",
              ["family_id", "n_clusters", "n_genomes", "n_species", "breadth",
               "representative_cluster"],
              [{"family_id": f,
                "n_clusters": len(clusters_of_family[f]),
                "n_genomes": len(genomes_of_family[f]),
                "n_species": len({species[a] for a in genomes_of_family[f]}),
                "breadth": "single_species"
                if len({species[a] for a in genomes_of_family[f]}) == 1 else "multi_species",
                "representative_cluster": min(clusters_of_family[f])}
               for f in sorted(clusters_of_family)])

    # --- the two matrices, compared --------------------------------------
    pairs = []
    for i, a in enumerate(accessions):
        for b in accessions[i + 1:]:
            pairs.append({
                "accession_a": a,
                "accession_b": b,
                "genomic_distance": round(jaccard_distance(sketches[a], sketches[b]), DECIMALS),
                "chemistry_distance": round(
                    1.0 - (len(repertoire[a] & repertoire[b]) / len(repertoire[a] | repertoire[b]))
                    if (repertoire[a] | repertoire[b]) else 1.0, DECIMALS),
                "same_species": "true" if species[a] == species[b] else "false",
            })
    print(f"{len(pairs)} genome pairs", flush=True)
    write_csv(args.out / "genome_pair_distances.csv",
              ["accession_a", "accession_b", "genomic_distance", "chemistry_distance",
               "same_species"], pairs)

    by_species: dict[str, list[dict]] = defaultdict(list)
    for pair in pairs:
        if pair["same_species"] == "true":
            by_species[species[pair["accession_a"]]].append(pair)
    counts = Counter(species.values())
    write_csv(args.out / "species_concordance.csv",
              ["species", "n_genomes", "n_pairs", "mean_genomic_distance",
               "mean_chemistry_distance"],
              [{"species": s, "n_genomes": counts[s], "n_pairs": len(group),
                "mean_genomic_distance": round(
                    float(np.mean([p["genomic_distance"] for p in group])), DECIMALS),
                "mean_chemistry_distance": round(
                    float(np.mean([p["chemistry_distance"] for p in group])), DECIMALS)}
               for s, group in sorted(by_species.items())])

    # Deciles of genomic distance, each holding an equal number of pairs. Ties
    # are broken by (accession_a, accession_b) so the split is deterministic.
    ordered = sorted(pairs, key=lambda p: (p["genomic_distance"], p["accession_a"],
                                           p["accession_b"]))
    per_bin = len(ordered) // N_BINS
    bins = []
    for index in range(N_BINS):
        start = index * per_bin
        stop = len(ordered) if index == N_BINS - 1 else start + per_bin
        group = ordered[start:stop]
        bins.append({
            "bin": index + 1,
            "n_pairs": len(group),
            "mean_genomic_distance": round(
                float(np.mean([p["genomic_distance"] for p in group])), DECIMALS),
            "mean_chemistry_distance": round(
                float(np.mean([p["chemistry_distance"] for p in group])), DECIMALS),
            "n_same_species_pairs": sum(1 for p in group if p["same_species"] == "true"),
        })
    write_csv(args.out / "distance_concordance_bins.csv",
              ["bin", "n_pairs", "mean_genomic_distance", "mean_chemistry_distance",
               "n_same_species_pairs"], bins)

    g = np.array([p["genomic_distance"] for p in pairs])
    c = np.array([p["chemistry_distance"] for p in pairs])
    print(f"concordance over all pairs: pearson={np.corrcoef(g, c)[0, 1]:.4f} "
          f"spearman={spearman(g, c):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
