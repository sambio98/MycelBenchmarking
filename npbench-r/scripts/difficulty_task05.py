"""Measure task 05's difficulty by grading plausible wrong pipelines.

Same approach and caveat as the task 02, 03 and 04 scripts. Task 05 was built
with shared upstream dependencies, like task 04; this says whether that
transfers to a second task or was a one-off.

The k-mer sketch is the expensive step, so sketches are cached per distinct
sketch configuration rather than per variant.

    python scripts/difficulty_task05.py --task tasks/05_salinispora_concordance \
        --truth tasks/05_salinispora_concordance/truth --scratch <dir>
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from npbench.grading import load_grading_config, score_task  # noqa: E402

MISTAKES = {
    "no_canonical_kmers": "uses forward k-mers only, ignoring contig orientation",
    "modulus_1024": "subsamples on % 1024 instead of % 1000",
    "kmer_multiset": "keeps duplicate k-mers instead of a set",
    "keeps_thin_clusters": "skips the fewer-than-3-distinct-domains exclusion",
    "no_transitive_closure": "treats each linked pair as a family",
    "families_numbered_by_rank": "names families GCF_0001.. by size rank",
    "species_keeps_strain": "groups by the whole organism string, strain included",
    "pairs_in_both_orders": "emits (a,b) and (b,a) for every pair",
    "equal_width_bins": "bins genomic distance into equal-width intervals, not deciles",
    "chemistry_over_clusters": "computes chemistry distance over clusters, not families",
}

_SKETCH_CACHE: dict[tuple, dict[str, np.ndarray]] = {}
_CLUSTERS: dict[str, list[dict]] = {}


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref05", task_dir / "solution" / "make_truth.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sketch(ref, archive: Path, mistakes: set[str]) -> np.ndarray:
    """The reference sketch, or the variant's wrong version of it."""
    modulus = 1024 if "modulus_1024" in mistakes else ref.MODULUS
    canonical = "no_canonical_kmers" not in mistakes
    dedupe = "kmer_multiset" not in mistakes
    kept = []
    K, MASK, BASE = ref.K, ref._MASK, ref._BASE
    with gzip.open(archive, "rb") as handle:
        for block in handle.read().split(b">")[1:]:
            sequence = block.split(b"\n", 1)[1].replace(b"\n", b"").upper()
            bases = BASE[np.frombuffer(sequence, dtype=np.uint8)]
            if len(bases) < K:
                continue
            usable = bases != 255
            n = len(bases) - K + 1
            forward = np.zeros(n, dtype=np.uint64)
            reverse = np.zeros(n, dtype=np.uint64)
            valid = np.ones(n, dtype=bool)
            for offset in range(K):
                base = bases[offset:n + offset].astype(np.uint64)
                forward = ((forward << np.uint64(2)) | base) & MASK
                reverse = (reverse >> np.uint64(2)) | (
                    (np.uint64(3) - base) << np.uint64(2 * (K - 1))
                )
                valid &= usable[offset:n + offset]
            chosen = np.minimum(forward, reverse) if canonical else forward
            kept.append(chosen[valid & ((chosen % np.uint64(modulus)) == 0)])
    if not kept:
        return np.zeros(0, dtype=np.uint64)
    joined = np.concatenate(kept)
    return np.unique(joined) if dedupe else np.sort(joined)


def sketches_for(ref, task_dir: Path, accessions: list[str], mistakes: set[str]):
    key = (
        "no_canonical_kmers" in mistakes,
        "modulus_1024" in mistakes,
        "kmer_multiset" in mistakes,
    )
    if key not in _SKETCH_CACHE:
        print(f"    (sketching for config {key})", flush=True)
        _SKETCH_CACHE[key] = {
            a: sketch(ref, task_dir / "data" / "genomes" / f"{a}.fna.gz", mistakes)
            for a in accessions
        }
    return _SKETCH_CACHE[key]


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    ua, ub = np.unique(a), np.unique(b)
    shared = np.intersect1d(ua, ub, assume_unique=True).size
    union = ua.size + ub.size - shared
    return 1.0 if union == 0 else 1.0 - shared / union


def build(task_dir: Path, work: Path, out: Path, mistakes: set[str]) -> None:
    ref = load_solution(task_dir)
    with (task_dir / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted(s["accession"] for s in strains)
    species = {
        s["accession"]: (s["organism"] if "species_keeps_strain" in mistakes
                         else " ".join(s["organism"].split()[:2]))
        for s in strains
    }

    if not _CLUSTERS:
        for accession in accessions:
            table = work / "gecco" / accession / f"{accession}.clusters.tsv"
            with table.open() as handle:
                _CLUSTERS[accession] = list(csv.DictReader(handle, delimiter="\t"))

    sk = sketches_for(ref, task_dir, accessions, mistakes)
    profiles = []
    for accession in accessions:
        contigs, total = ref.contig_stats(
            task_dir / "data" / "genomes" / f"{accession}.fna.gz"
        )
        profiles.append({"accession": accession, "n_contigs": contigs,
                         "assembly_bp": total, "n_kept_kmers": int(sk[accession].size)})
    w = ref.write_csv
    w(out / "genome_kmer_profiles.csv",
      ["accession", "n_contigs", "assembly_bp", "n_kept_kmers"], profiles)

    domains: dict[str, frozenset[str]] = {}
    accession_of: dict[str, str] = {}
    for accession in accessions:
        for row in _CLUSTERS[accession]:
            domains[row["cluster_id"]] = ref.domain_set(row["domains"])
            accession_of[row["cluster_id"]] = accession

    minimum = 0 if "keeps_thin_clusters" in mistakes else ref.MIN_DOMAINS
    eligible = sorted(c for c, d in domains.items() if len(d) >= minimum)
    index: dict[str, list[str]] = defaultdict(list)
    for cluster in eligible:
        for domain in domains[cluster]:
            index[domain].append(cluster)
    links = []
    seen: set[tuple[str, str]] = set()
    for cluster in eligible:
        for other in {o for d in domains[cluster] for o in index[d]}:
            if other == cluster:
                continue
            pair = (cluster, other) if cluster < other else (other, cluster)
            if pair in seen:
                continue
            seen.add(pair)
            a, b = domains[pair[0]], domains[pair[1]]
            if (a | b) and len(a & b) / len(a | b) >= ref.JACCARD_THRESHOLD:
                links.append(pair)

    if "no_transitive_closure" in mistakes:
        families = {c: c for c in eligible}
        for a, b in links:
            families[b] = min(families[b], a)
    else:
        parent = {c: c for c in eligible}

        def find(c: str) -> str:
            while parent[c] != c:
                parent[c] = parent[parent[c]]
                c = parent[c]
            return c

        for a, b in links:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
        families = {c: find(c) for c in eligible}

    if "families_numbered_by_rank" in mistakes:
        counts = Counter(families.values())
        order = sorted(counts, key=lambda f: (-counts[f], f))
        rename = {f: f"GCF_{i + 1:04d}" for i, f in enumerate(order)}
        families = {c: rename[f] for c, f in families.items()}

    size = Counter(families.values())
    w(out / "bgc_families.csv", ["cluster_id", "accession", "family_id", "family_size"],
      [{"cluster_id": c, "accession": accession_of[c], "family_id": families[c],
        "family_size": size[families[c]]} for c in sorted(families)])

    repertoire: dict[str, set[str]] = {a: set() for a in accessions}
    for cluster, family in families.items():
        repertoire[accession_of[cluster]].add(family)
    cluster_sets: dict[str, set[str]] = {a: set() for a in accessions}
    for cluster in domains:
        cluster_sets[accession_of[cluster]].add(cluster)
    clusters_per_genome = Counter(accession_of[c] for c in domains)

    genomes_of: dict[str, set[str]] = defaultdict(set)
    for accession, found in repertoire.items():
        for family in found:
            genomes_of[family].add(accession)
    w(out / "genome_repertoires.csv",
      ["accession", "n_clusters", "n_families", "n_unique_families"],
      [{"accession": a, "n_clusters": clusters_per_genome.get(a, 0),
        "n_families": len(repertoire[a]),
        "n_unique_families": sum(1 for f in repertoire[a] if len(genomes_of[f]) == 1)}
       for a in accessions])

    clusters_of: dict[str, list[str]] = defaultdict(list)
    for cluster, family in families.items():
        clusters_of[family].append(cluster)
    w(out / "family_taxonomic_breadth.csv",
      ["family_id", "n_clusters", "n_genomes", "n_species", "breadth",
       "representative_cluster"],
      [{"family_id": f, "n_clusters": len(clusters_of[f]), "n_genomes": len(genomes_of[f]),
        "n_species": len({species[a] for a in genomes_of[f]}),
        "breadth": "single_species" if len({species[a] for a in genomes_of[f]}) == 1
        else "multi_species",
        "representative_cluster": min(clusters_of[f])}
       for f in sorted(clusters_of)])

    chem_source = cluster_sets if "chemistry_over_clusters" in mistakes else repertoire
    pairs = []
    for i, a in enumerate(accessions):
        for b in accessions[i + 1:]:
            row = {
                "accession_a": a, "accession_b": b,
                "genomic_distance": round(jaccard(sk[a], sk[b]), ref.DECIMALS),
                "chemistry_distance": round(
                    1.0 - len(chem_source[a] & chem_source[b]) /
                    len(chem_source[a] | chem_source[b])
                    if (chem_source[a] | chem_source[b]) else 1.0, ref.DECIMALS),
                "same_species": "true" if species[a] == species[b] else "false",
            }
            pairs.append(row)
            if "pairs_in_both_orders" in mistakes:
                pairs.append({**row, "accession_a": b, "accession_b": a})
    w(out / "genome_pair_distances.csv",
      ["accession_a", "accession_b", "genomic_distance", "chemistry_distance",
       "same_species"], pairs)

    unique_pairs = [p for p in pairs if p["accession_a"] < p["accession_b"]]
    by_species: dict[str, list[dict]] = defaultdict(list)
    for pair in unique_pairs:
        if pair["same_species"] == "true":
            by_species[species[pair["accession_a"]]].append(pair)
    counts = Counter(species.values())
    w(out / "species_concordance.csv",
      ["species", "n_genomes", "n_pairs", "mean_genomic_distance",
       "mean_chemistry_distance"],
      [{"species": s, "n_genomes": counts[s], "n_pairs": len(g),
        "mean_genomic_distance": round(
            float(np.mean([p["genomic_distance"] for p in g])), ref.DECIMALS),
        "mean_chemistry_distance": round(
            float(np.mean([p["chemistry_distance"] for p in g])), ref.DECIMALS)}
       for s, g in sorted(by_species.items())])

    ordered = sorted(unique_pairs, key=lambda p: (p["genomic_distance"], p["accession_a"],
                                                  p["accession_b"]))
    groups = []
    if "equal_width_bins" in mistakes:
        values = np.array([p["genomic_distance"] for p in ordered])
        edges = np.linspace(values.min(), values.max(), ref.N_BINS + 1)
        for index in range(ref.N_BINS):
            lo, hi = edges[index], edges[index + 1]
            groups.append([p for p in ordered
                           if lo <= p["genomic_distance"] <= hi] if index == ref.N_BINS - 1
                          else [p for p in ordered if lo <= p["genomic_distance"] < hi])
    else:
        per_bin = len(ordered) // ref.N_BINS
        for index in range(ref.N_BINS):
            start = index * per_bin
            stop = len(ordered) if index == ref.N_BINS - 1 else start + per_bin
            groups.append(ordered[start:stop])
    w(out / "distance_concordance_bins.csv",
      ["bin", "mean_genomic_distance", "mean_chemistry_distance", "n_same_species_pairs"],
      [{"bin": i + 1,
        "mean_genomic_distance": round(
            float(np.mean([p["genomic_distance"] for p in g])), ref.DECIMALS) if g else 0.0,
        "mean_chemistry_distance": round(
            float(np.mean([p["chemistry_distance"] for p in g])), ref.DECIMALS) if g else 0.0,
        "n_same_species_pairs": sum(1 for p in g if p["same_species"] == "true")}
       for i, g in enumerate(groups)])


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: forward k-mers only": {"no_canonical_kmers"},
    "one slip: modulus 1024": {"modulus_1024"},
    "one slip: k-mer multiset": {"kmer_multiset"},
    "one slip: thin clusters kept": {"keeps_thin_clusters"},
    "one slip: no transitive closure": {"no_transitive_closure"},
    "one slip: families numbered by rank": {"families_numbered_by_rank"},
    "one slip: species keeps the strain": {"species_keeps_strain"},
    "one slip: pairs in both orders": {"pairs_in_both_orders"},
    "one slip: equal-width bins": {"equal_width_bins"},
    "one slip: chemistry over clusters": {"chemistry_over_clusters"},
    "hasty agent (3 slips)": {"no_canonical_kmers", "species_keeps_strain",
                              "equal_width_bins"},
    "careless agent (5 slips)": {"modulus_1024", "families_numbered_by_rank",
                                 "species_keeps_strain", "chemistry_over_clusters",
                                 "pairs_in_both_orders"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=None)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    work = args.work or args.task / "work"
    config = load_grading_config(args.task / "grading.yaml")

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':40}{'task':>6}  " + "".join(f"{n[:12]:>14}" for n in names)
    print(header)
    print("-" * len(header))
    passes = Counter()
    for label, mistakes in VARIANTS.items():
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace(",", "")
        build(args.task, work, out, mistakes)
        result = score_task(config, args.truth, out)
        row = f"{label:40}{result.task_score:6.2f}  "
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
