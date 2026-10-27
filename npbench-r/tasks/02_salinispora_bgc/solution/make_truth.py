"""Stage 2 of the reference solution: turn GECCO output into the graded artifacts.

Reads <work>/gecco/<accession>/<accession>.clusters.tsv for every assembly,
joins against the staged genomes (for contig lengths) and data/strains.tsv (for
species), and writes the five artifacts.

    python solution/make_truth.py --task tasks/02_salinispora_bgc --out <dir>

Writing to truth/ generates the ground truth; writing anywhere else produces a
submission, which is how `make selfgrade` checks the grader against the truth
it was built from.

Every output is sorted on its row key, so two runs produce byte-identical files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

# A cluster counts as sitting on a contig boundary if either end falls within
# this many bases of it. 119 of the 120 assemblies are scaffold-level, so a
# cluster that runs off the end of a scaffold is truncated rather than absent,
# and the inventory over-reports complete clusters if this is ignored.
EDGE_DISTANCE = 1000

# Gene cluster family assignment. PROJECT.md section 4, category 2.
#
# Families are connected components (single linkage) over pairs of clusters
# whose distinct-Pfam-domain sets reach JACCARD_THRESHOLD. Clusters carrying
# fewer than MIN_DOMAINS distinct domains are too thin to place and are left
# unassigned rather than forced into a family.
#
# A family is named for the lexicographically smallest cluster_id it contains.
# That is deliberate: naming a family by its rank would make every downstream
# row depend on a tie-break nobody can reproduce, so the identifier is derived
# from the family's own contents instead.
JACCARD_THRESHOLD = 0.50
MIN_DOMAINS = 3


def domain_set(field: str) -> frozenset[str]:
    """Distinct Pfam accessions in a GECCO `domains` field.

    The field is a multiset: a cluster with three copies of PF00165 lists it
    three times. Jaccard is defined over sets, so the duplicates come out here.
    """
    return frozenset(d for d in field.split(";") if d)


def assign_families(clusters: dict[str, frozenset[str]]) -> dict[str, str]:
    """Single-linkage families over Jaccard similarity. Returns cluster -> family."""
    eligible = sorted(c for c, d in clusters.items() if len(d) >= MIN_DOMAINS)

    parent = {c: c for c in eligible}

    def find(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            # Keep the smaller id as the root so the component name is stable.
            parent[max(ra, rb)] = min(ra, rb)

    # An inverted index over domains: two clusters can only reach the threshold
    # if they share at least one domain, which skips most of the n^2 pairs.
    by_domain: dict[str, list[str]] = {}
    for cluster in eligible:
        for domain in clusters[cluster]:
            by_domain.setdefault(domain, []).append(cluster)

    checked: set[tuple[str, str]] = set()
    for cluster in eligible:
        candidates = {other for d in clusters[cluster] for other in by_domain[d]}
        for other in candidates:
            if other == cluster:
                continue
            pair = (cluster, other) if cluster < other else (other, cluster)
            if pair in checked:
                continue
            checked.add(pair)
            a, b = clusters[pair[0]], clusters[pair[1]]
            if len(a & b) / len(a | b) >= JACCARD_THRESHOLD:
                union(*pair)

    return {c: find(c) for c in eligible}



def contig_lengths(archive: Path) -> dict[str, int]:
    """Length of every sequence in a gzipped FASTA, without loading it all."""
    lengths: dict[str, int] = {}
    name = None
    total = 0
    with gzip.open(archive, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths


def species_of(organism: str) -> str:
    """First two tokens of the organism string.

    NCBI records this inconsistently: some entries carry the strain
    ("Salinispora arenicola CNT849") and some do not. Both must collapse to one
    species, and "Salinispora sp." stays its own group.
    """
    return " ".join(organism.split()[:2])


def read_clusters(work: Path, accession: str) -> list[dict]:
    path = work / "gecco" / accession / f"{accession}.clusters.tsv"
    if not path.is_file():
        raise SystemExit(f"missing GECCO output for {accession}: {path}")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _median(values: list[int]) -> int:
    """Lower-of-the-two-middles for an even count, so the result stays an integer."""
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) // 2


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
    parser.add_argument("--work", type=Path, default=None, help="default: <task>/work")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    work = args.work or args.task / "work"

    with (args.task / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted(s["accession"] for s in strains)
    species = {s["accession"]: species_of(s["organism"]) for s in strains}

    inventory: list[dict] = []
    edge: list[dict] = []
    domains: dict[str, frozenset[str]] = {}
    assembly_bp: dict[str, int] = {}
    for accession in accessions:
        lengths = contig_lengths(args.task / "data" / "genomes" / f"{accession}.fna.gz")
        assembly_bp[accession] = sum(lengths.values())
        for row in read_clusters(work, accession):
            start, end = int(row["start"]), int(row["end"])
            sequence_id = row["sequence_id"]
            domains[row["cluster_id"]] = domain_set(row["domains"])
            inventory.append(
                {
                    "cluster_id": row["cluster_id"],
                    "accession": accession,
                    "sequence_id": sequence_id,
                    "start": start,
                    "end": end,
                    "length_bp": end - start + 1,
                    "n_distinct_domains": len(domain_set(row["domains"])),
                    "type": row["type"],
                }
            )
            contig_length = lengths[sequence_id]
            at_start = start <= EDGE_DISTANCE
            at_end = end >= contig_length - EDGE_DISTANCE + 1
            if at_start or at_end:
                edge.append(
                    {
                        "cluster_id": row["cluster_id"],
                        "accession": accession,
                        "sequence_id": sequence_id,
                        "contig_length_bp": contig_length,
                        "boundary": "both" if at_start and at_end
                        else ("start" if at_start else "end"),
                    }
                )

    inventory.sort(key=lambda r: r["cluster_id"])
    edge.sort(key=lambda r: r["cluster_id"])

    per_genome = Counter(r["accession"] for r in inventory)
    edge_per_genome = Counter(r["accession"] for r in edge)
    cluster_bp = Counter()
    for row in inventory:
        cluster_bp[row["accession"]] += row["length_bp"]
    per_type = Counter(r["type"] for r in inventory)
    genomes_per_type: dict[str, set[str]] = {}
    for row in inventory:
        genomes_per_type.setdefault(row["type"], set()).add(row["accession"])
    per_species_type = Counter((species[r["accession"]], r["type"]) for r in inventory)

    print(f"{len(inventory)} clusters across {len(accessions)} genomes")
    write_csv(args.out / "bgc_inventory.csv",
              ["cluster_id", "accession", "sequence_id", "start", "end", "length_bp",
               "n_distinct_domains", "type"],
              inventory)
    write_csv(args.out / "contig_edge_clusters.csv",
              ["cluster_id", "accession", "sequence_id", "contig_length_bp", "boundary"],
              edge)
    # Every staged assembly gets a row, including any with no clusters at all:
    # dropping them would hide a genome the agent failed to analyse.
    write_csv(
        args.out / "clusters_per_genome.csv",
        ["accession", "n_clusters", "n_edge_clusters", "assembly_bp", "pct_bp_in_clusters"],
        [
            {
                "accession": a,
                "n_clusters": per_genome.get(a, 0),
                "n_edge_clusters": edge_per_genome.get(a, 0),
                "assembly_bp": assembly_bp[a],
                "pct_bp_in_clusters": round(100 * cluster_bp.get(a, 0) / assembly_bp[a], 4),
            }
            for a in accessions
        ],
    )
    species_per_type: dict[str, set[str]] = {}
    for row in inventory:
        species_per_type.setdefault(row["type"], set()).add(species[row["accession"]])
    lengths_per_type: dict[str, list[int]] = {}
    for row in inventory:
        lengths_per_type.setdefault(row["type"], []).append(row["length_bp"])
    write_csv(
        args.out / "type_counts.csv",
        ["type", "n_clusters", "n_genomes", "n_species", "median_length_bp"],
        [
            {
                "type": t,
                "n_clusters": per_type[t],
                "n_genomes": len(genomes_per_type[t]),
                "n_species": len(species_per_type[t]),
                "median_length_bp": _median(lengths_per_type[t]),
            }
            for t in sorted(per_type)
        ],
    )
    write_csv(args.out / "species_type_counts.csv", ["species", "type", "n_clusters"],
              [{"species": s, "type": t, "n_clusters": n}
               for (s, t), n in sorted(per_species_type.items())])

    families = assign_families(domains)
    size = Counter(families.values())
    assignment = [
        {"cluster_id": c, "family_id": f, "family_size": size[f]}
        for c, f in sorted(families.items())
    ]
    accession_of = {r["cluster_id"]: r["accession"] for r in inventory}
    members: dict[str, list[str]] = {}
    for cluster, family in families.items():
        members.setdefault(family, []).append(cluster)
    summary = [
        {
            "family_id": family,
            "n_clusters": len(group),
            "n_genomes": len({accession_of[c] for c in group}),
            "n_species": len({species[accession_of[c]] for c in group}),
        }
        for family, group in sorted(members.items())
    ]
    print(f"{len(summary)} families over {len(assignment)} assignable clusters "
          f"({len(domains) - len(assignment)} below {MIN_DOMAINS} distinct domains)")
    write_csv(args.out / "gene_cluster_families.csv",
              ["cluster_id", "family_id", "family_size"], assignment)
    write_csv(args.out / "gcf_summary.csv",
              ["family_id", "n_clusters", "n_genomes", "n_species"], summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
