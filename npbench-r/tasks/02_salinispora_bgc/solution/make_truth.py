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
    for accession in accessions:
        lengths = contig_lengths(args.task / "data" / "genomes" / f"{accession}.fna.gz")
        for row in read_clusters(work, accession):
            start, end = int(row["start"]), int(row["end"])
            sequence_id = row["sequence_id"]
            inventory.append(
                {
                    "cluster_id": row["cluster_id"],
                    "accession": accession,
                    "sequence_id": sequence_id,
                    "start": start,
                    "end": end,
                    "length_bp": end - start + 1,
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
    per_type = Counter(r["type"] for r in inventory)
    genomes_per_type: dict[str, set[str]] = {}
    for row in inventory:
        genomes_per_type.setdefault(row["type"], set()).add(row["accession"])
    per_species_type = Counter((species[r["accession"]], r["type"]) for r in inventory)

    print(f"{len(inventory)} clusters across {len(accessions)} genomes")
    write_csv(args.out / "bgc_inventory.csv",
              ["cluster_id", "accession", "sequence_id", "start", "end", "length_bp", "type"],
              inventory)
    write_csv(args.out / "contig_edge_clusters.csv",
              ["cluster_id", "accession", "sequence_id", "contig_length_bp", "boundary"],
              edge)
    # Every staged assembly gets a row, including any with no clusters at all:
    # dropping them would hide a genome the agent failed to analyse.
    write_csv(args.out / "clusters_per_genome.csv", ["accession", "n_clusters"],
              [{"accession": a, "n_clusters": per_genome.get(a, 0)} for a in accessions])
    write_csv(args.out / "type_counts.csv", ["type", "n_clusters", "n_genomes"],
              [{"type": t, "n_clusters": per_type[t], "n_genomes": len(genomes_per_type[t])}
               for t in sorted(per_type)])
    write_csv(args.out / "species_type_counts.csv", ["species", "type", "n_clusters"],
              [{"species": s, "type": t, "n_clusters": n}
               for (s, t), n in sorted(per_species_type.items())])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
