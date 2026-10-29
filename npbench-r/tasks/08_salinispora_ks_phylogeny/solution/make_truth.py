"""Reference solution for tasks/08_salinispora_ks_phylogeny.

Ketosynthase phylogeny across the sequenced Salinispora, PROJECT.md §4
category 8: "phylogenetics of biosynthetic enzymes".

    python solution/make_truth.py --task tasks/08_salinispora_ks_phylogeny --out <dir>

The ketosynthase domain is the standard marker for polyketide synthase
phylogeny: it is present in every PKS module, long enough to align, and its
tree tracks the pathway rather than the organism, which is what makes the
taxonomy comparison in split_taxonomy.csv worth computing.

Why UPGMA from a distance matrix rather than a maximum-likelihood tree. An ML
tree needs RAxML or IQ-TREE, neither of which is pip-installable here, and
both of which are stochastic unless the seed is pinned -- PROJECT.md §9 wants
artifacts that are reproducible by construction. Average linkage over a
normalised alignment-score distance is UPGMA, a classical phylogenetic method,
and it is exact and deterministic given the matrix.

The tree is graded by its splits, not by a Newick string. A Newick string
depends on child order and on branch-length formatting; the set of splits is
the tree's actual topological content, and it is what Robinson-Foulds
compares. Each split is named by the lexicographically smallest domain below
it, so no identifier depends on an ordering the agent chose.

Alignment scores, never the traceback: a score is stable across Biopython
releases in a way the chosen alignment path is not.

The tree is built over **distinct domain sequences**, not over all 2,993
domains, and that is not an optimisation. 60 of the 120 assemblies are
S. arenicola, so the same ketosynthase recurs verbatim across strains: the
2,993 domains carry only 922 distinct sequences, one of them 53 times. Putting
all 2,993 in the tree leaves 2,552 of them at distance exactly 0 from their
nearest neighbour, and a UPGMA tree over that many zero-height ties has a
topology decided by the input order rather than by the sequences. It would
still be reproducible, and it would still be arbitrary. Dereplicating first is
also what a phylogeneticist does, and it makes the all-pairs step ten times
cheaper.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import gecco
import numpy as np
import pyhmmer
import pyrodigal
from Bio import Align
from Bio.Align import substitution_matrices
from pyhmmer.easel import Alphabet, TextSequence
from scipy.cluster.hierarchy import linkage

GENOMES = "data/genomes"
STRAINS = "strains.tsv"

PYRODIGAL_VERSION = "3.7.1"
PYHMMER_VERSION = "0.12.3"

# Pfam 35.0, as bundled with GECCO. PF00109 is the N-terminal ketoacyl
# synthase domain, the PKS marker.
KS_ACCESSION = "PF00109"
EVALUE = 1e-5

# A PF00109 alignment shorter than this is a fragment, not a domain: the model
# is 253 positions long, so 150 keeps partial but real domains and drops the
# tail ends of truncated genes at contig boundaries.
MIN_DOMAIN_AA = 150

GAP_OPEN = -11
GAP_EXTEND = -1
DECIMALS = 4


def check_versions() -> None:
    for module, pinned in ((pyrodigal, PYRODIGAL_VERSION), (pyhmmer, PYHMMER_VERSION)):
        if module.__version__ != pinned:
            raise SystemExit(
                f"this task's truth was produced by {module.__name__} {pinned}, found "
                f"{module.__version__}. Gene calls and domain boundaries move between "
                f"releases, so truth is only valid against the pinned build."
            )


def _hmm_path() -> Path:
    """The KS model, found through the installed gecco package.

    GECCO ships a 2,766-model subset of Pfam 35.0 with a recorded md5, so the
    task needs no database download and fetches nothing at runtime. Found
    through the package rather than an absolute path because the install
    location differs between environments.
    """
    path = Path(gecco.__file__).parent / "hmmer" / "Pfam.h3m"
    if not path.is_file():
        raise SystemExit(f"GECCO's bundled Pfam not found at {path}")
    return path


def load_strains(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_contigs(path: Path) -> list[str]:
    contigs, current = [], []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    contigs.append("".join(current))
                current = []
            else:
                current.append(line.strip())
    if current:
        contigs.append("".join(current))
    return contigs


def call_proteins(contigs: list[str]) -> list[tuple[str, str]]:
    """Gene calls for one assembly, trained on that assembly.

    Training on the whole assembly rather than per contig, and not using
    meta mode: both alternatives change the gene set, and PROJECT.md §9 asks
    for the choice to be written into the method rather than left open.
    """
    finder = pyrodigal.GeneFinder()
    finder.train(*contigs)
    proteins = []
    for index, contig in enumerate(contigs, start=1):
        for number, gene in enumerate(finder.find_genes(contig), start=1):
            proteins.append((f"c{index:04d}_g{number:05d}", gene.translate()))
    return proteins


def find_domains(models, proteins: list[tuple[str, str]], alphabet) -> list[dict]:
    """Included PF00109 domains, with their coordinates on the protein.

    Only hits and domains HMMER itself includes are kept, and the coordinates
    are the aligned region on the target protein, 1-based inclusive. The
    envelope would be a defensible alternative and gives slightly longer
    sequences; the method states which one is meant.
    """
    sequence_of = dict(proteins)
    digital = [
        TextSequence(name=name.encode(), sequence=sequence).digitize(alphabet)
        for name, sequence in proteins
    ]
    found = []
    for top in pyhmmer.hmmsearch(models, digital, E=EVALUE, cpus=4):
        for hit in top:
            if not hit.included:
                continue
            # pyhmmer returns str in some releases and bytes in others.
            protein_id = (hit.name.decode() if isinstance(hit.name, bytes)
                          else str(hit.name))
            for domain in hit.domains:
                if not domain.included:
                    continue
                alignment = domain.alignment
                start, end = alignment.target_from, alignment.target_to
                if end - start + 1 < MIN_DOMAIN_AA:
                    continue
                found.append({
                    "protein_id": protein_id,
                    "start": start,
                    "end": end,
                    "bitscore": round(domain.score, 1),
                    "sequence": sequence_of[protein_id][start - 1:end],
                })
    # Sorted so the domain numbering below does not depend on HMMER's hit order.
    found.sort(key=lambda d: (d["protein_id"], d["start"]))
    return found


def distances(sequences: list[str]) -> np.ndarray:
    """Condensed pairwise distance matrix over normalised local alignment score.

    distance = 1 - score(a,b) / min(score(a,a), score(b,b)), clamped to [0,1].
    Normalising by the smaller self-score rather than by length is what keeps a
    short domain from looking distant from the long one it sits inside.
    """
    aligner = Align.PairwiseAligner(scoring=None)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND
    aligner.mode = "local"

    self_score = [aligner.score(s, s) for s in sequences]
    n = len(sequences)
    condensed = np.empty(n * (n - 1) // 2, dtype=float)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            floor = min(self_score[i], self_score[j])
            ratio = aligner.score(sequences[i], sequences[j]) / floor if floor else 0.0
            condensed[k] = min(1.0, max(0.0, 1.0 - ratio))
            k += 1
        if (i + 1) % 250 == 0 or i + 1 == n:
            print(f"    distances: row {i + 1}/{n}", flush=True)
    print(flush=True)
    return condensed


def sequence_groups(domain_rows: list[dict], sequence_of: dict[str, str],
                    species_of: dict[str, str]) -> tuple[list[dict], dict[str, list[str]]]:
    """Distinct domain sequences, each named by its smallest member domain.

    This is the tree's leaf set. Naming a group by its smallest member rather
    than by a running number keeps the identifier independent of any ordering
    the agent chose, the same rule the splits use.
    """
    members: dict[str, list[str]] = defaultdict(list)
    for row in domain_rows:
        members[sequence_of[row["domain_id"]]].append(row["domain_id"])

    accession_of = {r["domain_id"]: r["accession"] for r in domain_rows}
    groups, by_group = [], {}
    for sequence, domain_ids in members.items():
        domain_ids = sorted(domain_ids)
        group_id = domain_ids[0]
        by_group[group_id] = domain_ids
        genomes = {accession_of[d] for d in domain_ids}
        groups.append({
            "group_id": group_id,
            "length_aa": len(sequence),
            "n_domains": len(domain_ids),
            "n_genomes": len(genomes),
            "n_species": len({species_of[a] for a in genomes}),
        })
    groups.sort(key=lambda g: g["group_id"])
    return groups, by_group


def nearest_neighbours(ids: list[str], condensed: np.ndarray) -> list[dict]:
    """Each group's closest other group, ties broken on the partner's id."""
    n = len(ids)
    square = np.zeros((n, n), dtype=float)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            square[i, j] = square[j, i] = condensed[k]
            k += 1
    rows = []
    for i in range(n):
        best, best_id = None, None
        for j in range(n):
            if i == j:
                continue
            value = round(float(square[i, j]), DECIMALS)
            if best is None or (value, ids[j]) < (best, best_id):
                best, best_id = value, ids[j]
        rows.append({"group_id": ids[i], "nearest_group_id": best_id,
                     "distance": best})
    return rows


def tree_splits(
    ids: list[str], condensed: np.ndarray
) -> tuple[list[dict], dict[tuple[str, int], list[str]]]:
    """UPGMA internal nodes, each as the set of sequence groups below it.

    scipy's 'average' linkage is UPGMA. Every merge is one internal node and so
    one split, named by the lexicographically smallest domain below it -- which
    is stable under any relabelling of the tree's internal nodes, unlike a node
    index.

    That name alone is not a key: a split and its parent can share a smallest
    leaf. The pair (split_id, n_leaves) is one, and provably so. Any two splits
    containing the same leaf are nested, because the ancestors of a leaf form a
    chain under inclusion, so two splits with the same smallest leaf and the
    same size are the same split.
    """
    merged = linkage(condensed, method="average")
    n = len(ids)
    leaves: dict[int, list[str]] = {i: [ids[i]] for i in range(n)}
    rows = []
    members: dict[tuple[str, int], list[str]] = {}
    for step, (left, right, height, _) in enumerate(merged):
        below = sorted(leaves[int(left)] + leaves[int(right)])
        leaves[n + step] = below
        key = (below[0], len(below))
        members[key] = below
        rows.append({
            "split_id": key[0],
            "n_leaves": key[1],
            "height": round(float(height), DECIMALS),
        })
    rows.sort(key=lambda r: (r["split_id"], r["n_leaves"]))
    return rows, members


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
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    check_versions()

    strains = load_strains(args.task / STRAINS)
    species_of = {s["accession"]: " ".join(s["organism"].split()[:2]) for s in strains}

    alphabet = Alphabet.amino()
    models = []
    with pyhmmer.plan7.HMMFile(_hmm_path()) as handle:
        for hmm in handle:
            if str(hmm.accession or "").split(".")[0] == KS_ACCESSION:
                models.append(hmm)
    if len(models) != 1:
        raise SystemExit(f"expected exactly one {KS_ACCESSION} model, found {len(models)}")

    domain_rows, protein_counts = [], {}
    sequence_of: dict[str, str] = {}
    for index, strain in enumerate(sorted(strains, key=lambda s: s["accession"]), start=1):
        accession = strain["accession"]
        contigs = read_contigs(args.task / strain["file"])
        proteins = call_proteins(contigs)
        protein_counts[accession] = len(proteins)
        domains = find_domains(models, proteins, alphabet)
        for number, domain in enumerate(domains, start=1):
            domain_id = f"{accession}_KS{number:03d}"
            sequence_of[domain_id] = domain["sequence"]
            domain_rows.append({
                "domain_id": domain_id,
                "accession": accession,
                "protein_id": domain["protein_id"],
                "start": domain["start"],
                "end": domain["end"],
                "length_aa": domain["end"] - domain["start"] + 1,
                "bitscore": domain["bitscore"],
            })
        print(f"[{index}/{len(strains)}] {accession}: {len(proteins)} proteins, "
              f"{len(domains)} KS domains", flush=True)

    write_csv(args.out / "ks_domains.csv",
              ["domain_id", "accession", "protein_id", "start", "end", "length_aa",
               "bitscore"], domain_rows)

    by_accession = defaultdict(list)
    for row in domain_rows:
        by_accession[row["accession"]].append(row)
    counts = []
    for strain in sorted(strains, key=lambda s: s["accession"]):
        accession = strain["accession"]
        rows = by_accession[accession]
        per_protein = Counter(r["protein_id"] for r in rows)
        counts.append({
            "accession": accession,
            "species": species_of[accession],
            "n_proteins": protein_counts[accession],
            "n_ks_proteins": len(per_protein),
            "n_ks_domains": len(rows),
            "max_domains_per_protein": max(per_protein.values()) if per_protein else 0,
        })
    write_csv(args.out / "genome_ks_counts.csv",
              ["accession", "species", "n_proteins", "n_ks_proteins", "n_ks_domains",
               "max_domains_per_protein"], counts)

    multidomain = []
    for accession, rows in sorted(by_accession.items()):
        per_protein = defaultdict(list)
        for row in rows:
            per_protein[row["protein_id"]].append(row)
        for protein_id, group in sorted(per_protein.items()):
            if len(group) < 2:
                continue
            multidomain.append({
                "accession": accession,
                "protein_id": protein_id,
                "n_domains": len(group),
                "first_domain_id": min(r["domain_id"] for r in group),
                "span_aa": max(r["end"] for r in group) - min(r["start"] for r in group) + 1,
            })
    write_csv(args.out / "multidomain_proteins.csv",
              ["accession", "protein_id", "n_domains", "first_domain_id", "span_aa"],
              multidomain)

    species_rows = defaultdict(list)
    for row in domain_rows:
        species_rows[species_of[row["accession"]]].append(row)
    genomes_per_species = Counter(species_of[s["accession"]] for s in strains)
    write_csv(args.out / "species_domain_counts.csv",
              ["species", "n_genomes", "n_domains", "mean_domains_per_genome"],
              [{"species": species,
                "n_genomes": genomes_per_species[species],
                "n_domains": len(group),
                "mean_domains_per_genome": round(
                    len(group) / genomes_per_species[species], DECIMALS)}
               for species, group in sorted(species_rows.items())])

    groups, domains_of_group = sequence_groups(domain_rows, sequence_of, species_of)
    write_csv(args.out / "sequence_groups.csv",
              ["group_id", "length_aa", "n_domains", "n_genomes", "n_species"], groups)

    ids = [group["group_id"] for group in groups]
    sequences = [sequence_of[i] for i in ids]
    print(f"aligning {len(ids)} distinct sequences from {len(domain_rows)} domains, "
          f"{len(ids) * (len(ids) - 1) // 2} pairs", flush=True)
    condensed = distances(sequences)

    write_csv(args.out / "group_nearest_neighbour.csv",
              ["group_id", "nearest_group_id", "distance"],
              nearest_neighbours(ids, condensed))

    splits, members = tree_splits(ids, condensed)
    write_csv(args.out / "tree_splits.csv", ["split_id", "n_leaves", "height"], splits)

    accession_of = {row["domain_id"]: row["accession"] for row in domain_rows}
    taxonomy = []
    for split in splits:
        below = members[(split["split_id"], split["n_leaves"])]
        # A split's leaves are sequence groups, so the genomes below it are the
        # genomes of every domain in every one of those groups.
        domains = [d for group_id in below for d in domains_of_group[group_id]]
        genomes = {accession_of[d] for d in domains}
        species = {species_of[a] for a in genomes}
        taxonomy.append({
            "split_id": split["split_id"],
            "n_leaves": split["n_leaves"],
            "n_domains": len(domains),
            "n_genomes": len(genomes),
            "n_species": len(species),
            "breadth": "single_species" if len(species) == 1 else "multi_species",
        })
    write_csv(args.out / "split_taxonomy.csv",
              ["split_id", "n_leaves", "n_domains", "n_genomes", "n_species", "breadth"],
              taxonomy)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
