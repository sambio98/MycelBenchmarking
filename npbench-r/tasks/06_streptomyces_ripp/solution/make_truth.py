"""Reference solution for tasks/06_streptomyces_ripp.

RODEO-style lasso peptide and RiPP locus mining across Streptomyces reference
genomes, PROJECT.md §4 category 6 and the RiPP half of §5's task 6.

    python solution/make_truth.py --task tasks/06_streptomyces_ripp --out <dir>

RODEO itself is not installable offline, so the pipeline is rebuilt from its
parts: pyrodigal for gene calling, pyhmmer against a fixed panel of RiPP Pfam
models, a neighbourhood rule anchored on the lasso cyclase, and a sequence
similarity network over the precursor candidates. The Pfam models ship inside
the pinned environment (GECCO bundles a 2,766-model Pfam 35.0 subset), so
nothing is fetched at runtime.

Artifacts share upstream dependencies on purpose, following tasks 04 and 05:
the gene call feeds the HMM hits, which feed the loci and the domain
prevalence, and the loci feed the precursors, which feed the families and the
per-genome counts. A mistake in the gene call reaches all seven.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import pyhmmer
import pyrodigal
from Bio import Align
from Bio.Align import substitution_matrices
from pyhmmer.easel import Alphabet, TextSequence

def _hmm_path() -> Path:
    """The Pfam models bundled inside the pinned GECCO install.

    Found through the package rather than an absolute path, so the task does
    not break on a different Python version or prefix. GECCO ships a
    2,766-model subset of Pfam 35.0 with a recorded md5, which is why this task
    needs no database download and nothing is fetched at runtime.
    """
    import gecco

    path = Path(gecco.__file__).parent / "hmmer" / "Pfam.h3m"
    if not path.is_file():
        raise SystemExit(f"the pinned GECCO install has no Pfam models at {path}")
    return path

# The RiPP marker panel, fixed. PF13471 is the anchor: it is the lasso peptide
# cyclase and the marker RODEO keys on. The rest describe what else the
# neighbourhood carries.
ANCHOR = "PF13471"
LASSO_SYNTHETASE = "PF00733"
RRE = "PF05402"
PANEL = {
    "PF13471": "Transglut_core3",
    "PF00733": "Asn_synthase",
    "PF05402": "PqqD",
    "PF04055": "Radical_SAM",
    "PF00881": "Nitroreductase",
    "PF02052": "Gallidermin",
    "PF08130": "Antimicrobial18",
    "PF05147": "LANC_like",
    "PF04738": "Lant_dehydr_N",
    "PF14028": "Lant_dehydr_C",
    "PF02624": "YcaO",
    "PF10439": "Bacteriocin_IIc",
    "PF07812": "TfuA",
}
OTHER_MARKERS = frozenset(PANEL) - {ANCHOR, LASSO_SYNTHETASE, RRE}

EVALUE = 1e-5
LOCUS_WINDOW_BP = 10_000
PRECURSOR_MIN_AA = 20
PRECURSOR_MAX_AA = 120

# Sequence similarity network over precursors. Local alignment scores only, no
# traceback: a score is deterministic across Biopython versions in a way the
# chosen alignment path is not.
SSN_THRESHOLD = 0.40
GAP_OPEN = -11
GAP_EXTEND = -1


def read_contigs(archive: Path) -> list[tuple[str, str]]:
    """(name, sequence) per contig, in file order. Order is part of the method."""
    contigs = []
    name = None
    buffer: list[str] = []
    with gzip.open(archive, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, "".join(buffer)))
                name = line[1:].split()[0]
                buffer = []
            else:
                buffer.append(line.strip())
    if name is not None:
        contigs.append((name, "".join(buffer)))
    return contigs


def call_genes(contigs: list[tuple[str, str]]) -> list[dict]:
    """pyrodigal in single-genome mode, trained on the whole assembly."""
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*[sequence for _, sequence in contigs])
    proteins = []
    for contig, sequence in contigs:
        for index, gene in enumerate(finder.find_genes(sequence), start=1):
            proteins.append(
                {
                    "protein_id": f"{contig}_{index}",
                    "contig": contig,
                    "start": gene.begin,
                    "end": gene.end,
                    "strand": "+" if gene.strand > 0 else "-",
                    "sequence": gene.translate().rstrip("*"),
                }
            )
    return proteins


def load_panel() -> list:
    models = []
    with pyhmmer.plan7.HMMFile(_hmm_path()) as handle:
        for hmm in handle:
            if str(hmm.accession or "").split(".")[0] in PANEL:
                models.append(hmm)
    missing = set(PANEL) - {str(m.accession or "").split(".")[0] for m in models}
    if missing:
        raise SystemExit(
            f"the pinned Pfam subset is missing {sorted(missing)}; truth was built "
            f"against a build that has all {len(PANEL)} panel families"
        )
    return models


def scan(models, proteins: list[dict]) -> list[dict]:
    """hmmsearch the panel over one genome's proteins. Returns included hits."""
    alphabet = Alphabet.amino()
    digital = [
        TextSequence(name=p["protein_id"].encode(), sequence=p["sequence"]).digitize(alphabet)
        for p in proteins
        if p["sequence"]
    ]
    by_id = {p["protein_id"]: p for p in proteins}
    hits = []
    for top in pyhmmer.hmmsearch(models, digital, E=EVALUE, cpus=4):
        pfam = str(top.query.accession or "").split(".")[0]
        for hit in top:
            if not hit.included:
                continue
            name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
            protein = by_id[name]
            hits.append(
                {
                    "protein_id": name,
                    "pfam": pfam,
                    "pfam_name": PANEL[pfam],
                    "contig": protein["contig"],
                    "start": protein["start"],
                    "end": protein["end"],
                    "strand": protein["strand"],
                    "bitscore": round(float(hit.score), 1),
                }
            )
    return hits


def midpoint(record: dict) -> float:
    return (record["start"] + record["end"]) / 2


def build_loci(proteins: list[dict], hits: list[dict]) -> list[dict]:
    """Loci anchored on the lasso cyclase.

    Anchors within one window of each other on the same contig merge into a
    single locus, named for the lexicographically smallest anchor protein, so
    a tandem pair does not produce two overlapping loci.
    """
    anchors = sorted(
        (h for h in hits if h["pfam"] == ANCHOR),
        key=lambda h: (h["contig"], h["start"], h["protein_id"]),
    )
    groups: list[list[dict]] = []
    for anchor in anchors:
        if groups and groups[-1][0]["contig"] == anchor["contig"] and \
                abs(midpoint(anchor) - midpoint(groups[-1][-1])) <= LOCUS_WINDOW_BP:
            groups[-1].append(anchor)
        else:
            groups.append([anchor])

    by_contig: dict[str, list[dict]] = defaultdict(list)
    for protein in proteins:
        by_contig[protein["contig"]].append(protein)
    pfams_of: dict[str, set[str]] = defaultdict(set)
    for hit in hits:
        pfams_of[hit["protein_id"]].add(hit["pfam"])

    loci = []
    for group in groups:
        contig = group[0]["contig"]
        centres = [midpoint(a) for a in group]
        members = [
            p for p in by_contig[contig]
            if any(abs(midpoint(p) - c) <= LOCUS_WINDOW_BP for c in centres)
        ]
        found = {pfam for p in members for pfam in pfams_of.get(p["protein_id"], ())}
        precursors = [
            p for p in members
            if PRECURSOR_MIN_AA <= len(p["sequence"]) <= PRECURSOR_MAX_AA
        ]
        if {ANCHOR, LASSO_SYNTHETASE, RRE} <= found:
            locus_class = "complete"
        elif {ANCHOR, LASSO_SYNTHETASE} <= found:
            locus_class = "core"
        else:
            locus_class = "minimal"
        loci.append(
            {
                "locus_id": min(a["protein_id"] for a in group),
                "contig": contig,
                "start": min(p["start"] for p in members),
                "end": max(p["end"] for p in members),
                "n_genes": len(members),
                "n_precursors": len(precursors),
                "locus_class": locus_class,
                "has_other_ripp_marker": "true" if found & OTHER_MARKERS else "false",
                "_precursors": precursors,
                "_pfams": found,
            }
        )
    return loci


def similarity_families(precursors: list[dict]) -> dict[str, str]:
    """Single-linkage families over normalised local alignment score."""
    aligner = Align.PairwiseAligner(scoring=None)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND
    aligner.mode = "local"

    labels = sorted(p["precursor_id"] for p in precursors)
    sequence_of = {p["precursor_id"]: p["sequence"] for p in precursors}
    self_score = {k: aligner.score(v, v) for k, v in sequence_of.items()}

    parent = {label: label for label in labels}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            floor = min(self_score[a], self_score[b])
            if floor <= 0:
                continue
            if aligner.score(sequence_of[a], sequence_of[b]) / floor >= SSN_THRESHOLD:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    return {label: find(label) for label in labels}


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
    args = parser.parse_args()

    with (args.task / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted(s["accession"] for s in strains)
    species = {s["accession"]: s["species"] for s in strains}
    models = load_panel()
    print(f"{len(accessions)} genomes, {len(models)} RiPP Pfam models", flush=True)

    genome_rows, hit_rows, locus_rows, precursor_rows = [], [], [], []
    for index, accession in enumerate(accessions, start=1):
        contigs = read_contigs(args.task / "data" / "genomes" / f"{accession}.fna.gz")
        proteins = call_genes(contigs)
        hits = scan(models, proteins)
        loci = build_loci(proteins, hits)

        genome_rows.append({
            "accession": accession, "species": species[accession],
            "n_contigs": len(contigs),
            "assembly_bp": sum(len(s) for _, s in contigs),
            "n_proteins": len(proteins),
        })
        hit_rows.extend({**h, "accession": accession} for h in hits)
        for locus in loci:
            locus_rows.append({**locus, "accession": accession})
            for precursor in locus["_precursors"]:
                precursor_rows.append({
                    "precursor_id": precursor["protein_id"],
                    "locus_id": locus["locus_id"],
                    "accession": accession,
                    "length_aa": len(precursor["sequence"]),
                    "sequence": precursor["sequence"],
                })
        if index % 25 == 0 or index == len(accessions):
            print(f"  {index}/{len(accessions)} "
                  f"({len(locus_rows)} loci, {len(precursor_rows)} precursors)", flush=True)

    hit_rows.sort(key=lambda r: (r["accession"], r["protein_id"], r["pfam"]))
    locus_rows.sort(key=lambda r: r["locus_id"])
    precursor_rows.sort(key=lambda r: r["precursor_id"])

    write_csv(args.out / "genome_proteins.csv",
              ["accession", "species", "n_contigs", "assembly_bp", "n_proteins"], genome_rows)
    write_csv(args.out / "hmm_hits.csv",
              ["accession", "protein_id", "pfam", "pfam_name", "contig", "start", "end",
               "strand", "bitscore"], hit_rows)
    write_csv(args.out / "ripp_loci.csv",
              ["locus_id", "accession", "contig", "start", "end", "n_genes",
               "n_precursors", "locus_class", "has_other_ripp_marker"], locus_rows)
    write_csv(args.out / "precursor_candidates.csv",
              ["precursor_id", "locus_id", "accession", "length_aa", "sequence"],
              precursor_rows)

    families = similarity_families(precursor_rows)
    size = Counter(families.values())
    print(f"{len(size)} precursor families over {len(families)} precursors", flush=True)
    write_csv(args.out / "precursor_families.csv",
              ["precursor_id", "family_id", "family_size"],
              [{"precursor_id": p, "family_id": f, "family_size": size[f]}
               for p, f in sorted(families.items())])

    loci_per_genome = Counter(r["accession"] for r in locus_rows)
    precursors_per_genome: dict[str, list[str]] = defaultdict(list)
    for row in precursor_rows:
        precursors_per_genome[row["accession"]].append(row["precursor_id"])
    write_csv(args.out / "genome_locus_counts.csv",
              ["accession", "species", "n_loci", "n_precursors", "n_families"],
              [{"accession": a, "species": species[a], "n_loci": loci_per_genome.get(a, 0),
                "n_precursors": len(precursors_per_genome.get(a, [])),
                "n_families": len({families[p] for p in precursors_per_genome.get(a, [])})}
               for a in accessions])

    hits_per_pfam = Counter(r["pfam"] for r in hit_rows)
    genomes_per_pfam: dict[str, set[str]] = defaultdict(set)
    for row in hit_rows:
        genomes_per_pfam[row["pfam"]].add(row["accession"])
    loci_per_pfam: dict[str, set[str]] = defaultdict(set)
    for locus in locus_rows:
        for pfam in locus["_pfams"]:
            loci_per_pfam[pfam].add(locus["locus_id"])
    write_csv(args.out / "domain_prevalence.csv",
              ["pfam", "pfam_name", "n_hits", "n_genomes", "n_loci"],
              [{"pfam": p, "pfam_name": PANEL[p], "n_hits": hits_per_pfam.get(p, 0),
                "n_genomes": len(genomes_per_pfam.get(p, ())),
                "n_loci": len(loci_per_pfam.get(p, ()))}
               for p in sorted(PANEL)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
