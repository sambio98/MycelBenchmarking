"""Measured difficulty of tasks/08_salinispora_ks_phylogeny.

Every variant is measured on all eight artifacts over the full dataset. That is
only affordable because the tree is built over the 922 distinct domain
sequences rather than all 2,993 domains: 425,000 alignments instead of 4.5
million, four minutes instead of thirty-nine. An earlier version of this script
had to leave the three tree artifacts unmeasured for any variant that changed
the domain set, and say so in the table. Dereplicating the leaves removed the
need.

Raw alignment scores are cached on disk per leaf set, and the extracted domain
tables per domain-layer configuration, so a second run is minutes. The three
distance-normalisation variants are three different functions of one cached
score matrix and cost nothing.

Same cheap-proxy caveat as the other tasks: synthetic variants, no substitute
for the model runs PROJECT.md section 8 asks for.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyhmmer  # noqa: E402
import pyrodigal  # noqa: E402
from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402
from pyhmmer.easel import Alphabet, TextSequence  # noqa: E402
from scipy.cluster.hierarchy import linkage  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

# Variants that leave the domain sequences alone, so the cached matrix holds.
CHEAP = {
    "normalise_by_length": "divides the score by the longer sequence, not the self-score",
    "normalise_by_max_self": "normalises by the larger self-score, not the smaller",
    "no_clamp": "leaves a negative distance negative instead of clamping to 0",
    "single_linkage": "builds the tree with single linkage instead of UPGMA",
    "complete_linkage": "builds the tree with complete linkage instead of UPGMA",
    "splits_numbered": "names splits N0001.. by merge order, not by smallest leaf",
    "species_keeps_strain": "groups by the whole organism string, strain included",
    "breadth_counts_genomes": "calls breadth on the genome count, not the species count",
}

# Variants that change the tree's leaf set rather than the domain set.
LEAF_LAYER = {
    "no_dereplication": "builds the tree over every domain, not over distinct sequences",
    "dereplicate_on_length": "groups domains by length instead of by exact sequence",
}

# Variants that change which domains exist or what their sequences are.
DOMAIN_LAYER = {
    "meta_mode": "calls genes with pyrodigal meta=True instead of training per assembly",
    "per_contig_training": "trains the gene finder per contig",
    "loose_evalue": "searches at E=1e-3 instead of 1e-5",
    "all_reported_hits": "keeps every reported hit and domain, not only the included",
    "envelope_coords": "takes the envelope region rather than the aligned region",
    "no_length_floor": "keeps domains shorter than 150 aa",
    "one_domain_per_protein": "reports only the best-scoring domain in each protein",
}

MISTAKES = {**CHEAP, **LEAF_LAYER, **DOMAIN_LAYER}

def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref_solution_08", task_dir / "solution" / "make_truth.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def extract_domains(
    ref, task_dir: Path, mistakes: set[str]
) -> tuple[list[dict], dict[str, str], dict[str, int]]:
    """The domain table, the sequences and the gene-call size per assembly."""
    strains = ref.load_strains(task_dir / ref.STRAINS)
    alphabet = Alphabet.amino()
    models = []
    with pyhmmer.plan7.HMMFile(ref._hmm_path()) as handle:
        for hmm in handle:
            if str(hmm.accession or "").split(".")[0] == ref.KS_ACCESSION:
                models.append(hmm)

    evalue = 1e-3 if "loose_evalue" in mistakes else ref.EVALUE
    require_included = "all_reported_hits" not in mistakes
    floor = 0 if "no_length_floor" in mistakes else ref.MIN_DOMAIN_AA

    rows, sequence_of, protein_counts = [], {}, {}
    for strain in sorted(strains, key=lambda s: s["accession"]):
        accession = strain["accession"]
        contigs = ref.read_contigs(task_dir / strain["file"])
        if "meta_mode" in mistakes:
            finder = pyrodigal.GeneFinder(meta=True)
            proteins = [
                (f"c{i:04d}_g{j:05d}", gene.translate())
                for i, contig in enumerate(contigs, start=1)
                for j, gene in enumerate(finder.find_genes(contig), start=1)
            ]
        elif "per_contig_training" in mistakes:
            proteins = []
            for i, contig in enumerate(contigs, start=1):
                finder = pyrodigal.GeneFinder()
                try:
                    finder.train(contig)
                except ValueError:       # a contig too short to train on
                    continue
                for j, gene in enumerate(finder.find_genes(contig), start=1):
                    proteins.append((f"c{i:04d}_g{j:05d}", gene.translate()))
        else:
            proteins = ref.call_proteins(contigs)

        protein_counts[accession] = len(proteins)
        sequence_by_id = dict(proteins)
        digital = [
            TextSequence(name=name.encode(), sequence=seq).digitize(alphabet)
            for name, seq in proteins
        ]
        found = []
        for top in pyhmmer.hmmsearch(models, digital, E=evalue, cpus=4):
            for hit in top:
                if require_included and not hit.included:
                    continue
                protein_id = (hit.name.decode() if isinstance(hit.name, bytes)
                              else str(hit.name))
                best = None
                for domain in hit.domains:
                    if require_included and not domain.included:
                        continue
                    alignment = domain.alignment
                    if "envelope_coords" in mistakes:
                        start, end = domain.env_from, domain.env_to
                    else:
                        start, end = alignment.target_from, alignment.target_to
                    if end - start + 1 < floor:
                        continue
                    entry = {
                        "protein_id": protein_id, "start": start, "end": end,
                        "bitscore": round(domain.score, 1),
                        "sequence": sequence_by_id[protein_id][start - 1:end],
                    }
                    if "one_domain_per_protein" in mistakes:
                        if best is None or domain.score > best[0]:
                            best = (domain.score, entry)
                    else:
                        found.append(entry)
                if "one_domain_per_protein" in mistakes and best is not None:
                    found.append(best[1])
        found.sort(key=lambda d: (d["protein_id"], d["start"]))
        for number, domain in enumerate(found, start=1):
            domain_id = f"{accession}_KS{number:03d}"
            sequence_of[domain_id] = domain["sequence"]
            rows.append({
                "domain_id": domain_id, "accession": accession,
                "protein_id": domain["protein_id"], "start": domain["start"],
                "end": domain["end"],
                "length_aa": domain["end"] - domain["start"] + 1,
                "bitscore": domain["bitscore"],
            })
    return rows, sequence_of, protein_counts


def raw_scores(sequences: list[str], cache: Path | None) -> tuple[np.ndarray, np.ndarray]:
    """Pairwise and self alignment scores, cached on disk.

    Separating the raw scores from the distance definition is what makes the
    normalisation variants free: they are three different functions of these
    same numbers.
    """
    if cache is not None and cache.is_file():
        stored = np.load(cache)
        if len(stored["self_scores"]) == len(sequences):
            print(f"    (raw scores from {cache.name})", flush=True)
            return stored["pairs"], stored["self_scores"]

    aligner = Align.PairwiseAligner(scoring=None)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"

    n = len(sequences)
    self_scores = np.array([aligner.score(s, s) for s in sequences], dtype=float)
    pairs = np.empty(n * (n - 1) // 2, dtype=float)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            pairs[k] = aligner.score(sequences[i], sequences[j])
            k += 1
        if (i + 1) % 250 == 0 or i + 1 == n:
            print(f"    aligning: row {i + 1}/{n}", flush=True)
    if cache is not None:
        np.savez(cache, pairs=pairs, self_scores=self_scores)
    return pairs, self_scores


def to_distances(pairs, self_scores, lengths, mistakes: set[str]) -> np.ndarray:
    n = len(self_scores)
    out = np.empty_like(pairs)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if "normalise_by_length" in mistakes:
                floor = float(max(lengths[i], lengths[j]))
            elif "normalise_by_max_self" in mistakes:
                floor = max(self_scores[i], self_scores[j])
            else:
                floor = min(self_scores[i], self_scores[j])
            value = 1.0 - (pairs[k] / floor if floor else 0.0)
            out[k] = value if "no_clamp" in mistakes else min(1.0, max(0.0, value))
            k += 1
    return out


def build(ref, task_dir: Path, out: Path, mistakes: set[str], scratch: Path) -> None:
    """Write one variant's eight artifacts."""
    strains = ref.load_strains(task_dir / ref.STRAINS)
    if "species_keeps_strain" in mistakes:
        species_of = {s["accession"]: s["organism"] for s in strains}
    else:
        species_of = {s["accession"]: " ".join(s["organism"].split()[:2]) for s in strains}

    domain_mistakes = mistakes & set(DOMAIN_LAYER)
    key = "reference" if not domain_mistakes else "+".join(sorted(domain_mistakes))
    cached_rows = scratch / f"domains_{key}.npz"
    if cached_rows.is_file():
        stored = np.load(cached_rows, allow_pickle=True)
        rows = list(stored["rows"])
        sequence_of = dict(stored["sequences"].item())
        protein_counts = dict(stored["protein_counts"].item())
    else:
        print(f"    (extracting domains: {key})", flush=True)
        rows, sequence_of, protein_counts = extract_domains(ref, task_dir, mistakes)
        np.savez(cached_rows, rows=np.array(rows, dtype=object),
                 sequences=np.array(sequence_of, dtype=object),
                 protein_counts=np.array(protein_counts, dtype=object))

    w = ref.write_csv
    w(out / "ks_domains.csv",
      ["domain_id", "accession", "protein_id", "start", "end", "length_aa", "bitscore"], rows)

    by_accession = defaultdict(list)
    for row in rows:
        by_accession[row["accession"]].append(row)
    counts = []
    for strain in sorted(strains, key=lambda s: s["accession"]):
        accession = strain["accession"]
        group = by_accession[accession]
        per_protein = Counter(r["protein_id"] for r in group)
        counts.append({
            "accession": accession, "species": species_of[accession],
            "n_proteins": protein_counts[accession],
            "n_ks_proteins": len(per_protein), "n_ks_domains": len(group),
            "max_domains_per_protein": max(per_protein.values()) if per_protein else 0,
        })
    w(out / "genome_ks_counts.csv",
      ["accession", "species", "n_proteins", "n_ks_proteins", "n_ks_domains",
       "max_domains_per_protein"], counts)

    multidomain = []
    for accession, group in sorted(by_accession.items()):
        per_protein = defaultdict(list)
        for row in group:
            per_protein[row["protein_id"]].append(row)
        for protein_id, members in sorted(per_protein.items()):
            if len(members) < 2:
                continue
            multidomain.append({
                "accession": accession, "protein_id": protein_id,
                "n_domains": len(members),
                "first_domain_id": min(r["domain_id"] for r in members),
                "span_aa": max(r["end"] for r in members) - min(r["start"] for r in members) + 1,
            })
    w(out / "multidomain_proteins.csv",
      ["accession", "protein_id", "n_domains", "first_domain_id", "span_aa"], multidomain)

    species_rows = defaultdict(list)
    for row in rows:
        species_rows[species_of[row["accession"]]].append(row)
    genomes_per_species = Counter(species_of[s["accession"]] for s in strains)
    w(out / "species_domain_counts.csv",
      ["species", "n_genomes", "n_domains", "mean_domains_per_genome"],
      [{"species": s, "n_genomes": genomes_per_species[s], "n_domains": len(g),
        "mean_domains_per_genome": round(len(g) / genomes_per_species[s], ref.DECIMALS)}
       for s, g in sorted(species_rows.items())])

    # The leaf set. Grouping on length instead of on the sequence is a modelled
    # mistake; not grouping at all is the one the prompt warns against, and it
    # is the expensive variant because it puts 2,993 leaves in the tree.
    accession_of = {r["domain_id"]: r["accession"] for r in rows}
    if "no_dereplication" in mistakes:
        by_group = {r["domain_id"]: [r["domain_id"]] for r in rows}
    elif "dereplicate_on_length" in mistakes:
        buckets = defaultdict(list)
        for row in rows:
            buckets[row["length_aa"]].append(row["domain_id"])
        by_group = {sorted(v)[0]: sorted(v) for v in buckets.values()}
    else:
        buckets = defaultdict(list)
        for row in rows:
            buckets[sequence_of[row["domain_id"]]].append(row["domain_id"])
        by_group = {sorted(v)[0]: sorted(v) for v in buckets.values()}

    groups = []
    for group_id, domain_ids in by_group.items():
        genomes = {accession_of[d] for d in domain_ids}
        groups.append({
            "group_id": group_id,
            "length_aa": len(sequence_of[group_id]),
            "n_domains": len(domain_ids),
            "n_genomes": len(genomes),
            "n_species": len({species_of[a] for a in genomes}),
        })
    groups.sort(key=lambda g: g["group_id"])
    w(out / "sequence_groups.csv",
      ["group_id", "length_aa", "n_domains", "n_genomes", "n_species"], groups)

    ids = [g["group_id"] for g in groups]
    sequences = [sequence_of[i] for i in ids]
    lengths = [len(s) for s in sequences]
    leaf_key = f"{key}_{'derep' if not (mistakes & set(LEAF_LAYER)) else '+'.join(sorted(mistakes & set(LEAF_LAYER)))}"
    print(f"    (scores for {leaf_key}: {len(ids)} leaves, "
          f"{len(ids) * (len(ids) - 1) // 2} pairs)", flush=True)
    pairs, self_scores = raw_scores(sequences, scratch / f"scores_{leaf_key}.npz")
    condensed = to_distances(pairs, self_scores, lengths, mistakes)

    w(out / "group_nearest_neighbour.csv",
      ["group_id", "nearest_group_id", "distance"],
      ref.nearest_neighbours(ids, condensed))

    method = ("single" if "single_linkage" in mistakes
              else "complete" if "complete_linkage" in mistakes else "average")
    merged = linkage(condensed, method=method)
    n = len(ids)
    leaves: dict[int, list[str]] = {i: [ids[i]] for i in range(n)}
    splits, members = [], {}
    for step, (left, right, height, _) in enumerate(merged):
        below = sorted(leaves[int(left)] + leaves[int(right)])
        leaves[n + step] = below
        name = f"N{step + 1:04d}" if "splits_numbered" in mistakes else below[0]
        members[(name, len(below))] = below
        splits.append({"split_id": name, "n_leaves": len(below),
                       "height": round(float(height), ref.DECIMALS)})
    splits.sort(key=lambda r: (r["split_id"], r["n_leaves"]))
    w(out / "tree_splits.csv", ["split_id", "n_leaves", "height"], splits)

    taxonomy = []
    for split in splits:
        below = members[(split["split_id"], split["n_leaves"])]
        domains = [d for group_id in below for d in by_group[group_id]]
        genomes = {accession_of[d] for d in domains}
        species = {species_of[a] for a in genomes}
        if "breadth_counts_genomes" in mistakes:
            breadth = "single_species" if len(genomes) == 1 else "multi_species"
        else:
            breadth = "single_species" if len(species) == 1 else "multi_species"
        taxonomy.append({"split_id": split["split_id"], "n_leaves": split["n_leaves"],
                         "n_domains": len(domains), "n_genomes": len(genomes),
                         "n_species": len(species), "breadth": breadth})
    w(out / "split_taxonomy.csv",
      ["split_id", "n_leaves", "n_domains", "n_genomes", "n_species", "breadth"],
      taxonomy)


VARIANTS = {
    "reference (no mistakes)": set(),
    # Cheap first, so a partial run still says something about every layer.
    "one slip: normalise by length": {"normalise_by_length"},
    "one slip: normalise by max self-score": {"normalise_by_max_self"},
    "one slip: no clamp on the distance": {"no_clamp"},
    "one slip: single linkage": {"single_linkage"},
    "one slip: complete linkage": {"complete_linkage"},
    "one slip: splits numbered": {"splits_numbered"},
    "one slip: species keeps the strain": {"species_keeps_strain"},
    "one slip: breadth counts genomes": {"breadth_counts_genomes"},
    "one slip: groups by length, not sequence": {"dereplicate_on_length"},
    "one slip: no dereplication at all": {"no_dereplication"},
    "hasty agent (3 slips)": {"complete_linkage", "splits_numbered",
                              "breadth_counts_genomes"},
    "one slip: meta-mode gene calling": {"meta_mode"},
    "one slip: trains per contig": {"per_contig_training"},
    "one slip: E=1e-3": {"loose_evalue"},
    "one slip: keeps all reported hits": {"all_reported_hits"},
    "one slip: envelope coordinates": {"envelope_coords"},
    "one slip: no 150 aa floor": {"no_length_floor"},
    "one slip: one domain per protein": {"one_domain_per_protein"},
    "careless agent (5 slips)": {"meta_mode", "loose_evalue", "no_length_floor",
                                 "splits_numbered", "species_keeps_strain"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--only", default=None,
                        help="substring filter over variant labels")
    args = parser.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)
    config = load_grading_config(args.task / "grading.yaml")
    ref = load_solution(args.task)

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':38}{'task':>6}  " + "".join(f"{n[:13]:>15}" for n in names)
    print(header)
    print("-" * len(header))

    passes, measured = Counter(), Counter()
    for label, mistakes in VARIANTS.items():
        if args.only and args.only not in label:
            continue
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace("/", "")
        build(ref, args.task, out, mistakes, args.scratch)
        result = score_task(config, args.truth, out)
        scored = {a.path: a for a in result.artifacts + result.soft_artifacts}
        row = f"{label:38}{result.task_score:6.2f}  "
        for spec in config.artifacts:
            artifact = scored[spec.path]
            row += f"{artifact.score:15.2f}"
            measured[spec.path] += 1
            if artifact.passed:
                passes[spec.path] += 1
        print(row, flush=True)

    print()
    total = max(measured.values()) if measured else 0
    print(f"artifact pass rate across the {total} variants "
          f"(PROJECT.md section 8 wants 0.1 to 0.8):")
    for artifact in config.artifacts:
        n = measured[artifact.path]
        rate = passes[artifact.path] / n if n else 0.0
        verdict = ("soft, out of the headline" if artifact.soft
                   else "ok" if 0.1 <= rate <= 0.8 else "OUTSIDE BAND")
        print(f"  {artifact.path:34} {rate:4.2f}  {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
