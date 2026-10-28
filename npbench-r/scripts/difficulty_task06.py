"""Measure task 06's difficulty by grading plausible wrong pipelines.

Same approach and caveat as the task 02 to 05 scripts: synthetic variants, a
cheap proxy for the model runs PROJECT.md §8 asks for, useful mainly for
finding artifacts nothing can fail.

Gene calling over 213 genomes takes about half an hour, and three of the
variants change how genes are called, so called proteins are cached on disk per
gene-calling configuration. Everything downstream of the proteins is seconds.

    python scripts/difficulty_task06.py --task tasks/06_streptomyces_ripp \
        --truth tasks/06_streptomyces_ripp/truth --scratch <dir>
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import warnings  # noqa: E402

import pyhmmer  # noqa: E402
import pyrodigal  # noqa: E402
from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402
from pyhmmer.easel import Alphabet, TextSequence  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

MISTAKES = {
    "meta_mode": "calls genes with pyrodigal meta=True instead of training on the assembly",
    "per_contig_training": "trains the gene finder per contig instead of per assembly",
    "loose_evalue": "searches at E=1e-3 instead of 1e-5",
    "all_reported_hits": "keeps every reported hit, not only the ones HMMER includes",
    "window_5kb": "uses a 5 kb neighbourhood instead of 10 kb",
    "no_anchor_merge": "treats every cyclase hit as its own locus, without merging",
    "no_precursor_lower_bound": "takes every ORF up to 120 aa, with no 20 aa floor",
    "class_ignores_rre": "calls a locus complete on cyclase plus synthetase alone",
    "ssn_global_mode": "scores the network with global instead of local alignment",
    "ssn_normalise_by_length": "normalises alignment score by length, not by self-score",
    "families_numbered_by_rank": "names families FAM_0001.. by size rank",
    "no_transitive_closure": "treats each linked pair as a family",
    "drops_zero_hit_families": "omits panel families with no hits from the prevalence table",
}

_PROTEIN_CACHE: dict[tuple, dict[str, list[dict]]] = {}


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref06", task_dir / "solution" / "make_truth.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def call_genes(ref, contigs, mistakes: set[str]) -> list[dict]:
    if "meta_mode" in mistakes:
        finder = pyrodigal.GeneFinder(meta=True)
        trainers = [None]
    elif "per_contig_training" in mistakes:
        finder = None
        trainers = None
    else:
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(*[s for _, s in contigs])
        trainers = [None]

    proteins = []
    for contig, sequence in contigs:
        if "per_contig_training" in mistakes:
            finder = pyrodigal.GeneFinder(meta=False)
            try:
                # Training on a short contig warns; that is the mistake behaving
                # as it would for an agent, not something to fix.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    finder.train(sequence)
            except Exception:
                finder = pyrodigal.GeneFinder(meta=True)
        for index, gene in enumerate(finder.find_genes(sequence), start=1):
            proteins.append({
                "protein_id": f"{contig}_{index}", "contig": contig,
                "start": gene.begin, "end": gene.end,
                "strand": "+" if gene.strand > 0 else "-",
                "sequence": gene.translate().rstrip("*"),
            })
    return proteins


def proteins_for(ref, task_dir: Path, accessions: list[str], mistakes: set[str],
                 cache_dir: Path | None):
    """Called proteins per gene-calling configuration, cached in memory and on disk."""
    key = ("meta_mode" in mistakes, "per_contig_training" in mistakes)
    if key in _PROTEIN_CACHE:
        return _PROTEIN_CACHE[key]

    name = "proteins_" + "".join("1" if k else "0" for k in key) + ".pkl.gz"
    on_disk = (cache_dir / name) if cache_dir else None
    if on_disk is not None and on_disk.is_file():
        print(f"    (reusing cached gene calls for config {key})", flush=True)
        with gzip.open(on_disk, "rb") as handle:
            _PROTEIN_CACHE[key] = pickle.load(handle)
        return _PROTEIN_CACHE[key]

    print(f"    (calling genes for config {key}; this is the slow step)", flush=True)
    built = {}
    for index, accession in enumerate(accessions, start=1):
        contigs = ref.read_contigs(task_dir / "data" / "genomes" / f"{accession}.fna.gz")
        built[accession] = (call_genes(ref, contigs, mistakes),
                            len(contigs), sum(len(s) for _, s in contigs))
        if index % 50 == 0:
            print(f"      {index}/{len(accessions)}", flush=True)
    if on_disk is not None:
        on_disk.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(on_disk, "wb") as handle:
            pickle.dump(built, handle, protocol=5)
    _PROTEIN_CACHE[key] = built
    return built


def scan(ref, models, proteins: list[dict], mistakes: set[str]) -> list[dict]:
    evalue = 1e-3 if "loose_evalue" in mistakes else ref.EVALUE
    alphabet = Alphabet.amino()
    digital = [
        TextSequence(name=p["protein_id"].encode(), sequence=p["sequence"]).digitize(alphabet)
        for p in proteins if p["sequence"]
    ]
    by_id = {p["protein_id"]: p for p in proteins}
    hits = []
    for top in pyhmmer.hmmsearch(models, digital, E=evalue, cpus=4):
        pfam = str(top.query.accession or "").split(".")[0]
        for hit in top:
            if "all_reported_hits" not in mistakes and not hit.included:
                continue
            name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
            protein = by_id[name]
            hits.append({
                "protein_id": name, "pfam": pfam, "pfam_name": ref.PANEL[pfam],
                "contig": protein["contig"], "start": protein["start"],
                "end": protein["end"], "strand": protein["strand"],
                "bitscore": round(float(hit.score), 1),
            })
    return hits


def build_loci(ref, proteins, hits, mistakes: set[str]) -> list[dict]:
    window = 5_000 if "window_5kb" in mistakes else ref.LOCUS_WINDOW_BP
    lower = 1 if "no_precursor_lower_bound" in mistakes else ref.PRECURSOR_MIN_AA
    anchors = sorted((h for h in hits if h["pfam"] == ref.ANCHOR),
                     key=lambda h: (h["contig"], h["start"], h["protein_id"]))
    if "no_anchor_merge" in mistakes:
        groups = [[a] for a in anchors]
    else:
        groups = []
        for anchor in anchors:
            if groups and groups[-1][0]["contig"] == anchor["contig"] and \
                    abs(ref.midpoint(anchor) - ref.midpoint(groups[-1][-1])) <= window:
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
        centres = [ref.midpoint(a) for a in group]
        members = [p for p in by_contig[contig]
                   if any(abs(ref.midpoint(p) - c) <= window for c in centres)]
        found = {pf for p in members for pf in pfams_of.get(p["protein_id"], ())}
        precursors = [p for p in members
                      if lower <= len(p["sequence"]) <= ref.PRECURSOR_MAX_AA]
        if "class_ignores_rre" in mistakes:
            locus_class = "complete" if {ref.ANCHOR, ref.LASSO_SYNTHETASE} <= found \
                else "minimal"
        elif {ref.ANCHOR, ref.LASSO_SYNTHETASE, ref.RRE} <= found:
            locus_class = "complete"
        elif {ref.ANCHOR, ref.LASSO_SYNTHETASE} <= found:
            locus_class = "core"
        else:
            locus_class = "minimal"
        loci.append({
            "locus_id": min(a["protein_id"] for a in group), "contig": contig,
            "start": min(p["start"] for p in members),
            "end": max(p["end"] for p in members),
            "n_genes": len(members), "n_precursors": len(precursors),
            "locus_class": locus_class,
            "has_other_ripp_marker": "true" if found & ref.OTHER_MARKERS else "false",
            "_precursors": precursors, "_pfams": found,
        })
    return loci


def families(ref, precursors: list[dict], mistakes: set[str]) -> dict[str, str]:
    aligner = Align.PairwiseAligner(scoring=None)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = ref.GAP_OPEN
    aligner.extend_gap_score = ref.GAP_EXTEND
    aligner.mode = "global" if "ssn_global_mode" in mistakes else "local"

    labels = sorted(p["precursor_id"] for p in precursors)
    seq = {p["precursor_id"]: p["sequence"] for p in precursors}
    self_score = {k: aligner.score(v, v) for k, v in seq.items()}
    parent = {label: label for label in labels}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    links = []
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            if "ssn_normalise_by_length" in mistakes:
                floor = min(len(seq[a]), len(seq[b]))
            else:
                floor = min(self_score[a], self_score[b])
            if floor <= 0:
                continue
            if aligner.score(seq[a], seq[b]) / floor >= ref.SSN_THRESHOLD:
                links.append((a, b))

    if "no_transitive_closure" in mistakes:
        assignment = {label: label for label in labels}
        for a, b in links:
            assignment[b] = min(assignment[b], a)
    else:
        for a, b in links:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
        assignment = {label: find(label) for label in labels}

    if "families_numbered_by_rank" in mistakes:
        counts = Counter(assignment.values())
        order = sorted(counts, key=lambda f: (-counts[f], f))
        rename = {f: f"FAM_{i + 1:04d}" for i, f in enumerate(order)}
        assignment = {p: rename[f] for p, f in assignment.items()}
    return assignment


def build(task_dir: Path, out: Path, mistakes: set[str], cache_dir: Path | None) -> None:
    ref = load_solution(task_dir)
    with (task_dir / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted(s["accession"] for s in strains)
    species = {s["accession"]: s["species"] for s in strains}
    models = ref.load_panel()
    called = proteins_for(ref, task_dir, accessions, mistakes, cache_dir)

    genome_rows, hit_rows, locus_rows, precursor_rows = [], [], [], []
    for accession in accessions:
        proteins, n_contigs, bases = called[accession]
        hits = scan(ref, models, proteins, mistakes)
        loci = build_loci(ref, proteins, hits, mistakes)
        genome_rows.append({"accession": accession, "species": species[accession],
                            "n_contigs": n_contigs, "assembly_bp": bases,
                            "n_proteins": len(proteins)})
        hit_rows.extend({**h, "accession": accession} for h in hits)
        for locus in loci:
            locus_rows.append({**locus, "accession": accession})
            for precursor in locus["_precursors"]:
                precursor_rows.append({
                    "precursor_id": precursor["protein_id"], "locus_id": locus["locus_id"],
                    "accession": accession, "length_aa": len(precursor["sequence"]),
                    "sequence": precursor["sequence"]})

    hit_rows.sort(key=lambda r: (r["accession"], r["protein_id"], r["pfam"]))
    locus_rows.sort(key=lambda r: r["locus_id"])
    precursor_rows.sort(key=lambda r: r["precursor_id"])

    w = ref.write_csv
    w(out / "genome_proteins.csv",
      ["accession", "species", "n_contigs", "assembly_bp", "n_proteins"], genome_rows)
    w(out / "hmm_hits.csv",
      ["accession", "protein_id", "pfam", "pfam_name", "contig", "start", "end",
       "strand", "bitscore"], hit_rows)
    w(out / "ripp_loci.csv",
      ["locus_id", "accession", "contig", "start", "end", "n_genes", "n_precursors",
       "locus_class", "has_other_ripp_marker"], locus_rows)
    w(out / "precursor_candidates.csv",
      ["precursor_id", "locus_id", "accession", "length_aa", "sequence"], precursor_rows)

    assignment = families(ref, precursor_rows, mistakes)
    size = Counter(assignment.values())
    w(out / "precursor_families.csv", ["precursor_id", "family_id", "family_size"],
      [{"precursor_id": p, "family_id": f, "family_size": size[f]}
       for p, f in sorted(assignment.items())])

    loci_per = Counter(r["accession"] for r in locus_rows)
    precursors_per: dict[str, list[str]] = defaultdict(list)
    for row in precursor_rows:
        precursors_per[row["accession"]].append(row["precursor_id"])
    w(out / "genome_locus_counts.csv",
      ["accession", "species", "n_loci", "n_precursors", "n_families"],
      [{"accession": a, "species": species[a], "n_loci": loci_per.get(a, 0),
        "n_precursors": len(precursors_per.get(a, [])),
        "n_families": len({assignment[p] for p in precursors_per.get(a, [])})}
       for a in accessions])

    hits_per = Counter(r["pfam"] for r in hit_rows)
    genomes_per: dict[str, set[str]] = defaultdict(set)
    for row in hit_rows:
        genomes_per[row["pfam"]].add(row["accession"])
    loci_per_pfam: dict[str, set[str]] = defaultdict(set)
    for locus in locus_rows:
        for pfam in locus["_pfams"]:
            loci_per_pfam[pfam].add(locus["locus_id"])
    panel = sorted(ref.PANEL)
    if "drops_zero_hit_families" in mistakes:
        panel = [p for p in panel if hits_per.get(p, 0)]
    w(out / "domain_prevalence.csv", ["pfam", "pfam_name", "n_hits", "n_genomes", "n_loci"],
      [{"pfam": p, "pfam_name": ref.PANEL[p], "n_hits": hits_per.get(p, 0),
        "n_genomes": len(genomes_per.get(p, ())), "n_loci": len(loci_per_pfam.get(p, ()))}
       for p in panel])


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: pyrodigal meta mode": {"meta_mode"},
    "one slip: trains per contig": {"per_contig_training"},
    "one slip: E=1e-3": {"loose_evalue"},
    "one slip: keeps all reported hits": {"all_reported_hits"},
    "one slip: 5 kb window": {"window_5kb"},
    "one slip: no anchor merge": {"no_anchor_merge"},
    "one slip: no 20 aa precursor floor": {"no_precursor_lower_bound"},
    "one slip: class ignores the RRE": {"class_ignores_rre"},
    "one slip: global alignment in the SSN": {"ssn_global_mode"},
    "one slip: SSN normalised by length": {"ssn_normalise_by_length"},
    "one slip: families numbered by rank": {"families_numbered_by_rank"},
    "one slip: no transitive closure": {"no_transitive_closure"},
    "one slip: drops zero-hit families": {"drops_zero_hit_families"},
    "hasty agent (3 slips)": {"window_5kb", "no_precursor_lower_bound",
                              "class_ignores_rre"},
    "careless agent (5 slips)": {"meta_mode", "loose_evalue", "no_anchor_merge",
                                 "ssn_normalise_by_length", "families_numbered_by_rank"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    config = load_grading_config(args.task / "grading.yaml")

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':40}{'task':>6}  " + "".join(f"{n[:12]:>14}" for n in names)
    print(header)
    print("-" * len(header))
    passes = Counter()
    for label, mistakes in VARIANTS.items():
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace(",", "")
        build(args.task, out, mistakes, args.scratch / "_protein_cache")
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
