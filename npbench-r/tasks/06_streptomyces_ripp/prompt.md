# Lasso peptide and RiPP locus mining across 213 *Streptomyces* reference genomes

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget two hours;
gene calling is most of it.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tools named below. Any other resource
  must go through the web access request tool, where it is adjudicated. **No
  database download is needed**: the HMM models ship inside the pinned
  environment.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  versions, and anything you had to decide for yourself.

There are 213 genomes and about 1.5 million proteins. Everything here has to be
scripted.

## Research objective

Lasso peptides are ribosomally synthesised and post-translationally modified
peptides (RiPPs) whose defining biosynthetic marker is a transglutaminase-like
cyclase. A locus is recognised by that cyclase together with a macrolactam
synthetase and, usually, a RiPP recognition element; the peptide itself is a
short precursor open reading frame in the same neighbourhood, too small and too
variable to be found by similarity alone.

Given 213 *Streptomyces* reference genomes, mine the lasso peptide and
associated RiPP loci: call the genes, find the marker domains, delimit the
loci, extract the precursor candidates, and group those precursors into a
sequence similarity network.

### Staged data

| Path | Contents |
|---|---|
| `data/genomes/*.fna.gz` | 213 gzipped nucleotide FASTA assemblies, named by RefSeq assembly accession |
| `strains.tsv` | one row per assembly: `accession`, `organism`, `species`, `strain`, `file` |
| `data.manifest.json` | the provenance record for every staged file |

The set is every NCBI RefSeq **reference genome** of a named *Streptomyces*
species with a complete, current assembly: exactly one genome per species, 213
species. `species` in `strains.tsv` is already the two-token species name; use
it as given rather than re-deriving it.

## Method guidance

Use **pyrodigal 3.7.1**, **pyhmmer 0.12.3** and **Biopython 1.88**. The Pfam
models come from the **GECCO 0.11.0** install, which bundles a 2,766-model
subset of **Pfam 35.0** at `<gecco package>/hmmer/Pfam.h3m`. Pin all four:
gene boundaries, bit scores and alignment scores are all version-sensitive.

Every rule below is part of the method; none of it is open to choice.

**Gene calling.** One `pyrodigal.GeneFinder(meta=False)` per assembly, trained
on **all** of that assembly's contigs in file order, then genes found per
contig in file order. A protein's `protein_id` is `{contig}_{index}` where
`index` is 1-based within its contig, in the order pyrodigal returns genes.
Translate with pyrodigal and strip a trailing `*`.

**The HMM panel.** These thirteen Pfam families, by accession:

| | | |
|---|---|---|
| `PF13471` Transglut_core3 | `PF00733` Asn_synthase | `PF05402` PqqD |
| `PF04055` Radical_SAM | `PF00881` Nitroreductase | `PF02052` Gallidermin |
| `PF08130` Antimicrobial18 | `PF05147` LANC_like | `PF04738` Lant_dehydr_N |
| `PF14028` Lant_dehydr_C | `PF02624` YcaO | `PF10439` Bacteriocin_IIc |
| `PF07812` TfuA | | |

Search them against each genome's proteins with `pyhmmer.hmmsearch` at
`E=1e-5`, and keep only hits HMMER itself includes (`hit.included`). Report the
per-sequence bit score **rounded to 1 decimal place**.

**Loci.** `PF13471` is the anchor. For each anchor hit, the locus is every
called gene on the same contig whose **midpoint** lies within **10,000 bp** of
the anchor gene's midpoint, where a midpoint is `(start + end) / 2`.

Two anchors on the same contig whose midpoints are within 10,000 bp of each
other **merge into one locus**, and the merged locus takes every gene within
10,000 bp of **either** anchor. A locus's `locus_id` is the lexicographically
smallest anchor `protein_id` it contains. Its `start` and `end` are the minimum
start and maximum end over its member genes.

**`locus_class`** is `complete` when the locus carries all three of `PF13471`,
`PF00733` and `PF05402`; `core` when it carries `PF13471` and `PF00733` but not
`PF05402`; and `minimal` otherwise.

**`has_other_ripp_marker`** is `true` when the locus carries any panel family
other than those three.

**Precursor candidates** are the locus's member genes whose translated length
is **between 20 and 120 amino acids inclusive**.

**The similarity network**, over all precursor candidates from all genomes:

1. Score pairs with a Biopython `PairwiseAligner` in **local** mode, substitution
   matrix **BLOSUM62**, `open_gap_score = -11`, `extend_gap_score = -1`.
2. Normalise: `score(a, b) / min(score(a, a), score(b, b))`. Skip a pair whose
   smaller self-score is not positive.
3. Two precursors are **linked** when that normalised score is **0.40 or
   greater**.
4. A **family** is a connected component of that graph — single linkage.
5. A family's **`family_id` is the lexicographically smallest `precursor_id`
   among its members.** Do not number families.

Use alignment **scores** only. Do not depend on the alignment path: scores are
reproducible across Biopython builds in a way the chosen traceback is not.

## Required outputs

Seven CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded. Booleans are the lowercase strings `true` and `false`. Coordinates are
pyrodigal's own, unchanged.

### `genome_proteins.csv`
One row per assembly.

`accession`, `species`, `n_contigs`, `assembly_bp`, `n_proteins`

### `hmm_hits.csv`
One row per (assembly, protein, Pfam family) hit that HMMER includes.

`accession`, `protein_id`, `pfam`, `pfam_name`, `contig`, `start`, `end`,
`strand`, `bitscore`

`pfam` is the bare accession without its version suffix, e.g. `PF13471`.
`strand` is `+` or `-`.

### `ripp_loci.csv`
One row per locus.

`locus_id`, `accession`, `contig`, `start`, `end`, `n_genes`, `n_precursors`,
`locus_class`, `has_other_ripp_marker`

### `precursor_candidates.csv`
One row per precursor candidate.

`precursor_id`, `locus_id`, `accession`, `length_aa`, `sequence`

`precursor_id` is the gene's `protein_id`.

### `precursor_families.csv`
One row per precursor candidate.

`precursor_id`, `family_id`, `family_size`

### `genome_locus_counts.csv`
One row per assembly, including any with no loci.

`accession`, `species`, `n_loci`, `n_precursors`, `n_families`

`n_families` counts the distinct families that assembly's precursors fall into.

### `domain_prevalence.csv`
One row per panel family, **all thirteen**, including any with no hits.

`pfam`, `pfam_name`, `n_hits`, `n_genomes`, `n_loci`

`n_loci` counts the loci carrying at least one hit of that family.
