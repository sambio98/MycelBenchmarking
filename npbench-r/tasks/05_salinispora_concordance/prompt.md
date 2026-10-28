# Taxonomy versus chemistry concordance across 120 *Salinispora* genomes

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget about three
hours; cluster detection is most of it.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tools named below. Any other resource
  must go through the web access request tool, where it is adjudicated.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  commands and versions, and anything you had to decide for yourself.

The inputs are far too large to read by eye. Expect to script everything, and
expect the pairwise step to need care: there are 7,140 genome pairs.

## Research objective

*Salinispora* species differ in their specialized metabolism, and the standing
question is how tightly that chemistry tracks the phylogeny: do strains that
are genomically close carry the same biosynthetic gene cluster families?

Given 120 *Salinispora* assemblies, build two independent distance matrices over
the same strains — one genomic, one chemical — and characterise how far they
agree, overall, within species, and across the range of genomic distance.

Two properties of the data matter. 119 of the 120 assemblies are
scaffold-level, so cluster detection is working on fragmented input. And
`strains.tsv` records the organism name as NCBI holds it, inconsistently: some
rows carry the strain designation inside the organism string and some do not,
so species-level grouping has to account for that.

### Staged data

| Path | Contents |
|---|---|
| `data/genomes/*.fna.gz` | 120 gzipped nucleotide FASTA assemblies, named by RefSeq assembly accession |
| `strains.tsv` | one row per assembly: `accession`, `organism`, `strain`, `assembly_level`, `file` |
| `data.manifest.json` | the provenance record for every staged file |

## Method guidance

Every rule below is part of the method; none of it is open to choice.

**Species** is the first two whitespace-separated tokens of the `organism`
field in `strains.tsv`.

### Cluster detection and families

Detect clusters with **GECCO 0.11.0** (the `gecco` distribution on PyPI; the
`gecco-tool` distribution is a renamed shim). Pin the version: cluster counts
and type assignments differ between releases.

```
gecco --jobs 4 run --genome <assembly>.fna --output-dir <out>/<accession>
```

Take all defaults; decompress each assembly first.

Group the clusters into families exactly as follows:

1. A cluster's **domain set** is the distinct Pfam accessions in its GECCO
   `domains` field, which lists a domain once per occurrence.
2. Clusters with **fewer than 3 distinct domains** are not assignable and
   appear in no family.
3. Two assignable clusters are **linked** when the Jaccard similarity of their
   domain sets is **0.50 or greater**.
4. A **family** is a connected component of that graph — single linkage.
5. A family's **`family_id` is the lexicographically smallest `cluster_id`
   among its members.** Do not number families.

A genome's **repertoire** is the set of families its assignable clusters fall
into.

### The genomic sketch

Do **not** use a MinHash: its distances depend on the hash function, and this
task is graded against an exact rule instead.

1. Read each assembly's sequences, uppercased. Encode bases `A`,`C`,`G`,`T` as
   `0`,`1`,`2`,`3`.
2. Consider every 21-mer containing no character outside `ACGT`.
3. A 21-mer's **encoding** is its bases read left to right as base-4 digits,
   most significant first. Its **canonical encoding** is the smaller of its own
   encoding and that of its reverse complement.
4. **Keep** a 21-mer when its canonical encoding is divisible by **1000**.
5. A genome's sketch is the **set** of kept canonical encodings — duplicates
   collapse.

### The two distances

For an unordered pair of genomes, both distances are one minus a Jaccard
index, **rounded to 6 decimal places**:

- **`genomic_distance`** over the two k-mer sketches.
- **`chemistry_distance`** over the two family repertoires. A pair whose
  repertoires are both empty has distance 1.0.

Every unordered pair of distinct assemblies appears once, with
`accession_a` < `accession_b` lexicographically.

### The bins

Sort all pairs by `genomic_distance` ascending, breaking ties by
`accession_a` then `accession_b`. Split the sorted list into **10 bins of equal
size**, numbered 1 to 10, bin 1 holding the smallest distances. Any remainder
goes into bin 10.

## Required outputs

Seven CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded. Booleans are the lowercase strings `true` and `false`. Means are
rounded to 6 decimal places.

### `genome_kmer_profiles.csv`
One row per assembly.

`accession`, `n_contigs`, `assembly_bp`, `n_kept_kmers`

### `bgc_families.csv`
One row per assignable cluster.

`cluster_id`, `accession`, `family_id`, `family_size`

### `genome_repertoires.csv`
One row per assembly, including any with no clusters.

`accession`, `n_clusters`, `n_families`, `n_unique_families`,
`largest_family_id`

`n_clusters` counts every detected cluster, assignable or not.
`n_unique_families` counts the families in that genome's repertoire that occur
in no other genome. `largest_family_id` is the `family_id` of the most populous
family in that genome's repertoire, counting members across all genomes, with
ties going to the lexicographically smallest `family_id`; it is the empty
string for a genome with an empty repertoire.

### `family_taxonomic_breadth.csv`
One row per family.

`family_id`, `n_clusters`, `n_genomes`, `n_species`, `breadth`,
`representative_cluster`

`breadth` is `single_species` when the family's genomes are all one species and
`multi_species` otherwise. `representative_cluster` is the lexicographically
smallest member `cluster_id`.

### `genome_pair_distances.csv`
One row per unordered pair of distinct assemblies: 7,140 rows.

`accession_a`, `accession_b`, `genomic_distance`, `chemistry_distance`,
`same_species`

### `species_concordance.csv`
One row per species having **two or more** assemblies.

`species`, `n_genomes`, `n_pairs`, `mean_genomic_distance`,
`mean_chemistry_distance`

Means are over that species' within-species pairs only.

### `distance_concordance_bins.csv`
One row per bin.

`bin`, `mean_genomic_distance`, `mean_chemistry_distance`,
`n_same_species_pairs`
