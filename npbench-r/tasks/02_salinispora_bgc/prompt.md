# Biosynthetic gene cluster inventory across 120 *Salinispora* genomes

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget about three
hours.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tools named below. Any other resource
  must go through the web access request tool, where it is adjudicated. Requests
  that would reveal published results for this dataset are denied.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  commands and versions, and anything you had to decide for yourself.

The input files are too large to read by eye. Expect to script everything.

## Research objective

*Salinispora* is a marine actinomycete genus known for prolific and
species-structured specialized metabolism. Given 120 *Salinispora* assemblies,
detect the biosynthetic gene clusters in each, and characterise how the cluster
inventory is distributed across contigs, genomes, cluster types and species.

Two properties of the staged data matter for the result.

First, 119 of the 120 assemblies are scaffold-level and only one is a complete
genome. A cluster lying at the edge of a scaffold is truncated, not absent, so
the inventory distinguishes clusters that sit on a contig boundary from those
that do not.

Second, `strains.tsv` records the organism name as NCBI holds it, which is
not consistent: some rows carry the strain designation inside the organism
string and some do not. Species-level grouping has to account for that.

### Staged data

| Path | Contents |
|---|---|
| `data/genomes/*.fna.gz` | 120 gzipped nucleotide FASTA assemblies, named by RefSeq assembly accession |
| `strains.tsv` | one row per assembly: `accession`, `organism`, `strain`, `assembly_level`, `file` |
| `data.manifest.json` | the provenance record for every staged file |

## Method guidance

Detect clusters with **GECCO 0.11.0** (the `gecco` distribution on PyPI; note
that the `gecco-tool` distribution is a renamed shim). Pin the version: cluster
counts and type assignments differ between GECCO releases, and this task is
graded against 0.11.0 output.

Run GECCO's full pipeline once per assembly, with de-novo gene calling:

```
gecco --jobs 4 run --genome <assembly>.fna --output-dir <out>/<accession>
```

Take all defaults. Do not change the probability threshold, the edge distance,
the CRF model, or the postprocessing, and do not pass `--no-pad`, `--no-trim`
or `--cds-feature`. GECCO reads uncompressed FASTA, so decompress each
assembly first. Its per-assembly `*.clusters.tsv` carries one row per detected
cluster with `sequence_id`, `cluster_id`, `start`, `end` and `type`.

Derive the artifacts from that output as follows.

- **Coordinates** are GECCO's own, unchanged, and 1-based inclusive. Cluster
  length is `end - start + 1`.
- **Contig boundary.** A cluster sits on a boundary if its start is at or below
  base 1000 of its contig, or its end is within 1000 bases of the contig's last
  base. Contig lengths come from the FASTA, not from the GECCO table.
- **Cluster type** is GECCO's `type` verbatim, including `Unknown` and
  multi-label values such as `NRP;Polyketide`. Do not split, re-map or drop
  them.
- **Species** is the first two whitespace-separated tokens of the `organism`
  field in `strains.tsv`.
- Every staged assembly appears in `clusters_per_genome.csv`, including any
  assembly in which no cluster is detected.

## Required outputs

Five CSV files in `outputs/`, each with a header row and the exact column names
below. Column order is not graded; column names are. Row order is not graded.

### `bgc_inventory.csv`
One row per detected cluster, over all 120 assemblies.

`cluster_id`, `accession`, `sequence_id`, `start`, `end`, `length_bp`, `type`

`cluster_id` is GECCO's own cluster identifier. `accession` is the RefSeq
assembly accession the cluster came from, as named in `strains.tsv`.

### `contig_edge_clusters.csv`
One row per cluster that sits on a contig boundary, by the rule above.

`cluster_id`, `accession`, `sequence_id`, `contig_length_bp`, `boundary`

`boundary` is `start`, `end`, or `both`.

### `clusters_per_genome.csv`
One row per staged assembly.

`accession`, `n_clusters`

### `type_counts.csv`
One row per distinct cluster type observed.

`type`, `n_clusters`, `n_genomes`

`n_genomes` is the number of distinct assemblies in which that type occurs.

### `species_type_counts.csv`
One row per observed species and cluster type combination.

`species`, `type`, `n_clusters`
