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

### Gene cluster families

Group the clusters into families by domain composition, as follows. Every rule
here is part of the method; none of it is open to choice.

1. A cluster's **domain set** is the distinct Pfam accessions in its GECCO
   `domains` field. That field lists a domain once per occurrence, so a cluster
   with three copies of `PF00165` contributes `PF00165` once.
2. Clusters with **fewer than 3 distinct domains** are not assignable. They
   appear in no family and in neither family artifact.
3. Two assignable clusters are **linked** when the Jaccard similarity of their
   domain sets, `|A ∩ B| / |A ∪ B|`, is **0.50 or greater**.
4. A **family** is a connected component of that link graph — single linkage.
   Two clusters that are not similar to each other belong to the same family if
   a chain of links connects them.
5. A family's **`family_id` is the lexicographically smallest `cluster_id`
   among its members.** Do not number families; derive the name from the
   contents, so the result does not depend on any ordering you chose.
6. A cluster linked to nothing is a family of one, named after itself.

## Required outputs

Seven CSV files in `outputs/`, each with a header row and the exact column names
below. Column order is not graded; column names are. Row order is not graded.

### `bgc_inventory.csv`
One row per detected cluster, over all 120 assemblies.

`cluster_id`, `accession`, `sequence_id`, `start`, `end`, `length_bp`,
`n_distinct_domains`, `type`

`cluster_id` is GECCO's own cluster identifier. `accession` is the RefSeq
assembly accession the cluster came from, as named in `strains.tsv`.
`n_distinct_domains` is the size of the cluster's domain set as defined under
*Gene cluster families* below.

### `contig_edge_clusters.csv`
One row per cluster that sits on a contig boundary, by the rule above.

`cluster_id`, `accession`, `sequence_id`, `contig_length_bp`, `boundary`

`boundary` is `start`, `end`, or `both`.

### `clusters_per_genome.csv`
One row per staged assembly.

`accession`, `n_clusters`, `n_edge_clusters`, `assembly_bp`, `pct_bp_in_clusters`

`assembly_bp` is the total length of all contigs in that assembly.
`n_edge_clusters` counts that assembly's clusters on a contig boundary, by the
rule above. `pct_bp_in_clusters` is the summed cluster length as a percentage
of `assembly_bp`, rounded to 4 decimal places.

### `type_counts.csv`
One row per distinct cluster type observed.

`type`, `n_clusters`, `n_genomes`, `n_species`, `median_length_bp`

`n_genomes` and `n_species` are the distinct assemblies and distinct species in
which that type occurs. `median_length_bp` is the median cluster length for that
type, as an integer; for an even number of clusters take the mean of the two
middle values, rounded down.

### `species_type_counts.csv`
One row per observed species and cluster type combination.

`species`, `type`, `n_clusters`

### `gene_cluster_families.csv`
One row per assignable cluster. Clusters with fewer than 3 distinct domains are
excluded entirely.

`cluster_id`, `family_id`, `family_size`

`family_size` is the number of clusters in that cluster's family.

### `gcf_summary.csv`
One row per family.

`family_id`, `n_clusters`, `n_genomes`, `n_species`

`n_genomes` and `n_species` count the distinct assemblies and distinct species
represented among the family's members.
