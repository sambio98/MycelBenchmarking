# Metabologenomic link ranking across 111 paired *Salinispora* strains

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget six hours.
Most of it is two heavy stages: BGC detection over 111 genomes, and parsing 111
LC-MS/MS files.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- You may write intermediates anywhere else in the working directory. The
  spectra are large; process them one file at a time rather than holding them
  all in memory.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tools named below. Any other resource
  must go through the web access request tool, where it is adjudicated. GNPS,
  MassIVE, the Paired Omics Data Platform and NCBI are denied: the staged files
  are the only inputs this task needs.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  versions, and anything you had to decide for yourself.

## Research objective

Metabologenomics asks which biosynthetic gene cluster makes which molecule, by
looking for gene clusters and metabolites that appear in the same strains and
are absent from the same strains. With enough strains, co-occurrence alone
ranks candidate links well enough to prioritise experiments.

You have 111 *Salinispora* strains with both a sequenced genome and an LC-MS/MS
run of the same isolate. Detect the gene clusters, work out which molecules
each strain produces, and rank the candidate gene-cluster-to-molecule links.

Four properties of the inputs matter, and all four are the data's own:

- **The spectra are strain-resolved; the clustered file is not.** Each strain
  has its own `.mzXML`. The deposit's GNPS-clustered consensus spectra, also
  staged, define the molecular families, but every one of its records says
  `FILENAME=specs_ms.pklbin` — it cannot tell you which strain a consensus
  spectrum came from. Recovering that is the first half of the task.
- **The raw MS2 spectra are profile-mode and mostly noise**, a median of about
  2,700 peaks each. They must be filtered before any similarity is meaningful.
- **Draft and finished assemblies are mixed**, so some gene clusters sit on
  contig boundaries and are truncated.
- **Absence is informative, and asymmetrically so.** A gene cluster present in
  a strain whose metabolome lacks the molecule is evidence against a link; a
  molecule present without the cluster is not equally strong evidence, because
  a detector can miss a cluster. The scoring below encodes that.

### Staged data

| file | what it is |
|---|---|
| `data/genomes/*.fna.gz` | 111 *Salinispora* assemblies from NCBI RefSeq, gzipped nucleotide FASTA, named by assembly accession. |
| `data/spectra/*.mzXML` | One LC-MS/MS run per strain, from MassIVE MSV000078836: the ethyl acetate extract of each isolate, all from the same screening batch and instrument. |
| `data/MSV000078836_specs_ms.mgf` | The deposit's GNPS-clustered consensus MS/MS spectra, 100,263 records. |
| `strains.tsv` | Accession, organism, strain code, assembly level, and the genome and spectrum file for each strain. |

Everything is pinned by sha256 in `data.manifest.json`.

## Method guidance

Use GECCO 0.11.0 for cluster detection, `pyteomics` 5.0.2 to read the mzXML,
`matchms` 0.33.1 for spectral similarity and `networkx` for the graphs.

### The genomic half

Run GECCO on each assembly with its defaults. For every cluster it reports,
take its `type` as the product type, using `Unknown` where GECCO reports none,
and take the **distinct** Pfam accessions in its `domains` field — the field is
a multiset and lists a repeated domain more than once.

Assign **gene cluster families** by single linkage: join two clusters when the
Jaccard similarity of their distinct-domain sets is at least **0.50**. A
cluster with fewer than **3** distinct domains is too thin to place and is left
out of every family rather than forced into one. Name a family for the
**lexicographically smallest cluster identifier** it contains.

### The metabolomic half

**Filter the consensus spectra** from the `.mgf`: discard those with no peaks,
normalise intensities, remove peaks within **17 Da** of the precursor, then
require at least **6** surviving peaks.

**Filter each strain's MS2 scans** the same way, with one extra step. For every
MS2 scan in a strain's mzXML: normalise intensities, remove peaks within 17 Da
of the precursor, keep the **50 most intense** surviving peaks, and require at
least 6. Break ties in that cut by m/z ascending. Ignore MS1 scans.

**Decide presence.** A consensus spectrum is **present** in a strain when that
strain has at least one filtered MS2 scan whose precursor m/z is within
**0.5 Da** of the consensus precursor and whose **modified cosine** against it,
at a fragment tolerance of **0.5 Da**, is at least **0.70**.

**Build the molecular families** over the consensus spectra, exactly as a GNPS
network: modified cosine at 0.5 Da tolerance, compare only pairs whose
precursor m/z differ by at most **2.0 Da**, keep a pair when the score is at
least **0.70** and at least **6** peaks matched, then apply the top-**10** rule
— an edge survives only if it is among the 10 highest-scoring edges of **both**
its nodes, ranking by descending cosine then ascending partner scan. Every
connected component is a molecular family, including a single spectrum on its
own. Name a family `MF` followed by its smallest member scan.

A molecular family is present in a strain when **any** of its spectra is.

### The links

Score every (gene cluster family, molecular family) pair with the
**standardised Metcalf score** over the 111 strains. For each strain, add:

| | molecular family present | absent |
|---|---|---|
| **gene cluster family present** | +10 | −10 |
| **absent** | 0 | +1 |

Score only molecular families present in at least one strain: a family nobody
detected scores identically against every cluster family and says nothing.

Then rank. Sort descending by score, breaking ties by `gcf_id` then `mf_id`
ascending, and keep the top **500** rows. Rank is **dense**: tied scores share
a rank and the next distinct score takes the next integer.

## Required outputs

Seven CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded.

### `strain_inventory.csv`
One row per strain in `strains.tsv`.

`accession`, `species`, `strain`, `n_ms2_scans`, `n_usable_ms2`,
`n_spectra_present`

`species` is the first two words of the organism string. `n_ms2_scans` counts
every MS2 scan in the file; `n_usable_ms2` counts those surviving the filter;
`n_spectra_present` counts the consensus spectra present in that strain.

### `bgc_inventory.csv`
One row per cluster GECCO reports.

`cluster_id`, `accession`, `product_type`, `n_domains`

`n_domains` is the number of **distinct** Pfam accessions.

### `gene_cluster_families.csv`
One row per gene cluster family.

`gcf_id`, `n_clusters`, `n_strains`, `n_species`, `representative_cluster`

### `spectrum_occurrence.csv`
One row per (consensus spectrum, strain) pair where the spectrum is present.
Spectra present in no strain are absent from this file.

`scan`, `accession`

### `molecular_families.csv`
One row per molecular family, including those present in no strain.

`mf_id`, `n_spectra`, `n_strains`, `n_species`

### `metcalf_scores.csv`
One row per scored (gene cluster family, molecular family) pair.

`gcf_id`, `mf_id`, `metcalf_score`, `n_shared_strains`

`n_shared_strains` counts strains where both are present.

### `top_links.csv`
The top 500 rows of `metcalf_scores.csv` by the ranking above.

`rank`, `gcf_id`, `mf_id`, `metcalf_score`, `n_shared_strains`
