# Ketosynthase phylogeny across the sequenced *Salinispora*

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget two hours:
the gene calling is about ten minutes over 120 assemblies and the all-pairs
alignment is the rest.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tools named below. Any other resource
  must go through the web access request tool, where it is adjudicated. NCBI,
  Pfam and InterPro are denied: the staged assemblies and the HMM that ships
  inside the installed `gecco` package are the only inputs this task needs.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  versions, and anything you had to decide for yourself.

There are 120 assemblies and several thousand domains. Expect to script
everything, and do not attempt a web tool for the phylogeny.

## Research objective

The ketosynthase domain is the standard marker for polyketide synthase
phylogeny: every PKS module carries one, it is long enough to align, and its
tree follows the biosynthetic pathway rather than the organism that hosts it.
That last property is the interesting one, and the final deliverable measures
it: a ketosynthase clade that contains domains from several *Salinispora*
species is evidence of a pathway older than the species boundary, while one
confined to a single species is evidence of a recent or lineage-specific
pathway.

Given the staged assemblies, find every ketosynthase domain, build a distance
tree over them, and characterise its clades taxonomically.

Two properties of the input matter:

- **These are draft assemblies as well as finished ones.** A gene at a contig
  boundary is truncated, which yields a partial domain. The method below sets a
  minimum domain length for exactly this reason.
- **Many ketosynthase domains sit in the same protein.** A modular PKS is one
  very long polypeptide carrying several modules, so reporting one domain per
  protein loses most of the data.

### Staged data

| file | what it is |
|---|---|
| `data/genomes/*.fna.gz` | 120 *Salinispora* genome assemblies from NCBI RefSeq, gzipped nucleotide FASTA, named by assembly accession. |
| `strains.tsv` | Accession, organism string, strain and assembly level for each. |

Every assembly is pinned by sha256 in `data.manifest.json`.

## Method guidance

Use `pyrodigal` 3.7.1 for gene calling, `pyhmmer` 0.12.3 for the domain search,
Biopython for the alignment scores and `scipy` for the tree. The ketosynthase
HMM does not need downloading: the installed `gecco` package ships a subset of
Pfam 35.0 at `gecco/hmmer/Pfam.h3m`, and the model you want is in it.

**Call the genes.** One `pyrodigal.GeneFinder` per assembly, **trained on that
whole assembly** — not in `meta` mode, and not trained per contig. Translate
every gene. Name proteins so that an identifier is stable across runs.

**Find the domains.** Search the single Pfam model whose accession is
**PF00109** (`ketoacyl-synt`) over each assembly's proteins at **E = 1e-5**.
Keep only hits and domains that HMMER itself marks as *included*. For each
included domain, take the **aligned region on the target protein** — the
alignment's target coordinates, 1-based and inclusive — not the envelope.
Discard any domain whose aligned region is shorter than **150 amino acids**.

Within one assembly, number the surviving domains `KS001`, `KS002`, … in order
of protein identifier then start coordinate, and give each the identifier
`<accession>_KS<nnn>`.

**Derive the species** from `strains.tsv` as the first two whitespace-separated
words of the organism string, so the strain suffix is dropped.

**Collapse identical sequences.** Group the domains by their exact amino acid
sequence, and name each group by the **lexicographically smallest domain
identifier** in it. These groups, not the individual domains, are the leaves of
the tree. This is not an optimisation: the same ketosynthase recurs verbatim
across strains of one species, so a tree built over every domain would be
mostly zero-length ties whose arrangement is decided by input order rather than
by the sequences.

**Compute distances.** Use a Biopython `PairwiseAligner` in **local** mode with
the **BLOSUM62** matrix, gap open **-11** and gap extend **-1**. Use alignment
**scores only, never the traceback**. For two domains a and b:

```
distance = 1 - score(a, b) / min(score(a, a), score(b, b))
```

clamped to the range 0 to 1. Normalising by the smaller self-score, rather than
by length, is what stops a short domain looking distant from the longer one it
sits inside.

**Build the tree.** Average linkage over the distance matrix of the groups —
this is UPGMA.
`scipy.cluster.hierarchy.linkage(condensed, method="average")` is the intended
call.

**Name the splits.** Every merge in the tree is one internal node, and so one
split: the set of sequence groups below it. Name a split by the
**lexicographically smallest group identifier below it**, and record how many leaves it has and
the height at which it was merged. A split's name is not unique on its own — a
split and its parent can share a smallest leaf — so a split is identified by
its name together with its leaf count.

## Required outputs

Eight CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded. Round every distance and height to 4 decimal places, and bit scores to 1.

### `ks_domains.csv`
One row per surviving ketosynthase domain.

`domain_id`, `accession`, `protein_id`, `start`, `end`, `length_aa`,
`bitscore`

`start` and `end` are on the protein, 1-based inclusive; `length_aa` is
`end - start + 1`; `bitscore` is the domain's score, not the hit's.

### `genome_ks_counts.csv`
One row per staged assembly, including any with no surviving domain.

`accession`, `species`, `n_proteins`, `n_ks_proteins`, `n_ks_domains`,
`max_domains_per_protein`

`n_proteins` is the whole gene call. `max_domains_per_protein` is `0` for an
assembly with no domain.

### `multidomain_proteins.csv`
One row per protein carrying **two or more** surviving domains.

`accession`, `protein_id`, `n_domains`, `first_domain_id`, `span_aa`

`first_domain_id` is the smallest domain identifier in that protein. `span_aa`
runs from the first domain's start to the last domain's end, inclusive.

### `sequence_groups.csv`
One row per distinct domain amino acid sequence.

`group_id`, `length_aa`, `n_domains`, `n_genomes`, `n_species`

`n_domains` is how many domains carry that exact sequence; `n_genomes` and
`n_species` count the assemblies and species those domains come from.

### `species_domain_counts.csv`
One row per species present in `strains.tsv` that has at least one domain.

`species`, `n_genomes`, `n_domains`, `mean_domains_per_genome`

`n_genomes` counts that species' staged assemblies, including any with no
domain, so the mean is over all of them.

### `group_nearest_neighbour.csv`
One row per sequence group: its closest other group.

`group_id`, `nearest_group_id`, `distance`

Ties are broken by the partner's identifier ascending.

### `tree_splits.csv`
One row per internal node of the tree.

`split_id`, `n_leaves`, `height`

### `split_taxonomy.csv`
One row per internal node, with the taxonomy of everything below it.

`split_id`, `n_leaves`, `n_domains`, `n_genomes`, `n_species`, `breadth`

`n_leaves` counts the sequence groups below the split; `n_domains` counts the
domains those groups stand for. `breadth` is `single_species` when every one of
those domains comes from one species, and `multi_species` otherwise.
