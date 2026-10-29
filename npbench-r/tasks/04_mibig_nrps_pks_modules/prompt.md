# NRPS and PKS module architecture and adenylation-domain substrate specificity in MIBiG 4.0

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. The analysis is
cheap; budget an hour, almost none of it compute.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tool named below. Any other resource
  must go through the web access request tool, where it is adjudicated. The
  MIBiG website is denied: the staged release is the only input this task needs.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  versions, and anything you had to decide for yourself.

## Research objective

Modular non-ribosomal peptide synthetases and polyketide synthases build their
products one module at a time, and the substrate each adenylation domain
selects largely determines the peptide's sequence. MIBiG curates this
architecture per cluster: an ordered list of modules, each with its domain
complement, and for adenylation domains the substrates they are believed to
activate together with the evidence behind each belief.

Given the frozen MIBiG 4.0 release, reconstruct that architecture and
characterise the substrate repertoire: what the modules are, which substrates
are claimed with what strength of evidence, which substrates are not
proteinogenic, and which clusters share an identical module architecture.

Five properties of the input matter, and all five are MIBiG's own data.

First, **module names repeat within a cluster.** Twelve such pairs exist, so a
module is identified by its **1-based position in the `modules` list**, not by
its name. Preserve that order; several outputs depend on it.

Second, **over half the modules are marked `active: false`**. They are curated
modules and all of them are in scope. Do not filter on `active`; report it.

Third, **almost every adenylation domain's `location` is the sentinel
`{"from": -1, "to": -1}`.** Coordinates are not usable here and are not asked
for. Do not filter modules on them.

Fourth, **an adenylation domain carries between zero and nine substrates**, and
295 carry no evidence entries at all. The substrate table is not one row per
module.

Fifth, **nine substrate names map to more than one distinct structure**,
`alanine` and `valine` among them, because MIBiG records D and L forms under
one name. One output asks for exactly that count.

Sixth, **some substrates are spelled two ways**, differing only in
capitalisation — `Leucine` and `leucine`, `Serine` and `serine` — so the
verbatim names overcount the substrates. Fold case before counting them.

### Staged data

| Path | Contents |
|---|---|
| `data/mibig_json_4.0.tar.gz` | the MIBiG 4.0 JSON release, one file per BGC |
| `data.manifest.json` | the provenance record for the staged file |

Each record carries an `accession`, a `biosynthesis.classes` list whose entries
have a `class`, and a `biosynthesis.modules` list. Only records with a
non-empty `modules` list are in scope. Each module has `name`, `type`,
`active`, `genes`, `carriers`, optionally `modification_domains`, and
optionally the domain slots `a_domain`, `c_domain`, `ks_domain`, `at_domain`.
An `a_domain` may carry `substrates` (each with `name`, `structure` and
`proteinogenic`) and `evidence` (each with a `method`).

## Method guidance

Use **RDKit 2026.9.1** (the `rdkit` wheel on PyPI; `rdkit.__version__` prints
it as `2026.09.1`) for the substrate descriptors. Pin the version: formula and
exact mass are version-sensitive.

Every rule below is part of the method; none of it is open to choice.

**Scope.** Only BGCs whose `biosynthesis.modules` list is non-empty. Within
those, every module, regardless of `active`.

**Module identity.** `module_index` is the module's 1-based position in the
`modules` list as it appears in the JSON.

**A module's domains** are the ones present in the slots `a_domain`,
`c_domain`, `ks_domain`, `at_domain` — taken **in that order** — followed by
every entry of `modification_domains` in list order. Each contributes its
`type` string. `n_domains` is how many that is.

**Substrate names.** Compare names case-insensitively after collapsing
whitespace, and report the **lexicographically smallest** spelling among a
substrate's variants. Use that spelling consistently in every output, including
`a_domain_substrates.csv`. Do not lowercase names: capitalisation carries
meaning in `N-methyl-L-phenylalanine`.

**Substrate scope.** Only substrates with a non-empty `name`.
`substrate_formula` is `rdMolDescriptors.CalcMolFormula` of the parsed
`structure`, verbatim, or the **empty string** where `structure` is missing or
RDKit cannot parse it. Two substrate entries in the release fall in that case.

**Best evidence.** An adenylation domain's `best_evidence` is the strongest
method among its `evidence` entries, by this ranking, strongest first:

```
ATP-PPi exchange assay  >  ACVS assay  >  Enzyme-coupled assay  >
Activity assay  >  Heterologous expression  >  Feeding study  >
Mass spectrometry  >  Structure-based inference  >  Homology  >
Sequence-based prediction
```

In vitro assays that measure the domain directly outrank in vivo work, which
outranks observing the product, which outranks inferring from the compound
structure, which outranks homology and bare sequence prediction. An
adenylation domain with no evidence entries has `best_evidence` of `none`. A
method outside this list is reported verbatim; none occur in this release.

**`proteinogenic` per substrate name.** MIBiG is not self-consistent about
this flag across entries sharing a name. Report the value carried by the
majority of that name's entries, with `true` winning a tie.

**Architecture.** A cluster's `architecture` is its module `type` values in
list order joined by `>`. Its `domain_composition` is every domain type from
every one of its modules, **sorted alphabetically** and joined by `;`,
keeping duplicates.

**Class cross-tabulation.** A cluster may carry several biosynthesis classes.
A substrate assignment counts once for each distinct class of its cluster.

## Required outputs

Eight CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded. Booleans are the lowercase strings `true` and `false`.

### `modules.csv`
One row per module in scope.

`bgc_accession`, `module_index`, `module_name`, `module_type`, `active`,
`n_genes`, `has_a_domain`, `n_modification_domains`, `n_domains`

### `a_domain_substrates.csv`
One row per (module, substrate name) where the module has an adenylation domain.

`bgc_accession`, `module_index`, `substrate_name`, `proteinogenic`,
`substrate_formula`, `best_evidence`

### `substrate_counts.csv`
One row per distinct substrate name.

`substrate_name`, `n_modules`, `n_bgcs`, `proteinogenic`,
`n_distinct_structures`

`n_distinct_structures` counts the distinct canonical SMILES recorded under
that name, ignoring entries RDKit cannot parse.

### `substrate_evidence_summary.csv`
One row per distinct `best_evidence` value.

`best_evidence`, `n_assignments`, `n_distinct_substrates`, `n_bgcs`

### `nonproteinogenic_substrates.csv`
One row per (substrate name, canonical SMILES) pair whose `proteinogenic` flag
on that assignment is false and whose structure RDKit can parse.

`substrate_name`, `canonical_smiles`, `formula`, `exact_mass`,
`n_heavy_atoms`, `n_modules`

`exact_mass` is `Descriptors.ExactMolWt` rounded to 4 decimal places.

### `class_substrate_matrix.csv`
One row per (biosynthesis class, substrate name) pair.

`biosynthesis_class`, `substrate_name`, `n_modules`, `n_bgcs`

### `module_architecture.csv`
One row per BGC in scope.

`bgc_accession`, `n_modules`, `n_a_domains`, `architecture`,
`domain_composition`

### `architecture_counts.csv`
One row per distinct `architecture` string.

`architecture`, `n_bgcs`, `n_modules`, `n_a_domains`, `representative_bgc`

`n_modules` is the module count of that architecture, which is the same for
every cluster sharing it. `n_a_domains` is summed over those clusters.
`representative_bgc` is the lexicographically smallest accession among them.
