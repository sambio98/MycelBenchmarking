# Exact-mass dereplication of a marine actinomycete molecular network against the NP Atlas

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget an hour.
Most of the compute is one all-pairs spectral comparison; the dereplication
itself is arithmetic.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tools named below. Any other resource
  must go through the web access request tool, where it is adjudicated. The NP
  Atlas website, GNPS and MassIVE are denied: the two staged files are the only
  inputs this task needs.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  versions, and anything you had to decide for yourself.

The reference file is 475 MB of JSON holding 36,454 records. Do not try to read
it by hand, and do not load it more than once.

## Research objective

Dereplication is the step that asks, of every feature in a metabolomics
dataset, whether it is already a known compound. Here you will dereplicate the
MS/MS features of a marine actinomycete survey against a frozen release of the
NP Atlas, the curated database of microbial natural products, and then ask
which parts of the molecular network remain unexplained.

Mass alone cannot identify a compound, and the deliverables are built so that
you report that honestly rather than collapse it. Three properties of the
inputs matter, and all three are the data's own:

- **96% of the spectrum file carries no fragment peaks.** Of its 100,263
  spectra only a small minority have any peaks at all. A pipeline that does not
  notice will search noise.
- **The precursor charge is often not stated.** Many spectra carry `CHARGE=0+`,
  which means the file does not say. Others state `2+`, where a
  singly-charged assumption would give the wrong neutral mass.
- **The reference has duplicate and near-duplicate entries.** A few records are
  a second name for a structure already present, and thousands of compounds
  differ from another only in stereochemistry, which an InChIKey's first block
  cannot see. Several hundred molecular formulae are shared by more than one
  compound.

### Staged data

| file | what it is |
|---|---|
| `data/MSV000078836_specs_ms.mgf` | GNPS-clustered consensus MS/MS spectra from MassIVE MSV000078836, the Crüsemann 146-strain marine actinomycete survey. 100,263 spectra, qTOF CID. |
| `data/NPAtlas_download.json` | The NP Atlas release 2024_09, full JSON export, 36,454 records. Each carries `npaid`, `mol_formula`, `exact_mass`, `inchikey`, `smiles`, an `origin_organism` block and an `npclassifier` block. |

Both are pinned by sha256 in `data.manifest.json`.

## Method guidance

Use `matchms` 0.33.1 for the spectra and `networkx` for the graph. The
dereplication needs no cheminformatics toolkit: it is a mass search.

**Build the reference table.** Deduplicate the release on the **full
InChIKey**, keeping the record with the lowest `id`. Do not deduplicate on the
InChIKey's first block: that block ignores stereochemistry and would discard
distinct stereoisomers. Read `exact_mass` as a number and round it to 4
decimals. Take `origin_type` and `genus` from `origin_organism`. For the
pathway label, take `npclassifier.pathway_results`, sort it, remove duplicates
and join it with `;`; where the block is absent or the list is empty, the label
is the literal string `Unclassified` — do not drop those compounds.

**Filter the spectra.** Discard spectra with no peaks. Normalise intensities.
Remove peaks within **17 Da** of the precursor m/z. Then require at least
**6** surviving peaks. Keep the surviving spectra sorted by scan number.

**Choose the query set.** Keep only spectra whose `CHARGE` is `0` or `1`,
treating both as singly charged. Exclude charge 2 and above rather than
guessing.

**Compute neutral masses.** For each query, from the precursor m/z:

- `[M+H]+`: subtract a proton, **1.007276**
- `[M+Na]+`: subtract **22.989218** (sodium minus one electron)

**Search.** For each query and each of the two adducts, report every
deduplicated compound whose `exact_mass` falls inside a window of **10 ppm of
the measured precursor m/z** — that is, `|neutral_mass - exact_mass| <=
precursor_mz * 10e-6`. Note the window is set by the measured m/z, not by the
candidate's mass. Report the error relative to the candidate:
`ppm_error = (neutral_mass - exact_mass) / exact_mass * 1e6`, rounded to 2
decimals, keeping its sign.

**Pick a best match per query.** Rank that query's candidates by **absolute**
ppm error ascending, then by NPAID ascending, and keep the first. The second
key matters: masses are rounded, so exact ties are common.

**Rebuild the molecular network**, so the per-family deliverable can be
computed. Use modified cosine with a fragment tolerance of **0.5 Da**, compare
only pairs whose precursor m/z differ by at most **2.0 Da**, and keep a pair as
an edge when the score is at least **0.70** and at least **6** peaks matched.
Then apply the top-**10** rule: an edge survives only if it is among the 10
highest-scoring edges of **both** of its nodes. Rank by descending cosine, then
ascending partner scan. A connected component with two or more nodes is a
molecular family; name it `C` followed by its smallest member scan.

Note the network is built from **all** retained spectra, including those
excluded from the query set by charge.

## Required outputs

Eight CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded.

### `atlas_compounds.csv`
One row per compound surviving the InChIKey dedupe.

`npaid`, `mol_formula`, `exact_mass`, `origin_type`, `genus`, `pathway`

### `queries.csv`
One row per spectrum eligible for mass search.

`scan`, `precursor_mz`, `charge`, `neutral_mass_h`, `neutral_mass_na`

`charge` is the integer from `CHARGE`, so `0` where the file states none.

### `annotations.csv`
One row per candidate: a (query, compound, adduct) triple inside the tolerance.
A compound matching one query under both adducts gives two rows.

`scan`, `npaid`, `adduct`, `ppm_error`, `mol_formula`

`adduct` is the literal string `[M+H]+` or `[M+Na]+`.

### `best_annotation.csv`
One row per query that has at least one candidate. Queries with none are
absent.

`scan`, `npaid`, `adduct`, `ppm_error`, `mol_formula`, `genus`, `pathway`

### `formula_candidates.csv`
One row per distinct molecular formula appearing in `annotations.csv`.

`mol_formula`, `n_compounds`, `n_scans`

`n_compounds` counts the compounds carrying that formula **in the whole
deduplicated release**, not only the matched ones. `n_scans` counts the
distinct queries the formula was offered to.

### `genus_annotations.csv`
One row per origin genus with at least one candidate.

`genus`, `n_compounds_matched`, `n_scans`

### `pathway_counts.csv`
One row per distinct pathway label appearing among the matched compounds. The
label is the combined string, so a compound with two pathways belongs to one
row, not two.

`pathway`, `n_compounds_matched`, `n_scans`

### `component_annotations.csv`
One row per molecular family of two or more nodes.

`component_id`, `n_nodes`, `n_annotated`, `n_distinct_formulas`, `top_pathway`

`n_annotated` counts the family's nodes that appear in `best_annotation.csv`,
and `0` is a normal answer. `n_distinct_formulas` counts the distinct formulae
offered to any of its nodes in `annotations.csv`. `top_pathway` is the most
common pathway label among the family's annotated nodes, ties broken
alphabetically; for a family with no annotated node it is `Unclassified`.
