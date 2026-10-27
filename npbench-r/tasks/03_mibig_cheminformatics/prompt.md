# Scaffold, isomer and similarity structure of the MIBiG 4.0 compound set

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

There are several thousand structures. Expect to script everything.

## Research objective

MIBiG is the reference repository of experimentally characterised biosynthetic
gene clusters, and each entry may report the compounds its cluster produces as
SMILES strings. Given the frozen MIBiG 4.0 release, characterise the chemical
space it covers: the per-compound descriptors, the Bemis-Murcko scaffolds, the
molecules that cannot be told apart by molecular formula alone, the ones that
cannot be told apart once stereochemistry is dropped, and the similarity
families the set falls into.

Four properties of the input matter for the result, and all four are MIBiG's
own data rather than anything added for this task.

First, **MIBiG's stated `formula` and `mass` fields disagree with the
structure** for a few hundred compound entries. Every descriptor asked for
below is computed from the SMILES. Do not copy either field.

Second, **a BGC can list the same compound name twice** with two different
structures. Six such pairs exist in the release — four identical names and two
differing only in capitalisation — so `(bgc_accession, compound_name)` is not
unique as given and needs the tie-break specified below.

Third, **the same molecule is reported by many different BGCs**. Artifacts keyed
on a structure are keyed on the distinct canonical SMILES, not on the MIBiG
entry.

Fourth, **hundreds of distinct structures share a molecular formula, and
hundreds collapse onto each other once stereochemistry is removed.** Both are
asked for explicitly; neither is an error to clean up.

### Staged data

| Path | Contents |
|---|---|
| `data/mibig_json_4.0.tar.gz` | the MIBiG 4.0 JSON release, one file per BGC |
| `data.manifest.json` | the provenance record for the staged file |

Each BGC record carries an `accession`, a `compounds` list whose entries may
have `name`, `structure` (SMILES), `formula` and `mass`, and a
`biosynthesis.classes` list whose entries carry a `class`.

## Method guidance

Use **RDKit 2026.9.1** (the `rdkit` wheel on PyPI; `rdkit.__version__` prints
it as `2026.09.1`). Pin the version: canonical
SMILES output, scaffold perception and fingerprint bits are all
version-sensitive, and this task is graded against that build.

Every rule below is part of the method; none of it is open to choice.

**Selection.** Consider only compound entries that have both a non-empty
`name` and a non-empty `structure`, and whose structure RDKit can parse with
`Chem.MolFromSmiles` without sanitisation errors. Discard the rest silently.

**Canonical form.** A structure's canonical SMILES is
`Chem.MolToSmiles(mol)` with RDKit's defaults, which keeps stereochemistry.

**Duplicate names.** Where one BGC lists the same `name` more than once, keep
the single entry whose canonical SMILES sorts lexicographically smallest, and
drop the others. Compare names **case-insensitively** after trimming
whitespace: two of the six collisions in the release differ only in
capitalisation.

**Descriptors**, all computed from the parsed molecule:
`formula` is `rdMolDescriptors.CalcMolFormula` verbatim, including any charge
suffix; `exact_mass` is `Descriptors.ExactMolWt` rounded to 4 decimal places;
`n_heavy_atoms` is `GetNumHeavyAtoms()`; `n_rings` is
`rdMolDescriptors.CalcNumRings`.

**Scaffold.** `MurckoScaffold.GetScaffoldForMol`, reported as its canonical
SMILES. An acyclic molecule yields an empty scaffold; keep it as an empty
string with zero rings and zero atoms rather than dropping the compound.

**Flat form.** The stereochemistry-free form is the canonical SMILES of the
molecule re-parsed from `Chem.MolToSmiles(mol, isomericSmiles=False)`.

**Similarity families.** Over the distinct structures only:

1. Fingerprint each with `rdFingerprintGenerator.GetMorganGenerator(radius=2,
   fpSize=2048)`.
2. Two structures are **linked** when their Tanimoto similarity is **0.50 or
   greater**.
3. A **family** is a connected component of that link graph — single linkage. A
   structure linked to nothing is a family of one.
4. A family's **`family_id` is the lexicographically smallest canonical SMILES
   among its members.** Do not number families; derive the name from the
   contents.

**Class cross-tabulation.** A BGC may carry several biosynthesis classes, and a
structure may be reported by several BGCs. A structure counts once for each
distinct class across all BGCs that report it.

## Required outputs

Eight CSV files in `outputs/`, each with a header row and the exact column
names below. Column order is not graded; column names are. Row order is not
graded. **SMILES are compared case-sensitively**, because case encodes
aromaticity; names and formulas are not.

### `compounds.csv`
One row per retained `(bgc_accession, compound_name)` pair.

`bgc_accession`, `compound_name`, `canonical_smiles`, `formula`, `exact_mass`,
`n_heavy_atoms`, `n_rings`

`compound_name` is MIBiG's `name` verbatim.

### `formula_disagreements.csv`
One row per retained compound whose MIBiG `formula` field is non-empty and
differs from the computed formula.

`bgc_accession`, `compound_name`, `stated_formula`, `computed_formula`

### `scaffolds.csv`
One row per distinct canonical structure.

`canonical_smiles`, `scaffold_smiles`, `n_scaffold_rings`, `n_scaffold_atoms`

### `scaffold_counts.csv`
One row per distinct scaffold.

`scaffold_smiles`, `n_structures`, `n_bgcs`, `n_classes`, `representative_smiles`

Counts are over distinct structures carrying that scaffold, the distinct BGCs
reporting any of them, and the distinct biosynthesis classes of those BGCs.
`representative_smiles` is the lexicographically smallest canonical SMILES
among the structures carrying that scaffold.

### `formula_collisions.csv`
One row per molecular formula carried by **two or more** distinct structures.

`formula`, `n_structures`, `n_bgcs`, `representative_smiles`

`representative_smiles` is the lexicographically smallest canonical SMILES
among the structures sharing that formula.

### `stereo_collapse_groups.csv`
One row per (flat form, structure) pair, for flat forms reached by **two or
more** distinct structures.

`flat_smiles`, `canonical_smiles`, `group_size`

`group_size` is the number of distinct structures sharing that flat form.

### `similarity_families.csv`
One row per distinct canonical structure.

`canonical_smiles`, `family_id`, `family_size`

### `class_scaffold_summary.csv`
One row per observed biosynthesis class and scaffold ring count.

`biosynthesis_class`, `n_scaffold_rings`, `n_structures`, `n_bgcs`,
`representative_smiles`

`representative_smiles` is the lexicographically smallest canonical SMILES
among that cell's structures.
