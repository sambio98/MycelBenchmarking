# Task 07 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

## Why mass search and not the tools this category names

PROJECT.md §4 category 4 names SIRIUS, CSI:FingerID and NP Atlas. §9 then rules
the first two out of any graded artifact: both call a web service whose model
version determines the output, and SIRIUS 4, 5.8.x and 6.x disagree with each
other. A GNPS library hit list is out for the same reason — the library is
continuously recurated, so a 2021 hit list is not reproducible today.

What §9 names as safe is matching against a database snapshot frozen inside the
bundle. So this task is an NP Atlas release pinned by sha256, and every number
in `truth/` is arithmetic over it. No network call happens at grading time or
at run time.

The cost is honest and is the point of the task: **mass alone cannot identify a
compound.** The artifacts are built to measure whether an agent reports that
ambiguity or collapses it into one confident answer.

## Pinning an unversioned download

NP Atlas serves its current release from `NPAtlas_download.json`, with no
version in the URL, so a new release replaces it in place and the release
cannot be recovered from the file. Two defences:

- The manifest pins bytes and sha256, so drift fails `npbench fetch` loudly
  rather than silently producing different truth.
- `solution/make_truth.py` asserts the record count before doing any work.
  v2024_09 has exactly 36,454 records.

Measured, and worth correcting against the estimate in PROJECT.md §9: the JSON
export is **475 MB and 36,454 records**, where §9 says 676 MB and 36,545
compounds. The licence note in §9 stands — NP Atlas is CC BY-NC 4.0, and
whether that is compatible with the benchmark's own licence is still an open
question for the release, not something this task can settle.

## Scale

| | |
|---|---|
| NP Atlas release | v2024_09, 36,454 records, 475 MB |
| after the InChIKey dedupe | 36,422 compounds |
| duplicate InChIKeys removed | 32 |
| compounds sharing a stereo-stripped skeleton | 2,925 |
| compounds with no npclassifier block | 3,923, reported as `Unclassified` |
| origin genera | 1,246 (*Streptomyces* 5,758 compounds, *Salinispora* 126) |
| spectra retained from the mgf | 4,081 of 100,263 |
| eligible queries | 3,811 — charge 2 excluded |
| candidate annotations | 10,953 |
| queries with at least one candidate | 2,058 of 3,811, so 46% stay unknown |
| distinct formulae matched | 1,538 |
| molecular families | 325, from 3,208 edges |
| artifacts / metric instances | 8 / 35 |
| staged | 493 MB |
| runtime | ~3 min |

That 46% of queries get no candidate at all is the normal and correct outcome
of mass-only dereplication against a curated database, not a failure.

## Reproducibility

Two independent runs produce byte-identical output for all eight artifacts.

## Measured difficulty

`scripts/difficulty_task07.py`, 21 variants spanning three layers — how the
reference release is read, how the query set is chosen, and how the search and
aggregates are computed. Same cheap-proxy caveat as the other tasks.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: dedupe on InChIKey block | 1.00 |
| one slip: no dedupe at all | 1.00 |
| one slip: release's own adduct m/z | 1.00 |
| one slip: searches charge 2 too | 1.00 |
| one slip: formula from the SMILES | 1.00 |
| one slip: window from neutral mass | 1.00 |
| one slip: drops the unclassified | 0.88 |
| one slip: average H, not a proton | 0.88 |
| one slip: splits pathway labels | 0.88 |
| one slip: families numbered | 0.88 |
| one slip: best match on signed ppm | 0.75 |
| one slip: first pathway only | 0.50 |
| one slip: 20 ppm | 0.38 |
| one slip: 5 ppm | 0.38 |
| one slip: average mass, not exact | 0.25 |
| one slip: keeps the empty spectra | 0.25 |
| one slip: forgets `[M+Na]+` | 0.25 |
| one slip: charge 1 only | 0.12 |
| hasty agent (3 slips) | 0.12 |
| careless agent (5 slips) | 0.25 |

| artifact | pass rate | |
|---|---|---|
| `best_annotation.csv` | 0.52 | in band |
| `pathway_counts.csv` | 0.52 | in band |
| `component_annotations.csv` | 0.52 | in band |
| `annotations.csv` | 0.62 | in band |
| `formula_candidates.csv` | 0.62 | in band |
| `genus_annotations.csv` | 0.62 | in band |
| `queries.csv` | 0.86 | **above band** |
| `atlas_compounds.csv` | 0.95 | **above band** |

Six of eight land in band and the deep artifacts all discriminate. The single
most expensive mistake is using the average molecular weight instead of the
monoisotopic mass, which takes five artifacts to 0.05 or below while leaving
the two shallow ones untouched — a good shape, because it means the task
distinguishes reading the inputs from doing the analysis.

### The 5-slip agent scores *above* the 3-slip agent

0.25 against 0.12, and this is not noise. The headline is the fraction of
artifacts that clear 0.80, so what matters is how many artifacts a mistake set
touches, not how many mistakes it contains. The careless set wrecks six
artifacts and leaves `atlas_compounds.csv` and `queries.csv` passing at 0.96
and 0.97; the hasty set includes `charge_1_only`, which drags `queries.csv` to
0.70 and fails it, so only one artifact survives.

The practical consequence is a caution about the 3-slip and 5-slip columns in
every one of these notes, this one included: **they are not monotone in the
number of mistakes, and they are not comparable across tasks** unless the
mistake sets happen to cover the dependency graph the same way. The
per-artifact pass rate is the figure to trust.

### A large reference table is a forgiving artifact

`atlas_compounds.csv` at 0.95 is the clearest instance of something worth
stating generally: it has 36,422 rows, and identifier F1 over 36,422 rows
cannot fail on a mistake that changes a few thousand of them. Deduplicating on
the InChIKey's first block discards **2,925 distinct stereoisomers** — a
serious error, and exactly the one the method warns against — and still scores
0.96, because 33,497 of the rows are untouched.

So the artifact is not easy in the sense of being obvious; it is forgiving in
the arithmetic sense. Three mistakes that a reviewer would call serious all
pass it. The only one that fails is `first_pathway_only` (0.24), which changes
a label on most rows rather than removing a minority of them.

This is the counterpart to task 01's finding. There, an artifact looked easy
because the variant set never touched its layer. Here the variants do touch it,
and the artifact absorbs them because row-set F1 over tens of thousands of rows
is a blunt instrument. **A reference-table artifact should be graded on a
column that moves with the mistake, not on its row set** — which is why the
pathway label, the genus and the formula are all graded here, and why the three
mistakes that change a per-row label are the ones that register.

`queries.csv` at 0.86 is the simpler case: it depends only on the spectrum
filter and the charge rule, and of the 21 variants only three touch either.

### Five variants that cost nothing, and why

Measured, not assumed. Three are properties of the data worth knowing:

- **Deduplicating the release is nearly a no-op.** Only 32 of 36,454 records
  are a second name for a structure already present, so skipping the dedupe
  entirely scores 1.00 everywhere — identifier F1 of 36,422/36,454 is 0.9996.
  The *wrong* dedupe (on the InChIKey's first block) is punished, at 0.96, so
  the rule is exercised in one direction only.
- **NP Atlas's stated formulae are internally consistent.** Recomputing every
  formula from the record's own SMILES with RDKit reproduces the `mol_formula`
  field for all 36,422 compounds, so that variant is a complete no-op. Worth
  contrasting with task 03, which found **56** formula disagreements in MIBiG's
  compound set. The two databases are not equally clean, and a task built on
  NP Atlas cannot borrow task 03's formula-disagreement artifact.
- **The release's precomputed adduct masses agree with computing them.**
  `m_plus_h` and `m_plus_na` match `exact_mass` plus the shifts to four
  decimals, so using them instead scores 1.00 on seven of eight artifacts.
- **Setting the 10 ppm window from the neutral mass rather than the measured
  m/z** changes the window by about 0.25% of its width, which moves four
  annotations out of 10,953.

The fifth was not a property of the data but a flaw in an artifact, and it was
fixed rather than recorded as difficulty. `formula_candidates.csv` originally
carried a fourth column, `n_compounds_matched`, and no mistake could get it
wrong. The reason is arithmetic: **a molecular formula fixes the monoisotopic
mass**, so every compound sharing a matched formula sits at the same mass and
falls inside the same window. It was provably equal to `n_compounds`, confirmed
on all 1,538 rows, and a column that cannot differ from its neighbour is not a
test. It is gone, and the variant that got it wrong is gone with it.

### Across the eight tasks

See the README for the current table. Task 07 behaves like 01, 04 and 05: the
artifacts share dependencies and the deep ones all discriminate.
