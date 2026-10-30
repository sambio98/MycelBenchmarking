# Task 09 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

**Read the scope section before using this task.** It stages 33.5 GB, which is
1.7x the 20 GB figure PROJECT.md §10 calls ideal, and that is the one decision
here most worth a second opinion.

## How this task became possible

PROJECT.md §4 category 5 — paired omics and metabologenomic link ranking — was
the last category with no task, and §11 calls it the differentiator worth
building. The README recorded it as unbuildable, with three dead ends:

- **The clustered spectra carry no per-strain provenance.** All 100,263 records
  in the staged `.mgf` say `FILENAME=specs_ms.pklbin`, so the file cannot say
  which strain a consensus spectrum came from.
- **MassIVE's own endpoints are unreachable from a current certificate store.**
  Its chain roots in an emSign CA present in neither the system trust store nor
  `certifi`. The GNPS2 proxy that works instead serves only exact paths and
  returns no directory listing, and the GNPS cluster-info tables — which would
  carry the per-file mapping — 404 at every documented name.
- **MSV000098188**, the other paired-omics candidate in §5, has peak files that
  404 at every layout tried, and its 101 "microscale" strains have zero overlap
  with the sequenced set.

All three are real and anyone rebuilding this will hit them. The conclusion
drawn from them was still wrong, because it was one source short: the **Paired
Omics Data Platform**, which §5 names for this very dataset, curates exactly
the missing mapping. Its project record for MSV000078836 holds **360
genome-to-file links** — 120 strains across three extraction methods — each
naming the deposit path of one strain's `.mzXML`. 111 of those strains already
had an NCBI assembly staged for task 02.

The generalisable lesson: *"unavailable" was a claim about the sources I had
tried, not about the data.* Every dead end above concerns file discovery on one
archive, while the mapping lived in a second archive built for precisely this
purpose.

## Scope, and what it costs

**One extraction method, every paired strain.** PoDP records three extraction
methods (ethyl acetate, butanol, methanol) for all 120 strains, all from one
`Screening` batch with identical sample preparation and instrumentation. The
task uses the ethyl acetate extracts only, which holds the analytical batch
constant — comparing metabolomes across extraction methods would confound every
occurrence call downstream. That rule, not an arbitrary cut, is what fixes the
strain count at 111.

**33.5 GB staged**, the largest fetch in the suite by a factor of sixty-seven.
The mzXML files run 113 to 790 MB, mean 300 MB. PROJECT.md §10 says to keep
every task well under 100 GB and ideally under 20 GB, so this satisfies the
limit and is **1.7× the ideal**.

Both of my earlier estimates were wrong and in the same direction, which is
worth recording as a method note: 17.8 GB came from extrapolating a single
160 MB file, and 27.9 GB from the first eight, whose mean was 252 MB. The true
mean is 300 MB and the distribution has a long right tail — the largest file is
790 MB, seven times the smallest. **Extrapolating a staging budget from a
sample of one, or of eight, is not estimating; the only honest number came
from fingerprinting all 111.**

I kept all 111 strains rather than cutting to fit, because "every paired strain
at one extraction method" is a rule and any smaller set needs an arbitrary one.
Staying under 20 GB would mean **66 strains**, dropping 45 of them. Metcalf
scoring degrades gracefully with fewer strains, so that reduction is available
to a reviewer who wants it; it is a curation decision rather than a technical
one, and it is the decision in this suite I would most like a second opinion
on.

**The molecular families are the depositors', not ours.** They come from the
GNPS clustering frozen in the `.mgf` rather than from re-clustering the raw
scans. §9 names network topology from a fixed `.mgf` as the reproducible
artifact, and reimplementing MS/MS clustering deterministically is a research
project rather than a task. What this task recovers is which *strain* each
consensus spectrum belongs to, which is the part the deposit does not record.

**One duplicate assembly dropped.** RefSeq holds two assemblies of
*Salinispora tropica* CNS416, both at scaffold level. The strain is the unit of
every occurrence and co-occurrence count, so one had to go or that strain would
be counted twice in the linkage matrix. The rule is recorded rather than
silently applied: better assembly level, then smallest accession.

## The link score rewards ubiquity

The score is the raw Metcalf weighted co-occurrence count from Doroghazi et al.
2014, which NPLinker implements as its `metcalf` function. It is **not**
NPLinker's *standardised* Metcalf, which subtracts the score's expected value
under a null model and divides by its standard deviation. An earlier draft of
this task called it standardised; a hand-check of the arithmetic caught that.

The raw score is what a graded artifact wants — a closed-form weighted sum with
nothing left to a judgement call. But it has a known bias, and it is visible
immediately: on a 26-strain smoke run the top-ranked links were a gene cluster
family present in all 26 strains against molecular families present in all 26,
scoring the maximum 260 whether or not they are related. Standardising against
a null model is the correction, and the obvious extension for a later revision.
Grading the raw score does not depend on the bias being absent; it depends on
the arithmetic being reproducible, which it is.

## Measured per-stage cost

Everything here was measured on this container rather than estimated.

| stage | cost |
|---|---|
| fetching one mzXML | 7 s cold, ~30 s sustained — the GNPS2 proxy throttles and returns 429 on a third rapid request, so the stager pauses 8 s between files and backs off to 240 s |
| reading and filtering one mzXML | ~2 s, 1,778 MS2 scans in the file probed |
| GECCO per assembly | 49 s on 2 of 4 vCPUs under contention, against task 02's measured 36 s on 4 idle ones |
| presence calls per strain | ~5 s, about 35,000 modified-cosine comparisons |
| consensus network | ~100 s, as task 01 |

The filter keeps almost everything: a typical file gives 1,627 usable MS2 scans
from 1,627, because after the precursor window and the top-50 cut nearly every
profile-mode scan still clears the six-peak floor. So the **six-peak floor is
very nearly a no-op on these spectra**, the same honest negative task 01
records for its consensus set. It stays because the artifact needs a defined
rule, not because it is hard.

## Staging had to stream

33.5 GB does not fit on this container beside everything else, and neither did
either of the figures it was estimated at. `scripts/stage_task09_spectra.py` therefore
fetches a batch, records bytes and sha256, runs the solution's own stage-2
extraction over it, and deletes the mzXML before fetching the next, keeping
peak disk at a few hundred megabytes. Deleting is safe because the manifest is
the record: `npbench fetch` restages any file and verifies its checksum.

That is a build-time convenience, not part of the contract. An agent attempting
the task stages all 111 at once.

One consequence is recorded rather than hidden: the difficulty sweep cannot
model the top-50 peak cut, because changing it means re-extracting from mzXML
this build no longer holds. Two variants for that lever were written and then
removed — a variant that silently measures nothing is worse than a documented
gap, which is the lesson task 07's provably-redundant column taught. An auditor
with the full staging can add it back.

## Reproducibility

Verified. Two independent runs of stage 3 produce byte-identical output for all
seven artifacts, including the 275,988-row link table.

## Scale

| | |
|---|---|
| paired strains | 111, from PoDP's 360 genome-to-file links |
| staged | 33.52 GB — 111 spectrum files (33.33 GB), 111 genomes, one consensus `.mgf` |
| MS2 scans per strain | median 1,510, of which 99.9% survive filtering |
| consensus spectra | 4,081 retained of 100,263 |
| present per strain | 226 to 908, median 605 |
| occurrence rows | 65,886 |
| BGCs | 2,730 over 111 strains |
| gene cluster families | 218 — 115 confined to one strain, one present in all 111 |
| molecular families | 2,581, of which 1,266 are detected in at least one strain |
| candidate links | 275,988; the ranked cut keeps 500, spanning 143 distinct ranks |
| artifacts / metric instances | 7 / 27 |
| runtime | ~90 min GECCO, ~2 s per mzXML, ~5 s per strain for presence, ~100 s network |

## Reproducibility

Verified. Two independent runs of stage 3 produce byte-identical output for all
seven artifacts, including the 275,988-row link table.

## Measured difficulty

`scripts/difficulty_task09.py`, 19 variants across the genomic, metabolomic and
linking layers. The presence calls and the consensus network are cached on
disk per configuration, and finished variants are appended to a results file
and skipped on a re-run, so a second pass is minutes. Same cheap-proxy caveat
as the other tasks.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: keeps thin clusters | 1.00 |
| one slip: plain cosine | 1.00 |
| one slip: drops untyped clusters | 1.00 |
| one slip: domain count with duplicates | 0.86 |
| one slip: splits hybrid types | 0.86 |
| one slip: family needs every spectrum | 0.86 |
| one slip: no +1 for neither | 0.86 |
| one slip: scores undetected families | 0.86 |
| one slip: ordinal rank, not dense | 0.86 |
| one slip: 2.0 Da presence window | 0.86 |
| one slip: symmetric Metcalf | 0.71 |
| one slip: domains as a multiset | 0.57 |
| one slip: families numbered by rank | 0.57 |
| one slip: Jaccard 0.30 | 0.57 |
| one slip: pairs, not single linkage | 0.57 |
| one slip: presence at cosine 0.50 | 0.57 |
| hasty agent (3 slips) | 0.57 |
| careless agent (5 slips) | 0.43 |

| artifact | pass rate | |
|---|---|---|
| `top_links.csv` | 0.37 | in band |
| `metcalf_scores.csv` | 0.53 | in band |
| `gene_cluster_families.csv` | 0.68 | in band |
| `strain_inventory.csv` | 0.89 | **above band** |
| `bgc_inventory.csv` | 0.89 | **above band** |
| `spectrum_occurrence.csv` | 1.00 | **above band** |
| `molecular_families.csv` | 1.00 | **above band** |

**Three of seven in band**, which puts this task with task 03 rather than with
04, 05 and 08. The deep artifacts work exactly as the design intended — the
ranked links are the most discriminating artifact in the whole suite at 0.37,
the raw score table next at 0.53 — but four shallow artifacts do not fail.

What the sweep establishes:

- **The Metcalf asymmetry is the method.** Scoring a strain that has the gene
  cluster but not the molecule as 0 instead of −10 leaves all five upstream
  artifacts at 1.00 and takes the score table to 0.13 and the ranked links to
  0.16. Nothing else in the linking layer comes close.
- **The ranking is far more sensitive than the table it ranks.** Dropping the
  +1 for a strain with neither barely moves `metcalf_scores.csv` (0.91) and
  takes `top_links.csv` to 0.28; ranking ordinally instead of densely leaves
  the table at 1.00 and the ranking at 0.24. A 500-row ordered cut is a much
  sharper instrument than a 275,988-row table.
- **The genomic grouping rule reaches both halves.** Naming families by rank,
  using the raw domain multiset, or loosening Jaccard to 0.30 each floor
  `gene_cluster_families.csv` at 0.00 and carry through to the links.

### The artifact that was unmeasured, not unfailable

`bgc_inventory.csv` came out at **1.00 across every variant** in the first
16-variant sweep, and the tempting read was that GECCO's output is simply not
something an agent gets wrong. It was not: **no variant in that set touched the
layer the artifact lives on.** Every genomic mistake modelled was about
*grouping* clusters, not about *reading* them, because re-running GECCO is 90
minutes and the sweep reuses one cached run.

Three parsing variants fixed that without re-running anything — dropping
untyped clusters, reporting the domain field's length rather than its distinct
size, splitting hybrid product types — and the artifact moved to 0.89, failing
outright on the hybrid-type split. Still above band, but no longer an artifact
nothing could reach.

This is the third time the same diagnosis has applied: task 01's `spectra.csv`,
task 07's reference table, and now this. **An above-band artifact is a claim
about the variant set before it is a claim about the artifact**, and the
question to ask first is which layer the artifact lives on and whether anything
in the sweep touches it.

### Why the other three stay above band

- **`spectrum_occurrence.csv` (1.00) and `molecular_families.csv` (1.00)** do
  move — the looser presence threshold takes them to 0.86 and 0.81 — but never
  below 0.80. They have 65,886 and 2,581 rows, and identifier F1 over row sets
  that large is a blunt instrument: a mistake has to change a large fraction of
  the rows to fail, and these mistakes change a few percent. This is the same
  finding task 07 recorded for its 36,422-row reference table, now confirmed a
  third time, and it is a property of the metric and the threshold rather than
  of the artifact.
- **`strain_inventory.csv` (0.89)** fails only on the presence threshold, which
  is the one modelled mistake that reaches its counts.

I did not add columns to any of them. These three carry unique and central
information — the occurrence matrix is the thing this task exists to produce —
so the honest reading is that the 0.80 threshold is blunt at this row count,
which is a question for the calibration PROJECT.md §8 prescribes rather than
for a contrived column.

### Three variants that cost nothing

- **The thin-cluster exclusion is a no-op**, scoring 1.00 everywhere. Tasks 02
  and 05 found the same thing — 2 of 2,935 clusters affected — and this is the
  third confirmation on a third cluster set. The rule stays because the
  artifact needs a defined answer.
- **Plain cosine in place of modified cosine is a no-op**, and the reason is
  specific and worth knowing: the presence test only compares spectra whose
  precursors are within 0.5 Da, and the modified cosine's extra term exists to
  match fragments across a *precursor mass shift*. Inside a 0.5 Da window there
  is essentially no shift to model, so the two scores agree. **The "modified"
  in modified cosine earns its keep in the networking step, where precursors
  may differ by up to 2 Da, and not in the presence test.** An agent that uses
  the plain cosine here is not making a mistake that matters.
- **Dropping untyped clusters is nearly a no-op** (0.82 on the inventory,
  passing), because GECCO types almost everything it reports in this set.

### Across the nine tasks

| | artifacts in band | 3 slips | 5 slips |
|---|---|---|---|
| 01, molecular networking | 5 of 5 hard | 0.20 | 0.00 |
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| 05, shared joins, two inputs | 7 of 7 | 0.29 | 0.00 |
| 06, RiPP mining | 4 of 7 | 0.43 | 0.29 |
| 07, dereplication | 6 of 8 | 0.12 | 0.25 |
| 08, KS phylogeny | 8 of 8 | 0.75 | 0.00 |
| **09, this task** | **3 of 7** | **0.57** | 0.43 |

Read the last two columns with the caution the README states: they count
artifacts clearing 0.80, so they are not monotone in the number of mistakes and
not comparable across tasks.
