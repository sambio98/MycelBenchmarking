# NPBench-R

A computational reproducibility benchmark for natural products and biosynthetic
gene clusters. Agents are given raw deposited data and prescribed methods, and
must produce structured artifacts that are graded against the ones the original
authors published.

Read [`PROJECT.md`](PROJECT.md) for the design and the reasoning behind it, and
[`CLAUDE.md`](CLAUDE.md) for how to work in this directory.

## State

Build order is set out in `CLAUDE.md`. Steps 1 and 2 are done and verified.

| Step | | |
|---|---|---|
| 1 | Grading library | done |
| 2 | Bundle validator | done; `npbench audit` adds the checks for a human curation pass |
| 3 | One task end to end | **done**: task 02 validates clean and self-grades 1.00 |
| 4 | Inspect AI harness | built; docker sandbox path unverified |
| 5 | Remaining tasks | all 6 built (01–06) |
| 6 | Difficulty calibration | machinery built (`npbench calibrate`); needs model runs |

Build order was reordered: **task 02 (Salinispora, BGC detection) was built
before task 01 (Crüsemann)**. MassIVE MSV000078836 turned out to be 311 GB
across 1,125 mzXML files, 15x the 20 GB target in PROJECT.md §10, so task 01
waited on a scope decision while task 02 staged from NCBI RefSeq in 191 MB.
Task 01 is now built against the GNPS-clustered consensus spectra from the same
deposit — 17 MB, and the artifact PROJECT.md §9 names — rather than the raw
mzXML. That reduction is the one decision in the suite most worth a second
opinion; its `NOTES.md` states what it costs and what it keeps.

## Tasks

| | |
|---|---|
| `01_crusemann_networking` | Molecular networking over the Crüsemann *Salinispora*/*Streptomyces* MS/MS deposit (MassIVE MSV000078836), from modified-cosine scoring through the top-10 cut to component structure. 17 MB staged, 4,081 spectra, 3,208 edges, 325 components, 6 artifacts, 27 metric instances, depth 2–5. All 5 hard artifacts in band; grades the GNPS-clustered spectra rather than the 311 GB of raw mzXML — see its `NOTES.md`. |
| `02_salinispora_bgc` | BGC inventory and gene cluster families across 120 *Salinispora* assemblies with GECCO 0.11.0. 191 MB staged, 2,935 clusters, 7 artifacts, 30 metric instances, depth 2–4. |
| `03_mibig_cheminformatics` | Scaffolds, isomers, stereochemistry collapse and similarity families over the MIBiG 4.0 compound set with RDKit. 936 KB staged, 4,397 compounds, 8 artifacts, 31 metric instances, depth 2–4. **Discriminates less well than the others** — see its `NOTES.md`. |
| `04_mibig_nrps_pks_modules` | NRPS/PKS module architecture and adenylation-domain substrate specificity in MIBiG 4.0. 936 KB staged, 4,034 modules, 8 artifacts, 39 metric instances, depth 2–4. Artifacts share upstream dependencies by design; all 8 land in band, but one scope mistake floors the task at 0.00. |
| `05_salinispora_concordance` | Taxonomy versus chemistry concordance over the same 120 *Salinispora* assemblies as task 02: a k-mer genomic distance against a BGC-family chemistry distance over 7,140 pairs. 7 artifacts, 32 metric instances, depth 2–5. **Best-behaved task** — all 7 in band and no floor. |
| `06_streptomyces_ripp` | RODEO-style lasso peptide and RiPP locus mining across 213 *Streptomyces* reference genomes, with a sequence similarity network over the precursors. 517 MB staged, 171 loci, 798 precursors, 7 artifacts, 37 metric instances, depth 2–5. 4 of 7 in band; two artifacts are drop candidates — see its `NOTES.md`. |
| `07_npatlas_dereplication` | Exact-mass dereplication of task 01's molecular network against a frozen NP Atlas v2024_09, with adducts, isomer ambiguity and per-family annotation. 493 MB staged, 36,422 compounds, 3,811 queries, 10,953 candidate annotations, 8 artifacts, 35 metric instances, depth 2–5. 6 of 8 in band. |
| `08_salinispora_ks_phylogeny` | Ketosynthase phylogeny across the 120 *Salinispora* assemblies: pyrodigal and pyhmmer to the domains, normalised local alignment to a distance matrix, UPGMA to a tree graded by its splits. 182 MB staged, 2,993 domains collapsing to 922 distinct sequences, 8 artifacts, 36 metric instances, depth 2–5. All 8 in band. |
| `09_salinispora_paired_omics` | Metabologenomic link ranking across 111 *Salinispora* strains with both a genome and an LC-MS/MS run, paired by the Paired Omics Data Platform. GECCO clusters and gene cluster families on one side, strain-resolved metabolite occurrence recovered from each strain's MS2 against the deposit's consensus spectra on the other, joined by Metcalf co-occurrence scoring. 33.5 GB staged, 7 artifacts, 27 metric instances, depth 2–5. 3 of 7 in band; its ranked link table is the most discriminating artifact in the suite. **The largest fetch in the suite at 1.7× §10's ideal — see its `NOTES.md`.** |

Tasks 01–08 validate clean, self-grade 1.00 and report no open `npbench audit`
findings. Task 09 is the newest and its verification is in progress; its
`task.yaml` records determinism as pending rather than claiming it. Each task's `NOTES.md` records
measured scale, reproducibility and difficulty, including what does *not*
discriminate.

## Coverage of the eight analysis categories

PROJECT.md §4 asks for six to eight tasks and forty to sixty graded artifacts
covering eight categories. The suite is at **nine tasks, 66 artifacts and 294
metric instances**, one task and six artifacts past the top of that range, and
**all eight categories are covered**.

| category | covered by |
|---|---|
| 1. BGC detection and annotation | 02 (GECCO), 06 (pyrodigal + pyhmmer), 09 |
| 2. Gene cluster family clustering | 02, 05, 06, 09 |
| 3. MS/MS molecular networking | 01, 07, 09 |
| 4. Dereplication and structure annotation | 07 |
| 5. Paired omics and link ranking | 09 |
| 6. NRPS and PKS domain parsing | 04, 08 |
| 7. Cheminformatics on NP structures | 03 |
| 8. Phylogenetics and taxonomy/chemistry concordance | 05 (concordance), 08 (phylogeny) |

Being over the task count is a curation decision, not an accident: §4's range
was written before the categories were costed, and category 5 needed a task of
its own. Dropping one is a reasonable call for the release and each task's
`NOTES.md` carries the evidence for making it.

Within category 2, *novelty scoring* specifically (BiG-SLiCE-style) is still
not covered; the clustering half is.

### How category 5 got unblocked

This section previously recorded category 5 as unbuildable, and the reasoning
was sound but one source short. Worth keeping as a record of the dead ends,
because they are real and anyone rebuilding this will hit them:

- **MSV000078836's clustered spectra carry no per-strain provenance.** Every
  one of the 100,263 spectra in the staged `.mgf` has `FILENAME=specs_ms.pklbin`,
  so the file cannot say which strain a consensus spectrum came from.
- **MassIVE's own endpoints are unreachable from this container.** Its HTTPS
  chain roots in an emSign CA present in neither the system trust store nor
  `certifi`, and the GNPS2 proxy that works instead serves only exact paths and
  returns no directory listing. The GNPS cluster-info tables, which would carry
  the per-file mapping, 404 at every documented name tried.
- **MSV000098188's peak files 404 at every layout tried**, and its 101
  "microscale" strains have zero overlap with the sequenced set.

What cleared it was the **Paired Omics Data Platform**, which PROJECT.md §5
names for this very dataset and which I had not tried. Its project record for
MSV000078836 curates exactly the missing mapping: 360 genome-to-file links, 120
strains across three extraction methods, each link naming the deposit path of
one strain's `.mzXML`. 111 of those strains already had an NCBI assembly staged
for task 02. That is task 09.

The lesson generalises past this task: **"unavailable" was a claim about the
sources I had tried, not about the data.** The three dead ends above are all
about file *discovery* on one archive, and the mapping existed in a different
archive built for precisely this purpose.

## What the difficulty measurements show

### The main design finding

A task discriminates when its artifacts **share upstream dependencies**, not
when it has many artifacts. Measured across the suite:

| | artifacts in the 0.1–0.8 band | 3 slips | 5 slips |
|---|---|---|---|
| 01, molecular networking | 5 of 5 hard | 0.20 | 0.00 |
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| 05, shared joins over two inputs | 7 of 7 | 0.29 | 0.00 |
| 06, RiPP mining | 4 of 7 | 0.43 | 0.29 |
| 07, dereplication | 6 of 8 | 0.12 | 0.25 |
| 08, KS phylogeny | 8 of 8 | 0.75 | 0.00 |
| 09, paired omics | 3 of 7 | 0.57 | 0.43 |

**Read the last two columns with care.** Task 07 scores *higher* at five
slips than at three, and that is not noise: the headline is the fraction of
artifacts clearing 0.80, so what matters is how many artifacts a mistake set
touches, not how many mistakes it contains. The composite columns are
therefore not monotone in the number of mistakes and not comparable across
tasks unless the mistake sets happen to cover each dependency graph the same
way. The per-artifact pass rate is the figure to trust.

Task 04 showed the opposite failure to task 03: one scope mistake floors it at
0.00, and a task that bottoms out loses information as surely as one that tops
out. Task 05 spreads the dependencies across two independent inputs — a k-mer
parse and a cluster parse — and gets every artifact informative with a
three-mistake agent at 0.29, nearest the 0.3–0.5 target in PROJECT.md §8.

A second finding, from trying the same fix on two tasks: a content-derived
column (a group's smallest member, its most populous family) only improves
discrimination when some mistake actually gets the grouping or naming wrong. It
moved task 05's weakest artifact into band and did nothing at all for task
03's. Adding one on principle buys nothing.

A third finding, from task 06: some rules cannot discriminate no matter what,
because the tools make them unreachable. Its 20-amino-acid precursor floor can
never bite, since pyrodigal's own minimum gene length puts the shortest
precursor at 29; and its E-value is dominated by HMMER's inclusion thresholds.
Those rules stay in the method because the artifacts need a defined answer, not
because they are hard, and its `NOTES.md` says so rather than counting them as
difficulty.

Task 08 independently replicated the E-value half of that on different genomes
and a different Pfam model: searching at 1e-3 instead of 1e-5 changes nothing
at all, and keeping every reported hit rather than only the included ones
changes almost nothing. HMMER's own inclusion thresholds dominate. Two
measurements agreeing is worth more than either alone, and it means an
HMMER-based task should not count its E-value as a difficulty lever.

A fourth finding, from task 01: **an above-band verdict is a claim about the
variant set at least as much as about the artifact.** Its `spectra.csv` first
measured 0.87, above band, and the tempting read was that it needed a harder
column. It did not. Every variant in that sweep had been drawn from the
pipeline's parameters, while that artifact lives almost entirely on the parsing
layer in front of them. Two realistic misreadings — a missing `CHARGE` read as
1+, a base peak read as highest m/z rather than most intense — take it to 0.57
and 0.17, and moved it into band with no change to the bundle. When that check was first made, the other five sweeps did not have the hole.
Tasks 07 and 09 later did: task 07's 36,422-row reference table absorbed
mistakes that a per-row label would have caught, and task 09's BGC inventory
read 1.00 across sixteen variants purely because none of them touched the layer
it lives on — adding three GECCO-parsing variants moved it to 0.89. **An
above-band artifact is a claim about the variant set before it is a claim about
the artifact**, and the first question is which layer the artifact lives on and
whether anything in the sweep reaches it.

All six `NOTES.md` files record what does *not* discriminate as well as what
does.

## Usage

```bash
pip install -e ".[dev]"                       # graders, validator, harness
T=02_salinispora_bgc

make test                                      # 282 tests
make validate TASK=$T                          # bundle contract check
make fetch TASK=$T                             # stage and checksum 191 MB of genomes
make truth TASK=$T                             # run the reference solution into truth/
make selfgrade TASK=$T                         # expect 1.00
make grade TASK=$T OUTPUTS=<dir>               # grade any outputs directory
make eval TASK=$T MODEL=anthropic/claude-opus-5

make digest                                    # refingerprint truth/ after regenerating it
make audit                                     # the human curation pass, all tasks
make audit TASK=$T                             # one task's dossier
```

Grading needs only numpy, pandas and pyyaml. Inspect AI is an extra
(`pip install -e ".[harness]"`) so the graders stay installable without
pulling in a hundred packages.

## Grading model

Four metrics, selected per artifact, multiplied into one artifact score. An
artifact passes at 0.80. The task score is the fraction of hard artifacts that
pass; soft artifacts are scored and reported but excluded from the headline.

| Metric | Used for |
|---|---|
| `identifier_f1` | row recovery: gene lists, BGC inventories, feature IDs |
| `lins_ccc` | numerical agreement, optionally log transformed |
| `macro_f1` | categorical agreement over aligned rows |
| `interval_overlap_f1` | BGC coordinates, under a reciprocal overlap rule |

Three rules are worth knowing before you write a `grading.yaml`.

**Nothing has a default.** Every parameter the spec leaves open — the
reciprocal overlap threshold, the log transform, the pseudocount, the
normalisation rules, the alias maps — must be stated per artifact or the
validator fails the bundle. A score computed under a threshold nobody chose is
the failure mode this is built to prevent. The 0.80 pass threshold is the one
exception: it is a property of the benchmark, not of a task, and lives in
`npbench/grading/config.py`.

**Aligned metrics need `identifier_f1` beside them.** CCC and macro F1 only see
the rows the submission and truth share, so on their own they would not notice
a submission that dropped half the rows. The validator rejects an artifact that
declares one without the other, and the product combiner carries the penalty.

**Duplicate row keys fail the artifact.** Row alignment is undefined when a key
appears twice, and keeping the first row silently would reward the bad-join
inflation that PROJECT.md §8 wants to catch. A duplicate key in `truth/` raises
instead, because that is a bundle bug.

Identifier normalisation runs before every comparison. Format compliance —
whether a column was found under its declared name or through an alias — is
reported as a separate axis and never costs artifact score, so a model that
produces correct science with wrong column names is visible as exactly that.

## The answer key is fingerprinted, not committed

`tasks/*/truth/` is gitignored. Truth is generated by a committed script from
checksummed inputs, and CLAUDE.md keeps raw and derived payloads out of git.

The gap that leaves is that the answer key becomes unverifiable: a dependency
bump, an edited solution script or a silently different reference release would
all produce a different `truth/`, and nothing would say so — while every score
the benchmark had reported was computed against the old one.

So the payload stays out and the fingerprint goes in. `tasks/*/truth.digest.json`
holds one sha256 and one row count per artifact, `make validate` checks it
whenever both the digest and a generated `truth/` are present, and `make digest`
rewrites it when truth is regenerated on purpose. Regenerating truth and finding
the digest unchanged is the evidence that a published score still means what it
meant.

**Whether to commit `truth/` itself is a release decision, not a technical one**,
and it is open. Shipping it would let anyone grade without first staging up to
500 MB per task and spending minutes to an hour regenerating; keeping it out
keeps the repo small and makes the generating script the single source of truth.
The digest makes either choice safe.

## Layout

```
npbench/
  grading/     metrics, normalisation, grading.yaml parsing, aggregation
  build/       bundle contract validation
  harness/     Inspect AI integration (step 4, empty)
  cli.py
tasks/         one directory per task (step 3, empty)
tests/
  fixtures/minimal_task/   hand-built bundle the validator is tested against
```

`tasks/*/data/`, `tasks/*/reference/` and `tasks/*/truth/` are never committed.
Git holds the manifests, not the payloads.
