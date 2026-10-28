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
| 2 | Bundle validator | done |
| 3 | One task end to end | **done**: task 02 validates clean and self-grades 1.00 |
| 4 | Inspect AI harness | built; docker sandbox path unverified |
| 5 | Remaining tasks | 4 of 6 built (02, 03, 04, 05) |
| 6 | Difficulty calibration | machinery built (`npbench calibrate`); needs model runs |

Build order was reordered: **task 02 (Salinispora, BGC detection) was built
before task 01 (Crüsemann)**. MassIVE MSV000078836 turned out to be 311 GB
across 1,125 mzXML files, which is 15x the 20 GB target in PROJECT.md §10, so
task 01 cannot be staged as specified and needs a subsetting decision. Task 02
stages from NCBI RefSeq in 191 MB.

## Tasks

| | |
|---|---|
| `02_salinispora_bgc` | BGC inventory and gene cluster families across 120 *Salinispora* assemblies with GECCO 0.11.0. 191 MB staged, 2,935 clusters, 7 artifacts, 30 metric instances, depth 2–4. |
| `03_mibig_cheminformatics` | Scaffolds, isomers, stereochemistry collapse and similarity families over the MIBiG 4.0 compound set with RDKit. 936 KB staged, 4,397 compounds, 8 artifacts, 31 metric instances, depth 2–4. **Discriminates less well than the others** — see its `NOTES.md`. |
| `04_mibig_nrps_pks_modules` | NRPS/PKS module architecture and adenylation-domain substrate specificity in MIBiG 4.0. 936 KB staged, 4,034 modules, 8 artifacts, 39 metric instances, depth 2–4. Artifacts share upstream dependencies by design; all 8 land in band, but one scope mistake floors the task at 0.00. |
| `05_salinispora_concordance` | Taxonomy versus chemistry concordance over the same 120 *Salinispora* assemblies as task 02: a k-mer genomic distance against a BGC-family chemistry distance over 7,140 pairs. 7 artifacts, 32 metric instances, depth 2–5. **Best-behaved task so far** — all 7 in band and no floor. |

All three validate clean and self-grade 1.00. Each task's `NOTES.md` records
measured scale, reproducibility and difficulty, including what does *not*
discriminate.

### The main design finding so far

A task discriminates when its artifacts **share upstream dependencies**, not
when it has many artifacts. Measured across the three:

| | artifacts in the 0.1–0.8 band | 3 slips | 5 slips |
|---|---|---|---|
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| 05, shared joins over two inputs | 7 of 7 | 0.29 | 0.00 |

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

All four `NOTES.md` files record what does *not* discriminate as well as what
does.

Task 01 (Crüsemann, MassIVE MSV000078836) is not built: the deposit is 311 GB
across 1,125 mzXML files, 15x the 20 GB target in PROJECT.md §10, so it needs a
subsetting decision before it can be staged.

## Usage

```bash
pip install -e ".[dev]"                       # graders, validator, harness
T=02_salinispora_bgc

make test                                      # 198 tests
make validate TASK=$T                          # bundle contract check
make fetch TASK=$T                             # stage and checksum 191 MB of genomes
make truth TASK=$T                             # run the reference solution into truth/
make selfgrade TASK=$T                         # expect 1.00
make grade TASK=$T OUTPUTS=<dir>               # grade any outputs directory
make eval TASK=$T MODEL=anthropic/claude-opus-5
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
