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
| 3 | One task end to end | task 02 bundle built; truth generation in progress |
| 4 | Inspect AI harness | built; docker sandbox path unverified |
| 5 | Remaining tasks | not started |
| 6 | Difficulty calibration | not started |

Build order was reordered: **task 02 (Salinispora, BGC detection) was built
before task 01 (Crüsemann)**. MassIVE MSV000078836 turned out to be 311 GB
across 1,125 mzXML files, which is 15x the 20 GB target in PROJECT.md §10, so
task 01 cannot be staged as specified and needs a subsetting decision. Task 02
stages from NCBI RefSeq in 191 MB.

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
