# CLAUDE.md

Operating instructions for this repository. Read `PROJECT.md` first; it explains what the benchmark is and why each design decision was made. This file tells you how to work in the repo.

---

## How to work

### 1. Think before coding

Don't assume. Don't hide confusion. Surface tradeoffs.

Before implementing:

- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them. Don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity first

Minimum code that solves the problem. Nothing speculative.

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical changes

Touch only what you must. Clean up only your own mess.

When editing existing code:

- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it. Don't delete it.

When your changes create orphans:

- Remove imports, variables and functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: every changed line should trace directly to the user's request.

### 4. Goal-driven execution

Define success criteria. Loop until verified.

Transform tasks into verifiable goals:

- "Add validation" becomes "Write tests for invalid inputs, then make them pass"
- "Fix the bug" becomes "Write a test that reproduces it, then make it pass"
- "Refactor X" becomes "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:

```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

These guidelines are working if: fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

### What these mean here specifically

A benchmark is graded, so a plausible-looking grader that is subtly wrong is worse than no grader. When you implement a metric, write the fixture first: a hand-computed input and expected score. When a source paper's method section is ambiguous about a parameter, that is a question for Saman, not a default you pick. When a download fails or an accession does not resolve, say so and stop; do not substitute a different dataset.

---

## What this repo is

Task bundles and a grading harness for NPBench-R, a benchmark in which AI agents reproduce published natural products and BGC analyses from raw data. Each task stages raw data, prescribes the method, and specifies output artifacts that are compared to the published artifacts with deterministic metrics.

The benchmark is not a model, not a pipeline, and not a Mycel component. It must remain runnable and credible independent of any system it evaluates.

## Layout

```
npbench-r/
  PROJECT.md              # design spec, read first
  CLAUDE.md               # this file
  tasks/<task_id>/        # one directory per task, contract in PROJECT.md §6
  npbench/
    grading/              # metric implementations
    harness/              # Inspect AI integration, scorer, adjudicator
    build/                # data fetch, manifest, bundle validation
  environments/           # Dockerfiles and conda locks, one per tool stack
  tests/                  # fixtures and unit tests for graders
  scripts/                # one-off operational scripts, not library code
```

`tasks/*/data/`, `tasks/*/reference/` and `tasks/*/truth/` are never committed. They are fetched by manifest and verified by checksum. Git holds the manifests, not the payloads.

## Build order

Work in this order. Each step has a verification that must pass before the next begins.

1. **Grading library.** Implement the four metrics in `npbench/grading/`: identifier F1, Lin's CCC, macro F1, interval overlap F1. Plus the product combiner and the 0.80 pass rule.
   Verify: unit tests with hand-computed fixtures, including degenerate cases (empty file, all rows missing, perfect match, constant-offset numeric column that should score well on Pearson and poorly on CCC).

2. **Bundle validator.** A function that takes a task directory and checks the contract: every file in `grading.yaml` exists in `truth/`, every manifest entry has a checksum, every required output has a declared metric and row key.
   Verify: validator passes on a hand-built minimal fixture task and fails with a clear message on each of five deliberately broken variants.

3. **One task end to end.** Build task 1 (Crüsemann, MassIVE MSV000078836). Fetch data, write the reference solution script, generate `truth/`, write the prompt, grade the reference solution against its own truth.
   Verify: the reference solution scores 1.00 on every artifact. If it does not, the grader or the truth is wrong, and that must be resolved before any other task is built.

4. **Harness.** Inspect AI ReAct loop, the five tools (bash session, python, text editor, web access request, submit), the network adjudicator, the 24-hour and message caps, and the scorer that reads `outputs/` off the VM.
   Verify: a deliberately failing agent run produces a scored result of 0.00 with no crash; the reference solution run through the harness scores 1.00.

5. **Remaining tasks.** Build tasks 2 through 6 in the order given in `PROJECT.md` §5.
   Verify per task: reference solution scores 1.00; bundle validator passes; total staged data under 20 GB; reference solution completes in under 8 hours on 32 vCPUs.

6. **Calibration.** Run two or three frontier models twice per task. Drop artifacts whose pass rate is outside roughly 0.1 to 0.8.
   Verify: top model task score lands between 0.3 and 0.5.

Do not start step N+1 while step N's verification is failing. Say so and stop.

## Hard rules

**Environments are pinned or they are wrong.** Every tool gets an exact version in a conda lock file or a Dockerfile with a digest-pinned base image. `antismash==7.1.0`, not `antismash`. This is not pedantry: antiSMASH detects different cluster type counts across major versions, and BiG-SCAPE GCF counts depend on the distance cutoff. A ground truth artifact is only valid against the exact version that produced it. Record the version in `task.yaml` and in the method guidance.

**Reference databases are frozen inside the bundle.** MIBiG, the antiSMASH database, GNPS library exports, NP Atlas. Snapshot, checksum, record the release. Never let a task resolve a database over the network at runtime.

**Never grade a live web service.** GNPS spectral library matching and SIRIUS/CSI:FingerID structure search both depend on server-side state that changes. If a task needs them, freeze the library export and grade against that, or mark the artifact soft and exclude it from the headline score. Network topology from a fixed `.mgf` with fixed parameters is fine.

**Raw data never enters git.** Manifests with URLs, accessions, byte counts and sha256. A fetch script that verifies. If an accession stops resolving, the manifest is the record of what we had.

**Truth is generated by a committed script.** Every file in `truth/` is produced by a script in `solution/`. If you cannot regenerate an artifact from a committed script, it is not usable as ground truth. Where a source paper's artifact involved manual curation, do not use it; grade only the tool-generated artifacts.

**Normalise identifiers before comparison.** Case, whitespace, accession formatting, and a declared per-artifact alias map. Report format compliance as a separate axis from the artifact score. A model that produces correct science with wrong column names should be visible as exactly that, not scored as a failure.

**Do not leak answers.** The adjudicator denies the source paper, its supplementary files, its processed matrices, and search queries containing paper-specific terms or accessions. It allows package registries, tool source, documentation and canonical reference data. When you add a task, check whether its truth is reachable through any allowed path and widen the denylist if so.

## Commands

```bash
make lock              # regenerate conda locks from environments/*.yml
make test              # unit tests for graders and validator
make validate TASK=01  # bundle contract check
make fetch TASK=01     # download and checksum staged data
make truth TASK=01     # run reference solution, write truth/
make selfgrade TASK=01 # grade reference solution against truth, expect 1.00
make eval TASK=01 MODEL=...  # full harness run
```

Keep `scripts/` for one-off operational work. Anything called more than twice belongs in `npbench/` with a test.

## What to ask about rather than decide

- Which artifact to drop when a source paper's supplementary table turns out to be figure-only or PDF.
- Any parameter the source paper did not state. Reconstructing it from context is a judgement call with downstream consequences for the whole task.
- Whether an artifact that only one tool version can reproduce is worth keeping.
- Substituting a dataset or accession when the one in `PROJECT.md` §5 does not resolve.
- The benchmark's licence, and anything touching the NP Atlas CC BY-NC term.

## Known unresolved items

Carried from `PROJECT.md` §5 and §13. These are open, not forgotten:

- Crüsemann supplementary file format and total MassIVE download size, unverified.
- mBio vertical inheritance Table S1 and S2 file types, and the OSF deposit contents, unverified.
- Leão cyanobacteria accession discrepancy: MSV000085038 in the NPLinker paper versus MSV000084950 in the NPClassScore paper.
