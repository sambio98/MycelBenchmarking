# NPBench-R: what is left

Everything the build produced is verified: nine tasks, 66 artifacts, 294 metric
instances, all eight PROJECT.md §4 categories, every bundle validating clean,
grading 1.00, byte-identical across two runs, zero open `npbench audit`
findings, 288 tests passing.

This is the list of what is *not* done, written so a reviewer can act on it
without rereading the history. Each item says who it is waiting on and where
the evidence lives.

A fresh clone was tested: 288 tests pass with no data, `make validate` reports
0 errors on an unstaged bundle, and `make fetch` stages and checksum-verifies.
See the README's "Starting from a fresh clone" for the sequence and the time and
disk it needs. Total staging across all nine tasks is 35.1 GB, of which task 09
is 33.5 GB.

---

## 1. Blocked on access or hardware

Nothing here needs a decision. Each is doable as written once the resource
exists.

### 1.1 Run the docker sandbox path of the Inspect harness
**Waiting on:** a machine with a docker daemon.
**Evidence:** `npbench/harness/task.py` docstring, README build-order step 4.

The harness is verified end to end against Inspect's mock model with
`--sandbox local`: bundle loads, sample built from `prompt.md`, ReAct loop runs
with its tools, scorer returns per-artifact reasons and metadata. Running it
found two defects that no other test could reach — see
`tests/test_harness.py::TestTheDocumentedEntrypointLoads`.

What is still unexecuted is `sandbox_spec`'s compose and Dockerfile branches,
which only unit tests cover. On a provisioned machine:

```bash
make eval TASK=03_mibig_cheminformatics MODEL=mockllm/model   # local sandbox
inspect eval npbench/harness/task.py@npbench_r \
  -T task_dir=tasks/03_mibig_cheminformatics --model mockllm/model
```

The second command exercises the real docker sandbox. Expect to debug the
image build before the loop runs; every task's Dockerfile installs from its
pinned lock and smoke-tests its own tools, but none has been built.

### 1.2 Real difficulty calibration
**Waiting on:** model API access.
**Evidence:** `npbench/calibrate.py`, every task's `NOTES.md`.

`npbench calibrate` is built and withholds a verdict below four runs. Every
difficulty number in the repository is a **synthetic-variant proxy** — I
perturbed the reference solution and scored the result. That is not what
PROJECT.md §8 asks for and no `NOTES.md` claims otherwise.

Two cautions for whoever runs the real thing:

- **The 3-slip and 5-slip columns are not monotone and not comparable across
  tasks.** Task 07 scores *higher* at five slips than three, because the
  headline counts artifacts clearing 0.80 and what matters is how many
  artifacts a mistake set touches. Trust the per-artifact pass rate.
- **An above-band artifact is a claim about the variant set before it is a
  claim about the artifact.** This caught me three times — task 01's
  `spectra.csv`, task 07's reference table, task 09's `bgc_inventory.csv`. Ask
  which layer the artifact lives on and whether anything in the sweep reaches
  it.

### 1.3 antiSMASH and BiG-SCAPE tasks
**Waiting on:** docker or a conda/mamba channel.
**Evidence:** task 02's `NOTES.md` on why GECCO was chosen instead.

PROJECT.md §4 names both. Neither is pip-installable and the binary download
for micromamba was refused in this environment, so GECCO carries categories 1
and 2. antiSMASH also needs a ~9 GB database. Version sensitivity is the real
hazard: §9 records 58 cluster types in v5, 71 in v6, 101 in v8, so any
antiSMASH task is only reproducible against one pinned version.

---

## 2. Decisions only you can make

### 2.1 Re-derive the 0.80 pass threshold
**Evidence:** PROJECT.md §7, `npbench/grading/config.py` (`PASS_THRESHOLD`).

The constant is mine, not calibrated. §7 says to repeat BixBench3's procedure
with NP domain experts on our own artifacts, and warns that reviewers will ask.
PROJECT.md §13 leaves open whether the panel is lab-internal or external.

One finding from this build should feed into it: **identifier F1 over a large
row set is a blunt instrument at any fixed threshold.** Task 07's 36,422-row
reference table absorbs the loss of 2,925 distinct stereoisomers and still
scores 0.96; task 09's 65,886-row occurrence matrix never drops below 0.80
under any modelled mistake. A row-count-aware threshold, or a rule that large
artifacts must be graded on a per-row label rather than a row set, may be a
better answer than one global number.

### 2.2 Task 09 stages 33.5 GB
**Evidence:** `tasks/09_salinispora_paired_omics/NOTES.md`, its `task.yaml`
`data_scope`.

1.7× the 20 GB §10 calls ideal, inside its 100 GB limit, and sixty-seven times
the next largest task. I kept all 111 paired strains because "every paired
strain at one extraction method" is a rule and any smaller set needs an
arbitrary one. **Dropping to 66 strains fits under 20 GB** and Metcalf scoring
degrades gracefully. This is the decision in the suite I would most like a
second opinion on.

### 2.3 Nine tasks against §4's six-to-eight
**Evidence:** README coverage section.

§4's range was written before the categories were costed, and category 5
needed a task of its own. Dropping one is reasonable; each `NOTES.md` carries
the evidence. The weakest by measured discrimination are **task 03 (3 of 8
artifacts in band)** and **task 09 (3 of 7)**, and task 03 is the one with no
category of its own to defend.

### 2.4 Whether to ship `truth/`
**Evidence:** README "The answer key is fingerprinted, not committed".

`tasks/*/truth/` is gitignored; `tasks/*/truth.digest.json` holds one sha256
and row count per artifact and `make validate` checks it. Shipping truth would
let anyone grade without staging up to 33.5 GB and regenerating; keeping it out
keeps the repo small and makes the generating script the single source. The
digest makes either choice safe.

### 2.5 Benchmark licence, and the NP Atlas CC BY-NC term
**Evidence:** PROJECT.md §9 and §13, task 07's `NOTES.md`.

NP Atlas v2024_09 is CC BY-NC 4.0 and task 07 stages it. Whether that
constrains the benchmark's own licence is §13's open question and is unresolved.
Only the manifest ships, not the payload, which may or may not matter.

### 2.6 Where to host staged data
**Evidence:** PROJECT.md §13.

Hugging Face (BixBench3's choice) or Zenodo. The manifests make either work;
nothing in the code assumes a host. Note the GNPS2 proxy that task 01, 07 and
09 fetch through throttles to ~1 MB/s sustained and returns 429 on a third
rapid request, which is a reason to rehost rather than depend on it.

### 2.7 How many models to evaluate at release
**Evidence:** PROJECT.md §13.

Thirteen is thorough and expensive; six probably shows the depth and data-size
effects.

---

## 3. Known weaknesses, documented rather than fixed

These are measured, recorded in the relevant `NOTES.md`, and left alone
deliberately. PROJECT.md §8's prescribed remedy is to build more candidate
artifacts than needed and drop the uninformative ones at calibration, so they
are drop candidates rather than bugs.

| task | measured | first drop candidates |
|---|---|---|
| 03 | 3 of 8 in band | its three aggregate artifacts |
| 06 | 4 of 7 in band | `genome_proteins.csv`, then `domain_prevalence.csv` |
| 07 | 6 of 8 in band | `atlas_compounds.csv`, `queries.csv` |
| 09 | 3 of 7 in band | `spectrum_occurrence.csv`, `molecular_families.csv` |

**Task 04 has the opposite problem:** all 8 artifacts in band, but one scope
mistake floors the whole task at 0.00. A task that bottoms out loses
information as surely as one that tops out.

**Task 01's `singletons.csv` is soft** — it cannot discriminate and was handing
the 3-slip and 5-slip agents a free pass. Kept because the 2,256-of-4,081
split is a reported result.

---

## 4. Technical follow-ups

### 4.1 Task 09: the peak-cut lever is unmeasured
Changing the top-50 MS2 peak cut needs re-extraction from the raw mzXML, and
the build streams those away — fetch, fingerprint, reduce, delete, because 111
files at a mean 300 MB do not fit beside everything else. Two variants for it
were written and removed; a variant that silently measures nothing is worse
than a documented gap. An auditor with the full 33.5 GB staged can add it back
to `scripts/difficulty_task09.py`.

### 4.2 Task 09: standardised Metcalf
The link score is the raw weighted co-occurrence count, which **rewards
ubiquity**: the top-ranked links are a gene cluster family present in all 111
strains against molecular families present in all 111, scoring near the maximum
whether or not they are related. NPLinker's standardised variant subtracts the
expected value under a null model and divides by the standard deviation.
Implementing it is the obvious revision; I kept the raw score because it is a
closed-form sum with nothing left to a judgement call, which is what a graded
artifact needs.

### 4.3 Category 2's novelty scoring is uncovered
BiG-SLiCE-style novelty scoring against a reference set has no task. The
clustering half of category 2 is covered four times over (02, 05, 06, 09).
Comparing the Salinispora GCFs against MIBiG by domain content would do it and
needs no new data — MIBiG is already staged for tasks 03 and 04.

### 4.4 Honest negatives worth not re-discovering
Measured no-ops, each recorded in its task's `NOTES.md`. Do not count any of
them as difficulty:

- **HMMER's inclusion thresholds dominate the E-value.** Searching at 1e-3
  instead of 1e-5 changes nothing, in task 06 *and* task 09, on different
  genomes and different Pfam models. An HMMER task should not treat its E-value
  as a lever.
- **The thin-cluster exclusion is a no-op** — 2 of 2,935 clusters — confirmed
  on three cluster sets (02, 05, 09).
- **Plain cosine equals modified cosine inside a tight precursor window.** Task
  09's presence test only compares spectra within 0.5 Da, and the modified
  cosine's extra term exists to match fragments across a *mass shift*. It earns
  its keep in the networking step, not there.
- **Task 08's distance clamp cannot bite**, provably: a local alignment score
  is bounded below by 0 and above by the shorter sequence's self-score, so the
  distance is always in [0,1]. Checked on all 424,581 pairs.
- **Task 01's ±17 Da precursor window is invisible to its metrics** — it
  changes 81 of 4,081 spectra's peak counts, pushes none below the six-peak
  floor and changes zero base peaks.
- **Six-peak floors are nearly no-ops** on both task 01's consensus spectra and
  task 09's profile-mode scans, where 99.9% of filtered scans clear them.

### 4.5 Keep task.yaml's pins honest
`task.yaml`'s `tool_versions` is what a reader installs from; the lock is what
truth was actually produced under. Task 01 claimed `networkx: 3.7` while its
lock pinned 3.4.2, so anyone following the documentation would have installed a
version the answer key was never built against. Fixed, and
`tests/test_validate.py::test_every_task_yaml_pin_matches_its_lock` now fails on
any such drift for every pip-installable tool. Database releases
(`mibig: 4.0`, `npatlas: 2024_09`, `pfam: 35.0`) have no lock line and are not
covered — check those by hand when adding a task.

### 4.6 A reproducibility hazard to watch
NP Atlas serves its current release from an **unversioned URL**. Task 07 pins
the bytes by sha256 and its solution asserts the 36,454-record count before
doing any work, so drift fails loudly rather than silently producing different
truth. Any future task built on a rolling database needs the same treatment.
