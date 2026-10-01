# NPBench-C — build notes (Phase 0 / 0.5)

Benchmark for natural-product and BGC computational agents. 32 campaigns from
24 templates, four-rung ladders, deterministic grading, no LLM anywhere in the
scoring path.

## What exists

| Component | Path | State |
|---|---|---|
| Scoring primitives | `src/npbench_c/grading/primitives.py` | 11 primitives, stdlib-only |
| Composition policy | `src/npbench_c/grading/compose.py` | product / geometric+floor / min |
| Rung ladders | `src/npbench_c/grading/ladder.py` | depth, monotonicity, chance floor |
| Grader | `src/npbench_c/grading/grade.py` | pure declarative comparator |
| Grader identity | `src/npbench_c/grading/version.py` | semver + content hash |
| Readiness gate | `src/npbench_c/readiness/gate.py` | 16 mechanical checks, 5 pending-on-agent-runs |
| Campaign L1 | `campaigns/construct-ecoli-rebh-01/` | oracle 1.0, all 4 rungs G1 |
| Tests | `tests/unit/` | 55 passing |

```bash
PYTHONPATH=src python3 -m pytest tests -q
PYTHONPATH=src python3 -m npbench_c.readiness.gate campaigns/construct-ecoli-rebh-01
```

## The rung ladder

The graded unit is the rung, not the campaign. 32 campaigns x 4 rungs = 128
graded (prompt, gold) pairs.

| Rung | Measures | Target clear rate |
|---|---|---|
| R1 Execute | toolchain runs, output contract honoured | 85-95% |
| R2 Correct | primary artifact matches gold within measured tolerance | 55-70% |
| R3 Commit | derived claim that depends on R2 | 30-45% |
| R4 Counterfactual | answer under a perturbation that makes the published answer wrong | 10-25% |

Campaign score = highest **contiguous** rung cleared / rung count.
`rungs_cleared` is reported alongside so monotonicity violations surface.

Rungs are runnable **cold** (agent attempts everything) or **warm** (agent
starts from rung *k-1* gold and attempts only rung *k*). Warm mode is what makes
a rung a standalone pair and what yields a per-rung difficulty profile — in cold
mode a failure at R2 censors everything above it, so cold runs alone cannot
calibrate a ladder.

A 0.50 headline under this design means: SOTA agents execute NP pipelines
correctly but do not reliably convert artifacts into scientific commitments or
predict perturbation outcomes. That is a pre-registered, falsifiable claim about
where agents fail, not a tuning target.

### Ladder validity rules (enforced, not aspirational)

1. Monotonicity is **measured** against internal agent runs, never assumed.
2. R1 must clear at >= 80% for a generic bash-only agent, or the campaign is
   measuring environment setup rather than science.
3. R4 must be parametric-proof: its gold depends on a perturbation.
4. Each rung grades independently of later rungs.
5. `chance_floor` <= 0.10, computed and stored. This is the rule that most
   protects calibration — a four-rung ladder whose R3 is a three-way enum is
   not a hard campaign however hard R1 and R2 are.
6. Warm-start separability: rung *k* attemptable from rung *k-1* gold.
7. At most one G2 rung per campaign.

## Architecture: measure, then grade

```
submission/ --> measure.py --> measurement.json --> grade() --> grade.json
```

R1/R2 cannot compare the agent's sequence to the oracle's, because many
sequences satisfy the constraints — grading must check the invariant, not the
representation. So a campaign-specific deterministic **measurer** computes
properties from the submission, and the generic **grader** compares those to
gold and thresholds. The agent cannot self-report, and the grader stays
declarative. This generalises to every campaign.

## Purity

The grader is a pure function of (submission, gold, grader_version). Enforced:

- stdlib-only scoring path, asserted by an AST test over `grading/*.py`
- no network, no LLM, no wall-clock, no RNG
- all float reductions in fixed (sorted) order; every score rounded to 6 dp
- 100 repeat gradings byte-identical; inputs never mutated; key insertion order
  irrelevant
- `grader_version` = semver + content hash over the scoring path

Purity is not validity. A pure grader can be reproducibly wrong, which is what
the separate construct-validity study measures (see below).

## Campaign L1: construct-ecoli-rebh-01

Design a synthetic ORF for RebH (UniProt Q8KHZ8, SV 1, 530 aa), then analyse
which design constraint actually binds, and how that changes with host.

- **R1** ORF present and schema-valid
- **R2** protein identity exact, all six constraints satisfied, CAI >= floor
- **R3** unconstrained-optimum violation counts (6 exact integers), binding
  constraint, non-binding set — all for *E. coli*
- **R4** same under a *Streptomyces* host swap, plus the flipped-constraint set

All four rungs are G1. Measured discrimination on built submissions:

| Submission | Depth | Fails at |
|---|---|---|
| oracle | 4 | — |
| naive reverse translation | 1 | R2 |
| correct construct, no analysis | 2 | R3 |
| correct construct + R3, no R4 | 3 | R4 |
| everything right, assumes binding constraint is host-invariant | 3 | R4 |

The last row is the parametric trap: an agent that carries its R3 answer across
to the new host instead of recomputing is caught exactly at R4. The binding
constraint genuinely flips (gc_window -> direct_repeat), which is what makes R4
unguessable.

### Two fields deliberately excluded from grading

Spec 8.5 says tolerances are measured and an artifact whose spread exceeds what
the science tolerates is redesigned, never given a wider tolerance. That fired
here. Running four legitimate optimiser policies:

- **binding constraint: stable** across all four policies, both hosts → graded
- **`binding_cost`: spread 0.0356** for gc_window, wider than the
  between-constraint signal for four of six constraints → **not graded**
- **`cai_full`: spread 0.0528** across policies → graded as a **floor**
  (0.902 E. coli, 0.929 *Streptomyces*), set below the worst legitimate policy,
  so a different valid optimiser is never penalised

Both exclusions are recorded in `task.yaml` with the measurement that justified
them. A reviewer can check the reasoning without rerunning anything.

Note this is a **second tolerance class** the spec did not anticipate: the
determinism audit measures run-to-run spread, which for a deterministic oracle
is zero. What matters for a quantity an agent computes with its own method is
*cross-method* spread. Campaigns must measure both.

## Readiness gate

A campaign is ready when `npbench_c.readiness.gate` returns all-pass — not when
someone believes it is. Checks requiring internal agent runs report **PENDING**,
never PASS: the gate does not launder an unmeasured property into a green tick.

Current state for L1: **16 pass / 0 fail / 5 pending**.

The five pending checks (monotonicity, R1 clear rate, no-tool leakage,
difficulty gate, discrimination) need the internal sweep. They are satisfied by
dropping an `internal_sweep.json` into the campaign directory.

**Never audit before calibrating.** A campaign dropped for non-discrimination
after a domain expert spent three hours on it has burned the scarcest resource
in the project. Order is: build → internal sweep → drop/fix → audit → ship.
Author ~44 to ship 32.

## A note on what the gate caught

Two real defects surfaced during this build, both by machinery rather than by
reading the code:

1. **Grader float boundary.** `abs_tol_match(1.0, 1.1, 0.1)` returned 0.0,
   because `|1.0 - 1.1|` is 0.10000000000000009 in binary floating point. An
   answer sitting exactly on the stated tolerance was marked wrong — a brittle
   match contradicting the documented contract, and the kind of defect that
   reads to an auditor as grader noise. The difference is now rounded to the
   declared precision before comparison. Regression test in
   `tests/unit/test_primitives.py`.
2. **Stale `grader_version` pin.** `task.yaml` carried a hash from before the
   above fix. The gate checked that a version was *recorded*, not that it
   *matched* what just ran, so a manual comparison caught what the gate should
   have. `grader_version_matches_pin` now closes it.

Both are worth repeating because they are the failure class the benchmark is
built to detect in others: reproducible, confidently reported, and wrong.

## Open items

1. **CAI table provenance.** The relative-adaptiveness values in
   `campaigns/*/oracle/tables.py` are representative values for highly-expressed
   genes and are flagged in-file. They need a cited source release, or a table
   computed over a declared gene set, with a hash — before the audit packet
   ships. The gold does not depend on them being *the* correct measurement, only
   on their being frozen, hashed and published; but an auditor will and should
   ask.
2. **Release packaging must exclude `gold/`, `oracle/` and `oracle_submission/`.**
   They are committed here because this is the build repo. Gold is never
   published for either split.
3. Internal sweep harness (the five PENDING checks).
4. Construct-validity study: inter-rater agreement first, then expert-grader
   agreement with Gwet's AC1 / Krippendorff's alpha alongside kappa, gate on
   Spearman against the continuous rating. This is the one place humans are
   involved, and it sits outside the grading pipeline by design.
