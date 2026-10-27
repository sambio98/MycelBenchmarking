# NPBench-C — build notes (Phase 0 / 0.5)

Benchmark for natural-product and BGC computational agents. **25 campaigns from
20 templates**, four-rung ladders, deterministic grading, no LLM anywhere in the
scoring path.

The suite is frozen at 25: that is the owner's decision of 2026-10-27 and it
closes open item 7. The catalogue describes 32 campaigns from 24 templates; three
rows were retired on measurement with their numbers recorded below, and the two
remaining instantiations (T-L4-1 x2, T-L4-2 x2) are deliberately not built. 25
campaigns x 4 rungs = **100 graded rungs** and **560 graded components**.

## What exists

| Component | Path | State |
|---|---|---|
| Scoring primitives | `src/npbench_c/grading/primitives.py` | 12 primitives, stdlib-only |
| Composition policy | `src/npbench_c/grading/compose.py` | product / geometric+floor / min |
| Rung ladders | `src/npbench_c/grading/ladder.py` | depth, monotonicity, chance floor |
| Grader | `src/npbench_c/grading/grade.py` | pure declarative comparator |
| Grader identity | `src/npbench_c/grading/version.py` | semver + content hash |
| Readiness gate | `src/npbench_c/readiness/gate.py` | 27 mechanical checks, 5 pending-on-agent-runs |
| Template: construct design | `src/npbench_c/templates/construct_design/` | shared engine, 2 ladders |
| Campaign: construct design | `campaigns/construct-ecoli-rebh-01/` | oracle 1.0, gate 25/0/5 |
| Campaign: constraint conflict | `campaigns/constraint-gc62-ecoli-rebh-01/` | oracle 1.0, gate 25/0/5 |
| Campaign: feasible control | `campaigns/constraint-gc55-ecoli-rebh-01/` | oracle 1.0, gate 25/0/5 |
| Template: NRPS mass balance | `src/npbench_c/templates/mass_balance_nrps/` | shared oracle, 2 instantiations |
| Campaign: malleobactin | `campaigns/massbalance-nrps-malleobactin-01/` | oracle 1.0, gate 25/0/5 |
| Campaign: sevadicin | `campaigns/massbalance-nrps-sevadicin-01/` | oracle 1.0, gate 25/0/5 |
| Internal sweep | `src/npbench_c/sweep/` | runner, stats, gates, stub fixtures |
| Refusal contract | `tests/unit/test_build_refusals.py` | every one of 149 build refusals across 12 templates named by a test |
| Release packager | `src/npbench_c/release/package.py` | builds the shippable tree; refuses gold, the oracle, `docs/`, the grading spec or any graded answer; assembles attribution |
| Cold-sandbox scan | `src/npbench_c/readiness/cold_sandbox.py` | gate check 22 and a project-wide test on one implementation; declared collisions in `audit/`, outside every sandbox, checked for staleness |
| Catalog | `docs/catalog.md` | the design catalogue: 32 campaigns from 24 templates, grounded. 25 of them are built and the suite is frozen there |
| Provisioning | `docs/environment.md` | tools, DBs, sandbox contract |
| Template: MIBiG release diff | `src/npbench_c/templates/mibig_diff/` | shared oracle |
| Campaign: MIBiG 3.1→4.0 diff | `campaigns/mibig-diff-3_1-to-4_0-01/` | oracle 1.0, gate 25/0/5 |
| Template: annotation audit | `src/npbench_c/templates/mibig_annotation_audit/` | shared oracle |
| Campaign: PrnA construct | `campaigns/construct-ecoli-prna-01/` | oracle 1.0, gate 25/0/5 |
| Campaign: gene-function evidence | `campaigns/mibig-gene-function-evidence-01/` | oracle 1.0, gate 25/0/5 |
| Template: structure features | `src/npbench_c/templates/structure_features/` | shared engine, 2 ladders |
| Campaign: PrnA residue evidence | `campaigns/residues-prna-01/` | oracle 1.0, gate 25/0/5 |
| Campaign: RebH pocket geometry | `campaigns/pocket-rebh-01/` | oracle 1.0, gate 25/0/5 |
| Template: chemical space | `src/npbench_c/templates/chemical_space/` | shared oracle, RDKit pinned at instantiation |
| Campaign: MIBiG 4.0 chemical space | `campaigns/chemspace-mibig-4_0-01/` | oracle 1.0, gate 25/0/5 |
| Template: kinetics consistency | `src/npbench_c/templates/kinetics_consistency/` | shared oracle, BRENDA 2026.1 |
| Campaign: BRENDA EC 1.1 kinetics | `campaigns/kinetics-brenda-ec1_1-01/` | oracle 1.0, gate 25/0/5 |
| Template: ChEMBL selectivity | `src/npbench_c/templates/chembl_selectivity/` | shared oracle, ChEMBL_37, first `share_alike` |
| Campaign: S. aureus topoisomerases | `campaigns/selectivity-saureus-topoisomerase-01/` | oracle 1.0, gate 25/0/5 |
| Template: EC domain audit | `src/npbench_c/templates/ec_domain_audit/` | HMMER 3.4 + Pfam 35.0, first tool-running template |
| Campaign: EC 1.1.3.15 FMN_dh audit | `campaigns/ecaudit-fmn-dh-1_1_3_15-01/` | oracle 1.0, gate 25/0/5, **first S1** |
| Template: annotation transfer | `src/npbench_c/templates/annotation_transfer/` | DIAMOND 2.2.8, closed reference set |
| Campaign: EC 1.14 transfer | `campaigns/transfer-ec1_14-uniprot-01/` | oracle 1.0, gate 25/0/5 |
| Template: BGC detection | `src/npbench_c/templates/bgc_detection/` | antiSMASH 8.0.4, whole chromosome, **first S2** |
| Campaign: S. coelicolor regions | `campaigns/bgcdetect-scoelicolor-01/` | oracle 1.0, gate 25/0/5 |
| Template: GCF cutoff | `src/npbench_c/templates/gcf_cutoff/` | BiG-SCAPE 2.0.3, one run over a declared grid |
| Campaign: MIBiG 32 cutoff sweep | `campaigns/gcfcutoff-mibig-32-01/` | oracle 1.0, gate 25/0/5 |
| Template: RiPP precursor audit | `src/npbench_c/templates/ripp_precursor/` | HMMER 3.4 + declared Pfam panel, elected coordinate convention |
| Campaign: MIBiG RiPP precursors | `campaigns/ripp-precursor-mibig-4_0-01/` | oracle 1.0, gate 25/0/5 |
| Template: domain architecture | `src/npbench_c/templates/domain_architecture/` | HMMER 3.4 + declared 65-model Pfam panel, declared overlap resolution |
| Campaign: MIBiG NRPS architectures | `campaigns/arch-nrps-mibig-4_0-01/` | oracle 1.0, gate 25/0/5 |
| Campaign: MIBiG PKS architectures | `campaigns/arch-pks-mibig-4_0-01/` | oracle 1.0, gate 25/0/5, cis-AT vs trans-AT strata |
| Template: gene calling | `src/npbench_c/templates/gene_calling/` | Prodigal 2.6.3, declared fragmentation rule, no external resource |
| Campaign: S. coelicolor gene calls | `campaigns/genecall-scoelicolor-01/` | oracle 1.0, gate 25/0/5 |
| Template: remote homology | `src/npbench_c/templates/remote_homology/` | MMseqs2 18.8cc5c over one superfamily, InterPro 110.0 as the label, declared shuffle ladder |
| Campaign: IPR036188 twilight zone | `campaigns/homology-twilight-ipr036188-01/` | oracle 1.0, gate 25/0/5, recovery from 5% to 97% across five queries |
| Template: frameshift response | `src/npbench_c/templates/frameshift_response/` | Prodigal 2.6.3 in the opposite mode, planted single-base deletion, coordinate lift past the cut |
| Campaign: S. coelicolor frameshifts | `campaigns/frameshift-scoelicolor-01/` | oracle 1.0, gate 25/0/5, ladder samples the measured turnover band |
| External resources | `src/npbench_c/tools/resources.py` | one declaration for a pinned file read from the image, shared by 3 campaigns |
| Pinned-tool image | `image/Dockerfile`, `image/environment.lock.json` | 9 tools, 226-package closure, lock generated from the prefix |
| Tool memo | `src/npbench_c/tools/cache.py` | stub sweeps only; 29m32s to 3m44s |
| One-pass archive reads | `src/npbench_c/tools/archive.py` | a .tar.gz read in order; ripp's sweep 7m51s to 45s |
| Phase 1 image | `image/Dockerfile`, `image/environment.lock.json` | 8 tools pinned to `version=build`, 185-package closure hashed, antiSMASH databases pinned at 9.4 GB |
| Tool registry | `src/npbench_c/tools/registry.py` | pins, controls, invocations, declared normalisations, canonicalisations, projections, enforced bans |
| Thread-invariance suite | `src/npbench_c/tools/invariance.py` | 8/8 tools invariant at 1 and 8 threads |
| Tests | `tests/` | 685 passing, 6 skipped |

```bash
PYTHONPATH=src python3 -m pytest tests -q
PYTHONPATH=src python3 -m npbench_c.sweep.cli campaigns/construct-ecoli-rebh-01 --stubs
PYTHONPATH=src python3 -m npbench_c.readiness.gate campaigns/construct-ecoli-rebh-01
```

## The rung ladder

The graded unit is the rung, not the campaign. 25 campaigns x 4 rungs = **100
graded (prompt, gold) pairs**, carrying 560 graded components in all -- a median
of 23 per campaign, from 11 to 38.

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
- **`binding_cost`**: the cross-optimiser spread is wider than the gap from that
  constraint's mean cost to the nearest other constraint's, for 4 of 6
  constraints in RebH's E. coli host (widest 0.035363, on gc_window — the
  constraint that binds) and 1 of 6 in its *Streptomyces* host (0.084266, on
  direct_repeat, which also binds). PrnA: 2 of 6 and 1 of 6. → **not graded**
- **`cai_full`**: the four policies disagree by 0.068207 / 0.084266 (RebH, two
  hosts) and 0.072642 / 0.108103 (PrnA), all wider than the 0.025 floor margin →
  graded as a **floor** (RebH 0.888 E. coli, 0.880 *Streptomyces*; PrnA 0.883
  and 0.848), each set that margin below the worst legitimate policy, so a
  different valid optimiser is never penalised

The spreads above are wider than the ones this section first reported, because
the codon tables they are measured over were re-derived from pinned genome
records rather than asserted — see **Where the codon tables come from** below.

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

Current state, measured over all 25 campaigns on 2026-10-27, one full gate run
each: **27 pass / 0 fail / 5 pending, every campaign, no exceptions.**

The five pending checks (monotonicity, R1 clear rate, no-tool leakage,
difficulty gate, discrimination) are satisfied by an `internal_sweep.json` from
a sweep over **real** systems. A stub-based sweep leaves them pending by design.

The twentieth check, `external_resources_pinned_and_present`, arrived with the
first S1 campaign. A campaign may need reference data it cannot ship — too large,
or under a licence the benchmark cannot pass on — and such a resource is declared
instead of copied: path, provider, version, sha256. The check resolves each one
and verifies the fingerprint, because an unverifiable external resource is the
solvability defect one step out: the sandbox looks complete and the answer depends
on bytes nobody pinned.

The twenty-seventh check, `contract_names_what_is_graded`, closes the gap
between what a campaign publishes and what it grades. The oracle never goes
through the published contract: `solve.py` renders a submission from gold's own
key tuples, so a contract that named a field the measurer does not read — or
failed to name one it does — still scored 1.0 while being unusable by any system
that followed it. The check projects the oracle's own submission onto the field
names `output_contract` declares and requires the projection to grade exactly as
the oracle does.

It was written because that projection, run for the first time over the frozen
suite, failed on **5 of 25 campaigns**:

| Campaign | Depth from the published contract | Cause |
| --- | --- | --- |
| `arch-nrps-mibig-4_0-01` | 0 | field names never updated when the template was generalised for PKS: `a_domain_verdicts`, `c_domain_verdicts`, `a_delta_census`, `module_type` against the measurer's `primary_verdicts`, `secondary_verdicts`, `primary_delta_census`, `primary_domain`/`secondary_domain`/`module_family`, plus three graded fields (`strata_census`, `domains_by_module_type`, `secondary_verdicts_by_stratum`) never declared at all |
| `arch-pks-mibig-4_0-01` | 0 | the same contract, the same way |
| `construct-ecoli-prna-01` | 2 | R3 and R4 answers live under `r3`/`r4` keys the contract never mentioned, so every analytical field was declared one level too high |
| `construct-ecoli-rebh-01` | 2 | the same |
| `residues-prna-01` | 0 | `mutagenesis_original_residue_agreements / _disagreements` — the second name abbreviated to a suffix no system can write |

Two of those scored **zero at R1**: a perfect system that read the contract and
reported exactly what it asked for would have been graded as having failed to
execute. Four more campaigns wrote a field nothing declared, or declared one
nothing read, without changing any score — cosmetic, and now fixed and checked
for too.

Contract keys that are not fields — a format, a closed vocabulary, a property of
a FASTA record — live under a reserved `conventions` key, which is how the check
tells a declaration about an artefact from a field of it, and how the prompt
generator knows to render it as a convention rather than as something to report.
The check reaches as deep as the contract names keys: where a contract stops
naming keys and describes the shape below in prose (a table's columns, a map
keyed by an identifier), a misnamed column is not something it can see.

**Never audit before calibrating.** A campaign dropped for non-discrimination
after a domain expert spent three hours on it has burned the scarcest resource
in the project. Order is: build → internal sweep → drop/fix → audit → ship.
The catalogue budgeted ~44 authored to ship 32; what it took to ship 25 was 28
designs, of which 3 ground out on measurement and 2 x2 instantiations were left
unbuilt by the freeze.

## Campaign solvability

**A campaign must be solvable from its own sandbox.** Every table, vocabulary,
counting rule and convention its gold depends on must be in `reference/` or
`task.yaml`.

This campaign shipped briefly with its codon weight table in `oracle/`, which
the sandbox excludes. R3's gold is a pure function of those weights, so the
campaign was unanswerable — and no run would have reported it: an agent would
simply have failed and looked incapable. The stub systems hid it because they
were handed `--oracle` on the command line.

Guards: `pinned_tables_reachable_by_agent` in the readiness gate, and a test
that recomputes R3's gold from sandbox files alone. The oracle loads the same
`reference/` JSON the agent receives, so the two cannot drift — verified by
regenerating gold after the extraction and diffing byte-for-byte.

See `docs/environment.md` for the full sandbox contract, the tool surfaces, and
the per-template database and tool requirements.

## A note on what the gate caught

Three real defects surfaced during this build, all by machinery rather than by
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

3. **An unsolvable campaign.** Covered above: the pinned tables the task
   requires were outside every path the sandbox ships.

All three are worth repeating because they are the failure class the benchmark
is built to detect in others: reproducible, confidently reported, and wrong.

## The internal sweep

Turns the readiness gate's five PENDING checks into numbers. Agent-agnostic: a
system is a `SystemSpec` (a command, a mode, a timeout), so a bash-only coding
agent, a general biomedical agent and a tool-equipped system all plug in without
the harness knowing anything about them.

```
cold runs   agent gets the objective, attempts everything  -> depth distribution
warm runs   agent starts from rung k-1 gold, attempts k     -> per-rung profile
```

Both are needed. Cold runs alone cannot calibrate a ladder: a failure at R2
censors every rung above it, so cold R3/R4 clear rates conflate "could not do
it" with "never got there".

### Gold isolation is the load-bearing property

If gold reaches a sandbox, every number the sweep produces is silently invalid
and nothing downstream notices — the runs simply look successful. So sandboxes
are built from an allowlist (`inputs/`, `task.yaml`) and then **audited**; a
violation raises.

Warm starts are **pre-rendered bundles**, never whole gold files. Declaring
`gold/gold.json#hosts.ecoli` as an R4 warm input and copying the file would hand
the run the *Streptomyces* answer it is being asked to compute. `gold/warm/r4/`
contains exactly the E. coli construct and the E. coli analysis, nothing more,
and a test asserts the *Streptomyces* gold values are absent.

### Statistics at n=3

Three repeats cannot support a claim if the unit of analysis is the campaign.
With a ladder the unit is the rung-run — 25 campaigns x 4 rungs x 3 repeats =
300 observations — and CIs come from a **cluster bootstrap resampling
campaigns**, which absorbs correlation between rungs of one campaign and
between repeats of one run. Once campaigns are instantiated from templates,
resample at the **template** level: four instantiations share an oracle and a
failure mode and are not four independent campaigns.

Pre-register the resulting resolution: **NPBench-C separates systems differing
by roughly 8 points or more at n=3, and makes no claim about finer ordering.**
Pass@1 and Pass^3 are both reported; the gap between them is the only variance
story available without paying for more repeats.

Infra failures (timeout, crash) are excluded from scores and reported
separately in `run_status_counts`, so the headline is not partly measuring our
own container.

### Two chance floors

The sweep exposed a contradiction in the original design. R1 is *meant* to be
easy — its gate is a >= 80% clear rate — and for a construct-design campaign a
no-tool system clears it by reverse-translating the protein, which needs the
genetic code and nothing else. The ablation measured R1's no-tool clear rate at
1.0 against a declared estimate of 0.30. So any four-rung ladder gives a no-tool
agent 0.25 for free, and a 0.10 cap on the **full** floor is unsatisfiable by
construction.

Resolved by splitting them:

- `chance_floor` — all rungs (~0.26 here). The baseline the **leakage ablation**
  is compared against.
- `discriminating_chance_floor` — R2 upward (~0.007 here). What the **0.10 cap**
  applies to, since R1 is an execution precondition, not a measurement.

A test asserts the split still catches a guessable R3, so excluding R1 has not
widened the hole the cap was for.

### Stub fixtures, and why they cannot turn a gate green

`src/npbench_c/sweep/stubs/` holds deterministic stub systems that reach preset
depths. They exist so the harness's arithmetic and gate logic can be validated
before any compute is spent — a sweep whose statistics are wrong is worse than
no sweep, because its numbers look authoritative.

A stub-based sweep therefore marks itself `stub_based: true`, and the readiness
gate keeps all five checks **PENDING** with the harness result shown but not
promoted. `ready_to_audit` stays `False` until real systems run. With the stub
sweep on record the gate reports **27 pass / 0 fail / 5 pending** — the count
that line carried when it was written was 16, which was the number of mechanical
checks at the time.

Harness self-test result (fixtures, not evidence): clear rates 1.00 / 0.75 /
0.50 / 0.25 across R1-R4, cold depths 1/2/3/4, all five gates green.

## Campaign 2: NRPS mass balance

`massbalance-nrps-malleobactin-01` — MIBiG BGC0000386, a four-module NRPS line
in *Burkholderia thailandensis* with four annotated congeners including a dimer.
Compute the assembly the modules would produce, quantify the gap to each
congener, predict every single-module deletion. All four rungs G1, solvable at
tool surface S0, stub discrimination 1/2/3/4, all five sweep gates green.

Three design decisions worth carrying to the other mass-balance instantiations:

- **Every annotated congener is analysed, not `compounds[0]`.** 17 of 43
  admissible entries annotate more than one compound, so "the" product is not
  well defined and picking the first listed would bake an undeclared convention
  into gold and mark a correct answer about another congener wrong.
- **Monomer formulas are resolved once with RDKit and pinned** into
  `reference/monomer_formulas.json` with their SMILES source. The oracle needs
  no chemistry toolkit at run time, so a toolkit version bump cannot move a
  score and the campaign stays S0-solvable.
- **The residual is quantified, never attributed.** Naming the chemistry that
  closes the gap is expert inference; its gold would be judgement. Excluded with
  the reason recorded in `task.yaml`.

Also: this campaign's stubs read **only the sandbox**, with no privileged oracle
path, precisely because the construct-design stubs' oracle access hid an
unsolvable campaign.

## Abstention gold that is provable

The constraint-conflict campaign is the first where the correct answer is to
**produce nothing** and say why. That makes its gold the most dangerous in the
catalog: if the constraints are merely hard rather than impossible, an agent
with a better optimiser finds a design and we mark it wrong.

So infeasibility is established by arithmetic, never by search failure. For each
residue independently, the highest-GC admissible codon is determined; because
codon choices are independent, the resulting whole-ORF bound is **attained**, not
merely valid — a test constructs the attaining sequence and confirms it encodes
the protein. A requirement above that bound is unreachable for every sequence,
whatever other constraints apply, because adding constraints only shrinks the
feasible set.

```
unrestricted achievable GC max        0.677338
after forbidding the 8 GC3 codons     0.584432
required minimum                      0.620000   -> impossible
```

Each constraint alone is satisfiable, so the conflict is genuinely pairwise, and
the build **refuses to emit gold** otherwise: it raises if the requirement is
still reachable, and separately if it exceeds even the unrestricted maximum.
Both refusals are tested.

Fabricating a design is a graded component, not advice: `no_design_claimed`
must be 0, which is the Output-Fabrication failure mode made measurable.

**The control, now built.** Abstention scored on its own rewards reflexive
abstention — a system that always answers "infeasible" would score full marks on
the conflict campaign. `constraint-gc55-ecoli-rebh-01` differs from it by
**one number** (GC floor 0.55 rather than 0.62); everything else is identical.
Precision and recall are reported separately across the pair, never F1 alone.

Both halves of the property are tested **at the grader**, not inferred from the
sweep: a submission with every bound correct but a reflexive verdict fails R3 at
depth 2 on each campaign. That check matters because both ablation stubs happen
to die at R2, so stub depths alone would have told us nothing about R3.

**The asymmetry.** Infeasibility is decidable by arithmetic. Feasibility of the
full set is only provable **constructively** — the achievable-GC bound settles
the GC axis, but the composition constraints could still exclude every
candidate. So the feasible half rests on an exhibited design, verified to have
zero violations and to encode the declared protein (GC 0.569994), recorded in
gold for the audit packet and **not graded**: many designs qualify, so grading
one would grade our optimiser rather than the agent's reasoning.

Building that witness needed a new engine mode. The greedy repair used by the
design ladder cannot climb a *global* GC floor — a single synonymous
substitution rarely crosses the threshold, so the violation count never strictly
improves and the search stalls. `design_under_conflict_set` instead starts from
the GC-maximising admissible assignment and repairs composition constraints
while refusing any substitution that would drop GC back below the floor, with
the GC delta computed in O(1) and a bounded rotating candidate window. First
attempt was quadratic in ORF length and too slow for the gate's determinism
re-run; second was too narrow and stalled on 50-nt GC windows.

**The build refuses to contradict its declared expectation**, in either
direction, so a campaign cannot silently become the opposite of what it was
authored to be. Tested both ways.

## A leakage ablation is not a competence probe

The MIBiG campaign made a distinction explicit that the earlier ones let slide.
Its `stub-nonormalise` loads both corpora, computes the raw diff correctly, and
reports it as the classified result — clearing R2 and failing R3. Labelling that
the no-tool ablation made the leakage gate **fail**, and the gate was right to
complain: it asks "is this answerable without doing the work?", and that system
did the work.

So the two are now separate roles:

- **competence probe** — a tooled system missing one skill. Useful for
  discrimination, meaningless for leakage.
- **leakage ablation** — withholds *computation itself*, which for an S0
  campaign is the only thing there is to withhold. `stub-noncompute` never opens
  the archives and answers with round numbers, so it fails R1 at depth 0.

This also corrected a claim carried in the earlier campaigns' chance levels. In
those, R1 is free — reverse-translating a protein or echoing a table needs no
computation — and the measured no-computation clear rate was 1.0. Here it is
**not**: the entry counts cannot be guessed, so the ablation fails R1. R1's
declared chance level stays 1.0 because any system that reads the inputs clears
it, and the capped floor is computed over R2 upward regardless, but the note in
`task.yaml` records that it is nominal here rather than measured.

## Template architecture

A campaign of a factored template contains **no code**: only `inputs/`,
`reference/`, `task.yaml`, `grading.json` and `gold/`. The oracle is shared, and
an instantiation is a data change plus two commands:

```bash
python3 -m npbench_c.templates.mass_balance_nrps.build_reference campaigns/<id>  # RDKit, once
python3 -m npbench_c.templates.mass_balance_nrps.build           campaigns/<id>  # stdlib, re-runnable
```

Reference generation is split from gold generation deliberately: the readiness
gate re-runs gold to check determinism, and it must be able to do that without a
chemistry toolkit and without silently re-deriving the pinned monomer table.
`grading.json` is **generated from gold**, so the spec cannot drift from the data
it grades.

`task.yaml` declares the three harness commands (`gold_command`, `solve_command`,
`measure_command`) and the stub fixtures, so neither the gate nor the sweep knows
anything campaign-specific. The construct-design campaign keeps a local oracle
and declares it the same way, so there is one code path rather than a fallback.

Both template families now work this way. `construct_design` holds one engine
and two ladders — `design/` (T-L6-1) and `conflict/` (T-L6-2) — which share the
constraint arithmetic and must not diverge. Migrating the construct campaign
onto it reproduced its gold **byte-identically**, which is how the port was
verified rather than assumed.

The payoff is measurable: the sevadicin instantiation is 2 JSON files and a
`task.yaml`. Copying six oracle files per campaign would have made drift between
instantiations invisible, and the catalog expects eight doubled templates.

## Catalogued designs corrected by grounding

Each was checked against real data before any oracle was written, and each would
have produced a plausible-looking campaign that could not be built or whose gold
would have been wrong.

1. **PKS mass balance** — `at_domain.substrates` is empty throughout MIBiG, so
   extender-unit identity is absent. The template is NRPS-only.
2. **"Validate the route"** — MIBiG annotates the enzymatic module count, not the
   number of monomer incorporations; only 1 of 43 admissible entries is exactly
   balanced. The campaign quantifies the gap instead.
3. **"What changed and in which field"** between MIBiG releases — the 4.0
   re-annotation restructured the schema, so a raw field diff is true and
   useless. The campaign reconciles on a declared mapping instead, and my own
   first pass reported 2442 taxon re-annotations where there are 8.
4. **Self-resistance target identification** — the resistance gene is a field
   lookup and the target is not in MIBiG; inferring it needs homology search.
   Re-tiered to `build: image`, with the S0-feasible evidence audit built in its
   place.

5. **Chemical space / scaffold analysis** — Bemis-Murcko reduction is far coarser
   than the design assumed: the most shared scaffold over the corpus is plain
   benzene, in 76 entries. The catalogued "distinguishing substitution" rung is
   therefore both near-vacuous and prose-only, and the catalogued substituent-swap
   counterfactual needs a chemistry toolkit at solve time. R3 grades the
   reconciliation between the naive and filtered sharing statistics instead, and
   R4 varies the declared filter and the evidence gate.

6. **Selectivity counterfactual (T-L5-1)** — the catalogued R4 excluded a
   low-confidence assay. Every one of the 105 assays behind this target pair's IC50
   values scores 7, so the exclusion rejects nothing. The rule stays in the
   contract; R4 varies the axes that move instead.
7. **Resistance mutation → target (T-L5-2)** — the structured
   `assay_variant_mutation` field is populated across 34,921 IC50 activities and
   carries almost nothing on theme: 45 variant records on the *S. aureus* gyrase
   complex, 0 on topoisomerase IV, and all 244 *P. falciparum* DHFR records reading
   "UNDEFINED MUTATION". Re-tiered rather than built on an anecdote.

8. **A-domain substrate specificity (T-L2-4)** — the position of the A-domain in
   its protein is `-1/-1` for 2,852 of 2,901 domains, the substrate labels are
   71.6% inference, and the cross-resource comparison gives 0.000, 0.151 or 0.995
   agreement depending on which key is used. Retired, not built.

The pattern is consistent enough to be a rule: a catalogue entry is a hypothesis
about data until the fields are inspected.

The fifth case sharpens it, because the fields *were* inspected and the design
still needed correcting. Scaffold SMILES parse, scaffolds are computable, and
every count in the re-grounding probe was real — what the probe did not show was
that the computed scaffolds are mostly bare rings. So the rule has a second half:
**a field being present does not make the quantity derived from it informative,
and only computing the distribution shows which.** The campaign's own
`notes_for_audit` names the benzene result first for that reason.

Cases 6 and 7 are that second half twice more, and the pair separates its two
outcomes. For T-L5-1 the uninformative field cost a counterfactual, which was
redesigned around axes that move; for T-L5-2 it cost the campaign. **The rung that
survives a vacuous field is the one whose claim can be re-pointed at something
else; the row that does not is the one whose claim *was* the field.**


## Campaigns 9 and 10: structure features, two ladders on one engine

`residues-prna-01` and `pocket-rebh-01` are the two ladders of
`structure_features`, the third factored template. They exist as one template
because they read the same two sources for the same proteins — UniProt's
per-residue features with their ECO evidence codes, and the mmCIF coordinates
those codes cite — and building them apart would have duplicated an mmCIF
parser, a numbering mapper and an ECO classifier, which would then have drifted.

The proteins are the two already in the benchmark. RebH and PrnA arrived as the
two construct-design instantiations, so no new provenance is introduced; what is
new is that both carry evidence-coded residue annotations and seven to ten
deposited structures each.

### What makes the annotation checkable at all

UniProt does not merely assert that position 348 of PrnA binds chloride. Under
`ECO:0007744` it names the deposited entry the assertion rests on:

```json
{"type": "Binding site", "location": {"start": {"value": 348}},
 "ligand": {"name": "chloride", "id": "ChEBI:CHEBI:17996"},
 "evidences": [{"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "16195462"},
               {"evidenceCode": "ECO:0007744", "source": "PDB", "id": "2AQJ"}]}
```

That field turns "is this annotation supported?" into arithmetic: load 2AQJ,
find the chloride, compute its contact shell, check whether 348 is in it. No
curator is consulted and no judgement is exercised, which is the same doctrine
the MIBiG evidence audit applies — turned on a second database, from the
structural side.

### What the two ladders found

Both are real results, and both shaped the rungs rather than being discovered
afterwards.

**The annotation is a strict subset of the geometry, in one direction only.**
Over 2AQJ, every annotated position is inside its ligand's contact shell — 2 of
2 for chloride, 10 of 10 for FAD, 5 of 5 for tryptophan — and this holds at
4.0 Å and 3.5 Å as well, so it is not an artefact of a generous cutoff. The
disagreement runs entirely the other way: 29 of FAD's 39 contacts, 9 of
tryptophan's 14 and 3 of chloride's 5 are not annotated. So `unconfirmed` is
zero on both structures, which means **R3's discriminating content is the contact
set and the unannotated remainder, not the confirmation verdict** — a system
could guess "all confirmed" and still fail R3 on the sets. The campaign's audit
notes say so plainly rather than letting an auditor discover it.

**Truncating a pocket splits it almost evenly.** Deleting each of the 14
residues lining RebH's substrate site back to alanine, one at a time, removes
the contact for exactly 7 of them and leaves 7 still touching through backbone
or CB. Neither verdict is the safe guess, which is what makes the rung worth
grading; the build refuses to emit gold if the split collapses.

### Three design decisions worth recording

**The pocket volume is a declared lattice count, not a published descriptor.**
Every published pocket-volume definition carries parameter choices, and a
benchmark that adopts one grades those choices. So `geometry_rules.json` pins a
1 Å lattice anchored on integers, a 1.4 Å probe, Bondi radii and three stated
admission conditions, and the graded quantity is the admitted point count — 57
for this site. Solvent-accessible surface area and any druggability score are
excluded for the same reason, with the reason recorded in the campaign.

**No toolkit, on purpose.** Neither ladder uses Biopython, gemmi, RDKit or any
alignment tool; the parser and the geometry are stdlib, and both campaigns
declare `min_tool_surface: S0`. The moment pocket geometry runs through a
toolkit, the toolkit's defaults for hydrogens, alternate locations, symmetry
expansion and radii become part of the answer.

**R4's perturbation scope is declared by rule, never by list.** The truncation
sweep covers the focal contact shell, which is R3's answer. Naming those
positions in `reference/campaign_keys.json` would have handed it over, so the
sandbox carries `truncation_scope: focal_shell` and a test asserts that no number
in the declared keys is a shell position. This is the same gold-isolation
reasoning that replaced whole gold files with pre-rendered warm bundles.

### What the gates caught this time

Three things, all before any agent ran.

1. **An mmCIF category written two ways.** `_struct_ref_seq` appears as a `loop_`
   in 2E4G and as bare `_category.item value` lines in 2AQJ, because it has two
   rows in one and one in the other. A parser handling only the loop form returns
   nothing for the second — and the numbering offset then gets *assumed* rather
   than read, which produces a shell that is internally consistent and entirely
   wrong. The parser handles both and a test pins both.

2. **An R2 component that depended on an R3 answer.** The graded nearest-contact
   distance was being looked up by the copy the *submission* nominated as focal,
   which lives in its R3 block. The stub that computed the shells correctly and
   stopped there therefore failed R2, and the difficulty gate failed with it. The
   measurer now identifies the focal copy from the campaign's own declared keys.
   Worth stating as a rule: **a rung's components must be answerable from that
   rung's work alone**, and the sweep is what surfaces a violation.

3. **Two components that restated an earlier rung.** The pocket ladder graded the
   focal shell in R3 when R2 already graded it per copy, and the residues ladder
   graded the alternate structure's numbering offsets in R4 when both entries of
   the protein carry the same offsets as the focal one. Both were caught by a
   test asserting that no warm bundle contains the answer to the rung it starts
   — the bundle legitimately carried them, because they were the *previous*
   rung's deliverable. Both components were dropped and the exclusions recorded
   in `excluded_from_grading`.

That third check was worth keeping, so it is now the gate's eighteenth
mechanical check: **if a rung's warm-start bundle contains the answer to one of
that rung's components, the component belongs to an earlier rung.** It is a
mechanical test for a design error that is otherwise invisible until a real
system clears a rung it should not have — the oracle scores 1.0 either way.

Only **composite** gold values are checked. A scalar drawn from a declared
vocabulary legitimately appears all over an earlier rung's census: the construct
ladder's R4 asks which constraint binds for the swapped host, gold is
`"direct_repeat"`, and the R4 bundle contains that string as a key of the first
host's violation counts. Flagging it would be a false positive there and in every
other enum-scored component in the benchmark. A test pins the exemption so it is
not "tightened" later.

Run across the whole benchmark, it caught one more thing: **the mass-balance
ladders were grading the compound list at R3, and every warm bundle handed it
over.** They had to — neither the assembly rung nor the reconciliation rung can
be posed without knowing which compounds the entry names. The list is also just a
field of the shipped MIBiG entry, so the component was measuring a parse the
sandbox had already answered. It is gone from both campaigns, with the reason in
their `excluded_from_grading`.

It is worth being precise about what this check is and is not. It does not detect
gold leaking into a sandbox — `_audit_sandbox` and `GoldLeak` do that, and they
catch a file in the wrong place. This catches something subtler and entirely of
the author's making: a rung graded on work that belongs to the rung below it. The
three cases it found were all mine, and all of them had passed every other gate.


## Campaign 11: chemical space, and the question with no single answer

`chemspace-mibig-4_0-01` asks how many distinct natural products MIBiG 4.0
contains. There is no single answer, and that is the campaign:

```
compound records                     3,449
distinct full InChIKeys              3,115     334 records restate a structure
distinct InChIKey block 1            2,996     119 further distinctions merged
distinct Bemis-Murcko scaffolds      1,843     plus 182 compounds with no ring
```

Three declared identity keys over the same records, three different numbers, and
every statistic downstream — duplicates, per-class spread, scaffold sharing —
inherits the choice. The connectivity block is not a rounding of the full key: it
merges 111 groups of stereoisomers, charge states and labelled variants, and
where stereochemistry is not determinable it is the honest granularity. The rung
that distinguishes them is the campaign's spine.

Chemistry is perceived once with RDKit at instantiation and pinned into
`reference/chemistry_table.json`, the same move that put monomer formulas in
mass balance's reference directory. The table supplies what needs a toolkit —
formula, InChIKey, scaffold, atom and ring counts, per record — and nothing
aggregated. Everything the ladder grades is bookkeeping over it, so the campaign
is **S0** and a toolkit upgrade cannot move a published score.

### Bemis-Murcko is much coarser than the catalogue assumed

This is the fourth catalogued design that grounding corrected, and the first
where the correction only appeared *after* the oracle ran:

```
most shared scaffold, unfiltered       benzene            76 entries
second                                 tetrahydropyran    26 entries
```

Bemis-Murcko keeps ring systems and the linkers between them and discards
everything else, so any natural product carrying one aromatic ring reduces to
benzene. "These 76 compounds share a scaffold" is true and says nothing. The
catalogued R3 asked for *the distinguishing substitution* between two compounds
sharing a scaffold — which on this data is almost the whole molecule, and is in
any case a set difference between molecular graphs that can only be named in
prose. So R3 grades the reconciliation instead, under a declared filter of ≥2
rings and ≥50% heavy-atom coverage:

```
shared scaffolds           279  ->  190 informative
cross-class scaffolds       36  ->   14 informative
```

Both sides are graded. A system that reports only the unfiltered number is not
wrong, it is incomplete, and R3 is where that shows — which is what the
`stub-unfiltered` competence probe exists to confirm: it does the whole
bookkeeping, reports the unfiltered statistic in the informative slots, clears R2
and fails R3.

The catalogued R4 — predict the scaffold after a substituent change — would need
scaffold perception at solve time and would make the campaign S1 for the sake of
one rung. R4 varies the declared thresholds instead (four settings, record counts
3,267 / 2,704 / 2,497 / 1,954) and then gates the corpus on experimental locus
evidence, which cuts it to 601 entries and the cross-class informative set from 14
scaffolds to 3.

### One counting rule the data forced

A scaffold "crosses biosynthetic classes" only when **two of its entries have
disjoint class sets**. The obvious rule — the union of its entries' classes has
more than one name — gives 115 where the declared rule gives 36, because 456
active entries carry more than one class and a single hybrid PKS/NRPS cluster
therefore crosses classes on its own, with nothing to compare against. Both
numbers are reported, because the naive one is what most analyses compute and the
gap is the point.

### Two smaller things worth keeping

**The record id is `(accession, ordinal)`, not `(accession, name)`.** Compound
names are not unique within an entry: BGC0002024 lists "nargenicin A1" twice and
BGC0002072 lists "linearmycin C" twice. Keying on the name drops two records and
shifts every count below it, silently.

**`scaffold_smiles` ships but is never a key.** The existing rule from the
mass-balance build — RDKit for formula, never canonical SMILES, InChIKey as the
comparison key — applies here too, so scaffold identity is the scaffold's
InChIKey. The two do not agree exactly: 1,843 InChIKeys against 1,850 canonical
SMILES, because InChI normalises tautomers and charge states the SMILES writer
keeps apart. That disagreement is the argument for declaring the key rather than
leaving it to the solver, and it is why both numbers appear in a test.

### A correction to the re-grounding report

The re-grounding pass recorded 1,655 distinct Murcko scaffolds, and that number
does not reproduce under any definition — 1,843 by scaffold InChIKey over active
entries, 1,850 by canonical SMILES, 2,254 including retired entries. It also
called the 334-record gap "cross-entry duplicates", which it is not: most of
those records repeat a structure inside one entry, and the cross-entry figure is
217 at the full key and 259 at the connectivity block. `docs/regrounding-2026-10-02.md`
carries the correction. The lesson is narrower than the earlier one about
catalogue entries being hypotheses: **a probe's summary statistic is a hypothesis
until the oracle recomputes it**, because a probe has no gate behind it.


## Phase 1: the pinned-tool image, and the two controls that were wrong

The provisioning notes set the rule: every tool runs at 1 and 8 threads on fixed
input and must produce identical output, and anything failing which cannot be
pinned into determinism is made single-threaded or kept off gold-producing paths.
`image/` now holds the image and `src/npbench_c/tools/` the suite that enforces
the rule.

One honest caveat first: **the Dockerfile has not been built end to end**, because
the host had the Docker CLI and no daemon. Every step in it was run directly —
the micromamba install of these exact specs, `registry --verify` against the
resulting prefix, and the full invariance suite — so the pins and the measured
determinism results are real, and the layer sequence is not yet proven.

Eight tools, each pinned to an exact `version=build` string with the resolved
185-package closure and every artifact's sha256 in `image/environment.lock.json`:
HMMER 3.4, DIAMOND 2.2.8, MMseqs2 18.8cc5c, Prodigal 2.6.3, MAFFT 7.526, BLAST+
2.17.0, antiSMASH 8.0.4 and BiG-SCAPE 2.0.3. **Ubuntu's archive cannot satisfy
two of the pins** —
it ships DIAMOND 2.1.9 against a ≥ 2.2.7 floor and BLAST+ 2.12.0 against 2.17.0 —
which settled the question of whether `apt` would do. A measured dry-run then
confirmed that adding antiSMASH moves **none** of the six sequence-tool pins, which
is why there is one environment rather than two.

**This unblocks 10 of the 12 `build: image` templates**: T-L1-1, T-L1-2, T-L1-3,
T-L1-4, T-L2-1, T-L2-3, T-L2-4, T-L3-1, T-L3-2 and T-L3-4, plus the catalogued
T-L3-5 target inference that the S0 build re-tiered to `image`. The remaining two
need matchms, deferred with a recorded reason rather than left as a gap in a
table.

### What measuring found

Four tools are invariant with nothing asked of them: BLAST+, MAFFT, Prodigal and —
after one declared normalisation — HMMER, which stamps its own command line, its
working directory and the wall-clock date into both `--tblout` and the HMM file.
The other two were not, and both answers were worth the trouble of getting.

**DIAMOND: the control this project had declared was backwards.** The provisioning
table listed `--no-reorder` as a required determinism control. Measured, *with*
the flag 2.2.8 emits the same 192 hits in three different query orders across four
multithreaded runs; *without* it every run at 1, 4 and 8 threads is byte-identical.
DIAMOND's default restores query order and `--no-reorder` documents itself as
switching that off for speed — so the flag was the cause, not the cure. It is now
forbidden on any gold-producing path. A standing **witness invocation** keeps the
failure under test and is reported as `xfail`, so a future release that fixes it
shows up as an `XPASS` to re-read rather than as a silently changed assumption.

**MMseqs2: no flag fixes it, so the contract does.** At one thread every run is
byte-identical; at eight the hit *multiset* is identical and only its order moves.
Nothing in the tool controls that, so the tabular output is sorted by line before
anything reads it.

### antiSMASH: the pin is a pair, and three things followed from adding it

**Half the pin is the databases.** The binary is 8.0.4; what decides which regions
come out is that plus the reference data `download-antismash-databases` fetched —
knownclusterblast 4.0, Pfam 35.0, MITE 1.3, as-js 0.16 and six more directories,
**9.4 GB on disk**, now recorded in the lock with their versions. A pin that names
only the binary leaves out the half that moves the answer.

On top of the version the lock records a **fingerprint over the detection rule
files** — `04add3eb0e86a816`, a sha256 over `strict.txt`, `relaxed.txt` and
`loose.txt`. The rule count is **103** at 8.0.4 (90 strict, 7 relaxed, 6 loose),
against the 88 at v7.1 and 58 at v5 this project had already written down. "Never
compare antiSMASH results across versions" was a standing instruction; the
fingerprint is its mechanical form, so a campaign pins what it was built against
and a changed rule set is detectable even when the version string is not what
moved.

**A banned binary arrives with it.** antiSMASH depends on the `fasttree` package,
which ships `FastTreeMP` — banned outright, since thread order affects its
neighbour-joining heuristic. The ban had been written down and nothing enforced it.
Now the image deletes the binary and `registry --verify` fails if it reappears,
which also sharpens a distinction worth keeping: `BANNED_BINARIES` is checked,
`NOT_IN_IMAGE` is a list of deferrals with reasons, and FastTreeMP belongs in
exactly one of them.

**It brings its own Python, and that nearly became the benchmark's.** With the tool
prefix first on `PATH`, `python3` silently resolves to antiSMASH's conda
interpreter, and the stdlib-only grading core would be running on an interpreter
nobody chose. It surfaced as `python3 -m pytest` failing to find pytest. antiSMASH's
console scripts carry an absolute shebang to their own interpreter, so the prefix
never needed to come first; the image now puts it last.

### BiG-SCAPE: it does not run as published

The sharpest result of the whole Phase 1 effort, and the one that justifies
running tools rather than installing them.

BiG-SCAPE 2.0.3's bioconda recipe asks for `sqlalchemy >= 2.0.2` with no upper
bound. A fresh solve therefore installs 2.1.x, and BiG-SCAPE raises
`ObjectNotExecutableError` **before it reads a single input file** — it passes a
compiled statement object to `Connection.execute`, which 2.0 tolerated and 2.1
rejects. The image carries `sqlalchemy=2.0.54`, a bound the upstream recipe does
not. **Installing is not the same as working**, and nothing short of running it
tells them apart.

Two smaller things the same session surfaced. BiG-SCAPE filters its input
directory by filename, defaulting to `cluster,region` to match antiSMASH's
`*.region001.gbk` convention; MIBiG reference files are named `BGC0000852.gbk`, so
without an explicit `--include-gbk` the run fails with *no valid input GBKs*
rather than with anything about names. And it needs `Pfam-A.hmm`, which the
antiSMASH databases already ship pressed — so the image reuses that copy, saving
1.5 GB and one separately-pinned release.

**The partition is the result; the labels are not.** `FAM_00001` is a name. Two
runs can agree completely on which clusters belong together and disagree on what
the groups are called, and a comparison of labels would score that as a
difference. The projection keeps, per class, the groups as sorted member sets,
sorted among themselves, and drops the family label and the connected-component
number. The cutoff stays in the *invocation* rather than the projection, because a
partition at a different cutoff is a different answer, not a different rendering
of one — which is exactly the distinction T-L3-2 ("GCF cutoff") is about.

**The fixture was built to have families.** 32 antiSMASH-processed MIBiG clusters,
selected by a rule computed from tables this repository already carries: the four
largest groups of 4–10 entries sharing an InChIKey connectivity block — ectoine
(10), ochratoxin A (7), kanamycin (5), aflatoxin B1 (5) — plus one singleton per
biosynthetic class as a negative control. BiG-SCAPE recovers those families at a
0.3 cutoff, which is what makes this a test of clustering rather than a test of
whether 32 singletons stay 32 singletons.

The one judgement in that selection is named rather than buried: groups whose
MIBiG compound name is a **class placeholder** — *capsular polysaccharide*,
*lipopolysaccharide*, *melanin*, *carotenoid*, *exopolysaccharide* — are excluded,
because MIBiG gives those generic names a representative structure and entries
sharing one need not be homologous. Two different molecular formulas appear under
*capsular polysaccharide* alone. Including them would have put apparent family
structure in the fixture that the biology does not support.

**One coupling to settle before T-L3-2 is authored.** The MIBiG reference set
ships as `mibig_antismash_4.0_gbk_as8b1.tar.bz2` — processed with antiSMASH **8.0
beta 1**, while the image pins **8.0.4**. Harmless for a determinism fixture,
since both arms of every comparison read the same files. Not harmless for a
campaign: comparing its own 8.0.4 regions against that reference set is precisely
the cross-version comparison antiSMASH's rules forbid. Either reprocess the
reference set with the pinned antiSMASH, or state that both sides come from the
published set.

*Settled.* T-L3-2 was authored as `gcfcutoff-mibig-32-01` and takes both sides
from the published set; so do the three other campaigns that read it, none of
which runs antiSMASH. The gate's check 25 now enforces it per campaign — see
"Item 5" below for the count and the field it reads.

### Projections, and why they are a third kind of licence

antiSMASH's output is one JSON object holding the input path, the tool version, a
record timestamp and the full HMM hit table in the same structure as the detected
regions. A line-drop normalisation cannot reach inside that. So the suite gained a
**projection**: a declared map from the raw artifact to the bytes that are actually
results.

| | What it says | Where it applies |
|---|---|---|
| **Normalisation** | this line is not a result, ignore it when comparing | inside the suite only |
| **Canonicalisation** | the content is right but its order is not guaranteed — fix the order before anybody reads it | every gold-producing use |
| **Projection** | these fields are the answer and the rest is not | every gold-producing use |

The two antiSMASH invocations project different things on purpose:
`minimal_detection` projects the regions, which is what a BGC-detection campaign
grades; `default_modules_domains` runs the analysis modules and projects the
ordered NRPS/PKS domain architecture per CDS, which is what a domain-architecture
or substrate-specificity campaign grades. The second drops the e-values and
bitscores — a float's last digits are a reduction-order artefact, and nothing
downstream needs them to state an architecture.

Two safety properties matter here more than they look:

- A projected invocation is reported as **`projection-identical`, never
  `raw-identical`**, because the raw artifact demonstrably is not identical. This
  is the same rule as the incidental-match one, one layer up: a claim must say what
  it is a claim about.
- A projection that finds nothing **raises**. A projection yielding empty bytes
  would turn the determinism test into a tautology — every run agreeing on the same
  emptiness — so an antiSMASH structure change or a fixture with no detectable
  cluster fails loudly instead of passing quietly. A test asserts the refusal.

### The antiSMASH fixture

Three complete deposited records, concatenated: *Vibrio anguillarum* 775 plasmid
pJM1 (anguibactin, NRPS, 65 kb), a *Kamptonema* landornamide cluster (ribosomal,
16.5 kb) and a *Streptomyces sampsonii* julichrome cluster (PKS, 16 kb). The run
detects three regions across four product names — `NRP-metallophore`+`NRPS`,
`lanthipeptide-class-ii`+`proteusin`, `T2PKS` — so the rules exercised are not all
of a kind, which a single-cluster fixture would not have achieved.

Three properties were chosen rather than stumbled into. Each record is **complete**
and its MIBiG locus starts at position 1, so the provenance is three accessions
with no coordinates to get wrong. Each is **annotated**, so antiSMASH's verdict is
about antiSMASH and not about a gene caller feeding it — which is why the
provisioning controls now forbid FASTA input on a gold-producing path. And all
three are **bacterial**, so one `--taxon` setting covers the file.

### Normalisation and canonicalisation are not the same licence

That second case forced a distinction the suite now carries explicitly, because
collapsing it is how a determinism claim becomes cheap:

- a **normalisation** says *this line is not a result, ignore it when comparing*.
  It lives inside the suite, and each rule carries a reason.
- a **canonicalisation** says *the content is right but its order is not
  guaranteed, so fix the order before anybody reads it*. It is part of the
  invocation contract for every gold-producing use — oracle and agent alike — and
  a tool that needs one and does not get it is not admissible as a gold source.

Sorting lines is sound for MMseqs2's tabular output only because it has no header
and each line is an independent record, which the declaration states. It would not
be sound for an aligned or sectioned format, and saying so is the difference
between a declared transform and a convenient one.

The suite reports the raw, normalised and canonical verdicts separately and names
the rules that fired. Reporting only the last would let any tool pass by stripping
whatever differs; reporting only the first would fail every tool that stamps its
own command line, which says nothing about determinism.

One smaller rule, which looks pedantic until it bites: **a raw byte match while a
normalisation rule fired is reported as incidental, not as raw-identical.**
`hmmbuild` writes a build date, and four runs inside the same second agree on it.
The first suite run reported `hmmbuild` as raw-identical for exactly that reason
and the next one did not, which is how the problem surfaced. Calling a timing
accident a property is the same error as laundering an unmeasured check into a
green tick, one layer down.

### The fixtures are committed, hashed, and one of them is synthetic

A determinism suite whose input moves measures nothing, so the fixtures are files
in the repo with recorded provenance. `halogenases.faa` is 24 reviewed UniProt
entries carrying Pfam PF04820: real homologs across the similarity range,
including RebH and PrnA — already in the benchmark, so no new provenance — plus
one fragment, because a short record is where a search tool's and an aligner's
tie-breaks show. `synthetic.fna` is reverse-translated from that set under a
declared codon cycle and spacer, giving Prodigal a contig with exactly 24 known
ORFs that regenerates from the repo with no download; a test rebuilds it in pure
Python and compares bytes. `halogenases.afa` is a frozen MAFFT alignment, and that
it came from a tool the suite also tests is not circular: once written it is a
file, so a later MAFFT regression cannot change the HMMER results computed from
it.


## Campaign 12: BRENDA kinetics, a corpus that checks itself

`kinetics-brenda-ec1_1-01` reconciles three numbers BRENDA stores separately and
arithmetic links: the Michaelis constant, the catalytic constant, and the
catalytic efficiency that is their quotient. No external standard, no curator —
the corpus is its own gold, and the shape of the disagreement is the finding.

Over EC subclass 1.1.- (437 EC numbers, 35,188 kinetic records, 2,528 complete
triples):

```
agree, within 5%                1598   63.2%
other                            708   28.0%
near, 5-25%                      166    6.6%
factor ~1000  (mM read as uM)     31    1.2%
factor ~60    (per-s as per-min)  25    1.0%
```

### The units are the task, and BRENDA records none of them

A value is a string — `0.05 {benzyl alcohol}` — a number and the substrate in
braces. The unit belongs to the **field**: Km in mM, kcat in s⁻¹, kcat/Km in
mM⁻¹s⁻¹. A system that does not know that cannot relate the three at all, and one
that assumes micromolar is wrong by exactly a factor of a thousand. R4 makes the
dependence explicit by applying the mistakes deliberately:

```
identity (control)      63.2%
Km read in micromolar    0.36%
kcat read per minute     0.32%
both                     0.04%
```

A 177-fold collapse from a unit assumption is as direct a demonstration as the
benchmark has that an implicit convention is load-bearing.

### Two traps that are in the data rather than in the task

**BRENDA writes -999 for a missing measurement** — large, negative, and present
in 545 records of this subclass alone. A pass that treats it as a number reports a
negative mean kinetic constant and silently loses every triple it belongs to. The
counting rules exclude it by numeric comparison, not string prefix, so `-999.0`
and `-9.99e2` are the sentinel and `-9990` is a measurement.

**Some values are ranges**, like `0.05 - 0.1`. They are counted as unparsed, never
averaged into a number nobody wrote down. 269 records in this subclass.

### The key has to be tight, and a test proves it

A triple is keyed on `(EC, proteins, substrate, references)`. Dropping the
reference pairs values from different papers; dropping the protein pairs different
enzymes. Either way the disagreement rate stops meaning anything — so a test
builds the loose key as well and asserts the two give different answers. A
counting rule nobody can tell the difference from its alternative is not doing
work.

### What this row cost the catalogue: two corrections at once

**The access was not blocked.** The catalogue recorded T-L2-5 as blocked on BRENDA
credentials. It is not: `download.php` serves a form, a licence-acceptance
checkbox enables the download buttons, and a POST with
`dlfile=dl-json&accept-license=1` returns 83 MB of JSON with no account and no API
key. The licence is plain CC BY 4.0 and the page says so. The earlier verdict came
from guessing an archive URL and reading its 404 as a wall.

**And the catalogued task was not buildable as written.** "Kinetics harmonisation
across unit conventions" presumes competing conventions in the source; BRENDA has
already normalised units per field and records none. The catalogued *intent*
survives — the unit conventions really are what the campaign turns on — but the
mechanism is the corpus checking itself rather than two sources being brought into
line. That is now five catalogued designs corrected by grounding, and the first
where the correction and the access correction arrived together.

## The ChEMBL decision

ChEMBL is CC BY-SA 3.0, and the project had been carrying it as an open licensing
question. **Decided: use it, under a new `license_class: share_alike`.** The
reasoning, recorded so it can be overturned on its merits rather than re-derived:

- ShareAlike attaches to an **Adaptation**, not to a **Collection**. A filtered
  ChEMBL subset is an adaptation, and so is gold computed from ChEMBL values, so
  both carry CC BY-SA 3.0 onward. A benchmark made of separable campaigns is a
  collection, not an adaptation of any one of them — so the code, the grader and
  every other campaign are untouched.
- The obligation is therefore satisfiable and **local**: attribution, the same
  licence on those data files, and a notice. A share-alike campaign can be dropped
  without touching anything else, which is exactly why the class is per-campaign.
- The alternatives do not help. IUPHAR/Guide to Pharmacology is CC BY-SA 4.0 — the
  same condition. PubChem BioAssay largely mirrors ChEMBL, so routing through it
  would be laundering rather than compliance. And ChEMBL is the source a reviewer
  expects for measured IC50s.

Made enforceable rather than asserted: `license_class` now admits
`open | nc | share_alike`, and the gate's nineteenth check requires a
`share_alike` campaign to carry a `license_notice` naming the source, the licence,
the attribution and what a redistributor must do. A licence condition recorded
only in a design document is a condition the person redistributing the files will
never see.

**What is left for the owner**: whether an evaluator's legal position rules out
shipping ShareAlike content at all. If so, T-L5-1 comes out and nothing else
changes — which is the property the per-campaign class was chosen to give.
**Exercised since**: `selectivity-saureus-topoisomerase-01` is the first
`share_alike` campaign and the gate's `share_alike_notice_present` check passes on
it, so the nineteenth check is confirmed against a real campaign rather than a
fixture. T-L5-2, the other row this decision was taken for, turned out not to be
blocked on licensing at all — see below.

## Campaign 13: ChEMBL selectivity, where the licence travels with the files

`selectivity-saureus-topoisomerase-01` is T-L5-1 and the first campaign in the
benchmark whose data carries an onward obligation. *S. aureus* DNA gyrase
(CHEMBL3038482) against topoisomerase IV (CHEMBL3038508), ChEMBL_37: the two type
II topoisomerases an antibacterial of this class can hit, so whether a compound
hits one or both is the real selectivity question rather than an arbitrary pair.

The admission rules **are** the task. Every one reads a structured ChEMBL field,
and two of them are where a potency table goes wrong in practice:

- **`confidence_score` lives on the assay, not the activity.** Admission needs an
  activity-to-assay join. A system that looks for the field on the activity record
  finds nothing and admits everything, and nothing in the data says it went wrong.
- **A `standard_relation` of `>` or `<` is a bound, not a measurement.** 235 of the
  1,117 IC50 records on this pair are exactly that. Admitting one as a value is the
  commonest way a potency table acquires numbers nobody measured.

Measured, with the rejections partitioning exactly — 739 of 1,117 admitted; 235
censored, 81 in mass-per-volume units that cannot be converted without a molecular
weight the campaign refuses to import, 33 flagged by ChEMBL, 29 potential
duplicates. 668 (molecule, target) pairs carry an aggregated value and **173
molecules are measured against both targets**. At twofold, 122 molecules favour
gyrase and 21 favour topoisomerase IV, so neither direction is a foregone
conclusion; the build refuses to emit gold if either side is empty. The most
selective compound is CHEMBL4555272 at 16,666.7×.

**The rejection breakdown is first in `notes_for_audit` for a mechanical reason**:
the reasons are tested in a declared order, each activity takes the first that
applies, so admitted plus rejections must equal the IC50 total per target. If that
sum fails, two rules are double-counting and every number below it is suspect. That
is a check an auditor can run in one line, which is the point.

### R4 had to be redesigned, and the reason generalises

The catalogued R4 was "a low-confidence assay is excluded — predict the new ratio".
It is a **no-op on this pair**: all 105 assays behind these IC50 values score 7,
"Direct protein complex subunits assigned". So the rule stays declared and graded —
a rule that happens to be a no-op on one instantiation is not a rule that can be
dropped from the contract — but a counterfactual has to vary an axis that moves.
R4 relaxes the censored-relation rule, the ChEMBL flag and the duplicate rule, and
swaps the aggregator. All five variants move something: admitting censored
relations takes eligibility 173 → 193, and **minimum instead of median leaves the
admitted count identical at 739 while moving the focal counts 122/74/18 →
139/90/25** — a reduction choice, not an admission choice, and it changes the
answer. That variant is the one worth keeping in mind: an agent can get every
admission decision right and still report a different selectivity profile.

This is the second half of the grounding rule again, in its sharpest form yet:
`confidence_score` is present on every assay, and computing its distribution is the
only thing that shows it carries no information here.

## T-L5-2: a field that is populated and says nothing

T-L5-2 (resistance mutation → target assignment) was the obvious companion build —
same source, same licence decision, and the re-grounding probe had confirmed the
structured `assay_variant_mutation` field it needs. **It does not survive counting.**
ChEMBL_37 has 34,921 IC50 activities carrying that field, and per on-theme target:

| Target | Variant records | Usable |
|---|---|---|
| *S. aureus* gyrase complex | 45 (D83N 41) | 1 clean WT/mutant pair |
| *S. aureus* GyrA | 39 (S84L 32, D83N 6) | mutant values mostly `None` |
| *S. aureus* ParE / topoisomerase IV complex | 0 | — |
| *P. falciparum* DHFR | 244 | 0 — all "UNDEFINED MUTATION" |
| HIV-1 RT | 5,423 (Y181C, P236L, K103N) | off-theme |

44 wild-type/mutant key pairs exist across the bacterial targets, but most mutant
`standard_value`s are `None`, so they cannot be turned into a ratio. The one clean
pair in the set is ciprofloxacin against gyrase: WT 5.0 nM versus D83N 580,000 nM.
A campaign on that is an anecdote, and a campaign on HIV-1 RT grades resistance
reasoning on an antiviral target whose chemistry is nothing like this benchmark's.

So the row is **re-tiered, and not for the reason it was previously blocked**:
ChEMBL is accepted, the data is thin. The alternatives are named in
`docs/catalog.md` and left to the owner — instantiate on HIV-1 RT and accept the
theme drift, re-scope the ladder to the census itself (how much resistance data
exists per target class is a real and gradeable question), or wait for a release
that fills the bacterial values. **The counts are recorded so the next person does
not re-probe.** This is also a correction to the re-grounding note, which called
T-L5-2 "data-verified": the probe verified the field, not the data behind it.

## Campaign 14: the first one that runs a tool

`ecaudit-fmn-dh-1_1_3_15-01` is T-L2-1 and the benchmark's first **S1** campaign.
Everything before it was pure computation over shipped data; this one needs HMMER
on the path and the pinned Pfam release in the image, and no S0 system can answer
any rung of it. That is declared in `task.yaml` rather than left for the results
matrix to discover, so the empty S0 column reads as tooling and not as capability.

UniProt 2026_03, every entry annotated EC 1.1.3.15 — an FMN-dependent
(S)-2-hydroxy-acid oxidase, so the annotation implies the FMN_dh domain.
**7,675 sequences: 29 reviewed and 7,646 unreviewed, both sections in full.** No
sampling, because the reviewed/unreviewed split is the campaign's control and a
sampled control is not one.

| | FMN_dh present | absent | rate |
|---|---|---|---|
| reviewed | 29 / 29 | **0** | 0.0000 |
| unreviewed | 6,300 / 7,646 | **1,346** | **0.1760** |

**The published figure does not reproduce, which is exactly why the G1 rule
exists.** Rembeza & Engqvist found 78% of EC 1.1.3.15 proteins in BRENDA 2017.1
lacked the FMN-dh domain, and that paper is why this row is in the catalog. On
UniProt 2026_03 the measured figure is 17.6%. Had the campaign graded the paper's
number it would have been a literature-recall test, and it would have been wrong.
Gold here is what the pinned tool says about the shipped sequences, and nothing
else; the literature figure is in `excluded_from_grading` with that reason.

### The two findings that make it a ladder rather than a tool invocation

**The absent sequences are not fragments.** This was the explanation that would
have made the whole audit an artefact: a sequence can miss a 348-column model
simply by being short, in which case "lacks the domain" is a statement about
sequencing completeness. Measured, the canonical-absent side is **longer** than
the present side — median 466 residues against 367 — and 2.4% of it falls under
300 residues against 12.6% of the present side. The length profile is a graded R2
component and the first item in `notes_for_audit`, because a refutation nobody can
check is just an assertion.

**They are a coherent enzyme family.** 967 of the 1,346 (71.8%) carry exactly
`FAD-oxidase_C + FAD_binding_4`: the FAD-linked glycolate oxidase subunit
architecture, which oxidises the same substrate through a different cofactor. The
rest are mostly `DAO` combinations, and 11 carry nothing from the panel. So the
absent set is not junk, and **absence of the canonical domain is not
misannotation**. The campaign therefore refuses to report a misannotation rate and
grades the count under three published policies instead:

| policy | accepts | misannotated |
|---|---|---|
| `p0_canonical_only` | FMN_dh | **1,346** |
| `p1_accept_fad_oxidase` | + the FAD-oxidase architecture | **319** |
| `p2_accept_any_flavin` | + any panel flavin family | **22** |

The campaign does not rule between p0 and p1. The gap between them *is* the
finding, and publishing the policy rather than asserting a rate is what keeps a
human judgement out of gold.

### R4 had to be rebuilt, and one of its axes is deliberately inert

The catalogued R4 was "a sequence is mutated to disrupt the FMN-dh domain —
predict the reclassification". It cannot be built honestly: mutating a sequence to
destroy an HMM match means choosing which residues to break, any choice large
enough to drop a curated gathering threshold is a choice about the answer, and
gold would then describe a sequence nobody deposited. R4 varies the two axes
already in the data — the declared cutoff and the declared policy.

The cutoff axis **barely moves**: 1,342 to 1,371 absent across six cutoffs, with
`--cut_ga`, `--cut_nc` and `--cut_tc` giving an identical answer. The baseline
cutoff is among the six, which makes it that axis's **identity control** — it
re-runs the audit under the flag that produced the headline and has to land on the
same number, the way the kinetics ladder's `v0_identity` variant does, and the
build refuses to emit gold if it drifts. Normally a
counterfactual that does not move is the defect that re-scoped the selectivity
campaign. Here it is a robustness measurement — the headline rate does not depend
on where the threshold sits — and it is admissible **only because the policy axis
in the same rung moves by a factor of sixty**. The build asserts exactly that: it
refuses to emit gold if both axes go flat, and refuses if the strict policy stops
reproducing the absence count.

The catalogued discrimination ECs are not used either, and measurement is why:
1.13.12.4 has 204 entries, 1.1.99.31 has 234, and **EC 1.1.2.3 is itself an FMN_dh
family member** — all three reviewed L-lactate dehydrogenase (cytochrome) entries
hit the model at the gathering threshold — so the canonical domain does not
discriminate it from EC 1.1.3.15 at all. Three more sequence sets, no signal.

### Declared, not discovered: how the tool stays affordable

The model panel is **published in `reference/`**, 13 Pfam models. It was fixed
once at build time by scanning the canonical-absent set against the whole of Pfam
35.0 and listing every family that fired above the gathering threshold. An agent
then scans with exactly those models: **11 seconds**, against the **2m36s** a full
Pfam-A pass costs. Making the discovery part of the task would have hidden the
answer behind an open-ended search whose cost is multiplied by every rung, every
repeat and every stub level in the sweep.

It also **partitions the set**: all 6,329 canonical-present sequences carry FMN_dh
and nothing else from the panel, so the architecture census is a real partition
rather than a ranking of overlapping families. That is a test, because a panel
that overlapped the canonical model would make the dominant-architecture claim
depend on which family happened to score higher.

### A resource the campaign may read but not ship

Pfam-A.hmm is 1.5 GB, so it is not redistributed. What replaces the bytes is a
declaration: the path, the provider (the antiSMASH 8.0.4 database layer), the
version and a **sha256**, recorded in `reference/audit_rules.json`. The build
refuses to emit gold if the file it finds does not hash to what it recorded,
because a domain verdict against a different Pfam release is not this campaign's
gold.

This is the mechanism the T-L2-4 grounding said every S2 campaign would need, and
the gate now has a **twentieth check**, `external_resources_pinned_and_present`,
which resolves each declared resource and verifies its fingerprint. Here the
reason for not shipping is **size** and Pfam is CC0 — recorded explicitly, because
the next campaign to use this mechanism will be withholding for licence instead,
and "we did not ship it" means different things in the two cases.

### The sweep, and what it cost

All five sweep gates pass, with cold depths exactly where the stub ladder predicts:

| stub | mode | cold depths |
|---|---|---|
| `stub-reader` | tooled | 1, 1, 1 |
| `stub-presence` | tooled | 2, 2, 2 |
| `stub-architecture` | tooled | 3, 3, 3 |
| `stub-complete` | tooled | 4, 4, 4 |
| `stub-inverted` | tooled | **1, 1, 1** |
| `stub-noncompute` | no_tool | 0, 0, 0 |

Per-rung clear rates 1.000 / 0.600 / 0.400 / 0.200, mean cold score 0.55, and the
no-tool ablation scores 0.0000 against a 0.2505 chance floor. `stub-inverted`
landing at depth 1 is the result worth reading: it runs the tool correctly with
the right panel and reads the hit table backwards, so it clears the inventory rung
and fails at presence — the rung whose claim the column swap actually corrupts.

The first run took **29m32s of wall clock and 47m of CPU**, against seconds for a
stdlib campaign, because every stub level re-ran the same HMMER call for every
rung and every repeat over an input that never changes — roughly a hundred
identical tool runs. Fixed below, and now **3m44s** with every gate result and
every depth bit-identical.

### A memo for tool output, and where it is not allowed

`npbench_c.tools.cache` memoises an expensive tool call to a JSON file. The
mechanism is three lines of hashing; the design is entirely about where it may
apply, because a memo in the wrong place turns a measurement into a tautology.

- **Off unless `NPBENCH_TOOL_CACHE` names a directory.** A build, a test, a human
  invocation and a real agent run all compute from scratch.
- **Only a stub, and only in tooled mode.** The sweep runner sets the variable for
  stub systems alone and **raises `ToolCacheMisuse`** if handed a memo for anything
  else. A sweep over real systems is partly a measurement of whether they can
  drive the tool at all, and a `no_tool` ablation's entire claim is that it did not
  compute — handing either one a memo would quietly void the leakage check.
- **The key is the whole identity of the computation**: the input file's digest,
  the library's fingerprint, the tool version, the panel and the flags. A key
  missing any of those would survive a change it must not survive and return an
  answer from the wrong world. Key parts are validated *before* the
  caching-disabled early return, so a malformed key fails everywhere rather than
  only on the one run where caching is on — a key defect that surfaces only inside
  a sweep is a key defect nobody sees.
- **The memo lives beside the sandboxes, never inside one**, so the gold-isolation
  audit still sees a clean sandbox, and it only ever holds tool output derived
  from shipped inputs.

Measured on this campaign: **29m32s → 3m44s**, CPU **47m → 3m12s**, six memo
entries (one per declared cutoff) in place of about a hundred identical runs, and
the gate block of `internal_sweep.json` byte-identical to the uncached run.

The cost of the mechanism is that a cached sweep's wall clock is no longer a cost
measurement of the campaign, so `_meta.tool_cache_note` says exactly that in the
report and a test asserts the sentence is there. That is the honest trade: the
sweep exists to check the grader and the ladder, and it now does that in minutes;
anyone who wants the campaign's real tool cost runs it without the variable.

### The parsing mistake that does not fail

`hmmsearch --domtblout` puts the **target** (a sequence) in column 1 and the
**query** (a model) in column 4. Reading them the other way round reports model
names as sequence accessions and accessions as models: no error, a
plausible-looking table, every number wrong. I made exactly that mistake while
grounding this row and it produced a convincing census of families named
`tr|Q7NQA5|Q7NQA5_CHRVO`. There is now a test that the returned mapping is keyed
by the shipped accessions and valued by the declared panel, and the `inverted`
stub level is the graded form of the same error — it runs the tool correctly, with
the right panel, and reads the result backwards.

## Campaign 16: the first S2, and two artefacts it had to refuse

`bgcdetect-scoelicolor-01` is T-L3-1 and the benchmark's **first S2 campaign**: it
needs antiSMASH 8.0.4 with its database layer, and it reads the three detection
rule files from the image. The complete *Streptomyces coelicolor* A3(2)
chromosome (AL645882.2, 8,667,507 bp) goes in; 29 regions over 22 distinct products
come out, covering **13.0% of the chromosome**, median span 29.8 kb, and **every
region `contig_edge=False`**.

### The slice that looked obvious and was wrong twice

A whole chromosome is a 6.8 MB input and the obvious move is to slice windows
around the clusters of interest. That fails on two counts, and both were measured
rather than reasoned about:

1. **A slice truncates CDS features at its edges**, which arrive carrying
   GenBank's `<`/`>` partiality markers, and antiSMASH rejects the record outright:
   `feature translation extends out of record`. Dropping those features makes it
   acceptable, which is how the second problem stays hidden.
2. **A region that reaches the slice edge has its boundary set by the window.** On
   an 81 kb window around the actinorhodin cluster — a 30 kb flank on each side —
   antiSMASH returned a single T2PKS region running to the record end with
   `contig_edge=True`. The boundary is what R2 and R3 grade, so that campaign would
   have graded my flank choice.

On the whole chromosome neither happens, and a `--minimal` pass takes **75
seconds**. The slice bought nothing at all. The build refuses to emit gold if any
region touches an edge, so the artefact cannot creep back in through a
re-instantiation.

### The reconciliation, which is the actual campaign

Running a tool and reporting its output is an R1. What makes this a ladder is the
comparison with the 15 cluster boundaries MIBiG 4.0 maps to the same accession:

| | |
|---|---|
| curated loci detected | **15 of 15** |
| fully inside a called region | 14 |
| partially covered | 1 (CDA, 95.3%) |
| median Jaccard | **0.295** |
| Jaccard range | 0.023 – 0.934 |
| called regions with no curated locus | **15 of 29** |

So detection is not the question — extent is. antiSMASH extends a region from the
protocluster core by a class-dependent distance, and this campaign measures by how
much rather than calling it an error. The extremes are instructive on their own: a
1.7 kb signalling cluster sits inside a 75.6 kb region (Jaccard 0.023) while CDA
agrees to 0.934. And **half of what the tool finds has no curated counterpart**,
which is a statement about the coverage of the curated record rather than a false
positive — MIBiG records what has been characterised.

One region holds two curated loci, because the SCB1 butyrolactone genes lie inside
the coelimycin cluster. A one-to-one join would have silently dropped one of them,
so the reconciliation reports which regions hold several.

**The two coordinate conventions differ, and the difference is verified rather than
assumed.** antiSMASH writes 0-based half-open; MIBiG's `from`/`to` are 1-based
inclusive. That is established against the source annotation: BGC0000194's `to` is
5,535,091, and SCO5092 ends at 5,535,091 in AL645882.2. A MIBiG locus therefore
spans `to - from + 1` bases, one more than its own numbers suggest at a glance. I
had assumed half-open before checking, which would have shifted every interval by
one base.

### R4 is the catalogued counterfactual, plus the control it needed

Delete the core biosynthetic genes of a declared region from the annotation — the
sequence untouched, so every other coordinate stays put — and run again:

| variant | outcome |
|---|---|
| `v0_identity` | nothing deleted; reproduces the baseline exactly |
| `v1_t2pks_core_both` | actinorhodin region **abolished** |
| `v2_t2pks_core_one` | **also abolished** — the rule is conjunctive |
| `v3_nrps_core` | CDA region abolished |
| `v4_lanthipeptide_core` | lanthipeptide region abolished |
| `v5_additional_not_core` | **nothing changes** |
| `v6_multiproduct_one_product` | region survives with one product fewer |

My first attempt had only v0–v4, and **the build refused it**: every variant lost
exactly one region, so the signature was constant and the rung was testing one
dependency four times. v5 and v6 are what the refusal bought. v5 deletes the
actinorhodin acyl carrier protein, which antiSMASH marks
`biosynthetic-additional` rather than core — the negative control, without which
"a deletion costs a region" is a free answer. v6 deletes a gene core to two nested
protoclusters, so the region is replaced rather than removed. The build now
requires three distinct outcome signatures, one of them a no-op and one of them an
outright loss.

### The sweep, and the cost the memo does not cover

All five sweep gates pass, with cold depths 1/2/3/4/1/0. `stub-offbyone` landing
at depth 1 is the result worth reading: it runs antiSMASH over the whole
chromosome correctly, projects all 29 regions correctly, and reports the 1-based
inclusive start without the conversion. It clears the inventory and fails at the
boundaries — the rung whose claim one arithmetic operation corrupts.

**It took 15m35s**, against 55 seconds for the S1 transfer campaign. The tool memo
is working — 7 antiSMASH runs serve all 72 rung-repeats — so what remains is not
the tool: it is reading and re-parsing a 26 MB GenBank record once per run, plus
copying a 6.8 MB input into each sandbox. **The memo caches the tool call, not the
work around it**, which is a limitation worth stating now rather than discovering
at the third S2 campaign. If S2 sweeps become routine, the next move is to cache
the parsed record alongside the tool output, or to let a stub level reuse one
sandbox across repeats.

### What else the gate caught here

`product_census` was graded at R2 **and** handed over by R2's warm bundle, because I
had put it in the inventory that R2's bundle publishes. Two scalars beside it,
`region_bases` and `fraction_of_record`, had the same defect and escaped only
because the warm-bundle check examines composites. The fix was to settle what each
rung is *about*: R1 is the inventory — run the tool and summarise — and R2 is the
per-region boundaries. The rung split now follows that, and the bundle contains
exactly R1's answers.

## Campaign 17: the tool's default cutoff recovers none of them

`gcfcutoff-mibig-32-01` is T-L3-2 and the second S2 campaign. BiG-SCAPE partitions
BGCs into families, and the partition is a function of a distance cutoff the caller
supplies — so the campaign runs the pinned tool **once** over a declared eight-point
grid and asks what the partition does across it.

The corpus is 32 antiSMASH-processed MIBiG clusters, and the thing they are
reconciled against comes from **chemistry rather than from clustering**: four
groups of entries sharing an InChIKey connectivity block in the chemical-space
campaign's table, plus five negative controls, one per biosynthetic class,
belonging to no such group. That independence is what makes R3 a comparison rather
than a restatement.

| cutoff | groups | recovered pairs | missed | false joins | pair Jaccard | families exact |
|---|---|---|---|---|---|---|
| 0.1 | 25 | 13 | 73 | 0 | 0.151 | 0 |
| 0.2 | 21 | 19 | 67 | 0 | 0.221 | 0 |
| **0.3** | 17 | 33 | 53 | 0 | **0.384** | **0** |
| 0.4 | 15 | 44 | 42 | 0 | 0.512 | 0 |
| **0.6** | 12 | 71 | 15 | 0 | **0.826** | 2 |
| 0.7 | 11 | 71 | 15 | **5** | 0.780 | 2 |

**At 0.3 — BiG-SCAPE's own documented default — not one of the four
chemistry-derived families comes back as exactly one group.** Every one is split;
one across five groups. Pair agreement is 0.384, less than half what the corpus
allows. The best cutoff is 0.6 at 0.826, and at 0.7 a negative control is absorbed
into a family, producing the first five false joins. So agreement is **not
monotone**, and both ends of the curve are inside the grid — the build refuses to
emit gold if the optimum sits at the edge, because then the claim would be an
artefact of where the sweep stopped. Two families are never recovered exactly at
any declared cutoff, which reports as `null` rather than a number.

### Three ways to get a plausible wrong answer

**The clustering table is not the partition.** It lists only clusters that joined
a family — at cutoff 0.1, ten of thirty-two records. The partition has to be
reconstructed with every omitted record as a singleton. Read directly, the table
gives a smaller corpus and a tidier family structure than the tool produced, and
nothing fails. That is the `tableonly` stub level, and it lands at depth 1.

**The join key is `GBK`, not `Record`.** `Record` carries a region suffix
(`BGC0000854.gbk_region_1`). The engine checks that each GBK appears exactly once
per cutoff and that every row is a `region` record, because a multi-region input
would need a different reconstruction rule than the declared one — checked rather
than assumed, after the antiSMASH campaign's column-order lesson.

**A family label is a name.** `FAM_00007` is not a result: two runs can agree
completely on which clusters belong together and disagree on what the groups are
called. So the partition is graded and agreement is counted over **pairs**.

### `interval_score` is the right primitive for the wrong campaign

The catalogue specified R4 "scored with `interval_score`". It is not used here, and
the reason is in that primitive's own docstring: `max_width` is "3x the observed
oracle spread from the determinism audit", so that both thresholds are *measured
rather than chosen*. This campaign has no spread to measure — BiG-SCAPE is
deterministic, every cutoff comes from one run, and the answer is exact on a
declared grid. Any `max_width` would therefore be a number chosen by me, which is
precisely what the primitive exists to prevent. R4 grades the per-cutoff table and
three boundary claims exactly instead.

That is worth separating from the other catalogue corrections in this document.
The earlier ones were cases where the *data* did not support a design. This one is
a case where **the campaign is too deterministic for the scorer**: `interval_score`
exists to make a system commit to a bracket when the quantity is uncertain, and
here it is not.

### Two pins the catalogue got backwards

`--mibig-version` was in the catalogued pin list. It is **forbidden** on a
gold-producing path: it downloads a reference set at run time, making the run
network-dependent and the reference set unpinned. The corpus is shipped instead.

And the open item about the MIBiG reference set's antiSMASH version closes here, in
the "declare" branch — correctly, for a reason worth stating. Every input in this
corpus was processed by the same antiSMASH version, and the campaign never compares
them against output from the image's antiSMASH 8.0.4. The cross-version rule
forbids *mixing* versions in one comparison; it does not require that every file in
the benchmark come from one version. No run here mixes.

The sweep passes all five gates with cold depths 1/2/3/4/1/0 in **6m27s** — against
15m35s for the whole-chromosome campaign, because this corpus is 32 small files
rather than one 26 MB record, which is the other half of the sweep-cost finding
above.

## The leak the sandbox was handing over in prose

Writing the first S2 campaign's rules file, I explained the decision to ship a
whole chromosome by saying "on the whole chromosome all 29 regions come back
contig_edge=False". **29 is R1's graded answer**, and `reference/` is a file the
agent reads before it has done anything. My own gold-isolation test caught it.

Checking the rest of the project found the same leak in **seven of sixteen
campaigns**, four of them written in earlier sessions:

| campaign | stated in the cold sandbox |
|---|---|
| `bgcdetect-scoelicolor-01` | the record length, and the region count in prose |
| `ecaudit-fmn-dh-1_1_3_15-01` | `sequences_total`, in an exclusion's reason |
| `selectivity-saureus-topoisomerase-01` | `activities_of_declared_type`, in the scope note |
| `transfer-ec1_14-uniprot-01` | `partial_ec_mentions` and `queries_with_a_hit`, twice each |
| `constraint-gc62-ecoli-rebh-01` | `unrestricted_gc_max`, `restricted_gc_max`, and most of R4 |
| `constraint-gc55-ecoli-rebh-01` | `restricted_gc_max`, and most of R4 |
| `mibig-diff-3_1-to-4_0-01` | `shared_entries` and both raw difference counts |

The constraint pair is the worst of them: the task description contains a worked
feasibility table quoting both GC maxima, which are R2's answers and also the
value most of R4's `relax_*_gc_max` components take. A real agent could clear R2
and most of R4 by reading the task description.

**Every one of these is invisible in everything the project already checks.** The
oracle still scores 1.0, because the oracle reads gold. The no-tool ablation still
scores zero, because a stub does not read prose. The warm-bundle check guards each
rung against its own bundle and says nothing about `task.yaml`. Only a real system
benefits — and then the per-rung difficulty profile is fiction in exactly the
direction that flatters the benchmark.

### Why it happens, and the rule

The mechanism is specific and worth naming: **a note explaining a design decision
drifts into quoting the measurement that motivated the decision.** Every leaked
number above sits in a sentence that was true, useful and well-intentioned — "the
reason we ship the whole chromosome is that all 29 regions then come back clean".
The explanation needs the phenomenon. It does not need the number.

So: `task.yaml` and `reference/` publish the objective, the contract, the
conventions and the vocabularies. **The measurements that motivated them belong in
`docs/`, which the sandbox never sees.** Where an audit note genuinely wants a
number, it says where to recompute it from rather than printing it: an auditor has
gold, and a system does not.

### Made mechanical, because a rule applied by hand drifts

The readiness gate's **twenty-first check**,
`cold_sandbox_withholds_answers`, builds the cold sandbox, strips `inputs/`, and
searches `task.yaml` and `reference/` for the string form of every graded value.
All sixteen campaigns pass it now.

Three design points in it, each one a false positive it had to stop producing:

- **Only distinctive scalars.** A graded count of 2 or a rate of 0.5 occurs in any
  prose by coincidence, so an integer needs four or more digits and a float four or
  more decimals. That leaves small graded values unchecked, which the docstring
  says plainly rather than implying coverage it does not have.
- **Word boundaries on integers.** `informative_records` is 2359 and the shipped
  chemistry table contains `BGC0002359`. Two campaigns failed on that substring
  before the match was anchored.
- **`inputs/` is excluded, and empty composites are skipped.** A number in `inputs/`
  is the data the agent was asked to read; a value readable straight out of it
  without doing the work is a *free component*, which is a different defect. And
  `[]` serialises to a string that occurs in any YAML file.

That last distinction found the only two remaining failures, and both were free
components rather than disclosures. `mibig-diff`'s `retired_class_terms` is
published in the class mapping — which the campaign *must* publish to be solvable
— so it is now dropped from grading with the exclusion recorded. That makes **five**
free components caught across the project, and the rule from campaign 15 holds
again: before grading a component, ask what would have to change in the world for
its value to change.

## Campaign 15: the closest hit is often the wrong donor

`transfer-ec1_14-uniprot-01` is T-L1-3, and it is the first catalogued design this
session that survived grounding almost intact — including the R4 that looked
expensive. The task is the oldest shortcut in functional genomics: take a
protein's closest database hit and copy its annotation across.

UniProt 2026_03, all **4,632 reviewed entries under EC 1.14** — oxidoreductases
acting on paired donors, the class holding the cytochromes P450 and the
flavin-dependent halogenases this benchmark already uses as subjects. The set is
searched against itself with pinned DIAMOND 2.2.8, the self hit is dropped while
the hit table is parsed, and each query's closest remaining hit donates its EC
number. **The reference set is closed and shipped**: transfer against a live
database would make gold a function of the day it ran.

| EC level | decidable | agree | rate |
|---|---|---|---|
| 3 (sub-subclass) | 3,650 | 3,556 | **0.9742** |
| 4 (serial) | 2,743 | 2,531 | **0.9227** |

And the curve the ladder exists to expose — level-4 agreement against the identity
of the closest hit:

| identity | n | agreement |
|---|---|---|
| 30–40% | 78 | **0.462** |
| 40–50% | 146 | 0.740 |
| 50–60% | 219 | 0.767 |
| 60–70% | 254 | 0.890 |
| 70–80% | 360 | 0.950 |
| 80–90% | 525 | 0.962 |
| 90–100% | 1,144 | **0.996** |

That is the twilight zone in one column, computed rather than cited. The build
refuses to emit gold if agreement does not rise with identity, because a campaign
whose central claim fails on its own instantiation is not a campaign.

R4 is the catalogued counterfactual verbatim — take the closest hit out of the
reference set, transfer from the next — and **it needs no second search**: removing
a hit means reading further down a ranking already computed. 501 of 4,552 verdicts
change, 111 of them from `supported` to `wrong_at_serial`; removing two moves 722.
The zero-removal variant is the axis's identity control and must reproduce the
headline census exactly.

### Three declared rules, each one load-bearing and each one measured

**A tie-break, because 602 queries have one.** One query in eight shares its top
bitscore with another hit, so "the closest hit" is not a function of the data until
the rule is published. The declared rule is bitscore descending, then subject
accession ascending. DIAMOND's own emission order is stable for this version — the
exact invocation came out byte-identical across 1, 2 and 8 threads over two repeats
each — but gold resting on an undeclared property of a particular build is gold
that breaks silently on a tool upgrade. Grounding this campaign with an ad-hoc
tie-break and then with the declared one shifted level-3 agreement from 3,558/3,652
to 3,556/3,650: small, and exactly the kind of small that nobody can reproduce.

**A dash is not a value.** 1,850 of this set's EC mentions stop early. An entry
annotated `1.14.-.-` makes no claim at level 3, so it can neither agree nor
disagree there, and a comparison is **decidable** only where both sides specify the
level. 902 queries land in `undecidable` rather than being scored either way.
Getting this wrong does not produce an error — it produces a report that looks
*more* complete than gold, which is why the `dashvalue` stub exists: it runs the
search correctly and compares EC strings with the dashes intact, clears R2 and
fails R3.

**Multi-EC is any-vs-any.** 446 entries carry more than one code, and a
multifunctional enzyme's annotation is a set rather than a ranking.

### The free component, caught a fourth time

EC levels 1 and 2 agree **1.0** and that is not a result: the query set is defined
by its EC prefix, so a hit inside the set shares the first two levels with every
query by construction. Both are computed and reported as `constant_depth_agreement`
for the record and **excluded from grading**, because a component whose value is
fixed adds weight to a conjunctive product and no information. That is the same
defect as `corpus_release`, `compounds_enumerated` and the pocket-geometry rung
redundancy — four times now, which makes it worth stating as a rule: **before
grading a component, ask what would have to change in the world for its value to
change. If the answer is "the campaign's own definition", it is not a measurement.**

The sweep passes all five gates with cold depths 1/2/3/4/2/0 in **55 seconds**,
which is what the tool memo bought: the same campaign under the previous runner
would have re-run the all-against-all search about a hundred times.

## T-L2-4: four comparison keys, four different answers

T-L2-4 (NRPS A-domain substrate specificity) was the next build after T-L5-1 — the
only `build: image` row with verified data, so the one that would have opened the
image tier. Grounded against MIBiG 4.0 and the Stachelhaus signature table
antiSMASH 8.0.4 ships, **it does not carry a ladder**, and the way it fails is
worth more than the campaign would have been.

The resources join on protein accession: MIBiG has 2,901 A-domains over 1,228 gene
accessions, the table has 2,319 rows over 1,637, and 714 accessions are shared. Of
those, **304 genes carry exactly one A-domain on both sides** — the only pairs that
can be matched without assuming MIBiG's module order equals the table's `.A1/.A2`
index, which nothing in either resource confirms (180 genes with equal counts above
one, holding 473 A-domains, hang entirely on that assumption).

On those 304 pairs, each comparison key tells a different story:

| key | decidable | agreement |
|---|---|---|
| raw substrate names | 304 | **0.000** |
| `aaSMILES` comment field read as a name table | 193 | 0.912 |
| stereo-aware InChIKey | 219 | **0.151** |
| connectivity InChIKey (block 1) | 219 | **0.995** |

Read top to bottom that is one resource compared with itself four times. Raw names
score zero because the vocabularies differ (`Ile` against `isoleucine`).
Stereo-aware structures score 0.151 because `aaSMILES` writes serine flat and
MIBiG writes the L-enantiomer — an annotation convention, not a disagreement about
chemistry. At connectivity level, the key this project already pinned for the
chemical-space campaign, the two resources agree on **218 of 219** comparable
A-domains; the single substantive disagreement in the whole set is ABA59548.1,
`Val` against pyruvic acid.

**The second row is the one to remember.** `aaSMILES.txt` ends each line with a `#`
comment that usually holds a substrate name, and using it as a name table yields 18
disagreements of which exactly one is real. The others are a hyphen
(`4-hydroxy-phenylglycine` against `4-hydroxyphenylglycine`) or a note the
maintainer left themselves: `salicylic acid (not in norine yet)`, `actually aile,
allo-isoleucine`. A field that usually parses is the most dangerous kind, because
nothing fails — it just quietly reports 17 findings that are not there. My own
first pass through this row did exactly that, and reported a stratified agreement
difference (0.875 experimental against 0.930 sequence-predicted) that **does not
survive controlling for substrate class**: the experimental stratum simply holds
more non-proteinogenic substrates, which are the ones the bridge cannot resolve.
That reading is withdrawn, not recorded as a finding.

Everything else about the row fails for an adjacent reason. The A-domain's position
in its protein — which the catalogued R3 needs to extract specificity residues — is
`-1/-1` for **2,852 of 2,901** domains, and MIBiG ships no sequences. The substrate
labels themselves are **71.6% inference-only and 10.2% unsupported**, with 1,049
naming "Sequence-based prediction", so grading an agent against them grades
agreement with a predictor. The catalogued proteinogenic/non-proteinogenic
instantiation axis gives **1.000 agreement with zero disagreements** on one side and
**83 of 112 unmappable** on the other. And the bridge-free route — group A-domains
by identical signature, then ask whether MIBiG's own labels agree — collapses as
well: the 34-residue extended signature is **unique for all 304**, grouping
nothing, while the 10-residue signature gives 33 groups of two or more whose 10
"inconsistencies" are mostly granularity (`threonine` against `allo-threonine`,
`glutamic acid` against `d-glutamic acid branched`).

So every rung available here is near-constant or decided by a naming artefact.
That is a vocabulary quiz with a hidden key, and the standing rule — no human
judgement anywhere, and a question-answer pair that depends on one comes out —
retires the row rather than dressing it up. Counts recorded in `docs/catalog.md`.

### The rule this adds, which binds every S2 campaign

`aaSMILES.txt` is antiSMASH package data under **AGPL-3.0-or-later**.
`stachelhaus/1.1/signatures.tsv` arrives from the antiSMASH database download with
**no stated licence at all**. So neither can be shipped in a campaign:
**you cannot redistribute data whose licence you cannot name**, and pulling AGPL
material into `inputs/` is not a thing to do by accident. An S2 campaign reads such
a resource from the pinned image at a declared path with a recorded sha256, and
redistributes nothing — which also means the tool surfaces need reading as
*installation* requirements, not tool-execution requirements. A campaign can need
S2 for a pinned table and never invoke a tool.

## The measurement environment was reading the ambient shell

Refreshing `image/thread_invariance.json` for this commit turned antiSMASH from
PASS to FAIL — both invocations, 4/4 runs, `RuntimeError: Modules failing
prerequisites`. The databases were present and the pin was right. The cause:
**the invariance runner invoked `{bin}/tool` by absolute path but never put the
tool prefix on the subprocess `PATH`.** antiSMASH shells out to nine helpers —
hmmscan, hmmsearch, hmmpress, hmmpfam2, blastp, makeblastdb, diamond, prodigal and
FastTree — and resolves every one of them on `PATH`, so the measurement was reading
whatever the measuring shell happened to export. On the shell that produced the
original 8/8 it was exported; on a clean one it was not, and
`antismash --check-prereqs` reports **44 prerequisite failures** across those nine,
every one of which is sitting in the prefix.

The verdicts were right and the way they were obtained was not, which is the worse
of the two failures: a determinism measurement that depends on the ambient
environment is not a measurement of the image. `_image_path()` now builds the
subprocess `PATH` as the ambient one with the prefix **appended** — last, matching
the image's own `PATH`, so the tool environment still cannot shadow the system
interpreter this suite runs under, which is the other half of a lesson this project
has already paid for once. A test pins both properties.

### And a witness assertion that was a 3% coin flip

The same pass caught `test_diamond_no_reorder_is_still_the_witness_it_is_declared_to_be`
failing on `assert not witness["repeatable_at_fixed_threads"]`. The witness
documents a per-run reordering, so every claim about it is a claim about a sample:
over 30 measurements the two single-thread repeats agreed **30/30** and the two
eight-thread repeats agreed **0/30** — but one four-run draw can still land on two
matching eight-thread runs, and that is what happened. The property under test is
that the eight-thread order is unspecified, which needs two distinct outputs
*somewhere* in the sample, not in one particular pair. The test now draws again
rather than let a coin decide a verdict. **A test that fails 3% of the time is not
a strict test; it is a test whose result is partly noise**, and on a suite whose
whole purpose is to certify determinism that is the one defect that cannot stand.

## Campaign 18: three records of one peptide, and a field with four shapes

`ripp-precursor-mibig-4_0-01` is T-L3-4, the third tool-running campaign and the
second S1. A RiPP is made from a precursor protein carrying a leader, a core and
sometimes a follower; the core becomes the product. MIBiG records the precursor's
gene, the core's sequence and — for some entries — a `leader_cleavage_location`.
The translation itself is in the antiSMASH-processed reference record for the same
cluster. So three independent records describe one peptide, and whether they agree
is entirely checkable.

**Half the catalogued ladder does not survive, and the reason is gradability
rather than cost.** The catalogued R3 was "the core peptide sequence and the
cleavage motif". The core peptide sequence *is* a MIBiG field: grading it grades
transcription, and the campaign must ship the field for anything else to be
solvable, which makes it a free component twice over. The cleavage motif is not in
MIBiG at all. The catalogued R4 was "the leader is altered — predict the new
core", and that has **no oracle**: no pinned tool computes it, the answer is a
wet-lab fact about a peptide nobody has made, and scoring it against a curator's
intuition is the human judgement this benchmark excludes from grading. R4 keeps
the counterfactual's intent — what does the leader carry — and makes it
measurable, by removing each declared leader and asking the tool what survives.

### The field has four JSON shapes and one of them is not a peptide

| shape of `core_sequence` | records |
|---|---|
| list of one string | 35 |
| bare string | 28 |
| list of several strings | 2 |
| **the printed form of a list, as a string** | **2** |

The last is literally `"['GGAGHVPEYFVGIGTPISFYG']"` — brackets and quotes inside
the value, from a serialiser that printed a list instead of emitting one. Three
such records exist across the whole class; two are in the corpus.

`core_sequence[0]` — the natural thing to write against the commonest shape — is
correct for the list shapes and takes **one residue** where the field is a bare
string. It raises no error on any record, and a one-residue core still localises;
it just localises many times over, so the report comes back saying *ambiguous*
rather than *I mis-parsed this*. That is the campaign's competence probe
(`stub-firstshape`), and it clears R1 and fails R2.

Case is the other half of the same problem and is handled differently on purpose.
Seven corpus cores are recorded in lower case against upper-case translations.
Folding case cannot invent a match that is not there, so it is a **normalisation**
and belongs in the baseline policy; a character that is not an amino-acid letter
*in either case* is a malformed core. Conflating the two reports a case difference
as a data defect — which is what the first draft of this campaign did, and the
build's own policy-sweep assertion caught it: with the fold hard-coded into the
comparison, the strict policy and the baseline produced identical censuses and the
build refused gold, saying two of three declared policies were the same policy
under two names.

### Localisation, and the six entries that pasted in the whole precursor

Under the baseline policy, over 67 precursor records from 54 clusters:

| verdict | records |
|---|---|
| core at the C terminus | 49 |
| core internal | 9 |
| **core is the whole precursor** | **6** |
| malformed core | 2 |
| core ambiguous | 1 |

Leader lengths run **1 to 85**, median 23. The minimum of 1 is real: bottromycin's
core sits at the N terminus and its leader is the initiator methionine alone.

Six entries put the **entire precursor**, leader included, into `core_sequence`.
That localises perfectly at offset zero, so a reader that only asks "is the core
findable" reports a clean record with a zero-length leader. The campaign gives it
its own verdict and nulls the derived coordinates, and the same is true of the
tandem-repeat cyanobactin precursor whose core occurs twice: *where is the core*
has no answer, and filling it with the first match would be a declared choice the
campaign has not declared.

### The coordinate convention is elected, not assumed

Nothing in MIBiG says whether `leader_cleavage_location` is 0-based or 1-based, or
whether `from` or `to` marks the core's first residue. All four readings are scored
against the measured position:

| reading | agrees | disagrees | underivable |
|---|---|---|---|
| **`to` as the 0-based core start** | **27** | **0** | 3 |
| `from` as the 0-based core start | 4 | 23 | 3 |
| `to` as the 1-based core start | 0 | 27 | 3 |
| `from` as the 1-based core start | 0 | 27 | 3 |

So `to` is the 0-based index of the first core residue, and on **every** record
where both a coordinate and a position exist, the two agree. The three that do not
are the three with no position: the tandem repeat, and two whose `core_sequence` is
the whole precursor while a cleavage site is *still* declared, with `from == to` —
which under no convention describes a bond. The four records where `from` also
agrees are exactly those `from == to` cases, which is why `from` scores 4 rather
than 0.

This is the same move T-L3-1 had to make for MIBiG's locus bounds, and the build
refuses gold if the top two readings tie, if none agrees anywhere, or if none
disagrees anywhere — any of which would make "elected" a tie-break dressed as a
finding.

### No Pfam family marks the boundary, and the leader is why

HMMER 3.4 against a declared eight-model Pfam 35.0 panel, fetched from the image
with `hmmfetch`. The panel was fixed at build time by scanning the corpus against
the whole of Pfam-A once — 27 seconds for 67 sequences of fifty-odd residues — and
listing every family that came back at the gathering threshold: `DUF5837`,
`DUF5840`, `DUF5973`, `Gallidermin`, `Lantibiotic_a`, `Nif11`, `RamS`,
`Thiopep_pre`.

**The panel hits 12 of 67 precursors.** Of the 10 that can be placed against a
measured boundary, **8 span the cleavage site** and 2 sit in the leader alone.
**None is core-specific.** The offset between a hit's envelope end and the leader's
length runs −21 to +30 and is **never zero** — so no family stops where the leader
does, and "use the domain call as the boundary" is not available. Five of the 12
envelopes cover the precursor end to end.

Then R4 asks the tool directly. Remove each record's declared leader, re-scan with
the same panel at the same threshold:

| outcome | records |
|---|---|
| **call lost entirely** | **7** |
| call retained | 3 |
| never recognised, still not (control) | 48 |
| no located core, not applicable | 9 |

Seven of the ten calls are carried by the leader. The three that hold are
`DUF5973` on cinA and `Lantibiotic_a` twice — families whose envelope lies mostly
in the core. Nothing is gained, which the vocabulary allows for and the corpus
declines to supply.

### Two measurements about the resources themselves

**The two exclusion clauses coincide exactly.** Of 447 ribosomal entries, 100 name
both a precursor gene and a core sequence. 46 of those are `retired`; 46 are absent
from the antiSMASH reference set; and they are **the same 46**, with no exceptions
in either direction. Applied in sequence the second clause would have reported
zero and the coincidence would have been invisible — a fact about the evaluation
order, not about the resources — so the campaign counts the two independently and
grades their overlap. The audit runs on what survives both: 67 records, 54
clusters, 30 declaring a cleavage coordinate.

**The reference records shipped are those for every entry of the class**, 363
files and 3.4 MB compressed, not the 54 the audit uses. Shipping only the entries
that pass the reference clause would make the count of those that do not readable
off a directory listing instead of computed from the two resources. That is a
**free component**, and it is the sixth one caught in this project.

Five echoes are reported and not graded for the same reason — `mibig_release`
comes from the input provenance, and `biosynthetic_class`, `baseline_policy`,
`threshold` and `panel_size` from the pinned rules table. Nothing in the world
would have to change for their values to change, only the campaign's own
declaration.

### What broke on the way

**The GenBank feature table does not start at column six.** The first parser ended
the feature table at the first line that did not begin with six spaces — and a
feature key begins at column *five*, so the table ended at the first feature and
the campaign found zero translations in every record. The stop condition is now
the first line that starts in column 0, which is what GenBank actually guarantees,
and there is a test that writes a minimal record and asserts the one translation
comes back.

**Factoring, not duplicating.** `Resource` — the declaration for a file read from
the image and not redistributed — existed twice already, in `ec_domain_audit` and
in `gcf_cutoff`, with the second supporting a globbed path and the first a
`verify()`. A third copy would have been the smell this project keeps finding in
other people's data, so both moved to `npbench_c.tools.resources`, which carries
the union: globbing, with the fingerprint deciding whether the match was the
declared one. `resolve()` now raises `ResourceError` rather than a template's own
error type, which is one line in the EC-audit tests, and both campaigns' suites
pass unchanged otherwise.

### The sweep, and what it cost

Gold builds in **7.7 seconds**, including both HMMER scans. The stub sweep —
six systems, four rungs, three repeats, 72 runs — takes **7m51s** with all five
gates passing and cold depths **1 / 2 / 3 / 4 / 1 / 0**:

| system | mode | cold depth |
|---|---|---|
| stub-reader | tooled | 1 |
| stub-localise | tooled | 2 |
| stub-evidence | tooled | 3 |
| stub-complete | tooled | 4 |
| **stub-firstshape** | tooled | **1** |
| stub-noncompute | no_tool | 0 |

Per-rung clear rates 1.00 / 0.60 / 0.40 / 0.20, cold mean score 0.55, and the
no-tool ablation at 0.0000 against a chance floor of 0.2502.

The residual cost is the same one open item 8 names, in its cheaper form: the tool
memo works, but each of the 72 runs copies 4.3 MB of inputs into a fresh sandbox
and re-reads a 3,013-entry annotation archive and several hundred GenBank records
to rebuild the corpus. That is about 6.5 seconds per run of pure re-parsing, and
it is the whole of the 7m51s.

## Campaign 19: the tool's output is not a domain list yet

`arch-nrps-mibig-4_0-01` is T-L1-2, the fourth tool-running campaign and the
third S1. An NRPS is one protein carrying a chain of catalytic domains, and its
architecture is the ordered list of them. Computing that from a pinned HMMER
against a pinned Pfam release sounds like a lookup. It is not: the answer is a
function of three choices nobody writes down, and this campaign declares all
three and measures what each one costs.

**Two of the four catalogued rungs do not survive, and a third is rewritten.**
The catalogued R3 claimed "which domain is catalytically essential, closed enum",
which is a literature claim with no computable oracle — the same refusal T-L3-4's
"predict the new core" earned. The catalogued R4's "predict the lost function"
goes with it. The catalogued R1 asked for a "non-overlapping" domain table, and
that is not what the tool produces.

### Hits overlap, and the resolution is the campaign

Measured over the corpus: **12,975 raw hits, 2,124 overlapping pairs, 315 of
1,136 genes carrying at least one overlap.** Nested Pfam families hit the same
residues — a ketosynthase core inside a thiolase model, five methyltransferase
families against each other, `adh_short` against `KR`. Resolving them drops
**1,031 hits** and changes the architecture on **314 genes**. A report that lists
every hit as a consecutive domain has invented a protein, and nothing in it says
so. That is the campaign's competence probe, and it is the mistake the catalogued
design made.

**Pfam's own rule is not available, which is worth recording for every future
campaign.** Pfam resolves intra-clan overlaps by keeping the best-scoring hit.
The pinned library — antiSMASH's copy of Pfam 35.0, 19,632 models — carries
**zero `CL` lines**. So the rule here is score and geometry, declared as this
campaign's own rather than borrowed and misattributed.

**An overlap *threshold* is not a subject.** Swept from 0.8 down to 0.0 the
architecture moves by at most a handful of genes, because the overlaps these
proteins produce are near-total rather than partial. A grid would be five points
reporting one answer, so the fraction is declared once at 0.5 and the axis R4
varies is whether overlaps are resolved at all.

### What the curated side can and cannot support

MIBiG 4.0 carries **4,034 curated modules over 667 entries**, with per-module gene
assignment and typed domain blocks. Two of its fields are unusable and the
campaign says so with a number instead of working around them:

- **Domain coordinates are placeholders**: 5,352 of 5,446 typed blocks carry
  `{-1, -1}`. A coordinate-level reconciliation would be computed over a
  sixtieth of the record and reported as though it covered all of it.
- **The `active` flag is unstated** on 2,025 of 4,034 modules (1,945 true, 64
  false). So "active modules" is not a quantity this database supports.

So counts are compared, not coordinates. The selection takes genes whose curated
modules are **all** of the declared type — a gene carrying both NRPS and PKS
modules draws its architecture from two vocabularies and belongs to neither
instantiation — on an entry the status rule admits, with a shipped translation:
**1,136 `nrps-type1` genes over 512 clusters**, 3,440,119 residues, median 2,549
aa, longest 16,367. **11,944 domains, 458 distinct architectures**, 0 to 61
domains per gene.

### Curation and the tool agree on one domain type and not the other

| comparison | agrees | computed more | computed fewer | not curated |
|---|---|---|---|---|
| adenylation (`AMP-binding`) | **1,059** | 70 | 7 | 0 |
| condensation (`Condensation`) | **157** | 122 | 7 | **850** |

The adenylation delta is a distribution, not an offset: −6, −2, −1 (×5), 0
(×1,059), +1 (×60), +2 (×9), +4. A curator records one adenylation domain per
module by construction, so that column is reliable and agreement is 93%. They
record a condensation domain only sometimes — 286 of 1,136 genes — and when they
do, agreement is **55%**. The campaign reports `not_curated` for the other 850
rather than a disagreement with zero, because there is nothing there to disagree
with.

The module decomposition is a third, different claim: counting the declared
pattern's non-overlapping occurrences in the architecture finds **1,881 complete
modules**, agreeing with the curated module count on 587 genes and falling short
on **541**. A gene can carry the right number of adenylation domains and not
assemble into complete modules at all.

### Which Pfam family *is* the adenylation domain

| mapping | agrees |
|---|---|
| `AMP-binding` | **1,059** |
| `AMP-binding_C` | 814 |
| both | **153** |

Pfam splits the adenylation domain into an N-terminal and a C-terminal family, so
counting both double-counts every module and agreement collapses by a factor of
seven. This is the axis the catalogued design did not have and the data demands,
and the build refuses to emit gold unless the best and worst declared mappings
differ by at least a factor of two.

### The cutoff, and the one flag that is not a substitute

Sixteen cells, one per (cutoff, resolution, coordinate span), against the
baseline:

| cutoff | genes moved, resolved | genes moved, unresolved | agreement, resolved | agreement, unresolved |
|---|---|---|---|---|
| `--cut_ga` (baseline) | 0 | 314 | 1,059 | 1,059 |
| `--cut_nc` | 4 | 317 | 1,058 | 1,058 |
| `--cut_tc` | 4 | 314 | 1,059 | 1,059 |
| **`-E 1e-5`** | **372** | **549** | **1,032** | **1,026** |

The three Pfam-curated cutoffs are near-interchangeable here — the same result
the EC-domain audit got on a different set — and the bare E-value is not. It adds
816 domains, moves a third of the corpus and costs 27 agreements. That contrast
is why it is in the grid: a per-model curated threshold and a global E-value are
not substitutes, and the E-value is the familiar flag.

The bottom-right cell is also where a defect in this campaign's own code showed
up. `reconcile` originally took the architecture table as an argument and then
**re-resolved the hits at the baseline**, so the sweep's agreement column ignored
two of its three axes. On the curated cutoffs the numbers coincided, which is why
it survived a reading; the fix moved the two unresolved E-value cells from 1,032
to **1,026**, because at that cutoff duplicate overlapping adenylation hits
inflate the count on six genes and the old code was quietly reporting the
resolved figure instead. A graded value computed by a function that ignores two
of its inputs is wrong even where the output happens to be right.

### The catalogued R4 was built, measured, and dropped

Removing each gene's first adenylation domain and re-scanning leaves the
architecture minus exactly that domain on **1,112 of the 1,126 genes it applies
to**, and the adenylation count falls by exactly one on 1,125. So "one domain
fewer" is a free answer and the rung would distinguish almost nothing. It is
recorded in `excluded_from_grading` with those numbers rather than silently
omitted, because the measurement is the reason.

### What broke on the way, and it broke twice

**A Pfam family name can contain the architecture separator.** The architecture
is a hyphen-joined family string, and `AMP-binding` holds a hyphen — so splitting
it back shreds every family name. The module decomposition did exactly that, and
the build's own refusal caught it: *"the declared module pattern is found on no
gene, so the decomposition is empty everywhere."* The decomposition now reads the
resolved hits; the string is for reading, not for parsing back. Then the test
written to check the fix made the same mistake on its first draft, asserting the
domain count against `architecture.split("-")`. Both have tests now, and the
second one asserts on real gold that splitting over-counts.

**And the declared module pattern is four families long for a module with three
domains.** `Condensation-AMP-binding-AMP-binding_C-PP-binding` occurs 1,881 times
and is by a wide margin the commonest four-gram in the corpus; the three-family
pattern a domain-level reading would write down occurs on a tenth as many genes.
The same fact the mapping axis measures, arriving a second time through a
different door.

**`Resource` was factored before this campaign, not during it.** It now carries
globbed paths and `verify()` in `npbench_c.tools.resources`, which this campaign
is the third to use.

### The sweep, and what it cost

Gold builds in about **two minutes**, which is four panel scans over 3.44 M
residues plus a sixteen-cell sweep in Python. The stub sweep — six systems, four
rungs, three repeats, 72 runs — passes all five gates with cold depths
**1 / 2 / 3 / 4 / 1 / 0**:

| system | mode | cold depth |
|---|---|---|
| stub-reader | tooled | 1 |
| stub-architect | tooled | 2 |
| stub-reconcile | tooled | 3 |
| stub-complete | tooled | 4 |
| **stub-unresolved** | tooled | **1** |
| stub-noncompute | no_tool | 0 |

Per-rung clear rates 1.00 / 0.60 / 0.40 / 0.20, cold mean score 0.55, no-tool
ablation 0.0000 against a chance floor of 0.2502. `stub-unresolved` is the probe
that matters: it runs the tool correctly and reports every hit as a domain, which
clears R1 and fails R2 — the catalogued design's own mistake, scored.

The readiness gate is slower here than on any earlier campaign, because two of
its checks regenerate gold and gold costs four tool scans. That is the gate
working as intended rather than a cost to optimise away.

### Why HMMER and not antiSMASH's own domain calls

The tool registry already declares an antiSMASH invocation,
`default_modules_domains`, whose projection is the ordered NRPS/PKS domain
architecture per CDS — which looks like exactly this campaign's answer, computed
by one call instead of four. It is not used, for two reasons worth stating.

The catalogued row is a **Pfam** row, and antiSMASH's NRPS/PKS domain calls come
from its own curated profile library and its own module-assembly logic, not from
Pfam. Grading those would be grading antiSMASH's opinion of a module, which is a
different and much more derived claim than "which Pfam families does this
sequence carry, and where" — and it would make the reconciliation against MIBiG
circular in part, because MIBiG's own processed records were produced by
antiSMASH.

And it would hide the campaign's subject. antiSMASH resolves overlaps and
assembles modules internally and reports the result; the three choices this
campaign declares and measures would all be made inside the tool, unstated. The
whole finding — that the architecture is a function of a resolution rule, a
coordinate span and a cutoff — is only visible from the raw hit table.

The antiSMASH route remains the right one for a campaign that wants to grade
antiSMASH's module assembly. That is a different row.

## Campaign 20: the second instantiation, and what it is for

`arch-pks-mibig-4_0-01` completes T-L1-2's catalogued ×2 axis. It is the same
template as the NRPS campaign and differs only in its declared rules file — the
module family, the Pfam panel, the two reconciled domain types, the module
pattern and the strata. Oracle 1.0 at depth 4, gate 22/0/5.

**464 `pks` genes over 162 clusters**, 1,460,845 residues, median 2,347 aa.
**5,954 domains, 276 distinct architectures**, 1 to 50 per gene.

An instantiation axis is only worth having if the second instance measures
something the first cannot, and this one does.

### The curated record encodes cis-AT versus trans-AT, exactly

| module type | modules | ketosynthase recorded | acyltransferase recorded |
|---|---|---|---|
| `pks-modular` | 612 | 612 | **612 (100%)** |
| `pks-modular-starter` | 17 | 0 | 17 |
| `pks-trans-at` | 471 | 471 | **0 (0%)** |
| `pks-trans-at-starter` | 33 | 0 | 0 |

A trans-AT polyketide synthase does not encode its own acyltransferase — the
enzyme acts in trans from a separate protein — and MIBiG records that in the
module type with no exceptions either way. So on the trans-AT side a curated
count of **zero** is a claim the computed architecture can agree or disagree
with, while elsewhere a zero is a silence. The campaign declares two strata over
the accepted types and grades the secondary reconciliation per stratum:

| stratum | genes | agrees | computed more | computed fewer | acyltransferases found |
|---|---|---|---|---|---|
| `cis_at` | 318 | 298 | 10 | 10 | 546 |
| **`trans_at`** | **128** | **128** | **0** | **0** | **0** |
| `mixed` | 18 | 9 | 6 | 3 | 26 |

**The tool finds no acyltransferase in any of the 128 trans-AT genes.** That is a
negative control the data supplies rather than one the campaign plants, and it is
a measurement: if HMMER had found some, the finding would have been about the
tool.

And which strata assert absence is **read from the per-module-type census rather
than asserted in the rules**. A MIBiG release that started recording trans-AT
acyltransferases would change the verdict instead of contradicting a hardcoded
fact.

The module decomposition separates the same strata through a second, independent
measure. The declared cis-AT module core is
`ketoacyl-synt-Ketoacyl-synt_C-KAsynt_C_assoc-Acyl_transf_1` — four Pfam families
for two catalytic domains — and it occurs 540 times, agreeing with the curated
module count on 292 genes. In the trans-AT stratum it occurs **zero** times,
because the pattern includes an acyltransferase those genes do not have.

### The sweeps

Both instantiations pass all five sweep gates with cold depths
**1 / 2 / 3 / 4 / 1 / 0** and identical per-rung clear rates (1.00 / 0.60 / 0.40
/ 0.20), cold mean score 0.55, no-tool ablation 0.0000 against a chance floor of
0.2502. One stub module serves both campaigns: every level reads the family, the
panel, the strata and the reconciled domain types out of the rules file, so there
is nothing family-specific in it — including `noncompute`, the ablation that reads
nothing and therefore must not name either family's vocabulary.

### The split family problem is worse here

| mapping for the ketosynthase | agrees |
|---|---|
| `ketoacyl-synt` | **398** |
| `Ketoacyl-synt_C` | 393 |
| `KAsynt_C_assoc` | 389 |
| **all three** | **13** |

Pfam splits the ketosynthase across **three** consecutive families, so counting
all of them triples every module and agreement falls by a factor of thirty — against
seven for the adenylation domain's two halves. The three single-family readings
agree closely but not identically, which is itself worth having: they are three
different models of the same region and they disagree on 9 to 55 genes each.

### And the overlaps are much heavier

| | genes | overlapping pairs | genes with an overlap |
|---|---|---|---|
| NRPS | 1,136 | 2,124 | 315 (**28%**) |
| PKS | 464 | 3,002 | 399 (**86%**) |

The reductase families a PKS carries — `KR`, `adh_short`, `adh_short_C2`,
`ADH_zinc_N`, `3Beta_HSD`, `Epimerase` — overlap each other far more densely than
the NRPS panel's do. So the resolution policy, which the catalogued design assumed
away, matters on **86% of this corpus**. The policy sweep says the same thing:
not resolving moves 399 of 464 genes, and `--cut_tc` moves 58 here against 4 in
the NRPS corpus.

### What the template had to grow, and the three defects it exposed

The instantiation axis is a declared rules file rather than a forked template, so
five fields became declarations: the module **family** (a set of accepted types,
because a modular PKS gene carries a loading module of its own type), the two
**reconciled domain types** (`ks_domain`/`at_domain` here, `a_domain`/`c_domain`
there), the **strata**, the module **pattern** and the **panel**. Both campaigns
name the same commands, the same stub module and the same pinned-table path; only
`template_params.rules_source` differs.

Generalising a working template broke it three times, and each break is worth
recording because none of them would have shown up in an oracle run:

- **The secondary domain's mapping key drifted from its name.** `reconcile` looks
  the mapping up by the declared secondary domain — `c_domain` — and the NRPS
  rules file still called it `condensation_domain`. The lookup returned nothing,
  so the secondary reconciliation computed **zero for every gene**: 157 agreements
  became 0 and 286 genes became `computed_fewer`. Every build assertion still
  passed, because a uniformly wrong column is still a populated column. The keying
  is now checked at load time, along with the baseline mapping not being the
  secondary domain's own.

- **`"mixed"` was a stratum the census had no room for.** A gene spanning two
  strata belongs to neither, so `stratum()` returns `"mixed"` — which the
  by-stratum census, initialised from the declared strata alone, raised a
  `KeyError` on. The NRPS family has one stratum and never produces it, so this
  surfaced on the PKS build's first run.

- **A build assertion was disabled by the fix for the one above.** The check "a
  single-stratum corpus whose secondary domain is always curated has one kind of
  row" counted the keys of the by-stratum census — which now always carries an
  empty `"mixed"` row, so the count was never below two and the check never fired.
  It counts **populated** strata now, and a test pins that.

Two shape mismatches between gold and the measurer also had to be settled, and
the honest fix in both cases was to close a vocabulary rather than teach the
measurer to drop things: the selection census now names all three of its drop
clauses whether or not each removed anything, so the census's shape does not
depend on the data.

## Campaign 21: the flag matters more than the assembly

`genecall-scoelicolor-01` is T-L1-1, the fifth tool-running campaign and the
fourth S1. A gene caller reads a contig and predicts coding sequences; cut the
same contig into pieces and run it again and the answer changes. How much of the
change is geometry and how much is the tool reacting to its input is the
campaign, and on this chromosome the answer is not the one the catalogued design
implies.

Pinned Prodigal 2.6.3 over the complete *Streptomyces coelicolor* A3(2)
chromosome (AL645882.2, 8,667,507 bp, GC 0.721). Oracle 1.0 at depth 4, gate
22/0/5. **7,738 calls** on the whole contig — 3,969 forward, 3,769 reverse —
covering **88.5%** of the chromosome, median call 855 bp.

**This is the first tool-running campaign here with no external resource to
pin.** Every other one reads a reference database from the image and declares it
by path, provider, version and digest. An ab initio caller reads only the
sequence it is given, so the version and the mode are the whole of the pin, and
`constraints.external_resources` is an empty list the gate reports as such.

### Nothing annotated is shipped, and that is the design

The detection campaign uses this same accession and ships the annotated GenBank,
because there the annotation is what antiSMASH reads. Here it is what the tool is
asked to reproduce, so the input is the bases alone: a 2.44 MB gzipped FASTA, and
the **7,825 CDS features of the source record are discarded**, with the count
recorded in the provenance so what was thrown away is stated rather than implied.

The curated annotation is also refused as a reference, with the reason in
`excluded_from_grading`: using it would make the campaign about calling
*accuracy* against an annotation that is itself largely pipeline output, which is
a much more derived claim than "what does this pinned caller do to this contig,
and what does cutting the contig do to that".

The provenance also does **not** state the sequence's length or digest, both of
which R1 grades. Writing them there would have made two free components, and the
digest is over the bases rather than the file precisely so that computing it is
where the declared normalisation gets checked — a FASTA rewrapped at another
width must fingerprint the same.

### The fragmentation is a rule, and the coordinates are the trap

The breakpoints are arithmetic: for a contig of length L cut into n fragments,
cut at `round(k·L/n)`. The campaign ships one unfragmented file and the agent
cuts it, because shipping pre-cut fragments would put each fragment's offset in
the input — and the offsets *are* the coordinate contract.

The caller reports 1-based inclusive positions **on the fragment**. A parent
coordinate is that plus the fragment's 0-based offset. Skip the lift and almost
every whole-contig call looks lost and almost every fragment call looks new, with
counts in the right order of magnitude and tables the right shape. Measured, the
unlifted comparison reports **more than twenty times** the real losses and gains,
with an identical set of keys and the same vocabulary. That is the competence
probe, and it clears R1 and R2 — which concern the whole contig and need no lift
— before failing at R3, where the two coordinate systems first have to meet.

The two conventions also meet once inside the engine, in one line: a cut at
0-based offset `c` falls strictly inside a 1-based inclusive call when
`start ≤ c < end`. Off by one either way and a gene spanning a breakpoint is
classified as merely near one.

### Fragmentation at the declared baseline

Cut into 50: **7,738 → 7,779 calls, 7,686 retained, 52 lost, 93 new.** The whole
contig has **no partial call at all** — a complete record has no edge for a call
to run off — and the fragmented run has 42 running off the left of their fragment
and 44 off the right.

The 52 losses, classified against the breakpoints:

| verdict | calls |
|---|---|
| crosses a cut | **43** |
| near a cut (≤ 5 kb) | 4 |
| **distant from every cut** | **5** |

The first is geometry and it dominates: a gene spanning a breakpoint cannot be
called whole. The third is the one that matters, because no cut can explain it.

### Only the mode that trains loses calls a cut cannot explain

Prodigal's single mode **trains on the input it is given**; its meta mode scores
against pre-trained models. So fragmenting the contig changes the training set in
one mode and not the other, and the grid separates the two effects:

| mode | fragments | lost | gained | crosses | near | **distant** |
|---|---|---|---|---|---|---|
| single | 2 | 3 | 3 | 0 | 1 | **2** |
| single | 10 | 9 | 17 | 7 | 2 | 0 |
| single | 50 | 52 | 93 | 43 | 4 | **5** |
| meta | 2 | 2 | 2 | 0 | 2 | **0** |
| meta | 10 | 10 | 16 | 7 | 3 | **0** |
| meta | 50 | 48 | 88 | 43 | 5 | **0** |

**The distant column is empty in meta mode at every fragmentation level and
populated in single mode.** The campaign does not assert that; it reports it per
cell, and the grid compares each mode against its *own* unfragmented run — not
against the baseline mode's, which would add the mode difference to the
fragmentation difference and report the sum as fragmentation.

### And the headline

The two modes on **identical input bytes**: 7,475 shared calls, 263 only in
single, 283 only in meta. **546 calls move because of one command-line flag,
against 145 that move from cutting the chromosome into fifty pieces.**

So a campaign that looks like it is about assembly quality is, on this genome,
mostly about a flag — which is why `environment.md` already warned that
Prodigal's mode must be pinned, and why the mode is a declared axis here rather
than a default.

### Two things worth recording about the tool

**The GFF is byte-identical on re-run.** Prodigal stamps no date and no command
line into its output, unlike HMMER, which needs a declared normalisation for
exactly that. So this campaign declares none, and that is a measurement rather
than an assumption.

**`retained_but_truncated` is zero everywhere, and is kept.** A call that crosses
a cut gets clipped, so it changes coordinates and appears as a loss plus a gain
rather than as a retained call with a changed flag. The quantity is not
structurally impossible — a call ending exactly at a fragment boundary could keep
its coordinates and gain a partial flag — so it is computed and reported inside
the fragmentation claim rather than dropped. On this contig it is zero at every
cell of the grid.

### What the catalogued design had to give up

The catalogued R2 was to be scored by overlap-F1 within a measured tolerance.
Scored exactly instead, for the reason the domain-architecture campaign found:
overlap-F1 absorbs boundary disagreement between two *independent* predictions,
and both sides here come from the same pinned caller over the same bytes. With a
byte-identical GFF there is no spread to take a tolerance from, so any tolerance
would absorb a difference that cannot occur while hiding one that can.

The catalogued R3 — "which genes are truncated at a contig boundary, from a
closed enum" — is built, and the enum is the caller's **own** `partial` flag
rather than one the campaign invents. What is added is the classification of the
losses, because that is what turns the rung from a lookup into a measurement.

### The sweep, and a probe that gets further than the others

Gold builds in about **two minutes** — eight caller runs over the chromosome at
roughly fifteen seconds each, with no database to load. The stub sweep passes all
five gates with cold depths **1 / 2 / 3 / 4 / 2 / 0** and per-rung clear rates
1.00 / 0.80 / 0.40 / 0.20, cold mean score 0.60, no-tool ablation 0.0000 against
a chance floor of 0.2502.

| system | mode | cold depth |
|---|---|---|
| stub-reader | tooled | 1 |
| stub-caller | tooled | 2 |
| stub-reconcile | tooled | 3 |
| stub-complete | tooled | 4 |
| **stub-unlifted** | tooled | **2** |
| stub-noncompute | no_tool | 0 |

The probe reaching **depth 2** rather than 1 is the point, and it is the first
one here that does. Every earlier competence probe failed at the rung after R1,
because its defect was visible in the first thing it reported. This one runs the
tool correctly, reports the whole contig correctly, and clears two rungs — the
coordinate lift it skipped does not exist yet at R1 or R2, because both are about
the unfragmented contig. It fails at R3, the first rung where fragment
coordinates and parent coordinates have to be compared.

That is worth having in the benchmark: a defect whose depth tells you *where* the
convention was missed rather than only that something was wrong.

## The guards nobody was guarding

The previous campaign ended with a build assertion that had stopped asserting: a
check for a single-stratum corpus counted the keys of a census that always
carried one extra, so the condition could never hold. Its test caught it — but
only because a test existed. That raised the obvious question, which turned out
to have an uncomfortable answer.

**Of 97 refusals across the nine templates that have a `_assert_has_content`, 43
were named by no test at all.** Three templates — `chembl_selectivity`,
`chemical_space`, `kinetics_consistency` — looked like they had none, though two
of those turned out to be an artefact of how the check first read the tests.

A `_assert_has_content` is the only thing between a campaign whose rungs measure
nothing and a campaign that ships. Every refusal in it is a claim about what the
data has to look like, and a claim no test exercises is a claim nobody has
checked — the same defect this project keeps finding in other people's
declarations, turned on its own.

### Two failure modes, and they need different mechanisms

- **A refusal is added without a test.** Nothing notices. The build still passes
  on real gold, the oracle still scores 1.0, and the refusal may be unreachable
  or inverted. Nothing in the project checked for this.
- **A refusal's condition becomes unreachable while its test still exists.** The
  test notices, because `pytest.raises` fails with "DID NOT RAISE". This is what
  caught the strata check.

`tests/unit/test_build_refusals.py` is the check for the first, and the two
together cover both directions. For every `raise` inside a template's
`_assert_has_content`, it requires the template's campaign test file to carry a
pattern that matches that refusal's message — which is why every refusal message
in this project ends up with a distinctive phrase. **The phrase is the contract
between the refusal and its test.**

**It is static and order-independent on purpose.** Tracing the suite and
asserting every `raise` line executed would be more precise, and it was the first
design. It was dropped: it would depend on the whole suite running in one process
and on test order, and this project runs its tests under `pytest-randomly`.

Three things the check needed before it was honest about its own subject:

- **It reads the syntax, not a regex over `match=`.** Two idioms are in use here:
  `pytest.raises(Err, match="...")` directly, and a parametrised refusal test
  whose `(mutate, message)` pairs put the pattern in the decorator and a variable
  at the call site. A regex sees only the first, and on the first run it reported
  `chemical_space` and `kinetics_consistency` as having zero coverage when they
  were largely covered. Two of the three "templates with nothing" were my own
  measurement error, and the correction is in the numbers above.
- **It resolves a message held in a variable.** Several templates loop over
  `(field, complaint)` pairs and raise the complaint, so the message is a `Name`
  at the raise site. The literals are still in the function and a refusal written
  that way is no less a refusal.
- **It takes only the tuple position the name unpacks from.** Collecting every
  literal in such a loop's iterable also collects the field names sitting beside
  the messages, and the check went looking for tests naming
  `entries_excluded_by_status`.

And it guards itself: one test asserts every template with a
`_assert_has_content` is listed, another that each listed template actually has
refusals, so emptying a template's assertions cannot make the coverage check pass
vacuously.

### It found a dead refusal immediately

`ec_domain_audit` refused gold when "neither the threshold nor the policy axis
moves". That branch was **unreachable**: a flat policy axis means every declared
policy gives the same count, which an earlier refusal already rejects
unconditionally, so the conjunction could never hold. The comment above it
explained the intent clearly and the code did nothing.

It is removed rather than tested, because a refusal that cannot fire is not a
guard, and the asymmetry it was reaching for is real and worth stating where the
code says it: a flat **policy** axis is refused outright, while an inert
**threshold** axis is allowed, because on this set the threshold barely mattering
is a measured robustness result. A test now pins both halves of that — the
refusal fires on a flat policy axis, and gold with a flattened threshold axis is
accepted.

### What closing the other 42 cost

42 new mutation tests across seven templates, each small: doctor one field of
real gold, assert the refusal fires. Two were worth the trouble on their own:

- `bgc_detection`'s extent-vocabulary refusal was first written against a verdict
  term called `contained`, which that campaign's vocabulary does not have — the
  terms are declared in its rules file and are `congruent`,
  `contained_generous`, `contained_much_larger`, `marginal`, `not_detected`,
  `partial`. The test now takes whichever term the vocabulary leads with, so
  renaming one in the rules file cannot quietly turn the test into a different
  test.
- Several refusals sit behind an earlier one that fires first on the obvious
  mutation, so reaching them needs a mutation that satisfies everything above.
  `bgc_detection`'s "no deletion abolishes a region outright" needs one inert
  variant to keep the negative-control refusal happy and the rest gaining as well
  as losing; `ec_domain_audit`'s policy asymmetry needs the threshold axis
  flattened alone. Writing those is where the ordering of a template's refusals
  stops being incidental and becomes part of its contract.

## Campaign 22: breaking the gene on purpose

T-L1-4 is frameshift and pseudogene detection, and its catalogued R3 is "which
locus is a pseudogene and why, closed enum". That claim cannot be graded. The
literature disagrees about which *S. coelicolor* loci are pseudogenes, the
annotation pipelines disagree with the literature and with each other, and there
is no oracle to appeal to. Scoring a system against a curator's reading is the
human judgement this benchmark keeps out of grading, so the catalogued claim is
excluded by name.

What replaces it: the campaign deletes a base itself. Take Prodigal's own call
table on the chromosome, pick the call closest to each declared gene length,
delete one base at its midpoint, run the caller again, lift the coordinates back
past the deletion, and report what became of the gene. The gene is broken by
construction, so the graded question is never whether it is broken. It is what
the caller says about a gene that demonstrably is.

`campaigns/frameshift-scoelicolor-01`: AL645882.2 again, 8,667,507 bp, Prodigal
2.6.3, 7,758 baseline calls, 21 targets from 99 to 22,392 bp. Oracle 1.0 at depth
4, gate 22/0/5. The build is 32 caller runs and takes 12m47s.

### The same tool, the opposite flag

Campaign 21 pins Prodigal's `single` mode on this exact chromosome. This one
pins `meta`, and the two reasons are the same measurement read in opposite
directions.

Single mode trains on the input it is given. In campaign 21 that is the finding:
fragmenting the contig changes the training set, which moves 546 calls against
145 that the cuts themselves explain, and the headline is that the flag matters
more than the assembly. Here the same property is the noise floor. One deleted
base out of 8.67 megabases moves 43 to 57 calls elsewhere on the chromosome in
single mode, and the signature this campaign is trying to see is one call lost
and two gained. The signal does not survive that.

So the campaign reports a ratio rather than an argument. The noise arm re-runs
four of the targets under single mode and gives, per target, how many moved calls
the broken gene itself explains against how many moved elsewhere: 53.0, 19.0,
18.0, 14.3. In meta mode the collateral columns are 0 on 19 of the 21 targets.
The build refuses to emit gold if the retraining mode comes out quieter than the
signal, which is the one case where the choice of baseline would need
re-measuring.

Two campaigns on one tool with opposite pins, each stated as a number, is the
closest this project has come to a result about a tool rather than about data.

### The lift, and a probe that surprised me

Delete a base at position *p* and every coordinate in the mutated run at or after
*p* stands for a parent coordinate one higher. The lift is one addition per call.
I expected skipping it to be catastrophic and obvious, since it displaces roughly
seven thousand calls by one base.

It is catastrophic and it is not obvious. The verdict usually survives. Overlap
with a gene's own span is robust to a one-base shift, so a system that never
lifted anything still reports that the 22 kb gene split in two and that the
189 bp gene truncated to one ORF. What breaks is every column that compares two
runs. On the 22 kb target
the truth is two new calls of 11,233 and 11,187 bases with no collateral at all;
the unlifted run reports three new calls (11,232, 11,187 and 11,013, the third
being a pre-existing call inside the span that moved and so counts as new) with
4,571 calls lost and 4,570 gained.

Measured across the whole ladder, the unlifted run gets **19 of the 21 verdicts
right** and flips two: the 240 bp gene reads as truncated instead of abolished
and the 249 bp one as split instead of truncated, because the displacement turns
a call that was already inside each span into a spurious new one. Its collateral
columns run from 1,539 to 7,745 against gold's 0 on 19 of 21 rows.

So `stub-unlifted` clears R1, fails R2 in every row, and fails R3 by two rows out
of twenty-one. That ratio is the probe's value: a system can be right about what
happened to nearly every gene in the experiment and still not have run the
experiment, and the rung that notices is the one where two coordinate systems
have to meet.

### The ladder was a threshold in disguise

The first draft declared the lengths anyone would declare: 100, 150, 200, 250,
300, 450, 700, 1200, 3000, 10000, 20000. The build passed, the oracle scored 1.0
at depth 4, every assertion held, and the verdict list came out in three blocks.
Abolished through 201 bp. Truncated at 249 and 300. Split from 450 up.

A list in contiguous blocks is its own cut points. Two numbers, 249 and 450,
regenerate all eleven entries, which makes R3's per-length verdict list one
threshold restated once per target. It is the free-component defect again, in a
shape I had not seen before: not a value readable out of the inputs, but a list
whose entries are a function of two of its own. Nothing in the build noticed,
because nothing was looking for it.

What made it visible was a stale sentence. The template docstring claimed, from a
grounding pass run before the ladder was declared, that a 189 bp gene survives as
one truncated ORF while 219 and 261 bp genes are abolished outright. The shipped
ladder contained none of those lengths and reported a clean threshold, so either
the docstring or the ladder was wrong. Re-measuring at a ten-base step from 170
to 280 settled it:

| length | 171 | 180 | 189 | 201 | 210 | 219 | 231 | 240 | 249 | 261 | 270 | 279 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| verdict | trunc | abol | trunc | abol | abol | abol | abol | abol | trunc | abol | abol | trunc |

There is no threshold. Whether a frameshifted gene leaves a callable ORF behind
depends on a start codon and a ribosome-binding site falling in the right place
downstream of the break, not only on how much room the halves have, and the round
numbers step straight over the band where that shows.

The temptation was to fix the docstring. What the ladder's own declared purpose
said, in the rules file, was that it "has to cross wherever that turns over", and
it did not. So the ladder changed instead: the wide anchors stay, and the band
from 170 to 280 is sampled at a uniform ten-base step. Uniform matters. Picking
189 and 219 because they show the inversion would be selecting the finding;
declaring every length in a band whose width is itself the measurement does not
select among them. The provenance of the band, the twelve verdicts above, is
recorded in the rules file next to the lengths.

The shipped campaign: abolished 10, truncated_to_one 5, split_into_two 6, and
`verdict_is_monotone_in_length` is False. R3 is now 21 verdicts that no two
numbers reproduce. And the build refuses gold when that boolean comes out True,
which is the guard the first draft needed and did not have — not a claim that the
data must be non-monotone, but a refusal to ship a rung that grades one cut point
21 times.

The ladder went from 11 targets to 21, the build from 22 caller runs to 32, and
the declared budget from 30 minutes to 60.

### R4 grades a coordinate, not a verdict

The catalogued R4 asks for the new pseudogene call as an indel is introduced
elsewhere. Walked across the 22 kb target at 5, 25, 50, 75 and 95 per cent, the
verdict does not move at all: five offsets, five times `split_into_two`. Grading
it would grade one answer five times.

What moves is where the gene splits, and it tracks the deletion exactly. First
fragment: 1,141 → 6,145 → 11,233 → 17,869 → 21,322 bases, against deletions
1,120 → 5,598 → 11,196 → 16,794 → 21,272 bases into the gene. So R4 grades the
split point, the verdict-as-a-function-of-position is excluded by name, and the
build refuses gold if those first-fragment lengths repeat or come out of order.

### Two refusals the contract caught

Registering this template in `test_build_refusals.py` reported 17 refusals with
no test, which is what that file is for, and writing them found two real defects.

The length-proximity rule read `abs(chosen.length - want) > want`, which is
one-sided. The gap between a 300 bp call and a declared 100 kb target is bounded
by 100 kb, so a declared length far larger than anything on the contig can never
trip it: a ladder asking for 20 kb on a contig whose longest gene is 300 bp would
have audited that gene as its 20 kb target without complaint. It is now a factor
of two in both directions.

And the first mutation I wrote for that refusal reached a different one. Five
copies of a declared length exhaust the single available call, so "no call is
left to stand for" fires before the proximity check ever runs. Same lesson as
last campaign, one layer down: a template's refusals have an order, the order
decides which mutations can reach which refusal, and nothing but these tests
records it.

### The sweep

Six stub levels, three repeats, 72 runs in 15m11s with the memo effective, and
the depth separation is exact: reader 1, respond 2, verdict 3, complete 4,
unlifted 1, noncompute 0. All five gates pass — R1 clear rate 1.0, clear rates
non-increasing at 1.0 / 0.6 / 0.4 / 0.2, cold mean 0.55, no-tool score 0.0000
against a chance floor of 0.2502.

The memo's value shows in the shape of the cost here rather than its size. The
whole campaign is 32 unique caller runs; 72 sweep runs need those same 32 and
nothing more, so the sweep costs about what one build costs plus the sandbox
overhead, instead of about what nine builds cost.

### The same leak, in a shape check 21 cannot see

Writing the ladder's own justification into the rules file, I explained why one
length band is sampled densely by saying what the verdicts in that band are. The
band is twelve of the twenty-one declared targets, so that sentence handed over
twelve of R3's answers in a file the agent reads before doing anything, and
`task.yaml` carried the same list in its deviation note.

Check 21 passed it. It searches `task.yaml` and `reference/` for the string form
of every graded value and looks only at distinctive scalars, by design: an
integer needs four digits and a float four decimals, because a graded count of 2
occurs in any prose by coincidence. Here neither half of the pair is a
distinctive scalar. `170` is a declared target length that has to be in the rules
file, and `truncated_to_one` is published vocabulary that also has to be. The
answer is the two of them in one sentence.

So the mechanism is the one the project already named — a note explaining a
design decision drifts into quoting the measurement that motivated the decision —
and the shape is new: not a value in prose, but a PAIRING in prose.

`tests/unit/test_cold_sandbox_prose.py` checks it across every built campaign:
for each graded value that is a map from a declared label to a term, it fails if
the label and its term occur within 160 characters of each other anywhere in the
cold sandbox. It is bounded to terms the campaign itself declares as a closed
vocabulary, and that bound is what makes it precise. A graded string that is not
a declared term — a molecular formula, a unit, an InChIKey — is published because
the campaign cannot be solved without it, and whether grading it back is a free
component is a different question with a different answer per campaign. Run
unbounded, the scan fired on six campaigns; triaged, one was this leak, three
were free components (open item 10), and two were the aggregate-over-a-published-
table case, which is not a defect at all.

Fixing it meant deleting the verdict list from both files and saying, where it
stood, that the band's verdicts are what R3 grades and are therefore not recorded
in the sandbox. The numbers are in this file instead, which the sandbox never
sees. The audit notes went the same way: they now say what to recompute and from
where, rather than printing what it comes to.

## Ten answers the sandbox was still handing over

Writing the categorical check for the frameshift leak, I generalised it one step
further than the leak needed: for every graded component whose gold value is a
scalar, does that exact value appear, on a word boundary, anywhere the agent can
read before computing? Twenty-five string hits and twenty-one numeric ones came
back across sixteen campaigns. Triaged, five of the numeric hits are genuine
coincidences and the rest split into ten disclosures and twenty free components.

The project has been here before and fixed it: check 21 exists because seven of
sixteen campaigns stated a graded answer in prose. What it does not do is reach
small values. It requires an integer to carry four digits and a float four
decimals, because a graded count of 2 occurs in any prose by coincidence, and
its own docstring says plainly that this leaves small graded values unchecked.
That is the gap, and these are what was in it.

### The construct pair, which is the worst leak so far

Both construct-design campaigns carried this sentence in their own `task.yaml`:

> every unconstrained-violation count (E. coli 20/1/0/29/0/0 against RebH's
> 11/3/0/50/1/0, Streptomyces direct_repeat 23 against 15) [...] So R3's six
> integer components and R4's flipped set remain real work while R3's
> binding-constraint component does not.

The sentence states the answer and then asserts it is real work, in that order.
`construct-ecoli-prna-01`'s R3 is the six E. coli counts — the same six integers
as `20/1/0/29/0/0`, in a different order — and its R4 swapped direct-repeat count
is the `23`. `construct-ecoli-rebh-01`'s are the other two groups in the same
sentence. So each campaign disclosed six of R3's components and one of R4's, and
each also disclosed the other campaign's, in a note whose purpose was release
management: it was arguing that the two instantiations must not be split across
the public and private sets. The argument is right. Printing the counts to make
it was not necessary, and it is the same mechanism check 21 was built for — an
explanation that needs the phenomenon quoting the measurement instead.

The other four, each in a sentence that was true and useful:

| campaign | component | stated as |
|---|---|---|
| `bgcdetect-scoelicolor-01` | `loci_total` | "a 15-row extract of MIBiG's locus table", and again as "the 15 regions with no curated counterpart" |
| `chemspace-mibig-4_0-01` | `multi_class_entries` | "456 active entries carry more than one class", in an exclusion's reason, in an audit note, and a third time in `reference/space_rules.json` |
| `mibig-gene-function-evidence-01` | `backed_Transport`, `unbacked_Transport` | "the starkest category in the backed census (23 backed against 153 unbacked)", in an audit note |
| `residues-prna-01` | `active_site_position` | "Position 79 carries four records at once", in an audit note |
| `selectivity-saureus-topoisomerase-01` | `distinct_assays`, `eligible_molecules` | "every one of the 105 assays", in a deviation note; "a dictionary of 173 floats", in an exclusion's reason |

Every one is now phrased without its number, each saying where the number lives
instead: in gold, which an auditor has and a system does not. The chemspace one
needed fixing in the template's constants as well as the campaign's copy, or the
next rebuild would have restored it.

### Three of the free components, as the pattern for the rest

The twenty remaining hits are a different defect. A graded component whose value
is the campaign's own accession, EC number, snapshot release or declared cutoff
discloses nothing an agent could not read off the task — that IS the task — but
grading it measures a copy. The newer campaigns already exclude exactly these by
name; the older ones grade them, which is drift rather than a decision.

Three are fixed here, chosen because their reason is the most interesting one.
Both mass-balance campaigns publish a per-monomer formula table and graded seven
components straight out of it, and `kinetics` published its unit-per-field table
and graded it back. In both cases the table has to stay: mass-balance is S0, the
agent has no chemistry toolkit, and the table exists precisely so that an RDKit
version bump cannot move a published score; the kinetics units exist so that a
correct number in the wrong convention is not scored as an error. So the table
stays and the component goes, which is the `mibig-diff` precedent exactly.
Mass-balance R2 goes from 6 and 7 components to 3, kinetics R1 from 8 to 5, and
what is left in those rungs is the work the table enables rather than the table:
`all_formulas_parse`, the summed assembly formula, the assembly mass.

### Check 22, and where an exemption has to live

The scan is now the gate's twenty-second mechanical check,
`cold_sandbox_withholds_small_answers`, running the same code as the project-wide
test: `src/npbench_c/readiness/cold_sandbox.py`. It builds the real cold sandbox,
as check 21 does, rather than reading `task.yaml` directly — the question is what
an agent gets, not what a file says. Every campaign was at 22/0/5 when
this landed; check 23 below took them to 23/0/5.

The design question worth recording is where a declared collision can live.

The first version put the exemptions in the test's source, which works but puts
a per-campaign fact in a file nobody reads when auditing a campaign. The obvious
fix is to declare them in `task.yaml`, next to `excluded_from_grading`, where the
project puts every other declaration — and that is exactly wrong here. An
exemption says which declared value a component's answer equals: "uv_gc_window
collides with the declared 50-nt window length" states that `uv_gc_window` is 50.
Declaring it in the sandbox would ship the disclosure the check exists to
prevent.

So they live in `audit/cold_sandbox_exemptions.json`, outside `campaigns/`, where
no sandbox can reach them, and a test asserts both that the file is outside every
sandbox and that none of its reasons has been copied into one. Eight declarations
across five campaigns, each naming its collision:

- four genuine coincidences — a GC-window violation count equal to the declared
  window length, a repeat count equal to the declared k-mer length, two
  corpus-level answers equal to heavy-atom counts in a shipped chemistry table,
  and a miss count equal to an identity-band edge;
- two deliberate groupings, both `environment` on the two Prodigal campaigns,
  where the tool identity is graded as one claim and the declared mode names
  inside it must be published. The composite also carries the caller's version,
  which the sandbox does not state, so it is not obtainable by reading.

Three things check 22 does **not** reach, stated in its docstring rather than
left to be assumed, because check 21's own honesty about its threshold is what
made this gap findable in the first place: a number inside a graded composite
keyed by a label (a census cell's count near its own term — tried on paper and
dropped, since small integers near published vocabulary collide constantly and
there is no measured instance to calibrate against); a number spelled out in
words, which nothing mechanical will ever find; and a value an agent can derive
from the sandbox without doing the work, which is a free component with a
different fix.

A ninth exemption was declared and removed within the minute. `chemspace`'s
`most_shared_informative_scaffold` collides under an unbounded scan — its
InChIKey is in the shipped chemistry table, because that is where the scaffolds
are — but the pair scan is bounded to terms the campaign declares as a closed
vocabulary, and an InChIKey is not one, so the component never fires and the
exemption was protecting nothing. The staleness check found it, which is twice
that check has paid for itself: the first time, removing a leaked sentence from
`bgcdetect` also removed the collision a neighbouring waiver was resting on. The
reason it was dropped is recorded in the file, because the next person to widen
the scan will be tempted to re-add it.

### The copies that were being graded

The twenty hits that were not disclosures were copies. A graded component whose
value is the campaign's own accession, EC number, snapshot release, declared
convention, invocation parameter or published cutoff tells an auditor nothing
about whether the work was done: the sandbox states it, so a system reports it
back. Nothing is leaked that the task does not already say — that IS the task —
and in a conjunctive rung a free component cannot make the rung passable for
nothing. What it does is flatter the per-rung difficulty profile, which is the
one direction that flatters the benchmark.

The newer campaigns already handle these by name. Both Prodigal campaigns
exclude the accession and the mode and grade a single `environment` claim that
includes the tool version, which cannot be read out of the sandbox. The older
campaigns graded them individually. That is drift rather than a decision, and
all twenty are now dropped with the exclusion recorded:

| campaign | dropped from grading | rungs after |
|---|---|---|
| `bgcdetect-scoelicolor-01` | `accession`, `coordinate_convention` | 10 / 2 / 6 / 1 |
| `ecaudit-fmn-dh-1_1_3_15-01` | `release`, `ec_number` | 7 / 4 / 5 / 2 |
| `gcfcutoff-mibig-32-01` | `corpus_size`, `baseline_cutoff`, `include_gbk` | 6 / 2 / 5 / 2 |
| `constraint-gc62-ecoli-rebh-01` | `witness_required` | 4 / 5 / 4 / 14 |
| `constraint-gc55-ecoli-rebh-01` | `witness_required` | 4 / 5 / 3 / 14 |
| `pocket-rebh-01` | `pdb_id`, `accession`, `contact_cutoff` | 8 / 6 / 8 / 6 |
| `residues-prna-01` | `accession`, `swap_structure` | 8 / 8 / 9 / 6 |
| `selectivity-saureus-topoisomerase-01` | `release`, `focal_target`, `reference_target` | 5 / 5 / 5 / 4 |
| `transfer-ec1_14-uniprot-01` | `release`, `ec_prefix`, `tie_break` | 7 / 4 / 4 / 2 |

Three distinctions the work turned on, each of which could have gone wrong
quietly:

- **The grid is not the parameter.** `gcfcutoff` loses `baseline_cutoff`, which
  is one declared number, and keeps `cutoffs`, the declared grid the run swept.
  Which cutoffs were run is a claim about the run; which one is the baseline is a
  copy of a line in the rules file.
- **The witness's requirement is not its achievement.** The constraint pair's
  `witness_required` is the campaign's own `gc_global_min` echoed back; what the
  constructive witness actually achieves is a different component and stays
  graded beside it. The same holds for `residues`: `swap_structure` names the
  structure the R4 swap uses, which `task.yaml` declares, while everything the
  swap does remains graded.
- **Not graded is not withheld.** `selectivity`'s three echoes stay in
  `CENSUS_KEYS`, because R2's warm bundle hands over the whole R1 census and
  removing them from the tuple would change what an agent is given. What to
  grade and what to hand over are separate decisions, and conflating them would
  have made a warm bundle quietly smaller.

Gold is byte-identical in all nine. Dropping a component changes `grading.json`
and nothing else, which is the useful property of generating the spec from gold:
`bgcdetect`'s gold came back at md5 `a1a0541…` for the third time today, through
five antiSMASH invocations, and `gcfcutoff`'s through BiG-SCAPE. Every one is
oracle 1.0 at depth 4 and gate 22/0/5.

### What the check looks like now

`tests/unit/test_cold_sandbox_prose.py` holds both scans — the categorical one
for a term stated against its label, and the scalar one above — and the scalar
one is anchored on the component rather than on a digit count. That produces
real coincidences, so each is listed with the collision named: a GC-window
violation count that equals the declared 50-nt window length, a direct-repeat
count that equals the declared 11-mer length, two corpus-level answers that
equal heavy-atom counts in the shipped chemistry table, a histogram edge.

Two further tests keep the waiver lists honest, and one of them earned its place
inside a minute. `test_every_waiver_still_describes_a_real_collision` fails if an
exempted component no longer collides — and it fired as soon as I removed the
bgcdetect sentence, because `regions_without_a_locus` had been waived for
matching the same `15`. A waiver nobody rechecks is a check nobody is getting.

All nine touched campaigns are back at gate 22/0/5. The three whose grading spec
changed regenerate byte-identical gold, because dropping a component changes
`grading.json` and nothing else. bgcdetect's regeneration was the one I nearly
skipped, on the assumption that antiSMASH over a chromosome is an afternoon: its
own measured budget says five invocations at about seventy-five seconds each, and
the whole gate run came back in minutes with gold byte-identical. The lesson is
smaller than the leak but worth keeping — the cost of a check is in the campaign's
budget note, which is there to be read rather than guessed at.

## Campaign 23: grading the search, not the homology

T-L2-3 is the remote-homology twilight zone, and its catalogued R3 is
"homologous or not per candidate, with the evidence class". The first half has no
oracle. At thirty per cent identity two proteins may share an ancestor or a fold,
the resources disagree, the thresholds are conventions, and a system asked to
decide is being asked to agree with whoever curated the answer. Same shape as
T-L1-4's pseudogene call, and the same resolution: don't ask it.

What replaces it is a measurement of the search. Take InterPro's own published
assertion — two proteins carry the same `family` entry, or only the same
`homologous_superfamily` entry — and grade what a pinned pairwise searcher
reports against it, band by band. The claim is never about homology. It is about
where a pairwise searcher stops seeing a relationship a profile-based resource
asserts, and what it reports instead.

`campaigns/homology-twilight-ipr036188-01`: every reviewed UniProtKB entry
(release 2026_03) cross-referencing the InterPro homologous superfamily
IPR036188, which is 4,105 sequences; labels from InterPro 110.0; MMseqs2
18.8cc5c at sensitivity 4.0. Oracle 1.0 at depth 4, gate 22/0/5, gold
byte-identical on regeneration, stub sweep 5/5 with depth separation
[1, 2, 3, 4]. Four search runs, each rebuilding and reindexing the corpus.

### Why the corpus is one superfamily

Because it makes the question unambiguous. Every pair in the corpus is a pair
InterPro already places in one homologous superfamily, so the campaign never has
to claim that two members are unrelated — a claim nothing published supports. The
question is only ever whether the search can see a relationship asserted at
family level.

That also fixes the vocabulary's third value. A pair sharing no entry is
`unclassified_pair`, not "non-homologous": InterPro records what has been
classified, which is not the same as what is related. The distinction is the same
one `bgcdetect` needed for regions with no curated counterpart, and it is
excluded from grading by name.

### The grounding, which decided the design

An all-vs-all run over the corpus, 1,045,338 hits in 7m28s single-threaded,
against the label:

| | |
|---|---|
| pairs InterPro places in one family | 590,116 |
| of those, the search misses at s=4.0 | **89,353** (15.1%) |
| pairs the search finds that InterPro does not call same-family | 43,319 |
| in the 20–30% identity band: same-family vs superfamily-only | 64,474 vs 29,801 |

The last row is the campaign's name. In that band the two labels occur in
comparable numbers, so identity alone cannot separate them — which is what makes
it a twilight zone rather than a threshold.

The shipped campaign is five declared queries rather than all-vs-all, and the
spread across them is the finding:

| query | declared family | hits | asserted relatives | found | recovery |
|---|---|---|---|---|---|
| A0A084GDU5 | IPR050097 | 18 | 312 | 17 | **5.4%** |
| A1AIE2 | IPR001100 | 188 | 242 | 159 | 65.7% |
| A0LDC6 | IPR004417 | 684 | 887 | 684 | 77.1% |
| A0AI79 | IPR002218 | 771 | 887 | 771 | 86.9% |
| A0AMD1 | IPR004416 | 856 | 887 | 856 | **96.5%** |

One searcher, one sensitivity, one corpus, and recovery runs from five per cent
to ninety-seven. That is the whole campaign in one column, and it is why the rung
grades a per-query census rather than a single rate.

### Two numbers that look like one

`family_members` and `same_family_asserted` differ on every declared query, in
both directions: 238 against 242, 888 against 887, 630 against 887. The first is
the declared family's size; the second counts every corpus member sharing ANY
family with the query, and a protein commonly carries several. Reading one for
the other makes the recovery rate look wrong, and "fixing" either to match would
be the bug, so both the output contract and a test say so.

A related limit, recorded rather than fixed: three of the five queries assert the
same relative count, because all three carry IPR002218 and InterPro families are
not a partition. The rows are still per query — different proteins, different hit
lists — but they are not independent samples, and a reading that treated them as
such would overstate what five queries show.

### The catalogued R4 was degenerate, measurably

The catalogue asked for "the query region carrying the shared motif" to be
shuffled. Measured, an InterPro family match on these proteins covers almost the
whole sequence — 7..482 of 518 on the halogenase this project uses elsewhere — so
that perturbation is a whole-sequence shuffle and every hit dies. One answer, per
query.

So R4 is a ladder of declared window lengths centred on the query's midpoint,
with a zero-width control, and the measured response is better than the
catalogued one in two ways. It is not uniform: `A0AI79` and `A0LDC6` lose nothing
at 60 or 120 residues and then lose 497 and 423 hits at 240, while `A0AMD1`
degrades from 60 upward. And it is ordered. **In all twenty cells the loss rate
falls from the distant band to the close one** — a local shuffle costs the
weakest hits first, at 100% of the under-25% band before it touches a third of
the over-35% band.

The rate is what makes that sayable. The bands differ in size by two orders of
magnitude on some queries (`A0AMD1` has 5 hits under 25% and 794 over 35%), so
the absolute counts say the opposite of the rates: at window 60 it loses 4 of the
weak band and 53 of the strong one. Gold reports the baseline band sizes beside
every rate for that reason, and a band with no baseline hit carries a null rather
than a zero — having nothing to lose is not the same claim as losing nothing.

### The probe is the self-hit

Every query is a member of the corpus it is searched against, so every query
returns itself at 100% identity. That hit is the query, not a finding. Keeping it
inflates each per-query hit count by exactly one, adds one to the top identity
band, and — since a protein shares every one of its own families with itself —
adds one same-family hit and one detected relative.

The symptom is an off-by-one that is internally consistent everywhere: the
recovery rate moves by a fraction of a per cent, the census keeps its shape,
every arithmetic identity in the report still closes, and nothing looks wrong on
its own.

There are two distinguishable places to make that mistake, and the first version
of the probe picked the wrong one. Keeping the self-hit in the HIT COUNTS fails
R1; keeping it in the LABEL JOIN leaves the counts right and fails R2. I built
the R1 version first, and the stub sweep rejected it — not on the campaign, on
the fixture:

    [ FAIL ] r1_clear_rate   R1 pass@1 for generic tooled systems = 0.75;
                             minimum 0.8

That gate exists to catch a campaign whose first rung is unclearable, on the
grounds that such a campaign measures environment setup rather than science. The
stub set is a fixture, not a population, so one stub built to fail R1 drags the
measured rate below the floor and makes the gate report a fixture artefact as a
campaign defect. Every other campaign here has its probe above R1 for the same
reason, which I had not noticed was a convention until this one broke it.

So the fixture does the label-join version and clears R1. The hit-count version
is measured and kept — there is a test pinning the off-by-one, and it is a real
error a real system could make — it just is not what the harness fixture is for.
Worth recording as the general point: a competence probe has to fail *above* the
rung the harness uses to decide whether the campaign is runnable at all.

## The MS2 fixture grounding, and what it says about the tolerance

Open item 4 has stood since the image was built: T-L4-1 and T-L4-2 need matchms,
which is pip-installable, but the lock file records that those rows need "an MS2
fixture set from MassBank with its own grounding pass first: for a spectral match
the invariance question becomes a tolerance question." This is that pass. No
campaign is built from it yet, and the reason is at the end.

### The corpus

MassBank's data repository at release **2026.03.1**, which is the release the
catalogue names, cloned shallow: 141,017 records, 817 MB, CC BY. Each record
carries an accession, an InChIKey, a formula, the precursor m/z, the instrument
type, the ion mode and a peak list, which is everything a dereplication fixture
needs and all of it published metadata rather than anyone's judgement.

One false start worth recording, because it is the kind of definition that looks
right and is not. The catalogue asks for "authentic-standard libraries only", and
MassBank records carry a `CONFIDENCE` comment whose Level 1 value means exactly
that — *Parent Substance with Reference Standard*. Across 141,017 records, **95
carry it**. It is one contributor's convention, not a corpus-wide field, so it
cannot define the fixture. What can is the structure of the corpus itself:

| | |
|---|---|
| MS2 records with an InChIKey, a precursor and ≥5 peaks | 86,735 |
| distinct compounds (InChIKey first block) | 9,182 |
| compounds with ≥2 spectra | 7,249 |
| compounds spanning ≥2 instrument types **and** ≥2 contributors | **2,394** |

That last row is the fixture. Dereplication is only gradeable when the query
compound is in the reference library via a *different* spectrum, and it is only
non-trivial when that spectrum comes from another instrument in another lab. So
the declared split is: one spectrum per compound as the query, taken from the
lexicographically first contributor, and every spectrum from any other
contributor as the library.

A second false start, caught by its own output: my first metadata extraction
matched `ION_MODE` as a substring and so counted `FRAGMENTATION_MODE` lines too,
which reported CID and HCD as ion modes — 39,818 and 39,097 of them. Polarity
cannot be mixed in a reference library, so the fixture would have been built on a
field that was a third fragmentation modes. The corrected pass gives 99,249
positive and 41,768 negative.

### The tolerance question, measured

matchms 0.30.2 with numpy 2.2.6, `CosineGreedy`, 150 queries against 2,034
library spectra, positive mode only:

| m/z tolerance | top-1 correct | rate | mean top score |
|---|---|---|---|
| 0.005 | 120 | 0.800 | 0.7563 |
| 0.01 | 121 | **0.807** | 0.7747 |
| 0.02 | 119 | 0.793 | 0.7803 |
| 0.05 | 115 | 0.767 | 0.8001 |
| 0.1 | 114 | 0.760 | **0.8239** |

Two things in that table, and the second is the one that matters.

**Accuracy and score move in opposite directions.** Loosening the tolerance
raises the mean top score monotonically and lowers the top-1 rate. A loose
tolerance makes everything look more similar, including the wrong things, so a
system reporting a higher similarity is not reporting a better identification.
That is a rung's worth of content on its own.

**The aggregate is stable and the per-query answer is not.** The rate moves by
less than five points across a twentyfold range of tolerance, which would suggest
pinning any value in it. But per query: **20 of 150 verdicts change somewhere**
across the sweep — 3 between 0.01 and 0.005, 4 at 0.02, 10 at 0.05, 17 at 0.1.
Thirteen per cent of the gold answers depend on a parameter whose aggregate
effect is invisible.

So the lock file's phrasing was right and the answer is not the comfortable one.
The tool itself is deterministic: repeated runs and a shuffled library order give
bit-identical score multisets (checked by digest, matchms 0.30.2 / numpy 2.2.6).
What is not invariant is the ANSWER, and the parameter it depends on is the
tolerance. A campaign here must pin the tolerance as part of gold's definition
and grade the sweep as an axis rather than picking one setting and calling the
rest a detail — the same shape as the GCF-cutoff campaign, where the cutoff grid
is the claim and the baseline cutoff is not graded at all.

### Why no campaign yet

Because gold has to come from the pinned image, and matchms is deliberately not
in it. I installed matchms 0.30.2 in this container to run the grounding, which
is fine for a measurement whose version only needs recording — but a campaign
built that way would produce gold no other container could reproduce, and the
gate's determinism check would pass here and nowhere else. The honest state is:
the fixture question is answered, the pin is measured, and the image needs the
three lines recorded in `image/environment.lock.json` before the row can be
built.

## T-L6-3 grounds out: a capped budget whose answer is mostly luck

The last unbuilt row with no provisioning blocker, and the reason it does not
survive is measurable rather than aesthetic.

**Candidate prioritisation under a hard container query cap.** R1 a ranked table,
R2 the ranking matching gold by Spearman within tolerance, R3 the top candidate
and the discriminating criterion, R4 the budget halved. The appeal is obvious: it
is the only catalogued row where the *resource limit* is the subject, and the
note that the budget is a container limit rather than a counted metric is a good
instinct — a system with coarser tool granularity should not be penalised for it.

Two problems, and the first was visible without measuring anything.

**R3 is a judgement or a copy, with nothing in between.** "The discriminating
criterion" is either declared by the campaign, in which case naming it back is a
copy out of the sandbox and check 22 would catch it, or it is left to the system,
in which case there is no oracle and it is the human judgement this benchmark
excludes from grading. Every other row that hit this fork had a third option —
the frameshift campaign replaced "which locus is a pseudogene" with what the
caller reports, the homology campaign replaced "homologous or not" with what the
search finds. Here the third option would be "the criterion the declared rule
used", which is the copy again.

**And R2 has no reachable oracle, which is the measurement.** The row needs the
cap to bind: if the budget is loose enough to measure every candidate, the task
collapses to sorting a table, which is arithmetic several built campaigns already
do. If it binds, gold is the best selection under full information, and the
question is whether a system can get there from what it can see for free.

Measured on the homology corpus, where the all-vs-all run is already paid for:
4,105 candidates, 3,776 of them with at least ten asserted relatives. The
EXPENSIVE quantity is what a search tells you — the fraction of a protein's
asserted relatives it recovers — at one query per candidate. The FREE features are
in the shipped label table.

| free feature | Spearman against the expensive quantity |
|---|---|
| sequence length | +0.239 |
| number of InterPro entries | +0.028 |
| largest declared family | +0.022 |
| asserted relatives | +0.023 |

There is no signal to prioritise on. One feature reaches 0.24 and the rest are
indistinguishable from noise, so a ranking by the expensive quantity cannot be
predicted from anything free — which is exactly the case where a capped budget
decides the answer and no system can reach it except by luck.

The selection version makes the problem concrete and adds the finding that
settles it:

| cap | best possible | best free ranking | random choice |
|---|---|---|---|
| k=10 | 10.0 | 10.0 (100%) | 7.7 (**77%**) |
| k=25 | 25.0 | 22.5 (90%) | 19.0 (**76%**) |
| k=50 | 50.0 | 44.8 (90%) | 38.0 (**76%**) |

**Random choice captures three quarters of the ceiling.** The whole measurable
band between guessing and the best informed strategy is about fourteen points,
and the strategy that achieves it ranks by a feature correlating at +0.022 — it
is not informed, it is lucky in a way that does not generalise. The 100% at k=10
is worse than it looks: thousands of candidates sit at recovery 1.0, so hitting a
ten-candidate ceiling is arithmetic, not skill.

That is a chance floor nowhere near this project's own limit. The readiness gate
requires `chance_floor <= 0.10` and computes it; a rung where guessing scores 76%
of the attainable value cannot pass it under any tolerance-based scorer. So the
row would be rejected by a check that already exists — and it is better to know
that from a measurement now than from a gate failure after a template is written.

**Retired on grounding**, the third row to go that way after T-L2-4 and T-L5-2.
The pattern across all three is worth naming: each died where its R2 or R3 needed
an answer that is either stated by the campaign or known to nobody. That is the
same fork the surviving campaigns resolved by grading what a pinned tool reports
rather than what is true, and the rows that ground out are the ones where no tool
output stands in for the claim.

If the budget idea is worth keeping, the version that would work is not this row:
cap a resource whose optimal use is *computable in advance* from the shipped
inputs, so gold is a declared optimum rather than a full-information one, and
grade the plan rather than the outcome. That is a different campaign and it needs
its own grounding.

## The release packager, and three readings of my own check

Open item 2 has said since Phase 1 that "release packaging must exclude `gold/`,
`oracle/` and `oracle_submission/`", and nothing existed to do it. The exclusion
list was a sentence in a document, which is the shape of defect this project
keeps finding in other people's declarations.

`src/npbench_c/release/package.py` builds the tree and refuses to build one that
carries the answers. It is the other side of `build_sandbox`: that assembles what
one system sees for one campaign and audits it, this assembles what leaves the
building, and a release that ships gold hands over every answer at once and
cannot take it back.

Three exclusions are less obvious than `gold/`, and two of them now have tests
whose job is to prove the exclusion is load-bearing rather than tidy.

**`docs/` is the sharpest one, and it is this file.** The project's own rule for a
measurement that motivated a design decision is that it belongs in `docs/`,
"which the sandbox never sees" — which is true of a sandbox and false of a
release that shipped `docs/`. Everything I have written here all session, the
verdict censuses and recovery rates and call counts, is the answer key. A test
asserts that `docs/*.md` really does state five or more graded four-digit values,
so if anyone later argues the exclusion is paranoid, the measurement answers.

**`src/` is the second.** Each template's `analyse.py` and `build.py` ARE the
oracle: they compute gold from the shipped inputs. A participant given the source
can regenerate every answer, tool runs included. No source ships at all, and the
test asserts the release contains no `.py` file.

**`grading.json` is the third**, and it is not an answer — it is a map of where
the marks are: which fields are graded, under which scorer, at what tolerance.
Organiser material. The manifest is checked for the same reason not to carry a
digest of gold: a digest is not gold, but it lets a system confirm a guess, which
is the same thing one bit at a time.

### The attribution check, and what it took to get right

A release that cannot attribute its data cannot ship — several inputs here are
CC BY, where attribution is a condition rather than a courtesy. So the packager
assembles `ATTRIBUTION.md` from each campaign's own statements and refuses if a
campaign states none.

I got that check wrong twice before it was right, in both directions.

**Too loose.** The first version accepted any key whose name contained "source",
and `construct-ecoli-prna-01` passed on the strength of `"source": "PubMed"`
inside an evidence record — a field about where a claim came from, not about what
may be redistributed. A rule that accepts that accepts almost any JSON.

**Then too shallow.** Tightening it to licence-naming keys only, read at the top
level, reported `pocket-rebh-01` and `residues-prna-01` as unattributed. Both are
properly attributed; they state their licences inside `annotation_source` and
`structure_sources`, two levels down. A licence nested in a provenance record is
still a licence, so the check walks the document.

With both fixed, the real answer is one campaign, not twelve and not three:
**`construct-ecoli-rebh-01` shipped UniProtKB data with no provenance file at
all**, while the second instantiation of the same template shipped one. The pair
diverged and nothing compared them. Fixed by adding the record, with the values
computed from the shipped FASTA rather than copied from its sibling — 530
residues, `sha256_sequence` 2853460230442a06, Q8KHZ8, CC BY 4.0 — and the
campaign re-gated at 22/0/5 with gold byte-identical, since adding a provenance
file changes the sandbox but no computation.

Both of my wrong readings are now tests rather than anecdotes: one asserting that
a `source` field is not a licence, one asserting that a nested licence counts.

### What the release is

23 campaigns, 161 files, 36 MB: per campaign the `inputs/`, the `reference/`
tables, `task.yaml` and an empty `submission/`, byte-identical to the build repo
— a release is a copy, not a rendering, because anything that transformed a
pinned table in transit would make the shipped campaign a different campaign from
the one the gate passed. Plus a manifest with a digest per file and the grader
version, an attribution file, and a README that says what is deliberately absent
so nobody hunts for it.

## Where the codon tables come from

Item 1 of the open list, closed. The construct-design campaigns grade a CAI floor
and six violation counts over a 64-codon relative-adaptiveness table per host,
and until now those tables were *representative values for highly-expressed
genes*: no release, no gene set, nothing an auditor could recompute. The JSON
said so, which is better than pretending otherwise, and the honest options were
to cite a source or to derive one.

Citing was not available. The shipped values are rounded to two decimals and
match no release I can verify; asserting a citation I cannot check would be worse
than the gap. So the tables are now derived, by a rule declared in
`constants/codon_table_sources.json` and run by
`npbench_c.templates.construct_design.prepare_tables`:

> Every CDS feature of a pinned INSDC record whose `/product` names a ribosomal
> protein, that is not `/pseudo`, whose location is one contiguous range, whose
> `/codon_start` is 1 and `/transl_table` is 11, whose bases are a whole number
> of ACGT codons ending in exactly one stop, and whose translation under the
> standard code equals the record's own `/translation`. Relative adaptiveness is
> count ÷ the most-used codon of the synonymous family, over all codons of those
> genes including each terminal stop, with an unobserved codon given 0.5
> observations first (Sharp & Li 1987).

Ribosomal proteins rather than a genome-wide table, because CAI wants a reference
set under translational selection and a genome-wide count measures mutational
bias as much as selection — which is the thing the weights are supposed to
separate. They are also the one set every host annotates under a name a regex can
match, and the regex had to stretch: *S. coelicolor* writes "putative 50S
ribosomal protein L28", BL21 writes "30S ribosomal subunit protein S12", and
both spellings had to match while "ribosomal protein S12 methylthiotransferase"
and "50S ribosomal protein L16 arginine hydroxylase" had to miss. What it
selects: **54 genes / 7,159 codons for `CP001509.3`** (E. coli BL21(DE3) — the
host the campaigns actually name, not K-12) and **60 genes / 7,912 codons for
`AL645882.2`** (*S. coelicolor* A3(2), the record already shipped as the BGC
detection campaign's input, so that host's table is recomputable from the
repository alone; the BL21 record is committed under `audit/records/` for the
same reason).

Two cross-checks fell out. The 54 BL21 genes are exactly the 54 CDS features
carrying an `rps*`/`rpl*`/`rpm*` gene name, which the rule never looks at. And
the derived E. coli table looks like a highly-expressed table rather than a
genome-average one — GCC drops from 0.72 to 0.188, CTG stays at 1.0 with CTC at
0.030 — which is what the reference set is for.

**What is pinned is the derived set, not the file.** The manifest records each
record's SHA-256 but does not enforce it: the same accession from another INSDC
mirror is a different byte stream and the same biology, and a refusal on the
transport encoding would fire on a correct re-run and teach the next person to
pass a flag. What it does enforce is the gene count and a digest over the
reference set itself (each gene's location, locus tag and length), so a changed
annotation release or a pattern that started matching something new refuses to
write a table. `--adopt-digests` re-pins deliberately.

### What it moved, measured before it was adopted

The reason to measure first: re-deriving a table that gold depends on is a gold
rebuild, and the size of the rebuild decides whether this is a provenance note or
a four-campaign job. Built against scratch copies, then again in place:

| campaign | gold keys moved |
|---|---|
| `construct-ecoli-rebh-01` | 20 of 41 |
| `construct-ecoli-prna-01` | 17 of 41 |
| `constraint-gc62-ecoli-rebh-01` | **0 of 62** |
| `constraint-gc55-ecoli-rebh-01` | **0 of 71** |

The two constraint campaigns are byte-identical because they read the table for
which codons encode which residue and what GC each can reach, and never read the
weights — now stated in both campaigns' audit notes, with the measurement behind
it. In the two design campaigns what moved was every CAI quantity (floors 0.902 →
0.888 and 0.929 → 0.880 for RebH, 0.908 → 0.883 and 0.914 → 0.848 for PrnA),
every cost spread, and **five of the twenty-four unconstrained-violation
counts** — three in RebH (E. coli `direct_repeat` 11 → 12 and `gc_window` 50 →
30, *Streptomyces* `direct_repeat` 15 → 14) and two in PrnA. What did *not* move:
the binding constraint in any host, the non-binding sets, and the flipped-
constraint sets. R3's and R4's categorical answers are stable under the table;
their integers are not, which is the right way round — the integers are the part
that depends on the weights.

Two graded integers changing also made two cold-sandbox exemptions stale, and the
gate failed the campaign until they were dropped. Both were genuine coincidences
(`uv_gc_window` equalled the declared 50-nt window length; `uv_direct_repeat`
equalled the 11-mer repeat length) and both stopped colliding the moment the
optimum moved. That is the staleness half of check 22 doing exactly what it was
written for: an exemption outlives the number it was written for, and the next
reader takes it for a standing allowance.

### One thing the release packager should not learn from this

The derived tables carry a `license` per record, and the release's attribution is
assembled by walking each campaign's provenance for licence-named keys — so the
obvious next move is to point that walk at `reference/` too. Measured before
doing it: over all 23 campaigns it would add the two INSDC records per
construct-design campaign, which is right, and a block from
`pocket-rebh-01`/`residues-prna-01`'s `eco_policy.json` whose only licence-shaped
key is `evidence_citation_unit` — a counting rule about PubMed citations, not an
attribution. That is the same false reading the attribution check has now
produced three times, so the walk stays on `task.yaml` and `inputs/provenance.json`
where it is calibrated. The records themselves are INSDC public domain, and the
provenance ships inside the reference file that travels with the weights, which
is where an auditor reading the table will look.

### Three more answers the sandbox was stating, and the shape they were in

Rewriting the two design campaigns' `task.yaml` meant re-reading the notes around
the numbers, and the split-constraint note had been naming the answers in words
since before the free-component audit:

> "3 of 6 graded sub-answers are identical: the binding constraint for both hosts
> (gc_window for E. coli, direct_repeat for Streptomyces)"

That is R3's `binding_constraint` and R4's `swap_binding_constraint`, both of
them, in both instantiations. The note's own tail says the earlier version had
printed R3's integers and that they were removed — and then it kept naming the
constraints, which is the same disclosure in words. Check 22 could not see it:
its pair scan looks for a label and a term *inside the graded value*, and here
the value is a bare string whose label is in the gold **path**
(`hosts.ecoli.binding_constraint`). Nor would a literal search have found the
pair, because the path writes the host as `ecoli` and the prose writes it as
"E. coli".

So the scan grew a third shape: for a graded value that is one term of published
vocabulary, take the innermost segment of its gold path as the label, and look
for the two within a window of each other with **both sides normalised to
lower-case alphanumerics** — which collapses "E. coli", "E.coli" and `ecoli` onto
one string, and `gc_window` onto "GC window" too. Two refinements the measurement
forced, in this order:

- **The label is the innermost path segment, not every segment.** With all
  segments, `hosts.ecoli.binding_constraint` also pairs on `hosts`, which matches
  the host block of `codon_tables.json` and fires on every construct campaign for
  nothing.
- **Shape three reads sentences, not structure.** A nested document states a
  label against a field name by construction: the block keyed `ecoli_bl21`
  declares a field called `gc_window`, which is not a disclosure. So shape three
  runs over the sandbox's prose with its structure dropped — JSON strings long
  enough to be sentences, YAML and Markdown lines with their leading key removed
  and their value kept. Shapes one and two still read everything.

Across all 23 campaigns the widened scan finds three answers, and it found the
third one in a campaign nobody was looking at: `massbalance-nrps-malleobactin-01`
stated "All four verdicts are additional_chemistry_required", which is four
graded components — *and* named the sevadicin instantiation's verdict, which is a
different campaign's graded answer. For the record, where no sandbox reads it:
malleobactin's four congeners are all `additional_chemistry_required` and
sevadicin is `balanced`, which is why the field is graded at all even though it
discriminates nothing within either campaign.

That cross-campaign case is a gap this check does not close. The scan compares
each sandbox against its own gold; a sandbox stating a *neighbouring* campaign's
answer passes 23 times out of 23. Both instantiations of a template already
declare a `split_constraint` for the related problem of shared answers, and the
arithmetic for a cross-campaign pass is 23 × 23 runs of a scan that takes two
seconds for all 23 — so it is affordable and it is now an open item rather than a
finished one.

## Check 23: the leak that passes a per-campaign scan every time

Checks 21 and 22 compare a cold sandbox against **its own** gold. A sentence in
one campaign's notes that names a *neighbouring* campaign's graded answer passes
both, 23 times out of 23 — and two of them had shipped. The malleobactin
campaign's audit notes named the sevadicin instantiation's verdict while
explaining that the field varies across instantiations (fixed in the previous
section), and the codon-table provenance block added in that same commit stated
`record_gc_fraction: 0.721187` for AL645882.2 in four campaigns' reference
tables — which is exactly what `genecall-scoelicolor-01` grades over that
chromosome. Four campaigns were stating a fifth's R2 answer, three days old, my
own, and found by the check written to look for it.

### Which shapes to run is a measurement, not a preference

All 506 ordered pairs of the 23 built campaigns, each shape counted separately,
before any of it was wired into the gate:

| shape | hits over 506 pairs |
|---|---|
| distinctive scalar (check 21's class: 4+ digits, 4+ decimals, whole composite) | 5 |
| small scalar (check 22's shape one: a graded 15 or 0.62 in a sentence) | 77 |
| label → term inside the graded value (check 22's shape two) | 1 |
| label from the gold path (check 22's shape three) | 0 |

So check 23 runs the distinctive scalars and both pair shapes, and **not** the
small scalars. What makes a two-digit integer checkable within one campaign is
anchoring on the component — 22 graded values, each searched once — and across 22
neighbours that anchor is gone: 77 hits, every one a graded count of 12 or 30
appearing somewhere in somebody else's prose. Running it anyway would have
produced a check nobody could keep green, which is how waiver lists start.

Of the six non-small hits, four were the provenance block, one was the
malleobactin case, and the sixth was a false positive worth keeping: the
chemical-space campaign ships a compound named **`S2502`** and the MIBiG diff
campaign grades a source-entry count of **2502**. Check 21's boundary excludes a
neighbouring *digit* — its own comment says so, about `BGC0002359` — but not a
neighbouring letter, so the compound name read as the count. Both checks now
share one `states_scalar`, with the boundary widened to any word character, and
one `distinctive_forms`. The duplication was the real defect: two
implementations of "is this number stated" drift, and only one of them had the
accession-table comment.

### Scope is the set shipped together, so it is passed in rather than discovered

A disclosure only matters where both campaigns can reach the same system. The
gate passes every built campaign beside the one under test; the release packager
passes **its own** subset, so dropping a campaign cannot refuse a release over a
sentence about it, and adding one forces a re-check. `audit_release` now runs
checks 22 and 23 against the assembled tree, which is the point at which
"shipped together" stops being a figure of speech.

### What made it affordable

A sandbox copies the campaign's inputs — tens of megabytes apiece — and the scan
needs 22 neighbours' declared vocabularies per run. It does not need their
sandboxes: a sandbox's JSON outside `inputs/` is exactly `reference/`, so
`declared_terms(campaign)` reads the vocabularies straight from the repository.
That shortcut is only valid while the two agree, so a test asserts
`declared_terms(campaign) == vocabularies(sandbox)` for all 23. Measured cost of
check 23 inside a gate run: **0.4s** for a construct-design campaign and **6.0s**
for the chemical-space one, whose own prose is a 3,000-row chemistry table.
Against a gate run of minutes that is free, and the whole-suite version in
`tests/unit/test_cold_sandbox_prose.py` — every campaign read against every other,
sandbox builds included — runs in about fifteen seconds, because there the
sandboxes are built once and shared.

### One declared collision, and it is the composite gap again

`frameshift-scoelicolor-01` states `genecall-scoelicolor-01`'s
`environment[baseline_mode]`. Both run the same ab initio caller on the same
chromosome, both must publish the mode names their own finding is about, and the
frameshift campaign's note explains that its baseline is the metagenome mode
rather than the single-genome one — which names the other campaign's mode cell.
But `environment` is graded as ONE claim: caller, version, mode and flags
together, and the version in it is the caller's self-reported string, obtainable
only by running the caller. One cell of a composite is not the composite, which
is the gap check 22's docstring already names. Declared under
`cross_campaign_exemptions` in `audit/cold_sandbox_exemptions.json`, keyed by the
reader and naming the campaign and component it covers, with the same staleness
rule as the within-campaign waivers: an exemption that no longer collides fails
the gate.

## Item 9: the field was decorative, and the census is why it is gone

`structural_elements` was a list in every campaign's `task.yaml` —
`verification`, `counterfactual`, `control`, `planted`, `abstention`,
`abstention_free`, `recovery` — that no code read and no check validated. The
open item said it needed either a definition and a gate check or removal. I
tried to write the definition first, and the census is what decided it:

| term | declared by | what an artefact could decide |
|---|---|---|
| `verification` | 23 of 23 | nothing to decide — every campaign recomputes from the submission; that is what measure-then-grade *is* |
| `counterfactual` | 23 of 23 | nothing to decide — all 23 ladders have exactly four rungs, and R4 is the counterfactual |
| `control` | 14 | no common signature: a negative subset (`ripp`), two chains that must agree (`pocket-rebh`), a clean comparison section (`ecaudit`), declared control entries in gold (`gcfcutoff`) |
| `planted` | 11 | **no campaign ships an altered input** |
| `abstention` | 1 | a graded gold value that means "there is no answer" — true for `constraint-conflict` |
| `abstention_free` | 2 | the absence of the above, which the absence of the term already says |
| `recovery` | 1 | what that campaign measures, not a property of its structure |

Two of the seven terms were universal, so they carried no information. Three
appeared once or twice and described what the campaign measures rather than how
it is built. And the two load-bearing ones did not survive contact with the
artefacts.

### `planted` was the one I expected to keep

The previous revision of this item said `frameshift-scoelicolor-01` "deletes a
base from the chromosome itself, and that planting is the whole reason its R3 is
gradeable", and offered it as the worked example that would make the definition
writable. It is not what the campaign ships. Its own provenance says the contig
is "the same extraction as the gene-calling campaign on this accession, so the
two are comparable: identical bases, identical digest" — and the mutated copies
are deliberately **not** shipped, because shipping them would hand over both the
target positions and the coordinate lift. The planting is a rule the campaign
applies at run time, not an artefact in `inputs/`. Same for
`homology-twilight-ipr036188-01`, the other campaign that explained its claim:
the perturbation is a declared shuffle with a declared seed, and its corpus
provenance says "Nothing is sampled or filtered: the whole reviewed set is the
corpus".

So in 23 campaigns, zero planted artefacts and eleven `planted` declarations.
The pattern in where those eleven sit is the tell: **every campaign that wrote
out its reasoning above the field concluded "nothing is planted"** — eight of
them, in those words — and the nine that claimed `planted` without explanation
are all from before that convention existed. A label nobody has to justify
drifts to whatever sounded good.

### And one of them was answer-adjacent

This is what settled it rather than merely argued it. The constraint control
pair grades one verdict each: `constraint-gc62-ecoli-rebh-01` answers
`infeasible_gc_unreachable` and `constraint-gc55-ecoli-rebh-01` answers
`feasible`. The conflict half declared `abstention` in its sandbox-facing
`task.yaml`; the feasible half declared `control`. The label told a reader which
half of a two-verdict pair they were holding — in a field added for tidiness,
which none of checks 21, 22 or 23 can see, because "abstention" is not the
verdict string and appears in no declared vocabulary.

### What replaced it

Nothing, in `task.yaml`. Check 24, `design_claims_are_checkable`, keeps the
field out under any of five names — the defect was unverifiability, not the
spelling — and validates the design claims that *do* resolve:

- a campaign name mentioned anywhere in `task.yaml` must be a campaign that
  exists, so a renamed or retired campaign cannot leave a sentence pointing at
  nothing;
- `control_pair_with`, where present, must name a campaign that names this one
  back. The constraint pair is symmetric today; a one-way pairing is a control
  pair that has quietly become two unrelated campaigns.

What the retired terms gestured at all still exists, in places something reads:
the four-rung ladder in `grading.json`, `gcfcutoff`'s five control entries in
gold where a component grades them, the abstention verdict in the graded
vocabulary, the deletion rule in `frameshift`'s pinned rules file. The two tests
that asserted `"planted" in structural_elements` now assert the artefacts
instead — that the frameshift contig's digest matches the gene-calling
campaign's over the same accession, and that the homology corpus is the whole
reviewed set — which is the difference between a label and a claim.

## Item 8: the sweep was slow for three reasons, none of them the one on record

The open item said the tool memo covers the tool and not the work around it, and
named two candidate fixes: cache the parsed inputs, or let a stub level reuse one
sandbox across repeats. Both were guesses. Profiling the four phases of a single
run — build the sandbox, run the stub, run the measurer, grade — put the cost
somewhere else entirely:

| phase | ripp | architecture | bgcdetect |
|---|---|---|---|
| `build_sandbox` | 0.17s | 0.01s | 0.01s |
| the stub subprocess | **10.06s** | **5.27s** | **4.19s** |
| the measurer | 0.07s | 0.07s | 0.07s |
| grading | 0.20s | 0.01s | 0.00s |

So the sandbox copy, which this item blamed in writing, costs a hundredth of a
second — the inputs are in the page cache after the first run and `shutil.copytree`
of 7 MB is nothing. Nor was parsing the whole story. Inside those stub
subprocesses were three separate defects:

**A `.tar.gz` read out of order.** The RiPP corpus opened its reference archive
once per record, through a per-record reader closure whose docstring explained
that holding twelve megabytes of GenBank in memory to read fifty of them was a
cost worth avoiding. The architecture corpus extracted members in *name* order,
which is the same defect wearing a sort. A gzip stream cannot seek backwards, so
each out-of-order read re-inflates from the beginning: 5.3 of ripp's 7.6 second
corpus load, 133,873 calls to `zlib.decompress`, against 0.2 seconds for one
pass. The fix is `npbench_c.tools.archive.read_members` — one open, one pass,
members in archive order, and the caller sorts the result. Four call sites moved
onto it, including `mibig_diff`, which had the same sort-before-reading pattern
and had not been measured.

**A 1.57 GB file hashed once per run.** Three stubs reported `pfam_sha256` by
hashing `Pfam-A.hmm` — 4.5 seconds a run, six minutes a sweep, to recompute a
number the campaign's own `reference/` already pins under
`external_resources.pfam_a_hmm.sha256`. The stubs now report the declared digest.
That is sound because the pin is not taken on trust where it matters: the oracle
calls `Resource.verify()`, which hashes the file and refuses to build gold on a
mismatch, and the readiness gate runs the oracle. A real system is handed no
shortcut; a stub is a fixture.

**A tool asked its own version twice per run.** This was the whole of bgcdetect's
4.19 seconds, and it is not parsing at all: `antismash --version` costs **2.12s**
because it imports the entire package to answer, and the stub asked twice.
Measured across the nine tools this benchmark drives, three behave that way —
antiSMASH at 2.12s, BiG-SCAPE at 2.38s and, measured when the dereplication
campaign was built, matchms at **9.6 to 11.7s**, which is the worst of them
because printing its version pulls numpy, scipy, numba and llvmlite — while
hmmscan, hmmsearch, prodigal, mmseqs and diamond all answer in 10ms or less. So
`cache.memoize_version` memoises a version string per process always (a tool
cannot change version inside one process, so this needs no opt-in) and on disk
where the tool memo is already allowed, and only those two call sites use it. A
memo around a 10ms call is a mechanism nobody needs.

### What it bought, measured from a cold memo on both sides

| campaign | recorded before | now | what remains |
|---|---|---|---|
| RiPP precursors (Campaign 18) | 7m51s | **45s** | HMMER, once per distinct scan |
| MIBiG 32 cutoff sweep (Campaign 17) | 6m27s | **83s** | one BiG-SCAPE run |
| S. coelicolor regions (Campaign 16) | 15m35s | **8m48s** | antiSMASH, ~80s per distinct run |
| NRPS architectures | not recorded | 3m44s | HMMER |
| PKS architectures | not recorded | 2m34s | HMMER |
| MIBiG 3.1→4.0 diff | not recorded | 1m03s | no tool |

Per-run overhead, which is what the fixes actually changed: ripp 10.06s → 0.61s,
architecture 5.27s → 0.52s, bgcdetect 4.19s → 0.07s. Every campaign's
`internal_sweep.json` is byte-identical to the one it had before, and every gate
still passes — the sweep measures the same thing, faster.

What is left is the tool itself. bgcdetect's remaining 8m48s is about six
distinct antiSMASH invocations at 80 seconds each, and the frameshift campaign's
sweep is roughly one build because its ladder calls Prodigal once per target.
Neither is overhead: it is the one computation the memo exists to pay for once.
Shaving it would mean memoising across campaigns or across sweeps, which trades a
property — that a sweep is reproducible from a clean checkout — for minutes, and
that trade is not worth making.

And the second candidate fix this item proposed, reusing one sandbox across
repeats, is now dead rather than deferred: it would have saved the 0.01 seconds
`build_sandbox` costs, at the price of three repeats no longer being independent.

## Item 5: the antiSMASH version coupling, settled by counting

MIBiG publishes its reference set as `mibig_antismash_4.0_gbk_as8b1.tar.bz2`,
processed with antiSMASH **8.0 beta 1**. The image pins **8.0.4**. antiSMASH's own
rules forbid comparing regions produced by two of its versions, so this was
recorded in three places as a coupling to settle before shipping, with two ways
out: reprocess the reference set with the pinned antiSMASH, or state that both
sides come from the published set.

Counting first, because the choice depends on who actually touches what. Over all
23 campaigns, from their own declarations and provenance:

| | runs antiSMASH | ships antiSMASH-processed records |
|---|---|---|
| `bgcdetect-scoelicolor-01` | **8.0.4** | no — a raw INSDC chromosome and MIBiG's curated loci |
| `ripp-precursor-mibig-4_0-01` | no — HMMER | yes, the published set |
| `arch-nrps-mibig-4_0-01` | no — HMMER | yes |
| `arch-pks-mibig-4_0-01` | no — HMMER | yes |
| `gcfcutoff-mibig-32-01` | no — BiG-SCAPE | yes |
| the other 18 | no | no |

**The two groups do not overlap.** The campaign that runs antiSMASH reads a
chromosome nobody processed, and the four that read processed records drive a
different tool entirely — they use the records as data: CDS translations to scan
with HMMER, GenBank files to cluster with BiG-SCAPE. Not one comparison in the
suite puts 8.0.4 output next to 8.0-beta-1 output.

So reprocessing buys no campaign anything, and it would cost something real: the
records are a published, citable artefact with a URL and a digest, and
reprocessing them locally makes them ours. The resolution is the second option,
and it was already written in prose in all four campaigns' provenance — ripp's
says the records "are NOT re-run through the image's antiSMASH: this campaign
never compares output from two antiSMASH versions". What was missing is the thing
this project keeps finding missing: anything that notices when a sentence stops
being true.

### Check 25, and the field it reads

`antismash_versions_are_not_mixed` is the gate's twenty-fifth mechanical check:

- a campaign that ships antiSMASH-processed records must **record which
  antiSMASH processed them**, in `inputs/provenance.json` under
  `antismash_version_that_processed_them` — a field, where a check can read it,
  rather than a version left to be inferred from an archive's filename.
  `gcfcutoff` already had the field; ripp and the two architecture campaigns said
  it in prose only, and now carry it;
- a campaign that *also* runs antiSMASH must run the version that processed them,
  and a failure names both ways out, because either is legitimate: reprocess with
  the pinned tool, or take both sides from the published set and stop running it;
- a declared `antismash@` pin must match the registry's, so a campaign cannot
  claim a tool the image does not carry.

Two details are deliberate. The scan reads only `source` and `url`, never a
free-text note: the one campaign that runs antiSMASH explains in its own note
that it ships the complete chromosome because a slice truncates CDS features
"which antiSMASH refuses", and a scan over prose would read that as shipping
processed records and fail the campaign doing it right. And the reprocessing
route has its own test, so the check cannot forbid the fix it recommends.

## Item 4: the image carries matchms, and two things were wrong on the way in

The pin was decided two days earlier and the grounding was done; what was left
was provisioning. That turned out to be the easy half.

**The install.** `micromamba install` into `/opt/npbench-tools` with the
Dockerfile's full spec list, dry-run first and read before it ran: **46 packages
added, nothing upgraded, nothing downgraded, nothing removed**. numpy, scipy and
pandas were already in the closure — antiSMASH brought them — so what arrived was
numba, llvmlite, rdkit, pyteomics, sparsestack, pickydict and their dependencies.
`registry --verify --prefix` then reports **9/9 tools at their pinned version**
where it reported 8/8, the banned `FastTreeMP` is still absent, and the full
invariance suite re-run after the install is **9 pass / 0 fail / 0 skip** with
every previous verdict unchanged.

Three things that were not just bookkeeping:

**The lock file is generated now, because the hand-assembled one went stale
immediately.** Its own prose said so: "THE NEXT BUILD MUST RE-SOLVE AND
REGENERATE THIS FILE". `python -m npbench_c.tools.lockfile` reads the prefix's
`conda-meta` — one JSON per package carrying exactly the fields the lock wants —
and rewrites the mechanical half while carrying the editorial half over, because
a per-pin rationale is writing and regenerating it would delete it. 185 packages
became **226**, and the file gained a `build_history`: five solves between
2026-10-02 and 2026-10-21, which is the honest shape of this prefix rather than
the single `micromamba create` the Dockerfile presents. A closure assembled over
five solves is not guaranteed to equal what one solve produces, and that is
exactly why the artifact hashes are there.

**A library pin cannot be probed through PATH.** `MATCHMS.binary` is
`python3.11`, and `shutil.which("python3.11")` finds the *system* interpreter —
which has no matchms, and whose being a different interpreter is the entire point
of the pin. So `--verify` reported the pin missing on an image that has it. The
probe now resolves the tool prefix explicitly. A binary like `hmmscan` never had
this problem because nothing else on PATH answers to that name.

**TOOLS membership means the suite measured it**, and matchms arrived with
`invocations=()`, which the suite's own contract forbids — "a tool with no
invocation is a tool not measured". So it got one: the prefix interpreter running
a committed driver (`fixtures/matchms_cosine.py`) over a committed fixture
(`fixtures/ms2_massbank.json.gz`), writing each query's top match at full float
precision, because a rounded score would hide the bit-level difference the suite
exists to catch. Verdict: **raw-identical** across repeats. The tool has no
thread flag, so the suite runs it at one thread — but the threading that could
bite a numpy-backed tool is BLAS, so that was measured separately:
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and
`NUMEXPR_NUM_THREADS` at 1 and at 8 give the byte-identical artifact.

### The release the grounding named does not exist

The grounding recorded its corpus as **MassBank 2026.03.1**. Cloning that tag
fails: `fatal: Remote branch 2026.03.1 not found`. The repository's tags run
2023.06 … 2025.10, **2026.03**, and that is the latest release. What does carry
the string is the *dev branch's* `VERSION` file:

```
main: version = 2026.03
dev:  version = 2026.03.1-SNAPSHOT
```

So the fixture was cut from a moving branch and its development version string
was written down as a release. Nobody — including me, four days later — can
reproduce that corpus. Re-run against the released tag, every census number
shifts by about half a per cent, exactly as a branch a few hundred records ahead
would:

| | tag 2026.03 | the dev snapshot |
|---|---|---|
| records | 139,240 | 141,017 |
| MS2 with an InChIKey, a precursor and ≥5 peaks | 85,146 | 86,735 |
| distinct compounds (InChIKey first block) | 9,078 | 9,182 |
| compounds with ≥2 spectra | 7,142 | 7,249 |
| spanning ≥2 instruments **and** ≥2 contributors | **2,384** | 2,394 |
| positive / negative | 97,803 / 41,437 | 99,249 / 41,768 |

The structure survives, which is the useful part: the selection rule reproduces
on a release anyone can fetch, and the fixture is now cut from the tag with the
rule, the census and the correction recorded in `fixtures/provenance.json`.

### The tolerance table, re-measured on the stack that will produce gold

The grounding ran on numpy 2.2.6 and scipy 1.16.3. This image has **numpy 2.2.5
and scipy 1.15.2** — so the table could not be inherited, because the whole
reason to pin a numerical stack is that the scores are floats. Re-measured here,
150 queries against 2,321 library spectra from the tag:

| m/z tolerance | top-1 correct | rate | mean top score | the grounding's rate |
|---|---|---|---|---|
| 0.005 | 117 | 0.780 | 0.7528 | 0.800 |
| 0.01 | 120 | **0.800** | 0.7727 | 0.807 |
| 0.02 | 118 | 0.787 | 0.7798 | 0.793 |
| 0.05 | 115 | 0.767 | 0.7980 | 0.767 |
| 0.1 | 118 | 0.787 | **0.8268** | 0.760 |

Both findings a campaign depends on survive. **Mean top score rises
monotonically with tolerance while accuracy does not** — a loose tolerance makes
everything look more similar, including the wrong things. And **the aggregate
hides the per-query answer**: the rate moves by 3.3 points across a twentyfold
range of tolerance while **19 of 150 verdicts change somewhere** in it (3 at
0.005, 4 at 0.02, 9 at 0.05, 16 at 0.1 against the 0.01 baseline). Thirteen per
cent of the gold answers depend on a parameter whose aggregate effect is nearly
invisible, so a campaign here pins the tolerance as part of gold and grades the
sweep as an axis.

One difference from the grounding worth stating rather than smoothing: its table
showed accuracy falling monotonically to 0.760 at tolerance 0.1, and this one
dips at 0.05 and recovers to 0.787. The monotone *fall* was an artefact of that
corpus; the monotone *rise* of the mean score, and the per-query instability, are
not.

## T-L4-1 grounded: three of its four rungs needed redesigning

*Every number in this section is the grounding's, measured on the split the
grounding declared. Building the campaign found that split wrong and the numbers
moved; the built campaign's are in the next section, and the two should not be
read as the same measurement.*

The catalogued ladder was: R1 an annotation table keyed by InChIKey · R2 matches
gold at the declared InChIKey block · R3 the confidence level per annotation and
the discriminating fragment · R4 an adduct shift is applied, predict which
annotations survive. Grounded against MassBank 2026.03 before any oracle was
written, two of those rungs do not survive contact with the data and one needed
a different shape.

### The pipeline the grounding declares

Candidates are library spectra whose precursor m/z is within a declared window
of the query's; they are scored with matchms `CosineGreedy` at a declared peak
tolerance; the top match is the highest score, then the most matched peaks, then
the lowest accession. The precursor filter is not an optimisation — it is what
dereplication practice does, and it is the thing an adduct misassignment breaks.
Over the eligible corpus it leaves a median of **18 candidates** out of 30,749,
with a maximum of 159.

The corpus: MassBank 2026.03, **[M+H]+ positive-mode MS2 only**, which holds the
adduct fixed so the R4 perturbation is a perturbation of one thing. 139,240
records → 54,599 usable → 6,173 compounds → **1,693 eligible**, where eligible
means the compound's spectra span two instrument types *and* two contributors.
The query is the lexicographically first contributor's spectrum and every other
spectrum of that compound joins the library, so the right answer is always
present and always from another lab.

### R2 holds: 90 of 100

On the 100 graded queries, the declared pipeline annotates 98 (two have no
candidate in the window at all) and **90 are the right compound**. Eight are
confidently wrong, which is what makes R2 worth grading: the set of correct
annotations is a 90-member answer nobody guesses.

### R3's "discriminating fragment" is undefined for most queries

The catalogued R3 asks, per annotation, for the fragment that discriminates it.
That presumes a competitor to discriminate against. Measured over 300 queries, a
runner-up **of a different compound** exists for only **129** — for the other 171
the precursor window holds nothing but that compound's own spectra, and the
question has no answer. On the 100 graded queries it is 28.

So the fragment claim is dropped and R3 grades two things that are defined
everywhere: the **confidence census** under a declared mechanical rule (score
≥0.80 and ≥6 matched peaks → confirmed; ≥0.50 and ≥3 → tentative; else
unsupported), which comes out 62 / 23 / 13 and discriminates properly; and the
**ambiguous set** — annotations whose nearest different-compound competitor sits
within a declared 0.05 margin — which is 7 of the 28 that have a competitor at
all. Where the catalogue wanted a per-annotation fragment, the data supports a
set.

MassBank's own `CONFIDENCE` field cannot do this job either, and the catalogue
already knew why: 95 of 139,240 records carry it. One contributor's convention,
not a corpus-wide field.

### R4's adduct shift makes the search silent, not wrong

Re-reading a protonated precursor as sodiated moves the candidate window 21.98 Da
down. Measured on 300 queries: **246 have no candidate at all**, 54 get a wrong
top-1, and **0 survive**. On the 100 graded queries: 73 / 27 / 0.

"Predict which annotations survive" therefore grades one answer — none — which is
the same degeneracy that retired T-L6-3 and the catalogued homology R4. But the
*census* is a real answer requiring the whole perturbed search, and the finding
in it is worth a rung: a misassigned adduct does not produce wrong
identifications, it produces silence, and a pipeline that reports "no match"
without reporting that its precursor assumption was never tested has told you
nothing.

So R4 grades the adduct census **and** a second axis the grounding found to be
the gradeable one: the tolerance grid. Across 0.005 to 0.1 the aggregate top-1
rate moves one point — 0.877 to 0.890 over 300 queries, which would suggest
pinning anything in the range — while **26 of 300 individual answers move
somewhere in it**. On the graded hundred the movement sets are 9, 3, 6 and 6
queries at tolerances 0.005, 0.02, 0.05 and 0.1 against the 0.01 baseline. Four
non-empty sets of named queries: a counterfactual with per-query content, which
is what R4 is for.

### What is committed, and what a half-built campaign may not do

The campaign directory holds its inputs, its pinned rules and nothing else yet:
no gold, no grading spec, no stub ladder. That state turned out to be one the
release packager handled wrongly — it ships every directory with a `task.yaml`,
so it packaged a campaign with no rungs to clear and no score to earn. The
criterion is mechanical now and needs no status flag anyone has to remember:
**no `grading.json`, not a campaign yet.** The manifest names what it left out,
`unfinished()` is the other half of `campaigns()`, and a test asserts both
directions. Being unfinished is not a licence exemption either — the MassBank
spectra are in the repository already, so the `license_notice` travels with the
data rather than with the grading spec, and a second test holds that.

| what the catalogue asked for | what the data supports |
|---|---|
| R3: confidence level per annotation | the confidence census under a declared rule (62 / 23 / 13) |
| R3: the discriminating fragment | **dropped** — no competing compound exists for 72 of 100 queries |
| R4: which annotations survive an adduct shift | **the census** (73 no candidate / 27 lost / 0 survive); the per-query bit has one answer |
| — | **added**: the tolerance-grid movement sets, 9 / 3 / 6 / 6 of 100 |

## T-L4-1 built, and its one-sentence claim was false

The grounding above left the pipeline, the corpus and all four rungs decided, so
this was meant to be assembly: gold, the warm bundles, the grading spec, the
measurer, the oracle, the stub ladder, the rungs in `task.yaml`. Five things came
out of it, and the first one is the campaign's whole construct.

### The split did not do what the campaign said it did

What makes a dereplication campaign gradeable without anybody's judgement is
that the query's compound is in the library through *another lab's* measurement.
That sentence was in the template's docstring. Audited against the shipped
hundred, it was false for a tenth of its own answers.

The grounding's rule was: the query is the lexicographically first contributor's
spectrum, and **every other spectrum** of that compound joins the library. The
second half is the defect. A lab that contributed three spectra of a compound
keeps two of them in the library, so the query can be answered from its own
lab's replicate. Measured: **19 of the 100 graded compounds** had a
same-contributor spectrum in the library, and **10 of the 90 correct answers**
came from one — at the same contributor *and* the same instrument type. Not a
retrieval. A lookup with one extra accession number.

The rule is now: the library takes each eligible compound's spectra from **every
contributor other than its own query's**. Eligibility already requires two
contributors per compound, so no query loses its answer, and the library goes
from 30,749 spectra to **24,794**. Of the ten replicate-won answers, **nine still
retrieve correctly from another lab** and one does not: the whole cost of the fix
is 90 → 89, which is also the measure of how little the replicates were doing
except flattering the campaign.

The catalogue had it right, which is the uncomfortable part. `docs/catalog.md`
describes the split as "every spectrum **from other contributors** as the
library"; the template's rule drifted to "every other spectrum of it" somewhere
between the catalogue row and the code, and the docstring kept asserting the
catalogue's version. Neither the gate nor the oracle can see a drift like that —
both sides are internally consistent and gold is exactly what the code computes.
What found it was reading the shipped records: a query and a library spectrum one
accession apart, same contributor, same instrument.

It is checked in two places, because the version that shipped was a sentence:
`prepare_input` refuses a split where any graded query is answerable from its own
contributor, and `Corpus.load` refuses a corpus where it holds — so a hand-built
or re-cut corpus cannot reintroduce it. A test builds exactly that corpus and
asserts the refusal.

And the claim that is *not* made is the instrument. Contributor exclusion does
not imply a different instrument type, so the campaign stops claiming one and
reports it instead: `inputs/provenance.json` carries an `instrument_support`
census, which comes out **100 cross-instrument / 0 same-instrument-only**. The
strong version of the claim is true on this corpus — and it is a fact about the
corpus, recorded as a census, not a rule the builder enforces.

### The numbers as built

| | grounding's split | as shipped |
|---|---|---|
| library spectra | 30,749 | **24,794** |
| candidate set, median / max / total | 18 / 159 / — | **16.5 / 159 / 2,778** |
| annotated / no candidate / no match | 98 / 2 / 0 | **98 / 2 / 0** |
| correct / misannotated | 90 / 8 | **89 / 9** |
| confirmed / tentative / unsupported | 62 / 23 / 13 | **61 / 23 / 14** |
| queries with a competing compound | 28 | **27** |
| inside the 0.05 margin | 7 | **7** |
| tolerance movement, 0.005 / 0.02 / 0.05 / 0.1 | 9 / 3 / 6 / 6 | **10 / 3 / 6 / 6** |
| adduct: no candidate / lost / survives | 73 / 27 / 0 | **74 / 26 / 0** |

Every shape the grounding found survived the recut, which is the useful part of
the table: the split fix cost one correct answer and moved nothing about what the
rungs measure.

### R1 was a free component in its first draft

The first grading spec graded `corpus_census` — graded queries, library spectra,
library compounds. All three are stated in `inputs/provenance.json`, which the
sandbox ships because the selection rule and the per-file digests belong with the
data. A component for them measures a copy, not a parse.

What R1 grades instead is the **candidate census**: the median, maximum and total
candidate-set size the declared precursor window leaves over this library. It
cannot be read anywhere — it is the filter applied to both files — and it is
answerable **without scoring a single peak**, which is the right shape for a
first rung and is why the stub ladder's `reader` level never calls the
similarity. `corpus_census` is still reported, and named in
`excluded_from_grading` with that reason.

### Two competence probes, and a third retired on measurement

- **`peaksfirst`** ranks candidates by matched peaks before score. A defensible
  reading of a spectral comparison and not the declared one: **65 of the 98
  annotated queries get a different top match**, only 2 of which change
  correctness. R1 is untouched — the filter is upstream of any ranking — so it
  clears R1 and fails R2, where the search is.
- **`scorefloor`** reads the declared confidence rule as a score threshold and
  drops its matched-peak floor. The single-condition slip: **17 of 98 move
  class**, 15 from `tentative` and 2 from `unsupported`, **every one of them
  upward** into `confirmed` — the direction that makes a system look more certain
  than its evidence supports. Clears R1 and R2, fails R3.
- **Dropping the declared ranking's tie-break keys** moves **nothing**. Eight
  queries have tied candidate scores and none of them at the top, so ordering by
  score alone reproduces gold exactly. It is not a fixture — a probe that fails
  nothing clears every rung and measures the harness — and the rule stays, since
  it is what makes the answer independent of the tool's iteration order on a
  corpus where a tie does reach the top. The measurement is recorded on
  `core.declared_ranking`, where the rule lives.

### One refusal was unreachable, and the test is what said so

`_assert_has_content` first refused gold where "every annotation survives the
misassigned adduct". That cannot happen: surviving everything means emptying no
candidate window, and the refusal above it already fires on that. The
parametrised refusal test could not construct a gold that reached it — which is
the only way an unreachable guard announces itself, since a guard that never
fires looks exactly like a guard nothing has triggered yet. It now refuses the
reachable degeneracy: **as many right answers after the shift as before it**,
which would make R4's census a restatement of R2's.

### Where it stands

Gate: **25 pass / 0 fail / 5 pending**, the pending five being the agent-run
checks every campaign holds open. Stub sweep: **0 gates failing**, cold depths
1 / 2 / 3 / 4 for the four ladder levels, 1 for `peaksfirst`, 2 for
`scorefloor`, 0 for the no-tool ablation, R1 clear rate 1.0 and a cold mean score
of 0.54. The suite gained **55 tests** and the cross-campaign scan is clean in
both directions over all 24 campaigns. **24 campaigns are built**, and the
release packager ships all of them.

## The library pin imported off the host, and `--verify` said 9/9

Found while starting T-L4-2, by running a grounding script with `python3 -I`:

```
ModuleNotFoundError: No module named 'requests'
  matchms/importing/load_from_usi.py, line 4
```

matchms 0.30.2 declares `requests >= 2.31.0` in its conda recipe and imports it
at import time, through `matchms/__init__` → `importing` → `load_from_usi`. The
pinned prefix did not carry it. Nor `urllib3`, `certifi`, `idna` or
`charset-normalizer`: of matchms's twenty declared run dependencies, nineteen
were in the prefix and `requests` was not.

**It imported anyway, and that is the part worth keeping.** The session container
has `requests` in `/root/.local/lib/python3.11/site-packages` — dated 2026-03-31,
installed for `conan`, nothing to do with this project — and the prefix
interpreter is the same Python minor version, so a plain `python3.11 -c "import
matchms"` finds it on the user site. Which means:

- every matchms score in T-L4-1's gold was produced through an import resolved
  partly outside the pin;
- `registry --verify` reported `matchms 0.30.2` and `9/9 tools at their pinned
  version`, because its version probe is exactly that import;
- `image/environment.lock.json`, which is generated **from the prefix**, recorded
  a 226-package closure that cannot import its own library pin;
- a container built from `image/Dockerfile` — which installs Ubuntu's `python3`
  and `python3-yaml` and nothing else from PyPI — would have had no `requests`
  anywhere, so the driver would have died on the first campaign that ran it.

A pin that only imports next to the right host is not a pin.

### Where the gap came from, and where it did not

A fresh `micromamba create` of the same specs links the whole chain: measured, in
a throwaway prefix, five packages linked and `import requests` clean. So the
single `create` in the Dockerfile is the complete path and the **staged** prefix
is the incomplete one — the one assembled over five solves between 2026-10-02 and
10-21, which the Dockerfile's own comment already warned "is not guaranteed to
equal what one solve produces". How the 10-21 matchms step came to link matchms
without `requests` is not established: that install's log is gone, and my own
first `install requests` into the prefix reported "Install: 5 packages" and left
four of them unlinked, which is the same shape but one occurrence rather than a
reproduction.

That is the argument for a check rather than a note. The cause is unknown; the
condition is one line of JSON away from being checkable.

### Three changes, each one measured

**The prefix carries the closure.** `requests 2.34.2` plus `urllib3`, `certifi`,
`idna` and `charset-normalizer`, installed with exact `version=build` specs.
Nothing was upgraded, downgraded or removed — the dry run said so and the result
confirms it: numpy stays 2.2.5 and scipy 1.15.2, which is what makes the next
paragraph possible. The lock is regenerated from the prefix: **231 packages**,
seven build steps, and the five new `requested` entries carry their rationale.

**The probes no longer read the host.** The library version probe and the
invariance suite's `cosine_top1` invocation both run `python3.11 -s`, which drops
the user site. `-s` and not `-I`: `-I` implies `-E`, and the image sets
`PYTHONHASHSEED=0` on purpose. The dereplication campaign's own two calls — the
version query and the scoring driver — carry `-s` for the same reason, so a score
comes from the pinned prefix or not at all. **The invariance digest is unchanged
by all of it**: `9793db03e82f42f2` before and after, so `-s` costs nothing
numerically and `image/thread_invariance.json` stands as measured.

**`--verify` checks the closure.** `dependency_gaps()` reads each library pin's
own `depends` list out of the prefix's `conda-meta` and reports any package the
prefix lacks — the check that would have caught this on the day it happened,
and the one the Dockerfile's build-time `--verify` now runs. A binary pin needs
no equivalent: a binary either runs or it does not, and `--verify` runs every
one. A library is imported, and an import searches paths nobody pinned.

### What it did not change

T-L4-1's gold is byte-identical after all of it — same gate, 25 pass / 0 fail /
5 pending, including the regeneration check that rebuilds gold and compares it
byte for byte against what is committed. The numbers in the campaign were right;
what was wrong was the claim that only the pinned image could produce them.

## The queries carried their own answer

Found while designing T-L4-2's inputs, which must not ship the identity of the
spectra they ask about: T-L4-1's already did.

`inputs/queries.json.gz` carried `compound_block` and `formula` on every query
record, while `inputs/provenance.json` said, in the campaign's own words:

> No annotation of the queries: the compound each query IS is the answer.

The first block of an InChIKey **is** the compound, and the molecular formula
names it almost uniquely inside a precursor window of ten millidaltons. So the
answer to "which compound is this spectrum?" was in the sandbox, on the record
being asked about.

### What it cost, exactly

Not a cleared rung, which is why no gate check saw it. R2 grades the per-query
top match, and a system exploiting the leak — score only the library spectra of
the query's own compound, report the best one — reproduces gold for the 89
queries whose retrieval is correct and differs on the 9 that are not, so the
exact-match component still fails. The leak handed over 89 of 98 answers without
clearing the rung that asks for them.

It also made one component unaskable-but-asked. `correct_queries` — the set of
queries whose top match is a spectrum of the query's own compound — cannot be
computed without knowing what the query is. Grading it **required** the
disclosure. That is the tell worth remembering: a component whose answer needs
the truth to be in the sandbox is not a component, it is a reason the sandbox
holds the truth.

### The fix, and the knock-on it forced

A query record now carries what an experimental spectrum carries — accession,
contributor, instrument type, precursor m/z, precursor type, peaks — and nothing
else. The identities go to `gold/queries_truth.json`, inside a secret campaign
path that no sandbox copies and no release ships. `library.json.gz` is unchanged
and byte-identical: a reference library without identities is not a reference
library.

Three refusals keep it that way, in the grading core rather than in prose:

- a query record carrying `compound_block`, `inchikey`, `formula`, `smiles` or
  `compound_name` is refused outright;
- a truth file that does not cover every shipped query is refused, so a stale one
  cannot silently make an answer undefined;
- the split's two refusals — the query's compound is in the library, and never
  from the query's own contributor — now run where the truth is readable, which
  is the build and the oracle. That is exactly where they matter: a system cannot
  reintroduce a split defect it cannot see.

R2 lost a component with the disclosure (12 graded components, now 11) and grades
what the search answered: the retrieval census and the per-query top match. Gold
still records `correct_queries` and `misannotated_queries` for the audit, and no
component reads them.

**And R4's census had to be redefined, which improved it.** The adduct census
counted `survives` = "still names the right compound", which needs the truth and
was therefore only computable because of the leak. It now compares the shifted
search against the pipeline's **own unshifted answer**: `no_candidate`,
`changed`, `unchanged`. Measured, the two readings coincide on this corpus —

| | graded, against the pipeline | gold-side, against the truth |
|---|---|---|
| window emptied | 74 | 74 silenced |
| kept a candidate, different answer | 26 | 26 wrong |
| kept a candidate, same answer | 0 | 0 still correct |

— so the campaign's claim is unchanged and is now one a system can check: every
query that keeps a candidate under the wrong adduct names a different compound,
and none of them is right. Gold keeps the truth-side census as
`adduct_against_the_truth`, graded nowhere, because the prose rests on it and a
claim in prose should have a number an auditor can find.

### Two defects of one shape, in one campaign

The split's claim and the identity claim were both sentences in `provenance.json`
that the shipped artefact contradicted, and both were written by the same hand in
the same week. Neither the gate nor the oracle can see that class of defect: the
prose is not an input to anything, and gold is whatever the code computes. What
found them was reading the shipped records — two accessions one apart with the
same contributor, and a `compound_block` key where the provenance said there was
none.

So both are refusals now. A declaration that no code reads is a decoration, and
this campaign shipped two.

## T-L4-2: the network was buildable, two of its four rungs were not

The catalogued ladder: R1 a network edge table · R2 matches gold · R3 which
annotation is unsupported and why, from a closed confound vocabulary · R4 the
planted annotation is removed, predict the network change. Grounded against
MassBank 2026.03 before any oracle was written, R1 and R2 hold as catalogued and
the other two do not.

### Choosing a node set is choosing the campaign

A molecular network is one experiment: spectra acquired together, most of them
unidentified. So the node-set rule is one contributor on one instrument type, and
the choice of which was a survey of all 38 contributors in the release.

| node set | spectra / compounds | edges at 0.7/6 | analogue edges | largest component | verdict |
|---|---|---|---|---|---|
| **AAFC / LC-ESI-ITFT** | 478 / 121 | 676 | 383 | 61 | **built** — aflatoxins, AAL toxin, mellein: fungal secondary metabolism |
| Eawag / LC-ESI-ITFT | 3,256 / 299 | 10,123 | 3,667 | 110 | the second instantiation: cyanobacterial peptides, and four components over 90 nodes, which needs the neighbour and component-size caps this template does not declare |
| Antwerp_Univ | 688 / 169 | 1,840 | 1,060 | 153 | rejected: citrate and trimellitate plasticizers, and one 153-node hairball |
| BGC_Munich | 557 / 195 | 455 | 343 | 38 | rejected: free amino acids — primary metabolism, not natural products |
| Washington_State_Univ | 544 / 228 | 150 | 70 | 10 | rejected: on-theme alkaloids, too few edges |
| MPI_for_Chemical_Ecology | 249 / 62 | 266 | **22** | 12 | rejected: almost every edge joins replicates of one compound |
| RIKEN | 4,658 / 828 | — | — | — | rejected: pantothenate, glycerol, shikimate — primary metabolites |

The library is every OTHER contributor's spectra, which is what makes an
annotation cross-laboratory, filtered to the declared precursor window: 1,195 of
54,101 cross-contributor spectra, because a spectrum outside every node's window
cannot be a candidate for any node by the definition of the filter. That keeps 13
MB of unreachable spectra out of the repository and changes no answer.

### R3's confound vocabulary does not survive the data

The catalogued question is which annotation is unsupported and why. Five
mechanical clauses were written and measured over the accepted annotations:

| clause | fires | of which wrong | precision |
|---|---|---|---|
| `single_lab_support` | 112 | 29 | 0.26 |
| `unopposed_in_the_window` | 53 | 9 | 0.17 |
| `contradicted_by_a_neighbour` | 12 | 4 | 0.33 |
| `contradicted_across_a_near_identical_edge` | 4 | 2 | 0.50 |
| nothing fires | 10 | 2 | 0.20 |

The base rate of wrong annotations is 0.24, and the union of the clauses flags
**121 of 131** — recall 0.94 at precision 0.25. A label that applies to
ninety-two per cent of the annotations is not a finding, and the two clauses that
beat the base rate cover a dozen cases between them.

**And the one sound-looking inference is false.** "Two near-identical spectra at
the same precursor must be the same compound" is the argument a network adds to a
library search. There are four such disagreeing pairs on this node set, and **two
of them are genuinely different compounds** — one pair at cosine 0.909 with 26
matched peaks, both annotations correct, different InChIKey skeletons. Isomers
fragment alike. So the inference that would have made the confound vocabulary
work is the one the data refuses.

### What replaced it: what the network CLAIMS

Not whether an annotation is right — the nodes ship without their identities, so
that is gold's business — but what the declared rule says about the nodes the
library could not name. Four cases, and the distinction between the middle two is
the campaign:

| claim | nodes | what the truth says |
|---|---|---|
| `annotated_directly` | 121 | 94 right, 27 wrong (not graded) |
| `no_annotated_neighbour` | 304 | — |
| `analogue_with_mass_shift` | 44 | **44 of 44** are indeed a different compound, as the claim says |
| `same_compound` | 7 | 3 right, **4 wrong** |
| `conflicting_neighbours` | 2 | — |

That asymmetry is the campaign's subject. An analogue claim — "this node is a
relative of that compound, shifted by this mass" — is borne out every time it is
made here. A same-compound claim is wrong more often than right, and the naive
rule that networking practice invites, inheriting identity across every edge
regardless of the precursor difference, is wrong on **51 of the 53** nodes it
reaches. The declared rule separates the two cases; the probe in the stub ladder
conflates them.

### R4's planted annotation is retired twice over

Nothing in this benchmark ships a planted artefact — the census that retired
`structural_elements` found eleven campaigns claiming `planted` and none doing it
— and **removing an annotation cannot change a network**: the graph is built from
spectra. What a removal moves is the annotation layer. So the counterfactual is
two declared re-cuts, and both move:

| threshold | edges | components | largest | nodes the annotations reach |
|---|---|---|---|---|
| 0.6 | 902 | 160 | 97 | 55 |
| **0.7** | **676** | **196** | **61** | **51** |
| 0.8 | 410 | 259 | 26 | 30 |
| 0.9 | 190 | 344 | 21 | 18 |

and dropping the lower accepted confidence tier removes 48 annotations, changes
every claim count and leaves **28 components with no annotation at all**.

The pinned threshold's own row is the grid's identity rung: it must reproduce the
network the rungs above report, and the build refuses gold otherwise.

### The edge pairs are graded and the edge scores are not

A cosine is a float out of a BLAS-backed library. What the campaign claims is
which pairs clear the declared threshold, which is a comparison — and measured on
this node set, no pair score sits within 1e-4 of 0.6, 0.7 or 0.8, and exactly one
sits 5.6e-05 below 0.9. That is the closest call in the campaign and still three
orders of magnitude above the noise this image's pinned numpy and scipy produce,
which is why the pair set is safe to grade and the scores are reported and read
by nothing.

### Where it stands

Gate: **25 pass / 0 fail / 5 pending**. Twelve graded components over four rungs.
Two similarity calls serve every rung, including the whole R4 grid: the pairwise
matrix is computed once at the lowest declared threshold and every higher row is
a filter over it.

Stub sweep: **0 gates failing**, 84 runs and no crash, cold depths 1 / 2 / 3 / 4
for the four ladder levels, 1 for `seedall` (seeds from every library match
rather than the declared tiers, so R2's component census fails), 2 for
`propagate` (inherits identity across every edge, so R3 fails), 0 for the no-tool
ablation. R1 clear rate 1.0, cold mean score 0.54.

A third plausible slip is recorded rather than fixtured: dropping the declared
matched-peak floor from the edge rule takes the network from 676 edges to 822. It
fails R1, and a fixture that fails the first rung drags the harness's own R1
clear-rate gate below its floor and reports a fixture artefact as a campaign
defect.

## The sweep passed every gate while twelve runs crashed

Worth recording because it was my own typo that found it, and because the hole it
found is one a real campaign could fall into.

A stub level called a template function whose signature had changed one commit
earlier. Twelve of eighty-four runs crashed, the `complete` fixture produced no
submission at all, R4 was cleared by nobody — and the sweep reported **0 gates
failing**. Each gate was behaving as designed:

- a crash is an INFRA failure, which the headline deliberately excludes so that
  the score is not measuring our own container;
- a rung no system clears is a legitimate finding about a hard campaign, so
  `discrimination` named it in its detail text (`rungs never cleared by any
  system: ['r4']`) and passed.

Both of those are right for a sweep over real systems and wrong for a sweep over
fixtures. A stub is deterministic code in this repository: it cannot crash for an
environmental reason, and the top fixture exists precisely to show that the
ladder is clearable. So a stub sweep now carries two gates of its own:

- **`stub_systems_do_not_crash`** — any crashed or timed-out fixture run fails
  the sweep, with the fixture named. A fixture that crashes is a harness bug, not
  a finding.
- **`every_rung_is_cleared_by_some_fixture`** — a rung no fixture ever clears
  fails the sweep. With real systems that would be a difficulty measurement; with
  fixtures built to reach preset depths it is a defect in the fixture or the
  campaign.

Neither applies to a sweep without stubs, which is the sweep whose numbers get
published.

## Check 26: three campaigns were shipping their verdict in one word

Found while cleaning up metadata before the freeze, by asking a question none of
checks 21, 22 or 23 asks: **where a graded component's whole answer is one term
of a published vocabulary, is that term anywhere in the cold sandbox?**

Checks 21 and 22 pair a graded value with a label — the label inside the value,
or the label in the gold PATH. Neither reaches a term stated against some
*other* label. So this passed every scan, in `template_params`, for as long as
the campaign had existed:

```yaml
conflict:
  forbidden_codons: [CCC, CCG, CGC, CGG, GCC, GCG, GGC, GGG]
  gc_global_min: 0.55
  expect: feasible          # <- R3's entire answer
```

### What the scan found, across all 25 campaigns

Eighteen components whose whole answer is one declared term. Three campaigns were
genuinely stating it:

- **the constraint control pair** — `expect:` in `template_params` on both
  halves, plus prose on each (`"The correct answer to this campaign is that no
  design exists"`, `"This is the feasible half of a control pair"`, a chance-level
  note naming the most guessable verdict, and an `excluded_from_grading` reason
  copied from the wrong half);
- **`massbalance-nrps-sevadicin-01`** — two audit notes and a chance-level note
  stating both the verdict and the residual it follows from.

And one campaign was grading two **echoes**: `ecaudit`'s `canonical_model` and
`threshold` are declared verbatim in `reference/`, because the campaign cannot be
posed without naming the domain model it audits for and the cutoff to run it at.
Both are now excluded from grading with the reason recorded, which is the same
treatment `release` and `ec_number` already had in that campaign.

The remaining twelve are collisions, declared in
`audit/cold_sandbox_exemptions.json` with the collision named: a verdict
PRECEDENCE rule has to name the terms it orders (`mibig-diff`, both mass-balance
instantiations), the six construct constraint NAMES are the published contract
and one of them is also R3's answer, and a corpus statistic that one of 43
admissible entries is exactly balanced identifies none of them.

### The campaign names were part of the leak

`constraint-conflict-ecoli-rebh-01` and `constraint-feasible-ecoli-rebh-01`. The
campaign id is in `task.yaml`, which every sandbox gets, so the directory name
answered R3 before the agent read a single constraint. They are
**`constraint-gc62-ecoli-rebh-01`** and **`constraint-gc55-ecoli-rebh-01`** now —
named for the parameter under test, which each objective already states.

### And the pair could not say that its halves disagree

The subtler one. Both halves wanted to explain why a control pair exists, and
every honest phrasing of it ("a system that always answers one way scores full
marks on one half and zero on the other") states that the two verdicts differ.
A higher GC floor is never easier, so *that sentence answers both campaigns at
once* — no computation required. Both halves now carry the same
`control_pair_note`, which says the pair exists so abstention can be scored as
precision and recall, and says nothing about the relationship. Identical prose in
both halves is the structural fix: identical text cannot tell you which half you
are holding.

### What the fix cost, and what it did not

The verdict is now derived from the arithmetic rather than from the authored
expectation (`"feasible" if reachable else "infeasible_gc_unreachable"`), the
expectation moved to `gold/expectation.json` outside every sandbox, and the build
refuses a `task.yaml` that declares it or a campaign that has no expectation file
at all. **Both halves' gold is byte-identical after all of it** — the answers were
right, and what was wrong was that the sandbox told you them.

## Release 1.0.0

Generated, not committed: `npbench_c.release.package` builds the shippable tree
and refuses to write it inside this repository, so that the next packaging run
cannot package its own output. What is committed is the record --
`release/MANIFEST-1.0.0.json` (every shipped file with its sha256, the grader
version, the campaigns included, the paths withheld),
`release/ATTRIBUTION-1.0.0.md` (assembled from each campaign's own
`license_notice`) and the tarball's digest. `release/README.md` carries the
command that regenerates it and a script that verifies a tarball against its
manifest.

| | |
|---|---|
| campaigns | 25 |
| files | 171 |
| grader | `2.0.0+9fc95658` |
| size | 43 MB unpacked, 37 MB packed |
| withheld per campaign | `gold/`, `oracle/`, `oracle_submission/`, `grading.json`, `internal_sweep.json` |
| excluded trees | `docs/`, `src/`, `tests/`, `image/`, `audit/`, `oracle/` |

A campaign in the release is `inputs/`, `reference/`, `task.yaml` and an empty
`submission/` -- the cold sandbox and nothing else. Verified on this build: zero
files matching `*truth*`, `expectation.json`, `gold*`, `grading.json`,
`internal_sweep.json` or `measurement.json` anywhere in the tree, and the two
MS2 campaigns' shipped query and node records carry only the six fields an
experimental spectrum carries.

**What the release is NOT yet.** Every campaign reports `ready_to_audit: False`,
for one reason: the five agent-run checks are satisfied only by a sweep over real
systems, and every sweep on record is stub-based. The gate will not promote a
fixture's number to a measurement. So 1.0.0 is a complete, mechanically verified
benchmark whose difficulty and discrimination are declared and unmeasured -- see
open item 3, which is the whole of what remains.

## The missing half of the sweep: a prompt, and a declared system set

The readiness gate has kept five checks PENDING on every campaign since the
harness existed, and the reason was recorded as "needs real systems". That was
true and incomplete. The runner has always been agent-agnostic -- a system is a
command, `SystemSpec` is a name, a command, a mode and a timeout -- so what was
missing was never the runner. It was two things, and neither existed:

**A stub reads `task.yaml` itself.** It is code in this repository, so it needs
no instructions. A real agent does, and nothing in the harness produced any. A
sweep over real systems would have handed each one an empty working directory and
measured what it guessed.

**A set of systems to run.** `--stubs` was the only path the CLI had; the other
branch printed "no real systems configured" and exited 2.

### `npbench_c.sweep.prompt`

One generator for all 25 campaigns, with three properties that are the whole
point of generating rather than writing 25 prompts:

- **built from the SANDBOX, not the campaign** -- the only inputs are the files a
  system already has, so it cannot leak what `task.yaml` does not. Verified the
  same way the campaigns are: all 25 cold prompts pass checks 22 and 26's scans,
  as a test, with the prompt rendered into the sandbox being scanned;
- **identical for every system** -- a difficulty comparison between two systems
  given different instructions measures the instructions;
- **it states the contract and declines to state the method.** The objective, the
  scope note, the output contract field by field, the pinned tables, the budget
  and the tool surface are rendered. `notes_for_audit` is not: it is written for
  someone holding gold, and in front of a system it reads as a procedure.

A warm prompt additionally names the files handed over and says to take them as
correct. 100 prompts render across the suite -- 25 cold and 75 warm -- and the
first draft raised on two of them, because it assumed a warm bundle is JSON and
the construct campaigns' R2 bundle is a FASTA.

### `npbench_c.sweep.systems`

The systems come from a config file, and the module's value is what it refuses:

| refusal | why it exists |
|---|---|
| a command that never mentions `{submission}` | the run is scored on an empty directory and reads as a capability result |
| no `no_tool` ablation, or two | the leakage check is computed against it; without it the gate fails rather than passing quietly |
| no reference system, or two | the difficulty gate is judged on its per-rung profile |
| every tooled system is the reference | the R1 clear rate is measured over the tooled systems that are NOT the reference, so the gate becomes uncomputable and reports `None` as a campaign defect |
| a `no_tool` system that can still see the tool prefix | the ablation's whole claim is that it did not compute; the loader strips the prefix from its PATH and refuses a config that puts it back |
| `is_stub` | a fixture in a real sweep is laundered into evidence |

The **timeout is not configurable upward**: it is derived from each campaign's
own `budgets.wall_clock_minutes`, so no system is silently given more wall clock
than the campaign budgeted. A config may ask for less.

`sweep_systems.example.json` is a template; `--systems <config>` is the CLI path.

### And a harness defect only a real command could find

`_format` substituted the harness's placeholders with `str.format`. A stub's argv
is written in this repository and never contains a brace; a real agent's command
plausibly contains several, and the smoke test's `printf '{}' > {submission}/r.json`
raised before it reached a shell. Substitution is explicit now -- four names,
every other brace left alone -- with a test for `${HOME}` and a JSON literal in
one command line.

**Smoke-tested end to end**: trivial systems against one campaign, cold plus
three warm, prompts written, submissions measured, graded, and three sweep gates
correctly failing because a system that writes `{}` clears nothing. The path
works; what it needs now is an agent.

### The near-miss that flag exists for

`sweep_campaign` records into the campaign by default, which is right for a real
sweep and wrong for a trial. The first smoke run wrote its result into
`constraint-gc62-ecoli-rebh-01/internal_sweep.json` marked `stub_based: false` --
and the readiness gate READS that file and promotes a non-stub sweep to a
measurement. Two throwaway systems that write `{}` had become the campaign's
evidence, with three gates failing, and it was committed before it was noticed.

Restored from the previous commit, and `--no-record` exists now so a trial run
cannot touch the record at all. A test holds the flag by asserting what the CLI
passes through, because the failure mode is silent: the sweep prints the same
summary either way, and the damage is in a file nobody re-reads until the gate
does. Hardened on the way past: the summary printer raised on a report with no
gates.

## Open items

1. ~~CAI table provenance.~~ **Closed.** The tables are derived from two pinned
   INSDC records by a declared ribosomal-protein rule, the reference set is
   hashed and enforced, the per-codon counts and gene lists ship in
   `audit/codon_table_provenance.json`, and the rebuild moved 37 gold keys across
   the two design campaigns and none in the two constraint campaigns. See "Where
   the codon tables come from" above, including the two stale exemptions and the
   three disclosures it surfaced on the way.
2. ~~Release packaging must exclude `gold/`, `oracle/` and
   `oracle_submission/`.~~ **Done**, and it is a program rather than a sentence:
   `npbench_c.release.package` builds the tree, refuses one that carries gold,
   the oracle, the grading spec, `docs/` or any graded answer (re-running the
   gate's check-22 scan against what would ship), and assembles the attribution
   from each campaign's own statements. See the section above for the three
   exclusions that are less obvious than `gold/` and for the one campaign whose
   missing licence record it found.
3. **Real systems for the sweep — the only thing between 1.0.0 and an auditable
   benchmark.** The harness side is now complete: the prompt generator, the
   validated system-set loader, the `--systems` CLI path and an end-to-end smoke
   test are all in (see "The missing half of the sweep" above), and
   `sweep_systems.example.json` is a template with three roles filled in. What
   remains is not code: it is **which agents to run and the credential and
   compute to run them**. The arithmetic, at the frozen 25 campaigns: 4
   rung-modes x 25 campaigns = 100 runs per system per repeat, so three systems
   at two repeats is 600 agent sessions. Expect some campaigns to FAIL
   difficulty or discrimination on real agents and budget a rescope pass for
   them -- that failure is the benchmark working.
4. ~~matchms in the image.~~ **Done.** The prefix carries
   matchms 0.30.2 (46 packages added, nothing moved), `--verify` reports 9/9,
   the invariance suite measures it as an invocation like every other tool
   (raw-identical, plus a separate BLAS-thread measurement), the lock is
   generated from the prefix at 226 packages with a build history, and the
   fixture is cut from MassBank's released tag with the tolerance table
   re-measured on this image's numpy. Two defects came out of it: the grounding's
   "MassBank 2026.03.1" was a dev-branch snapshot string that no clone can
   reproduce, and a library pin probed through PATH reports itself missing on an
   image that has it. See "Item 4" above. **T-L4-1's first instantiation is
   now built** — see "T-L4-1 built" above, including the split defect the
   assembly found — and **T-L4-2 is built too**. T-L4-1's second instantiation
   is the only matchms row left, and it needs an instantiation axis decided
   (item 14).

5. ~~Reprocess the MIBiG reference set with the pinned antiSMASH, or declare
   that a T-L3-2 campaign takes both sides from the published `as8b1` set.~~
   **Declared, and now enforced.** Counted over all 23 campaigns, the group that
   runs antiSMASH (one campaign, on a raw INSDC chromosome) and the group that
   ships antiSMASH-processed records (four, all driving HMMER or BiG-SCAPE) do
   not overlap, so no comparison in the suite mixes versions and reprocessing
   would only cost the provenance of a published set. The four now record
   `antismash_version_that_processed_them` as a field, and check 25 refuses a
   campaign that ships processed records without it or runs a different antiSMASH
   over them. See "Item 5" above.

6. **T-L5-2 is a decision, not a build task.** T-L5-1 is built. T-L5-2 has no
   on-theme data (counts above); the owner picks between HIV-1 RT with the theme
   drift, a re-scoped census ladder, and waiting for a release that fills the
   bacterial variant values. Nothing is blocked on licensing.
7. ~~The published count.~~ **Decided: the suite is 25 campaigns from 20
   templates, frozen there** (owner, 2026-10-27). Restated in this file's
   headline, in `docs/catalog.md`'s title and in `docs/environment.md`, with the
   difference from the catalogue recorded rather than dropped: three rows ground
   out on measurement (T-L2-4 ×2, T-L5-2, T-L6-3, numbers above), and the two
   remaining second instantiations — T-L4-1 ×2 and T-L4-2 ×2 — are buildable and
   deliberately not built (items 14 and 15). T-L3-5a counts as one of the 25: it
   is a built campaign whatever it replaced. 25 campaigns × 4 rungs = 100 graded
   rungs carrying 560 graded components.
8. ~~Sweep cost: the tool memo covers the tool, not the work around it.~~
   **Closed, and the diagnosis in this item was wrong.** Profiled per phase, the
   sandbox copy it blamed costs 0.01s and the parse was only one of three causes:
   a `.tar.gz` read out of order (5.3s of a 7.6s corpus load), a 1.57 GB Pfam
   file hashed once per run to report a digest `reference/` already pins (4.5s),
   and a tool asked its own version twice per run where that tool imports a
   Python package to answer (2.12s and 2.38s for the only two of nine that do).
   RiPP 7m51s → **45s**, cutoff sweep 6m27s → **83s**, bgcdetect 15m35s →
   **8m48s**, with every `internal_sweep.json` byte-identical and every gate
   still green. What remains in a sweep is the distinct tool invocations
   themselves. See "Item 8" above.

9. ~~`structural_elements` is decorative and nothing validates it.~~
   **Closed by removal**, with the census that justified it: two of its seven
   terms were declared by 23 of 23 campaigns, `planted` was declared by 11 and no
   campaign ships an altered input, `control` meant three unrelated things, and on
   the constraint control pair the label stated which half of a two-verdict pair
   the reader was holding. The field is gone from every `task.yaml`, check 24
   keeps it out under five names and validates the design references that do
   resolve, and the two tests that asserted the label now assert the artefacts.
   See "Item 9" above.

10. ~~Twenty graded components restate the campaign's own definition.~~
   **Closed.** All twenty are dropped from grading across nine campaigns, each
   with its exclusion recorded in the campaign's `task.yaml`, and
   `SCALAR_FREE_COMPONENTS` in `tests/unit/test_cold_sandbox_prose.py` is now
   empty with a test asserting it stays that way. See "The copies that were
   being graded" above for what they were and what each rung kept.

11. ~~The categorical-leak check belongs in the readiness gate, as check 22.~~
   **Done.** `cold_sandbox_withholds_small_answers` is the gate's twenty-second
   mechanical check, which took the gate to 22/0/5 at the time. The scan lives in
   `src/npbench_c/readiness/cold_sandbox.py`, the gate and
   `tests/unit/test_cold_sandbox_prose.py` run the same implementation, and the
   test covers all 22 campaigns in two seconds where the gate covers one at a
   time. The migration cost was an afternoon's worth of gate runs, as expected;
   what made it affordable was discovering that a gate run costs minutes rather
   than the hours assumed when the item was written.

12. Construct-validity study: inter-rater agreement first, then expert-grader
   agreement with Gwet's AC1 / Krippendorff's alpha alongside kappa, gate on
   Spearman against the continuous rating. This is the one place humans are
   involved, and it sits outside the grading pipeline by design.

13. ~~A cross-campaign disclosure scan.~~ **Done**, as the gate's
   twenty-third mechanical check, which took the gate to 23/0/5 at the time
   (item 9 below added a twenty-fourth). Which shapes
   it runs was measured over all 506 ordered pairs first, and the small-scalar
   shape is deliberately excluded with the number that excluded it. It found
   four campaigns stating the gene-calling campaign's graded GC fraction — a
   provenance block three days old — and one false positive that tightened the
   scalar boundary both checks now share. Scope is the set shipped together:
   the gate passes every built campaign, `audit_release` passes the release's
   own subset. See "Check 23" above.
14. **Not built, by the freeze.** T-L4-1's second instantiation needs an axis,
   and it is the owner's call.
   The catalogue names no instantiation axis for the dereplication row, and the
   alternative sources it lists are not interchangeable with MassBank: NIST14 is
   unshippable, the GNPS reference-standard sets and Bruker MetaboBASE need
   licence checks, and CASMI 2022 is a difficulty anchor rather than a library.
   Within MassBank the honest axes are a different assumed adduct -- which also
   re-points R4 -- a negative-mode corpus, or a different release tag for a
   reproducibility-style pair. Nothing is blocked on code: the template takes an
   entry file and two build commands.
15. **Not built, by the freeze.** T-L4-2's second instantiation is a build
   rather than a decision, and the node set is chosen: Eawag on LC-ESI-ITFT, 3,256 spectra of cyanobacterial peptides.
   It needs the two caps this template deliberately does not declare -- a
   neighbour cap and a maximum component size -- because four of its components
   hold more than ninety nodes, and the build refuses a hairball on purpose. So
   the work is declaring both caps with the measurement that sets them, not
   another corpus hunt.
