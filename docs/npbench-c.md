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
| Template: construct design | `src/npbench_c/templates/construct_design/` | shared engine, 2 ladders |
| Campaign: construct design | `campaigns/construct-ecoli-rebh-01/` | oracle 1.0, gate 17/0/5 |
| Campaign: constraint conflict | `campaigns/constraint-conflict-ecoli-rebh-01/` | oracle 1.0, gate 17/0/5 |
| Campaign: feasible control | `campaigns/constraint-feasible-ecoli-rebh-01/` | oracle 1.0, gate 17/0/5 |
| Template: NRPS mass balance | `src/npbench_c/templates/mass_balance_nrps/` | shared oracle, 2 instantiations |
| Campaign: malleobactin | `campaigns/massbalance-nrps-malleobactin-01/` | oracle 1.0, gate 17/0/5 |
| Campaign: sevadicin | `campaigns/massbalance-nrps-sevadicin-01/` | oracle 1.0, gate 17/0/5 |
| Internal sweep | `src/npbench_c/sweep/` | runner, stats, gates, stub fixtures |
| Catalog | `docs/catalog.md` | 24 templates / 32 campaigns, grounded |
| Provisioning | `docs/environment.md` | tools, DBs, sandbox contract |
| Template: MIBiG release diff | `src/npbench_c/templates/mibig_diff/` | shared oracle |
| Campaign: MIBiG 3.1→4.0 diff | `campaigns/mibig-diff-3_1-to-4_0-01/` | oracle 1.0, gate 17/0/5 |
| Template: annotation audit | `src/npbench_c/templates/mibig_annotation_audit/` | shared oracle |
| Campaign: PrnA construct | `campaigns/construct-ecoli-prna-01/` | oracle 1.0, gate 17/0/5 |
| Campaign: gene-function evidence | `campaigns/mibig-gene-function-evidence-01/` | oracle 1.0, gate 17/0/5 |
| Tests | `tests/unit/` | 153 passing |

```bash
PYTHONPATH=src python3 -m pytest tests -q
PYTHONPATH=src python3 -m npbench_c.sweep.cli campaigns/construct-ecoli-rebh-01 --stubs
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
difficulty gate, discrimination) are satisfied by an `internal_sweep.json` from
a sweep over **real** systems. A stub-based sweep leaves them pending by design.

**Never audit before calibrating.** A campaign dropped for non-discrimination
after a domain expert spent three hours on it has burned the scarcest resource
in the project. Order is: build → internal sweep → drop/fix → audit → ship.
Author ~44 to ship 32.

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
With a ladder the unit is the rung-run — 32 campaigns x 4 rungs x 3 repeats =
384 observations — and CIs come from a **cluster bootstrap resampling
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
sweep on record the gate reports **16 pass / 0 fail / 5 pending**.

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
the conflict campaign. `constraint-feasible-ecoli-rebh-01` differs from it by
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

## Four catalogued designs corrected by grounding

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

The pattern is consistent enough to be a rule: a catalogue entry is a hypothesis
about data until the fields are inspected.

## Open items

1. **CAI table provenance.** The relative-adaptiveness values now live in
   `campaigns/*/reference/codon_tables.json` (shipped to agents; the oracle
   loads the same file). They are representative values for highly-expressed
   genes, flagged in the JSON, and need a cited source release — a named
   Kazusa / HIVE-CUT release, or a table computed over a declared
   ribosomal-protein gene set — with a hash, before the audit packet ships. Gold
   does not depend on them being *the* correct measurement, only on their being
   frozen, hashed and published; but an auditor will and should ask.
2. **Release packaging must exclude `gold/`, `oracle/` and `oracle_submission/`.**
   They are committed here because this is the build repo. Gold is never
   published for either split.
3. Real systems for the sweep. The harness is built and self-tested; it needs
   `SystemSpec` entries for the actual agents (bash-only, general biomedical,
   tool-equipped) before the five checks can go green.
4. Construct-validity study: inter-rater agreement first, then expert-grader
   agreement with Gwet's AC1 / Krippendorff's alpha alongside kappa, gate on
   Spearman against the continuous rating. This is the one place humans are
   involved, and it sits outside the grading pipeline by design.
