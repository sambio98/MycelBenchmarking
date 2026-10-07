# NPBench-C — build notes (Phase 0 / 0.5)

Benchmark for natural-product and BGC computational agents. 32 campaigns from
24 templates, four-rung ladders, deterministic grading, no LLM anywhere in the
scoring path.

## What exists

| Component | Path | State |
|---|---|---|
| Scoring primitives | `src/npbench_c/grading/primitives.py` | 12 primitives, stdlib-only |
| Composition policy | `src/npbench_c/grading/compose.py` | product / geometric+floor / min |
| Rung ladders | `src/npbench_c/grading/ladder.py` | depth, monotonicity, chance floor |
| Grader | `src/npbench_c/grading/grade.py` | pure declarative comparator |
| Grader identity | `src/npbench_c/grading/version.py` | semver + content hash |
| Readiness gate | `src/npbench_c/readiness/gate.py` | 20 mechanical checks, 5 pending-on-agent-runs |
| Template: construct design | `src/npbench_c/templates/construct_design/` | shared engine, 2 ladders |
| Campaign: construct design | `campaigns/construct-ecoli-rebh-01/` | oracle 1.0, gate 20/0/5 |
| Campaign: constraint conflict | `campaigns/constraint-conflict-ecoli-rebh-01/` | oracle 1.0, gate 20/0/5 |
| Campaign: feasible control | `campaigns/constraint-feasible-ecoli-rebh-01/` | oracle 1.0, gate 20/0/5 |
| Template: NRPS mass balance | `src/npbench_c/templates/mass_balance_nrps/` | shared oracle, 2 instantiations |
| Campaign: malleobactin | `campaigns/massbalance-nrps-malleobactin-01/` | oracle 1.0, gate 20/0/5 |
| Campaign: sevadicin | `campaigns/massbalance-nrps-sevadicin-01/` | oracle 1.0, gate 20/0/5 |
| Internal sweep | `src/npbench_c/sweep/` | runner, stats, gates, stub fixtures |
| Catalog | `docs/catalog.md` | 24 templates / 32 campaigns, grounded |
| Provisioning | `docs/environment.md` | tools, DBs, sandbox contract |
| Template: MIBiG release diff | `src/npbench_c/templates/mibig_diff/` | shared oracle |
| Campaign: MIBiG 3.1→4.0 diff | `campaigns/mibig-diff-3_1-to-4_0-01/` | oracle 1.0, gate 20/0/5 |
| Template: annotation audit | `src/npbench_c/templates/mibig_annotation_audit/` | shared oracle |
| Campaign: PrnA construct | `campaigns/construct-ecoli-prna-01/` | oracle 1.0, gate 20/0/5 |
| Campaign: gene-function evidence | `campaigns/mibig-gene-function-evidence-01/` | oracle 1.0, gate 20/0/5 |
| Template: structure features | `src/npbench_c/templates/structure_features/` | shared engine, 2 ladders |
| Campaign: PrnA residue evidence | `campaigns/residues-prna-01/` | oracle 1.0, gate 20/0/5 |
| Campaign: RebH pocket geometry | `campaigns/pocket-rebh-01/` | oracle 1.0, gate 20/0/5 |
| Template: chemical space | `src/npbench_c/templates/chemical_space/` | shared oracle, RDKit pinned at instantiation |
| Campaign: MIBiG 4.0 chemical space | `campaigns/chemspace-mibig-4_0-01/` | oracle 1.0, gate 20/0/5 |
| Template: kinetics consistency | `src/npbench_c/templates/kinetics_consistency/` | shared oracle, BRENDA 2026.1 |
| Campaign: BRENDA EC 1.1 kinetics | `campaigns/kinetics-brenda-ec1_1-01/` | oracle 1.0, gate 20/0/5 |
| Template: ChEMBL selectivity | `src/npbench_c/templates/chembl_selectivity/` | shared oracle, ChEMBL_37, first `share_alike` |
| Campaign: S. aureus topoisomerases | `campaigns/selectivity-saureus-topoisomerase-01/` | oracle 1.0, gate 20/0/5 |
| Template: EC domain audit | `src/npbench_c/templates/ec_domain_audit/` | HMMER 3.4 + Pfam 35.0, first tool-running template |
| Campaign: EC 1.1.3.15 FMN_dh audit | `campaigns/ecaudit-fmn-dh-1_1_3_15-01/` | oracle 1.0, gate 20/0/5, **first S1** |
| Phase 1 image | `image/Dockerfile`, `image/environment.lock.json` | 8 tools pinned to `version=build`, 185-package closure hashed, antiSMASH databases pinned at 9.4 GB |
| Tool registry | `src/npbench_c/tools/registry.py` | pins, controls, invocations, declared normalisations, canonicalisations, projections, enforced bans |
| Thread-invariance suite | `src/npbench_c/tools/invariance.py` | 8/8 tools invariant at 1 and 8 threads |
| Tests | `tests/` | 353 passing, 2 skipped |

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

Current state, every campaign: **20 pass / 0 fail / 5 pending**.

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

**It took 29m32s of wall clock and 47m of CPU**, against seconds for a stdlib
campaign. That is the whole sweep cost of one S1 campaign with a 13-model panel;
a full-Pfam panel would have been roughly fourteen times worse, and an S2 campaign
running antiSMASH per rung-repeat will be worse again. The open items carry this
as a decision to take before the first S2 build, not after.

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
4. **matchms fixtures.** The last provisioning gap, for T-L4-1 and T-L4-2. The
   tools are pip-installable; what the rows need first is an MS² fixture set from
   MassBank with its own grounding pass, because for a spectral match the
   invariance question becomes a tolerance question and that is a different suite
   from this one. (This item appeared twice in earlier revisions, as 5 and 6; the
   duplicate is merged here.)
5. **Reprocess the MIBiG reference set with the pinned antiSMASH**, or declare
   that a T-L3-2 campaign takes both sides from the published `as8b1` set. As
   shipped, the reference clusters were processed with antiSMASH 8.0 beta 1 and
   the image pins 8.0.4.
6. **T-L5-2 is a decision, not a build task.** T-L5-1 is built. T-L5-2 has no
   on-theme data (counts above); the owner picks between HIV-1 RT with the theme
   drift, a re-scoped census ladder, and waiting for a release that fills the
   bacterial variant values. Nothing is blocked on licensing.
7. **Two L2/L5 rows are retired on measurement, so the catalogue's 32-campaign
   target is now 30 at most.** T-L5-2 and T-L2-4 both ground out. The count in this
   file's header and in `docs/catalog.md` is the *catalogued* target and has not
   been restated downward, because the replacement question is the owner's call:
   re-scope those rows, or accept a smaller benchmark. Fourteen campaigns are
   built. This should be settled before the audit packet quotes a number.
8. **The S1 sweep is slow enough to need a decision — measured: 29m32s.** Every
   stub level of a tool-running campaign re-runs the tool, so
   `ecaudit-fmn-dh-1_1_3_15-01` swept in **29m32s** (47m of CPU) where a stdlib
   campaign takes seconds, and its `complete` level alone runs seven HMMER passes
   per rung-repeat. The declared panel keeps
   that affordable, but the real-systems sweep over several S1 and S2 campaigns
   will not fit the current serial runner. Either the runner caches tool output per
   (input, threshold) across stub levels, which is honest for stubs but must not
   leak into real runs, or the sweep gets parallelised. Worth settling before the
   first S2 campaign, not after.
9. **`structural_elements` is decorative and nothing validates it.** Nine of the
   fourteen built campaigns declare `planted` and several of them plant nothing —
   `mibig-diff-3_1-to-4_0-01`'s own notes say its sharpest case is "supplied by the
   corpus rather than planted". The new campaign declares only what it has
   (`verification, counterfactual, control`), but the field needs either a
   definition and a gate check or removal, and fixing the other nine is the owner's
   taxonomy call rather than a silent edit. A declaration no check reads is the
   thing this project keeps finding in other people's data.
10. Construct-validity study: inter-rater agreement first, then expert-grader
   agreement with Gwet's AC1 / Krippendorff's alpha alongside kappa, gate on
   Spearman against the continuous rating. This is the one place humans are
   involved, and it sits outside the grading pipeline by design.
