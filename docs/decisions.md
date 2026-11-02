# Decisions the owner has to make, and the recommendation for each

NPBench-C is mechanically finished: 25 campaigns at 28 pass / 0 fail / 5 pending,
1,539 tests, release 1.0.0 verified, the tool image built end to end with a
closure that reproduces. Nothing further is blocked on code. What is left is nine
decisions, and they are listed here rather than in the open-items index because
that index records what the PROJECT owes; this records what the OWNER owes.

Each entry states the constraints already fixed by the data -- those are not open
-- and then a recommendation. The recommendations are not neutral: where one
option is better the entry says so.

Status: **none of the nine is decided yet.** Dates go in as they are.

---

## A. Blocking the internal sweep

### 1. Which systems to run, and what each one is for

- [ ] **Decided:**

The gate needs three roles, and they are scientific controls rather than
contestants:

| role | why the gate needs it |
|---|---|
| **reference** (exactly one) | the difficulty gate is judged on its per-rung profile |
| **another tooled system** (at least one) | the R1 clear rate is measured over tooled NON-reference systems; one tooled system alone makes it uncomputable and the gate reports it as a campaign defect |
| **`no_tool` ablation** (exactly one) | the leakage check compares its score against the chance floor; the loader strips the tool prefix from its PATH and refuses a config that puts it back |

`npbench_c.sweep.systems` refuses a config that breaks any of this, so the
decision is which agents fill the roles, not whether the roles exist.

**Recommendation.** Make the three differ in **tool access, not model family**:
the strongest tooled agent as reference, the same agent with domain tools
removed, and the same agent again as the ablation. Varying the model instead
measures capability and tooling at once, and the ablation stops being a control
of anything.

### 2. Pilot before committing the 600 sessions

- [ ] **Decided:**

The full sweep at the frozen suite is 4 rung modes x 25 campaigns x 3 systems x 2
repeats = **600 agent sessions**, worst case 336 h of declared budget serial, 84 h
at four concurrent lanes, peak 2 CPU and 8 GB per run.

**Recommendation.** Run **one campaign x 3 systems x 2 repeats = 24 runs** first,
on an S1 campaign with a long budget (`arch-pks-mibig-4_0-01`, 45 min). For 4% of
the cost it answers three things nobody has measured: real wall clock per run,
whether an agent survives the declared budgets at all, and whether the credential
and concurrency wiring holds. Then commit to the rest.

### 3. The drop / rescope rule -- decide BEFORE the numbers arrive

- [ ] **Decided:**

Some campaigns will fail difficulty or discrimination on real agents. That is the
sweep working. The risk is not the failure; it is choosing the remedy after
seeing WHICH campaigns failed, because the suite then gets tuned toward whatever
the reference agent happens to do well.

The thresholds already exist in the gate: R1 pass@1 >= 0.8 for tooled
non-reference systems; the reference clears both lower rungs at >= 0.5 with the
upper rungs' minimum <= 0.5; discrimination needs >= 2 distinct cold depths and a
mean score strictly between 0 and 1; the no-tool ablation must sit within z = 2.0
of the chance floor; per-rung clear rates must be non-increasing.

**Recommendation.** Commit now, in writing, to those thresholds and to one fixed
remedy: **a campaign may be rescoped once at its rung boundaries, then dropped.**
Publish the failing numbers either way. A dropped campaign with its measurements
on the record is a stronger artefact than a rescued one.

---

## B. Blocking publication

### 4. The public / private split

- [ ] **Decided:**

Constraints already in the data, not open to choice:

- both `construct-design` campaigns carry an explicit `split_constraint`: 3 of 6
  graded sub-answers coincide across them, so a system that solved one transfers
  those to the other without recomputing;
- `constraint-gc62-ecoli-rebh-01` names a `required_sibling` -- it is half of a
  two-verdict control pair and is meaningless alone;
- **five templates have two instantiations**: domain-architecture,
  construct-conflict, construct-design, mass-balance-nrps, structure-features.

**Recommendation.** Split **by template, never by campaign** -- both
instantiations of a template go the same way -- and hold back roughly a third.
That keeps every declared pair intact and stops a public instantiation leaking
its private twin.

### 5. The one share-alike campaign

- [ ] **Decided:**

24 campaigns are `license_class: open`. One is `share_alike`:
`selectivity-saureus-topoisomerase-01`. The obligation attaches to that
campaign's own data files, not to the suite: a collection of separable campaigns
is not an adaptation of any one of them.

**Recommendation.** Ship it, under the notice `ATTRIBUTION-1.0.0.md` already
assembles from the campaign's own `license_notice`. It is separable by
construction -- an evaluator who cannot accept the terms drops that one campaign
and the other 24 are unaffected, which is why the licence class is per-campaign
and the gate enforces the notice.

### 6. How the image and the grader reach evaluators

- [ ] **Decided:**

The release ships no `gold/` and no `grading.json` by design, so someone has to
grade. And the image needs about 30 GB free to export, though its layers now
build and its closure reproduces set-identically from the lock.

**Recommendation.** Publish **both** the build inputs (`Dockerfile`,
`environment.lock.json`, `environment.explicit.txt`) **and** a registry digest: the
lock makes a rebuild reproducible, and a digest is what makes two scores
comparable without trusting whoever built it. Grade **centrally** -- evaluators
submit, you measure. It is the only arrangement in which a private set stays
private.

---

## C. Can wait for calibration

### 7. Lift the 25-campaign freeze for items 14 and 15?

- [ ] **Decided:**

**Recommendation: no, not yet.** The sweep may retire campaigns, and building two
more now risks building into a suite that is about to shrink. Note also that item
15 is real work rather than a copy: T-L4-2's second instantiation needs a
neighbour cap and a maximum component size declared with the measurement that
sets them, because four of its components hold more than ninety nodes and the
build refuses a hairball on purpose. Revisit both after calibration.

### 8. T-L5-2 (open item 6)

- [ ] **Decided:**

**Recommendation: do not build it.** All three honest options cost something --
HIV-1 RT drifts off theme, a re-scoped census ladder weakens the ladder, waiting
for a release that fills the bacterial variant values is indefinite. Record it as
deliberately absent with the counts that made it so. A gap with a reason is a
better artefact than a campaign with a compromise.

### 9. The construct-validity study (open item 12)

- [ ] **Decided:**

**Recommendation: after calibration, and scoped to three or four campaigns rather
than 25** -- the ones whose gold is most judgement-adjacent
(`kinetics-brenda-ec1_1-01`, `mibig-gene-function-evidence-01`,
`chemspace-mibig-4_0-01`). Inter-rater agreement on those buys the most per
expert-hour, and expert attention is the scarcest input in this project.

---

## Not a decision: one thing the project still owes

- [ ] The gate's **mechanical check count** in `docs/npbench-c.md` is maintained
  by hand. Every other structural figure in the documents is now counted from the
  artefacts and asserted by `tests/unit/test_readiness.py`, but deriving this one
  means running the gate and rebuilding gold, so it has no test behind it -- and
  it has lagged the code at 25, 26 and 27 checks already. Either accept it as
  prose or give the gate a cheap way to report its own check names without
  running them.
