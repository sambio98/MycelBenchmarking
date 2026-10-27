# Task 04 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

## Why this task is shaped the way it is

Task 03's notes record a finding: **a task discriminates when its artifacts
share upstream dependencies, not when it has many artifacts.** Task 03 has
eight artifacts that are largely independent derivations of one object, and
almost every single mistake costs exactly one of them, so a five-mistake agent
still scored 0.62.

Task 04 was built to test that claim. Its dependency graph is deliberate:

```
modules.csv ──┬─> a_domain_substrates.csv ──┬─> substrate_counts.csv
              │                             ├─> substrate_evidence_summary.csv
              │                             ├─> nonproteinogenic_substrates.csv
              │                             └─> class_substrate_matrix.csv
              └─> module_architecture.csv ────> architecture_counts.csv
```

Nothing stands alone. A mistake in the module parse reaches all eight.

## Reproducibility

Two independent runs produce byte-identical output for all eight artifacts.

## Scale and mess

| | |
|---|---|
| BGCs with curated modules | 667 of 3,013 |
| modules | 4,034, of which **2,089 are `active: false`** |
| A-domain substrate assignments | 3,182 over 2,879 A-domains |
| substrates per A-domain | 0 to 9 |
| A-domains with no evidence | 295 |
| A-domain locations that are the `-1` sentinel | 2,852 of 2,901 |
| module names repeating within a BGC | 12 |
| substrate names spelled two ways (case only) | 12, so 281 names are 269 substrates |
| substrate names mapping to >1 structure | 9, incl. `alanine`, `valine` (D/L forms) |
| distinct architectures | 131 over 667 BGCs, sized 1 to 88 |
| artifacts / metric instances | 8 / 39 |
| runtime | 3 s |

## Measured difficulty

`scripts/difficulty_task04.py`, 13 variants. Same cheap-proxy caveat as the
other two: synthetic variants, multiples of 1/13, no substitute for model runs.

| variant | task score | artifacts damaged |
|---|---|---|
| reference (no mistakes) | 1.00 | — |
| one slip: `proteinogenic` from first entry | 1.00 | none |
| one slip: counts rows, not modules | 1.00 | none |
| one slip: first evidence, not strongest | 0.88 | 2 |
| one slip: class counted once | 0.88 | 1 |
| one slip: unsorted domain composition | 0.88 | 1 |
| one slip: keys modules on name | 0.75 | 2 |
| one slip: sorts the module list | 0.62 | 4 |
| one slip: raw substrate names | 0.62 | 3 |
| one slip: filters on the `-1` location | 0.38 | 5 |
| one slip: drops inactive modules | **0.00** | **all 8** |
| hasty agent (3 slips) | **0.00** | all 8 |
| careless agent (5 slips) | **0.00** | all 8 |

Every artifact lands in the 0.1–0.8 band:

| artifact | pass rate |
|---|---|
| `a_domain_substrates.csv` | 0.54 |
| `class_substrate_matrix.csv` | 0.54 |
| `modules.csv` | 0.62 |
| `substrate_counts.csv` | 0.62 |
| `nonproteinogenic_substrates.csv` | 0.62 |
| `module_architecture.csv` | 0.62 |
| `substrate_evidence_summary.csv` | 0.69 |
| `architecture_counts.csv` | 0.69 |

### The design claim, tested

| | artifacts in band | hasty agent (3 slips) |
|---|---|---|
| task 02 (shared joins) | 6 of 7 | 0.14 |
| task 03 (independent artifacts) | 3 of 8 | 0.62 |
| task 04 (shared joins by design) | **8 of 8** | **0.00** |

The claim holds. It is the strongest design lever found so far and should drive
tasks 05 and 06.

### The opposite risk, which this task now has

Three variants bottom out at exactly 0.00, and the single `drops inactive
modules` slip is enough on its own. That mistake removes 2,089 of 4,034 rows at
the root of the graph, so everything downstream fails at once. A task that
floors at 0.00 loses information just as surely as one that ceilings at 1.00 —
§8 says so about artifacts and it applies to task scores too.

So the honest reading is that task 04 is **bimodal**: get the scope right and
most artifacts follow; get it wrong and nothing does. Whether that is
acceptable depends on how often real models make a scope error of that size,
which is a question for the real calibration and not for this proxy. If it
turns out that frontier models commonly drop the inactive modules, the lever is
too large and the fix is to make `active` a reported column on a smaller set
rather than a filter over the whole task.

### Two prescribed rules that do not discriminate

Honest negatives, both measured:

- **`proteinogenic` by majority vote.** Taking the first entry's flag instead
  scores 1.00 everywhere. MIBiG is self-consistent enough on this flag that the
  majority rule almost never changes an answer. The rule stays because it makes
  the output defined, but it is not a difficulty lever.
- **`n_modules` as distinct modules rather than rows.** Counting rows instead
  also scores 1.00, because an A-domain's substrates carry distinct names, so
  the two counts coincide almost everywhere.

Neither should be described as contributing difficulty. Both are in the method
because the artifact needs a defined answer, not because they are hard.

### One design choice, flagged as such

The A-domain evidence strength ranking in `solution/make_truth.py` is this
task's rule, not a standard. It puts in vitro assays above in vivo work, above
product observation, above structural inference, above homology and sequence
prediction. CLAUDE.md says a parameter a source paper did not state is a
question rather than a default; this is not from a source paper, but it is the
same kind of judgement and is recorded here for review.

An earlier version ranked only six of the ten methods the release uses, which
let `Homology` outrank a direct enzyme assay through the unranked fallback. The
ranking is now total over the data. Changing it changes
`substrate_evidence_summary.csv` and one column of
`a_domain_substrates.csv`, and would require regenerating truth.
