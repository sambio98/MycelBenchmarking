# Task 03 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

**Read the difficulty section before using this task.** By the proxy available
here it is materially weaker than task 02, and the reason is structural rather
than fixable by adding columns.

## Tool choice

RDKit ships manylinux wheels, so this is the one PROJECT.md §4 category that
runs end to end with no conda channel, no compiler and no reference-database
download. The whole reference solution takes 14 s.

That makes it deliberately light on tool weight. §8 lists installation pain as
a confound rather than a difficulty lever, so pairing a 14-second task with
task 02's 75-minute GECCO run separates analytical capability from
infrastructure luck. It is not an attempt to make the task easy.

## Reproducibility

Two independent runs of the reference solution produce byte-identical output
for all eight artifacts. RDKit canonical SMILES, Bemis-Murcko scaffolds, Morgan
fingerprints and the single-linkage assignment are all deterministic, and every
artifact is sorted on its row key before writing.

## Scale and mess

| | |
|---|---|
| MIBiG entries | 3,013 |
| staged data | 936 KB |
| compound entries with a usable structure | 4,397 |
| distinct canonical structures | 3,923 |
| distinct Bemis-Murcko scaffolds | 2,251, of which 219 are empty (acyclic) |
| formulas shared by 2+ distinct structures | 533, covering 1,270 structures |
| structures collapsing on stereochemistry removal | 345, onto 167 flat forms |
| compounds whose MIBiG `formula` contradicts the structure | 56 |
| duplicate `(accession, name)` pairs | 6: four identical, two differing only in case |
| similarity families | 1,480 |
| artifacts / metric instances | 8 / 31 |

All of that mess is MIBiG's own. PROJECT.md §8 names isomers sharing a formula
and stereochemistry collapsed by canonical SMILES as things to leave in; both
are graded here directly rather than merely tolerated.

MIBiG's `mass` field disagrees with the structure for several hundred entries
and its `formula` field for 56 of the retained compounds, which is why every
descriptor is specified as computed from the SMILES.

## Measured difficulty

`scripts/difficulty_task03.py`, 13 variants. Same cheap-proxy caveat as task
02: the variants are synthetic rather than drawn from observed model behaviour,
with 13 of them a pass rate only lands on multiples of 1/13, and this does not
replace the model runs §8 asks for.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: drops acyclic compounds | 1.00 |
| one slip: class counted per BGC | 1.00 |
| one slip: trusts MIBiG's formula field | 0.88 |
| one slip: no dedupe at all | 0.88 |
| one slip: dedupes on the raw name | 0.88 |
| one slip: strips stereo characters textually | 0.88 |
| one slip: no transitive closure | 0.88 |
| one slip: families numbered by rank | 0.88 |
| one slip: Morgan radius 3 | 0.88 |
| one slip: keyed on entries, not structures | 0.75 |
| hasty agent (3 slips) | 0.62 |
| careless agent (5 slips) | 0.62 |

| artifact | pass rate | |
|---|---|---|
| `compounds.csv` | 0.69 | in band |
| `formula_disagreements.csv` | 0.77 | in band |
| `similarity_families.csv` | 0.62 | in band |
| `scaffolds.csv` | 0.92 | **above band** |
| `stereo_collapse_groups.csv` | 0.85 | **above band** |
| `scaffold_counts.csv` | 1.00 | **above band** |
| `formula_collisions.csv` | 1.00 | **above band** |
| `class_scaffold_summary.csv` | 1.00 | **above band** |

### Why this task discriminates less than task 02

Task 02's mistakes propagate: getting the species rule wrong damages two
artifacts, getting the contig-edge rule wrong damages two others, because its
artifacts share upstream joins. A three-mistake agent there scores 0.14.

Task 03's artifacts are largely independent derivations of a single upstream
object, the set of distinct structures. A mistake lands on one artifact and
stops, so almost every single slip costs exactly 0.12 — one artifact of eight —
and a five-mistake agent still scores 0.62.

The lesson generalises and is worth carrying into tasks 04 onward: **a task
discriminates when its artifacts share upstream dependencies, not when it has
many artifacts.** Eight independent artifacts measure roughly as much as one.

### An attempted fix that did not work

`scaffold_counts.csv`, `formula_collisions.csv` and `class_scaffold_summary.csv`
originally reported counts only, and passed every variant. Each now also
carries `representative_smiles`, the lexicographically smallest canonical SMILES
among its members — content-derived, like `family_id`, so it can only be right
if the membership is.

It did not move the proxy at all: those three still pass 13 of 13. None of the
variants gets group *membership* wrong in a way that changes a group's minimum
member, so the new column has nothing to catch. It is kept because it closes a
real hole — counts can be right while membership is wrong — but it should not
be described as having improved discrimination, because it did not.

### What to do with this at calibration

The five above-band artifacts are the first candidates to drop, in the order
given by their pass rate. §8 prescribes exactly this: build more candidate
artifacts than needed, measure, and keep the informative ones. That decision
needs real model runs, not this proxy — a real agent can fail
`formula_collisions.csv` by running out of budget or mis-parsing the MIBiG JSON,
failure modes no variant here models.

If the whole task still fails to discriminate after real runs, the honest
options are to fold its three in-band artifacts into another task or to drop it
from the headline set and report it separately.
