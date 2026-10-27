# Task 02 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

## Tool choice

PROJECT.md §4 lists antiSMASH, GECCO and DeepBGC for BGC detection. GECCO was
chosen for the first detection task because it pins reproducibly from PyPI
wheels with no reference-database download and no conda channel, which makes
the environment a 60-line lock file instead of a multi-gigabyte image. antiSMASH
remains the right tool for a later task that needs its cluster typology; it
needs a provisioned machine with docker or conda.

## Determinism

Two independent GECCO 0.11.0 runs of `GCF_000016425.1` produced byte-identical
`clusters.tsv`. Gene calling (pyrodigal) and the CRF prediction are both
deterministic. Truth is therefore stable across re-runs.

## Scale

Measured on 21 of the 120 assemblies: 516 clusters, 24.6 per genome, which
matches the per-strain cluster counts reported for *Salinispora*. Extrapolated,
the full inventory is roughly 2,900 rows — past the point where an agent can
eyeball it, as PROJECT.md §8 asks for.

142 of those 516 clusters (27.5%) sit on a contig boundary. That is the
scaffold-assembly mess the task deliberately leaves in: an agent that ignores
contig edges over-reports complete clusters by more than a quarter.

Runtime: about 37 s per assembly on 4 vCPUs, so roughly 75 minutes for 120.

## Grader behaviour on real output

Each variant perturbs `bgc_inventory.csv` only, so a task score of 0.80 means
four of the five artifacts were untouched and passed.

| Perturbation | Artifact score | Passes | Format axis |
|---|---|---|---|
| none (truth vs itself) | 1.00 | yes | 1.00 |
| drop 10% of clusters | 0.90 | **yes** | 1.00 |
| shift every coordinate +50 kb | 0.06 | no | 1.00 |
| relabel every type as `Unknown` | 0.06 | no | 1.00 |
| rename `cluster_id` to `bgc_id` | 1.00 | yes | 0.80 |
| duplicate 20 rows | 0.00 | no | 1.00 |

Three of these are the graders working as intended. Coordinate drift and label
collapse both fail hard, which is what interval overlap F1 and macro F1 exist
for. A renamed column costs nothing on the artifact score and shows up only on
the format axis, which is the BixBench3 confound PROJECT.md §3 sets out to
remove. Duplicated rows fail the artifact outright.

**One finding for the calibration panel.** Dropping 10% of clusters scores 0.90
and passes. The aligned metrics only see the intersection, so the whole penalty
comes from identifier F1, and at 90% recall with perfect precision that is
0.947. Whether losing a tenth of a BGC inventory should count as reproducing it
is exactly the judgement PROJECT.md §7 says to put to domain experts rather than
inherit from BixBench3. If 0.80 survives that review, this artifact is more
forgiving of silent row loss than its other three metrics suggest.
