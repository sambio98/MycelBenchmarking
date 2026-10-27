# Task 02 build notes

Measurements taken while building this task, on the full 120-assembly set.
They are not part of the contract; they are the evidence behind the choices in
`task.yaml` and the input to the difficulty calibration in PROJECT.md §8.

## Tool choice

PROJECT.md §4 lists antiSMASH, GECCO and DeepBGC for BGC detection. GECCO was
chosen for the first detection task because it pins reproducibly from PyPI
wheels with no reference-database download and no conda channel, which makes
the environment a 28-line lock file instead of a multi-gigabyte image.
antiSMASH remains the right tool for a later task that needs its cluster
typology; it needs a provisioned machine with docker or conda.

## Reproducibility

- Two independent GECCO 0.11.0 runs of `GCF_000016425.1` produced byte-identical
  `clusters.tsv`.
- Two independent runs of the whole truth pipeline produced byte-identical
  output for all seven artifacts. Gene calling (pyrodigal), the CRF prediction
  and the family assignment are all deterministic.
- The family assignment uses an inverted domain index and union-find. On the
  real 1,791-cluster intermediate it produced exactly the same 180 families as
  explicit all-pairs single linkage, with identical naming;
  `tests/test_task02_solution.py` repeats that check against brute force on
  randomised inputs.

## Scale

| | |
|---|---|
| assemblies | 120 (1 complete, 119 scaffold) |
| staged data | 191 MB |
| clusters | 2,935 (24.5 per genome) |
| clusters on a contig boundary | 929 (31.7%) |
| cluster types | 8, including `Unknown` and hybrids such as `NRP;Polyketide` |
| gene cluster families | 226 over 2,933 assignable clusters |
| graded artifacts / metric instances | 7 / 30 |
| GECCO runtime | ~36 s per assembly on 4 vCPUs, ~75 min for 120 |

The 31.7% contig-boundary figure is the scaffold-assembly mess the task
deliberately leaves in: an agent that ignores contig edges over-reports
complete clusters by nearly a third.

## Measured difficulty

`scripts/measure_difficulty.py` grades 15 variants, each a mistake a competent
agent actually makes, built from the same GECCO output the reference solution
reads. This is a **cheap proxy, not the calibration**. PROJECT.md §8 asks for
two runs of each of two or three frontier models, which this does not replace:
the variants are synthetic, they are not drawn from observed model behaviour,
and with 15 of them a pass rate can only land on multiples of 1/15.

What it is good for is finding artifacts that nothing can fail. It found two.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: length off by one | 1.00 |
| one slip: coordinates 0-based | 1.00 |
| one slip: thin clusters kept | 1.00 |
| one slip: thin clusters dropped everywhere | 1.00 |
| one slip: domains used as a multiset | 0.86 |
| one slip: domain count counts duplicates | 0.86 |
| one slip: families numbered by rank | 0.71 |
| one slip: no transitive closure | 0.71 |
| one slip: edge rule too strict | 0.71 |
| one slip: species keeps the strain | 0.57 |
| one slip: hybrid types split into rows | 0.57 |
| one slip: drops `Unknown`-type clusters | 0.43 |
| hasty agent (3 slips) | 0.14 |
| careless agent (5 slips) | 0.14 |

Artifact pass rate across the 15 variants, against the 0.1–0.8 band in §8:

| artifact | pass rate | |
|---|---|---|
| `bgc_inventory.csv` | 0.73 | in band |
| `contig_edge_clusters.csv` | 0.87 | **above band** |
| `clusters_per_genome.csv` | 0.67 | in band |
| `type_counts.csv` | 0.73 | in band |
| `species_type_counts.csv` | 0.80 | in band |
| `gene_cluster_families.csv` | 0.60 | in band |
| `gcf_summary.csv` | 0.60 | in band |

### What this says, and what it does not

A single mistake costs between 0 and 0.57 of the task score, and three
simultaneous mistakes cost 0.86. A model that executes the method cleanly
scores near 1.00, so the top-model score depends entirely on how many of these
slips a real model makes — which is the thing only a real run can tell us. The
single-slip band of 0.43–1.00 is consistent with the 0.3–0.5 target for a model
that makes two or three, but that is an inference, not a measurement.

### Three things that do not discriminate, and why

- **A 1 bp coordinate convention is invisible**, whether as an off-by-one length
  or a 0-based re-basing. Interval overlap F1 still matches at 0.9999 reciprocal
  overlap, and CCC on a 20–100 kb spread cannot see it. This is correct
  behaviour, not a gap: a 1 bp difference is not a failure to reproduce.
- **The `MIN_DOMAINS` filter is nearly a no-op.** Only 2 of 2,935 clusters carry
  fewer than 3 distinct domains, so applying it wrongly, or everywhere, changes
  almost nothing. It stays in the method because it is needed for a defined
  partition, but it is not a difficulty lever and should not be described as one.
- **`contig_edge_clusters.csv` sits above the band at 0.87.** Only two of the
  fifteen variants touch it, because its row set follows from a single rule.
  That does not make it a free mark for a real agent — it needs contig lengths
  from the FASTA, which the GECCO table does not contain, and its
  `contig_length_bp`, `accession` and `sequence_id` columns are all graded — but
  this proxy cannot distinguish it, and the real calibration should be expected
  to either move it into band or justify keeping it.

### An earlier finding, now superseded

On an earlier 53-genome snapshot, dropping 10% of clusters scored 0.90 and
passed, because the aligned metrics only see the intersection and the whole
penalty fell on identifier F1. Grading every required column (30 metric
instances rather than 11) makes a 10% row loss considerably more expensive,
since each lost row also perturbs the per-genome, per-type and per-family
aggregates. The underlying question still belongs to the calibration panel:
whether losing a tenth of a BGC inventory counts as reproducing it is the
judgement PROJECT.md §7 says to put to domain experts rather than inherit from
BixBench3.
