# Task 08 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

## Why UPGMA and not a maximum-likelihood tree

PROJECT.md §4 category 8 asks for phylogenetics of biosynthetic enzymes. The
tools that field actually uses — RAxML, IQ-TREE — are not pip-installable in
this environment, and both are stochastic unless the seed is pinned, which §9
excludes in favour of artifacts reproducible by construction.

Average linkage over a distance matrix is UPGMA, a classical phylogenetic
method, and it is exact and deterministic given the matrix. The cost is real
and is stated in `task.yaml`: UPGMA assumes a molecular clock and will not
recover a topology an ML method would. So what this task grades is the domain
extraction, the dereplication, the distance definition and the tree's derived
structure — not the quality of phylogenetic inference.

The tree is graded by its **splits**, never by a Newick string. A Newick string
varies with child order and branch-length formatting; the set of splits is the
tree's topological content and is what Robinson-Foulds compares. A split is
named by the lexicographically smallest leaf below it, so no identifier depends
on an ordering the agent chose.

That name alone is not a key — a split and its parent can share a smallest leaf
— so the row key is `(split_id, n_leaves)`. That pair is provably unique: the
splits containing a given leaf are its ancestors, which form a chain under
inclusion, so two splits with the same smallest leaf and the same size are the
same split.

## The dereplication step is the whole reason the task works

The first version of this task built the tree over all 2,993 domains. The
numbers that forced a rebuild:

| | over all domains | over distinct sequences |
|---|---|---|
| leaves | 2,993 | 922 |
| pairwise alignments | 4,477,528 | 424,581 |
| alignment runtime | ~39 min | ~4 min |
| leaves at distance 0 from their nearest neighbour | **2,552** | 74 |
| median nearest-neighbour distance | **0.0000** | 0.0040 |
| merges at height exactly 0 | hundreds | 37 |

60 of the 120 assemblies are *S. arenicola*, so the same ketosynthase recurs
verbatim across strains: one sequence appears 53 times. With 85% of leaves tied
at distance zero, the UPGMA topology below those ties is decided by the input
order rather than by the sequences. It was still reproducible — two runs agreed
byte for byte — and it was still **arbitrary**, which is worse than slow,
because an agent that ordered its leaves differently would have been marked
wrong for an equally valid tree.

Collapsing identical sequences first is also simply what a phylogeneticist
does. It made the task ten times cheaper, made the tree meaningful, and made
the difficulty sweep affordable across every layer instead of leaving the three
tree artifacts unmeasured.

`domain_length_bins.csv` was dropped in the same rebuild and `sequence_groups.csv`
took its place. Six rows of binned lengths were a coarse artifact with little
to get wrong; the 922-row group table is where the collapse above is graded,
and it is a prerequisite for the tree being interpretable at all.

## Scale

| | |
|---|---|
| assemblies | 120, 182 MB staged |
| called proteins | ~590,000, 4,682 to 5,398 per assembly |
| PF00109 domains kept | 2,993, after the 150 aa floor |
| domains per assembly | 8 to 46 |
| proteins carrying ≥2 domains | 342, up to 6 domains in one protein |
| distinct domain sequences | 922, of which 491 appear once and one appears 53 times |
| domain length | 150 to 287 aa, median 248 |
| species | 9 (*S. arenicola* 60 genomes, *S. pacifica* 30, *S. tropica* 13) |
| tree splits | 921 |
| artifacts / metric instances | 8 / 36 |
| runtime | ~9 min gene calling and domain search, ~4 min alignment |

## What the tree says

Two numbers worth reading together:

- **Only 18 of the 922 sequence groups span more than one species.** An
  *identical* ketosynthase is almost always confined to one species.
- **436 of the 921 splits are multi-species, and of the 217 splits with ten or
  more leaves, 162 are.** Clades cross the species boundary even though exact
  sequences do not.

That is the shape the vertical-inheritance literature describes, and it is why
`split_taxonomy.csv` is the deepest artifact: the interesting claim is not how
many domains there are but whether a ketosynthase lineage predates the species
carrying it.

## Reproducibility

Two independent runs produce byte-identical output for all eight artifacts.

## Measured difficulty

`scripts/difficulty_task08.py`, 20 variants, all eight artifacts measured for
every variant over the full dataset. That is only affordable because the tree
is built over the 922 distinct sequences rather than all 2,993 domains. Same
cheap-proxy caveat as the other tasks.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: no clamp on the distance | 1.00 |
| one slip: complete linkage | 1.00 |
| one slip: E=1e-3 | 1.00 |
| one slip: keeps all reported hits | 1.00 |
| one slip: breadth counts genomes | 0.88 |
| one slip: single linkage | 0.75 |
| one slip: splits numbered | 0.75 |
| one slip: meta-mode gene calling | 0.75 |
| one slip: normalise by length | 0.62 |
| one slip: normalise by max self-score | 0.62 |
| one slip: species keeps the strain | 0.50 |
| one slip: groups by length, not sequence | 0.50 |
| one slip: no dereplication at all | 0.50 |
| one slip: envelope coordinates | 0.25 |
| one slip: trains per contig | 0.00 |
| one slip: no 150 aa floor | 0.00 |
| one slip: one domain per protein | 0.00 |
| hasty agent (3 slips) | 0.75 |
| careless agent (5 slips) | 0.00 |

| artifact | pass rate | |
|---|---|---|
| `split_taxonomy.csv` | 0.30 | in band |
| `tree_splits.csv` | 0.40 | in band |
| `group_nearest_neighbour.csv` | 0.55 | in band |
| `sequence_groups.csv` | 0.60 | in band |
| `ks_domains.csv` | 0.70 | in band |
| `multidomain_proteins.csv` | 0.70 | in band |
| `genome_ks_counts.csv` | 0.75 | in band |
| `species_domain_counts.csv` | 0.75 | in band |

**All eight artifacts land in band**, 0.30 to 0.75, and the deepest artifact is
the most discriminating — which is the shape the dependency argument predicts
and the one tasks 03 and 06 did not achieve. The task also has floors: three
single mistakes take it to 0.00, all of them in the domain layer.

What the sweep establishes:

- **The dereplication step is load-bearing, not an optimisation.** Grouping by
  length instead of exact sequence, or not grouping at all, both score 0.50 and
  drive `sequence_groups.csv` and `split_taxonomy.csv` to 0.00. Building the
  tree over all 2,993 domains — the thing the first version of this task did —
  is now a modelled mistake that fails half the task.
- **The domain layer is where the floors are.** Training the gene finder per
  contig, dropping the 150 aa floor, or keeping one domain per protein each
  score 0.00. The last is the mistake the prompt warns about in its own words,
  and it is the most expensive single slip here: a modular PKS is one
  polypeptide carrying up to six ketosynthases, so one-per-protein discards
  most of the data.
- **The distance normaliser matters more than the linkage.** Normalising by
  length fails three artifacts outright (0.00, 0.03, 0.04); normalising by the
  larger self-score lands at 0.78–0.79, just under threshold. Complete linkage
  scores 0.83 on both tree artifacts and passes.
- **The 3-slip composite scores 0.75, above eight of the single slips.** Same
  non-monotonicity task 07 showed, and for the same reason: its three mistakes
  are all in the tree and taxonomy layers, so they fail two artifacts between
  them. The composite columns count artifacts clearing 0.80, not severity.

### Four variants that cost nothing

- **The clamp cannot bite.** Proven below.
- **E = 1e-3 instead of 1e-5 changes nothing at all**, and **keeping every
  reported hit rather than only the included ones changes almost nothing**
  (0.98 on two artifacts, passing everywhere). This replicates task 06's
  finding exactly, on different genomes and a different Pfam model: HMMER's own
  inclusion thresholds dominate the E-value, so prescribing one is honest
  documentation of what was run rather than a difficulty lever. Two independent
  tasks now agree on that, which is worth more than either measurement alone.
- **Complete linkage barely differs from average** over these distances, at
  0.83 on both tree artifacts.

### The clamp cannot bite, and that is provable

`no_clamp` scores 1.00 on all eight artifacts, and this is not a sampling
accident. The distance is `1 - score(a,b) / min(self(a), self(b))`, and for
local alignment:

- `score(a,b) <= min(self(a), self(b))`, because a local alignment of a against
  b cannot score higher than the best local alignment of the shorter sequence
  against itself, so the distance is never negative;
- `score(a,b) >= 0`, because local alignment permits the empty alignment, so the
  distance never exceeds 1.

Measured over all 424,581 pairs: zero pairs exceed the floor, zero have a
negative raw score, and the largest ratio seen is exactly 1.000000 — reached by
the identical and substring pairs that give distance 0. The clamp is defensive,
not a difficulty lever.

Unlike task 07's `n_compounds_matched`, which was a provably redundant *column*
and was removed, this is a line of method guidance that asks nothing extra of
the agent. It stays, documented as a measured negative rather than counted as
difficulty.

### Across the eight tasks

| | artifacts in band | 3 slips | 5 slips |
|---|---|---|---|
| 01, molecular networking | 5 of 5 hard | 0.20 | 0.00 |
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| 05, shared joins, two inputs | 7 of 7 | 0.29 | 0.00 |
| 06, RiPP mining | 4 of 7 | 0.43 | 0.29 |
| 07, dereplication | 6 of 8 | 0.12 | 0.25 |
| **08, this task** | **8 of 8** | **0.75** | 0.00 |

Task 08 joins 04 and 05 as a task where every artifact is informative. Read the
last two columns with the caution the README states: they count artifacts
clearing 0.80, so they are not monotone in the number of mistakes and not
comparable across tasks.
