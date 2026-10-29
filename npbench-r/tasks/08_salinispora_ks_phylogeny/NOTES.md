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

Pending — `scripts/difficulty_task08.py` is written and runs the reference
variant to an exact reproduction of `truth/`, but the full 18-variant sweep had
not finished when these notes were written. The script measures all eight
artifacts for every variant over the full dataset, which only became affordable
after the dereplication rework above.

The variants it covers, by layer:

- **Gene call**: meta mode, per-contig training.
- **Domain search**: E = 1e-3, keeping all reported hits, envelope rather than
  aligned coordinates, no 150 aa floor, one domain per protein.
- **Leaf set**: no dereplication at all, grouping by length instead of by
  sequence.
- **Distance**: normalising by length, normalising by the larger self-score,
  leaving the distance unclamped.
- **Tree**: single linkage, complete linkage, splits numbered by merge order.
- **Taxonomy**: species keeping the strain suffix, breadth counted on genomes
  instead of species.

Run it with:

```bash
python scripts/difficulty_task08.py --task tasks/08_salinispora_ks_phylogeny \
  --truth tasks/08_salinispora_ks_phylogeny/truth --scratch <dir>
```

It caches extracted domain tables and raw alignment scores to `<dir>`, so a
second run is minutes.
