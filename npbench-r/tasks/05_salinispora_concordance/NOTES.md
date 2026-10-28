# Task 05 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

## Why this task is shaped the way it is

Task 03 found that a task discriminates when its artifacts share upstream
dependencies rather than when it has many; task 04 confirmed it but ended up
bimodal, flooring at 0.00 on a single scope mistake. Task 05 keeps the shared
dependencies and spreads the weight across two independent inputs so that no
single mistake takes everything:

```
k-mer parse ──> genome_kmer_profiles ──┐
                                        ├──> genome_pair_distances ──┬─> species_concordance
GECCO clusters ──> bgc_families ──┬─────┘                            └─> distance_concordance_bins
                                  ├──> genome_repertoires
                                  └──> family_taxonomic_breadth
```

## The result the task is built around

Pearson 0.924 and Spearman 0.711 between genomic and chemistry distance over
7,140 pairs, with mean chemistry distance rising monotonically across deciles
of genomic distance, 0.253 to 0.746. Specialized metabolism tracks phylogeny in
this genus, which is what the source study argued.

## The sketch: three candidates, measured

The genomic distance needed a relatedness measure that is cheap, discriminating
within one genus, and **reproducible without agreeing on an implementation**.

| candidate | within-species | between-species | 120 assemblies | verdict |
|---|---|---|---|---|
| tetranucleotide composition | — | — | 73 s | rejected: all distances fell in 0.0004–0.004, too coarse to carry signal within a genus |
| MinHash, 2,000 elements, k=21 | 0.646 | 0.920 | 833 s | rejected: distances depend on the hash function, so an agent that picked a different one fails every pair for a reason that is not science |
| **canonical 21-mers kept when the 2-bit encoding is divisible by 1000** | **0.447** | **0.905** | **186 s** | chosen: nothing to agree on but the encoding and the modulus, 4.5x faster, and better separated |

Canonicalisation matters and is not cosmetic: contig orientation is arbitrary,
so without taking the smaller of a k-mer and its reverse complement two
assemblies of the same strain can share almost nothing.

## Reproducibility

Two independent runs produce byte-identical output for all seven artifacts.
The family rule is identical to task 02's, and both tasks independently find
**226 families over 2,933 assignable clusters** — a cross-check on both.

## Scale

| | |
|---|---|
| assemblies | 120 (1 complete, 119 scaffold), 191 MB staged |
| kept k-mers per assembly | 5,961 to 7,695 |
| clusters | 2,935; assignable 2,933; families 226 |
| genome pairs | 7,140 |
| species with two or more assemblies | 7 |
| artifacts / metric instances | 7 / 32 |
| runtime | GECCO ~75 min, sketches ~3 min, the rest seconds |

## Measured difficulty

`scripts/difficulty_task05.py`, 13 variants. Same cheap-proxy caveat as the
other tasks: synthetic variants, rates in multiples of 1/13, no substitute for
the model runs §8 asks for.

| variant | task score | artifacts damaged |
|---|---|---|
| reference (no mistakes) | 1.00 | — |
| one slip: k-mer multiset | 1.00 | none |
| one slip: thin clusters kept | 1.00 | none |
| one slip: pairs in both orders | 0.86 | 1 |
| one slip: equal-width bins | 0.86 | 1 |
| one slip: modulus 1024 | 0.71 | 4 |
| one slip: families numbered by rank | 0.57 | 3 |
| one slip: chemistry over clusters | 0.57 | 3 |
| one slip: forward k-mers only | 0.43 | 4 |
| one slip: species keeps the strain | 0.43 | 4 |
| one slip: no transitive closure | 0.14 | 5 |
| hasty agent (3 slips) | 0.29 | 6 |
| careless agent (5 slips) | 0.00 | all 7 |

All seven artifacts land in the 0.1–0.8 band:

| artifact | pass rate |
|---|---|
| `genome_pair_distances.csv` | 0.46 |
| `species_concordance.csv` | 0.46 |
| `distance_concordance_bins.csv` | 0.46 |
| `family_taxonomic_breadth.csv` | 0.62 |
| `genome_kmer_profiles.csv` | 0.69 |
| `bgc_families.csv` | 0.77 |
| `genome_repertoires.csv` | 0.77 |

### Across the four tasks

| | artifacts in band | 3 slips | 5 slips |
|---|---|---|---|
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| **05, shared joins over two inputs** | **7 of 7** | **0.29** | 0.00 |

Task 05 has the most usable profile: every artifact informative, and a
three-mistake agent at 0.29 — nearest the 0.3–0.5 target in §8 — without the
whole task collapsing the way task 04 does on one scope error. Only a
five-mistake agent floors, which is the right behaviour.

### The content-derived column: two data points, one conclusion

`genome_repertoires.csv` first measured 0.92, above the band, because it
reported counts only: family **naming** could be wrong without it noticing, and
the rank-naming variant scored 1.00 on it. Adding `largest_family_id` — the
genome's most populous family, ties to the smallest id — moved it to **0.77**,
in band. The rank-naming variant went from 1.00 to 0.00 on that artifact.

The same device did **nothing** for task 03's three aggregate artifacts, and
task 03's notes say so. The difference is now clear and worth stating as a
rule: **a content-derived column only helps when something actually gets the
grouping or the naming wrong.** Task 03 had no variant that moved a group's
smallest member, so the column had nothing to catch; task 05 had one. Adding
such a column on principle, without a mistake that probes it, buys nothing.

### Three rules that do not discriminate

Honest negatives, all measured:

- **The k-mer set rule is nearly a no-op.** Keeping duplicates scores 1.00
  everywhere except 0.97 on `n_kept_kmers`, which still passes, because both
  distances are computed over sets and the duplicates never reach them.
- **The thin-cluster filter is a no-op**, exactly as task 02's notes record:
  only 2 of 2,935 clusters carry fewer than 3 distinct domains. Same rule, same
  non-effect, now confirmed on a second task.
- **Emitting each pair in both orders** costs only `genome_pair_distances.csv`
  (0.67, still a pass). The duplicated rows are extra identifiers rather than
  wrong ones, and identifier F1 at 7,140 true positives against 7,140 false
  positives is 0.67 — above the 0.80 threshold only because the other three
  metrics on that artifact are unaffected. Worth knowing before reading a near
  pass there as a clean run.
