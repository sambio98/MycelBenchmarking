# Task 06 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

**Read the difficulty section before using this task.** Three of its seven
artifacts cannot be failed by any mistake modelled here, and for two of them
that is a structural property of the tools rather than something a better
column would fix.

## Rebuilding RODEO from its parts

RODEO is a web tool and is not installable offline, so the pipeline is
reconstructed: pyrodigal for gene calling, pyhmmer against a fixed 13-family
RiPP Pfam panel, a neighbourhood rule anchored on the lasso cyclase PF13471,
and a sequence similarity network over the precursor candidates.

The paper's 1,419 lasso BGCs are **not** the truth here. That set came from a
different genome collection, a different tool and a curation step, and is not
reproducible from any input we can stage. This task measures execution of the
prescribed pipeline, as `task.yaml`'s `truth_provenance` says.

Two things made it feasible, both checked rather than assumed:

- **GECCO bundles the HMMs.** Its install ships a 2,766-model subset of Pfam
  35.0 with a recorded md5, and it contains all three defining lasso families —
  PF13471 (cyclase), PF00733 (macrolactam synthetase), PF05402 (RRE) — plus ten
  more RiPP markers. The task needs no database download and fetches nothing at
  runtime. The solution finds the file through the `gecco` package rather than
  an absolute path.
- **Biopython does ~30,000 local BLOSUM62 scores per second**, so an all-pairs
  network over 798 precursors costs seconds and needs no BLAST. The network uses
  alignment **scores only, never the traceback**: a score is stable across
  Biopython versions in a way the chosen alignment path is not.

## The genome set is a rule, not a cap

Every NCBI RefSeq **reference genome** of a named *Streptomyces* species with a
complete, current assembly. There are 1,543 complete *Streptomyces* assemblies
but only 213 reference genomes, exactly one per species, so NCBI's own curation
does the sampling. That avoids both an arbitrary cutoff and the
over-representation of heavily sequenced species, which would have distorted
the similarity network and any species-breadth claim.

## Reproducibility

Two independent runs produce byte-identical output for all seven artifacts.

## Scale

| | |
|---|---|
| genomes | 213 (one per named species), 517 MB staged |
| called proteins | ~1.5 million, 4,892 to 10,637 per genome |
| panel HMM hits | 8,330 |
| lasso loci | 171 — `complete` 125, `core` 38, `minimal` 8 |
| precursor candidates | 798, lengths 29 to 120 aa |
| precursor families | 483, of which 424 singletons, largest 106 |
| artifacts / metric instances | 7 / 37 |
| runtime | ~8 s per 8.7 Mb genome, so ~30 min for 213 |

The 424 singleton families out of 798 precursors are not a failure of the
network: lasso precursors are hypervariable, which is exactly why they cannot be
found by similarity alone and need the cyclase-anchored neighbourhood rule.

## Measured difficulty

`scripts/difficulty_task06.py`, 16 variants spanning the whole chain. Same
cheap-proxy caveat as the other tasks: synthetic variants, rates in multiples
of 1/16, no substitute for model runs.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: E=1e-3 | 1.00 |
| one slip: keeps all reported hits | 1.00 |
| one slip: no 20 aa precursor floor | 1.00 |
| one slip: global alignment in the SSN | 1.00 |
| one slip: drops zero-hit families | 1.00 |
| one slip: no anchor merge | 0.86 |
| one slip: class ignores the RRE | 0.86 |
| one slip: families numbered by rank | 0.86 |
| one slip: no transitive closure | 0.86 |
| one slip: SSN normalised by length | 0.71 |
| one slip: trains per contig | 0.43 |
| one slip: 5 kb window | 0.43 |
| hasty agent (3 slips) | 0.43 |
| careless agent (5 slips) | 0.29 |

| artifact | pass rate | |
|---|---|---|
| `precursor_families.csv` | 0.50 | in band |
| `ripp_loci.csv` | 0.62 | in band |
| `precursor_candidates.csv` | 0.62 | in band |
| `genome_locus_counts.csv` | 0.75 | in band |
| `hmm_hits.csv` | 0.81 | **borderline** |
| `genome_proteins.csv` | 1.00 | **above band** |
| `domain_prevalence.csv` | 1.00 | **above band** |

### Across the five tasks

| | artifacts in band | 3 slips | 5 slips |
|---|---|---|---|
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| 05, shared joins, two inputs | 7 of 7 | 0.29 | 0.00 |
| **06, this task** | **4 of 7** | **0.43** | 0.29 |

Task 06 sits in the middle: better than task 03, weaker than 02, 04 and 05.
The dependency chain is there and the deep artifacts do discriminate — the
locus, precursor and family artifacts all land in band — but the two ends of
the chain do not.

### Why five of the sixteen variants cost nothing

These are structural, not oversights, and no added column fixes them:

- **The 20 aa precursor floor can never bite.** pyrodigal has its own minimum
  gene length, so the shortest precursor in the truth is 29 aa. Removing the
  floor changes nothing. The rule stays because the artifact needs a defined
  answer, not because it is hard.
- **The E-value and the inclusion filter are both no-ops.** HMMER's own
  inclusion thresholds dominate the search: the lowest kept bit score is 23.1,
  and loosening E to 1e-3 or keeping every reported hit adds almost nothing
  `included` would have rejected. Prescribing E=1e-5 is honest documentation of
  what was run, not a difficulty lever.
- **Global instead of local alignment barely moves the network** (0.89 on
  `precursor_families.csv`, still a pass) because the precursors are short and
  of similar length, so the two modes score them alike.
- **Dropping the zero-hit panel families is one row short of biting.** Three of
  the thirteen families — PF02052, PF08130, PF10439 — have no hits in these
  genomes, so omitting them gives identifier F1 of 0.87, just above the 0.80
  threshold. The rule is exercised but not punished.

### The two above-band artifacts, and why I did not contrive a fix

- **`genome_proteins.csv` (1.00).** It summarises the gene call, and both
  plausible gene-calling choices produce near-identical protein counts: even
  per-contig training leaves its CCC at 1.00. There is no realistic mistake in
  this pipeline that moves a 4,892–10,637 protein count enough to fail.
- **`domain_prevalence.csv` (1.00).** Thirteen rows of aggregate counts over
  the hits; nothing modelled here gets the per-family totals wrong.

I considered adding content-derived columns to both, as was done for task 05's
weakest artifact. I did not, because task 03 already answered what that is
worth: **a content-derived column only helps when some mistake actually gets
the grouping or naming wrong**, and here no modelled mistake does. Adding one
anyway would manufacture a number without adding information — the exact error
task 03's notes record.

The spec-aligned action is the one PROJECT.md §8 prescribes: build more
candidate artifacts than needed, measure, and drop the uninformative ones at
calibration. These two are the drop candidates, `genome_proteins.csv` first.

`hmm_hits.csv` at 0.81 is one variant above the band and this proxy cannot
resolve it: 1/16 is 0.0625, so 0.81 and 0.75 are adjacent values. Treat it as
undecided until real runs say otherwise, which is what `npbench calibrate` does
below four runs.
