# Task 01 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

**Read the scope section before using this task.** It grades the GNPS-clustered
consensus spectra, not the 311 GB of raw mzXML, and that is a deliberate
reduction with a stated justification. It is the one decision here a reviewer
might reasonably overrule.

## Scope: the clustered spectra, not the raw files

MSV000078836 holds 1,125 mzXML files totalling 311 GB — fifteen times the 20 GB
ceiling in PROJECT.md §10 — so reprocessing the raw deposit is out of budget.
The same deposit carries the GNPS-clustered consensus spectra as a single 17 MB
`.mgf`, and PROJECT.md §9 names precisely that file as the reproducible
artifact. The task stages the `.mgf`.

What this costs: the agent does not perform the clustering step, so no artifact
here tests peak picking or consensus building. What it keeps: everything
downstream of clustering — filtering, pairwise spectral similarity, the
top-K cut, and the whole network topology — which is where the paper's
molecular networking claims actually live.

The paper's own networks, strain counts and prioritisation are **not** the truth
here, for the reason `task.yaml`'s `truth_provenance` states: they came from a
GNPS web job whose version is not recoverable. This task measures execution of
the prescribed pipeline against a committed script.

## Getting the file at all

MassIVE's own HTTPS endpoint is unreachable from a current certificate store —
its chain roots in an emSign CA present in neither the container trust store nor
`certifi`. The manifest therefore fetches through `massiveproxy.gnps2.org`,
which GNPS2 operates in front of the same deposit. The manifest pins the
accession, the byte count and the sha256, so the route the bytes travel does not
affect what the task grades.

## Scale

| | |
|---|---|
| spectra in the file | 100,263 |
| spectra with peaks | 4,081 — the other 96,182 carry none |
| retained after filtering | 4,081 |
| candidate edges | 15,427, of which 3,208 survive the top-10 cut |
| cosine range | 0.7002 to 0.9999, median 0.932 |
| components (≥2 nodes) | 325, sizes 2 to 115, median 3 |
| singletons | 2,256 of 4,081 |
| node degree | 0 to 10, mean 1.57 |
| precursor mass bins | 13, from 200 to 1500 Da |
| artifacts / metric instances | 6 / 27 |
| staged | 17 MB |
| runtime | 50 s end to end |

The 2,256 unnetworked spectra are the normal shape of a GNPS network at these
settings, not a failure. The degree ceiling of exactly 10 is the top-K cut
doing its job.

## Reproducibility

Two independent runs produce byte-identical output for all six artifacts.

## Measured difficulty

`scripts/difficulty_task01.py`, 17 variants. Same cheap-proxy caveat as the
other tasks: synthetic variants, rates in multiples of 1/17, no substitute for
the model runs PROJECT.md §8 asks for.

| variant | task score |
|---|---|
| reference (no mistakes) | 1.00 |
| one slip: skips the precursor window | 1.00 |
| one slip: plain cosine, not modified | 1.00 |
| one slip: 0.1 Da fragment tolerance | 0.80 |
| one slip: links at cosine 0.6 | 0.80 |
| one slip: 5 Da precursor window | 0.80 |
| one slip: singletons in components | 0.80 |
| one slip: bins by rounding | 0.80 |
| one slip: missing charge is 1+ | 0.80 |
| one slip: base peak is highest m/z | 0.80 |
| one slip: components numbered | 0.60 |
| one slip: keeps the empty spectra | 0.40 |
| one slip: no matched-peak minimum | 0.20 |
| one slip: no top-10 cut | 0.20 |
| one slip: top-10 of either end | 0.20 |
| hasty agent (3 slips) | 0.20 |
| careless agent (5 slips) | 0.00 |

| artifact | pass rate | |
|---|---|---|
| `components.csv` | 0.47 | in band |
| `precursor_mass_bins.csv` | 0.53 | in band |
| `network_nodes.csv` | 0.59 | in band |
| `edges.csv` | 0.71 | in band |
| `spectra.csv` | 0.76 | in band |
| `singletons.csv` | 0.88 | soft, out of the headline |

Every artifact that counts toward the headline lands in band. At three slips the
task scores 0.20, a little below the 0.3–0.5 the spec targets, so this is on the
hard side of the suite without bottoming out the way task 04 does.

The dependency chain is the reason: the retained spectrum set feeds the edges,
the edges feed the nodes, and the nodes feed the components, the singletons and
the mass bins. The three variants that cost the most — dropping the
matched-peak minimum, skipping the top-10 cut, applying it at either end — all
act on the edge list and are each felt by four artifacts.

### The measurement was wrong before it was right

The first sweep ran 15 variants and put `spectra.csv` at 0.87, above band. The
tempting read was that the artifact is too easy and needs a harder column. It
was not: **every variant in that set was drawn from the pipeline's parameters,
and `spectra.csv` is almost entirely a product of the parsing layer.** Only one
of the fifteen touched parsing at all.

Two realistic parsing misreadings, added as variants 16 and 17, fail it hard:

| parsing mistake | `spectra.csv` | which metric |
|---|---|---|
| missing `CHARGE` read as 1+ | 0.57 | macro F1 on `charge` falls to 0.568 |
| base peak read as highest m/z | 0.17 | Lin's CCC on `base_peak_mz` falls to 0.170 |

Neither is a typo; both are readings a competent agent could defend. 1,741 of
the 4,081 spectra carry no `CHARGE` line, and assuming singly-charged in
positive-mode ESI collapses two of the three classes macro F1 averages over.
"Base peak" means most intense, but it reads like highest.

That moved `spectra.csv` to 0.76 with **no change to the bundle**. The general
lesson is about measurement, not difficulty: a variant set built from a
pipeline's parameters measures the pipeline and is blind to the parsing in front
of it, so an artifact that lives on parsing looks easy when it is only
unmeasured. An above-band verdict is a claim about the variant set at least as
much as about the artifact.

I checked whether the other five sweeps have the same hole. They do not: 02
models zero-based coordinates, off-by-one lengths and strain-suffixed organism
strings; 03 models stated-versus-computed formulae, raw-name deduplication and
stereo stripping; 04 models module keying, ordering and evidence selection; 05
models k-mer canonicalisation and pair orientation; 06 models both gene-calling
modes and both hit-inclusion rules. Every one of those is a field-reading
mistake, and in task 06's case they sit on exactly the two layers its two
above-band artifacts live on — which is why its notes can say those artifacts
are structurally safe rather than merely untested.

Task 01 was the exception because it is the only task whose input is a raw
instrument format rather than a table or a JSON document, where the
field-reading decisions are obvious enough to model without thinking about
them. So the hole is real, it was worth finding, and it was local.

### Why `singletons.csv` is soft

It is the one artifact here that cannot discriminate, and that is structural.
Its three columns — `scan`, `precursor_mz`, `n_peaks` — are copies of
`spectra.csv` restricted to the degree-zero nodes, so once the node set is right
the columns follow automatically. Its only independent claim is the membership
of a 2,256-row set, and no scoring mistake reaches the interior of that set:
changing the cosine threshold or the top-K rule only flips the handful of nodes
sitting on the connectivity boundary, which is why it scores 0.96 to 0.99 under
variants that take `components.csv` down to 0.17.

It failed twice in seventeen variants, both times on the single mistake that
changes the node universe. Worse, while it counted toward the headline it handed
the 3-slip and 5-slip agents a free pass: marking it soft took the hasty agent
from 0.33 to 0.20 and sharpened the headline rather than merely tidying a band
number.

It stays in the bundle because the 2,256-of-4,081 split is a reported result and
the agent should have to produce it. It stays out of the headline because it
cannot tell a careful agent from a lucky one. I did not add a network-derived
column to rescue it — the near-miss cosine was the obvious candidate — because
task 03 already established what that is worth: a derived column only helps when
some modelled mistake gets the grouping wrong, and here the only mistake that
moves this artifact is one the other five artifacts already punish harder.

### Two prescribed filters that cannot bite

Measured, not assumed, and no added column fixes either:

- **The ±17 Da precursor-window removal is invisible to the metrics.** It is
  not a no-op on the data: it changes the peak count of 81 of the 4,081 spectra,
  dropping 83 peaks in total. But it pushes none below the six-peak floor, so
  identifier F1 stays at exactly 1.0000; it moves `n_peaks` so little that Lin's
  CCC is 1.0000 to four decimals; and it changes the base peak of **zero**
  spectra, because the base peak is essentially never within 17 Da of the
  precursor. Skipping it scores 1.00 on all six artifacts. The step stays
  because it is standard practice and the artifact needs a defined rule.
- **The six-peak minimum never removes anything.** The peak-count histogram
  jumps straight from 0 to 6 — nothing in this file has between 1 and 5 peaks —
  so the floor only ever re-rejects what the empty-spectrum check already
  rejected. Same reasoning: kept for a defined answer, not for difficulty.

Both belong to the family task 06 found with its 20 amino-acid precursor floor:
a rule the spec requires, that the tool or the data makes unreachable. Writing
it down is honest documentation of what was run, not a difficulty lever.

### Across the six tasks

| | artifacts in band | 3 slips | 5 slips |
|---|---|---|---|
| **01, this task** | **5 of 5 hard** | **0.20** | 0.00 |
| 02, shared joins | 6 of 7 | 0.14 | 0.14 |
| 03, independent artifacts | 3 of 8 | 0.62 | 0.62 |
| 04, shared joins by design | 8 of 8 | 0.00 | 0.00 |
| 05, shared joins, two inputs | 7 of 7 | 0.29 | 0.00 |
| 06, RiPP mining | 4 of 7 | 0.43 | 0.29 |

Task 01 behaves like 04 and 05 rather than like 03: the artifacts share
dependencies and every one of them is informative. Of the three well-behaved
tasks, 05 is still nearest the spec's 3-slip target, 01 is slightly harder, and
04 bottoms out.
