# Task 09 build notes

Measurements taken while building this task. They are not part of the contract;
they are the evidence behind the choices in `task.yaml` and the input to the
difficulty calibration in PROJECT.md §8.

**Read the scope section before using this task.** It stages 33.5 GB, which is
1.7x the 20 GB figure PROJECT.md §10 calls ideal, and that is the one decision
here most worth a second opinion.

## How this task became possible

PROJECT.md §4 category 5 — paired omics and metabologenomic link ranking — was
the last category with no task, and §11 calls it the differentiator worth
building. The README recorded it as unbuildable, with three dead ends:

- **The clustered spectra carry no per-strain provenance.** All 100,263 records
  in the staged `.mgf` say `FILENAME=specs_ms.pklbin`, so the file cannot say
  which strain a consensus spectrum came from.
- **MassIVE's own endpoints are unreachable from a current certificate store.**
  Its chain roots in an emSign CA present in neither the system trust store nor
  `certifi`. The GNPS2 proxy that works instead serves only exact paths and
  returns no directory listing, and the GNPS cluster-info tables — which would
  carry the per-file mapping — 404 at every documented name.
- **MSV000098188**, the other paired-omics candidate in §5, has peak files that
  404 at every layout tried, and its 101 "microscale" strains have zero overlap
  with the sequenced set.

All three are real and anyone rebuilding this will hit them. The conclusion
drawn from them was still wrong, because it was one source short: the **Paired
Omics Data Platform**, which §5 names for this very dataset, curates exactly
the missing mapping. Its project record for MSV000078836 holds **360
genome-to-file links** — 120 strains across three extraction methods — each
naming the deposit path of one strain's `.mzXML`. 111 of those strains already
had an NCBI assembly staged for task 02.

The generalisable lesson: *"unavailable" was a claim about the sources I had
tried, not about the data.* Every dead end above concerns file discovery on one
archive, while the mapping lived in a second archive built for precisely this
purpose.

## Scope, and what it costs

**One extraction method, every paired strain.** PoDP records three extraction
methods (ethyl acetate, butanol, methanol) for all 120 strains, all from one
`Screening` batch with identical sample preparation and instrumentation. The
task uses the ethyl acetate extracts only, which holds the analytical batch
constant — comparing metabolomes across extraction methods would confound every
occurrence call downstream. That rule, not an arbitrary cut, is what fixes the
strain count at 111.

**33.5 GB staged**, the largest fetch in the suite by a factor of sixty-seven.
The mzXML files run 113 to 790 MB, mean 300 MB. PROJECT.md §10 says to keep
every task well under 100 GB and ideally under 20 GB, so this satisfies the
limit and is **1.7× the ideal**.

Both of my earlier estimates were wrong and in the same direction, which is
worth recording as a method note: 17.8 GB came from extrapolating a single
160 MB file, and 27.9 GB from the first eight, whose mean was 252 MB. The true
mean is 300 MB and the distribution has a long right tail — the largest file is
790 MB, seven times the smallest. **Extrapolating a staging budget from a
sample of one, or of eight, is not estimating; the only honest number came
from fingerprinting all 111.**

I kept all 111 strains rather than cutting to fit, because "every paired strain
at one extraction method" is a rule and any smaller set needs an arbitrary one.
Staying under 20 GB would mean **66 strains**, dropping 45 of them. Metcalf
scoring degrades gracefully with fewer strains, so that reduction is available
to a reviewer who wants it; it is a curation decision rather than a technical
one, and it is the decision in this suite I would most like a second opinion
on.

**The molecular families are the depositors', not ours.** They come from the
GNPS clustering frozen in the `.mgf` rather than from re-clustering the raw
scans. §9 names network topology from a fixed `.mgf` as the reproducible
artifact, and reimplementing MS/MS clustering deterministically is a research
project rather than a task. What this task recovers is which *strain* each
consensus spectrum belongs to, which is the part the deposit does not record.

**One duplicate assembly dropped.** RefSeq holds two assemblies of
*Salinispora tropica* CNS416, both at scaffold level. The strain is the unit of
every occurrence and co-occurrence count, so one had to go or that strain would
be counted twice in the linkage matrix. The rule is recorded rather than
silently applied: better assembly level, then smallest accession.

## The link score rewards ubiquity

The score is the raw Metcalf weighted co-occurrence count from Doroghazi et al.
2014, which NPLinker implements as its `metcalf` function. It is **not**
NPLinker's *standardised* Metcalf, which subtracts the score's expected value
under a null model and divides by its standard deviation. An earlier draft of
this task called it standardised; a hand-check of the arithmetic caught that.

The raw score is what a graded artifact wants — a closed-form weighted sum with
nothing left to a judgement call. But it has a known bias, and it is visible
immediately: on a 26-strain smoke run the top-ranked links were a gene cluster
family present in all 26 strains against molecular families present in all 26,
scoring the maximum 260 whether or not they are related. Standardising against
a null model is the correction, and the obvious extension for a later revision.
Grading the raw score does not depend on the bias being absent; it depends on
the arithmetic being reproducible, which it is.

## Measured per-stage cost

Everything here was measured on this container rather than estimated.

| stage | cost |
|---|---|
| fetching one mzXML | 7 s cold, ~30 s sustained — the GNPS2 proxy throttles and returns 429 on a third rapid request, so the stager pauses 8 s between files and backs off to 240 s |
| reading and filtering one mzXML | ~2 s, 1,778 MS2 scans in the file probed |
| GECCO per assembly | 49 s on 2 of 4 vCPUs under contention, against task 02's measured 36 s on 4 idle ones |
| presence calls per strain | ~5 s, about 35,000 modified-cosine comparisons |
| consensus network | ~100 s, as task 01 |

The filter keeps almost everything: a typical file gives 1,627 usable MS2 scans
from 1,627, because after the precursor window and the top-50 cut nearly every
profile-mode scan still clears the six-peak floor. So the **six-peak floor is
very nearly a no-op on these spectra**, the same honest negative task 01
records for its consensus set. It stays because the artifact needs a defined
rule, not because it is hard.

## Staging had to stream

33.5 GB does not fit on this container beside everything else, and neither did
either of the figures it was estimated at. `scripts/stage_task09_spectra.py` therefore
fetches a batch, records bytes and sha256, runs the solution's own stage-2
extraction over it, and deletes the mzXML before fetching the next, keeping
peak disk at a few hundred megabytes. Deleting is safe because the manifest is
the record: `npbench fetch` restages any file and verifies its checksum.

That is a build-time convenience, not part of the contract. An agent attempting
the task stages all 111 at once.

One consequence is recorded rather than hidden: the difficulty sweep cannot
model the top-50 peak cut, because changing it means re-extracting from mzXML
this build no longer holds. Two variants for that lever were written and then
removed — a variant that silently measures nothing is worse than a documented
gap, which is the lesson task 07's provably-redundant column taught. An auditor
with the full staging can add it back.

## Reproducibility

Pending. Two independent runs must be compared before `task.yaml`'s
`determinism` field is changed from its current "pending verification".

## Scale

Pending the first full run.

## Measured difficulty

Pending. `scripts/difficulty_task09.py` is written: 16 variants across the
genomic, metabolomic and linking layers, ordered cheap-first, with the
expensive presence calls and the consensus network cached on disk per
configuration so a second run is minutes.
