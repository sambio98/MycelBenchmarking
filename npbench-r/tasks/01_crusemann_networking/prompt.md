# Molecular networking of the Crüsemann 146-strain marine actinomycete MS/MS survey

## General instructions

You are working in a container with 4 vCPUs and 8 GB of RAM. Budget an hour;
the analysis itself takes under a minute once it is written correctly.

- `data/` is read-only. Do not modify or move anything in it.
- Write every deliverable to `outputs/`. Nothing outside `outputs/` is collected.
- Direct network access is limited to an allowlist of package registries and OS
  mirrors, which is enough to install the tool named below. Any other resource
  must go through the web access request tool, where it is adjudicated. The
  GNPS hosts are denied: the staged spectrum file is the only input this task
  needs.
- Also write `outputs/METHODS.md`: a short account of what you ran, the exact
  versions, and anything you had to decide for yourself.

## Research objective

Molecular networking groups MS/MS spectra by fragmentation similarity, so that
related metabolites fall into connected families even when none of them is in a
spectral library. Given the clustered consensus spectra deposited with a survey
of 146 marine actinomycete strains, build the classic GNPS molecular network and
characterise its topology: which spectra survive filtering, which pairs are
similar enough to link, what families the links produce, and which spectra stay
unconnected.

**One property of the input dominates this task.** The file contains 100,263
spectra, and **about 96% of them carry no fragment peaks at all** — it is mostly
headers. Spectrum 2, for instance, has a precursor mass and no peak list.
Roughly four thousand spectra are usable. An agent that feeds the file straight
into a similarity calculation will network noise, or will attempt five billion
comparisons and run out of time.

A second property matters for reading the results rather than producing them:
the survey is described as 146 strains in its title, and the deposit and
downstream reprocessings have been reported with different counts. The clustered
file carries no per-strain attribution at all, so nothing here is reported per
strain.

### Staged data

| Path | Contents |
|---|---|
| `data/MSV000078836_specs_ms.mgf` | the GNPS-clustered consensus MS/MS spectra, 17 MB |
| `data.manifest.json` | the provenance record for the staged file |

Each spectrum carries `PEPMASS`, `CHARGE`, `MSLEVEL`, `ACTIVATION`,
`INSTRUMENT`, `TITLE` and `SCANS`. `SCANS` is the spectrum's identifier and is
unique.

## Method guidance

Use **matchms 0.33.1** and **networkx 3.7**. Pin matchms in particular: the
similarity class used here was called `ModifiedCosine` in earlier releases and
is `ModifiedCosineGreedy` in 0.33.1, and cosine scores are
implementation-dependent.

Every rule below is part of the method; none of it is open to choice.

**Filtering**, in this order:

1. Load with `matchms.importing.load_from_mgf`.
2. Discard any spectrum with no peaks.
3. `matchms.filtering.normalize_intensities`.
4. `matchms.filtering.remove_peaks_around_precursor_mz` with
   `mz_tolerance=17`. This removes fragment peaks within 17 Da of the
   precursor. It affects only a small number of spectra, but it does change
   their peak counts and base peaks, so do not skip it.
5. Keep spectra with **at least 6 peaks** remaining.

Sort the retained spectra by `SCANS` ascending and keep that order.

**Similarity.** `ModifiedCosineGreedy(tolerance=0.5)`. Compare a pair only when
the absolute difference of their precursor m/z is **at most 2.0 Da**; the rest
cannot form an edge and do not need scoring.

**Candidate edges.** A pair becomes a candidate edge when its score is
**0.70 or greater** *and* it has **at least 6 matched peaks**.

**The top-K cut.** Rank each spectrum's candidate edges by descending cosine,
breaking ties by ascending partner scan number. An edge survives only if it is
within the **top 10** of **both** of its endpoints. This is GNPS's rule, and
the tie-break matters because cosines rounded to 4 decimal places tie often.

**Components** are the connected components of the surviving graph, over all
retained spectra including unconnected ones. A component's `component_id` is
the **smallest scan number** it contains. Do not number components.

**Mass bins** are 100 Da wide, labelled by their lower edge: a spectrum with
precursor m/z 566.89 falls in bin 500. Report only bins containing at least one
retained spectrum.

Round every m/z, cosine, delta and mean to **4 decimal places**.

## Required outputs

Six CSV files in `outputs/`, each with a header row and the exact column names
below. Column order is not graded; column names are. Row order is not graded.

### `spectra.csv`
One row per **retained** spectrum.

`scan`, `precursor_mz`, `charge`, `n_peaks`, `base_peak_mz`

`charge` is the integer from `CHARGE`, which is `0` where the file does not
state one. `n_peaks` and `base_peak_mz` describe the spectrum **after**
filtering; `base_peak_mz` is the m/z of its most intense surviving peak.

### `edges.csv`
One row per surviving edge, with `scan_a` < `scan_b`.

`scan_a`, `scan_b`, `cosine`, `n_matched_peaks`, `precursor_delta`

`precursor_delta` is the absolute difference in precursor m/z.

### `network_nodes.csv`
One row per retained spectrum, including unconnected ones.

`scan`, `degree`, `component_id`, `component_size`

### `components.csv`
One row per component containing **two or more** nodes.

`component_id`, `n_nodes`, `n_edges`, `min_precursor_mz`, `max_precursor_mz`,
`mean_cosine`

`mean_cosine` is the mean over that component's edges.

### `singletons.csv`
One row per retained spectrum with degree 0.

`scan`, `precursor_mz`, `n_peaks`

### `precursor_mass_bins.csv`
One row per occupied 100 Da bin.

`bin_start`, `n_spectra`, `n_networked`, `n_singletons`, `mean_degree`

`n_networked` counts spectra of degree 1 or more; `mean_degree` is over all
retained spectra in the bin.
