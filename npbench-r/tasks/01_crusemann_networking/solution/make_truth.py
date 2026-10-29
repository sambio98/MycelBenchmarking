"""Reference solution for tasks/01_crusemann_networking.

Classic GNPS-style molecular networking over the clustered MS/MS spectra
deposited with the Crusemann marine actinomycete survey, PROJECT.md §4
category 3.

    python solution/make_truth.py --task tasks/01_crusemann_networking --out <dir>

Why the clustered .mgf and not the raw files. MSV000078836 holds 1,125 mzXML
files totalling 311 GB, fifteen times the 20 GB target in PROJECT.md §10, so
reprocessing it is out of budget. The deposit also contains the GNPS-clustered
consensus spectra as a single 17 MB .mgf, and §9 names exactly this as the
safer artifact: "Network topology from a fixed .mgf with fixed parameters is
reproducible". The clustering is the depositors' own, frozen and checksummed in
the manifest, so the task starts where the reproducible part starts.

The file's own mess is left in and is the main difficulty. Of its 100,263
spectra, 96% carry no fragment peaks at all -- 89,923 peak lines across 100,263
spectra, and spectrum 2 is literally empty. Only 4,081 survive a minimum-peak
filter. An agent that does not notice will network noise or time out.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matchms
import networkx as nx
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosineGreedy

MGF = "data/MSV000078836_specs_ms.mgf"

# matchms renamed ModifiedCosine to ModifiedCosineGreedy, so the pin is
# load-bearing in a way a version number alone does not convey.
MATCHMS_VERSION = "0.33.1"

# Classic GNPS molecular networking parameters.
PRECURSOR_WINDOW_DA = 17.0    # peaks within this of the precursor are removed
MIN_PEAKS = 6
FRAGMENT_TOLERANCE = 0.5
MIN_COSINE = 0.70
MIN_MATCHED_PEAKS = 6
MAX_PRECURSOR_DELTA = 2.0
TOP_K = 10
MASS_BIN_DA = 100
DECIMALS = 4


def check_version() -> None:
    if matchms.__version__ != MATCHMS_VERSION:
        raise SystemExit(
            f"this task's truth was produced by matchms {MATCHMS_VERSION}, found "
            f"{matchms.__version__}. The similarity class was renamed between "
            f"releases and cosine scores are implementation-dependent."
        )


def load_and_filter(path: Path) -> tuple[list, int]:
    """Retained spectra sorted by scan, plus how many were loaded.

    Two of the three filter steps are no-ops on this file, measured rather than
    assumed, and that is a property of the deposit worth knowing:

      * 96,182 of the 100,263 spectra carry no peaks at all, and the peak-count
        histogram jumps straight from 0 to 6 -- nothing has between 1 and 5
        peaks. So the MIN_PEAKS filter never removes anything the empty-spectrum
        check did not already remove. It stays because the artifact needs a
        defined rule, not because it is hard.
      * The +/-17 Da precursor-window removal changes only 81 of the 4,081
        spectra with peaks, dropping 83 peaks in total and pushing none below
        six. It is kept because it is standard practice and it does move
        n_peaks, base_peak_mz and a few cosines -- a small, easily missed
        difference rather than a no-op.
    """
    loaded = list(load_from_mgf(str(path)))
    retained = []
    for spectrum in loaded:
        if spectrum is None or len(spectrum.peaks.mz) == 0:
            continue
        normalised = normalize_intensities(spectrum)
        if normalised is None:
            continue
        trimmed = remove_peaks_around_precursor_mz(
            normalised, mz_tolerance=PRECURSOR_WINDOW_DA
        )
        if trimmed is None or len(trimmed.peaks.mz) < MIN_PEAKS:
            continue
        retained.append(trimmed)

    retained.sort(key=lambda s: int(s.get("scans")))
    return retained, len(loaded)


def candidate_edges(spectra: list) -> list[dict]:
    """Modified cosine over pairs inside the precursor-mass window."""
    similarity = ModifiedCosineGreedy(tolerance=FRAGMENT_TOLERANCE)
    scans = [int(s.get("scans")) for s in spectra]
    masses = [float(s.get("precursor_mz")) for s in spectra]
    edges = []
    for i in range(len(spectra)):
        for j in range(i + 1, len(spectra)):
            delta = masses[j] - masses[i]
            if abs(delta) > MAX_PRECURSOR_DELTA:
                continue
            result = similarity.pair(spectra[i], spectra[j])
            score = float(result["score"])
            matches = int(result["matches"])
            if score >= MIN_COSINE and matches >= MIN_MATCHED_PEAKS:
                a, b = (scans[i], scans[j]) if scans[i] < scans[j] else (scans[j], scans[i])
                edges.append({
                    "scan_a": a, "scan_b": b,
                    "cosine": round(score, DECIMALS),
                    "n_matched_peaks": matches,
                    "precursor_delta": round(abs(delta), DECIMALS),
                })
    return edges


def apply_top_k(edges: list[dict]) -> list[dict]:
    """GNPS's top-K rule: an edge survives only if it is in the top K of BOTH ends.

    Ranking is by descending cosine, then ascending partner scan so the cut is
    deterministic when scores tie -- and with rounded cosines they tie often.
    """
    neighbours: dict[int, list[tuple[float, int, dict]]] = defaultdict(list)
    for edge in edges:
        neighbours[edge["scan_a"]].append((-edge["cosine"], edge["scan_b"], edge))
        neighbours[edge["scan_b"]].append((-edge["cosine"], edge["scan_a"], edge))

    kept_by: dict[int, set[tuple[int, int]]] = {}
    for scan, group in neighbours.items():
        group.sort(key=lambda item: (item[0], item[1]))
        kept_by[scan] = {
            (e["scan_a"], e["scan_b"]) for _, _, e in group[:TOP_K]
        }
    return [
        edge for edge in edges
        if (edge["scan_a"], edge["scan_b"]) in kept_by[edge["scan_a"]]
        and (edge["scan_a"], edge["scan_b"]) in kept_by[edge["scan_b"]]
    ]


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    check_version()
    spectra, n_loaded = load_and_filter(args.task / MGF)
    print(f"{n_loaded} spectra loaded, {len(spectra)} retained", flush=True)

    scan_of = [int(s.get("scans")) for s in spectra]
    write_csv(
        args.out / "spectra.csv",
        ["scan", "precursor_mz", "charge", "n_peaks", "base_peak_mz"],
        [
            {
                "scan": int(s.get("scans")),
                "precursor_mz": round(float(s.get("precursor_mz")), DECIMALS),
                "charge": int(s.get("charge") or 0),
                "n_peaks": len(s.peaks.mz),
                # The most intense surviving fragment: a cheap check that the
                # normalise-then-trim order was followed.
                "base_peak_mz": round(
                    float(s.peaks.mz[int(s.peaks.intensities.argmax())]), DECIMALS),
            }
            for s in spectra
        ],
    )

    raw = candidate_edges(spectra)
    edges = apply_top_k(raw)
    print(f"{len(raw)} candidate edges, {len(edges)} after top-{TOP_K}", flush=True)
    edges.sort(key=lambda e: (e["scan_a"], e["scan_b"]))
    write_csv(args.out / "edges.csv",
              ["scan_a", "scan_b", "cosine", "n_matched_peaks", "precursor_delta"], edges)

    graph = nx.Graph()
    graph.add_nodes_from(scan_of)
    for edge in edges:
        graph.add_edge(edge["scan_a"], edge["scan_b"], cosine=edge["cosine"])

    component_of: dict[int, int] = {}
    members: dict[int, list[int]] = {}
    for component in nx.connected_components(graph):
        name = min(component)
        members[name] = sorted(component)
        for scan in component:
            component_of[scan] = name

    write_csv(
        args.out / "network_nodes.csv",
        ["scan", "degree", "component_id", "component_size"],
        [{"scan": scan, "degree": graph.degree(scan),
          "component_id": component_of[scan],
          "component_size": len(members[component_of[scan]])}
         for scan in scan_of],
    )

    mz_of = {int(s.get("scans")): float(s.get("precursor_mz")) for s in spectra}
    edges_by_component: dict[int, list[dict]] = defaultdict(list)
    for edge in edges:
        edges_by_component[component_of[edge["scan_a"]]].append(edge)
    # Only components with an edge are reported; a lone spectrum is a singleton,
    # which has its own artifact.
    write_csv(
        args.out / "components.csv",
        ["component_id", "n_nodes", "n_edges", "min_precursor_mz", "max_precursor_mz",
         "mean_cosine"],
        [
            {
                "component_id": name,
                "n_nodes": len(group),
                "n_edges": len(edges_by_component[name]),
                "min_precursor_mz": round(min(mz_of[s] for s in group), DECIMALS),
                "max_precursor_mz": round(max(mz_of[s] for s in group), DECIMALS),
                "mean_cosine": round(
                    sum(e["cosine"] for e in edges_by_component[name])
                    / len(edges_by_component[name]), DECIMALS),
            }
            for name, group in sorted(members.items()) if len(group) > 1
        ],
    )

    singletons = [s for s in scan_of if graph.degree(s) == 0]
    write_csv(
        args.out / "singletons.csv", ["scan", "precursor_mz", "n_peaks"],
        [{"scan": int(s.get("scans")),
          "precursor_mz": round(float(s.get("precursor_mz")), DECIMALS),
          "n_peaks": len(s.peaks.mz)}
         for s in spectra if graph.degree(int(s.get("scans"))) == 0],
    )
    print(f"{len(members)} components, {len(singletons)} singletons", flush=True)

    bins: dict[int, list[int]] = defaultdict(list)
    for scan in scan_of:
        bins[int(mz_of[scan] // MASS_BIN_DA) * MASS_BIN_DA].append(scan)
    write_csv(
        args.out / "precursor_mass_bins.csv",
        ["bin_start", "n_spectra", "n_networked", "n_singletons", "mean_degree"],
        [
            {
                "bin_start": start,
                "n_spectra": len(group),
                "n_networked": sum(1 for s in group if graph.degree(s) > 0),
                "n_singletons": sum(1 for s in group if graph.degree(s) == 0),
                "mean_degree": round(
                    sum(graph.degree(s) for s in group) / len(group), DECIMALS),
            }
            for start, group in sorted(bins.items())
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
