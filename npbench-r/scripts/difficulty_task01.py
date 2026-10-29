"""Measure task 01's difficulty by grading plausible wrong pipelines.

Same approach and caveat as the task 02 to 06 scripts: synthetic variants, a
cheap proxy for the model runs PROJECT.md §8 asks for, useful mainly for
finding artifacts nothing can fail.

Candidate edges are cached per scoring configuration, since several variants
change only what happens after scoring.

One shortcut is worth stating because it is exact rather than an
approximation: the "keeps the empty spectra" variant does NOT need the 62
million comparisons it looks like it needs. A spectrum with no peaks can never
reach six matched peaks, so it can never form an edge. Keeping the empties
changes every node-level artifact and leaves the edge list untouched, which is
what the variant below computes.

    python scripts/difficulty_task01.py --task tasks/01_crusemann_networking \
        --truth tasks/01_crusemann_networking/truth --scratch <dir>
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import networkx as nx  # noqa: E402
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz  # noqa: E402
from matchms.importing import load_from_mgf  # noqa: E402
from matchms.similarity import CosineGreedy, ModifiedCosineGreedy  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

MISTAKES = {
    "keeps_empty_spectra": "never discards the 96% of spectra with no peaks",
    "skips_precursor_window": "omits the +/-17 Da precursor-window removal",
    "plain_cosine": "uses CosineGreedy instead of ModifiedCosineGreedy",
    "tolerance_0_1": "scores with a 0.1 Da fragment tolerance instead of 0.5",
    "threshold_0_6": "links at cosine 0.6 instead of 0.7",
    "no_min_matched_peaks": "drops the six-matched-peak requirement",
    "no_top_k": "skips the top-10 cut entirely",
    "top_k_either_end": "keeps an edge in the top 10 of either end, not both",
    "precursor_window_5da": "pairs spectra up to 5 Da apart instead of 2",
    "components_numbered": "numbers components C0001.. instead of naming them by scan",
    "components_include_singletons": "reports one-node components in components.csv",
    "bins_by_rounding": "bins precursor mass by rounding instead of flooring",
    "charge_defaults_to_one": "reads a missing CHARGE as 1+ rather than unknown",
    "base_peak_highest_mz": "reads the base peak as the highest m/z, not the most intense",
}

_SPECTRA: dict[bool, list] = {}
_EDGES: dict[tuple, list[dict]] = {}


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref01", task_dir / "solution" / "make_truth.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def spectra_for(ref, task_dir: Path, mistakes: set[str]):
    """Filtered spectra. Cached on whether the precursor window was applied."""
    windowed = "skips_precursor_window" not in mistakes
    if windowed in _SPECTRA:
        return _SPECTRA[windowed]
    print(f"    (parsing and filtering; window={windowed})", flush=True)
    kept = []
    for spectrum in load_from_mgf(str(task_dir / ref.MGF)):
        if spectrum is None or len(spectrum.peaks.mz) == 0:
            continue
        normalised = normalize_intensities(spectrum)
        if normalised is None:
            continue
        if windowed:
            normalised = remove_peaks_around_precursor_mz(
                normalised, mz_tolerance=ref.PRECURSOR_WINDOW_DA)
            if normalised is None:
                continue
        if len(normalised.peaks.mz) < ref.MIN_PEAKS:
            continue
        kept.append(normalised)
    kept.sort(key=lambda s: int(s.get("scans")))
    _SPECTRA[windowed] = kept
    return kept


def edges_for(ref, spectra, mistakes: set[str]) -> list[dict]:
    """Candidate edges, cached per scoring configuration."""
    key = (
        "skips_precursor_window" in mistakes,
        "plain_cosine" in mistakes,
        "tolerance_0_1" in mistakes,
        "precursor_window_5da" in mistakes,
        "threshold_0_6" in mistakes,
        "no_min_matched_peaks" in mistakes,
    )
    if key in _EDGES:
        return _EDGES[key]
    print(f"    (scoring for config {key})", flush=True)
    tolerance = 0.1 if "tolerance_0_1" in mistakes else ref.FRAGMENT_TOLERANCE
    similarity = (CosineGreedy(tolerance=tolerance) if "plain_cosine" in mistakes
                  else ModifiedCosineGreedy(tolerance=tolerance))
    max_delta = 5.0 if "precursor_window_5da" in mistakes else ref.MAX_PRECURSOR_DELTA
    min_cosine = 0.6 if "threshold_0_6" in mistakes else ref.MIN_COSINE
    min_matches = 0 if "no_min_matched_peaks" in mistakes else ref.MIN_MATCHED_PEAKS

    scans = [int(s.get("scans")) for s in spectra]
    masses = [float(s.get("precursor_mz")) for s in spectra]
    edges = []
    for i in range(len(spectra)):
        for j in range(i + 1, len(spectra)):
            delta = masses[j] - masses[i]
            if abs(delta) > max_delta:
                continue
            result = similarity.pair(spectra[i], spectra[j])
            score = float(result["score"])
            matches = int(result["matches"])
            if score >= min_cosine and matches >= min_matches:
                a, b = (scans[i], scans[j]) if scans[i] < scans[j] else (scans[j], scans[i])
                edges.append({"scan_a": a, "scan_b": b,
                              "cosine": round(score, ref.DECIMALS),
                              "n_matched_peaks": matches,
                              "precursor_delta": round(abs(delta), ref.DECIMALS)})
    _EDGES[key] = edges
    return edges


def top_k(ref, edges: list[dict], mistakes: set[str]) -> list[dict]:
    if "no_top_k" in mistakes:
        return list(edges)
    neighbours: dict[int, list] = defaultdict(list)
    for edge in edges:
        neighbours[edge["scan_a"]].append((-edge["cosine"], edge["scan_b"], edge))
        neighbours[edge["scan_b"]].append((-edge["cosine"], edge["scan_a"], edge))
    kept_by = {}
    for scan, group in neighbours.items():
        group.sort(key=lambda item: (item[0], item[1]))
        kept_by[scan] = {(e["scan_a"], e["scan_b"]) for _, _, e in group[:ref.TOP_K]}
    if "top_k_either_end" in mistakes:
        return [e for e in edges
                if (e["scan_a"], e["scan_b"]) in kept_by[e["scan_a"]]
                or (e["scan_a"], e["scan_b"]) in kept_by[e["scan_b"]]]
    return [e for e in edges
            if (e["scan_a"], e["scan_b"]) in kept_by[e["scan_a"]]
            and (e["scan_a"], e["scan_b"]) in kept_by[e["scan_b"]]]


def charge_of(spectrum, mistakes: set[str]) -> int:
    """Charge, where 1,741 of the 4,081 spectra carry no CHARGE line at all.

    matchms leaves those at 0. Defaulting them to 1+ -- the usual assumption for
    positive-mode ESI, and a reading rather than a typo -- collapses two of the
    three classes macro F1 is averaged over.
    """
    raw = spectrum.get("charge")
    if raw:
        return int(raw)
    return 1 if "charge_defaults_to_one" in mistakes else 0


def base_peak_of(spectrum, mistakes: set[str]) -> float:
    """The m/z of the most intense peak, or of the highest peak if misread."""
    mz = spectrum.peaks.mz
    if not len(mz):
        return 0.0
    if "base_peak_highest_mz" in mistakes:
        return float(mz.max())
    return float(mz[int(spectrum.peaks.intensities.argmax())])


def build(task_dir: Path, out: Path, mistakes: set[str]) -> None:
    ref = load_solution(task_dir)
    spectra = spectra_for(ref, task_dir, mistakes)
    edges = top_k(ref, edges_for(ref, spectra, mistakes), mistakes)
    edges = sorted(edges, key=lambda e: (e["scan_a"], e["scan_b"]))

    # The node universe. An empty spectrum can never reach six matched peaks,
    # so keeping the empties adds isolated nodes and no edges.
    node_spectra = list(spectra)
    if "keeps_empty_spectra" in mistakes:
        extra = []
        for spectrum in load_from_mgf(str(task_dir / ref.MGF)):
            if spectrum is not None and len(spectrum.peaks.mz) == 0:
                extra.append(spectrum)
        node_spectra = sorted(spectra + extra, key=lambda s: int(s.get("scans")))

    w = ref.write_csv
    w(out / "spectra.csv",
      ["scan", "precursor_mz", "charge", "n_peaks", "base_peak_mz"],
      [{"scan": int(s.get("scans")),
        "precursor_mz": round(float(s.get("precursor_mz")), ref.DECIMALS),
        "charge": charge_of(s, mistakes),
        "n_peaks": len(s.peaks.mz),
        "base_peak_mz": round(base_peak_of(s, mistakes), ref.DECIMALS)}
       for s in node_spectra])
    w(out / "edges.csv",
      ["scan_a", "scan_b", "cosine", "n_matched_peaks", "precursor_delta"], edges)

    scan_of = [int(s.get("scans")) for s in node_spectra]
    graph = nx.Graph()
    graph.add_nodes_from(scan_of)
    for edge in edges:
        graph.add_edge(edge["scan_a"], edge["scan_b"], cosine=edge["cosine"])

    members: dict[int, list[int]] = {}
    component_of: dict[int, int] = {}
    for component in nx.connected_components(graph):
        name = min(component)
        members[name] = sorted(component)
        for scan in component:
            component_of[scan] = name
    label = dict.fromkeys(members, None)
    if "components_numbered" in mistakes:
        order = sorted(members, key=lambda n: (-len(members[n]), n))
        label = {n: f"C{i + 1:04d}" for i, n in enumerate(order)}
    else:
        label = {n: n for n in members}

    w(out / "network_nodes.csv", ["scan", "degree", "component_id", "component_size"],
      [{"scan": s, "degree": graph.degree(s),
        "component_id": label[component_of[s]],
        "component_size": len(members[component_of[s]])} for s in scan_of])

    mz_of = {int(s.get("scans")): float(s.get("precursor_mz")) for s in node_spectra}
    edges_by: dict[int, list[dict]] = defaultdict(list)
    for edge in edges:
        edges_by[component_of[edge["scan_a"]]].append(edge)
    floor = 1 if "components_include_singletons" in mistakes else 2
    w(out / "components.csv",
      ["component_id", "n_nodes", "n_edges", "min_precursor_mz", "max_precursor_mz",
       "mean_cosine"],
      [{"component_id": label[name], "n_nodes": len(group),
        "n_edges": len(edges_by[name]),
        "min_precursor_mz": round(min(mz_of[s] for s in group), ref.DECIMALS),
        "max_precursor_mz": round(max(mz_of[s] for s in group), ref.DECIMALS),
        "mean_cosine": round(sum(e["cosine"] for e in edges_by[name])
                             / len(edges_by[name]), ref.DECIMALS) if edges_by[name] else 0.0}
       for name, group in sorted(members.items()) if len(group) >= floor])

    w(out / "singletons.csv", ["scan", "precursor_mz", "n_peaks"],
      [{"scan": int(s.get("scans")),
        "precursor_mz": round(float(s.get("precursor_mz")), ref.DECIMALS),
        "n_peaks": len(s.peaks.mz)}
       for s in node_spectra if graph.degree(int(s.get("scans"))) == 0])

    bins: dict[int, list[int]] = defaultdict(list)
    for scan in scan_of:
        if "bins_by_rounding" in mistakes:
            start = int(round(mz_of[scan] / ref.MASS_BIN_DA)) * ref.MASS_BIN_DA
        else:
            start = int(mz_of[scan] // ref.MASS_BIN_DA) * ref.MASS_BIN_DA
        bins[start].append(scan)
    w(out / "precursor_mass_bins.csv",
      ["bin_start", "n_spectra", "n_networked", "n_singletons", "mean_degree"],
      [{"bin_start": start, "n_spectra": len(group),
        "n_networked": sum(1 for s in group if graph.degree(s) > 0),
        "n_singletons": sum(1 for s in group if graph.degree(s) == 0),
        "mean_degree": round(sum(graph.degree(s) for s in group) / len(group),
                             ref.DECIMALS)}
       for start, group in sorted(bins.items())])


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: keeps the empty spectra": {"keeps_empty_spectra"},
    "one slip: skips the precursor window": {"skips_precursor_window"},
    "one slip: plain cosine, not modified": {"plain_cosine"},
    "one slip: 0.1 Da fragment tolerance": {"tolerance_0_1"},
    "one slip: links at cosine 0.6": {"threshold_0_6"},
    "one slip: no matched-peak minimum": {"no_min_matched_peaks"},
    "one slip: no top-10 cut": {"no_top_k"},
    "one slip: top-10 of either end": {"top_k_either_end"},
    "one slip: 5 Da precursor window": {"precursor_window_5da"},
    "one slip: components numbered": {"components_numbered"},
    "one slip: singletons in components": {"components_include_singletons"},
    "one slip: bins by rounding": {"bins_by_rounding"},
    "one slip: missing charge is 1+": {"charge_defaults_to_one"},
    "one slip: base peak is highest m/z": {"base_peak_highest_mz"},
    "hasty agent (3 slips)": {"skips_precursor_window", "top_k_either_end",
                              "bins_by_rounding"},
    "careless agent (5 slips)": {"keeps_empty_spectra", "plain_cosine", "no_top_k",
                                 "components_numbered", "bins_by_rounding"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    config = load_grading_config(args.task / "grading.yaml")

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':40}{'task':>6}  " + "".join(f"{n[:12]:>14}" for n in names)
    print(header)
    print("-" * len(header))
    passes = Counter()
    for sweep_label, mistakes in VARIANTS.items():
        out = args.scratch / sweep_label.replace(" ", "_").replace(":", "").replace(",", "")
        build(args.task, out, mistakes)
        result = score_task(config, args.truth, out)
        # A soft artifact is scored and reported but stays out of the headline,
        # so the columns are rebuilt in config order from both lists.
        scored = {a.path: a for a in result.artifacts + result.soft_artifacts}
        row = f"{sweep_label:40}{result.task_score:6.2f}  "
        for spec in config.artifacts:
            artifact = scored[spec.path]
            row += f"{artifact.score:14.2f}"
            if artifact.passed:
                passes[artifact.path] += 1
        print(row, flush=True)

    print()
    print(f"artifact pass rate across the {len(VARIANTS)} variants "
          f"(PROJECT.md section 8 wants 0.1 to 0.8):")
    for artifact in config.artifacts:
        rate = passes[artifact.path] / len(VARIANTS)
        # The band is a property of the headline, so a soft artifact is reported
        # without a verdict rather than judged against it.
        verdict = ("soft, out of the headline" if artifact.soft
                   else "ok" if 0.1 <= rate <= 0.8 else "OUTSIDE BAND")
        print(f"  {artifact.path:34} {rate:4.2f}  {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
