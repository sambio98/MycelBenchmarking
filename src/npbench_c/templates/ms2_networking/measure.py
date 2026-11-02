"""Measure a networking submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies in `reference/networking_rules.json` -- the confidence classes, the
propagation vocabulary, the threshold grid -- from this template's own closed
column lists, or from the shipped nodes. The node accessions are the one thing
the measurer cannot declare in advance, and it reads them from
`inputs/nodes.json.gz` rather than out of the submission, so a submission that
answered for a node set of its own scores zero instead of being read on its own
terms.

Usage:
  python -m npbench_c.templates.ms2_networking.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.ms2_networking.analyse import COMPONENT_BUCKETS, load_params
from npbench_c.templates.ms2_networking.core import Network

REPORT = "network_report.json"

ENVIRONMENT_KEYS = ("matchms_version", "similarity", "peak_tolerance",
                    "edge_score_threshold", "edge_min_matched_peaks")
CENSUS_KEYS = ("nodes", "edges", "components", "singletons",
               "largest_component", "nodes_with_an_edge", "max_degree")
GRID_KEYS = ("edges", "components", "largest_component")
DROP_KEYS = ("annotations_removed", "components_losing_all_annotation")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _str_list(value: object) -> list[str]:
    return sorted(str(v) for v in value) if isinstance(value, list) else []


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _edge_pairs(value: object) -> list[list[str]]:
    """The submission's edge table, projected onto the claim that is graded.

    An edge row may carry its score and matched-peak count -- the contract asks
    for them -- and neither is read here: a cosine is a float out of a
    BLAS-backed library. What is graded is which pairs clear the declared
    threshold, normalised so that a row's order within the pair and the table's
    order do not decide the answer.
    """
    if not isinstance(value, list):
        return []
    pairs = []
    for row in value:
        if isinstance(row, (list, tuple)) and len(row) >= 2:
            a, b = row[0], row[1]
            if isinstance(a, str) and isinstance(b, str):
                pairs.append(sorted((a, b)))
    return sorted(pairs)


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    network = Network.load(campaign, params["nodes_input"], params["library_input"])
    rules = network.rules
    accessions = sorted(node["accession"] for node in network.nodes)
    classes = sorted(rules.confidence_vocabulary)
    claims = sorted(rules.propagation_vocabulary)
    thresholds = [f"{threshold}" for threshold in rules.threshold_grid]

    path = pathlib.Path(submission_dir) / REPORT
    out: dict = {"present": {"report": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["report"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["report"] = 0
        return out

    # ---- R1: what built the graph, and the graph
    environment = sub.get("environment")
    environment = environment if isinstance(environment, dict) else {}
    out["environment"] = {
        "matchms_version": _str(environment.get("matchms_version")),
        "similarity": _str(environment.get("similarity")),
        "peak_tolerance": _float(environment.get("peak_tolerance")),
        "edge_score_threshold": _float(environment.get("edge_score_threshold")),
        "edge_min_matched_peaks": _int(environment.get("edge_min_matched_peaks")),
    }
    out["network_census"] = _int_map(sub.get("network_census"), CENSUS_KEYS)
    out["edges"] = _edge_pairs(sub.get("edges"))

    # ---- R2: the annotation layer
    out["annotation_census"] = _int_map(sub.get("annotation_census"), classes)
    matches = sub.get("top_matches")
    matches = matches if isinstance(matches, dict) else {}
    out["top_matches"] = {accession: _str(matches.get(accession))
                          for accession in accessions}
    out["component_annotation_census"] = _int_map(
        sub.get("component_annotation_census"), COMPONENT_BUCKETS)

    # ---- R3: what the network claims
    claimed = sub.get("propagation_claims")
    claimed = claimed if isinstance(claimed, dict) else {}
    out["propagation_claims"] = {accession: _str(claimed.get(accession))
                                 for accession in accessions}
    out["propagation_census"] = _int_map(sub.get("propagation_census"), claims)
    shifts = sub.get("mass_shifts")
    shifts = shifts if isinstance(shifts, dict) else {}
    # Only the nodes the submission names, because which nodes carry a shift is
    # part of the claim: a shift reported for a node gold does not list is a
    # wrong answer rather than an ignored one.
    out["mass_shifts"] = {accession: _float(value)
                          for accession, value in sorted(shifts.items())
                          if accession in set(accessions)}

    # ---- R4: the threshold grid and the dropped tier
    grid = sub.get("threshold_grid")
    grid = grid if isinstance(grid, dict) else {}
    out["threshold_grid"] = {}
    for threshold in thresholds:
        row = grid.get(threshold)
        row = row if isinstance(row, dict) else {}
        out["threshold_grid"][threshold] = {
            **_int_map(row, GRID_KEYS),
            "claims": _int_map(row.get("claims"), claims),
        }
    drop = sub.get("dropped_tier")
    drop = drop if isinstance(drop, dict) else {}
    out["dropped_tier"] = {
        "tier_dropped": _str_list(drop.get("tier_dropped")),
        **_int_map(drop, DROP_KEYS),
        "claims": _int_map(drop.get("claims"), claims),
    }
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>",
              file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
