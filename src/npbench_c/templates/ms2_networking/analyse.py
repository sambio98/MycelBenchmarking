"""Everything the networking campaign grades, computed from the pinned rules.

The oracle's analysis and the measurer's reading share no code: gold is this
file, and `measure.py` reads a submission. Nothing here reads gold.

Each rung's projection is a module-level function rather than a block inside
`analyse`, because the stub ladder needs exactly the first N of them and a stub
that reimplemented a census would drift from the one gold is built with.
"""

from __future__ import annotations

import collections
import json
import pathlib
from typing import Mapping, Sequence

import yaml

from .core import Network, NetworkError

#: How a component is counted in the per-component annotation census: by how
#: many DISTINCT compounds its accepted annotations name. Three or more is one
#: bucket, because beyond two the number says the same thing.
COMPONENT_BUCKETS = ("0", "1", "2", "3_or_more")


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


# ------------------------------------------------------------------ R1: execute

def environment(network: Network) -> dict:
    """What built the graph. The version is the half the sandbox cannot state."""
    rules = network.rules
    return {
        "matchms_version": network.matchms_version(),
        "similarity": rules.similarity,
        "peak_tolerance": rules.peak_tolerance,
        "edge_score_threshold": rules.edge_score_threshold,
        "edge_min_matched_peaks": rules.edge_min_matched_peaks,
    }


def edge_table(network: Network, threshold: float | None = None) -> list[list[str]]:
    """The declared edge rule as a set of PAIRS, sorted.

    The scores are not here on purpose: a cosine is a float out of a BLAS-backed
    library, and the claim this campaign makes is which pairs clear the declared
    threshold. The scores are reported in the submission and graded nowhere.
    """
    return [[a, b] for a, b, _score, _matches in network.edges(threshold)]


def network_census(network: Network, threshold: float | None = None) -> dict:
    edges = network.edges(threshold)
    components = network.components(threshold)
    degree: collections.Counter = collections.Counter()
    for a, b, *_ in edges:
        degree[a] += 1
        degree[b] += 1
    return {
        "nodes": len(network.nodes),
        "edges": len(edges),
        "components": len(components),
        "singletons": sum(1 for members in components if len(members) == 1),
        "largest_component": len(components[0]) if components else 0,
        "nodes_with_an_edge": len(degree),
        "max_degree": max(degree.values()) if degree else 0,
    }


# ------------------------------------------------------------------ R2: correct

def annotation_layer(network: Network, rows: Sequence[Mapping]) -> dict:
    """What the library says about each node, per node and in aggregate."""
    census = {label: 0 for label in network.rules.confidence_vocabulary}
    for row in rows:
        if row["confidence"] not in census:
            raise NetworkError(
                f"{row['node']} was classified {row['confidence']!r}, which is "
                "not in the declared confidence vocabulary")
        census[row["confidence"]] += 1
    return {
        "annotation_census": {label: census[label] for label in sorted(census)},
        "top_matches": {row["node"]: row["top_match"] for row in rows},
    }


def component_annotation_census(network: Network, rows: Sequence[Mapping],
                                threshold: float | None = None) -> dict:
    """Components by how many distinct compounds their annotations name.

    This is the number that says whether a network is annotated at all: a
    component with none is dark, one with two or more carries a disagreement the
    propagation rule has to resolve.
    """
    seeds = network.accepted(rows)
    buckets = {bucket: 0 for bucket in COMPONENT_BUCKETS}
    for members in network.components(threshold):
        blocks = {seeds[node] for node in members if node in seeds}
        buckets[COMPONENT_BUCKETS[min(len(blocks), 3)]] += 1
    return buckets


# ------------------------------------------------------------------- R3: commit

def propagation(network: Network, rows: Sequence[Mapping],
                threshold: float | None = None,
                seeds: Mapping[str, str] | None = None) -> dict:
    claims = network.claims(rows, threshold, seeds)
    census = {label: 0 for label in network.rules.propagation_vocabulary}
    for record in claims.values():
        census[record["claim"]] += 1
    return {
        "propagation_claims": {node: record["claim"]
                               for node, record in sorted(claims.items())},
        "propagation_census": {label: census[label] for label in sorted(census)},
        "mass_shifts": {node: record["mass_shift"]
                        for node, record in sorted(claims.items())
                        if record["mass_shift"] is not None},
    }


# ----------------------------------------------------------- R4: counterfactual

def threshold_grid(network: Network, rows: Sequence[Mapping]) -> dict:
    """The edge rule is a cut through a continuous similarity, so the campaign
    grades the whole cut. One tool call serves every row of this table."""
    out = {}
    for threshold in network.rules.threshold_grid:
        census = network_census(network, threshold)
        claims = propagation(network, rows, threshold)["propagation_census"]
        out[f"{threshold}"] = {
            "edges": census["edges"],
            "components": census["components"],
            "largest_component": census["largest_component"],
            "claims": claims,
        }
    return out


def dropped_tier(network: Network, rows: Sequence[Mapping]) -> dict:
    """What the network claims when only the top confidence tier seeds it.

    The catalogued counterfactual was "the planted annotation is removed --
    predict the network change". Nothing is planted here, and removing an
    annotation cannot move an edge: the graph is built from spectra. What a
    removal changes is the annotation layer, so the declared removal is a whole
    TIER and the change is measured where it happens.
    """
    rules = network.rules
    top = rules.accepted[0]
    seeds = {row["node"]: row["block"] for row in rows
             if row["confidence"] == top}
    kept = network.accepted(rows)
    dark = 0
    for members in network.components():
        if any(node in kept for node in members) and not any(
                node in seeds for node in members):
            dark += 1
    return {
        "tier_dropped": [label for label in rules.accepted if label != top],
        "annotations_removed": len(kept) - len(seeds),
        "claims": propagation(network, rows, None, seeds)["propagation_census"],
        "components_losing_all_annotation": dark,
    }


# ------------------------------------------------------- gold-side truth checks

def annotation_truth(network: Network, rows: Sequence[Mapping]) -> dict:
    """How many accepted annotations are right. Gold-side, graded nowhere: a
    system cannot read the nodes' identities and nothing asks it to."""
    accepted = [row for row in rows if row["confidence"] in network.rules.accepted]
    right = sum(1 for row in accepted if row["correct"])
    return {"accepted": len(accepted), "correct": right,
            "misannotated": len(accepted) - right}


def propagation_truth(network: Network, rows: Sequence[Mapping]) -> dict:
    """The campaign's finding, measured and graded nowhere.

    An `analogue_with_mass_shift` claim says the node is a RELATIVE of the named
    compound, so it is borne out when the node's own compound differs. A
    `same_compound` claim says the node IS that compound, so it is borne out when
    the compound matches. Measured on the shipped node set the two come out very
    differently, which is what the campaign is about.
    """
    claims = network.claims(rows)
    out = collections.Counter()
    for node, record in claims.items():
        if record["claim"] not in ("same_compound", "analogue_with_mass_shift"):
            continue
        same = record["block"] == network.truth[node]
        if record["claim"] == "same_compound":
            out["same_compound_right" if same else "same_compound_wrong"] += 1
        else:
            out["analogue_a_different_compound" if not same
                else "analogue_the_same_compound"] += 1
    return dict(sorted(out.items()))


# -------------------------------------------------------------------- the whole

def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    network = Network.load(campaign, params["nodes_input"], params["library_input"])
    rows = network.annotations()

    gold = {
        "environment": environment(network),
        "network_census": network_census(network),
        "edges": edge_table(network),
        **annotation_layer(network, rows),
        "component_annotation_census": component_annotation_census(network, rows),
        **propagation(network, rows),
        "threshold_grid": threshold_grid(network, rows),
        "dropped_tier": dropped_tier(network, rows),
    }
    if network.truth is None:
        raise NetworkError(
            "gold cannot be built without the nodes' identities; the template "
            "expects them at gold/nodes_truth.json, written by prepare_input")
    gold["annotation_truth"] = annotation_truth(network, rows)
    gold["propagation_truth"] = propagation_truth(network, rows)
    return gold


def main(argv: list[str] | None = None) -> int:
    import sys

    campaign = pathlib.Path((argv or sys.argv[1:])[0])
    gold = analyse(campaign)
    bulky = ("edges", "top_matches", "propagation_claims", "mass_shifts")
    print(json.dumps({k: v for k, v in gold.items() if k not in bulky}, indent=1))
    print("edges", len(gold["edges"]), "mass shifts", len(gold["mass_shifts"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
