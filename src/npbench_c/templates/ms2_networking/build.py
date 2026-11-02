"""Build gold, warm bundles and the grading spec for the networking ladder.

Usage: python -m npbench_c.templates.ms2_networking.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.ms2_networking.analyse import COMPONENT_BUCKETS, analyse
from npbench_c.templates.ms2_networking.core import (
    CONSTANTS,
    RULES,
    NetworkError,
    Rules,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's report, and therefore exactly what R2's warm bundle may hand over.
EXECUTE_KEYS = ("environment", "network_census", "edges")
#: R2's answers: what the reference library says about the nodes.
ANNOTATION_KEYS = ("annotation_census", "top_matches",
                   "component_annotation_census")
#: R3's answers: what the network claims about the nodes it cannot annotate.
PROPAGATION_KEYS = ("propagation_claims", "propagation_census", "mass_shifts")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("threshold_grid", "dropped_tier")

#: A network of fewer edges than this is a handful of pairs, and its component
#: census says nothing about structure.
MIN_EDGES = 100
#: Below this fraction of nodes in a component of their own, the threshold is so
#: low that the campaign is grading one hairball.
MAX_LARGEST_COMPONENT_FRACTION = 0.25
#: An annotation layer thinner than this cannot seed a propagation census.
MIN_ACCEPTED = 25


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(campaign: pathlib.Path, gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    rules = Rules.load(campaign)

    if gold["environment"]["matchms_version"] != rules.tool_version:
        raise NetworkError(
            f"the pin names matchms {rules.tool_version} and the image answered "
            f"{gold['environment']['matchms_version']}; every score in gold was "
            "produced by a library the campaign does not declare")

    census = gold["network_census"]
    if census["edges"] < MIN_EDGES:
        raise NetworkError(
            f"the network has {census['edges']} edges, fewer than the "
            f"{MIN_EDGES} this template requires; a handful of pairs has no "
            "structure to grade")
    if len(gold["edges"]) != census["edges"]:
        raise NetworkError(
            f"the edge table holds {len(gold['edges'])} pairs and the census "
            f"counts {census['edges']}")
    if census["largest_component"] > MAX_LARGEST_COMPONENT_FRACTION * census["nodes"]:
        raise NetworkError(
            f"the largest component holds {census['largest_component']} of "
            f"{census['nodes']} nodes, more than the declared "
            f"{MAX_LARGEST_COMPONENT_FRACTION:g} fraction; at that density the "
            "campaign grades one hairball and the edge rule needs a neighbour "
            "cap this template does not declare")
    if census["singletons"] == census["nodes"]:
        raise NetworkError(
            "every node is a singleton, so the edge rule found no similarity at "
            "all and there is no network")
    if census["components"] <= 1:
        raise NetworkError(
            "the network is one component, so the component census is a "
            "constant")

    annotations = gold["annotation_census"]
    if sum(annotations.values()) != census["nodes"]:
        raise NetworkError(
            f"the annotation census covers {sum(annotations.values())} nodes "
            f"and the network has {census['nodes']}")
    accepted = sum(annotations[label] for label in rules.accepted)
    if accepted < MIN_ACCEPTED:
        raise NetworkError(
            f"only {accepted} nodes carry an accepted annotation, fewer than "
            f"the {MIN_ACCEPTED} this template requires to seed a propagation "
            "census")
    if accepted == census["nodes"]:
        raise NetworkError(
            "every node carries its own accepted annotation, so nothing is left "
            "for the network to say and the propagation census is a constant")
    empty = sorted(label for label in rules.accepted if not annotations[label])
    if empty:
        raise NetworkError(
            f"the accepted confidence tiers {empty} are empty, so R4's tier drop "
            "removes nothing")
    if len(gold["top_matches"]) != census["nodes"]:
        raise NetworkError(
            f"{len(gold['top_matches'])} nodes carry a top match and the "
            f"network has {census['nodes']}")

    buckets = gold["component_annotation_census"]
    if sum(buckets.values()) != census["components"]:
        raise NetworkError(
            f"the per-component annotation census counts {sum(buckets.values())} "
            f"components and the network has {census['components']}")
    if not buckets["0"] or not buckets["1"]:
        raise NetworkError(
            "the components are either all annotated or all dark, so the "
            "per-component census separates nothing")

    claims = gold["propagation_census"]
    if sum(claims.values()) != census["nodes"]:
        raise NetworkError(
            f"the propagation census covers {sum(claims.values())} nodes and "
            f"the network has {census['nodes']}")
    reached = claims["same_compound"] + claims["analogue_with_mass_shift"]
    if not reached:
        raise NetworkError(
            "no unannotated node is adjacent to an annotated one, so the "
            "propagation this campaign exists to measure never happens")
    if not claims["analogue_with_mass_shift"]:
        raise NetworkError(
            "no node inherits an analogue claim, so every propagation on this "
            "network is a same-compound claim and the vocabulary's two halves "
            "are not both in use")
    if len(gold["mass_shifts"]) != claims["analogue_with_mass_shift"]:
        raise NetworkError(
            f"{len(gold['mass_shifts'])} mass shifts are reported for "
            f"{claims['analogue_with_mass_shift']} analogue claims")
    if len(set(gold["mass_shifts"].values())) < 2:
        raise NetworkError(
            "every analogue claim carries the same mass shift, so the shift "
            "column is a constant")

    grid = gold["threshold_grid"]
    if len(grid) < 3:
        raise NetworkError(
            f"a grid of {len(grid)} thresholds cannot show a response to the "
            "threshold")
    pinned = f"{rules.edge_score_threshold}"
    if grid[pinned]["edges"] != census["edges"]:
        raise NetworkError(
            f"the grid's own {pinned} row reports {grid[pinned]['edges']} edges "
            f"and the network has {census['edges']}, so the grid is not "
            "computed from the same pipeline as the rung above it")
    counts = [row["edges"] for _threshold, row in sorted(grid.items())]
    if len(set(counts)) < len(counts):
        raise NetworkError(
            f"two declared thresholds give the same edge count {counts}, so the "
            "grid has a rung that measures nothing")
    if sorted(counts, reverse=True) != counts:
        raise NetworkError(
            f"the edge count {counts} does not fall as the threshold rises, "
            "which means the grid is not a cut through one similarity matrix")

    dropped = gold["dropped_tier"]
    if not dropped["annotations_removed"]:
        raise NetworkError(
            "the declared tier drop removes no annotation, so R4's second axis "
            "restates R3")
    if dropped["claims"] == gold["propagation_census"]:
        raise NetworkError(
            "dropping a whole confidence tier changes no propagation claim, so "
            "the annotation layer's lower tier reaches nothing")
    if not dropped["components_losing_all_annotation"]:
        raise NetworkError(
            "no component depends on the dropped tier for its only annotation, "
            "so the coverage this axis measures does not move")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", EXECUTE_KEYS, "network.json"),
                             ("r3", ANNOTATION_KEYS, "annotations.json"),
                             ("r4", PROPAGATION_KEYS, "propagation.json")):
        (warm / rung).mkdir(parents=True)
        (warm / rung / name).write_text(
            json.dumps({k: gold[k] for k in keys}, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    grading = task.get("grading") or {}
    chance = {**DEFAULT_CHANCE_LEVELS, **(grading.get("chance_levels") or {})}
    notes = grading.get("chance_level_notes") or {}

    def ex(name, path=None):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{path or name}"}}

    # The edge SCORES are reported and not graded: a cosine is a float out of a
    # BLAS-backed library, and what this campaign claims is which pairs clear the
    # declared threshold. The tool identity is graded as ONE claim that also
    # carries the library's version, which the sandbox does not state.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("environment"), ex("network_census"), ex("edges")]
    r2 = [ex("annotation_census"), ex("top_matches"),
          ex("component_annotation_census")]
    r3 = [ex("propagation_claims"), ex("propagation_census"), ex("mass_shifts")]
    r4 = [ex("threshold_grid"), ex("dropped_tier")]

    rungs = []
    for rid, ordinal, label, components in (("r1", 1, "execute", r1),
                                            ("r2", 2, "correct", r2),
                                            ("r3", 3, "commit", r3),
                                            ("r4", 4, "counterfactual", r4)):
        rung = {"rung_id": rid, "ordinal": ordinal, "name": label, "gold_tier": "G1",
                "composition": "product", "pass_threshold": 1.0,
                "chance_level": chance[rid], "components": components}
        if rid in notes:
            rung["chance_level_note"] = notes[rid]
        rungs.append(rung)

    spec = {"campaign_id": task["campaign_id"], "template_id": task["template_id"],
            "grading_spec_version": "1.0.0", "rungs": rungs}
    (pathlib.Path(campaign) / "grading.json").write_text(
        json.dumps(spec, indent=2) + "\n")
    return spec


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)

    gold = analyse(campaign)
    _assert_has_content(campaign, gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    env, census = gold["environment"], gold["network_census"]
    print(f"{campaign.name}: {census['nodes']} nodes, {census['edges']} edges, "
          f"{census['components']} components (largest "
          f"{census['largest_component']}), matchms {env['matchms_version']} "
          f"{env['similarity']} at {env['edge_score_threshold']}/"
          f"{env['edge_min_matched_peaks']}")
    print(f"  annotations: {gold['annotation_census']}")
    print(f"  components by distinct annotated compounds: "
          f"{ {b: gold['component_annotation_census'][b] for b in COMPONENT_BUCKETS} }")
    print(f"  propagation: {gold['propagation_census']}")
    print(f"  threshold grid: "
          f"{ {t: row['edges'] for t, row in sorted(gold['threshold_grid'].items())} }")
    print(f"  tier drop: {gold['dropped_tier']['annotations_removed']} removed, "
          f"{gold['dropped_tier']['components_losing_all_annotation']} components "
          "go dark")
    print(f"  not graded: {gold['annotation_truth']}, {gold['propagation_truth']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
