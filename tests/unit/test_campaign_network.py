"""MS2 molecular networking: grading what the network claims, not what is true.

The catalogue asked for a planted annotation and a confound vocabulary. Neither
survived the grounding, and the measurements that retired them are the reason
this campaign is shaped the way it is -- so the tests below pin those
measurements alongside the arithmetic.

Four things here produce a wrong answer without producing an error.

An edge means two spectra fragment alike. At the same precursor that supports
"the same compound"; at a different precursor it supports "an analogue of it",
with a mass shift. Treating every edge as the first is the mistake molecular
networking invites, and it is one of the two competence probes.

Seeding the network from every library match rather than from the declared
confidence tiers is the other: it reads as diligence and it contaminates the
component census.

Two near-identical spectra at the same precursor are NOT necessarily the same
compound. Measured on this node set, two of the four disagreeing pairs of that
kind are genuinely different compounds -- which is why no rung asks which
annotation is unsupported.

And a node carries no identity. Whether an annotation is right is gold's
business; the sandbox cannot tell, and nothing graded asks it to.
"""

from __future__ import annotations

import copy
import gzip
import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.templates.ms2_dereplication.prepare_input import QUERY_FIELDS
from npbench_c.templates.ms2_networking import analyse as net
from npbench_c.templates.ms2_networking import core as core_module
from npbench_c.templates.ms2_networking.build import (
    ANNOTATION_KEYS,
    COUNTERFACTUAL_KEYS,
    EXECUTE_KEYS,
    PROPAGATION_KEYS,
    _assert_has_content,
)
from npbench_c.templates.ms2_networking.core import Network, NetworkError, Rules
from npbench_c.templates.ms2_networking.measure import measure
from npbench_c.templates.ms2_networking.solve import REPORT, render
from npbench_c.tools.resources import tool_prefix

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "network-massbank-2026_03-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


def _spec() -> dict:
    return json.loads((CAMPAIGN / "grading.json").read_text())


def _provenance() -> dict:
    return json.loads((CAMPAIGN / "inputs" / "provenance.json").read_text())


@pytest.fixture(scope="module")
def network() -> Network:
    params = net.load_params(CAMPAIGN)
    return Network.load(CAMPAIGN, params["nodes_input"], params["library_input"])


# ------------------------------------------------------------ what is shipped


def test_only_the_two_corpora_and_their_provenance_are_shipped():
    """The network is not shipped: it follows from the declared edge rule over
    the shipped nodes, and shipping it would hand over R1 entirely."""
    assert sorted(p.name for p in (CAMPAIGN / "inputs").iterdir()) == [
        "library.json.gz", "nodes.json.gz", "provenance.json"]


def test_a_shipped_node_carries_no_identity():
    params = net.load_params(CAMPAIGN)
    nodes = json.loads(gzip.open(CAMPAIGN / params["nodes_input"], "rt").read())
    for row in nodes:
        assert sorted(row) == sorted(QUERY_FIELDS)
        assert not [f for f in core_module.IDENTITY_FIELDS if f in row]
    library = json.loads(gzip.open(CAMPAIGN / params["library_input"], "rt").read())
    assert all("compound_block" in row for row in library)


def test_the_truth_is_in_gold_and_graded_nowhere(network):
    from npbench_c.sweep.runner import PUBLIC_CAMPAIGN_PATHS, SECRET_CAMPAIGN_PATHS

    truth = pathlib.PurePath(core_module.TRUTH)
    assert truth.parts[0] in SECRET_CAMPAIGN_PATHS
    assert truth.parts[0] not in PUBLIC_CAMPAIGN_PATHS
    assert network.truth and len(network.truth) == len(network.nodes)

    gold = _gold()
    assert {"annotation_truth", "propagation_truth"} <= set(gold)
    graded = {component["args"]["gold"] for rung in _spec()["rungs"]
              for component in rung["components"]
              if isinstance(component["args"]["gold"], str)}
    assert "gold:annotation_truth" not in graded
    assert "gold:propagation_truth" not in graded
    assert "annotation_truth" not in render(gold)


def test_the_node_set_is_one_laboratory_on_one_instrument(network):
    params = net.load_params(CAMPAIGN)
    assert {row["contributor"] for row in network.nodes} == {params["contributor"]}
    assert {row["instrument_type"] for row in network.nodes} == \
        {params["instrument"]}


def test_the_library_excludes_the_node_contributor(network):
    params = net.load_params(CAMPAIGN)
    assert params["contributor"] not in {row["contributor"]
                                         for row in network.library}
    assert len({row["contributor"] for row in network.library}) >= 5


def test_the_library_filter_is_declared_and_cannot_change_an_answer(network):
    """Every shipped library spectrum is inside some node's window, and the
    provenance records what the filter was applied to."""
    window = network.rules.precursor_window
    masses = sorted(row["precursor_mz"] for row in network.nodes)
    import bisect

    for row in network.library:
        index = bisect.bisect_left(masses, row["precursor_mz"] - window)
        assert index < len(masses)
        assert masses[index] <= row["precursor_mz"] + window, row["accession"]
    census = _provenance()["census"]
    assert census["library_spectra"] == len(network.library)
    assert (census["cross_contributor_spectra_before_the_window_filter"]
            > census["library_spectra"])


def test_the_loader_refuses_a_node_set_that_mixes_laboratories(network, tmp_path):
    params = net.load_params(CAMPAIGN)
    nodes = json.loads(gzip.open(CAMPAIGN / params["nodes_input"], "rt").read())
    nodes[0] = {**nodes[0], "contributor": "SOMEWHERE_ELSE"}
    campaign = _sandbox(tmp_path, nodes, [dict(row) for row in network.library])
    with pytest.raises(NetworkError, match="a network is one"):
        Network.load(campaign, "inputs/nodes.json.gz", "inputs/library.json.gz")


def test_the_loader_refuses_a_node_that_names_its_own_compound(network, tmp_path):
    params = net.load_params(CAMPAIGN)
    nodes = json.loads(gzip.open(CAMPAIGN / params["nodes_input"], "rt").read())
    nodes[0] = {**nodes[0], "compound_block": "WHATEVERTHISIS"}
    campaign = _sandbox(tmp_path, nodes, [dict(row) for row in network.library])
    with pytest.raises(NetworkError, match="which name the compound each"):
        Network.load(campaign, "inputs/nodes.json.gz", "inputs/library.json.gz")


def test_the_loader_refuses_a_library_from_the_node_contributor(network, tmp_path):
    params = net.load_params(CAMPAIGN)
    nodes = json.loads(gzip.open(CAMPAIGN / params["nodes_input"], "rt").read())
    library = [dict(row) for row in network.library]
    library[0] = {**library[0], "contributor": params["contributor"]}
    campaign = _sandbox(tmp_path, nodes, library)
    with pytest.raises(NetworkError, match="own laboratory's replicate"):
        Network.load(campaign, "inputs/nodes.json.gz", "inputs/library.json.gz")


def _sandbox(tmp_path, nodes, library) -> pathlib.Path:
    campaign = tmp_path / "campaign"
    (campaign / "inputs").mkdir(parents=True)
    (campaign / "reference").mkdir()
    (campaign / "reference" / "networking_rules.json").write_text(
        (CAMPAIGN / "reference" / "networking_rules.json").read_text())
    with gzip.open(campaign / "inputs" / "nodes.json.gz", "wt") as handle:
        json.dump(nodes, handle)
    with gzip.open(campaign / "inputs" / "library.json.gz", "wt") as handle:
        json.dump(library, handle)
    return campaign


# ---------------------------------------------------------------- the network


def test_the_edge_table_is_pairs_and_the_scores_are_not_graded():
    gold = _gold()
    assert all(len(edge) == 2 for edge in gold["edges"])
    assert all(a < b for a, b in gold["edges"])
    assert gold["edges"] == sorted(gold["edges"])
    assert len(gold["edges"]) == gold["network_census"]["edges"]
    excluded = " ".join(entry["field"] for entry
                        in _task()["excluded_from_grading"])
    assert "cosine score" in excluded


@image_present
def test_no_edge_score_sits_on_a_declared_threshold(network):
    """The measurement that makes a pair set safe to grade. One pair sits
    5.6e-05 below the highest declared threshold, which is the closest call in
    the campaign and still far above the float noise of a pinned stack."""
    margins = []
    for threshold in network.rules.threshold_grid:
        for _a, _b, score, _matches in network.pairs():
            margins.append((abs(score - threshold), threshold))
    closest, threshold = min(margins)
    assert closest > 1e-5, (closest, threshold)
    assert sum(1 for margin, _ in margins if margin < 1e-4) == 1


def test_the_component_census_covers_every_node(network):
    gold = _gold()
    components = network.components()
    assert sum(len(members) for members in components) == \
        gold["network_census"]["nodes"]
    assert len(components) == gold["network_census"]["components"]
    assert len(components[0]) == gold["network_census"]["largest_component"]
    assert sum(gold["component_annotation_census"].values()) == len(components)


def test_no_component_is_a_hairball():
    """A low threshold produces one giant component, which is why the classical
    convention caps neighbours and component size. The declared threshold needs
    neither here, and the build refuses gold if that stops being true."""
    from npbench_c.templates.ms2_networking.build import (
        MAX_LARGEST_COMPONENT_FRACTION,
    )

    census = _gold()["network_census"]
    assert census["largest_component"] <= (MAX_LARGEST_COMPONENT_FRACTION
                                           * census["nodes"])


# ------------------------------------------------------------- the annotation


def test_the_annotation_census_covers_every_node_and_both_silent_states():
    gold = _gold()
    census = gold["annotation_census"]
    assert sum(census.values()) == gold["network_census"]["nodes"]
    assert census["no_candidate"] and census["no_match"]
    assert len(gold["top_matches"]) == gold["network_census"]["nodes"]


def test_a_silent_node_has_no_top_match():
    gold = _gold()
    silent = [node for node, match in gold["top_matches"].items() if match is None]
    census = gold["annotation_census"]
    assert len(silent) == census["no_candidate"] + census["no_match"]


def test_only_the_declared_tiers_seed_the_network(network):
    rules = network.rules
    assert rules.accepted == ("confirmed", "tentative")
    assert rules.classify(0.95, 6) == "confirmed"
    assert rules.classify(0.95, 5) == "tentative"
    assert rules.classify(0.4, 20) == "unsupported"
    assert "unsupported" not in rules.accepted


# ------------------------------------------------------------- the propagation


def test_every_node_carries_exactly_one_declared_claim(network):
    gold = _gold()
    vocabulary = set(network.rules.propagation_vocabulary)
    assert set(gold["propagation_census"]) == vocabulary
    assert set(gold["propagation_claims"].values()) <= vocabulary
    assert len(gold["propagation_claims"]) == gold["network_census"]["nodes"]
    assert sum(gold["propagation_census"].values()) == \
        gold["network_census"]["nodes"]
    tally = {label: 0 for label in vocabulary}
    for claim in gold["propagation_claims"].values():
        tally[claim] += 1
    assert tally == gold["propagation_census"]


def test_a_mass_shift_is_reported_for_an_analogue_claim_and_nothing_else():
    gold = _gold()
    analogues = {node for node, claim in gold["propagation_claims"].items()
                 if claim == "analogue_with_mass_shift"}
    assert set(gold["mass_shifts"]) == analogues
    assert len(analogues) == gold["propagation_census"]["analogue_with_mass_shift"]


def test_a_mass_shift_is_exact_arithmetic_over_the_shipped_precursors(network):
    """Not tool output, so it is reproducible by hand: the shift is the node's
    precursor minus its best-scoring annotated neighbour's."""
    gold = _gold()
    precursor = {row["accession"]: row["precursor_mz"] for row in network.nodes}
    window = network.rules.precursor_window
    for node, shift in gold["mass_shifts"].items():
        assert abs(shift) > window, node
        candidates = [round(precursor[node] - precursor[other],
                            network.rules.mass_shift_precision)
                      for other in precursor]
        assert shift in candidates


def test_a_same_compound_claim_stays_inside_the_annotation_window(network):
    gold = _gold()
    assert gold["propagation_census"]["same_compound"]
    for node, claim in gold["propagation_claims"].items():
        if claim == "same_compound":
            assert node not in gold["mass_shifts"]


def test_the_claims_the_truth_bears_out_are_the_analogue_ones():
    """The campaign's finding, measured and graded nowhere: the analogue claim is
    borne out every time it is made, and the same-compound claim is wrong more
    often than it is right. That asymmetry is why the declared rule separates the
    two cases that networking practice conflates."""
    truth = _gold()["propagation_truth"]
    census = _gold()["propagation_census"]
    assert truth["analogue_a_different_compound"] == \
        census["analogue_with_mass_shift"]
    assert "analogue_the_same_compound" not in truth
    assert truth["same_compound_wrong"] > truth["same_compound_right"]


def test_the_misannotation_rate_is_real_and_not_graded():
    """A fifth of the accepted annotations name the wrong compound, which is the
    mis-seeding the catalogue asked about. No component reads it: a system cannot
    know, and the clause set that would have found it does not work -- see the
    exclusion."""
    truth = _gold()["annotation_truth"]
    assert truth["misannotated"] > 0.1 * truth["accepted"]
    excluded = " ".join(entry["field"] for entry
                        in _task()["excluded_from_grading"])
    assert "which annotation is unsupported" in excluded


# ------------------------------------------------------------ the counterfactual


def test_the_grid_holds_the_declared_thresholds_and_its_own_identity_rung(network):
    gold = _gold()
    grid = gold["threshold_grid"]
    assert sorted(grid) == sorted(f"{t}" for t in network.rules.threshold_grid)
    pinned = f"{network.rules.edge_score_threshold}"
    assert grid[pinned]["edges"] == gold["network_census"]["edges"]
    assert grid[pinned]["components"] == gold["network_census"]["components"]
    assert grid[pinned]["claims"] == gold["propagation_census"]


def test_the_grid_falls_monotonically_and_moves_every_number():
    grid = _gold()["threshold_grid"]
    rows = [grid[key] for key in sorted(grid)]
    edges = [row["edges"] for row in rows]
    assert edges == sorted(edges, reverse=True)
    assert len(set(edges)) == len(edges)
    components = [row["components"] for row in rows]
    assert components == sorted(components)
    reached = [row["claims"]["same_compound"]
               + row["claims"]["analogue_with_mass_shift"] for row in rows]
    assert reached == sorted(reached, reverse=True)


def test_the_tier_drop_removes_annotations_and_darkens_components():
    gold = _gold()
    drop = gold["dropped_tier"]
    assert drop["tier_dropped"] == ["tentative"]
    assert drop["annotations_removed"] == gold["annotation_census"]["tentative"]
    assert drop["claims"] != gold["propagation_census"]
    assert drop["claims"]["annotated_directly"] == \
        gold["annotation_census"]["confirmed"]
    assert drop["components_losing_all_annotation"] > 0
    assert drop["claims"]["no_annotated_neighbour"] > \
        gold["propagation_census"]["no_annotated_neighbour"]


# ------------------------------------------------------------ gold and grading


def test_the_rung_key_sets_partition_the_report():
    sets = [set(keys) for keys in (EXECUTE_KEYS, ANNOTATION_KEYS,
                                   PROPAGATION_KEYS, COUNTERFACTUAL_KEYS)]
    for i, left in enumerate(sets):
        for right in sets[i + 1:]:
            assert not left & right
    assert set().union(*sets) == set(render(_gold()))


def test_every_graded_component_resolves_in_gold():
    gold = _gold()
    for rung in _spec()["rungs"]:
        for component in rung["components"]:
            reference = component["args"]["gold"]
            if isinstance(reference, str) and reference.startswith("gold:"):
                assert reference.split(":", 1)[1] in gold


def test_the_oracle_submission_scores_one(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / REPORT).write_text(json.dumps(render(_gold())))
    result = grade(measure(CAMPAIGN, submission), _gold(), _spec())
    assert result["score"] == 1.0
    assert result["depth"] == 4


def test_a_missing_report_scores_nothing(tmp_path):
    submission = tmp_path / "submission"
    submission.mkdir()
    assert grade(measure(CAMPAIGN, submission), _gold(), _spec())["depth"] == 0


def test_the_edge_table_is_read_as_a_set_of_pairs(tmp_path):
    """Order inside a pair and order of the table are not the claim; the pair
    set is. A missing pair is."""
    gold = _gold()
    report = render(gold)
    report["edges"] = [[b, a, 0.0, 0] for a, b in reversed(gold["edges"])]
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / REPORT).write_text(json.dumps(report))
    assert grade(measure(CAMPAIGN, submission), gold, _spec())["score"] == 1.0

    report["edges"] = [[a, b] for a, b in gold["edges"][:-1]]
    (submission / REPORT).write_text(json.dumps(report))
    assert grade(measure(CAMPAIGN, submission), gold, _spec())["depth"] == 0


def test_the_measurer_reads_the_node_set_from_the_inputs(tmp_path):
    report = render(_gold())
    report["propagation_claims"] = {f"INVENTED-{i}": "same_compound"
                                    for i in range(478)}
    submission = tmp_path / "submission"
    submission.mkdir()
    (submission / REPORT).write_text(json.dumps(report))
    measurement = measure(CAMPAIGN, submission)
    assert set(measurement["propagation_claims"]) == set(_gold()["propagation_claims"])
    assert all(v is None for v in measurement["propagation_claims"].values())
    assert grade(measurement, _gold(), _spec())["depth"] == 2


# ----------------------------------------------------------- build refusals


def _a_different_library_version() -> dict:
    gold = copy.deepcopy(_gold())
    gold["environment"]["matchms_version"] = "0.0.1"
    return gold


def _a_thin_network() -> dict:
    gold = copy.deepcopy(_gold())
    gold["network_census"]["edges"] = 12
    gold["edges"] = gold["edges"][:12]
    return gold


def _a_miscounted_edge_table() -> dict:
    gold = copy.deepcopy(_gold())
    gold["edges"] = gold["edges"][:-1]
    return gold


def _a_hairball() -> dict:
    gold = copy.deepcopy(_gold())
    gold["network_census"]["largest_component"] = \
        gold["network_census"]["nodes"] // 2
    return gold


def _one_component() -> dict:
    gold = copy.deepcopy(_gold())
    gold["network_census"]["components"] = 1
    return gold


def _every_node_a_singleton() -> dict:
    gold = copy.deepcopy(_gold())
    gold["network_census"]["singletons"] = gold["network_census"]["nodes"]
    return gold


def _a_short_top_match_table() -> dict:
    gold = copy.deepcopy(_gold())
    gold["top_matches"].pop(sorted(gold["top_matches"])[0])
    return gold


def _a_miscounted_component_census() -> dict:
    gold = copy.deepcopy(_gold())
    gold["component_annotation_census"]["0"] += 1
    return gold


def _a_miscounted_annotation_census() -> dict:
    gold = copy.deepcopy(_gold())
    gold["annotation_census"]["confirmed"] += 1
    return gold


def _an_empty_accepted_tier() -> dict:
    gold = copy.deepcopy(_gold())
    census = gold["annotation_census"]
    census["unsupported"] += census["tentative"]
    census["tentative"] = 0
    return gold


def _a_fully_annotated_network() -> dict:
    gold = copy.deepcopy(_gold())
    nodes = gold["network_census"]["nodes"]
    gold["annotation_census"] = {"confirmed": nodes - 1, "tentative": 1,
                                 "unsupported": 0, "no_candidate": 0,
                                 "no_match": 0}
    return gold


def _a_dark_component_census() -> dict:
    gold = copy.deepcopy(_gold())
    buckets = gold["component_annotation_census"]
    buckets["0"] += buckets["1"]
    buckets["1"] = 0
    return gold


def _a_miscounted_propagation_census() -> dict:
    gold = copy.deepcopy(_gold())
    gold["propagation_census"]["same_compound"] += 1
    return gold


def _a_network_that_propagates_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    census = gold["propagation_census"]
    census["no_annotated_neighbour"] += (census["same_compound"]
                                         + census["analogue_with_mass_shift"])
    census["same_compound"] = 0
    census["analogue_with_mass_shift"] = 0
    gold["mass_shifts"] = {}
    return gold


def _no_analogue_claims() -> dict:
    gold = copy.deepcopy(_gold())
    census = gold["propagation_census"]
    census["same_compound"] += census["analogue_with_mass_shift"]
    census["analogue_with_mass_shift"] = 0
    gold["mass_shifts"] = {}
    return gold


def _a_miscounted_shift_table() -> dict:
    gold = copy.deepcopy(_gold())
    gold["mass_shifts"].popitem()
    return gold


def _one_mass_shift_everywhere() -> dict:
    gold = copy.deepcopy(_gold())
    gold["mass_shifts"] = {node: 14.0157 for node in gold["mass_shifts"]}
    return gold


def _a_two_point_grid() -> dict:
    gold = copy.deepcopy(_gold())
    keep = sorted(gold["threshold_grid"])[:2]
    gold["threshold_grid"] = {k: gold["threshold_grid"][k] for k in keep}
    return gold


def _a_grid_that_disagrees_with_its_own_rung() -> dict:
    gold = copy.deepcopy(_gold())
    gold["threshold_grid"]["0.7"]["edges"] += 1
    return gold


def _a_grid_that_does_not_fall() -> dict:
    gold = copy.deepcopy(_gold())
    grid = gold["threshold_grid"]
    grid["0.9"]["edges"] = grid["0.6"]["edges"] + 1
    return gold


def _a_grid_with_a_repeated_row() -> dict:
    gold = copy.deepcopy(_gold())
    grid = gold["threshold_grid"]
    grid["0.8"]["edges"] = grid["0.7"]["edges"]
    return gold


def _a_tier_drop_that_removes_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    gold["dropped_tier"]["annotations_removed"] = 0
    return gold


def _a_tier_drop_that_changes_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    gold["dropped_tier"]["claims"] = dict(gold["propagation_census"])
    return gold


def _a_tier_drop_that_darkens_nothing() -> dict:
    gold = copy.deepcopy(_gold())
    gold["dropped_tier"]["components_losing_all_annotation"] = 0
    return gold


@pytest.mark.parametrize("mutate,pattern", [
    (_a_different_library_version, "a library the campaign does not declare"),
    (_a_thin_network, "fewer than the"),
    (_a_miscounted_edge_table, "the edge table holds"),
    (_a_hairball, "more than the declared"),
    (_one_component, "the component census is a"),
    (_every_node_a_singleton, "every node is a singleton"),
    (_a_short_top_match_table, "nodes carry a top match"),
    (_a_miscounted_component_census, "per-component annotation census counts"),
    (_a_miscounted_annotation_census, "the annotation census covers"),
    (_an_empty_accepted_tier, "accepted confidence tiers"),
    (_a_fully_annotated_network, "nothing is left"),
    (_a_dark_component_census, "separates nothing"),
    (_a_miscounted_propagation_census, "the propagation census covers"),
    (_a_network_that_propagates_nothing, "never happens"),
    (_no_analogue_claims, "no node inherits an analogue claim"),
    (_a_miscounted_shift_table, "mass shifts are reported for"),
    (_one_mass_shift_everywhere, "the shift\n?.*column is a constant"),
    (_a_two_point_grid, "cannot show a response to the"),
    (_a_grid_that_disagrees_with_its_own_rung, "not\n?.*computed from the same"),
    (_a_grid_that_does_not_fall, "does not fall as the threshold rises"),
    (_a_grid_with_a_repeated_row, "give the same edge count"),
    (_a_tier_drop_that_removes_nothing, "removes no annotation"),
    (_a_tier_drop_that_changes_nothing, "changes no propagation claim"),
    (_a_tier_drop_that_darkens_nothing, "no component depends on the dropped"),
])
def test_the_build_refuses_gold_that_says_nothing(mutate, pattern):
    with pytest.raises(NetworkError, match=pattern):
        _assert_has_content(CAMPAIGN, mutate())


def test_the_real_gold_passes_every_refusal():
    _assert_has_content(CAMPAIGN, _gold())


def test_the_rules_refuse_a_grid_without_the_pinned_threshold(tmp_path):
    document = json.loads(
        (CAMPAIGN / "reference" / "networking_rules.json").read_text())
    document["network"]["threshold_grid"] = [0.5, 0.6, 0.8]
    campaign = tmp_path / "campaign"
    (campaign / "reference").mkdir(parents=True)
    (campaign / "reference" / "networking_rules.json").write_text(
        json.dumps(document))
    with pytest.raises(NetworkError, match="has no identity rung"):
        Rules.load(campaign)
