"""Everything the annotation-transfer ladder grades, from the sandbox alone.

One search, then arithmetic. The counterfactual needs no second search: removing
a hit from the reference set means reading further down a ranking already
computed, which is why this campaign's R4 is cheap where its catalogued form
looked expensive.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.annotation_transfer.core import Transfer


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    transfer = Transfer.load(campaign, params["sequence_input"],
                             params["annotation_input"])
    rules = transfer.rules

    hits = transfer.hits()
    verdicts = transfer.verdicts(hits)
    census = transfer.verdict_census(verdicts)

    removal = {}
    transitions = {}
    for name, remove in sorted(rules.removal_variants.items()):
        removal[name] = transfer.census_after_removal(hits, remove)
        transitions[name] = (transfer.verdict_transitions(hits, remove)
                             if remove else {})

    graded = transfer.agreement_by_depth(hits, rules.graded_depths)
    constant = transfer.agreement_by_depth(hits, rules.constant_depths)
    deepest = max(rules.graded_depths)

    errors = sorted(q for q, v in verdicts.items()
                    if v == "wrong_at_subsubclass")

    return {
        "release": transfer.release,
        "ec_prefix": rules.ec_prefix,
        "diamond_version": transfer.diamond_version(),
        "tie_break": rules.tie_break,
        "entries_total": len(transfer.sequences),
        **transfer.annotation_census(),
        "queries_with_a_hit": len(hits),
        "queries_without_a_hit": len(transfer.sequences) - len(hits),
        "tied_top_queries": transfer.tied_top_queries(hits),
        "top_hit_identity_census": transfer.identity_census(hits),
        "agreement_by_depth": graded,
        "agreement_by_identity": transfer.agreement_by_identity(hits, deepest),
        "verdict_census": census,
        "transfer_errors_at_subsubclass": errors,
        "removal_census": removal,
        "removal_changed": {
            name: sum(moves.values()) for name, moves in sorted(transitions.items())},
        # Audit material, not graded: see excluded_from_grading.
        "constant_depth_agreement": constant,
        "removal_transitions": transitions,
    }
