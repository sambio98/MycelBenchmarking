"""Build gold, warm bundles and the grading spec for the annotation-transfer ladder.

Usage: python -m npbench_c.templates.annotation_transfer.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.annotation_transfer.analyse import analyse, load_params
from npbench_c.templates.annotation_transfer.core import (
    CONSTANTS,
    RULES,
    TransferError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.002, "r3": 0.001, "r4": 0.001}

SET_KEYS = ("release", "ec_prefix", "diamond_version", "tie_break",
            "entries_total", "distinct_ec_codes", "ec_codes_per_entry",
            "partial_ec_mentions", "entries_with_multiple_ec")
SEARCH_KEYS = ("queries_with_a_hit", "queries_without_a_hit",
               "tied_top_queries", "top_hit_identity_census")
VERDICT_KEYS = ("agreement_by_depth", "agreement_by_identity", "verdict_census",
                "transfer_errors_at_subsubclass")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    census = gold["verdict_census"]
    populated = [v for v, n in census.items() if n]
    if len(populated) < 3:
        raise TransferError(
            f"only {populated} verdicts occur; a partition with two cells is a "
            "boolean dressed up as an enum")
    wrong = census["wrong_at_serial"] + census["wrong_at_subsubclass"]
    if wrong < 20:
        raise TransferError(
            f"only {wrong} transfer errors in the whole set, so the campaign's "
            "subject barely exists here")
    if not census["supported"]:
        raise TransferError(
            "no transfer is supported, which would mean the search or the "
            "annotation join is broken rather than that transfer is unreliable")

    depth = gold["agreement_by_depth"]
    deepest = max(depth, key=lambda k: int(k))
    rate = depth[deepest]["rate"]
    if rate is None or not 0.0 < rate < 1.0:
        raise TransferError(
            f"agreement at depth {deepest} is {rate}; a rung whose answer is 0 or "
            "1 measures nothing")

    # The claim the identity bands exist to support. If accuracy did not fall off
    # with identity, the bands would be decoration and R3 would be one number.
    bands = [b for b in gold["agreement_by_identity"].values()
             if b["decidable"] >= 10 and b["rate"] is not None]
    if len(bands) < 4:
        raise TransferError(
            "fewer than four identity bands carry ten decidable queries, so the "
            "accuracy-against-identity curve is not supported by this set")
    first, last = bands[0]["rate"], bands[-1]["rate"]
    if not last > first:
        raise TransferError(
            f"agreement does not rise with identity ({first} at the bottom band "
            f"against {last} at the top), so the campaign's central claim does "
            "not hold on this set")

    identity = gold["removal_census"].get("v0_identity")
    if identity != census:
        raise TransferError(
            "the zero-removal variant does not reproduce the headline census, so "
            "R4's identity control is not controlling anything")
    moved = [n for name, n in gold["removal_changed"].items()
             if name != "v0_identity"]
    if not moved or max(moved) < 50:
        raise TransferError(
            f"removing the closest hit changes at most {max(moved) if moved else 0} "
            "verdicts, so the counterfactual shows nothing")
    if gold["removal_changed"].get("v0_identity"):
        raise TransferError(
            "the zero-removal variant moved a verdict, which is impossible unless "
            "the removal logic is reading the wrong ranking")

    if not gold["tied_top_queries"]:
        raise TransferError(
            "no query has a tied top bitscore, so the declared tie-break is a "
            "no-op here; that may be fine, but it must be re-examined rather "
            "than inherited from another instantiation")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", SET_KEYS, "set_inventory.json"),
                             ("r3", SEARCH_KEYS, "search_inventory.json"),
                             ("r4", VERDICT_KEYS, "verdicts.json")):
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

    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("release"), ex("ec_prefix"), ex("tie_break"), ex("entries_total"),
          ex("distinct_ec_codes"), ex("ec_codes_per_entry"),
          ex("partial_ec_mentions"), ex("entries_with_multiple_ec"),
          ex("diamond_version")]

    r2 = [ex("queries_with_a_hit"), ex("queries_without_a_hit"),
          ex("tied_top_queries"), ex("top_hit_identity_census")]

    r3 = [ex("agreement_by_depth"), ex("agreement_by_identity"),
          ex("verdict_census"), ex("transfer_errors_at_subsubclass")]

    r4 = [ex("removal_census"), ex("removal_changed")]

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
    _assert_has_content(gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: UniProt {gold['release']}, EC {gold['ec_prefix']}.*, "
          f"{gold['entries_total']} reviewed entries")
    print(f"  {gold['diamond_version']}, {gold['queries_with_a_hit']} queries with "
          f"a non-self hit, {gold['tied_top_queries']} with a tied top bitscore")
    print(f"  {gold['partial_ec_mentions']} EC mentions stop early; "
          f"{gold['entries_with_multiple_ec']} entries carry more than one")
    print("  agreement by EC depth:")
    for depth, block in sorted(gold["agreement_by_depth"].items()):
        print(f"    level {depth}: {block['agree']}/{block['decidable']} "
              f"= {block['rate']}")
    print("  agreement by top-hit identity:")
    for band, block in gold["agreement_by_identity"].items():
        if block["decidable"]:
            print(f"    {band:>9}%  n={block['decidable']:<5} rate={block['rate']}")
    print(f"  verdicts: {gold['verdict_census']}")
    print(f"  removal counterfactual: changed {gold['removal_changed']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
