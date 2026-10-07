"""Build gold, warm bundles and the grading spec for the GCF-cutoff ladder.

Usage: python -m npbench_c.templates.gcf_cutoff.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.gcf_cutoff.analyse import analyse, load_params
from npbench_c.templates.gcf_cutoff.core import (
    CONSTANTS,
    RULES,
    CutoffError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("corpus_size", "designed_families", "controls",
                  "bigscape_version", "pfam_version", "pfam_sha256", "cutoffs",
                  "baseline_cutoff", "include_gbk", "class_census")
#: R2's answers: the partition at the declared cutoff.
PARTITION_KEYS = ("baseline_partition", "groups", "size_census",
                  "largest_group", "singletons", "clustered_records")
#: R3's answers: how that partition compares with the chemistry-derived families.
RECONCILIATION_KEYS = ("per_family", "verdict_census", "families_exact",
                       "controls_absorbed", "designed_pairs",
                       "co_clustered_pairs", "recovered_pairs", "missed_pairs",
                       "false_join_pairs", "pair_jaccard")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if gold["corpus_size"] < 10:
        raise CutoffError(
            f"a corpus of {gold['corpus_size']} clusters cannot carry a family "
            "structure worth partitioning")
    if len(gold["designed_families"]) < 3:
        raise CutoffError(
            "fewer than three chemistry-derived families, so the reconciliation "
            "is an anecdote rather than a comparison")
    if not gold["controls"]:
        raise CutoffError(
            "no negative controls, so nothing in the corpus can show the tool "
            "over-merging and the upper end of the sweep means nothing")

    # The reconstruction is the point of R2. If the table happened to list every
    # record, the trap would not exist here and the rule would be inherited from
    # another instantiation rather than exercised.
    if not gold["singletons"]:
        raise CutoffError(
            "every record joined a family at the baseline cutoff, so the "
            "singleton reconstruction R2 turns on is not exercised")
    if gold["groups"] >= gold["corpus_size"]:
        raise CutoffError(
            "the baseline partition is all singletons; the cutoff is too low for "
            "the campaign to be about clustering")

    sweep = gold["cutoff_sweep"]
    if len(sweep) < 4:
        raise CutoffError(
            f"a sweep over {len(sweep)} cutoffs cannot show a dependence")
    jaccards = [row["pair_jaccard"] for row in sweep.values()]
    if max(jaccards) - min(jaccards) < 0.1:
        raise CutoffError(
            f"the pair agreement moves by {max(jaccards) - min(jaccards):.3f} "
            "across the whole grid, so membership does not depend on the cutoff "
            "on this corpus and the campaign has no subject")
    groups = [row["groups"] for row in sweep.values()]
    if len(set(groups)) < 3:
        raise CutoffError(
            "the partition takes fewer than three distinct shapes over the grid")

    # The campaign's claim is that agreement is not monotone: too low fragments
    # the families, too high absorbs a control. Both ends have to be visible or
    # the "best cutoff" is just the largest one.
    if gold["first_cutoff_with_a_false_join"] is None:
        raise CutoffError(
            "no declared cutoff produces a false join, so the grid never reaches "
            "over-merging and the best cutoff is an artefact of where it stops")
    best = gold["best_cutoff"]
    if best == max(sweep, key=float):
        raise CutoffError(
            f"the best cutoff {best} is the highest declared one, so the optimum "
            "may lie outside the grid")
    if gold["best_pair_jaccard"] >= 1.0:
        raise CutoffError(
            "some cutoff reproduces the chemistry-derived families exactly, so "
            "the reconciliation has no disagreement to grade")

    recovered = gold["first_cutoff_recovering_each_family"]
    if all(v is None for v in recovered.values()):
        raise CutoffError(
            "no family is exactly recovered at any declared cutoff, so the "
            "per-family interval claim is empty everywhere")
    if all(v is not None for v in recovered.values()) and len(set(
            v for v in recovered.values())) == 1:
        raise CutoffError(
            "every family is recovered at the same single cutoff, so the "
            "per-family claim is one number repeated")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", PARTITION_KEYS, "partition.json"),
                             ("r4", RECONCILIATION_KEYS, "reconciliation.json")):
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
          ex("corpus_size"), ex("designed_families"), ex("controls"),
          ex("cutoffs"), ex("baseline_cutoff"), ex("include_gbk"),
          ex("class_census"), ex("environment")]

    r2 = [ex("baseline_partition"), ex("group_census")]

    r3 = [ex("per_family"), ex("verdict_census"), ex("families_exact"),
          ex("controls_absorbed"), ex("pair_statistics")]

    r4 = [ex("cutoff_sweep_summary"), ex("sweep_boundaries")]

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


def _flatten_for_grading(gold: dict) -> dict:
    """Group what is one claim into one graded composite.

    The tool identity is one claim -- which BiG-SCAPE and which Pfam produced this
    partition. The group census is one claim about the partition's shape. The pair
    statistics are one claim about agreement, and grading the five counts apart
    would let a system that miscounted lose one component instead of the claim.
    The sweep is graded as the per-cutoff table plus the boundary summary, which
    are different claims and stay apart.
    """
    gold["environment"] = {"bigscape_version": gold["bigscape_version"],
                           "pfam_version": gold["pfam_version"],
                           "pfam_sha256": gold["pfam_sha256"]}
    gold["group_census"] = {k: gold[k] for k in
                            ("groups", "size_census", "largest_group",
                             "singletons", "clustered_records")}
    gold["pair_statistics"] = {k: gold[k] for k in
                               ("designed_pairs", "co_clustered_pairs",
                                "recovered_pairs", "missed_pairs",
                                "false_join_pairs", "pair_jaccard")}
    # The sweep's per-cutoff rows, reduced to the quantities that make the curve.
    gold["cutoff_sweep_summary"] = {
        label: {k: row[k] for k in ("groups", "singletons", "recovered_pairs",
                                    "missed_pairs", "false_join_pairs",
                                    "pair_jaccard", "families_exact")}
        for label, row in sorted(gold["cutoff_sweep"].items(), key=lambda kv: float(kv[0]))}
    gold["sweep_boundaries"] = {
        "best_cutoff": gold["best_cutoff"],
        "best_pair_jaccard": gold["best_pair_jaccard"],
        "first_cutoff_with_a_false_join": gold["first_cutoff_with_a_false_join"],
        "first_cutoff_absorbing_a_control": gold["first_cutoff_absorbing_a_control"],
        "first_cutoff_recovering_each_family":
            gold["first_cutoff_recovering_each_family"],
    }
    return gold


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)

    gold = _flatten_for_grading(analyse(campaign))
    _assert_has_content(gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: {gold['corpus_size']} clusters, "
          f"{len(gold['designed_families'])} chemistry-derived families, "
          f"{len(gold['controls'])} controls")
    print(f"  {gold['bigscape_version']} against Pfam {gold['pfam_version']}, "
          f"cutoffs {gold['cutoffs']}")
    print(f"  baseline c={gold['baseline_cutoff']}: {gold['groups']} groups "
          f"({gold['singletons']} singletons), {gold['families_exact']} families "
          f"exactly recovered, pair Jaccard {gold['pair_jaccard']}")
    print("  per family at the baseline:")
    for name, block in sorted(gold["per_family"].items()):
        print(f"    {name:<16} {block['verdict']:<24} "
              f"spans={block['groups_spanned']} sizes={block['group_sizes']} "
              f"outsiders={block['outsiders']}")
    print("  cutoff sweep:")
    for label, row in sorted(gold["cutoff_sweep_summary"].items(),
                             key=lambda kv: float(kv[0])):
        print(f"    c={label:<4} groups={row['groups']:<3} "
              f"recovered={row['recovered_pairs']:<4} missed={row['missed_pairs']:<4} "
              f"false={row['false_join_pairs']:<3} jaccard={row['pair_jaccard']:<9} "
              f"exact={row['families_exact']}")
    print(f"  best cutoff {gold['best_cutoff']} at Jaccard "
          f"{gold['best_pair_jaccard']}; first false join at "
          f"{gold['first_cutoff_with_a_false_join']}; first control absorbed at "
          f"{gold['first_cutoff_absorbing_a_control']}")
    print(f"  first cutoff recovering each family: "
          f"{gold['first_cutoff_recovering_each_family']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
