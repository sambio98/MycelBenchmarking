"""Build gold, warm bundles and the grading spec for the gene-calling ladder.

Usage: python -m npbench_c.templates.gene_calling.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.gene_calling.analyse import analyse, load_params
from npbench_c.templates.gene_calling.core import (
    CONSTANTS,
    RULES,
    CallingError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("accession", "baseline_mode", "baseline_fragments",
                  "locality_window", "fixed_flags", "prodigal_version",
                  "contig_length", "contig_sha256", "gc_fraction",
                  "call_census", "partial_census", "start_type_census",
                  "breakpoint_census")
#: R2's answers: the calls themselves.
TABLE_KEYS = ("call_table", "call_span")
#: R3's answers: the same contig cut up.
FRAGMENTATION_KEYS = (
    "fragmented_baseline_calls", "fragmented_fragmented_calls",
    "fragmented_retained", "fragmented_lost", "fragmented_gained",
    "fragmented_retained_but_truncated", "fragmented_loss_verdicts",
    "fragmented_lost_table", "fragmented_gained_partial_census",
    "fragmented_partial_census", "fragmented_call_census")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("fragmentation_grid", "mode_difference")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    census = gold["call_census"]
    if census["calls"] < 500:
        raise CallingError(
            f"only {census['calls']} calls on the whole contig; a fragmentation "
            "experiment on a record this quiet is an anecdote")
    if not census["on_the_forward_strand"] or not census["on_the_reverse_strand"]:
        raise CallingError(
            "every call is on one strand, so the contig is not a genome and the "
            "call table has a free column")
    if not 0.0 < census["coding_fraction"] < 1.0:
        raise CallingError(
            f"the coding fraction is {census['coding_fraction']}; a contig that "
            "is all or none coding is not an assembly")

    # The whole point of giving the agent an unfragmented contig: a complete
    # record has no edge for a call to run off, so every call is whole. If any
    # were partial here, the input is already fragmented and R3's comparison
    # would be between two cut assemblies.
    partial = gold["partial_census"]
    if partial["00"] != census["calls"]:
        raise CallingError(
            f"{census['calls'] - partial['00']} calls on the unfragmented contig "
            "are already flagged partial, so the input is not the whole record "
            "this campaign declares")
    if len(gold["start_type_census"]) < 2:
        raise CallingError(
            "every call has the same start codon, so the start-type census is a "
            "constant dressed up as a distribution")

    # R3 has to move, and it has to move in a way the cuts explain.
    if not gold["fragmented_lost"]:
        raise CallingError(
            "fragmenting the contig loses no call, so the campaign's subject "
            "does not exist at the declared fragmentation")
    if not gold["fragmented_gained"]:
        raise CallingError(
            "fragmenting the contig gains no call, so the fragments' edges "
            "produce no new prediction and the comparison is one-sided")
    if not gold["fragmented_retained"]:
        raise CallingError(
            "no call survives fragmentation, which is a broken run rather than a "
            "finding about assembly quality")
    verdicts = gold["fragmented_loss_verdicts"]
    if not verdicts["crosses_a_cut"]:
        raise CallingError(
            "no lost call spans a breakpoint, which is the one loss the geometry "
            "guarantees; the coordinate lift or the breakpoint rule is wrong")
    populated = [v for v, n in verdicts.items() if n]
    if len(populated) < 2:
        raise CallingError(
            f"only {populated} loss verdicts occur, so the vocabulary that "
            "separates a cut's geometry from the tool's reaction to it is not "
            "separating anything here")
    if len(gold["fragmented_lost_table"]) != gold["fragmented_lost"]:
        raise CallingError(
            f"the lost table carries {len(gold['fragmented_lost_table'])} rows "
            f"and the count says {gold['fragmented_lost']}; one of them is wrong")

    # R4's two axes, each measured rather than assumed.
    grid = gold["fragmentation_grid"]
    if len(grid) < 6:
        raise CallingError(
            f"a grid of {len(grid)} cells cannot cross two modes with a "
            "fragmentation series")
    unfragmented = [cell for label, cell in grid.items() if label.endswith("|1")]
    if not unfragmented:
        raise CallingError(
            "the grid has no unfragmented cell, so nothing in it is a control")
    if any(cell["lost"] or cell["gained"] for cell in unfragmented):
        raise CallingError(
            "an unfragmented cell loses or gains a call against itself, which is "
            "impossible unless the comparison is reading two different runs")
    losses = [cell["lost"] for label, cell in grid.items()
              if not label.endswith("|1")]
    if len(set(losses)) < 3:
        raise CallingError(
            f"the grid's loss counts take {len(set(losses))} distinct values, so "
            "the fragmentation series is not a series")

    difference = gold["mode_difference"]
    if not difference:
        raise CallingError(
            "no mode comparison, so the campaign cannot separate the flag from "
            "the assembly")
    for label, block in difference.items():
        moved = sum(n for key, n in block.items() if key != "shared")
        if not moved:
            raise CallingError(
                f"{label}: the two modes agree on every call from identical "
                "input bytes, so the mode is not a parameter on this contig")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", TABLE_KEYS, "call_table.json"),
                             ("r4", FRAGMENTATION_KEYS, "fragmentation.json")):
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

    # The echoes of the rules file -- the mode, the fragmentation, the window and
    # the flags -- are reported and NOT graded: their value is fixed by the
    # campaign's own declaration.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("contig_length"), ex("contig_sha256"),
          ex("gc_fraction"), ex("call_census"), ex("partial_census"),
          ex("start_type_census"), ex("breakpoint_census"), ex("environment")]

    r2 = [ex("call_table"), ex("call_span")]

    r3 = [ex("fragmentation_summary"), ex("fragmented_loss_verdicts"),
          ex("fragmented_lost_table"), ex("fragmented_partial_census"),
          ex("fragmented_gained_partial_census"), ex("fragmented_call_census")]

    r4 = [ex("fragmentation_grid"), ex("mode_difference")]

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

    The tool identity is one claim -- which Prodigal, at which mode, under which
    flags. The fragmentation counts are one claim about what the cuts did, and
    grading the five apart would let a system that mis-lifted every coordinate
    lose one component instead of the claim.
    """
    gold["environment"] = {"prodigal_version": gold["prodigal_version"],
                           "baseline_mode": gold["baseline_mode"],
                           "fixed_flags": gold["fixed_flags"]}
    gold["fragmentation_summary"] = {
        k: gold[f"fragmented_{k}"] for k in
        ("baseline_calls", "fragmented_calls", "retained", "lost", "gained",
         "retained_but_truncated")}
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

    census = gold["call_census"]
    print(f"{campaign.name}: {gold['accession']}, {gold['contig_length']} bp, "
          f"GC {gold['gc_fraction']:.3f}")
    print(f"  {gold['prodigal_version']} in {gold['baseline_mode']} mode, flags "
          f"{' '.join(gold['fixed_flags'])}")
    print(f"  {census['calls']} calls ({census['on_the_forward_strand']}+ / "
          f"{census['on_the_reverse_strand']}-), "
          f"{census['coding_fraction']:.1%} of the contig coding, median length "
          f"{census['median_length']}")
    print(f"  partial flags on the whole contig: {gold['partial_census']}")
    print(f"  start types: {gold['start_type_census']}")
    print(f"  cut into {gold['baseline_fragments']}: "
          f"{gold['fragmentation_summary']}")
    print(f"  loss verdicts: {gold['fragmented_loss_verdicts']}")
    print(f"  partial flags after fragmentation: "
          f"{gold['fragmented_partial_census']}")
    print("  grid (mode|fragments -> lost / gained / truncated):")
    for label, cell in sorted(gold["fragmentation_grid"].items()):
        print(f"    {label:<14} lost={cell['lost']:<4} gained={cell['gained']:<4} "
              f"truncated={cell['retained_but_truncated']:<4} "
              f"{ {k: v for k, v in cell['loss_verdicts'].items() if v} }")
    print(f"  mode difference on identical bytes: {gold['mode_difference']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
