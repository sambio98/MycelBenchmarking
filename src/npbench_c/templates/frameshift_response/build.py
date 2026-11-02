"""Build gold, warm bundles and the grading spec for the frameshift ladder.

Usage: python -m npbench_c.templates.frameshift_response.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.frameshift_response.analyse import analyse, load_params
from npbench_c.templates.frameshift_response.core import (
    CONSTANTS,
    RULES,
    FrameshiftError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("accession", "baseline_mode", "noise_mode", "indel_fraction",
                  "fixed_flags", "prodigal_version", "contig_length",
                  "contig_sha256", "baseline_calls", "target_table",
                  "target_length_summary")
#: R2's answers: what the caller did to each broken gene.
RESPONSE_KEYS = ("response_table",)
#: R3's answers: the verdicts.
VERDICT_KEYS = ("verdict_census", "verdict_by_length",
                "verdict_is_monotone_in_length")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("offset_sweep", "noise_arm")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if gold["baseline_calls"] < 500:
        raise FrameshiftError(
            f"only {gold['baseline_calls']} calls on the contig; a frameshift "
            "experiment on a record this quiet is an anecdote")

    targets = gold["target_table"]
    summary = gold["target_length_summary"]
    if summary["targets"] < 5:
        raise FrameshiftError(
            f"{summary['targets']} targets cannot span a length range")
    if summary["longest"] < 20 * summary["shortest"]:
        raise FrameshiftError(
            f"the targets run {summary['shortest']}..{summary['longest']} bp, "
            "less than a factor of twenty; the ladder cannot cross wherever the "
            "caller stops recovering a broken gene")
    spans = {(t["start"], t["end"], t["strand"]) for t in targets.values()}
    if len(spans) != len(targets):
        raise FrameshiftError(
            "two declared lengths selected the same call, so one gene is "
            "audited twice and a length stands for nothing")
    for label, target in sorted(targets.items()):
        if not target["start"] <= target["deleted_at"] <= target["end"]:
            raise FrameshiftError(
                f"target {label}: the deleted base at {target['deleted_at']} is "
                f"outside the gene ({target['start']}..{target['end']}), so the "
                "frameshift is not in the gene it is supposed to break")

    responses = gold["response_table"]
    if set(responses) != set(targets):
        raise FrameshiftError(
            "the response table and the target table describe different gene "
            "sets, so a response belongs to a gene nobody declared")

    # A planted frameshift that leaves the gene intact means the deletion did
    # not land where the campaign thinks it did, or the lift is wrong.
    unchanged = [label for label, row in responses.items()
                 if row["verdict"] == "unchanged"]
    if len(unchanged) == len(responses):
        raise FrameshiftError(
            "every target survives its own frameshift unchanged, which is "
            "impossible unless the deletion or the coordinate lift is wrong")

    census = gold["verdict_census"]
    populated = sorted(v for v, n in census.items() if n)
    if len(populated) < 3:
        raise FrameshiftError(
            f"only {populated} verdicts occur; the vocabulary that separates a "
            "gene the caller recovers from one it loses is not separating "
            "anything on this ladder")
    if not census["split_into_two"]:
        raise FrameshiftError(
            "no target splits in two, which is the response a frameshift in a "
            "long gene guarantees; the classification or the lift is wrong")
    if not census["abolished"]:
        raise FrameshiftError(
            "no target is abolished, so the ladder never reaches a gene too "
            "short for the caller to recover either half of")
    if sum(census.values()) != len(responses):
        raise FrameshiftError(
            f"the verdict census counts {sum(census.values())} responses and "
            f"the table has {len(responses)}")

    # The ladder samples the measured turnover band at a uniform step precisely
    # so that R3 is not a threshold wearing a list's clothing. If every verdict
    # occupies one contiguous block, the per-length list is recoverable from two
    # cut points, and R3 grades one number restated once per target.
    if gold["verdict_is_monotone_in_length"]:
        raise FrameshiftError(
            f"the verdict list {gold['verdict_by_length']} is monotone in "
            "length, so it collapses to a threshold and R3 grades one cut point "
            "restated once per target; the ladder is not crossing the band where "
            "the caller recovers some broken genes and not others")

    sweep = gold["offset_sweep"]
    if len(sweep) < 4:
        raise FrameshiftError(
            f"an offset sweep over {len(sweep)} positions cannot show where the "
            "split tracks the indel")
    splits = [row["new_call_lengths"] for row in sweep.values()]
    if not all(splits):
        raise FrameshiftError(
            "some offset produces no overlapping call at all, so the sweep's "
            "target is not a gene the caller recovers and the split point it is "
            "supposed to track does not exist")
    firsts = [lengths[0] for lengths in splits]
    if len(set(firsts)) < len(firsts):
        raise FrameshiftError(
            f"the sweep's first-fragment lengths repeat ({firsts}), so the split "
            "point does not track the indel and the axis measures nothing")
    if sorted(firsts) != firsts:
        raise FrameshiftError(
            f"the first fragment does not grow as the indel moves along the gene "
            f"({firsts}); the offsets or the coordinates are not in order")

    noise = gold["noise_arm"]
    if not noise:
        raise FrameshiftError(
            "the noise arm is empty, so the campaign asserts its choice of mode "
            "rather than measuring it")
    ratios = [row["noise_to_signal"] for row in noise.values()
              if row["noise_to_signal"] is not None]
    if not ratios:
        raise FrameshiftError(
            "no target in the noise arm explains any call, so the ratio the arm "
            "exists to report is undefined everywhere")
    if max(ratios) < 1.0:
        raise FrameshiftError(
            f"the retraining mode moves fewer calls elsewhere than the target "
            f"explains ({ratios}), so it is not the noisy mode this campaign "
            "declares and the baseline choice needs re-measuring")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", RESPONSE_KEYS, "responses.json"),
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
          ex("contig_length"), ex("contig_sha256"), ex("baseline_calls"),
          ex("target_table"), ex("target_length_summary"), ex("environment")]

    r2 = [ex("response_table")]

    r3 = [ex("verdict_census"), ex("verdict_by_length"),
          ex("verdict_is_monotone_in_length")]

    r4 = [ex("offset_sweep"), ex("noise_arm")]

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
    """The tool identity is one claim: which caller, in which mode, under which
    flags, with which mode held back as the noise arm."""
    gold["environment"] = {"prodigal_version": gold["prodigal_version"],
                           "baseline_mode": gold["baseline_mode"],
                           "noise_mode": gold["noise_mode"],
                           "fixed_flags": gold["fixed_flags"]}
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

    print(f"{campaign.name}: {gold['accession']}, {gold['contig_length']} bp, "
          f"{gold['baseline_calls']} baseline calls")
    print(f"  {gold['prodigal_version']} in {gold['baseline_mode']} mode "
          f"(noise arm: {gold['noise_mode']}), indel at "
          f"{gold['indel_fraction']} of the gene")
    print(f"  {gold['target_length_summary']}")
    print("  response to a planted single-base deletion:")
    for label, row in sorted(gold["response_table"].items(), key=lambda kv: int(kv[0])):
        print(f"    len~{label:<6} actual={row['target_length']:<6} "
              f"{row['verdict']:<17} new={row['new_call_lengths']} "
              f"covered={row['bases_of_the_target_still_called']} "
              f"other_lost={row['other_calls_lost']} "
              f"other_gained={row['other_calls_gained']}")
    print(f"  verdicts: {gold['verdict_census']}")
    print(f"  monotone in length: {gold['verdict_is_monotone_in_length']}")
    print("  offset sweep:")
    for frac, row in sorted(gold["offset_sweep"].items(), key=lambda kv: float(kv[0])):
        print(f"    f={frac:<5} {row['bases_into_the_gene']:>6} bases in  "
              f"{row['verdict']:<17} {row['new_call_lengths']}")
    print("  noise arm (the retraining mode):")
    for label, row in sorted(gold["noise_arm"].items(), key=lambda kv: int(kv[0])):
        print(f"    len~{label:<6} {row['verdict']:<17} "
              f"target_explains={row['calls_explained_by_the_target']:<3} "
              f"other_moved={row['other_calls_moved']:<4} "
              f"noise/signal={row['noise_to_signal']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
