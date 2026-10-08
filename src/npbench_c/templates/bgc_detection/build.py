"""Build gold, warm bundles and the grading spec for the BGC-detection ladder.

Usage: python -m npbench_c.templates.bgc_detection.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.bgc_detection.analyse import analyse, load_params
from npbench_c.templates.bgc_detection.core import (
    CONSTANTS,
    RULES,
    DetectionError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.002, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("accession", "record_length", "antismash_version",
                  "detection_rule_fingerprint", "coordinate_convention",
                  "search_flags", "regions", "distinct_products",
                  "product_census", "contig_edge_regions", "region_bases",
                  "fraction_of_record")
#: R2's answers: the per-region boundaries themselves.
TABLE_KEYS = ("region_table", "region_span")
RECONCILIATION_KEYS = ("loci_total", "verdict_census", "per_locus",
                       "regions_matching_a_locus", "regions_without_a_locus",
                       "regions_holding_several_loci", "jaccard_median",
                       "jaccard_min", "jaccard_max")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if gold["regions"] < 5:
        raise DetectionError(
            f"only {gold['regions']} regions called; a detection campaign on a "
            "record this quiet is an anecdote")

    # The whole reason the input is an unsliced chromosome. A region touching the
    # record edge has its boundary set by the record, not by the tool, and the
    # boundary is what R2 and R3 grade.
    if gold["contig_edge_regions"]:
        raise DetectionError(
            f"{gold['contig_edge_regions']} regions report contig_edge=True, so "
            "their boundaries are set by the record rather than by the detection "
            "rules; the input is not whole or not the one this campaign declares")

    if gold["distinct_products"] < 5:
        raise DetectionError(
            f"only {gold['distinct_products']} distinct products, so the product "
            "census is nearly a constant")

    census = gold["verdict_census"]
    populated = [v for v, n in census.items() if n]
    if len(populated) < 3:
        raise DetectionError(
            f"only {populated} verdicts occur; an extent vocabulary with two "
            "cells is a boolean dressed up as an enum")
    if census["not_detected"] == gold["loci_total"]:
        raise DetectionError(
            "no curated locus is detected at all, which is a broken run rather "
            "than a finding about extent")

    if not 0.0 < gold["jaccard_median"] < 1.0:
        raise DetectionError(
            f"the median Jaccard is {gold['jaccard_median']}; a rung whose answer "
            "is 0 or 1 measures nothing")
    if gold["jaccard_max"] - gold["jaccard_min"] < 0.2:
        raise DetectionError(
            "every curated locus agrees with its region to the same degree, so "
            "the extent vocabulary is not distinguishing anything")

    outcomes = gold["deletion_outcomes"]
    identity = outcomes.get("v0_identity")
    if identity is None:
        raise DetectionError("no identity control among the deletion variants")
    if (identity["regions"] != gold["regions"] or identity["regions_lost"]
            or identity["regions_gained"]):
        raise DetectionError(
            f"the zero-deletion variant does not reproduce the baseline "
            f"({identity}), so R4's identity control is not controlling anything")
    variants = {name: o for name, o in outcomes.items() if name != "v0_identity"}
    if not variants:
        raise DetectionError("no deletion variant besides the identity control")
    if not any(o["regions_lost"] or o["regions_gained"] for o in variants.values()):
        raise DetectionError(
            "no declared deletion changes the region set, so R4 shows nothing")

    # The count alone is not the signal. Each variant targets one region, so
    # "one region lost" is true of most of them by construction; what
    # distinguishes them is WHICH product goes and whether the region survives at
    # all. The signature carries both.
    signatures = {
        (o["regions_lost"], o["regions_gained"], tuple(o["lost_products"]),
         tuple(o["gained_products"])) for o in variants.values()}
    if len(signatures) < 3:
        raise DetectionError(
            f"the deletion variants produce only {len(signatures)} distinct "
            "outcomes, so they are not probing different dependencies")

    # A negative control: at least one declared deletion must leave the region
    # set alone. Without it a system can score this rung by answering "one region
    # is lost" to every variant.
    if not any(not o["regions_lost"] and not o["regions_gained"]
               for o in variants.values()):
        raise DetectionError(
            "every declared deletion changes the region set, so R4 has no "
            "negative control and 'a deletion costs a region' is a free answer")
    if not any(o["regions_lost"] and not o["regions_gained"]
               for o in variants.values()):
        raise DetectionError(
            "no declared deletion abolishes a region outright, so the rung never "
            "shows detection failing")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", TABLE_KEYS, "region_table.json"),
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

    # R1 is the inventory: run the tool and summarise what came back. R2 is the
    # boundaries themselves. The split matters because R2's warm bundle hands over
    # R1's answers, so a quantity graded at R2 must not be part of the inventory --
    # product_census, region_bases and fraction_of_record all were, and the
    # warm-bundle check caught the first of them.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("accession"), ex("record_length"), ex("coordinate_convention"),
          ex("search_flags"), ex("regions"), ex("distinct_products"),
          ex("product_census"), ex("contig_edge_regions"), ex("region_bases"),
          ex("fraction_of_record"), ex("environment")]

    r2 = [ex("region_table"), ex("region_span")]

    r3 = [ex("verdict_census"), ex("per_locus"), ex("loci_total"),
          ex("regions_without_a_locus"), ex("regions_holding_several_loci"),
          ex("jaccard_summary")]

    r4 = [ex("deletion_summary")]

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
    """Group the pinned identity and the two summaries into graded composites.

    The tool version and the rule fingerprint are one claim -- which rules
    produced this region set -- and grading them apart would let a system that
    pinned neither lose one component instead of the claim.
    """
    gold["environment"] = {
        "antismash_version": gold["antismash_version"],
        "detection_rule_fingerprint": gold["detection_rule_fingerprint"],
    }
    gold["jaccard_summary"] = {"median": gold["jaccard_median"],
                               "min": gold["jaccard_min"],
                               "max": gold["jaccard_max"]}
    gold["deletion_summary"] = {
        name: {k: outcome[k] for k in
               ("regions", "regions_lost", "regions_gained", "regions_unchanged",
                "lost_products", "gained_products")}
        for name, outcome in sorted(gold["deletion_outcomes"].items())}
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

    print(f"{campaign.name}: {gold['accession']}, {gold['record_length']} bp, "
          f"{gold['antismash_version']} rules {gold['detection_rule_fingerprint']}")
    print(f"  {gold['regions']} regions, {gold['distinct_products']} distinct "
          f"products, {gold['contig_edge_regions']} at a record edge")
    print(f"  {gold['region_bases']} bp in regions = "
          f"{gold['fraction_of_record']:.1%} of the record; span median "
          f"{gold['region_span']['median']/1000:.1f} kb")
    print(f"  {gold['loci_total']} curated loci: {gold['verdict_census']}")
    print(f"  Jaccard median {gold['jaccard_median']}, "
          f"range {gold['jaccard_min']}-{gold['jaccard_max']}")
    print(f"  {gold['regions_without_a_locus']} of {gold['regions']} regions have "
          f"no curated counterpart; regions holding several loci: "
          f"{gold['regions_holding_several_loci']}")
    print("  deletion counterfactual:")
    for name, outcome in sorted(gold["deletion_summary"].items()):
        print(f"    {name:<24} regions={outcome['regions']:<3} "
              f"lost={outcome['regions_lost']} gained={outcome['regions_gained']} "
              f"lost_products={outcome['lost_products']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
