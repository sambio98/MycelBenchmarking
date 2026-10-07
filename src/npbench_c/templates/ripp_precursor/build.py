"""Build gold, warm bundles and the grading spec for the RiPP-precursor ladder.

Usage: python -m npbench_c.templates.ripp_precursor.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.ripp_precursor.analyse import analyse, load_params
from npbench_c.templates.ripp_precursor.core import (
    CONSTANTS,
    RULES,
    PrecursorError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("mibig_release", "biosynthetic_class", "baseline_policy",
                  "threshold", "panel_size", "hmmer_version", "pfam_version",
                  "pfam_sha256", "selection_census",
                  "core_sequence_shape_census", "core_sequence_case_census",
                  "records_declaring_a_cleavage_location", "ripp_type_census",
                  "status_census", "precursor_length_summary")
#: R2's answers: where each declared core sits in its declared precursor.
LOCALISATION_KEYS = ("localisation_table", "localisation_census",
                     "leader_length_summary")
#: R3's answers: the coordinate convention, and what the domain models say.
EVIDENCE_KEYS = ("convention_tally", "elected_convention", "cleavage_table",
                 "domain_table", "panel_census", "records_with_a_panel_hit",
                 "boundary_relation_census", "boundary_offset_summary")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("policy_sweep", "deletion_table", "deletion_census")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    census = gold["selection_census"]
    if census["corpus_records"] < 20:
        raise PrecursorError(
            f"{census['corpus_records']} precursor records resolved; an audit "
            "of that many annotations is an anecdote")
    if census["corpus_entries"] < 10:
        raise PrecursorError(
            f"the corpus spans {census['corpus_entries']} clusters, so the "
            "audit is about a handful of submissions rather than a database")

    # The selection rule has three clauses and each one has to remove something,
    # or the campaign inherited it from another instantiation rather than
    # exercising it here.
    for field, complaint in (
            ("entries_excluded_by_status",
             "no entry is excluded by status, so that clause of the selection "
             "rule is not exercised"),
            ("entries_absent_from_reference_set",
             "every selected entry has a reference record, so the campaign "
             "never has to notice that the two resources cover different sets")):
        if not census[field]:
            raise PrecursorError(complaint)

    shapes = gold["core_sequence_shape_census"]
    if len(shapes) < 3:
        raise PrecursorError(
            f"core_sequence arrives in {len(shapes)} shapes ({sorted(shapes)}); "
            "a declared ingestion policy is pointless where the field is "
            "consistent")
    if len(gold["core_sequence_case_census"]) < 2:
        raise PrecursorError(
            "every core is in the same case, so case folding is a transform "
            "with nothing to transform and the baseline policy is arbitrary")

    localisation = gold["localisation_census"]
    populated = sorted(v for v, n in localisation.items() if n)
    if len(populated) < 4:
        raise PrecursorError(
            f"only {populated} localisation verdicts occur; a vocabulary that "
            "collapses to two or three cells is an enum dressed up")
    located = sum(1 for row in gold["localisation_table"]
                  if row["core_start"] is not None)
    if located == len(gold["localisation_table"]):
        raise PrecursorError(
            "every declared core is cleanly located, so the verdict vocabulary "
            "never has to report a coordinate it cannot derive")
    if not located:
        raise PrecursorError("no core is located at all, which is a broken read")

    # The convention is elected, so the election has to be decisive on this data
    # and the losing readings have to actually lose. Otherwise 'elected' is a
    # tie-break dressed as a measurement.
    tally = gold["convention_tally"]
    elected = gold["elected_convention"]
    agreements = sorted((block["agrees"] for block in tally.values()),
                        reverse=True)
    if agreements[0] == 0:
        raise PrecursorError(
            "no declared convention agrees with the measured core position on "
            "any record, so the coordinate field cannot be reconciled at all")
    if len(agreements) > 1 and agreements[0] == agreements[1]:
        raise PrecursorError(
            f"two conventions agree on {agreements[0]} records each, so "
            f"{elected} is the declared tie-break rather than a finding")
    if not any(block["disagrees"] for block in tally.values()):
        raise PrecursorError(
            "no convention disagrees with the measured position on any record, "
            "so the reconciliation cannot tell the readings apart")

    relations = gold["boundary_relation_census"]
    if not gold["records_with_a_panel_hit"]:
        raise PrecursorError(
            "the declared panel matches nothing in the corpus, so R3's domain "
            "evidence is empty and R4's deletion has nothing to remove")
    if len([r for r, n in relations.items() if n]) < 2:
        raise PrecursorError(
            "every panel hit sits the same way across the leader/core boundary, "
            "so the relation vocabulary is not distinguishing anything")

    sweep = gold["policy_sweep"]
    if len(sweep) < 3:
        raise PrecursorError(
            f"a sweep over {len(sweep)} ingestion policies cannot show a "
            "dependence")
    censuses = {json.dumps(block["verdict_census"], sort_keys=True)
                for block in sweep.values()}
    if len(censuses) < len(sweep):
        raise PrecursorError(
            f"{len(sweep)} declared policies produce {len(censuses)} distinct "
            "censuses, so at least two of them are the same policy under two "
            "names on this corpus")

    deletion = gold["deletion_census"]
    if not deletion["unchanged_absent"]:
        raise PrecursorError(
            "every record in the corpus carries a panel hit, so the deletion "
            "counterfactual has no negative control and 'the call is lost' is "
            "the only answer the rung can take")
    moved = deletion["lost"] + deletion["changed"] + deletion["gained"]
    if not moved:
        raise PrecursorError(
            "removing the declared leader changes no domain call, so R4 shows "
            "nothing about which side of the boundary the signal lives on")
    if not deletion["retained"]:
        raise PrecursorError(
            "every panel call is lost with its leader, so 'the signal is in the "
            "leader' is a free answer and the rung distinguishes no record")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", LOCALISATION_KEYS, "localisation.json"),
                             ("r4", EVIDENCE_KEYS, "evidence.json")):
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

    # The echoes of the rules file -- biosynthetic_class, baseline_policy,
    # threshold, panel_size and mibig_release -- are reported and NOT graded.
    # Their value is fixed by the campaign's own declaration, so grading them
    # adds weight to a conjunctive product and no information.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("selection_census"), ex("core_sequence_shape_census"),
          ex("core_sequence_case_census"),
          ex("records_declaring_a_cleavage_location"), ex("ripp_type_census"),
          ex("status_census"), ex("precursor_length_summary"), ex("environment")]

    r2 = [ex("localisation_table"), ex("localisation_census"),
          ex("leader_length_summary")]

    r3 = [ex("convention_tally"), ex("elected_convention"), ex("cleavage_table"),
          ex("domain_table"), ex("panel_census"), ex("records_with_a_panel_hit"),
          ex("boundary_relation_census"), ex("boundary_offset_summary")]

    r4 = [ex("policy_sweep"), ex("deletion_table"), ex("deletion_census")]

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

    The tool identity is one claim -- which HMMER and which Pfam produced these
    domain calls -- and grading the three fields apart would let a system that
    pinned none of them lose one component instead of the claim.
    """
    gold["environment"] = {"hmmer_version": gold["hmmer_version"],
                           "pfam_version": gold["pfam_version"],
                           "pfam_sha256": gold["pfam_sha256"]}
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

    census = gold["selection_census"]
    print(f"{campaign.name}: MIBiG {gold['mibig_release']}, "
          f"{census['entries_in_class']} {gold['biosynthetic_class']} entries")
    print(f"  {census['entries_with_a_named_precursor']} name a precursor gene "
          f"and a core; {census['entries_excluded_by_status']} excluded by "
          f"status, {census['entries_absent_from_reference_set']} absent from "
          f"the reference set")
    print(f"  corpus: {census['corpus_records']} records over "
          f"{census['corpus_entries']} clusters, "
          f"{gold['records_declaring_a_cleavage_location']} declaring a cleavage "
          f"location")
    print(f"  core_sequence shapes: {gold['core_sequence_shape_census']}")
    print(f"  case: {gold['core_sequence_case_census']}")
    print(f"  localisation under {gold['baseline_policy']}: "
          f"{ {k: v for k, v in gold['localisation_census'].items() if v} }")
    print(f"  leader lengths: {gold['leader_length_summary']}")
    print("  cleavage conventions:")
    for name, block in gold["convention_tally"].items():
        mark = " <- elected" if name == gold["elected_convention"] else ""
        print(f"    {name:<28} {block}{mark}")
    print(f"  {gold['hmmer_version']} against {gold['pfam_version']}, "
          f"panel of {gold['panel_size']} at {gold['threshold']}")
    print(f"  {gold['records_with_a_panel_hit']} of {census['corpus_records']} "
          f"records carry a panel hit: "
          f"{ {k: v for k, v in gold['panel_census'].items() if v} }")
    print(f"  boundary relation: "
          f"{ {k: v for k, v in gold['boundary_relation_census'].items() if v} }")
    print(f"  boundary offset: {gold['boundary_offset_summary']}")
    print("  ingestion-policy sweep:")
    for name, block in sorted(gold["policy_sweep"].items()):
        print(f"    {name:<24} located={block['located']:<3} "
              f"{ {k: v for k, v in block['verdict_census'].items() if v} }")
    print(f"  leader removed: "
          f"{ {k: v for k, v in gold['deletion_census'].items() if v} }")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
