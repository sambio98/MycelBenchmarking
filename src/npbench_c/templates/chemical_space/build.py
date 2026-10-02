"""Build gold, warm bundles and the grading spec for the chemical-space ladder.

Usage: python -m npbench_c.templates.chemical_space.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

from npbench_c.templates.chemical_space.analyse import analyse, load_params
from npbench_c.templates.chemical_space.core import CONSTANTS, RULES, TABLE, SpaceError

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.003, "r3": 0.001, "r4": 0.001}

CENSUS_KEYS = ("corpus_release", "corpus_entries", "active_entries",
               "compound_records", "entries_with_compounds",
               "distinct_compound_keys", "distinct_connectivity_keys",
               "distinct_scaffold_keys", "acyclic_compounds",
               "duplicate_records", "class_vocabulary", "entries_per_class",
               "multi_class_entries")
GRANULARITY_KEYS = ("connectivity_groups_with_multiple_compound_keys",
                    "distinctions_lost_at_connectivity",
                    "cross_entry_duplicates", "records_per_class",
                    "distinct_compounds_per_class",
                    "distinct_scaffolds_per_class",
                    "scaffold_ring_count_histogram")
COMMITMENT_KEYS = ("informative_rule", "informative_records",
                   "informative_scaffold_keys", "shared_scaffolds",
                   "cross_class_scaffolds", "cross_class_informative_entries",
                   "most_shared_informative_scaffold")


def copy_constants(campaign: pathlib.Path) -> None:
    """Publish the rules into reference/. The chemistry table is built, not copied.

    Perception is expensive and version-dependent, so build_reference.py writes
    the table once with RDKit and it is thereafter campaign data; re-deriving it
    on every gold build would let a toolkit upgrade move a published score.
    """
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)
    if not (reference / TABLE).is_file():
        raise SpaceError(
            f"reference/{TABLE} is missing. Run "
            "`python -m npbench_c.templates.chemical_space.build_reference "
            "<campaign>` once, with RDKit available, before building gold.")


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if gold["distinctions_lost_at_connectivity"] <= 0:
        raise SpaceError(
            "the connectivity key merges nothing, so the two granularities are "
            "the same statistic and R2 has no content")
    cross = gold["cross_class_scaffolds"]
    if cross["naive_all"] <= cross["all"]:
        raise SpaceError(
            f"the naive cross-class rule does not exceed the declared one "
            f"({cross['naive_all']} vs {cross['all']}), so the distinction the "
            "campaign is built on does not bite on this corpus")
    shared = gold["shared_scaffolds"]
    if shared["informative"] >= shared["all"]:
        raise SpaceError(
            f"the informativeness filter removes no shared scaffold "
            f"({shared['informative']} of {shared['all']}), so R3's "
            "reconciliation restates R2")
    gate = gold["evidence_gate"]
    if gate["cross_class_scaffolds"]["informative"] >= cross["informative"]:
        raise SpaceError(
            "the evidence gate does not shrink the cross-class set, so R4's "
            "gate has no counterfactual content")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", CENSUS_KEYS, "corpus_census.json"),
                             ("r3", GRANULARITY_KEYS, "granularity.json"),
                             ("r4", COMMITMENT_KEYS, "sharing_claim.json")):
        (warm / rung).mkdir(parents=True)
        (warm / rung / name).write_text(
            json.dumps({k: gold[k] for k in keys}, indent=2, sort_keys=True) + "\n")


def _flatten_for_grading(gold: dict) -> dict:
    """Per-variant views the grading spec reads by a single dotted path."""
    variants = gold["threshold_variants"]
    gold["threshold_variant_records"] = {
        k: v["informative_records"] for k, v in sorted(variants.items())}
    gold["threshold_variant_scaffolds"] = {
        k: v["informative_scaffold_keys"] for k, v in sorted(variants.items())}
    gold["threshold_variant_sharing"] = {
        k: {"shared_informative": v["shared_informative"],
            "cross_class_informative": v["cross_class_informative"],
            "naive_cross_class_informative": v["naive_cross_class_informative"]}
        for k, v in sorted(variants.items())}
    gate = gold["evidence_gate"]
    gold["evidence_gate_counts"] = {
        k: gate[k] for k in ("entries", "compound_records",
                             "distinct_compound_keys", "distinct_scaffold_keys",
                             "informative_records", "informative_scaffold_keys")}
    gold["evidence_gate_sharing"] = {
        "shared_scaffolds": gate["shared_scaffolds"],
        "cross_class_scaffolds": gate["cross_class_scaffolds"]}
    return gold


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = __import__("yaml").safe_load((campaign / "task.yaml").read_text())
    grading = task.get("grading") or {}
    chance = {**DEFAULT_CHANCE_LEVELS, **(grading.get("chance_levels") or {})}
    notes = grading.get("chance_level_notes") or {}

    def ex(name, path):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{path}"}}

    def st(name, path):
        return {"name": name, "scorer": "set_exact",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{path}"}}

    r1 = [
        {"name": "report_present", "scorer": "exact",
         "args": {"pred": "pred:present.report", "gold": 1}},
        # corpus_release is NOT graded: the chemistry table's header states it,
        # so a component for it would score reading one field. It stays in the
        # report as provenance; see excluded_from_grading. compound_records is
        # also in that header, but it is one of four counts in this composite and
        # the other three need the archive.
        ex("corpus_counts", "corpus_counts"),
        ex("distinct_keys", "distinct_keys"),
        ex("acyclic_compounds", "acyclic_compounds"),
        ex("duplicate_records", "duplicate_records"),
        st("class_vocabulary", "class_vocabulary"),
        ex("entries_per_class", "entries_per_class"),
        ex("multi_class_entries", "multi_class_entries"),
    ]

    r2 = [
        ex("connectivity_groups_with_multiple_compound_keys",
           "connectivity_groups_with_multiple_compound_keys"),
        ex("distinctions_lost_at_connectivity", "distinctions_lost_at_connectivity"),
        ex("cross_entry_duplicates", "cross_entry_duplicates"),
        ex("records_per_class", "records_per_class"),
        ex("distinct_compounds_per_class", "distinct_compounds_per_class"),
        ex("distinct_scaffolds_per_class", "distinct_scaffolds_per_class"),
        ex("scaffold_ring_count_histogram", "scaffold_ring_count_histogram"),
    ]

    r3 = [
        ex("informative_rule", "informative_rule"),
        ex("informative_records", "informative_records"),
        ex("informative_scaffold_keys", "informative_scaffold_keys"),
        ex("shared_scaffolds", "shared_scaffolds"),
        ex("cross_class_scaffolds", "cross_class_scaffolds"),
        ex("cross_class_informative_entries", "cross_class_informative_entries"),
        ex("most_shared_informative_scaffold", "most_shared_informative_scaffold"),
    ]

    r4 = [
        ex("threshold_variant_records", "threshold_variant_records"),
        ex("threshold_variant_scaffolds", "threshold_variant_scaffolds"),
        ex("threshold_variant_sharing", "threshold_variant_sharing"),
        ex("evidence_gate_counts", "evidence_gate_counts"),
        ex("evidence_gate_sharing", "evidence_gate_sharing"),
        st("evidence_gate_cross_class_keys",
           "evidence_gate.cross_class_informative_keys"),
    ]

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
    # Two grouped views so R1 grades four corpus counts and three key counts as
    # two claims rather than seven, keeping the conjunctive product honest.
    gold["corpus_counts"] = {k: gold[k] for k in (
        "corpus_entries", "active_entries", "compound_records",
        "entries_with_compounds")}
    gold["distinct_keys"] = {
        "compound": gold["distinct_compound_keys"],
        "connectivity": gold["distinct_connectivity_keys"],
        "scaffold": gold["distinct_scaffold_keys"]}
    gold = _flatten_for_grading(gold)
    _assert_has_content(gold)

    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    params = load_params(campaign)
    print(f"{campaign.name}: MIBiG {gold['corpus_release']}, "
          f"{gold['active_entries']} active entries, "
          f"{gold['compound_records']} compound records from "
          f"{gold['entries_with_compounds']} entries")
    print(f"  how many compounds is that? "
          f"{gold['distinct_compound_keys']} full InChIKeys, "
          f"{gold['distinct_connectivity_keys']} connectivity blocks "
          f"({gold['distinctions_lost_at_connectivity']} distinctions merged by "
          f"{gold['connectivity_groups_with_multiple_compound_keys']} groups), "
          f"{gold['distinct_scaffold_keys']} scaffolds, "
          f"{gold['acyclic_compounds']} acyclic")
    print(f"  cross-entry duplicates {gold['cross_entry_duplicates']}")
    print(f"  shared scaffolds {gold['shared_scaffolds']}; "
          f"cross-class {gold['cross_class_scaffolds']}")
    print(f"  most shared informative scaffold: "
          f"{gold['most_shared_informative_scaffold']}")
    print("  threshold variants (rings, coverage -> records / shared / cross-class):")
    for name, v in sorted(gold["threshold_variants"].items()):
        print(f"    {name:<18} ({v['min_ring_count']}, "
              f"{v['min_heavy_atom_coverage']}) -> {v['informative_records']:>5} / "
              f"{v['shared_informative']:>4} / {v['cross_class_informative']:>3}")
    gate = gold["evidence_gate"]
    print(f"  evidence gate: {gate['entries']} entries, "
          f"{gate['compound_records']} records, cross-class informative "
          f"{gate['cross_class_scaffolds']['informative']} "
          f"({len(gate['cross_class_informative_keys'])} keys)")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          f"components; {len(params['threshold_variants'])} declared variants")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
