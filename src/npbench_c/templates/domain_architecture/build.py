"""Build gold, warm bundles and the grading spec for the domain-architecture ladder.

Usage: python -m npbench_c.templates.domain_architecture.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.domain_architecture.analyse import analyse, load_params
from npbench_c.templates.domain_architecture.core import (
    CONSTANTS,
    RULES,
    ArchitectureError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("mibig_release", "module_type", "baseline_cutoff",
                  "baseline_resolution", "baseline_span", "baseline_mapping",
                  "panel_size", "hmmer_version", "pfam_version", "pfam_sha256",
                  "selection_census", "curated_coordinate_usability",
                  "curated_active_flag", "curated_module_census",
                  "curated_domain_totals", "protein_length_summary",
                  "panel_census")
#: R2's answers: the architecture per gene.
ARCHITECTURE_KEYS = ("architecture_table", "architecture_census",
                     "overlap_census")
#: R3's answers: against the curated modules.
RECONCILIATION_KEYS = ("reconciliation_table", "a_domain_verdicts",
                       "c_domain_verdicts", "a_delta_census", "module_table",
                       "module_verdicts", "complete_modules_total")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("policy_sweep", "mapping_sweep")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    census = gold["selection_census"]
    if census["corpus_genes"] < 50:
        raise ArchitectureError(
            f"{census['corpus_genes']} genes survived the selection; an audit of "
            "that many architectures is an anecdote")
    if census["corpus_entries"] < 20:
        raise ArchitectureError(
            f"the corpus spans {census['corpus_entries']} clusters, so it is "
            "about a handful of submissions rather than a database")
    for field in ("drop_excluded_by_status", "drop_other_module_type"):
        if not census.get(field):
            raise ArchitectureError(
                f"{field} removes nothing, so that clause of the selection rule "
                "is not exercised here and was inherited from another "
                "instantiation")

    arch = gold["architecture_census"]
    if arch["distinct_architectures"] < 20:
        raise ArchitectureError(
            f"only {arch['distinct_architectures']} distinct architectures over "
            f"{arch['genes']} genes; the corpus is one protein repeated")
    if arch["domains_per_gene"]["max"] < 4:
        raise ArchitectureError(
            "no gene carries more than three domains, so nothing here is a "
            "multidomain protein and the architecture is a domain name")

    # The resolution policy exists because the raw hits overlap. If they did not,
    # the policy would be ceremony and R4's first axis would be a no-op.
    overlaps = gold["overlap_census"]
    if not overlaps["overlapping_pairs"]:
        raise ArchitectureError(
            "no two hits overlap anywhere in the corpus, so the resolution "
            "policy resolves nothing and the campaign's central complication is "
            "absent")
    if not any(row["resolved_away"] for row in gold["architecture_table"]):
        raise ArchitectureError(
            "the resolution policy drops no hit, so the baseline and the "
            "unresolved architecture are the same answer")

    # R3 needs real disagreement in both directions of its vocabulary, or the
    # reconciliation is a restatement of one side.
    a = gold["a_domain_verdicts"]
    if not a["agrees"]:
        raise ArchitectureError(
            "the computed adenylation count agrees with the curated one on no "
            "gene at all, which is a broken mapping rather than a finding")
    disagreeing = a["computed_more"] + a["computed_fewer"]
    if not disagreeing:
        raise ArchitectureError(
            "the computed and curated adenylation counts agree on every gene, so "
            "the reconciliation has nothing to grade and the curated record is "
            "recomputable from sequence")
    if len(gold["a_delta_census"]) < 3:
        raise ArchitectureError(
            f"the adenylation delta takes {len(gold['a_delta_census'])} values, "
            "so the disagreement is a single offset rather than a distribution")
    if not gold["c_domain_verdicts"]["not_curated"]:
        raise ArchitectureError(
            "every gene curates a condensation domain, so the not_curated "
            "outcome is never exercised")
    if not gold["module_verdicts"]["agrees"]:
        raise ArchitectureError(
            "the declared module pattern is found on no gene, so the "
            "decomposition is empty everywhere")

    # R4: each axis has to move something, measured here rather than assumed.
    sweep = gold["policy_sweep"]
    if len(sweep) < 8:
        raise ArchitectureError(
            f"a policy sweep over {len(sweep)} cells cannot show three axes")
    moved = [cell["genes_differing_from_baseline"] for cell in sweep.values()]
    if max(moved) == 0:
        raise ArchitectureError(
            "no declared reading changes any architecture, so R4 has no axis and "
            "the three choices the campaign declares are all cosmetic")
    distinct = {(cell["domains_total"], cell["distinct_architectures"],
                 cell["a_domain_agreement"]) for cell in sweep.values()}
    if len(distinct) < 3:
        raise ArchitectureError(
            f"the sweep's {len(sweep)} cells take {len(distinct)} distinct "
            "shapes, so the axes are not independent on this corpus")

    mapping = gold["mapping_sweep"]
    if len(mapping) < 3:
        raise ArchitectureError(
            f"a mapping sweep over {len(mapping)} readings cannot show that the "
            "family-to-domain mapping is load-bearing")
    agreements = {name: block["verdicts"]["agrees"]
                  for name, block in mapping.items()}
    if len(set(agreements.values())) < len(agreements):
        raise ArchitectureError(
            f"two declared mappings give the same agreement ({agreements}), so "
            "at least two of them are the same mapping under two names")
    if max(agreements.values()) < 2 * min(agreements.values()):
        raise ArchitectureError(
            f"the mapping choice moves agreement by less than a factor of two "
            f"({agreements}), so declaring it is not worth a rung")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", ARCHITECTURE_KEYS, "architecture.json"),
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

    # The echoes of the rules file and the input provenance are reported and NOT
    # graded: their value is fixed by the campaign's own declaration.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("selection_census"), ex("curated_coordinate_usability"),
          ex("curated_active_flag"), ex("curated_module_census"),
          ex("curated_domain_totals"), ex("protein_length_summary"),
          ex("panel_census"), ex("environment")]

    r2 = [ex("architecture_table"), ex("architecture_census"),
          ex("overlap_census")]

    r3 = [ex("reconciliation_table"), ex("a_domain_verdicts"),
          ex("c_domain_verdicts"), ex("a_delta_census"), ex("module_table"),
          ex("module_verdicts"), ex("complete_modules_total")]

    r4 = [ex("policy_sweep"), ex("mapping_sweep")]

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
    """The tool identity is one claim -- which HMMER and which Pfam produced these
    domain calls -- so grading the three fields apart would let a system that
    pinned none of them lose one component instead of the claim."""
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
    arch = gold["architecture_census"]
    print(f"{campaign.name}: MIBiG {gold['mibig_release']}, "
          f"{census['modules']} curated modules over "
          f"{census['entries_with_modules']} entries")
    print(f"  corpus: {census['corpus_genes']} {gold['module_type']} genes over "
          f"{census['corpus_entries']} clusters, from "
          f"{census['shipped_proteins']} shipped translations")
    print(f"  curated coordinate usability: {gold['curated_coordinate_usability']}")
    print(f"  curated `active` flag: {gold['curated_active_flag']}")
    print(f"  {gold['hmmer_version']} against {gold['pfam_version']}, panel of "
          f"{gold['panel_size']} at {gold['baseline_cutoff']}")
    print(f"  lengths: {gold['protein_length_summary']}")
    print(f"  {arch['domains_total']} domains, {arch['distinct_architectures']} "
          f"distinct architectures, per gene "
          f"{arch['domains_per_gene']}")
    print(f"  overlaps before resolution: {gold['overlap_census']}")
    print(f"  adenylation reconciliation: {gold['a_domain_verdicts']}")
    print(f"  delta distribution: {gold['a_delta_census']}")
    print(f"  condensation reconciliation: {gold['c_domain_verdicts']}")
    print(f"  complete declared modules: {gold['complete_modules_total']}, "
          f"against curated: {gold['module_verdicts']}")
    print("  policy sweep:")
    for cell, block in sorted(gold["policy_sweep"].items()):
        print(f"    {cell:<40} domains={block['domains_total']:<6} "
              f"distinct={block['distinct_architectures']:<5} "
              f"agree={block['a_domain_agreement']:<5} "
              f"moved={block['genes_differing_from_baseline']}")
    print("  mapping sweep:")
    for name, block in sorted(gold["mapping_sweep"].items()):
        print(f"    {name:<20} {'+'.join(block['families']):<30} "
              f"{block['verdicts']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
