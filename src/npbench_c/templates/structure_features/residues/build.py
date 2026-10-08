"""Build gold, warm bundles and the grading spec for the residues ladder.

Usage: python -m npbench_c.templates.structure_features.residues.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

from npbench_c.templates.structure_features.core import StructureError
from npbench_c.templates.structure_features.params import (
    StructureSpec,
    copy_constants,
    load_task,
)
from npbench_c.templates.structure_features.residues.analyse import analyse

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.003, "r3": 0.002, "r4": 0.001}


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    focal = gold["confirmation"]
    if not any(v for v in focal["contact_not_annotated_by_ligand"].values()):
        raise StructureError(
            "no ligand has a contact position outside its annotation set, so R3 "
            "would restate R2 and the geometry adds nothing")
    if not any(not v for v in focal["annotated_by_ligand"].values()):
        raise StructureError(
            "every ligand is annotated from the focal entry, so no row tests "
            "whether a system respects WHICH structure evidences an annotation")

    swap = gold["structure_swap"]
    deltas = swap["confirmed_delta_by_ligand"]
    if not any(d > 0 for d in deltas.values()) or not any(d < 0 for d in deltas.values()):
        raise StructureError(
            "the alternate structure does not move the confirmed counts in both "
            f"directions (deltas {deltas}), so R4 has no counterfactual content")

    totals = {name: v["admitted_total"] for name, v in gold["policy_variants"].items()}
    if len(set(totals.values())) < len(totals):
        raise StructureError(
            f"two evidence-policy variants admit the same record count ({totals}), "
            "so at least one R4 component is not discriminating")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    """Pre-rendered bundles, never whole gold files.

    Each bundle holds exactly the previous rung's deliverable, so a warm start
    removes the earlier work without handing over any part of the answer the
    rung itself is asked for.
    """
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "annotation_census.json").write_text(json.dumps({
        k: gold[k] for k in (
            "accession", "sequence_length", "functional_record_counts",
            "secondary_structure_counts", "annotated_positions",
            "ligand_vocabulary", "pdb_cross_references",
            "mutagenesis_original_residue_agreements",
            "mutagenesis_original_residue_disagreements")
    }, indent=2, sort_keys=True) + "\n")

    (warm / "r3").mkdir(parents=True)
    (warm / "r3" / "evidence_census.json").write_text(json.dumps({
        k: gold[k] for k in (
            "focal_structure", "eco_code_citation_totals", "record_class_counts",
            "binding_positions_by_ligand", "pdb_evidence_index",
            "unevidenced_records", "inferred_only_records")
    }, indent=2, sort_keys=True) + "\n")

    # R4 starts from the focal confirmation and learns nothing about the
    # alternate structure or the policy variants.
    focal = gold["confirmation"]
    (warm / "r4").mkdir(parents=True)
    (warm / "r4" / "focal_confirmation.json").write_text(json.dumps({
        "structure": focal["structure"],
        "chain": focal["chain"],
        "numbering": focal["numbering"],
        "annotated_by_ligand": focal["annotated_by_ligand"],
        "contact_by_ligand": focal["contact_by_ligand"],
        "confirmed_count_by_ligand": focal["confirmed_count_by_ligand"],
        "contact_not_annotated_by_ligand": focal["contact_not_annotated_by_ligand"],
        "totals": focal["totals"],
        "active_site_position": focal["active_site_position"],
        "active_site_residue": focal["active_site_residue"],
        "active_site_contacts": focal["active_site_contacts"],
    }, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = load_task(campaign)
    grading = task.get("grading") or {}
    chance = {**DEFAULT_CHANCE_LEVELS, **(grading.get("chance_levels") or {})}
    notes = grading.get("chance_level_notes") or {}

    def ex(name, pred, gold_path):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{gold_path}"}}

    def st(name, pred, gold_path):
        return {"name": name, "scorer": "set_exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{gold_path}"}}

    r1 = [
        {"name": "report_present", "scorer": "exact",
         "args": {"pred": "pred:present.report", "gold": 1}},
        ex("accession", "accession", "accession"),
        ex("sequence_length", "sequence_length", "sequence_length"),
        ex("functional_record_counts", "functional_record_counts",
           "functional_record_counts"),
        ex("secondary_structure_counts", "secondary_structure_counts",
           "secondary_structure_counts"),
        st("annotated_positions", "annotated_positions", "annotated_positions"),
        ex("ligand_vocabulary", "ligand_vocabulary", "ligand_vocabulary"),
        st("pdb_cross_references", "pdb_cross_references", "pdb_cross_references"),
        ex("mutagenesis_check", "mutagenesis_check", "mutagenesis_check"),
    ]

    r2 = [
        st("focal_ligand_copies", "focal_ligand_copies",
           "focal_structure.ligand_copies"),
        ex("focal_polymer_atoms", "focal_polymer_atom_records",
           "focal_structure.polymer_atom_records"),
        ex("eco_code_citation_totals", "eco_code_citation_totals",
           "eco_code_citation_totals"),
        ex("record_class_counts", "record_class_counts", "record_class_counts"),
        ex("binding_positions_by_ligand", "binding_positions_by_ligand",
           "binding_positions_by_ligand"),
        ex("pdb_evidence_index", "pdb_evidence_index", "pdb_evidence_index"),
        st("unevidenced_records", "unevidenced_records", "unevidenced_records"),
        st("inferred_only_records", "inferred_only_records", "inferred_only_records"),
    ]

    r3 = [
        ex("numbering", "numbering", "confirmation.numbering"),
        ex("annotated_citing_focal", "annotated_citing_focal",
           "confirmation.annotated_by_ligand"),
        ex("contact_positions_focal", "contact_positions_focal",
           "confirmation.contact_by_ligand"),
        ex("confirmed_focal", "confirmed_focal",
           "confirmation.confirmed_count_by_ligand"),
        ex("contact_not_annotated_focal", "contact_not_annotated_focal",
           "confirmation.contact_not_annotated_by_ligand"),
        ex("confirmation_totals", "confirmation_totals", "confirmation.totals"),
        ex("active_site_position", "active_site_position",
           "confirmation.active_site_position"),
        ex("active_site_residue", "active_site_residue",
           "confirmation.active_site_residue"),
        st("active_site_contacts", "active_site_contacts",
           "confirmation.active_site_contacts"),
    ]

    r4 = [
        ex("swap_structure", "swap_structure", "structure_swap.structure"),
        # The alternate entry's numbering offsets are NOT graded: both entries of
        # this protein carry the same offsets, so the component would restate
        # R3's and the warm R4 bundle would hand it over. It stays in the report
        # because a system should still show it read the second alignment record.
        ex("swap_contact_positions", "swap_contact_positions",
           "structure_swap.contact_by_ligand"),
        ex("swap_confirmed", "swap_confirmed",
           "structure_swap.confirmed_count_by_ligand"),
        ex("swap_confirmed_delta", "swap_confirmed_delta",
           "structure_swap.confirmed_delta_by_ligand"),
        ex("policy_admitted_totals", "policy_admitted_totals",
           "policy_admitted_totals"),
        ex("policy_admitted_counts", "policy_admitted_counts",
           "policy_admitted_counts"),
        ex("policy_admitted_ligands", "policy_admitted_ligands",
           "policy_admitted_ligands"),
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


def _flatten_for_grading(gold: dict) -> dict:
    """The derived views the grading spec reads by a single dotted path."""
    gold["mutagenesis_check"] = {
        "agreements": gold["mutagenesis_original_residue_agreements"],
        "disagreements": gold["mutagenesis_original_residue_disagreements"],
    }
    variants = gold["policy_variants"]
    gold["policy_admitted_totals"] = {k: v["admitted_total"]
                                      for k, v in sorted(variants.items())}
    gold["policy_admitted_counts"] = {k: v["admitted_record_counts"]
                                      for k, v in sorted(variants.items())}
    gold["policy_admitted_ligands"] = {
        k: v["admitted_binding_positions_by_ligand"]
        for k, v in sorted(variants.items())}
    return gold


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)
    StructureSpec.load(campaign).write_reference(campaign)

    gold = _flatten_for_grading(analyse(campaign))
    _assert_has_content(gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    focal = gold["confirmation"]
    print(f"{campaign.name}: {gold['accession']} "
          f"({gold['sequence_length']} aa), "
          f"{sum(gold['functional_record_counts'].values())} functional records "
          f"over {len(gold['annotated_positions'])} positions")
    print(f"  evidence classes: "
          f"{ {t: {k: v for k, v in c.items() if v} for t, c in gold['record_class_counts'].items()} }")
    print(f"  focal {focal['structure']} chain {focal['chain']} "
          f"numbering {focal['numbering']}")
    print("  ligand                       comp  annot  contact  confirmed  extra")
    for chebi, block in focal["per_ligand"].items():
        name = gold["ligand_vocabulary"].get(chebi, chebi)
        print(f"    {name:<26} {block['component_id']:<5} "
              f"{len(block['annotated_citing_this_structure']):>5} "
              f"{len(block['contact_positions']):>8} {block['confirmed']:>10} "
              f"{len(block['contact_not_annotated']):>6}")
    print(f"  totals {focal['totals']}")
    print(f"  swap to {gold['structure_swap']['structure']}: "
          f"deltas {gold['structure_swap']['confirmed_delta_by_ligand']}")
    print(f"  policy variants: {gold['policy_admitted_totals']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
