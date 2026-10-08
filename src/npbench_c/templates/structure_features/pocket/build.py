"""Build gold, warm bundles and the grading spec for the pocket ladder.

Usage: python -m npbench_c.templates.structure_features.pocket.build <campaign>
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
from npbench_c.templates.structure_features.pocket.analyse import analyse

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.002, "r3": 0.001, "r4": 0.001}
#: Tolerance on a reported descriptor: one unit in its last declared decimal.
DESCRIPTOR_TOLERANCE = 0.001


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if gold["numbering"]["label_to_uniprot"] == gold["numbering"]["auth_to_uniprot"]:
        raise StructureError(
            "the two structure numbering axes carry the same offset, so R1's "
            "numbering component does not distinguish a system that read the "
            "alignment record from one that guessed")
    sets = gold["descriptors"]["confirmation_sets"]
    if not sets["contact_not_annotated"]:
        raise StructureError(
            "every contact is annotated, so the shell adds nothing the annotation "
            "set does not already carry")
    if gold["descriptors"]["pocket_grid_points"] <= 0:
        raise StructureError("the pocket rule admits no lattice point")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "structure_inventory.json").write_text(json.dumps({
        k: gold[k] for k in (
            "pdb_id", "accession", "atom_inventory", "polymer_chains",
            "ligand_copies", "numbering", "resolved_residues_per_chain",
            "resolved_span_per_chain", "sequence_agreement")
    }, indent=2, sort_keys=True) + "\n")

    (warm / "r3").mkdir(parents=True)
    (warm / "r3" / "contact_shells.json").write_text(
        json.dumps(gold["shells"], indent=2, sort_keys=True) + "\n")

    # R4 starts from the shell and the descriptors, and learns nothing about the
    # cutoff sweep or which side chains matter.
    descriptors = gold["descriptors"]
    (warm / "r4").mkdir(parents=True)
    (warm / "r4" / "pocket_descriptors.json").write_text(json.dumps({
        k: descriptors[k] for k in (
            "focal_copy", "focal_ligand", "heavy_atoms", "centroid",
            "radius_of_gyration", "grid_census", "pocket_grid_points",
            "pocket_volume_a3", "contact_positions_uniprot",
            "annotated_positions_uniprot", "confirmation_sets",
            "annotation_precision", "annotation_recall")
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

    def tol(name, pred, gold_path, abs_tol=DESCRIPTOR_TOLERANCE):
        return {"name": name, "scorer": "abs_tol_match",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{gold_path}",
                         "abs_tol": abs_tol}}

    r1 = [
        {"name": "report_present", "scorer": "exact",
         "args": {"pred": "pred:present.report", "gold": 1}},
        ex("pdb_id", "pdb_id", "pdb_id"),
        ex("accession", "accession", "accession"),
        ex("atom_inventory", "atom_inventory", "atom_inventory"),
        st("polymer_chains", "polymer_chains", "polymer_chains"),
        ex("ligand_copies", "ligand_copies", "ligand_copies"),
        ex("numbering", "numbering", "numbering"),
        ex("resolved_residues_per_chain", "resolved_residues_per_chain",
           "resolved_residues_per_chain"),
        ex("resolved_span_per_chain", "resolved_span_per_chain",
           "resolved_span_per_chain"),
        ex("sequence_agreement", "sequence_agreement", "sequence_agreement"),
    ]

    r2 = [
        tol("contact_cutoff", "contact_cutoff_angstrom",
            "shells.contact_cutoff_angstrom"),
        ex("contact_counts_by_copy", "contact_counts_by_copy",
           "shells.contact_counts_by_copy"),
        ex("contacts_by_copy", "contacts_by_copy", "shells.contacts_by_copy"),
        ex("contact_positions_by_copy", "contact_positions_by_copy",
           "shells.contact_positions_by_copy"),
        ex("nearest_contact_by_copy", "nearest_contact_by_copy",
           "shells.nearest_contact_by_copy"),
        # One distance, under the declared tie-break, at one unit in the last
        # reported decimal. The full distance map is reported but not graded;
        # see excluded_from_grading.
        tol("nearest_contact_distance", "nearest_contact_distance_focal",
            "nearest_contact_distance_focal", 0.01),
        ex("chain_symmetry", "chain_symmetry", "shells.chain_symmetry"),
    ]

    r3 = [
        tol("centroid_x", "centroid_x", "descriptors.centroid.x"),
        tol("centroid_y", "centroid_y", "descriptors.centroid.y"),
        tol("centroid_z", "centroid_z", "descriptors.centroid.z"),
        tol("radius_of_gyration", "radius_of_gyration",
            "descriptors.radius_of_gyration"),
        ex("grid_census", "grid_census", "descriptors.grid_census"),
        ex("pocket_grid_points", "pocket_grid_points",
           "descriptors.pocket_grid_points"),
        # The focal shell itself is NOT graded here: R2 already grades it as
        # contact_positions_by_copy, so a component for it would score the same
        # computation twice and the warm R3 bundle would hand it over. What R3
        # adds is the reconciliation against the annotation set.
        st("annotated_positions_uniprot", "annotated_positions_uniprot",
           "descriptors.annotated_positions_uniprot"),
        ex("confirmation_sets", "confirmation_sets", "descriptors.confirmation_sets"),
    ]

    r4 = [
        ex("cutoff_sweep_counts", "cutoff_sweep_counts",
           "counterfactuals.cutoff_sweep_counts"),
        ex("cutoff_sweep_positions", "cutoff_sweep_positions",
           "counterfactuals.cutoff_sweep_positions"),
        ex("truncation_shell_size", "truncation_shell_size",
           "counterfactuals.truncation.shell_size"),
        ex("truncation_atoms_removed", "truncation_atoms_removed",
           "counterfactuals.truncation.atoms_removed_by_position"),
        st("truncation_contact_survives", "truncation_contact_survives",
           "counterfactuals.truncation.contact_survives"),
        st("truncation_contact_lost", "truncation_contact_lost",
           "counterfactuals.truncation.contact_lost"),
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
    """The one derived value the grading spec reads by a single dotted path."""
    focal = gold["descriptors"]["focal_copy"]
    distances = gold["shells"]["nearest_contact_distance_by_copy"]
    if focal not in distances:
        raise StructureError(
            f"focal copy {focal} has no nearest contact, so the shell is empty")
    gold["nearest_contact_distance_focal"] = distances[focal]
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

    shells, descriptors = gold["shells"], gold["descriptors"]
    truncation = gold["counterfactuals"]["truncation"]
    print(f"{campaign.name}: {gold['pdb_id']} / {gold['accession']}, "
          f"{gold['atom_inventory']['polymer_atom_records']} polymer atoms, "
          f"chains {gold['polymer_chains']}, numbering {gold['numbering']}")
    print(f"  ligand copies {gold['ligand_copies']}")
    print(f"  shells at {shells['contact_cutoff_angstrom']} A: "
          f"{shells['contact_counts_by_copy']}; symmetry "
          f"{ {k: v['position_sets_identical'] for k, v in shells['chain_symmetry'].items()} }")
    print(f"  focal {descriptors['focal_copy']}: rg "
          f"{descriptors['radius_of_gyration']}, "
          f"{descriptors['pocket_grid_points']} grid points "
          f"({descriptors['pocket_volume_a3']} A^3)")
    print(f"  annotation vs geometry: {len(descriptors['annotated_positions_uniprot'])} "
          f"annotated, {len(descriptors['contact_positions_uniprot'])} in contact, "
          f"{len(descriptors['confirmation_sets']['confirmed'])} confirmed, "
          f"{len(descriptors['confirmation_sets']['contact_not_annotated'])} "
          f"contacts unannotated")
    print(f"  cutoff sweep {gold['counterfactuals']['cutoff_sweep_counts']}")
    print(f"  truncation: {len(truncation['contact_lost'])} of "
          f"{truncation['shell_size']} shell residues lose the contact, "
          f"{len(truncation['contact_survives'])} keep it")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
