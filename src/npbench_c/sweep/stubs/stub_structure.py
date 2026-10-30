"""Deterministic stub systems for the two structure-features ladders.

Harness fixtures, never a published result: a stub-based sweep marks itself and
the readiness gate keeps every agent-run check PENDING.

Each ladder gets one stub per rung boundary, one competence probe and one
leakage ablation, following the separation the MIBiG diff campaign forced:

  residues   parser      clears R1: the annotation census
             evidence    clears R2: ECO classes and the PDB citation index
             confirm     clears R3: annotation reconciled against geometry
             complete    clears R4
             annotated   COMPETENCE probe, tooled: reads both sources but takes
                         the annotation set to BE the contact set, so it clears
                         R2 and fails R3
             noncompute  LEAKAGE ablation: never opens inputs/

  pocket     parser      clears R1: inventory and numbering offsets
             shells      clears R2: the contact shells
             descriptors clears R3: pocket descriptors and reconciliation
             complete    clears R4
             labelseq    COMPETENCE probe, tooled: reads the alignment record,
                         reports the offsets correctly, then applies the ENTITY
                         offset when mapping contacts to sequence positions, so
                         it clears R1 and fails R2
             noncompute  LEAKAGE ablation: never opens inputs/
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.structure_features.core import Rules, Structure, Annotations
from npbench_c.templates.structure_features.params import StructureSpec
from npbench_c.templates.structure_features.pocket import analyse as pocket_analyse
from npbench_c.templates.structure_features.residues import analyse as residues_analyse

LADDERS = ("residues", "pocket")
RESIDUE_LEVELS = ("parser", "evidence", "confirm", "complete", "annotated",
                  "noncompute")
POCKET_LEVELS = ("parser", "shells", "descriptors", "complete", "labelseq",
                 "noncompute")


def _write(path: pathlib.Path, payload: dict) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


# ------------------------------------------------------------------ residues


def _residues(sandbox: pathlib.Path, level: str, out: pathlib.Path) -> int:
    target = out / "residue_report.json"
    if level == "noncompute":
        # Shaped like a report, informed by nothing. Record counts and position
        # lists cannot be guessed, so this fails R1 and sets the leakage floor.
        return _write(target, {
            "accession": "P00000", "sequence_length": 500,
            "functional_record_counts": {"Active site": 1, "Binding site": 10,
                                         "Site": 1, "Mutagenesis": 2},
            "secondary_structure_counts": {}, "annotated_positions": [],
            "ligand_vocabulary": {}, "pdb_cross_references": [],
            "mutagenesis_original_residue_agreements": 2,
            "mutagenesis_original_residue_disagreements": 0,
            "eco_code_citation_totals": {}, "record_class_counts": {},
            "binding_positions_by_ligand": {}, "pdb_evidence_index": {},
            "unevidenced_records": [], "inferred_only_records": [],
            "focal_structure": {}, "confirmation": {}, "structure_swap": {},
            "policy_variants": {},
        })

    rules = Rules.load(sandbox)
    spec = StructureSpec.load(sandbox)
    annotations = Annotations.load(sandbox / spec.annotation_input, rules)
    report = residues_analyse.census(annotations, rules)
    if level == "parser":
        return _write(target, report)

    focal = Structure.load(sandbox / spec.focal_input, spec.focal_pdb_id)
    report.update(residues_analyse.evidence_census(annotations, rules))
    report["focal_structure"] = {
        "pdb_id": focal.pdb_id,
        "polymer_atom_records": len(focal.polymer()),
        "ligand_copies": sorted(f"{comp}:{ch}:{res}"
                                for comp, ch, res in focal.ligand_copies()),
    }
    if level == "evidence":
        return _write(target, report)

    confirmation = residues_analyse.confirm(
        annotations, focal, spec.focal_chain, spec, rules)
    if level == "annotated":
        # The error the rung exists to catch: "annotated" read as "confirmed",
        # with the contact shell never computed.
        collapsed = {}
        for chebi, block in confirmation["per_ligand"].items():
            annotated = block["annotated_citing_this_structure"]
            collapsed[chebi] = {
                "component_id": block["component_id"],
                "copies_in_chain": block["copies_in_chain"],
                "annotated_citing_this_structure": annotated,
                "contact_positions": annotated,
                "confirmed": len(annotated),
                "annotated_not_in_contact": [],
                "contact_not_annotated": [],
            }
        report["confirmation"] = {
            "structure": confirmation["structure"],
            "chain": confirmation["chain"],
            "numbering": confirmation["numbering"],
            "per_ligand": collapsed,
            "totals": confirmation["totals"],
            "active_site_position": confirmation["active_site_position"],
            "active_site_residue": confirmation["active_site_residue"],
            "active_site_contacts": [],
        }
        return _write(target, report)

    report["confirmation"] = confirmation
    if level == "confirm":
        return _write(target, report)

    alternate = Structure.load(sandbox / spec.alternate_input, spec.alternate_pdb_id)
    swap = residues_analyse.confirm(
        annotations, alternate, spec.focal_chain, spec, rules)
    swap["confirmed_delta_by_ligand"] = {
        chebi: (swap["per_ligand"][chebi]["confirmed"]
                - confirmation["per_ligand"][chebi]["confirmed"])
        for chebi in sorted(spec.ligand_components)}
    report["structure_swap"] = swap
    report["policy_variants"] = residues_analyse.policy_variants(
        annotations, spec, rules)
    return _write(target, report)


# -------------------------------------------------------------------- pocket


def _pocket(sandbox: pathlib.Path, level: str, out: pathlib.Path) -> int:
    target = out / "pocket_report.json"
    if level == "noncompute":
        return _write(target, {
            "pdb_id": "0XXX", "accession": "P00000",
            "atom_inventory": {"polymer_atom_records": 4000,
                               "ligand_atom_records": 20,
                               "water_atom_records": 300,
                               "excluded_hydrogen_records": 0,
                               "excluded_altloc_records": 0},
            "polymer_chains": ["A"], "ligand_copies": {},
            "numbering": {"auth_to_uniprot": 0, "label_to_uniprot": 0},
            "resolved_residues_per_chain": {}, "resolved_span_per_chain": {},
            "sequence_agreement": {"agreements": 0, "disagreements": 0},
            "shells": {}, "descriptors": {}, "counterfactuals": {},
        })

    rules = Rules.load(sandbox)
    spec = StructureSpec.load(sandbox)
    annotations = Annotations.load(sandbox / spec.annotation_input, rules)
    structure = Structure.load(sandbox / spec.focal_input, spec.focal_pdb_id)

    report = pocket_analyse.inventory(structure, annotations, rules)
    if level == "parser":
        return _write(target, report)

    shells = pocket_analyse.shells(structure, annotations, rules)
    if level == "labelseq":
        # Reads the alignment record and reports both offsets correctly, then
        # applies the ENTITY offset when turning contacts into sequence
        # positions. R1 is clean, R2 is uniformly shifted.
        offsets = report["numbering"]
        shift = offsets["label_to_uniprot"] - offsets["auth_to_uniprot"]
        shells["contact_positions_by_copy"] = {
            copy: [p + shift for p in positions]
            for copy, positions in shells["contact_positions_by_copy"].items()}
        shells["chain_symmetry"] = {
            comp: {**block,
                   "symmetric_difference": [p + shift
                                            for p in block["symmetric_difference"]]}
            for comp, block in shells["chain_symmetry"].items()}
        report["shells"] = shells
        return _write(target, report)

    report["shells"] = shells
    if level == "shells":
        return _write(target, report)

    report["descriptors"] = pocket_analyse.descriptors(
        structure, annotations, spec, rules)
    if level == "descriptors":
        return _write(target, report)

    report["counterfactuals"] = pocket_analyse.counterfactuals(
        structure, annotations, spec, rules)
    return _write(target, report)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ladder", required=True, choices=LADDERS)
    ap.add_argument("--level", required=True)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    levels = RESIDUE_LEVELS if args.ladder == "residues" else POCKET_LEVELS
    if args.level not in levels:
        ap.error(f"--level must be one of {levels} for ladder {args.ladder}")

    sandbox = pathlib.Path.cwd()
    out = pathlib.Path(args.submission)
    out.mkdir(parents=True, exist_ok=True)
    if args.ladder == "residues":
        return _residues(sandbox, args.level, out)
    return _pocket(sandbox, args.level, out)


if __name__ == "__main__":
    raise SystemExit(main())
