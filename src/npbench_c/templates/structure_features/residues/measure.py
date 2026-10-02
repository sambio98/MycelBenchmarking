"""Measure a residues-ladder submission. Deterministic, agent-independent.

The measurer restructures the submission into the flat per-ligand views the
grading spec compares, so the agent's contract can stay in the readable nested
shape while the spec keeps one component per claim instead of one per
(ligand, field) pair. Nothing here reads gold: the keys it iterates come from
the campaign's own declared vocabularies.

Usage:
  python -m npbench_c.templates.structure_features.residues.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.structure_features.core import Rules
from npbench_c.templates.structure_features.params import StructureSpec

REPORT = "residue_report.json"


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over whatever keys the submission used, zeros dropped.

    A submission that lists every code with a zero compares equal to gold,
    which omits the zeros: presentation is not graded.
    """
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def _positions(value: object) -> list[int]:
    if not isinstance(value, list):
        return []
    return sorted({n for v in value if (n := _int(v)) is not None})


def _position_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _positions(block.get(k)) for k in keys}


def _strings(value: object) -> list[str]:
    return sorted({str(v) for v in value}) if isinstance(value, list) else []


def _str_map(value: object) -> dict:
    block = value if isinstance(value, dict) else {}
    return {str(k): str(v) for k, v in sorted(block.items())}


def _nested(value: object, key: str) -> object:
    return value.get(key) if isinstance(value, dict) else None


def _per_ligand(block: object, field: str, ligands, as_int: bool = False) -> dict:
    """Pull one field out of a submission's per_ligand map, for every ligand."""
    per = _nested(block, "per_ligand")
    per = per if isinstance(per, dict) else {}
    out = {}
    for chebi in ligands:
        row = per.get(chebi)
        row = row if isinstance(row, dict) else {}
        out[chebi] = _int(row.get(field)) if as_int else _positions(row.get(field))
    return out


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    spec = StructureSpec.load(campaign)
    types = tuple(sorted(rules.functional_types))
    secondary = tuple(sorted(rules.secondary_types))
    ligands = tuple(sorted(spec.ligand_components))
    variants = tuple(sorted(spec.policy_variants))

    path = pathlib.Path(submission_dir) / REPORT
    out: dict = {"present": {"report": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["report"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["report"] = 0
        return out

    # ---- R1: the annotation census
    out["accession"] = str(sub.get("accession")) if sub.get("accession") is not None else None
    out["sequence_length"] = _int(sub.get("sequence_length"))
    out["functional_record_counts"] = _int_map(sub.get("functional_record_counts"), types)
    out["secondary_structure_counts"] = _int_map(
        sub.get("secondary_structure_counts"), secondary)
    out["annotated_positions"] = _positions(sub.get("annotated_positions"))
    out["ligand_vocabulary"] = _str_map(sub.get("ligand_vocabulary"))
    out["pdb_cross_references"] = _strings(sub.get("pdb_cross_references"))
    out["mutagenesis_check"] = {
        "agreements": _int(sub.get("mutagenesis_original_residue_agreements")),
        "disagreements": _int(sub.get("mutagenesis_original_residue_disagreements")),
    }

    # ---- R2: the evidence census and the structure inventory
    focal_structure = sub.get("focal_structure")
    out["focal_ligand_copies"] = _strings(_nested(focal_structure, "ligand_copies"))
    out["focal_polymer_atom_records"] = _int(
        _nested(focal_structure, "polymer_atom_records"))
    out["eco_code_citation_totals"] = _free_int_map(sub.get("eco_code_citation_totals"))
    classes = sub.get("record_class_counts")
    classes = classes if isinstance(classes, dict) else {}
    out["record_class_counts"] = {
        t: _int_map(classes.get(t), sorted(rules.class_vocabulary)) for t in types}
    out["binding_positions_by_ligand"] = _position_map(
        sub.get("binding_positions_by_ligand"), ligands)
    pdb_index = sub.get("pdb_evidence_index")
    pdb_index = pdb_index if isinstance(pdb_index, dict) else {}
    out["pdb_evidence_index"] = {str(k).upper(): _positions(v)
                                 for k, v in sorted(pdb_index.items())
                                 if _positions(v)}
    out["unevidenced_records"] = _strings(sub.get("unevidenced_records"))
    out["inferred_only_records"] = _strings(sub.get("inferred_only_records"))

    # ---- R3: confirmation against the focal structure
    confirmation = sub.get("confirmation")
    out["numbering"] = _int_map(_nested(confirmation, "numbering"),
                                ("auth_to_uniprot", "label_to_uniprot"))
    out["annotated_citing_focal"] = _per_ligand(
        confirmation, "annotated_citing_this_structure", ligands)
    out["contact_positions_focal"] = _per_ligand(
        confirmation, "contact_positions", ligands)
    out["confirmed_focal"] = _per_ligand(confirmation, "confirmed", ligands, as_int=True)
    out["contact_not_annotated_focal"] = _per_ligand(
        confirmation, "contact_not_annotated", ligands)
    out["confirmation_totals"] = _int_map(
        _nested(confirmation, "totals"),
        ("annotated_citing_this_structure", "confirmed", "unconfirmed"))
    out["active_site_position"] = _int(_nested(confirmation, "active_site_position"))
    residue = _nested(confirmation, "active_site_residue")
    out["active_site_residue"] = str(residue) if isinstance(residue, str) else None
    out["active_site_contacts"] = _strings(_nested(confirmation, "active_site_contacts"))

    # ---- R4: the structure swap and the evidence-policy variants
    swap = sub.get("structure_swap")
    structure = _nested(swap, "structure")
    out["swap_structure"] = str(structure).upper() if isinstance(structure, str) else None
    out["swap_numbering"] = _int_map(_nested(swap, "numbering"),
                                     ("auth_to_uniprot", "label_to_uniprot"))
    out["swap_contact_positions"] = _per_ligand(swap, "contact_positions", ligands)
    out["swap_confirmed"] = _per_ligand(swap, "confirmed", ligands, as_int=True)
    out["swap_confirmed_delta"] = _int_map(
        _nested(swap, "confirmed_delta_by_ligand"), ligands)

    policies = sub.get("policy_variants")
    policies = policies if isinstance(policies, dict) else {}
    totals, counts, per_ligand_positions = {}, {}, {}
    for name in variants:
        block = policies.get(name)
        block = block if isinstance(block, dict) else {}
        totals[name] = _int(block.get("admitted_total"))
        counts[name] = _int_map(block.get("admitted_record_counts"), types)
        # Gold omits a ligand with no admitted position, so an empty list from
        # the submission is dropped rather than compared against a missing key.
        positions = _position_map(
            block.get("admitted_binding_positions_by_ligand"), ligands)
        per_ligand_positions[name] = {k: v for k, v in positions.items() if v}
    out["policy_admitted_totals"] = totals
    out["policy_admitted_counts"] = counts
    out["policy_admitted_ligands"] = per_ligand_positions
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>",
              file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
