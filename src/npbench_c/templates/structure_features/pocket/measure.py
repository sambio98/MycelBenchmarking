"""Measure a pocket-ladder submission. Deterministic, agent-independent.

Nothing here reads gold. The keys it iterates come from the campaign's own
declared vocabularies -- the chains and ligand copies a submission reports, the
cutoff variants task.yaml declares -- so a submission cannot be rewarded for
omitting a block it found inconvenient.

Usage:
  python -m npbench_c.templates.structure_features.pocket.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.structure_features.params import StructureSpec

REPORT = "pocket_report.json"

INVENTORY_KEYS = ("polymer_atom_records", "ligand_atom_records",
                  "water_atom_records", "excluded_hydrogen_records",
                  "excluded_altloc_records")
CONFIRMATION_KEYS = ("confirmed", "annotated_not_in_contact",
                     "contact_not_annotated")
GRID_KEYS = ("grid_points_near_ligand", "grid_points_rejected_clash",
             "grid_points_rejected_bulk")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _int_map(value: object, keys=None) -> dict:
    block = value if isinstance(value, dict) else {}
    if keys is not None:
        return {k: _int(block.get(k)) for k in keys}
    return {str(k): n for k, v in sorted(block.items()) if (n := _int(v)) is not None}


def _positions(value: object) -> list[int]:
    if not isinstance(value, list):
        return []
    return sorted({n for v in value if (n := _int(v)) is not None})


def _position_map(value: object) -> dict:
    block = value if isinstance(value, dict) else {}
    return {str(k): _positions(v) for k, v in sorted(block.items())}


def _label_map(value: object) -> dict:
    block = value if isinstance(value, dict) else {}
    return {str(k): sorted({str(x) for x in v}) for k, v in sorted(block.items())
            if isinstance(v, list)}


def _str_map(value: object) -> dict:
    block = value if isinstance(value, dict) else {}
    return {str(k): str(v) for k, v in sorted(block.items()) if isinstance(v, str)}


def _nested(value: object, *keys: str) -> object:
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    spec = StructureSpec.load(pathlib.Path(campaign))

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

    # ---- R1: the structure inventory and the numbering
    pdb_id = sub.get("pdb_id")
    out["pdb_id"] = str(pdb_id).upper() if isinstance(pdb_id, str) else None
    accession = sub.get("accession")
    out["accession"] = str(accession) if isinstance(accession, str) else None
    out["atom_inventory"] = _int_map(sub.get("atom_inventory"), INVENTORY_KEYS)
    chains = sub.get("polymer_chains")
    out["polymer_chains"] = sorted({str(c) for c in chains}) if isinstance(chains, list) else []
    out["ligand_copies"] = _int_map(sub.get("ligand_copies"))
    out["numbering"] = _int_map(sub.get("numbering"),
                                ("auth_to_uniprot", "label_to_uniprot"))
    out["resolved_residues_per_chain"] = _int_map(
        sub.get("resolved_residues_per_chain"))
    out["resolved_span_per_chain"] = _position_map(sub.get("resolved_span_per_chain"))
    out["sequence_agreement"] = _int_map(sub.get("sequence_agreement"),
                                         ("agreements", "disagreements"))

    # ---- R2: the contact shells
    shells = sub.get("shells")
    out["contact_cutoff_angstrom"] = _float(
        _nested(shells, "contact_cutoff_angstrom"))
    out["contacts_by_copy"] = _label_map(_nested(shells, "contacts_by_copy"))
    out["contact_counts_by_copy"] = _int_map(
        _nested(shells, "contact_counts_by_copy"))
    out["contact_positions_by_copy"] = _position_map(
        _nested(shells, "contact_positions_by_copy"))
    out["nearest_contact_by_copy"] = _str_map(
        _nested(shells, "nearest_contact_by_copy"))

    symmetry = _nested(shells, "chain_symmetry")
    symmetry = symmetry if isinstance(symmetry, dict) else {}
    out["chain_symmetry"] = {}
    for component, block in sorted(symmetry.items()):
        block = block if isinstance(block, dict) else {}
        compared = block.get("copies_compared")
        identical = block.get("position_sets_identical")
        out["chain_symmetry"][str(component)] = {
            "copies_compared": sorted({str(c) for c in compared})
                               if isinstance(compared, list) else [],
            "position_sets_identical": identical if isinstance(identical, bool) else None,
            "symmetric_difference": _positions(block.get("symmetric_difference")),
        }

    # ---- R3: the pocket descriptors
    descriptors = sub.get("descriptors")
    focal = _nested(descriptors, "focal_copy")
    focal = str(focal) if isinstance(focal, str) else None
    out["focal_copy"] = focal
    for axis in ("x", "y", "z"):
        out[f"centroid_{axis}"] = _float(_nested(descriptors, "centroid", axis))
    out["radius_of_gyration"] = _float(_nested(descriptors, "radius_of_gyration"))
    out["grid_census"] = _int_map(_nested(descriptors, "grid_census"), GRID_KEYS)
    out["pocket_grid_points"] = _int(_nested(descriptors, "pocket_grid_points"))
    out["contact_positions_uniprot"] = _positions(
        _nested(descriptors, "contact_positions_uniprot"))
    out["annotated_positions_uniprot"] = _positions(
        _nested(descriptors, "annotated_positions_uniprot"))
    out["confirmation_sets"] = {
        k: _positions(_nested(descriptors, "confirmation_sets", k))
        for k in CONFIRMATION_KEYS}

    # The graded distance is the focal copy's, and the focal copy is identified
    # from the campaign's OWN declared keys -- the focal ligand's component and
    # the focal chain -- not from the submission's R3 block. Reading it from the
    # submission would make an R2 component depend on an R3 answer, so a system
    # that computed the shells correctly and stopped there would fail R2 for
    # reasons that have nothing to do with R2.
    distances = _nested(shells, "nearest_contact_distance_by_copy")
    distances = distances if isinstance(distances, dict) else {}
    component = spec.ligand_components.get(spec.focal_ligand_chebi or "")
    prefix = f"{component}:{spec.focal_chain}:" if component else None
    matching = ([k for k in sorted(distances) if str(k).startswith(prefix)]
                if prefix else [])
    out["nearest_contact_distance_focal"] = (
        _float(distances[matching[0]]) if len(matching) == 1 else None)

    # ---- R4: the cutoff sweep and the truncation sweep
    counterfactuals = sub.get("counterfactuals")
    labels = tuple(f"{c:.1f}" for c in spec.cutoff_variants)
    out["cutoff_sweep_counts"] = _int_map(
        _nested(counterfactuals, "cutoff_sweep_counts"), labels)
    sweep = _nested(counterfactuals, "cutoff_sweep_positions")
    sweep = sweep if isinstance(sweep, dict) else {}
    out["cutoff_sweep_positions"] = {k: _positions(sweep.get(k)) for k in labels}

    truncation = _nested(counterfactuals, "truncation")
    out["truncation_shell_size"] = _int(_nested(truncation, "shell_size"))
    out["truncation_atoms_removed"] = _int_map(
        _nested(truncation, "atoms_removed_by_position"))
    out["truncation_contact_survives"] = _positions(
        _nested(truncation, "contact_survives"))
    out["truncation_contact_lost"] = _positions(_nested(truncation, "contact_lost"))
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
