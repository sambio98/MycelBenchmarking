"""The pocket ladder: contact shells, pocket descriptors and a geometric
counterfactual, computed from deposited coordinates alone.

Everything here is arithmetic over the numbers in the file. Nothing is
superposed, no symmetry mate is generated, no side chain is rebuilt and no
energy is evaluated, so every reported quantity is a function of the file and
the pinned rules in geometry_rules.json -- which is the only way a geometry task
can be graded without grading the convention instead of the computation.

The ladder's third rung reconciles the computed shell against UniProt's own
binding-site annotations for this entry, so the two sources meet here exactly as
they do in the residues ladder, from the other side.

Runs from a sandbox: inputs/, reference/ and task.yaml, nothing else.
"""

from __future__ import annotations

import pathlib

from npbench_c.templates.structure_features.core import (
    DISTANCE_DECIMALS,
    Annotations,
    Rules,
    Structure,
    StructureError,
    centroid,
    contact_label,
    pocket_volume,
    radius_of_gyration,
    uniprot_positions,
)
from npbench_c.templates.structure_features.params import StructureSpec


def copy_id(key: tuple[str, str, str]) -> str:
    comp, chain, residue = key
    return f"{comp}:{chain}:{residue}"


def inventory(structure: Structure, annotations: Annotations, rules: Rules) -> dict:
    """R1: what the file contains, and how its numbering relates to UniProt."""
    polymer = structure.polymer()
    copies = structure.ligand_copies()
    offsets = structure.offsets(annotations.accession)
    auth = offsets["auth_to_uniprot"]

    per_chain_count: dict[str, int] = {}
    per_chain_span: dict[str, list[int]] = {}
    agree = disagree = 0
    for chain in structure.chains():
        resolved = structure.resolved_residues(chain)
        per_chain_count[chain] = len(resolved)
        per_chain_span[chain] = [min(resolved), max(resolved)] if resolved else []
        for number, comp in resolved.items():
            position = number + auth
            if not 1 <= position <= len(annotations.sequence):
                continue
            if rules.one_letter(comp) == annotations.sequence[position - 1]:
                agree += 1
            else:
                disagree += 1

    return {
        "pdb_id": structure.pdb_id,
        "accession": annotations.accession,
        "atom_inventory": {
            "polymer_atom_records": len(polymer),
            "ligand_atom_records": len(structure.atoms) - len(polymer),
            "water_atom_records": structure.water_atoms,
            "excluded_hydrogen_records": structure.excluded_hydrogens,
            "excluded_altloc_records": structure.excluded_altloc,
        },
        "polymer_chains": list(structure.chains()),
        "ligand_copies": {copy_id(k): len(v) for k, v in copies.items()},
        "numbering": offsets,
        "resolved_residues_per_chain": dict(sorted(per_chain_count.items())),
        "resolved_span_per_chain": dict(sorted(per_chain_span.items())),
        "sequence_agreement": {"agreements": agree, "disagreements": disagree},
    }


def shells(structure: Structure, annotations: Annotations, rules: Rules) -> dict:
    """R2: the contact shell of every ligand copy, and the symmetry check.

    A structure with two copies of the same ligand in two chains is a free
    internal control: the two shells must agree once each is expressed in
    sequence positions, and a system whose numbering is wrong usually fails
    that agreement before it fails anything else.
    """
    polymer = structure.polymer()
    auth = structure.offsets(annotations.accession)["auth_to_uniprot"]
    contacts_by_copy: dict[str, list[str]] = {}
    counts: dict[str, int] = {}
    nearest: dict[str, str] = {}
    nearest_distance: dict[str, float] = {}
    positions_by_copy: dict[str, list[int]] = {}

    for key, atoms in structure.ligand_copies().items():
        name = copy_id(key)
        hits = structure.contacts(atoms, rules.cutoff, polymer=polymer)
        labels = [contact_label(chain, residue, comp, rules)
                  for chain, residue, comp in hits]
        contacts_by_copy[name] = sorted(labels)
        counts[name] = len(hits)
        # Declared tie-break: the smallest distance, and among equal distances
        # the lexicographically smallest label. Without the second clause the
        # answer would depend on dict order.
        ranked = sorted(
            (round(distance, DISTANCE_DECIMALS),
             contact_label(chain, residue, comp, rules))
            for (chain, residue, comp), distance in hits.items())
        if ranked:
            nearest_distance[name], nearest[name] = ranked[0]
        positions_by_copy[name] = list(uniprot_positions(hits, key[1], auth))

    # Copies of the same chemical component, compared on sequence positions.
    groups: dict[str, list[str]] = {}
    for key in structure.ligand_copies():
        groups.setdefault(key[0], []).append(copy_id(key))
    symmetry = {}
    for comp, names in sorted(groups.items()):
        if len(names) < 2:
            continue
        sets = [set(positions_by_copy[n]) for n in sorted(names)]
        union, intersection = set(), set(sets[0])
        for s in sets:
            union |= s
            intersection &= s
        symmetry[comp] = {
            "copies_compared": sorted(names),
            "position_sets_identical": all(s == sets[0] for s in sets),
            "symmetric_difference": sorted(union - intersection),
        }

    return {
        "contact_cutoff_angstrom": rules.cutoff,
        "contacts_by_copy": contacts_by_copy,
        "contact_counts_by_copy": dict(sorted(counts.items())),
        "nearest_contact_by_copy": dict(sorted(nearest.items())),
        "nearest_contact_distance_by_copy": dict(sorted(nearest_distance.items())),
        "contact_positions_by_copy": dict(sorted(positions_by_copy.items())),
        "chain_symmetry": symmetry,
    }


def _focal_copy(structure: Structure, spec: StructureSpec) -> tuple[tuple, tuple]:
    if not spec.focal_ligand_chebi:
        raise StructureError("the pocket ladder needs focal_ligand_chebi")
    component = spec.ligand_components.get(spec.focal_ligand_chebi)
    if not component:
        raise StructureError(
            f"no component declared for focal ligand {spec.focal_ligand_chebi}")
    found = [(k, v) for k, v in structure.ligand_copies().items()
             if k[0] == component and k[1] == spec.focal_chain]
    if len(found) != 1:
        raise StructureError(
            f"{structure.pdb_id}: expected exactly one {component} copy in chain "
            f"{spec.focal_chain}, found {len(found)}")
    return found[0]


def descriptors(structure: Structure, annotations: Annotations, spec: StructureSpec,
                rules: Rules) -> dict:
    """R3: pocket descriptors for the focal copy, reconciled with the annotations."""
    key, atoms = _focal_copy(structure, spec)
    polymer = structure.polymer()
    auth = structure.offsets(annotations.accession)["auth_to_uniprot"]

    hits = structure.contacts(atoms, rules.cutoff, polymer=polymer)
    contact = set(uniprot_positions(hits, spec.focal_chain, auth))
    annotated = set(annotations.positions_for_ligand(
        spec.focal_ligand_chebi, structure.pdb_id))
    if not annotated:
        raise StructureError(
            f"{annotations.accession} has no binding-site annotation for "
            f"{spec.focal_ligand_chebi} citing {structure.pdb_id}, so R3 has "
            "nothing to reconcile")
    confirmed = sorted(annotated & contact)

    cx, cy, cz = centroid(atoms)
    volume = pocket_volume(atoms, polymer, rules)
    true_positive = len(confirmed)
    precision = round(true_positive / len(contact), 6) if contact else 0.0
    recall = round(true_positive / len(annotated), 6) if annotated else 0.0
    return {
        "focal_copy": copy_id(key),
        "focal_ligand": spec.focal_ligand_chebi,
        "heavy_atoms": len(atoms),
        "centroid": {"x": cx, "y": cy, "z": cz},
        "radius_of_gyration": radius_of_gyration(atoms),
        "grid_census": {
            "grid_points_near_ligand": volume["grid_points_near_ligand"],
            "grid_points_rejected_clash": volume["grid_points_rejected_clash"],
            "grid_points_rejected_bulk": volume["grid_points_rejected_bulk"],
        },
        "pocket_grid_points": volume["pocket_grid_points"],
        "pocket_volume_a3": volume["pocket_volume_a3"],
        "contact_positions_uniprot": sorted(contact),
        "annotated_positions_uniprot": sorted(annotated),
        "confirmation_sets": {
            "confirmed": confirmed,
            "annotated_not_in_contact": sorted(annotated - contact),
            "contact_not_annotated": sorted(contact - annotated),
        },
        "annotation_precision": precision,
        "annotation_recall": recall,
    }


def counterfactuals(structure: Structure, annotations: Annotations,
                    spec: StructureSpec, rules: Rules) -> dict:
    """R4: the two ways the shell is perturbed -- the cutoff, and a side chain."""
    key, atoms = _focal_copy(structure, spec)
    polymer = structure.polymer()
    auth = structure.offsets(annotations.accession)["auth_to_uniprot"]

    sweep_counts: dict[str, int] = {}
    sweep_positions: dict[str, list[int]] = {}
    for cutoff in spec.cutoff_variants:
        hits = structure.contacts(atoms, cutoff, polymer=polymer)
        label = f"{cutoff:.1f}"
        sweep_counts[label] = len(hits)
        sweep_positions[label] = list(uniprot_positions(hits, spec.focal_chain, auth))
    if len(set(map(tuple, sweep_positions.values()))) < 2:
        raise StructureError(
            "every declared cutoff gives the same contact set, so the sweep has "
            "no counterfactual content")

    if spec.truncation_scope != "focal_shell":
        raise StructureError(
            "the pocket ladder needs truncation_scope: focal_shell; the positions "
            "are the shell itself and must not be listed in the sandbox, because "
            "listing them would hand over the previous rung's answer")

    before = structure.contacts(atoms, rules.cutoff, polymer=polymer)
    shell = sorted(before)
    removed_by_position: dict[str, int] = {}
    survivors: list[int] = []
    lost: list[int] = []
    for chain, residue, comp in shell:
        if chain != spec.focal_chain or not residue.lstrip("-").isdigit():
            continue
        position = int(residue) + auth
        if comp in ("GLY", "ALA"):
            # Nothing to truncate, so the contact survives by definition. The
            # row is kept rather than skipped so the table covers the shell.
            removed_by_position[str(position)] = 0
            survivors.append(position)
            continue
        kept, removed = structure.truncated_polymer(chain, residue, rules)
        after = structure.contacts(atoms, rules.cutoff, polymer=kept)
        removed_by_position[str(position)] = removed
        if (chain, residue, comp) in after:
            survivors.append(position)
        else:
            lost.append(position)

    if not survivors or not lost:
        raise StructureError(
            "truncating the shell gives the same verdict everywhere "
            f"({len(survivors)} survive, {len(lost)} lost), so R4's in-silico "
            "substitution has a safe guess and no discriminating content")

    return {
        "cutoff_sweep_counts": sweep_counts,
        "cutoff_sweep_positions": sweep_positions,
        "truncation": {
            "scope": spec.truncation_scope,
            "shell_size": len(shell),
            "atoms_removed_by_position": dict(sorted(
                removed_by_position.items(), key=lambda kv: int(kv[0]))),
            "contact_survives": sorted(survivors),
            "contact_lost": sorted(lost),
        },
    }


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    spec = StructureSpec.load(campaign)
    annotations = Annotations.load(campaign / spec.annotation_input, rules)
    if annotations.accession != spec.accession:
        raise StructureError(
            f"annotation input is {annotations.accession}, task.yaml declares "
            f"{spec.accession}")
    structure = Structure.load(campaign / spec.focal_input, spec.focal_pdb_id)
    return {
        **inventory(structure, annotations, rules),
        "shells": shells(structure, annotations, rules),
        "descriptors": descriptors(structure, annotations, spec, rules),
        "counterfactuals": counterfactuals(structure, annotations, spec, rules),
    }
