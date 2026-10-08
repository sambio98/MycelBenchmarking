"""The residues ladder: annotation census, evidence classification, then
confirmation of the annotations against the geometry of the structure that
UniProt itself cites as their evidence.

The last step is the point of the ladder. UniProt does not merely assert that a
position binds a ligand; under ECO:0007744 it names the deposited structure the
assertion rests on. That turns "is this annotation supported?" into arithmetic
over coordinates, with no curator in the loop -- which is the same doctrine the
MIBiG evidence audit applies to a different database.

Runs from a sandbox: inputs/, reference/ and task.yaml, nothing else.
"""

from __future__ import annotations

import pathlib

from npbench_c.templates.structure_features.core import (
    Annotations,
    Rules,
    Structure,
    StructureError,
    uniprot_positions,
)
from npbench_c.templates.structure_features.params import StructureSpec

ACTIVE_SITE = "Active site"
BINDING_SITE = "Binding site"


def _label(record) -> str:
    return f"{record.type}@{record.start}"


def census(annotations: Annotations, rules: Rules) -> dict:
    """R1: what the two sources contain, before any evidence is weighed."""
    counts = {t: 0 for t in rules.functional_types}
    positions: set[int] = set()
    agree = disagree = 0
    for record in annotations.records:
        counts[record.type] += 1
        positions.add(record.start)
        # A mutagenesis feature records the residue it mutated; the sequence
        # records what is there. Comparing them is a cheap internal-consistency
        # check that an agent can only get right by actually indexing the
        # sequence at the annotated position.
        if (record.type == "Mutagenesis" and record.original
                and record.start <= len(annotations.sequence)):
            if record.original == annotations.sequence[record.start - 1]:
                agree += 1
            else:                           # pragma: no cover - defensive
                disagree += 1
    return {
        "accession": annotations.accession,
        "sequence_length": len(annotations.sequence),
        "functional_record_counts": dict(sorted(counts.items())),
        "secondary_structure_counts": dict(annotations.secondary_counts),
        "annotated_positions": sorted(positions),
        "ligand_vocabulary": annotations.ligands(),
        "pdb_cross_references": list(annotations.pdb_entries),
        "mutagenesis_original_residue_agreements": agree,
        "mutagenesis_original_residue_disagreements": disagree,
    }


def evidence_census(annotations: Annotations, rules: Rules) -> dict:
    """R2: classify every record by its ECO codes and index the PDB citations."""
    code_totals: dict[str, int] = {}
    class_counts = {t: {c: 0 for c in rules.class_vocabulary}
                    for t in rules.functional_types}
    binding_positions: dict[str, set[int]] = {}
    pdb_index: dict[str, set[int]] = {}
    unevidenced: list[str] = []
    inferred_only: list[str] = []

    for record in annotations.records:
        for code, source, ident in record.citations:
            code_totals[code] = code_totals.get(code, 0) + 1
            if code == "ECO:0007744" and source == "PDB":
                pdb_index.setdefault(ident.upper(), set()).add(record.start)
        klass = rules.eco_class(record.codes)
        class_counts[record.type][klass] += 1
        if klass == "none":
            unevidenced.append(_label(record))
        elif klass == "inferred":
            inferred_only.append(_label(record))
        if record.type == BINDING_SITE and record.ligand_id:
            binding_positions.setdefault(record.ligand_id, set()).add(record.start)

    return {
        "eco_code_citation_totals": dict(sorted(code_totals.items())),
        "record_class_counts": {t: dict(sorted(c.items()))
                                for t, c in sorted(class_counts.items())},
        "binding_positions_by_ligand": {k: sorted(v)
                                        for k, v in sorted(binding_positions.items())},
        "pdb_evidence_index": {k: sorted(v) for k, v in sorted(pdb_index.items())},
        "unevidenced_records": sorted(set(unevidenced)),
        "inferred_only_records": sorted(set(inferred_only)),
    }


def confirm(annotations: Annotations, structure: Structure, chain: str,
            spec: StructureSpec, rules: Rules) -> dict:
    """Reconcile the annotations with the geometry of one structure.

    An annotation is structure-confirmed when it names this entry as its
    structural evidence AND the ligand it names is present in the coordinates
    AND the annotated position is in that ligand's contact shell at the pinned
    cutoff. All three, in the campaign's declared chain.
    """
    offsets = structure.offsets(annotations.accession)
    auth = offsets["auth_to_uniprot"]
    copies = structure.ligand_copies()
    polymer = structure.polymer()

    per_ligand: dict[str, dict] = {}
    confirmed_positions: dict[str, list[int]] = {}
    contact_cache: dict[str, set[int]] = {}
    for chebi, comp in sorted(spec.ligand_components.items()):
        annotated = set(annotations.positions_for_ligand(chebi, structure.pdb_id))
        present = [atoms for (c, ch, _res), atoms in copies.items()
                   if c == comp and ch == chain]
        contact: set[int] = set()
        for atoms in present:
            hits = structure.contacts(atoms, rules.cutoff, polymer=polymer)
            contact |= set(uniprot_positions(hits, chain, auth))
        contact_cache[chebi] = contact
        confirmed = sorted(annotated & contact)
        confirmed_positions[chebi] = confirmed
        per_ligand[chebi] = {
            "component_id": comp,
            "copies_in_chain": len(present),
            "annotated_citing_this_structure": sorted(annotated),
            "contact_positions": sorted(contact),
            "confirmed": len(confirmed),
            "annotated_not_in_contact": sorted(annotated - contact),
            "contact_not_annotated": sorted(contact - annotated),
        }

    active = [r for r in annotations.records if r.type == ACTIVE_SITE]
    if len(active) != 1:
        raise StructureError(
            f"{annotations.accession}: expected exactly one Active site record, "
            f"found {len(active)}; the campaign's R3 claim is not well posed")
    position = active[0].start
    active_contacts = sorted(chebi for chebi, hit in contact_cache.items()
                             if position in hit)

    annotated_total = sum(len(v["annotated_citing_this_structure"])
                          for v in per_ligand.values())
    confirmed_total = sum(v["confirmed"] for v in per_ligand.values())
    return {
        "structure": structure.pdb_id,
        "chain": chain,
        "numbering": offsets,
        "per_ligand": per_ligand,
        # Flat per-ligand views. per_ligand is the readable form for an auditor;
        # these are the graded form, because a grading spec that walked into
        # per_ligand would need one component per ligand per field and the
        # conjunctive product over 24 terms would make the rung unclearable for
        # reasons that have nothing to do with the task.
        "annotated_by_ligand": {k: v["annotated_citing_this_structure"]
                                for k, v in sorted(per_ligand.items())},
        "contact_by_ligand": {k: v["contact_positions"]
                              for k, v in sorted(per_ligand.items())},
        "confirmed_count_by_ligand": {k: v["confirmed"]
                                      for k, v in sorted(per_ligand.items())},
        "contact_not_annotated_by_ligand": {k: v["contact_not_annotated"]
                                            for k, v in sorted(per_ligand.items())},
        "confirmed_positions": confirmed_positions,
        "totals": {
            "annotated_citing_this_structure": annotated_total,
            "confirmed": confirmed_total,
            "unconfirmed": annotated_total - confirmed_total,
        },
        "active_site_position": position,
        "active_site_residue": annotations.sequence[position - 1],
        "active_site_contacts": active_contacts,
    }


def policy_variants(annotations: Annotations, spec: StructureSpec,
                    rules: Rules) -> dict:
    """R4's second axis: which records survive a different evidence policy.

    Admission is per CITATION class, not per record class: a record citing both
    a published experiment and a structure is admitted by a policy accepting
    either. The census in R2 reports the strongest class per record instead, and
    the two rules are stated separately in eco_policy.json because computing one
    while reporting the other is how two blocks of the same report end up
    disagreeing.
    """
    out: dict[str, dict] = {}
    for name, variant in sorted(spec.policy_variants.items()):
        accepted = frozenset(variant["accepted_classes"])
        counts = {t: 0 for t in rules.functional_types}
        ligands: dict[str, set[int]] = {}
        for record in annotations.records:
            classes = {rules.eco_class([code]) for code in record.codes} or {"none"}
            if not (classes & accepted):
                continue
            counts[record.type] += 1
            if record.type == BINDING_SITE and record.ligand_id:
                ligands.setdefault(record.ligand_id, set()).add(record.start)
        out[name] = {
            "accepted_classes": sorted(accepted),
            "admitted_record_counts": dict(sorted(counts.items())),
            "admitted_binding_positions_by_ligand": {
                k: sorted(v) for k, v in sorted(ligands.items())},
            "admitted_total": sum(counts.values()),
        }
    return out


def analyse(campaign: pathlib.Path) -> dict:
    """Everything the ladder grades, computed from the sandbox alone."""
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    spec = StructureSpec.load(campaign)
    annotations = Annotations.load(campaign / spec.annotation_input, rules)
    if annotations.accession != spec.accession:
        raise StructureError(
            f"annotation input is {annotations.accession}, task.yaml declares "
            f"{spec.accession}")

    focal = Structure.load(campaign / spec.focal_input, spec.focal_pdb_id)
    result = {
        **census(annotations, rules),
        **evidence_census(annotations, rules),
        "focal_structure": {
            "pdb_id": focal.pdb_id,
            "polymer_atom_records": len(focal.polymer()),
            "ligand_copies": sorted(f"{comp}:{ch}:{res}"
                                    for comp, ch, res in focal.ligand_copies()),
        },
        "confirmation": confirm(annotations, focal, spec.focal_chain, spec, rules),
    }

    if not spec.alternate_pdb_id or not spec.alternate_input:
        raise StructureError("the residues ladder needs an alternate structure for R4")
    alternate = Structure.load(campaign / spec.alternate_input, spec.alternate_pdb_id)
    swap = confirm(annotations, alternate, spec.focal_chain, spec, rules)
    swap["confirmed_delta_by_ligand"] = {
        chebi: (swap["per_ligand"][chebi]["confirmed"]
                - result["confirmation"]["per_ligand"][chebi]["confirmed"])
        for chebi in sorted(spec.ligand_components)
    }
    result["structure_swap"] = swap
    result["policy_variants"] = policy_variants(annotations, spec, rules)
    return result
