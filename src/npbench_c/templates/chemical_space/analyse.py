"""Everything the chemical-space ladder grades, from the sandbox alone.

Runs over reference/chemistry_table.json joined to the release archive in
inputs/. Pure stdlib set arithmetic in sorted order; no chemistry toolkit, since
perception already happened once in build_reference.py.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.chemical_space.core import KEYS, ChemicalSpace, SpaceError


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def census(space: ChemicalSpace) -> dict:
    """R1: the corpus, and how many compounds it is under each declared key."""
    records = space.records
    multi_class = sum(1 for classes in space.classes.values() if len(classes) > 1)
    entries_per_class = {
        name: sum(1 for classes in space.classes.values() if name in classes)
        for name in sorted(space.rules.class_vocabulary)}
    return {
        "corpus_release": space.release,
        "corpus_entries": space.corpus_entries,
        "active_entries": space.active_entries,
        "compound_records": len(records),
        "entries_with_compounds": len({r["accession"] for r in records}),
        "distinct_compound_keys": space.distinct(records, "compound"),
        "distinct_connectivity_keys": space.distinct(records, "connectivity"),
        "distinct_scaffold_keys": space.distinct(records, "scaffold"),
        "acyclic_compounds": sum(1 for r in records
                                 if r["scaffold_inchikey"] is None),
        "duplicate_records": len(records) - space.distinct(records, "compound"),
        "class_vocabulary": sorted(space.rules.class_vocabulary),
        "entries_per_class": entries_per_class,
        "multi_class_entries": multi_class,
    }


def granularity(space: ChemicalSpace) -> dict:
    """R2: what each identity key merges, and how the corpus spreads over classes."""
    records = space.records
    return {
        "connectivity_groups_with_multiple_compound_keys":
            space.multi_key_connectivity_groups(records),
        "distinctions_lost_at_connectivity":
            space.distinct(records, "compound")
            - space.distinct(records, "connectivity"),
        "cross_entry_duplicates": {
            key: space.cross_entry_duplicates(records, key) for key in KEYS},
        "records_per_class": space.per_class(records, None),
        "distinct_compounds_per_class": space.per_class(records, "compound"),
        "distinct_scaffolds_per_class": space.per_class(records, "scaffold"),
        "scaffold_ring_count_histogram": space.ring_histogram(records),
    }


def sharing(space: ChemicalSpace, records, min_rings=None, min_coverage=None) -> dict:
    """The filtered-versus-unfiltered sharing statistic over one record set."""
    informative = space.informative(records, min_rings, min_coverage)
    all_cross = space.cross_class(records)
    informative_cross = space.cross_class(informative)
    return {
        "informative_records": len(informative),
        "informative_scaffold_keys": space.distinct(informative, "scaffold"),
        "shared_scaffolds": {
            "all": len(space.shared_scaffolds(records)),
            "informative": len(space.shared_scaffolds(informative)),
        },
        "cross_class_scaffolds": {
            "all": len(all_cross["declared"]),
            "informative": len(informative_cross["declared"]),
            "naive_all": all_cross["naive_count"],
            "naive_informative": informative_cross["naive_count"],
        },
        "_cross_class_informative_entries": informative_cross["declared"],
        "_informative": informative,
    }


def commitment(space: ChemicalSpace) -> dict:
    """R3: the sharing claim that survives the declared informativeness filter."""
    result = sharing(space, space.records)
    entries = result.pop("_cross_class_informative_entries")
    informative = result.pop("_informative")

    shared = space.shared_scaffolds(informative)
    if not shared:
        raise SpaceError("no informative scaffold is shared, so R3 has no claim")
    # Declared tie-break: most accessions, then the smallest key. Without the
    # second clause the answer would depend on iteration order.
    ranked = sorted(((-len(accessions), key) for key, accessions in shared.items()))
    count, key = ranked[0]
    return {
        "informative_rule": {
            "min_ring_count": space.rules.min_ring_count,
            "min_heavy_atom_coverage": space.rules.min_coverage,
        },
        **result,
        "cross_class_informative_entries": entries,
        "most_shared_informative_scaffold": {
            "scaffold_inchikey": key, "entries": -count},
    }


def counterfactuals(space: ChemicalSpace, params: dict) -> dict:
    """R4: vary the declared filter, then gate the corpus on locus evidence."""
    variants = {}
    for name, spec in sorted((params.get("threshold_variants") or {}).items()):
        block = sharing(space, space.records,
                        int(spec["min_ring_count"]),
                        float(spec["min_heavy_atom_coverage"]))
        # Only the numbers the filter moves. The unfiltered counts are the same
        # in every variant by construction, and a constant column in a
        # counterfactual table is a component that cannot discriminate.
        variants[name] = {
            "min_ring_count": int(spec["min_ring_count"]),
            "min_heavy_atom_coverage": float(spec["min_heavy_atom_coverage"]),
            "informative_records": block["informative_records"],
            "informative_scaffold_keys": block["informative_scaffold_keys"],
            "shared_informative": block["shared_scaffolds"]["informative"],
            "cross_class_informative": block["cross_class_scaffolds"]["informative"],
            "naive_cross_class_informative":
                block["cross_class_scaffolds"]["naive_informative"],
        }
    if len({v["informative_records"] for v in variants.values()}) < len(variants):
        raise SpaceError(
            "two threshold variants admit the same number of records, so at "
            f"least one R4 component is not discriminating: "
            f"{ {k: v['informative_records'] for k, v in variants.items()} }")

    gated_records = space.evidence_gated(space.records)
    if not gated_records:
        raise SpaceError("the evidence gate admits no record")
    gate = sharing(space, gated_records)
    gate_entries = gate.pop("_cross_class_informative_entries")
    gate.pop("_informative")
    return {
        "threshold_variants": variants,
        "evidence_gate": {
            "entries": len({r["accession"] for r in gated_records}),
            "compound_records": len(gated_records),
            "distinct_compound_keys": space.distinct(gated_records, "compound"),
            "distinct_scaffold_keys": space.distinct(gated_records, "scaffold"),
            **gate,
            "cross_class_informative_keys": sorted(gate_entries),
        },
    }


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    space = ChemicalSpace.load(campaign, params["corpus_archive"])
    if space.release != params["corpus_release"]:
        raise SpaceError(
            f"the chemistry table was built for release {space.release}, "
            f"task.yaml declares {params['corpus_release']}")
    return {
        **census(space),
        **granularity(space),
        **commitment(space),
        **counterfactuals(space, params),
    }
