"""Does each campaign's answer cohere with itself, and with the chemistry?

Check 27 and check 28 are about whether an answer can be WRITTEN. This file is
about whether the answer is RIGHT -- as far as arithmetic can say so without a
second implementation of each campaign's science.

Every assertion here is an identity the answer asserts about itself and that a
wrong answer would usually break: a census that partitions a set must sum to it,
a rate must equal its own numerator over its own denominator, a committed set
must be the size its census claims, a signed difference must be the difference
of the two numbers it sits beside, and a peptide's formula must be the sum of
its monomers less one water per bond. The constants used here (the monoisotopic
masses) are written out rather than imported, so the chemistry is checked
against an independent source rather than against the code that produced it.

A generic version of this file -- scan for negative counts, rates outside
[0, 1], partitions that do not sum -- was written first and reported zero
findings over all 25 golds. It was measured and found vacuous: four campaigns
evaluated no assertions at all, most evaluated fewer than five, and the
partition rule never matched a single key pair. Hence identities named one at a
time, per campaign.
"""

from __future__ import annotations

import json
import math
import pathlib
from collections import Counter

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2] / "campaigns"

#: Monoisotopic masses, from the standard atomic-mass tables rather than from
#: the template that computed the formulas under test.
MASS = {"C": 12.0, "H": 1.0078250319, "N": 14.0030740052, "O": 15.9949146221,
        "P": 30.97376151, "S": 31.97207069, "Cl": 34.96885271, "Br": 78.9183376,
        "F": 18.9984032, "Fe": 55.9349421, "Zn": 63.9291466}


def gold(name: str) -> dict:
    path = ROOT / name / "gold" / "gold.json"
    if not path.is_file():
        pytest.skip(f"{name} has no gold")
    return json.loads(path.read_text())


def formula(text: str) -> Counter:
    """Parse a Hill-order formula into element counts."""
    import re

    out: Counter = Counter()
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", text):
        if element:
            out[element] += int(count or 1)
    return out


def mass_of(counts: Counter) -> float:
    return sum(MASS[element] * n for element, n in counts.items())


# ---------------------------------------------------------------- arch (x2)


@pytest.mark.parametrize("name", ["arch-nrps-mibig-4_0-01", "arch-pks-mibig-4_0-01"])
def test_arch_censuses_partition_the_corpus(name):
    g = gold(name)
    genes = g["selection_census"]["corpus_genes"]
    assert genes > 0
    for table in ("architecture_table", "module_table", "reconciliation_table"):
        assert len(g[table]) == genes, table
    assert g["architecture_census"]["genes"] == genes
    assert g["protein_length_summary"]["count"] == genes
    for census in ("primary_verdicts", "secondary_verdicts", "module_verdicts",
                   "primary_delta_census", "curated_module_census",
                   "strata_census"):
        assert sum(g[census].values()) == genes, census
    # the per-stratum split of one census must rebuild that census
    pooled: Counter = Counter()
    for block in g["secondary_verdicts_by_stratum"].values():
        pooled.update(block)
    assert dict(pooled) == {k: v for k, v in g["secondary_verdicts"].items()}


@pytest.mark.parametrize("name", ["arch-nrps-mibig-4_0-01", "arch-pks-mibig-4_0-01"])
def test_arch_totals_are_the_tables_they_summarise(name):
    g = gold(name)
    rows = g["architecture_table"]
    assert g["architecture_census"]["domains_total"] == sum(r["domains"] for r in rows)
    assert g["architecture_census"]["genes_with_no_domain"] == sum(
        1 for r in rows if r["domains"] == 0)
    assert g["architecture_census"]["distinct_architectures"] == len(
        {r["architecture"] for r in rows})
    assert g["overlap_census"]["hits_total"] == sum(g["panel_census"].values())
    assert g["complete_modules_total"] == sum(
        r["complete_modules"] for r in g["module_table"])


@pytest.mark.parametrize("name", ["arch-nrps-mibig-4_0-01", "arch-pks-mibig-4_0-01"])
def test_arch_the_baseline_cell_is_the_baseline(name):
    """The sweep must contain the headline answer, not merely surround it."""
    g = gold(name)
    cell = g["policy_sweep"][
        f"{g['baseline_cutoff']}|{g['baseline_resolution']}|{g['baseline_span']}"]
    assert cell["domains_total"] == g["architecture_census"]["domains_total"]
    assert cell["distinct_architectures"] == \
        g["architecture_census"]["distinct_architectures"]
    assert cell["genes_differing_from_baseline"] == 0
    assert cell["primary_domain_agreement"] == g["primary_verdicts"]["agrees"]
    assert g["mapping_sweep"][g["baseline_mapping"]]["verdicts"] == \
        g["primary_verdicts"]


@pytest.mark.parametrize("name", ["arch-nrps-mibig-4_0-01", "arch-pks-mibig-4_0-01"])
def test_arch_each_row_agrees_with_its_own_verdict(name):
    g = gold(name)
    delta = Counter()
    for row in g["reconciliation_table"]:
        assert row["primary_delta"] == row["computed_primary"] - row["curated_primary"]
        expected = ("agrees" if row["primary_delta"] == 0
                    else "computed_more" if row["primary_delta"] > 0
                    else "computed_fewer")
        assert row["primary_verdict"] == expected, row
        delta[str(row["primary_delta"])] += 1
    assert dict(delta) == {k: v for k, v in g["primary_delta_census"].items()}


# -------------------------------------------------------- massbalance (x2)


@pytest.mark.parametrize("name", ["massbalance-nrps-malleobactin-01",
                                  "massbalance-nrps-sevadicin-01"])
def test_massbalance_assembly_is_the_monomers_less_one_water_per_bond(name):
    g = gold(name)
    assert g["peptide_bonds"] == g["n_modules"] - 1
    assert len(g["monomers"]) == g["n_modules"]
    total: Counter = Counter()
    for monomer in g["monomers"]:
        total.update(formula(g["monomer_formulas"][monomer]))
    total["H"] -= 2 * g["peptide_bonds"]
    total["O"] -= 1 * g["peptide_bonds"]
    assert +total == +formula(g["naive_assembly_formula"]), name
    assert math.isclose(g["naive_assembly_monoisotopic_mass"],
                        mass_of(total), abs_tol=5e-4)


@pytest.mark.parametrize("name", ["massbalance-nrps-malleobactin-01",
                                  "massbalance-nrps-sevadicin-01"])
def test_massbalance_residuals_are_product_minus_assembly(name):
    g = gold(name)
    naive = formula(g["naive_assembly_formula"])
    for compound, block in g["per_compound"].items():
        product = formula(block["product_formula"])
        residual = Counter(product)
        residual.subtract(naive)
        stated = Counter(block["residual"])
        assert {k: v for k, v in residual.items() if v} == dict(stated), compound
        assert math.isclose(block["product_monoisotopic_mass"],
                            mass_of(product), abs_tol=5e-4), compound


@pytest.mark.parametrize("name", ["massbalance-nrps-malleobactin-01",
                                  "massbalance-nrps-sevadicin-01"])
def test_massbalance_a_deletion_removes_one_monomer_and_one_bond(name):
    g = gold(name)
    naive = formula(g["naive_assembly_formula"])
    for index, block in g["single_module_deletions"].items():
        gone = formula(g["monomer_formulas"][block["deleted_monomer"]])
        expected = Counter(naive)
        expected.subtract(gone)
        expected["H"] += 2          # the bond that is no longer formed
        expected["O"] += 1
        assert +expected == +formula(block["formula"]), index
        assert math.isclose(block["mass_shift"],
                            mass_of(expected) - mass_of(naive), abs_tol=5e-4), index
        assert block["peptide_bonds"] == g["peptide_bonds"] - 1


# ---------------------------------------------------------------- transfer


def test_transfer_every_query_is_accounted_for():
    g = gold("transfer-ec1_14-uniprot-01")
    assert g["entries_total"] == g["queries_with_a_hit"] + g["queries_without_a_hit"]
    assert sum(g["ec_codes_per_entry"].values()) == g["entries_total"]
    assert g["entries_with_multiple_ec"] == \
        g["entries_total"] - g["ec_codes_per_entry"]["1"]
    assert sum(g["top_hit_identity_census"].values()) == g["queries_with_a_hit"]
    assert sum(g["verdict_census"].values()) == g["queries_with_a_hit"]
    assert len(g["transfer_errors_at_subsubclass"]) == \
        g["verdict_census"]["wrong_at_subsubclass"]


def test_transfer_rates_are_their_own_numerators_over_denominators():
    g = gold("transfer-ec1_14-uniprot-01")
    for level, block in g["agreement_by_depth"].items():
        assert block["agree"] <= block["decidable"], level
        assert math.isclose(block["rate"], block["agree"] / block["decidable"],
                            abs_tol=1e-6), level
    deepest = max(g["agreement_by_depth"], key=int)
    banded = g["agreement_by_identity"]
    assert sum(b["agree"] for b in banded.values()) == \
        g["agreement_by_depth"][deepest]["agree"]
    assert sum(b["decidable"] for b in banded.values()) == \
        g["agreement_by_depth"][deepest]["decidable"]
    for variant, census in g["removal_census"].items():
        assert sum(census.values()) == g["queries_with_a_hit"], variant


# ------------------------------------------------------------- selectivity


def test_selectivity_admission_partitions_the_activities():
    g = gold("selectivity-saureus-topoisomerase-01")
    assert sum(g["admitted_by_target"].values()) == g["variant_admitted"]["v0_baseline"]
    assert sum(g["aggregated_by_target"].values()) == g["aggregated_pairs"]
    for target, admitted in g["admitted_by_target"].items():
        rejected = sum(g["rejected_by_target"][target].values())
        assert admitted + rejected <= g["activities_of_declared_type"], target
    assert g["eligible_molecules"] <= min(g["aggregated_by_target"].values())
    assert g["distinct_assays"] == sum(g["assay_confidence_scores"].values())
    assert len(g["focal_selective_at_highest_threshold"]) == \
        g["focal_selective_counts"][max(g["focal_selective_counts"], key=float)]


def test_selectivity_counts_fall_as_the_threshold_rises():
    g = gold("selectivity-saureus-topoisomerase-01")
    for census in (g["focal_selective_counts"], g["reference_selective_counts"]):
        ordered = [census[k] for k in sorted(census, key=float)]
        assert ordered == sorted(ordered, reverse=True), census
    assert g["most_selective_ratio"] == max(
        -ratio for ratio, _ in g["ratio_extremes"]["most_focal_selective"])
    assert g["most_selective_molecule"] == \
        g["ratio_extremes"]["most_focal_selective"][0][1]


# ------------------------------------------------------------------ derep


def test_derep_every_query_gets_exactly_one_verdict():
    g = gold("derep-massbank-2026_03-01")
    queries = g["corpus_census"]["graded_queries"]
    assert sum(g["confidence_census"].values()) == queries
    assert sum(g["retrieval_census"].values()) == queries
    assert sum(g["adduct_misassignment_census"].values()) == queries
    assert len(g["top_matches"]) == queries
    assert len(g["confidence_by_query"]) == queries
    assert len(g["correct_queries"]) + len(g["misannotated_queries"]) == \
        g["retrieval_census"]["annotated"]
    assert g["candidate_census"]["median"] <= g["candidate_census"]["max"]
    assert g["confidence_census"]["no_candidate"] == \
        g["retrieval_census"]["no_candidate"]
    assert g["confidence_census"]["no_match"] == g["retrieval_census"]["no_match"]


# ---------------------------------------------------------------- network


def test_network_every_node_gets_one_annotation_and_one_claim():
    g = gold("network-massbank-2026_03-01")
    nodes = g["network_census"]["nodes"]
    assert sum(g["annotation_census"].values()) == nodes
    assert sum(g["propagation_census"].values()) == nodes
    assert len(g["propagation_claims"]) == nodes
    assert len(g["top_matches"]) == nodes
    assert g["network_census"]["edges"] == len(g["edges"])
    assert g["network_census"]["nodes_with_an_edge"] + \
        g["network_census"]["singletons"] == nodes
    assert g["network_census"]["largest_component"] <= nodes
    assert g["network_census"]["components"] <= nodes
    assert len(g["mass_shifts"]) == \
        g["propagation_census"]["analogue_with_mass_shift"]


def test_network_the_truth_census_covers_what_was_claimed():
    g = gold("network-massbank-2026_03-01")
    assert g["annotation_truth"]["accepted"] == \
        g["annotation_truth"]["correct"] + g["annotation_truth"]["misannotated"]
    assert g["annotation_truth"]["accepted"] == \
        g["propagation_census"]["annotated_directly"]
    truth = g["propagation_truth"]
    assert truth["analogue_a_different_compound"] == \
        g["propagation_census"]["analogue_with_mass_shift"]
    assert truth["same_compound_right"] + truth["same_compound_wrong"] == \
        g["propagation_census"]["same_compound"]


# ----------------------------------------------------------------- pocket


def test_pocket_the_confirmation_sets_partition_both_position_sets():
    g = gold("pocket-rebh-01")
    d = g["descriptors"]
    annotated = set(d["annotated_positions_uniprot"])
    contact = set(d["contact_positions_uniprot"])
    sets = d["confirmation_sets"]
    assert set(sets["confirmed"]) == annotated & contact
    assert set(sets["annotated_not_in_contact"]) == annotated - contact
    assert set(sets["contact_not_annotated"]) == contact - annotated
    assert math.isclose(d["annotation_precision"], len(sets["confirmed"]) / len(contact),
                        abs_tol=1e-6)
    assert math.isclose(d["annotation_recall"], len(sets["confirmed"]) / len(annotated),
                        abs_tol=1e-6)


def test_pocket_the_shell_is_the_union_of_its_copies():
    g = gold("pocket-rebh-01")
    shells = g["shells"]
    by_copy = shells["contact_positions_by_copy"]
    for copy, positions in by_copy.items():
        assert shells["contact_counts_by_copy"][copy] == len(positions), copy
    assert set(g["descriptors"]["contact_positions_uniprot"]) == {
        p for positions in by_copy.values() for p in positions}
    for component, block in shells["chain_symmetry"].items():
        first, second = (set(by_copy[c]) for c in block["copies_compared"])
        assert block["position_sets_identical"] == (first == second), component
        assert set(block["symmetric_difference"]) == first ^ second, component
    assert g["sequence_agreement"]["agreements"] == sum(
        g["resolved_residues_per_chain"].values())


def test_pocket_a_wider_cutoff_never_loses_a_contact():
    g = gold("pocket-rebh-01")
    sweep = g["counterfactuals"]["cutoff_sweep_positions"]
    counts = g["counterfactuals"]["cutoff_sweep_counts"]
    previous: set[int] = set()
    for cutoff in sorted(sweep, key=float):
        positions = set(sweep[cutoff])
        assert counts[cutoff] == len(positions), cutoff
        assert previous <= positions, cutoff
        previous = positions


# --------------------------------------------------------------- residues


def test_residues_record_classes_partition_the_records():
    g = gold("residues-prna-01")
    for feature, counts in g["record_class_counts"].items():
        assert sum(counts.values()) == g["functional_record_counts"][feature], feature
    assert g["mutagenesis_check"]["agreements"] == \
        g["mutagenesis_original_residue_agreements"]
    assert g["mutagenesis_check"]["disagreements"] == \
        g["mutagenesis_original_residue_disagreements"]
    assert g["mutagenesis_check"]["agreements"] + \
        g["mutagenesis_check"]["disagreements"] == \
        g["functional_record_counts"]["Mutagenesis"]
    for variant, total in g["policy_admitted_totals"].items():
        assert sum(g["policy_admitted_counts"][variant].values()) == total, variant


def test_residues_confirmation_is_annotated_intersected_with_contact():
    g = gold("residues-prna-01")
    for block in (g["confirmation"], g["structure_swap"]):
        annotated = block["annotated_by_ligand"]
        contact = block["contact_by_ligand"]
        for ligand, positions in annotated.items():
            confirmed = set(block["confirmed_positions"][ligand])
            assert confirmed == set(positions) & set(contact[ligand]), ligand
            assert block["confirmed_count_by_ligand"][ligand] == len(confirmed)
            assert set(block["contact_not_annotated_by_ligand"][ligand]) == \
                set(contact[ligand]) - set(positions), ligand
        assert block["totals"]["confirmed"] == sum(
            block["confirmed_count_by_ligand"].values())
        assert block["totals"]["annotated_citing_this_structure"] == \
            block["totals"]["confirmed"] + block["totals"]["unconfirmed"]


# ---------------------------------------------------------------- homology


def test_homology_recovery_splits_what_the_labels_assert():
    g = gold("homology-twilight-ipr036188-01")
    for query, block in g["recovery"].items():
        assert block["same_family_detected"] + block["same_family_missed"] == \
            block["same_family_asserted"], query
        assert math.isclose(block["recovery_rate"],
                            block["same_family_detected"] /
                            block["same_family_asserted"], abs_tol=1e-6), query
        assert block["hits"] == g["baseline_hits"][query], query
    for query, bands in g["baseline_band_sizes"].items():
        assert sum(bands.values()) == g["baseline_hits"][query], query
        for band, size in bands.items():
            assert sum(g["band_census"][query][band].values()) == size, (query, band)
    assert sum(g["verdict_census"].values()) == sum(
        len(bands) for bands in g["band_verdicts"].values())


# --------------------------------------------------------------- gcfcutoff


def test_gcfcutoff_pair_arithmetic_is_the_designed_families():
    g = gold("gcfcutoff-mibig-32-01")
    designed = sum(n * (n - 1) // 2 for n in g["designed_families"].values())
    assert g["designed_pairs"] == designed
    assert g["recovered_pairs"] + g["missed_pairs"] == g["designed_pairs"]
    assert g["recovered_pairs"] == g["co_clustered_pairs"] - g["false_join_pairs"]
    union = g["designed_pairs"] + g["co_clustered_pairs"] - g["recovered_pairs"]
    assert math.isclose(g["pair_jaccard"], g["recovered_pairs"] / union, abs_tol=1e-6)
    assert g["pair_statistics"] == {
        "designed_pairs": g["designed_pairs"], "co_clustered_pairs":
        g["co_clustered_pairs"], "recovered_pairs": g["recovered_pairs"],
        "missed_pairs": g["missed_pairs"], "false_join_pairs": g["false_join_pairs"],
        "pair_jaccard": g["pair_jaccard"]}


def test_gcfcutoff_groups_partition_the_clustered_records():
    g = gold("gcfcutoff-mibig-32-01")
    assert sum(g["size_census"].values()) == g["groups"]
    # the size census covers the whole corpus; `clustered_records` is the corpus
    # less the singletons, which is what "clustered" means here.
    assert sum(int(size) * n for size, n in g["size_census"].items()) == \
        g["corpus_size"]
    assert g["singletons"] == g["size_census"].get("1", 0)
    assert g["clustered_records"] == g["corpus_size"] - g["singletons"]
    assert g["largest_group"] == max(int(s) for s in g["size_census"])
    assert sum(g["verdict_census"].values()) == len(g["designed_families"])
    for family, block in g["per_family"].items():
        assert block["groups_spanned"] == len(block["group_sizes"]), family
        assert sum(block["group_sizes"]) == g["designed_families"][family], family


# ----------------------------------------------------------------- genecall


def test_genecall_the_fragmentation_summary_balances():
    g = gold("genecall-scoelicolor-01")
    s = g["fragmentation_summary"]
    assert s["retained"] + s["lost"] == s["baseline_calls"]
    assert s["retained"] + s["gained"] == s["fragmented_calls"]
    assert s["baseline_calls"] == g["call_census"]["calls"]
    assert s["fragmented_calls"] == g["fragmented_call_census"]["calls"]
    assert len(g["fragmented_lost_table"]) == s["lost"]
    assert sum(g["fragmented_loss_verdicts"].values()) == s["lost"]
    assert sum(g["fragmented_gained_partial_census"].values()) == s["gained"]


def test_genecall_the_call_census_is_the_call_table():
    g = gold("genecall-scoelicolor-01")
    rows = g["call_table"]
    c = g["call_census"]
    assert c["calls"] == len(rows)
    assert c["on_the_forward_strand"] + c["on_the_reverse_strand"] == c["calls"]
    assert sum(g["partial_census"].values()) == c["calls"]
    assert sum(g["start_type_census"].values()) == c["calls"]
    assert math.isclose(c["coding_fraction"], c["coding_bases"] / g["contig_length"],
                        abs_tol=1e-6)
    assert c["coding_bases"] <= g["contig_length"]
    modes = g["mode_difference"]["single_vs_meta"]
    assert modes["shared"] + modes["only_in_single"] == c["calls"]


# --------------------------------------------------------------- bgcdetect


def test_bgcdetect_regions_and_loci_are_each_accounted_for():
    g = gold("bgcdetect-scoelicolor-01")
    assert g["regions"] == len(g["region_table"])
    assert g["regions_matching_a_locus"] + g["regions_without_a_locus"] == g["regions"]
    assert g["loci_total"] == len(g["per_locus"])
    assert sum(g["verdict_census"].values()) == g["loci_total"]
    assert g["verdict_census"]["not_detected"] + g["regions_matching_a_locus"] >= \
        g["loci_total"] - len(g["regions_holding_several_loci"])
    assert g["distinct_products"] == len(g["product_census"])
    assert math.isclose(g["fraction_of_record"],
                        g["region_bases"] / g["record_length"], abs_tol=1e-6)
    assert g["region_bases"] <= g["record_length"]
    assert g["jaccard_summary"] == {"min": g["jaccard_min"],
                                    "median": g["jaccard_median"],
                                    "max": g["jaccard_max"]}
    assert g["jaccard_min"] <= g["jaccard_median"] <= g["jaccard_max"]


# ---------------------------------------------------------------- chemspace


def test_chemspace_keys_collapse_monotonically():
    """Compound keys refine connectivity keys, which refine scaffolds."""
    g = gold("chemspace-mibig-4_0-01")
    assert g["distinct_scaffold_keys"] <= g["distinct_connectivity_keys"] <= \
        g["distinct_compound_keys"] <= g["compound_records"]
    assert g["distinctions_lost_at_connectivity"] == \
        g["distinct_compound_keys"] - g["distinct_connectivity_keys"]
    assert g["duplicate_records"] == g["compound_records"] - g["distinct_compound_keys"]
    assert g["distinct_keys"] == {"compound": g["distinct_compound_keys"],
                                  "connectivity": g["distinct_connectivity_keys"],
                                  "scaffold": g["distinct_scaffold_keys"]}
    assert g["corpus_counts"] == {
        "corpus_entries": g["corpus_entries"], "active_entries": g["active_entries"],
        "entries_with_compounds": g["entries_with_compounds"],
        "compound_records": g["compound_records"]}
    assert g["active_entries"] <= g["corpus_entries"]
    assert g["entries_with_compounds"] <= g["corpus_entries"]
    # the histogram counts RECORDS by the ring count of their scaffold
    assert sum(g["scaffold_ring_count_histogram"].values()) == g["compound_records"]
    assert g["acyclic_compounds"] == g["scaffold_ring_count_histogram"]["0"]
    assert g["informative_records"] <= g["compound_records"]
    assert g["informative_scaffold_keys"] <= g["distinct_scaffold_keys"]


def test_chemspace_a_class_split_covers_at_least_the_corpus():
    g = gold("chemspace-mibig-4_0-01")
    # an entry may carry several classes, so a per-class split over-counts it
    # once per extra class -- at least one extra for each multi-class entry
    assigned = sum(g["entries_per_class"].values())
    assert assigned >= g["entries_with_compounds"]
    assert assigned - g["entries_with_compounds"] >= g["multi_class_entries"]
    assert sum(g["records_per_class"].values()) >= g["compound_records"]
    assert set(g["entries_per_class"]) == set(g["class_vocabulary"])


# ------------------------------------------------------- constraint (x2)


@pytest.mark.parametrize("name,feasible", [("constraint-gc55-ecoli-rebh-01", True),
                                           ("constraint-gc62-ecoli-rebh-01", False)])
def test_constraint_the_verdict_is_the_arithmetic(name, feasible):
    """The whole claim: is the required GC floor inside the reachable range?"""
    g = gold(name)
    reachable = g["gc_global_min"] <= g["restricted_gc_max"]
    assert reachable is feasible, (g["gc_global_min"], g["restricted_gc_max"])
    assert g["verdict"].startswith("feasible" if feasible else "infeasible")
    assert g["expect"] == ("feasible" if feasible else "infeasible")
    assert (g["constructive_witness"] is not None) is feasible
    assert bool(g["conflicting_constraints"]) is not feasible
    assert g["restricted_gc_max"] <= g["unrestricted_gc_max"]
    assert g["restricted_gc_min"] >= g["unrestricted_gc_min"]
    assert g["orf_length_nt"] == 3 * (g["protein_length_aa"] + 1)
    w = g["witness"]
    assert w["required_gc_min"] == g["gc_global_min"]
    assert w["max_achievable_gc_under_forbidden_codons"] == g["restricted_gc_max"]
    # the margin is signed the way the question is asked: how much GC the floor
    # demands beyond what the codon set can reach, so feasible means <= 0
    assert math.isclose(w["margin"], w["required_gc_min"]
                        - w["max_achievable_gc_under_forbidden_codons"],
                        abs_tol=1e-6)
    assert (w["margin"] <= 0) is feasible


def test_constraint_a_feasible_witness_satisfies_what_it_witnesses():
    g = gold("constraint-gc55-ecoli-rebh-01")
    witness = g["constructive_witness"]
    assert witness["encodes_declared_protein"] is True
    assert witness["gc_fraction"] >= g["gc_global_min"]
    assert sum(witness["violations"].values()) == 0


# --------------------------------------------------------- construct (x2)


@pytest.mark.parametrize("name", ["construct-ecoli-prna-01", "construct-ecoli-rebh-01"])
def test_construct_the_flip_is_the_symmetric_difference(name):
    g = gold(name)
    first, second = (set(block["non_binding_constraints"])
                     for block in g["hosts"].values())
    assert set(g["flipped_constraints"]) == first ^ second
    binding = [block["binding_constraint"] for block in g["hosts"].values()]
    assert g["binding_constraint_flips"] == int(binding[0] != binding[1])
    for host, block in g["hosts"].items():
        # a constraint is non-binding when relaxing it buys no CAI
        costs = block["cost_spread_across_policies"]
        assert block["binding_constraint"] not in block["non_binding_constraints"], host
        assert set(block["non_binding_constraints"]) <= set(costs), host
        assert block["cai_floor"] <= block["cai_min_across_policies"], host
        assert math.isclose(block["cai_spread_across_policies"],
                            block["cai_max_across_policies"]
                            - block["cai_min_across_policies"], abs_tol=1e-6), host


# ----------------------------------------------------------------- ecaudit


def test_ecaudit_present_and_absent_partition_every_section():
    g = gold("ecaudit-fmn-dh-1_1_3_15-01")
    assert sum(g["sequences_by_section"].values()) == g["sequences_total"]
    for section, total in g["sequences_by_section"].items():
        assert g["canonical_present_by_section"][section] + \
            g["canonical_absent_by_section"][section] == total, section
    assert sum(g["canonical_absent_by_section"].values()) == g["canonical_absent_total"]
    assert g["family_census"][g["canonical_model"]] == sum(
        g["canonical_present_by_section"].values())
    assert sum(g["architecture_census"].values()) == g["sequences_total"]
    assert sum(g["absent_architecture_census"].values()) == g["canonical_absent_total"]
    assert g["threshold_absent"][g["threshold"]] == g["canonical_absent_total"]
    assert g["misannotated_by_policy"]["p0_canonical_only"] == \
        g["canonical_absent_total"]
    assert g["length_profile"]["canonical_absent"]["count"] == \
        g["canonical_absent_total"]
    assert g["model_panel_size"] == len(g["family_census"])
    # a looser policy can only admit more, so it can only misannotate fewer
    policies = [g["misannotated_by_policy"][k]
                for k in sorted(g["misannotated_by_policy"])]
    assert policies == sorted(policies, reverse=True)


# --------------------------------------------------------------- frameshift


def test_frameshift_every_target_has_one_verdict_and_one_response():
    g = gold("frameshift-scoelicolor-01")
    targets = g["target_length_summary"]["targets"]
    assert len(g["target_table"]) == targets
    assert len(g["response_table"]) == targets
    assert len(g["verdict_by_length"]) == targets
    assert sum(g["verdict_census"].values()) == targets
    assert dict(Counter(g["verdict_by_length"].values())) == {
        k: v for k, v in g["verdict_census"].items() if v}
    for length, row in g["response_table"].items():
        target = g["target_table"][length]
        assert row["calls_before"] == g["baseline_calls"], length
        assert row["target_start"] == target["start"], length
        assert row["target_end"] == target["end"], length
        assert row["target_length"] == target["length"], length
        assert target["end"] - target["start"] + 1 == target["length"], length
        assert target["start"] <= target["deleted_at"] <= target["end"], length
        assert row["verdict"] == g["verdict_by_length"][length], length
    assert g["target_length_summary"]["shortest"] == min(
        r["length"] for r in g["target_table"].values())
    assert g["target_length_summary"]["longest"] == max(
        r["length"] for r in g["target_table"].values())


# ----------------------------------------------------------------- kinetics


def test_kinetics_bands_partition_the_checked_triples():
    g = gold("kinetics-brenda-ec1_1-01")
    assert g["checked"] == g["complete_triples"]
    assert g["completeness_counts"]["3"] == g["complete_triples"]
    assert sum(g["band_counts"].values()) == g["checked"]
    assert math.isclose(g["agreement_rate"], g["band_counts"]["agree"] / g["checked"],
                        abs_tol=1e-6)
    for field, records in g["records_by_field"].items():
        assert g["parsed_by_field"][field] + g["sentinel_by_field"][field] + \
            g["unparsed_by_field"][field] == records, field
    assert g["sentinel_and_unparsed"] == {"sentinel": g["sentinel_by_field"],
                                          "unparsed": g["unparsed_by_field"]}
    for variant, agree in g["unit_error_agree"].items():
        assert math.isclose(g["unit_error_rates"][variant], agree / g["checked"],
                            abs_tol=1e-6), variant
        assert g["unit_error_variants"][variant]["agree"] == agree, variant
        assert g["unit_error_variants"][variant]["checked"] == g["checked"], variant
    assert g["unit_error_agree"]["v0_identity"] == g["band_counts"]["agree"]


# --------------------------------------------------------------- mibig-diff


def test_mibig_diff_the_two_releases_partition_into_shared_and_only():
    g = gold("mibig-diff-3_1-to-4_0-01")
    assert g["shared_entries"] + g["only_in_source"] == g["source_entries"]
    assert g["shared_entries"] + g["only_in_target"] == g["target_entries"]
    assert sum(g["shared_status_in_target"].values()) == g["shared_entries"]
    for field, raw in g["raw_difference_counts"].items():
        classified = g["classified_verdict_counts"][field]
        assert sum(classified.values()) >= raw, field
        assert g["raw_versus_classified_gap"][field] == \
            raw - classified.get("substantive_change", 0), field
        assert raw <= g["shared_entries"], field
    assert set(g["raw_difference_counts"]) == set(g["fields"])


# ------------------------------------------------- mibig-gene-function


def test_gene_function_evidence_splits_every_annotation():
    g = gold("mibig-gene-function-evidence-01")
    assert sum(g["admissible_by_function"].values()) == g["admissible_total"]
    assert sum(g["backed_by_function"].values()) + \
        sum(g["unbacked_by_function"].values()) == g["total_function_annotations"]
    for function, backed in g["backed_by_function"].items():
        assert g["backed_minus_admissible"][function] == \
            backed - g["admissible_by_function"][function], function
        assert g["admissible_by_function"][function] <= backed, function
    assert set(g["backed_by_function"]) == set(g["function_vocabulary"])
    assert g["admissible_total"] <= sum(g["backed_by_function"].values())
    assert len(g["focal_admissible_entries"]) <= \
        g["admissible_by_function"][g["focal_function"]]
    for variant, block in g["policy_variants"].items():
        assert sum(block["admissible_by_function"].values()) == block["total"], variant
        assert set(block["admissible_by_function"]) == set(g["function_vocabulary"])


# --------------------------------------------------------------------- ripp


def test_ripp_every_record_is_localised_and_typed_exactly_once():
    g = gold("ripp-precursor-mibig-4_0-01")
    records = g["selection_census"]["corpus_records"]
    assert len(g["localisation_table"]) == records
    assert g["precursor_length_summary"]["count"] == records
    for census in ("localisation_census", "status_census", "ripp_type_census",
                   "core_sequence_shape_census"):
        assert sum(g[census].values()) == records, census
    assert sum(g["core_sequence_case_census"].values()) <= records
    assert g["records_declaring_a_cleavage_location"] == len(g["cleavage_table"])
    assert sum(g["boundary_relation_census"].values()) == len(g["domain_table"])
    assert g["boundary_offset_summary"]["count"] == len(g["deletion_table"])
    assert sum(g["deletion_census"].values()) == records
    assert g["leader_length_summary"]["count"] <= records
    for policy, block in g["policy_sweep"].items():
        assert sum(block["verdict_census"].values()) == records, policy
        assert block["located"] <= block["cores_declared"], policy


# ----------------------------------------------------------------- coverage


def test_every_campaign_has_closure_assertions():
    """A new campaign must not join the suite without identities of its own.

    This file is per-campaign by necessity -- the generic version of it was
    measured and found vacuous -- so nothing but a coverage test keeps it from
    falling behind the suite it is about.
    """
    import re

    text = pathlib.Path(__file__).read_text()
    named = set(re.findall(r'gold\("([^"]+)"\)', text))
    named |= set(re.findall(r'"([a-z0-9][a-z0-9_.-]*-\d+)"', text))
    uncovered = sorted(p.name for p in ROOT.iterdir() if p.name not in named)
    assert not uncovered, f"no closure assertions for: {uncovered}"
