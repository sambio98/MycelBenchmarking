"""Chemical space: how many compounds a corpus contains depends on the key.

The campaign's subject is the identity key, so these tests pin the three keys
apart, pin the counting rules that would otherwise drift, and pin the build's
refusal to emit gold where a distinction does not bite.

The scaffold tests matter most. Bemis-Murcko reduction discards substituents, so
the most widely shared scaffold in this corpus is plain benzene; a chemical-space
analysis that reports that number without saying what it did about it is
reporting an artefact of the reduction. The informativeness filter is the
declared response and the ladder grades both sides of it.
"""

from __future__ import annotations

import json
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.chemical_space import analyse as chem
from npbench_c.templates.chemical_space.build import _assert_has_content
from npbench_c.templates.chemical_space.core import (
    TABLE,
    ChemicalSpace,
    SpaceError,
    connectivity_key,
)

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "chemspace-mibig-4_0-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def space():
    params = chem.load_params(CAMPAIGN)
    return ChemicalSpace.load(CAMPAIGN, params["corpus_archive"])


@pytest.fixture(scope="module")
def table():
    return json.loads((CAMPAIGN / "reference" / TABLE).read_text())


# ---------------------------------------------------------------- the gates

def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_all_rungs_are_g1():
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


# -------------------------------------------------------- the identity keys

def test_the_three_keys_give_three_strictly_different_counts():
    """The campaign's premise. If any two agreed, the rung that distinguishes
    them would be grading the same statistic twice."""
    keys = _gold()["distinct_keys"]
    assert keys["compound"] > keys["connectivity"] > keys["scaffold"]


def test_the_connectivity_key_merges_real_variants(space):
    """Stereoisomers, charge states and labelled variants. If this were zero the
    connectivity granularity would be a distinction without a difference."""
    gold = _gold()
    assert gold["connectivity_groups_with_multiple_compound_keys"] > 0
    assert (gold["distinctions_lost_at_connectivity"]
            == gold["distinct_compound_keys"] - gold["distinct_connectivity_keys"])
    assert (gold["distinctions_lost_at_connectivity"]
            >= gold["connectivity_groups_with_multiple_compound_keys"])


def test_connectivity_key_is_the_first_inchikey_block(space):
    record = space.records[0]
    assert connectivity_key(record["inchikey"]) == record["inchikey"].split("-")[0]
    assert len(connectivity_key(record["inchikey"])) == 14


def test_scaffold_smiles_is_not_the_key(table):
    """Canonical SMILES is a toolkit output, so it cannot identify anything in a
    benchmark whose scores must stay comparable. The two counts differ, which is
    exactly why the key has to be declared rather than left to the solver."""
    records = table["records"]
    smiles = {r["scaffold_smiles"] for r in records if r["scaffold_smiles"]}
    inchikeys = {r["scaffold_inchikey"] for r in records if r["scaffold_inchikey"]}
    assert len(smiles) != len(inchikeys)
    assert len(inchikeys) == _gold()["distinct_scaffold_keys"]


def test_record_id_is_accession_and_ordinal_not_the_name(table):
    """Two entries in this release list the same compound name twice. Keying on
    the name would silently drop records and shift every count below it."""
    records = table["records"]
    assert len({r["record_id"] for r in records}) == len(records)
    by_name = {(r["accession"], r["name"]) for r in records}
    assert len(by_name) < len(records)


# --------------------------------------------------------- scaffold reduction

def test_acyclic_compounds_are_excluded_from_scaffold_statistics(space, table):
    gold = _gold()
    acyclic = [r for r in table["records"] if r["scaffold_inchikey"] is None]
    assert len(acyclic) == gold["acyclic_compounds"] > 0
    assert gold["scaffold_ring_count_histogram"]["0"] == gold["acyclic_compounds"]
    assert all(r["scaffold_smiles"] is None for r in acyclic)


def test_the_unfiltered_statistic_is_dominated_by_bare_rings(space):
    """The campaign's reason for existing, asserted rather than described.

    Over this corpus the most widely shared Bemis-Murcko scaffold is benzene. Any
    chemical-space analysis that reports 'these compounds share a scaffold'
    without saying what it did about that is reporting an artefact.
    """
    labels = {r["scaffold_inchikey"]: r["scaffold_smiles"] for r in space.records
              if r["scaffold_inchikey"]}
    shared = space.shared_scaffolds(space.records)
    top = min(shared.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    assert labels[top[0]] == "c1ccccc1"
    assert len(top[1]) > 50


def test_the_filter_removes_bare_rings_and_keeps_real_frameworks(space):
    informative = space.informative(space.records)
    labels = {r["scaffold_inchikey"]: r["scaffold_smiles"] for r in space.records
              if r["scaffold_inchikey"]}
    shared = space.shared_scaffolds(informative)
    top = min(shared.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    assert labels[top[0]] != "c1ccccc1"
    # Both floors must actually be doing work.
    assert all(r["scaffold_ring_count"] >= space.rules.min_ring_count
               for r in informative)
    assert all(r["scaffold_heavy_atoms"] / r["compound_heavy_atoms"]
               >= space.rules.min_coverage for r in informative)


def test_the_filter_shrinks_the_sharing_claim():
    shared = _gold()["shared_scaffolds"]
    cross = _gold()["cross_class_scaffolds"]
    assert shared["informative"] < shared["all"]
    assert cross["informative"] < cross["all"]


# ------------------------------------------------------------ the class join

def test_the_naive_cross_class_rule_is_badly_inflated():
    """456 active entries carry more than one class, so under the naive rule a
    single hybrid cluster makes its scaffold cross classes with nothing to
    compare it against."""
    cross = _gold()["cross_class_scaffolds"]
    assert cross["naive_all"] > 2 * cross["all"]
    assert cross["naive_informative"] > cross["informative"]
    assert _gold()["multi_class_entries"] > 0


def test_a_cross_class_scaffold_really_has_two_disjoint_entries(space):
    """The declared rule, checked against the classes rather than trusted."""
    for scaffold, accessions in _gold()["cross_class_informative_entries"].items():
        sets = [space.classes[a] for a in accessions]
        assert any(not (a & b) for i, a in enumerate(sets) for b in sets[i + 1:]), \
            scaffold


def test_per_class_counts_sum_above_the_corpus(space):
    """Deliberate: a record contributes to every class its entry carries, because
    MIBiG supplies no rule for choosing a primary class."""
    gold = _gold()
    assert sum(gold["records_per_class"].values()) > gold["compound_records"]
    assert set(gold["records_per_class"]) == set(gold["class_vocabulary"])


# ------------------------------------------------------- the counterfactuals

def test_every_threshold_variant_admits_a_different_record_count():
    records = _gold()["threshold_variant_records"]
    assert len(set(records.values())) == len(records)


def test_the_unfiltered_variant_reproduces_the_unfiltered_statistic():
    """t1 sets both floors to their no-op values, so it must agree with the
    all-scaffold numbers. An anchor against an off-by-one in the filter."""
    gold = _gold()
    t1 = gold["threshold_variants"]["t1_unfiltered"]
    assert t1["informative_records"] == gold["compound_records"] - gold["acyclic_compounds"]
    assert t1["shared_informative"] == gold["shared_scaffolds"]["all"]
    assert t1["cross_class_informative"] == gold["cross_class_scaffolds"]["all"]


def test_the_evidence_gate_shrinks_everything(space):
    gold, gate = _gold(), _gold()["evidence_gate"]
    assert gate["entries"] < gold["entries_with_compounds"]
    assert gate["compound_records"] < gold["compound_records"]
    assert (gate["cross_class_scaffolds"]["informative"]
            < gold["cross_class_scaffolds"]["informative"])
    assert len(gate["cross_class_informative_keys"]) == \
        gate["cross_class_scaffolds"]["informative"]


def test_the_gate_reads_locus_evidence_not_compound_evidence(space):
    """The same axis the MIBiG audit campaign pins apart. Reading
    compounds[].evidence here would gate on how a structure was determined rather
    than on how the cluster was linked to its product."""
    gated = {r["accession"] for r in space.evidence_gated(space.records)}
    assert gated
    for accession in sorted(gated)[:50]:
        assert space.locus_evidence[accession] & space.rules.accepted_evidence


# -------------------------------------------------------------- build refusal

@pytest.mark.parametrize("mutate,message", [
    (lambda g: g.__setitem__("distinctions_lost_at_connectivity", 0),
     "merges nothing"),
    (lambda g: g["cross_class_scaffolds"].__setitem__(
        "naive_all", g["cross_class_scaffolds"]["all"]), "naive cross-class"),
    (lambda g: g["shared_scaffolds"].__setitem__(
        "informative", g["shared_scaffolds"]["all"]), "removes no shared scaffold"),
    (lambda g: g["evidence_gate"]["cross_class_scaffolds"].__setitem__(
        "informative", g["cross_class_scaffolds"]["informative"]),
     "does not shrink"),
])
def test_build_refuses_gold_without_content(mutate, message):
    doctored = json.loads(json.dumps(_gold()))
    mutate(doctored)
    with pytest.raises(SpaceError, match=message):
        _assert_has_content(doctored)


def test_class_outside_the_vocabulary_is_refused(space):
    with pytest.raises(SpaceError, match="not a declared identity key"):
        space.key(space.records[0], "inchi")


# --------------------------------------------------- sandbox self-containment

def test_sandbox_carries_no_gold(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "cold")
    assert sorted(p.name for p in s.iterdir()) == ["inputs", "reference",
                                                   "submission", "task.yaml"]


def test_campaign_is_solvable_from_its_sandbox(tmp_path):
    """Gold is a pure function of the sandbox. The pinned chemistry table is an
    INPUT -- it supplies scaffold perception, which needs a toolkit -- and every
    graded number is bookkeeping over it that the agent must do itself."""
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    recomputed = chem.analyse(s)
    gold = _gold()
    for key in ("distinct_compound_keys", "cross_entry_duplicates",
                "shared_scaffolds", "cross_class_scaffolds",
                "cross_class_informative_entries", "threshold_variants",
                "evidence_gate"):
        assert recomputed[key] == gold[key], key


def test_the_chemistry_table_is_per_record_chemistry_and_provenance_only(tmp_path):
    """No aggregate the campaign grades as a claim.

    The table's only non-record fields are provenance -- which release, which
    RDKit, which scaffold definition, and a record count that is a self-check on
    the file. Everything the ladder grades is an aggregate, a set or a join that
    is not in here, which is what makes the table an input rather than an answer.
    """
    s = build_sandbox(CAMPAIGN, tmp_path / "table")
    table = json.loads((s / "reference" / TABLE).read_text())
    assert set(table) == {"description", "corpus_release", "corpus_filter",
                          "rdkit_version", "scaffold_definition", "record_count",
                          "records"}
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    assert not (graded & set(table)), graded & set(table)
    assert set(table["records"][0]) == {
        "record_id", "accession", "ordinal", "name", "formula", "inchikey",
        "compound_heavy_atoms", "compound_ring_count", "scaffold_smiles",
        "scaffold_inchikey", "scaffold_heavy_atoms", "scaffold_ring_count"}


def test_pinned_tables_are_all_reachable(tmp_path):
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    s = build_sandbox(CAMPAIGN, tmp_path / "tables")
    for pinned in task["constraints"]["pinned_tables"]:
        assert (s / pinned).is_file(), pinned


def test_catalog_deviation_is_recorded():
    """The catalogued R3 and R4 are not what this campaign does, and the artifact
    must say so rather than the difference living only in a commit message."""
    task = yaml.safe_load((CAMPAIGN / "task.yaml").read_text())
    assert "catalog_deviation" in task
    deviation = task["catalog_deviation"]
    assert "distinguishing substitution" in deviation and "S1" in deviation
    excluded = {e["field"] for e in task["excluded_from_grading"]}
    assert any("distinguishing substitution" in f for f in excluded)
    assert any("scaffold_smiles" in f for f in excluded)
