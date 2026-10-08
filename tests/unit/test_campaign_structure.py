"""Structure-features: UniProt per-residue annotations against mmCIF geometry.

Two ladders share one engine, so these tests pin the engine (mmCIF parsing,
numbering, selection, geometry) once and then pin what each ladder claims.

The engine tests matter more than usual here. Every number this template grades
is a function of the coordinate file and the pinned conventions, so a parser that
quietly drops the key-value form of a category, or a selection rule applied
inconsistently, produces a self-consistent and entirely wrong answer -- which no
amount of downstream checking would catch.
"""

from __future__ import annotations

import json
import math
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.structure_features.core import (
    Annotations,
    Rules,
    Structure,
    StructureError,
    centroid,
    pocket_volume,
    radius_of_gyration,
    read_cif,
    uniprot_positions,
    _tokens,
)
from npbench_c.templates.structure_features.params import StructureSpec
from npbench_c.templates.structure_features.pocket import analyse as pocket_analyse
from npbench_c.templates.structure_features.pocket.build import (
    _assert_has_content as pocket_assert,
)
from npbench_c.templates.structure_features.residues import analyse as residues_analyse
from npbench_c.templates.structure_features.residues.build import (
    _assert_has_content as residues_assert,
)

ROOT = pathlib.Path(__file__).resolve().parents[2] / "campaigns"
RESIDUES = ROOT / "residues-prna-01"
POCKET = ROOT / "pocket-rebh-01"

pytestmark = pytest.mark.skipif(
    not (RESIDUES / "gold" / "gold.json").is_file()
    or not (POCKET / "gold" / "gold.json").is_file(),
    reason="campaign gold not generated")


def _gold(campaign: pathlib.Path) -> dict:
    return json.loads((campaign / "gold" / "gold.json").read_text())


@pytest.fixture(scope="module")
def rules():
    return Rules.load(POCKET)


@pytest.fixture(scope="module")
def rebh():
    spec = StructureSpec.load(POCKET)
    return Structure.load(POCKET / spec.focal_input, spec.focal_pdb_id)


# ------------------------------------------------------------- the mmCIF parser

def test_tokeniser_keeps_quoted_values_whole():
    assert _tokens("A 'two words' B") == ["A", "two words", "B"]
    assert _tokens('X "a b" Y') == ["X", "a b", "Y"]
    # An apostrophe inside a token is not a quote: it does not open a string
    # because it is not at a token boundary, and C5' is a real atom name.
    assert _tokens("C5' N1 2.4") == ["C5'", "N1", "2.4"]
    assert _tokens("value # trailing comment") == ["value"]


def test_parser_reads_both_serialisations():
    """A category with one row is written as bare key-value lines, not a loop.

    2AQJ writes _struct_ref_seq that way and 2E4G writes it as a loop. A parser
    that handles only the loop form returns nothing for the first, and the
    numbering offset then gets assumed instead of read.
    """
    text = (
        "data_TEST\n"
        "_struct_ref_seq.pdbx_db_accession   P12345\n"
        "_struct_ref_seq.db_align_beg        1\n"
        "#\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.label_comp_id\n"
        "ATOM ALA\n"
        "ATOM GLY\n"
        "#\n"
    )
    blocks = read_cif(text, ("struct_ref_seq", "atom_site"))
    assert blocks["struct_ref_seq"] == [{"pdbx_db_accession": "P12345",
                                        "db_align_beg": "1"}]
    assert [r["label_comp_id"] for r in blocks["atom_site"]] == ["ALA", "GLY"]


def test_parser_reassembles_a_row_split_across_lines():
    text = ("loop_\n_t.a\n_t.b\n_t.c\n"
            "1 2\n3\n4 5 6\n#\n")
    rows = read_cif(text, ("t",))["t"]
    assert rows == [{"a": "1", "b": "2", "c": "3"}, {"a": "4", "b": "5", "c": "6"}]


# ------------------------------------------------------------------- numbering

def test_the_two_structure_numbering_axes_disagree(rebh):
    """The expression tag is the point of the pocket ladder's R1.

    label_seq_id starts at 21 where the UniProt sequence starts at 1, while the
    author numbering is aligned. A system that assumes either structure axis is
    the sequence position gets a coherent, wrong shell.
    """
    offsets = rebh.offsets("Q8KHZ8")
    assert offsets == {"label_to_uniprot": -20, "auth_to_uniprot": 0}
    assert offsets["label_to_uniprot"] != offsets["auth_to_uniprot"]


def test_numbering_refuses_a_per_chain_disagreement(rebh):
    """Two chains with different offsets would make "the" position ambiguous."""
    doctored = list(rebh.ref_seq)
    doctored[1] = {**doctored[1], "pdbx_auth_seq_align_beg": "7"}
    broken = Structure(pdb_id=rebh.pdb_id, atoms=rebh.atoms,
                       water_atoms=rebh.water_atoms,
                       excluded_hydrogens=rebh.excluded_hydrogens,
                       excluded_altloc=rebh.excluded_altloc,
                       ref_seq=tuple(doctored))
    with pytest.raises(StructureError, match="disagree"):
        broken.offsets("Q8KHZ8")


def test_numbering_refuses_an_accession_the_structure_does_not_cite(rebh):
    with pytest.raises(StructureError, match="no _struct_ref_seq record"):
        rebh.offsets("P95480")


# ------------------------------------------------------- the pinned conventions

def test_rules_refuse_an_unknown_eco_code(rules):
    with pytest.raises(StructureError, match="outside the pinned policy table"):
        rules.eco_class(["ECO:9999999"])


def test_rules_refuse_an_unknown_residue_and_element(rules):
    with pytest.raises(StructureError, match="residue alphabet"):
        rules.one_letter("ZZZ")
    with pytest.raises(StructureError, match="van der Waals"):
        rules.radius("Xx")


def test_evidence_class_is_the_strongest_citation(rules):
    """The census reports one class per record; the policy admits per citation.

    Both rules are declared, because computing one while reporting the other is
    how two blocks of the same report end up disagreeing.
    """
    assert rules.eco_class(["ECO:0000269", "ECO:0007744"]) == "experimental"
    assert rules.eco_class(["ECO:0007744"]) == "structural"
    assert rules.eco_class(["ECO:0000305"]) == "inferred"
    assert rules.eco_class([]) == "none"


def test_waters_hydrogens_and_alternate_locations_are_excluded(rebh):
    assert rebh.water_atoms > 0
    assert all(a.comp != "HOH" for a in rebh.atoms)
    assert all(a.element not in ("H", "D") for a in rebh.atoms)


# -------------------------------------------------------------------- geometry

def test_contacts_are_symmetric_between_the_two_chains():
    """The second chain is an internal control, not extra data."""
    gold = _gold(POCKET)
    symmetry = gold["shells"]["chain_symmetry"]["TRP"]
    assert symmetry["position_sets_identical"] is True
    assert symmetry["symmetric_difference"] == []
    assert len(symmetry["copies_compared"]) == 2


def test_contact_set_grows_monotonically_with_the_cutoff():
    sweep = _gold(POCKET)["counterfactuals"]["cutoff_sweep_positions"]
    ordered = [sweep[k] for k in sorted(sweep, key=float)]
    for tighter, looser in zip(ordered, ordered[1:]):
        assert set(tighter) <= set(looser)


def test_descriptors_are_reproducible_and_order_independent(rebh, rules):
    """The reductions are fixed-order, so reversing the atom list must not move
    a reported descriptor beyond its last declared decimal."""
    atoms = rebh.ligand_copies()[("TRP", "A", "1001")]
    forward, backward = centroid(atoms), centroid(tuple(reversed(atoms)))
    assert all(abs(a - b) <= 0.001 for a, b in zip(forward, backward))
    assert abs(radius_of_gyration(atoms)
               - radius_of_gyration(tuple(reversed(atoms)))) <= 0.001


def test_the_volume_lattice_accounts_for_every_point_it_visits(rebh, rules):
    atoms = rebh.ligand_copies()[("TRP", "A", "1001")]
    v = pocket_volume(atoms, rebh.polymer(), rules)
    assert (v["grid_points_near_ligand"]
            == v["grid_points_rejected_clash"] + v["grid_points_rejected_bulk"]
            + v["pocket_grid_points"])
    assert v["pocket_volume_a3"] == pytest.approx(
        v["pocket_grid_points"] * rules.grid_spacing ** 3, abs=1e-6)


def test_the_volume_lattice_is_anchored_on_integers(rebh, rules):
    """An integer-anchored lattice does not move when the box does, which is
    what makes the count reproducible rather than merely repeatable."""
    atoms = rebh.ligand_copies()[("TRP", "A", "1001")]
    lo = math.floor(min(a.x for a in atoms) - rules.cutoff)
    assert lo == float(int(lo))


def test_truncation_removes_only_the_named_side_chain(rebh, rules):
    polymer = rebh.polymer()
    kept, removed = rebh.truncated_polymer("A", "455", rules)
    assert removed == len(polymer) - len(kept) > 0
    survivors = {a.name for a in kept if a.chain == "A" and a.residue == "455"}
    assert survivors <= rules.retained_atoms
    # Nothing else moved or disappeared.
    assert {(a.chain, a.residue) for a in kept} == {(a.chain, a.residue)
                                                    for a in polymer}


def test_truncation_refuses_a_residue_with_no_side_chain(rebh, rules):
    glycine = next(a for a in rebh.polymer() if a.comp == "GLY")
    with pytest.raises(StructureError, match="no side chain"):
        rebh.truncated_polymer(glycine.chain, glycine.residue, rules)


def test_truncation_splits_the_shell_both_ways():
    """Half the shell reaches the ligand only through its side chain and half
    still reaches it through backbone or CB, so neither verdict is a safe guess
    and R4 cannot be cleared by a constant answer."""
    truncation = _gold(POCKET)["counterfactuals"]["truncation"]
    assert truncation["contact_survives"] and truncation["contact_lost"]
    assert (len(truncation["contact_survives"]) + len(truncation["contact_lost"])
            == truncation["shell_size"])


# ------------------------------------------------------- the residues ladder

def test_dedup_keeps_two_substitutions_at_one_position_apart():
    """Two mutagenesis records at position 272 differ only in the substitution,
    and the declared record identity includes it, so they stay two records."""
    rules = Rules.load(RESIDUES)
    spec = StructureSpec.load(RESIDUES)
    annotations = Annotations.load(RESIDUES / spec.annotation_input, rules)
    at_272 = [r for r in annotations.records
              if r.type == "Mutagenesis" and r.start == 272]
    assert len(at_272) == 2
    assert {r.substitution for r in at_272} == {"A", "F"}


def test_dedup_keeps_two_ligands_at_one_position_apart():
    rules = Rules.load(RESIDUES)
    spec = StructureSpec.load(RESIDUES)
    annotations = Annotations.load(RESIDUES / spec.annotation_input, rules)
    at_346 = [r for r in annotations.records
              if r.type == "Binding site" and r.start == 346]
    assert {r.ligand_id for r in at_346} == {"ChEBI:CHEBI:57912",
                                            "ChEBI:CHEBI:58713"}


def test_the_annotation_is_a_strict_subset_of_the_geometry():
    """The campaign's finding, and the reason R3 is not a restatement of R2:
    every annotated position lines the pocket, but most contacts are not
    annotated. A system that reports the annotation set as the contact set
    clears R2 and fails R3."""
    focal = _gold(RESIDUES)["confirmation"]
    extra = focal["contact_not_annotated_by_ligand"]
    assert any(extra.values())
    for chebi, annotated in focal["annotated_by_ligand"].items():
        contact = focal["contact_by_ligand"][chebi]
        assert set(annotated) <= set(contact)


def test_one_ligand_is_annotated_only_from_the_other_structure():
    """7-chloro-L-tryptophan cites 2AR8, not 2AQJ, and is absent from 2AQJ's
    coordinates. Its focal row is a true zero, and the R4 swap is what shows a
    system whether it respected WHICH entry evidences an annotation."""
    gold = _gold(RESIDUES)
    chebi = "ChEBI:CHEBI:58713"
    assert gold["confirmation"]["annotated_by_ligand"][chebi] == []
    assert gold["structure_swap"]["annotated_by_ligand"][chebi]
    deltas = gold["structure_swap"]["confirmed_delta_by_ligand"]
    assert deltas[chebi] > 0 and deltas["ChEBI:CHEBI:57912"] < 0


def test_evidence_axes_separate_the_unevidenced_from_the_inferred():
    gold = _gold(RESIDUES)
    assert gold["unevidenced_records"] == ["Site@79"]
    assert gold["inferred_only_records"] == ["Site@346"]
    assert gold["record_class_counts"]["Site"]["none"] == 1
    assert gold["record_class_counts"]["Site"]["inferred"] == 1


def test_every_policy_variant_admits_a_different_count():
    totals = _gold(RESIDUES)["policy_admitted_totals"]
    assert len(set(totals.values())) == len(totals)


# ---------------------------------------------------------------- the gates

@pytest.mark.parametrize("campaign", [RESIDUES, POCKET], ids=["residues", "pocket"])
def test_oracle_scores_exactly_one(campaign):
    m = json.loads((campaign / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((campaign / "grading.json").read_text())
    r = grade(m, _gold(campaign), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


@pytest.mark.parametrize("campaign", [RESIDUES, POCKET], ids=["residues", "pocket"])
def test_all_rungs_are_g1(campaign):
    spec = json.loads((campaign / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


def test_build_refuses_gold_when_the_geometry_adds_nothing():
    gold = _gold(RESIDUES)
    doctored = json.loads(json.dumps(gold))
    for chebi in doctored["confirmation"]["contact_not_annotated_by_ligand"]:
        doctored["confirmation"]["contact_not_annotated_by_ligand"][chebi] = []
    with pytest.raises(StructureError, match="outside its annotation set"):
        residues_assert(doctored)


def test_build_refuses_gold_when_the_swap_moves_nothing():
    doctored = json.loads(json.dumps(_gold(RESIDUES)))
    doctored["structure_swap"]["confirmed_delta_by_ligand"] = {
        k: 0 for k in doctored["structure_swap"]["confirmed_delta_by_ligand"]}
    with pytest.raises(StructureError, match="both"):
        residues_assert(doctored)


def test_build_refuses_gold_when_the_numbering_axes_agree():
    doctored = json.loads(json.dumps(_gold(POCKET)))
    doctored["numbering"]["label_to_uniprot"] = doctored["numbering"]["auth_to_uniprot"]
    with pytest.raises(StructureError, match="same offset"):
        pocket_assert(doctored)


# --------------------------------------------------- sandbox self-containment

@pytest.mark.parametrize("campaign", [RESIDUES, POCKET], ids=["residues", "pocket"])
def test_sandbox_carries_no_gold(campaign, tmp_path):
    s = build_sandbox(campaign, tmp_path / "cold")
    assert not (s / "gold").exists() and not (s / "grading.json").exists()
    assert sorted(p.name for p in s.iterdir()) == ["inputs", "reference",
                                                   "submission", "task.yaml"]


def test_residues_campaign_is_solvable_from_its_sandbox(tmp_path):
    """The gold is a pure function of the sandbox. Without this a campaign can
    depend on a table only the oracle can see, and every system fails for a
    reason that looks like incapacity."""
    s = build_sandbox(RESIDUES, tmp_path / "solvable")
    recomputed = residues_analyse.analyse(s)
    gold = _gold(RESIDUES)
    assert recomputed["confirmation"] == gold["confirmation"]
    assert recomputed["record_class_counts"] == gold["record_class_counts"]
    assert recomputed["policy_variants"] == gold["policy_variants"]


def test_pocket_campaign_is_solvable_from_its_sandbox(tmp_path):
    s = build_sandbox(POCKET, tmp_path / "solvable")
    recomputed = pocket_analyse.analyse(s)
    gold = _gold(POCKET)
    assert recomputed["shells"] == gold["shells"]
    assert recomputed["descriptors"] == gold["descriptors"]
    assert recomputed["counterfactuals"] == gold["counterfactuals"]


def test_pocket_r4_does_not_list_the_shell_it_perturbs(tmp_path):
    """The truncation sweep covers the focal shell, which is R3's answer. Naming
    the positions in the sandbox would hand that answer over, so the scope is
    declared by rule instead."""
    s = build_sandbox(POCKET, tmp_path / "scope")
    keys = json.loads((s / "reference" / "campaign_keys.json").read_text())
    assert keys["truncation_scope"] == "focal_shell"
    shell = set(_gold(POCKET)["descriptors"]["contact_positions_uniprot"])

    # Structural, not textual: a substring test would fire on "79" inside the
    # ChEBI identifier 57912. What must not be in the sandbox is a declared
    # NUMBER that is a shell position.
    def numbers(value):
        if isinstance(value, bool):
            return
        if isinstance(value, int):
            yield value
        elif isinstance(value, dict):
            for v in value.values():
                yield from numbers(v)
        elif isinstance(value, list):
            for v in value:
                yield from numbers(v)

    assert not (set(numbers(keys)) & shell)


@pytest.mark.parametrize("campaign", [RESIDUES, POCKET], ids=["residues", "pocket"])
def test_warm_bundles_hand_over_the_previous_rung_only(campaign, tmp_path):
    spec = json.loads((campaign / "grading.json").read_text())
    for rung in spec["rungs"][1:]:
        rid = rung["rung_id"]
        s = build_sandbox(campaign, tmp_path / f"warm_{rid}", f"gold/warm/{rid}")
        handed = [p for p in s.iterdir() if p.suffix == ".json"]
        assert handed, f"{rid} bundle is empty"
        blob = json.dumps([json.loads(p.read_text()) for p in handed])
        gold = _gold(campaign)

        # Value-based, not key-based: the R4 bundle legitimately carries the
        # FOCAL structure id under the key "structure" while R4 grades the
        # ALTERNATE one, so a key-name test would fire on a bundle that gives
        # nothing away. Scalars short enough to collide by accident -- a count
        # of 0, a one-letter residue -- are skipped; what matters is that no
        # composite answer is sitting in the bundle.
        for component in rung["components"]:
            gold_ref = component["args"].get("gold")
            if not isinstance(gold_ref, str) or not gold_ref.startswith("gold:"):
                continue
            value = gold
            for part in gold_ref.split(":", 1)[1].split("."):
                value = value[part]
            serialised = json.dumps(value, sort_keys=True)
            if len(serialised) <= 8:
                continue
            assert serialised not in blob, (
                f"{rid} bundle carries the answer to {component['name']}")


@pytest.mark.parametrize("campaign", [RESIDUES, POCKET], ids=["residues", "pocket"])
def test_pinned_tables_are_all_reachable(campaign, tmp_path):
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    s = build_sandbox(campaign, tmp_path / "tables")
    for table in task["constraints"]["pinned_tables"]:
        assert (s / table).is_file(), table
