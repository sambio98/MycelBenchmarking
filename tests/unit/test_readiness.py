"""The readiness gate's own checks, in particular the warm-bundle check.

A rung may not be handed the answer to one of its own components. The failure is
invisible in an oracle run -- the oracle scores 1.0 either way -- and only shows
up when a real system clears a rung it should not have, by which point the
per-rung difficulty profile is fiction. So it needs a mechanical check, and the
check needs tests of its own: one that it fires, and one that it does not fire on
the scalar vocabulary members that legitimately appear in every census.
"""

from __future__ import annotations

import json
import pathlib
import shutil

import pytest

from npbench_c.readiness.gate import (
    ANTISMASH_PROCESSING_KEY,
    LICENSE_CLASSES,
    PASS,
    UNVALIDATED_DESIGN_LABELS,
    _antismash_versions_are_not_mixed,
    _contract_and_provenance,
    _contract_names_what_is_graded,
    _design_claims_are_checkable,
    _resolve_gold,
    _warm_bundles_withhold_answers,
)

ROOT = pathlib.Path(__file__).resolve().parents[2] / "campaigns"
CAMPAIGNS = sorted(p for p in ROOT.iterdir()
                   if (p / "gold" / "gold.json").is_file())
#: Tiny inputs, so copying it for a doctored case costs nothing.
SMALL = ROOT / "construct-ecoli-rebh-01"

pytestmark = pytest.mark.skipif(not CAMPAIGNS, reason="no campaign gold generated")


def test_gold_paths_all_resolve():
    """An unresolvable gold path is a spec bug that would silently score 0."""
    for campaign in CAMPAIGNS:
        gold = json.loads((campaign / "gold" / "gold.json").read_text())
        spec = json.loads((campaign / "grading.json").read_text())
        for rung in spec["rungs"]:
            for component in rung["components"]:
                reference = component["args"].get("gold")
                if not isinstance(reference, str) or not reference.startswith("gold:"):
                    continue
                _resolve_gold(gold, reference.split(":", 1)[1])


@pytest.mark.parametrize("campaign", CAMPAIGNS, ids=lambda p: p.name)
def test_no_campaign_hands_a_rung_its_own_answer(campaign):
    check = _warm_bundles_withhold_answers(campaign)
    assert check.status == PASS, check.detail
    # The check must actually be looking at something on every campaign.
    assert int(check.detail.split(None, 1)[0]) > 0, check.detail


def test_the_check_fires_when_a_bundle_carries_a_composite_answer(tmp_path):
    copy = tmp_path / SMALL.name
    shutil.copytree(SMALL, copy)
    spec = json.loads((copy / "grading.json").read_text())
    gold = json.loads((copy / "gold" / "gold.json").read_text())

    rung = next(r for r in spec["rungs"] if r["rung_id"] == "r4")
    planted = None
    for component in rung["components"]:
        reference = component["args"].get("gold")
        if not isinstance(reference, str) or not reference.startswith("gold:"):
            continue
        value = _resolve_gold(gold, reference.split(":", 1)[1])
        if isinstance(value, (dict, list)):
            planted = (component["name"], value)
            break
    assert planted, "r4 has no composite component to plant"

    (copy / "gold" / "warm" / "r4" / "leaked.json").write_text(
        json.dumps({"handed_over": planted[1]}, sort_keys=True))
    check = _warm_bundles_withhold_answers(copy)
    assert check.status != PASS
    assert f"r4:{planted[0]}" in check.detail


def test_the_check_ignores_a_scalar_from_a_declared_vocabulary():
    """The deliberate exemption, pinned so nobody "tightens" it later.

    The construct ladder grades the binding constraint of the swapped host as an
    enum. Its value, 'direct_repeat', is a vocabulary member and appears in the
    R4 bundle as a key of the first host's violation census -- which gives
    nothing away, because the question is which constraint binds for the OTHER
    host. Flagging a scalar would make this a false positive, and would do the
    same in every enum-scored component in the benchmark.
    """
    import tempfile

    from npbench_c.sweep.runner import build_sandbox

    spec = json.loads((SMALL / "grading.json").read_text())
    gold = json.loads((SMALL / "gold" / "gold.json").read_text())
    component = next(c for r in spec["rungs"] if r["rung_id"] == "r4"
                     for c in r["components"] if c["scorer"] == "enum_exact")
    value = _resolve_gold(gold, component["args"]["gold"].split(":", 1)[1])
    assert isinstance(value, str) and value in component["args"]["vocabulary"]

    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(SMALL, pathlib.Path(tmp) / "r4", "gold/warm/r4")
        blob = json.dumps([json.loads(p.read_text())
                           for p in sorted(sandbox.iterdir())
                           if p.suffix == ".json"])
    assert value in blob, "the premise of the exemption no longer holds"
    assert _warm_bundles_withhold_answers(SMALL).status == PASS


def test_the_compound_list_is_no_longer_graded():
    """What the check caught when it was first run across the benchmark.

    The compound list is a field of the shipped entry and every warm bundle hands
    it over, because neither the assembly rung nor the reconciliation rung can be
    posed without it. Grading it measured a parse the sandbox had already
    answered.
    """
    import yaml

    for name in ("massbalance-nrps-malleobactin-01", "massbalance-nrps-sevadicin-01"):
        campaign = ROOT / name
        spec = json.loads((campaign / "grading.json").read_text())
        names = {c["name"] for r in spec["rungs"] for c in r["components"]}
        assert "compounds_enumerated" not in names
        task = yaml.safe_load((campaign / "task.yaml").read_text())
        reasons = {e["field"]: e["reason"] for e in task["excluded_from_grading"]}
        assert any("compound list" in f for f in reasons)


# ------------------------------------------------------------- licence classes

def test_share_alike_is_a_declarable_class():
    """ChEMBL is CC BY-SA 3.0. The decision taken was to allow it, because the
    ShareAlike obligation attaches to an ADAPTATION -- a derived subset and the
    gold computed from it -- and not to a COLLECTION of separable campaigns. So
    the benchmark's code, grader and other campaigns are unaffected, and a
    share-alike campaign can be dropped without touching anything else."""
    assert set(LICENSE_CLASSES) == {"open", "nc", "share_alike"}


@pytest.mark.parametrize("campaign", CAMPAIGNS, ids=lambda p: p.name)
def test_every_campaign_declares_an_allowed_licence_class(campaign):
    import yaml

    task = yaml.safe_load((campaign / "task.yaml").read_text())
    assert task["license_class"] in LICENSE_CLASSES


def test_a_share_alike_campaign_must_carry_its_notice(tmp_path):
    """A licence condition recorded only in a design document is a condition the
    person redistributing the files will never see, so the notice travels with the
    campaign and the gate checks it."""
    import shutil

    import yaml

    source = next(c for c in CAMPAIGNS if (c / "task.yaml").is_file())
    copy = tmp_path / source.name
    shutil.copytree(source, copy)
    task = yaml.safe_load((copy / "task.yaml").read_text())
    task["license_class"] = "share_alike"
    # The notice is removed explicitly rather than assumed absent: whichever
    # campaign sorts first is incidental, and one of them now carries a notice for
    # its own reasons, which quietly turned this into a test of nothing.
    task.pop("license_notice", None)
    (copy / "task.yaml").write_text(yaml.safe_dump(task, sort_keys=False))

    by_name = {c.name: c for c in _contract_and_provenance(copy)}
    assert by_name["license_class_declared"].status == PASS
    assert by_name["share_alike_notice_present"].status != PASS

    task["license_notice"] = {
        "source": "ChEMBL", "license": "CC BY-SA 3.0",
        "attribution": "ChEMBL, EMBL-EBI",
        "redistribution": "redistribute these files under CC BY-SA 3.0",
    }
    (copy / "task.yaml").write_text(yaml.safe_dump(task, sort_keys=False))
    by_name = {c.name: c for c in _contract_and_provenance(copy)}
    assert by_name["share_alike_notice_present"].status == PASS


# ----------------------------------------------- design claims a reader can check


@pytest.mark.parametrize("campaign", CAMPAIGNS, ids=lambda c: c.name)
def test_no_campaign_declares_an_unvalidated_design_label(campaign):
    """What retired `structural_elements`, asserted over the whole suite.

    It was a free-text list of design features that nothing read: 23 of 23
    campaigns declared `verification` and `counterfactual`, which is what the
    four-rung ladder means; 11 declared `planted` and no campaign ships an
    altered input; 14 declared `control` for three unrelated things. And on the
    constraint control pair the conflict half declared `abstention` while the
    feasible half did not, which states which half of a two-verdict pair the
    reader is holding.
    """
    result = _design_claims_are_checkable(campaign)
    assert result.status == PASS, result.detail


def test_the_check_fires_on_a_restored_design_label(tmp_path):
    """The check has to be able to fail, under any of the field's names."""
    for field in UNVALIDATED_DESIGN_LABELS:
        campaign = tmp_path / field
        shutil.copytree(SMALL, campaign)
        task = campaign / "task.yaml"
        task.write_text(task.read_text()
                        + f"\n{field}: [verification, counterfactual, planted]\n")
        result = _design_claims_are_checkable(campaign)
        assert result.status != PASS, field
        assert field in result.detail


def test_the_check_fires_on_a_dangling_campaign_reference(tmp_path):
    """A renamed or retired campaign leaves a sentence pointing at nothing."""
    campaign = tmp_path / SMALL.name
    shutil.copytree(SMALL, campaign)
    task = campaign / "task.yaml"
    task.write_text(task.read_text() + (
        "\nsplit_note: >\n  compare against construct-ecoli-retired-09, which "
        "shares this template.\n"))
    result = _design_claims_are_checkable(campaign)
    assert result.status != PASS
    assert "construct-ecoli-retired-09" in result.detail


def test_the_check_fires_on_a_one_way_control_pair(tmp_path):
    """A control pair whose halves do not name each other is two campaigns."""
    pair = ROOT / "constraint-gc62-ecoli-rebh-01"
    if not (pair / "task.yaml").is_file():
        pytest.skip("the constraint control pair is not built")
    campaign = tmp_path / "conflict"
    shutil.copytree(pair, campaign)
    # The copy keeps the pairing, but it is no longer the campaign its partner
    # names -- which is exactly the asymmetry the check looks for.
    (tmp_path / "constraint-gc55-ecoli-rebh-01").mkdir()
    shutil.copy(ROOT / "constraint-gc55-ecoli-rebh-01" / "task.yaml",
                tmp_path / "constraint-gc55-ecoli-rebh-01" / "task.yaml")
    result = _design_claims_are_checkable(campaign)
    assert result.status != PASS
    assert "names" in result.detail


def test_a_campaign_naming_its_real_partner_passes(tmp_path):
    """The other side of the pairing test: a symmetric pair is fine."""
    for name in ("constraint-gc62-ecoli-rebh-01",
                 "constraint-gc55-ecoli-rebh-01"):
        source = ROOT / name
        if not (source / "task.yaml").is_file():
            pytest.skip("the constraint control pair is not built")
        (tmp_path / name).mkdir()
        shutil.copy(source / "task.yaml", tmp_path / name / "task.yaml")
    result = _design_claims_are_checkable(
        tmp_path / "constraint-gc62-ecoli-rebh-01")
    assert result.status == PASS, result.detail


# ------------------------------------------------- one antiSMASH per campaign


def _fake_campaign(root: pathlib.Path, snapshots: dict,
                   provenance: dict | None = None) -> pathlib.Path:
    """The two files check 25 reads, and nothing else."""
    campaign = root / "campaign"
    (campaign / "inputs").mkdir(parents=True, exist_ok=True)
    (campaign / "task.yaml").write_text(
        json.dumps({"snapshots_required": snapshots}, indent=2) + "\n")
    if provenance is not None:
        (campaign / "inputs" / "provenance.json").write_text(
            json.dumps(provenance, indent=2) + "\n")
    return campaign


PUBLISHED = ("https://dl.secondarymetabolites.org/mibig/"
             "mibig_antismash_4.0_gbk_as8b1.tar.bz2")


@pytest.mark.parametrize("campaign", CAMPAIGNS, ids=lambda c: c.name)
def test_no_campaign_mixes_two_antismash_versions(campaign):
    """The measurement that resolved the open item, asserted over the suite.

    Four campaigns ship records out of MIBiG's antiSMASH 8.0-beta-1 set and one
    campaign runs the pinned antiSMASH 8.0.4. The two groups do not overlap, so
    nothing is mixed -- and this is what notices if that changes.
    """
    result = _antismash_versions_are_not_mixed(campaign)
    assert result.status == PASS, result.detail


def test_processed_records_without_a_version_record_are_refused(tmp_path):
    campaign = _fake_campaign(
        tmp_path, {"shipped": ["mibig_antismash@4.0_as8b1"], "image": []},
        {"records": {"source": "MIBiG antiSMASH-processed reference set",
                     "url": PUBLISHED}})
    result = _antismash_versions_are_not_mixed(campaign)
    assert result.status != PASS
    assert ANTISMASH_PROCESSING_KEY in result.detail


def test_running_one_antismash_over_another_s_records_is_refused(tmp_path):
    """The forbidden comparison, planted: the pinned tool against the published
    set. The refusal names both ways out, because either is legitimate."""
    campaign = _fake_campaign(
        tmp_path,
        {"shipped": ["mibig_antismash@4.0_as8b1"], "image": ["antismash@8.0.4"]},
        {"records": {"source": "MIBiG antiSMASH-processed reference set",
                     "url": PUBLISHED,
                     ANTISMASH_PROCESSING_KEY: "8.0 beta 1, as the archive name "
                                               "as8b1 states"}})
    result = _antismash_versions_are_not_mixed(campaign)
    assert result.status != PASS
    assert "cross-version" in result.detail
    assert "reprocess" in result.detail and "published set" in result.detail


def test_a_reprocessed_set_is_not_refused(tmp_path):
    """The other way out has to pass, or the check would forbid the fix.

    A campaign that reprocessed the records with the pinned antiSMASH and says so
    is running one version over its own output, which is what the rule asks for.
    """
    campaign = _fake_campaign(
        tmp_path,
        {"shipped": ["mibig_antismash@4.0_reprocessed"],
         "image": ["antismash@8.0.4"]},
        {"records": {"source": "MIBiG reference set, reprocessed here",
                     "url": PUBLISHED,
                     ANTISMASH_PROCESSING_KEY: "8.0.4, the version this image "
                                               "pins; reprocessed for this "
                                               "campaign"}})
    assert _antismash_versions_are_not_mixed(campaign).status == PASS


def test_a_pin_the_registry_does_not_carry_is_refused(tmp_path):
    campaign = _fake_campaign(tmp_path, {"shipped": [], "image": ["antismash@7.1.0"]})
    result = _antismash_versions_are_not_mixed(campaign)
    assert result.status != PASS
    assert "registry pins" in result.detail


def test_a_note_mentioning_antismash_is_not_shipped_processed_data(tmp_path):
    """Why the scan reads `source` and `url` and never a note.

    The campaign that RUNS antiSMASH explains in its provenance that it ships the
    complete chromosome because a slice truncates CDS features "which antiSMASH
    refuses". A scan over prose reads that as shipping processed records and
    fails the one campaign doing it right.
    """
    campaign = _fake_campaign(
        tmp_path, {"shipped": ["genbank@AL645882.2"], "image": ["antismash@8.0.4"]},
        {"sequence_source": "NCBI Nucleotide efetch, rettype=gbwithparts",
         "note": "The COMPLETE chromosome, unsliced. A slice truncates CDS "
                 "features at its edges, which antiSMASH refuses."})
    assert _antismash_versions_are_not_mixed(campaign).status == PASS


# --- check 27: the published contract against what the measurer reads --------
#
# The failure this covers was found by projecting each oracle submission onto the
# contract's own field names: five of twenty-five campaigns graded lower that way
# than the oracle did, two of them at depth 0. The oracle never goes through the
# contract -- `solve.py` renders from gold's key tuples -- so nothing in the gate
# had ever compared the two.


@pytest.mark.parametrize("campaign", CAMPAIGNS, ids=lambda p: p.name)
def test_every_campaign_can_be_answered_from_its_published_contract(campaign):
    check = _contract_names_what_is_graded(campaign)
    assert check.status == PASS, check.detail


def test_a_contract_naming_a_field_the_measurer_does_not_read_is_refused(tmp_path):
    """The arch failure mode: a template generalised, its contract left behind."""
    copy = tmp_path / SMALL.name
    shutil.copytree(SMALL, copy)
    task = (copy / "task.yaml").read_text()
    assert "      binding_constraint: >" in task
    (copy / "task.yaml").write_text(
        task.replace("      binding_constraint: >", "      which_constraint_binds: >"))

    check = _contract_names_what_is_graded(copy)
    assert check.status != PASS
    assert "r3.binding_constraint written, not declared" in check.detail
    assert "r3.which_constraint_binds declared, not written" in check.detail


def test_a_contract_silent_about_a_wrapper_key_is_refused(tmp_path):
    """The construct failure mode: the answer lives a level down from the names.

    Removing the `r3` and `r4` keys from the contract and promoting their fields
    is exactly how both construct campaigns shipped: every field name correct,
    and a system following them writes all of them one level too high.
    """
    import yaml

    copy = tmp_path / SMALL.name
    shutil.copytree(SMALL, copy)
    task = yaml.safe_load((copy / "task.yaml").read_text())
    analysis = task["output_contract"]["analysis.json"]
    task["output_contract"]["analysis.json"] = {
        "conventions": analysis["conventions"],
        **analysis["r3"], **analysis["r4"],
    }
    (copy / "task.yaml").write_text(yaml.safe_dump(task, sort_keys=False))

    check = _contract_names_what_is_graded(copy)
    assert check.status != PASS
    assert "analysis.json.r3 written, not declared" in check.detail
    # R1 and R2 are read off the ORF files, which the promotion leaves alone, so
    # the ladder stops exactly where the analytical claims start.
    assert "depth=2" in check.detail


def test_conventions_are_not_read_as_fields(tmp_path):
    """A format, a closed vocabulary and a FASTA's properties are not fields.

    Without the reserved key the check would demand that the oracle write a
    `formula_format` key, and the only way to satisfy it would be to grade one.
    """
    import yaml

    copy = tmp_path / SMALL.name
    shutil.copytree(SMALL, copy)
    task = yaml.safe_load((copy / "task.yaml").read_text())
    conventions = task["output_contract"]["analysis.json"].pop("conventions")
    task["output_contract"]["analysis.json"].update(conventions)
    (copy / "task.yaml").write_text(yaml.safe_dump(task, sort_keys=False))

    check = _contract_names_what_is_graded(copy)
    assert check.status != PASS
    assert "analysis.json.constraint_vocabulary declared, not written" in check.detail


# --- the gold-side inputs a build may read ----------------------------------


def test_every_declared_gold_side_input_exists():
    """A declaration naming a file nothing ships would waive a check silently."""
    from npbench_c.readiness.gate import GOLD_SIDE_INPUTS, gold_side_inputs

    assert GOLD_SIDE_INPUTS.is_file(), GOLD_SIDE_INPUTS
    declared = json.loads(GOLD_SIDE_INPUTS.read_text())["inputs"]
    for campaign, entries in declared.items():
        assert (ROOT / campaign).is_dir(), campaign
        for name, reason in gold_side_inputs(campaign).items():
            assert (ROOT / campaign / "gold" / name).is_file(), (campaign, name)
            assert len(reason) > 80, (campaign, name)


def test_no_graded_component_reads_a_gold_side_input():
    """The point of the declaration: these files are not part of any answer.

    A graded value derived from the queries' own identities would be a component
    no system could reach, and the determinism check restores these files before
    rebuilding -- so if one of them fed a graded path, that check would be
    covering for an unanswerable campaign rather than catching it.
    """
    from npbench_c.readiness.gate import GOLD_SIDE_INPUTS, gold_side_inputs

    declared = json.loads(GOLD_SIDE_INPUTS.read_text())["inputs"]
    for campaign in declared:
        spec = json.loads((ROOT / campaign / "grading.json").read_text())
        paths = [arg for rung in spec["rungs"] for component in rung["components"]
                 for arg in component["args"].values()
                 if isinstance(arg, str) and arg.startswith("gold:")]
        assert paths, campaign
        for name in gold_side_inputs(campaign):
            stem = pathlib.Path(name).stem            # queries_truth, expectation
            root = stem.split("_")[0]                 # queries, expectation
            for path in paths:
                assert stem not in path, (campaign, path)
                assert not path.startswith(f"gold:{root}"), (campaign, path)


def test_an_unneeded_gold_side_declaration_is_refused(monkeypatch):
    """The staleness half, without paying for two real gold rebuilds.

    A campaign that passes says so: "the sandbox plus expectation.json, each of
    them required" is only reachable when the rebuild WITHOUT the declared file
    failed. This covers the other branch -- a declared file the build does not
    need, which would waive the derivability question for that campaign.
    """
    from npbench_c.readiness import gate

    monkeypatch.setattr(gate, "gold_side_inputs",
                        lambda name, path=None: {"expectation.json": "planted"})
    monkeypatch.setattr(gate, "_rebuild_gold",
                        lambda campaign, into, restore: (True, "byte-identical"))
    check = gate._oracle_determinism(SMALL)
    assert check.status != PASS
    assert "are not needed" in check.detail
