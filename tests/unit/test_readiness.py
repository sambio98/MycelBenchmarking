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
    LICENSE_CLASSES,
    PASS,
    UNVALIDATED_DESIGN_LABELS,
    _contract_and_provenance,
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
    pair = ROOT / "constraint-conflict-ecoli-rebh-01"
    if not (pair / "task.yaml").is_file():
        pytest.skip("the constraint control pair is not built")
    campaign = tmp_path / "conflict"
    shutil.copytree(pair, campaign)
    # The copy keeps the pairing, but it is no longer the campaign its partner
    # names -- which is exactly the asymmetry the check looks for.
    (tmp_path / "constraint-feasible-ecoli-rebh-01").mkdir()
    shutil.copy(ROOT / "constraint-feasible-ecoli-rebh-01" / "task.yaml",
                tmp_path / "constraint-feasible-ecoli-rebh-01" / "task.yaml")
    result = _design_claims_are_checkable(campaign)
    assert result.status != PASS
    assert "names" in result.detail


def test_a_campaign_naming_its_real_partner_passes(tmp_path):
    """The other side of the pairing test: a symmetric pair is fine."""
    for name in ("constraint-conflict-ecoli-rebh-01",
                 "constraint-feasible-ecoli-rebh-01"):
        source = ROOT / name
        if not (source / "task.yaml").is_file():
            pytest.skip("the constraint control pair is not built")
        (tmp_path / name).mkdir()
        shutil.copy(source / "task.yaml", tmp_path / name / "task.yaml")
    result = _design_claims_are_checkable(
        tmp_path / "constraint-conflict-ecoli-rebh-01")
    assert result.status == PASS, result.detail
