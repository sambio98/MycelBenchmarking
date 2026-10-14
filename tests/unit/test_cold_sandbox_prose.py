"""The cold sandbox must not state a graded ANSWER, in any shape.

The scan itself lives in `npbench_c.readiness.cold_sandbox` and runs as the
readiness gate's check 22, so this file is deliberately thin: it exercises that
one implementation across every built campaign at once, which the gate does only
one campaign at a time and only when someone runs it.

Why the check exists is worth keeping here, because it is the second time this
project has found the same defect. Check 21 searches the cold sandbox for the
string form of every graded value and requires an integer to carry four digits
and a float four decimals, since a graded count of 2 appears in any prose by
coincidence. Two shapes live in that gap:

- a small scalar stated in a sentence -- `15` as a locus total, `456` as a corpus
  count, `79` as a position, `0.62` as a GC floor;
- a declared term stated against its own label, where the label is a declared
  parameter and the term is published vocabulary, so neither half is distinctive
  and the pair is the answer.

An audit found twelve answers across nine campaigns in those two shapes, ten of
them disclosures and the rest copies of the campaign's own parameters. Every one
is fixed; this is what keeps them fixed.
"""

from __future__ import annotations

import json
import pathlib
import tempfile

import pytest

from npbench_c.readiness import cold_sandbox
from npbench_c.readiness.gate import _resolve_gold
from npbench_c.sweep.runner import build_sandbox

ROOT = pathlib.Path(__file__).resolve().parents[2]
CAMPAIGNS = ROOT / "campaigns"


def _campaigns() -> list[pathlib.Path]:
    return [c for c in sorted(CAMPAIGNS.iterdir())
            if (c / "gold" / "gold.json").is_file()
            and (c / "grading.json").is_file()]


@pytest.fixture(scope="module")
def sandboxes(tmp_path_factory) -> dict[str, pathlib.Path]:
    """One cold sandbox per campaign, built once for the whole module."""
    root = tmp_path_factory.mktemp("cold")
    return {c.name: build_sandbox(c, root / c.name, None) for c in _campaigns()}


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_no_graded_answer_is_stated_in_the_cold_sandbox(campaign, sandboxes):
    hits = cold_sandbox.findings(campaign, sandboxes[campaign.name], _resolve_gold)
    assert not hits, (
        f"{campaign.name}: the cold sandbox states {len(hits)} graded answer(s) "
        "that check 21's digit threshold does not reach:\n"
        + "\n".join(f"  {hit}" for hit in sorted(hits))
        + "\n\nRemove the value from task.yaml/reference (the measurement belongs "
          "in docs/, which no sandbox sees), or -- only if the match is a genuine "
          "collision -- declare it in audit/cold_sandbox_exemptions.json with the "
          "collision named."
    )


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_every_declared_collision_still_collides(campaign, sandboxes):
    """A waiver list rots silently, and this one has rotted twice already.

    The first time, removing a leaked sentence from `bgcdetect` also removed the
    collision a neighbouring waiver rested on. The second, an exemption was
    written for a component that the bounded scan never flags at all. Both were
    caught here rather than by reading the list.
    """
    stale = cold_sandbox.stale_exemptions(
        campaign, sandboxes[campaign.name], _resolve_gold)
    assert not stale, (
        f"{campaign.name}: {len(stale)} declared collision(s) in "
        "audit/cold_sandbox_exemptions.json no longer describe anything:\n"
        + "\n".join(f"  {s}" for s in stale))


def test_the_scan_reaches_every_built_campaign():
    """A scan that silently covered nothing would pass forever."""
    assert len(_campaigns()) >= 22


def test_every_exemption_carries_a_reason():
    declared = json.loads(cold_sandbox.EXEMPTIONS.read_text())["exemptions"]
    for name, entries in sorted(declared.items()):
        assert (CAMPAIGNS / name).is_dir(), f"{name} is not a campaign"
        for entry in entries:
            assert entry.get("component"), name
            assert len(entry.get("reason", "")) > 40, (name, entry.get("component"))


def test_the_exemptions_are_outside_every_sandbox(sandboxes):
    """The file names which declared value a component's answer equals, which
    states the answer. Shipping it would be the disclosure the check prevents."""
    assert cold_sandbox.EXEMPTIONS.parent.name == "audit"
    assert CAMPAIGNS not in cold_sandbox.EXEMPTIONS.parents
    body = cold_sandbox.EXEMPTIONS.read_text()
    for name, sandbox in sandboxes.items():
        for path in sandbox.rglob("*"):
            if path.is_file() and path.name == cold_sandbox.EXEMPTIONS.name:
                raise AssertionError(f"{name} ships the exemptions file")
    # And the reasons themselves must not have been copied into a sandbox.
    for name, sandbox in sandboxes.items():
        prose = cold_sandbox.sandbox_prose(sandbox)
        for entry in json.loads(body)["exemptions"].get(name, []):
            assert entry["reason"] not in prose, (name, entry["component"])


def test_the_scan_finds_a_planted_pair(tmp_path):
    """The check has to be able to fail. Rebuilds one campaign's sandbox with a
    note that states a verdict against its own declared label, which is the
    sentence that shipped in `frameshift-scoelicolor-01` before it was caught."""
    source = CAMPAIGNS / "frameshift-scoelicolor-01"
    if not (source / "gold" / "gold.json").is_file():
        pytest.skip("frameshift campaign gold not generated")
    gold = json.loads((source / "gold" / "gold.json").read_text())
    label, verdict = sorted(gold["verdict_by_length"].items(),
                            key=lambda kv: int(kv[0]))[0]

    campaign = tmp_path / source.name
    (campaign / "gold").mkdir(parents=True)
    (campaign / "gold" / "gold.json").write_text(json.dumps(gold))
    (campaign / "grading.json").write_text((source / "grading.json").read_text())
    sandbox = tmp_path / "sandbox"
    (sandbox / "reference").mkdir(parents=True)
    (sandbox / "inputs").mkdir()
    # The real rules file, because it declares the verdict vocabulary the pair
    # scan is bounded to.
    (sandbox / "reference" / "frameshift_rules.json").write_text(
        (source / "reference" / "frameshift_rules.json").read_text())
    (sandbox / "task.yaml").write_text(
        f"objective: the gene selected for {label} bases comes back {verdict}\n")

    hits = cold_sandbox.findings(campaign, sandbox, _resolve_gold)
    assert any(f"[{label}]" in hit for hit in hits), hits


def test_a_small_scalar_in_the_sandbox_is_found(tmp_path):
    """The other shape, planted: a graded two-digit integer in a sentence, which
    check 21 skips because it carries fewer than four digits."""
    campaign = tmp_path / "campaign"
    (campaign / "gold").mkdir(parents=True)
    (campaign / "gold" / "gold.json").write_text(json.dumps({"loci_total": 15}))
    (campaign / "grading.json").write_text(json.dumps({"rungs": [
        {"rung_id": "r1", "components": [
            {"name": "loci_total", "scorer": "exact",
             "args": {"pred": "pred:loci_total", "gold": "gold:loci_total"}}]}]}))
    sandbox = tmp_path / "sandbox"
    (sandbox / "inputs").mkdir(parents=True)
    (sandbox / "task.yaml").write_text(
        "license_notice: inputs/ carries a 15-row extract of the locus table\n")
    assert cold_sandbox.findings(campaign, sandbox, _resolve_gold) == [
        "loci_total = 15"]


def test_a_number_in_inputs_is_not_a_disclosure(tmp_path):
    """`inputs/` is excluded for the reason check 21 gives: a value there is the
    data the campaign asked the agent to read. Whether grading it is a FREE
    COMPONENT is a different question with its own answer per campaign."""
    campaign = tmp_path / "campaign"
    (campaign / "gold").mkdir(parents=True)
    (campaign / "gold" / "gold.json").write_text(json.dumps({"loci_total": 15}))
    (campaign / "grading.json").write_text(json.dumps({"rungs": [
        {"rung_id": "r1", "components": [
            {"name": "loci_total", "scorer": "exact",
             "args": {"pred": "pred:loci_total", "gold": "gold:loci_total"}}]}]}))
    sandbox = tmp_path / "sandbox"
    (sandbox / "inputs").mkdir(parents=True)
    (sandbox / "task.yaml").write_text("objective: count the loci\n")
    (sandbox / "inputs" / "loci.json").write_text(json.dumps({"rows": 15}))
    assert cold_sandbox.findings(campaign, sandbox, _resolve_gold) == []
