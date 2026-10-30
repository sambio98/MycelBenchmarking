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
  and the pair is the answer;
- the same pair with the label in the gold PATH rather than in the graded value,
  and spelled as prose spells it: `hosts.ecoli.binding_constraint` against "for
  E. coli". Both sides are normalised before the windows are compared, and this
  shape reads the sandbox's sentences rather than its structure, because a block
  keyed `ecoli_bl21` declaring a field named `gc_window` is not a disclosure.

And the same scan runs ACROSS campaigns, as check 23: a sandbox stating a
neighbouring campaign's answer passes a per-campaign scan every time, which is
how the malleobactin notes came to state sevadicin's verdict and how a codon
table's provenance block came to state the gene-calling campaign's graded GC
fraction in four campaigns at once. That scan runs the distinctive scalars and
the two pair shapes but NOT the small scalars, because anchoring on the component
is what makes a two-digit integer checkable and across 22 neighbours the anchor
is gone -- measured at 77 hits over 506 ordered pairs, all coincidence.

An audit found twelve answers across nine campaigns in the first two shapes, ten
of them disclosures and the rest copies of the campaign's own parameters. The
third shape was added when the codon tables were re-derived, and found three
more: both construct-design campaigns were naming the binding constraint for
both hosts -- R3's answer and R4's -- and the malleobactin campaign was naming
its own four verdicts and the sevadicin instantiation's. Every one is fixed; this
is what keeps them fixed.
"""

from __future__ import annotations

import json
import pathlib
import shutil
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
    # And the reasons themselves must not have been copied into a sandbox,
    # from either list: a cross-campaign reason names another campaign's answer.
    declared = json.loads(body)
    for name, sandbox in sandboxes.items():
        prose = cold_sandbox.sandbox_prose(sandbox)
        for key in ("exemptions", "cross_campaign_exemptions"):
            for entry in (declared.get(key) or {}).get(name, []):
                assert entry["reason"] not in prose, (key, name, entry["component"])


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


def _one_component_campaign(root: pathlib.Path, gold: dict, path: str,
                            name: str = "binding_constraint") -> pathlib.Path:
    campaign = root / "campaign"
    (campaign / "gold").mkdir(parents=True)
    (campaign / "gold" / "gold.json").write_text(json.dumps(gold))
    (campaign / "grading.json").write_text(json.dumps({"rungs": [
        {"rung_id": "r3", "components": [
            {"name": name, "scorer": "exact",
             "args": {"pred": f"pred:{name}", "gold": f"gold:{path}"}}]}]}))
    return campaign


def _sandbox_with(root: pathlib.Path, task: str,
                  vocabulary: list[str] | None = None) -> pathlib.Path:
    sandbox = root / "sandbox"
    (sandbox / "reference").mkdir(parents=True)
    (sandbox / "inputs").mkdir(exist_ok=True)
    (sandbox / "task.yaml").write_text(task)
    (sandbox / "reference" / "limits.json").write_text(json.dumps(
        {"constraint_order": vocabulary or
         ["forbidden_sites", "gc_window", "gc_global", "homopolymer",
          "direct_repeat", "internal_rbs"]}))
    return sandbox


GOLD_PAIR = {"hosts": {"ecoli": {"binding_constraint": "gc_window"},
                       "streptomyces": {"binding_constraint": "direct_repeat"}}}


def test_a_term_keyed_by_its_gold_path_is_found(tmp_path):
    """The shape that shipped: the label is in the path, not in the value.

    This is the sentence both construct-design campaigns carried, down to the
    spelling -- the gold path says `ecoli` and the note says "E. coli", so no
    literal search for the pair finds it.
    """
    campaign = _one_component_campaign(
        tmp_path, GOLD_PAIR, "hosts.ecoli.binding_constraint")
    sandbox = _sandbox_with(tmp_path, (
        "split_constraint: >\n"
        "  3 of 6 graded sub-answers are identical: the binding constraint for\n"
        "  both hosts (gc_window for E. coli, direct_repeat for Streptomyces).\n"))
    assert cold_sandbox.findings(campaign, sandbox, _resolve_gold) == [
        "binding_constraint@ecoli = 'gc_window'"]


def test_the_same_sentence_without_the_label_is_not_a_hit(tmp_path):
    """A published vocabulary term on its own is not a disclosure.

    Which term applies is the answer, and the campaign has to publish the
    vocabulary for two correct answers to be written the same way.
    """
    campaign = _one_component_campaign(
        tmp_path, GOLD_PAIR, "hosts.ecoli.binding_constraint")
    sandbox = _sandbox_with(tmp_path, (
        "split_constraint: >\n"
        "  3 of 6 graded sub-answers are identical, including the binding\n"
        "  constraint for each host; gold/gold.json names them.\n"))
    assert cold_sandbox.findings(campaign, sandbox, _resolve_gold) == []


def test_a_nested_structure_is_not_a_pair(tmp_path):
    """Shape three reads sentences, not structure.

    `codon_tables.json` keys a host block `ecoli_bl21` and declares a field in it
    called `gc_window`. Scanning the structure pairs every label with every field
    under it, which fired on both construct campaigns for nothing. The real
    table is used, so the exclusion is tested against what ships.
    """
    campaign = _one_component_campaign(
        tmp_path, GOLD_PAIR, "hosts.ecoli.binding_constraint")
    sandbox = _sandbox_with(tmp_path, "objective: design the construct\n")
    (sandbox / "reference" / "codon_tables.json").write_text(
        (ROOT / "src" / "npbench_c" / "templates" / "construct_design"
         / "constants" / "codon_tables.json").read_text())
    assert cold_sandbox.findings(campaign, sandbox, _resolve_gold) == []


def test_the_window_bounds_the_pair(tmp_path):
    """Two unrelated sentences in one note are not a statement about the label."""
    campaign = _one_component_campaign(
        tmp_path, GOLD_PAIR, "hosts.ecoli.binding_constraint")
    filler = "  " + ("the measurement is recorded in gold rather than here. " * 12)
    sandbox = _sandbox_with(tmp_path, (
        "notes_for_audit: >\n"
        "  the E. coli host is the primary one.\n"
        + filler + "\n"
        "  gc_window is one of the six declared constraints.\n"))
    assert cold_sandbox.findings(campaign, sandbox, _resolve_gold) == []


def test_only_the_innermost_path_segment_is_a_label(tmp_path):
    """`hosts` names a shape, not a thing asked about.

    With every segment taken as a label, a container name pairs with any field
    declared under it and the check drowns in its own output.
    """
    assert cold_sandbox._path_labels("hosts.ecoli.binding_constraint") == ["ecoli"]
    assert cold_sandbox._path_labels("per_query.0.verdict") == []
    assert cold_sandbox._path_labels("verdict") == []
    assert cold_sandbox._path_labels("hosts.ec.binding_constraint") == []


def test_normalisation_collapses_the_spellings():
    assert cold_sandbox._normalise("E. coli") == "ecoli"
    assert cold_sandbox._normalise("E.coli") == "ecoli"
    assert cold_sandbox._normalise("GC window") == "gcwindow"
    assert cold_sandbox._normalise("malleobactin D") == "malleobactind"


def test_sentences_drop_structure_and_keep_prose(tmp_path):
    sandbox = tmp_path / "sandbox"
    (sandbox / "inputs").mkdir(parents=True)
    (sandbox / "task.yaml").write_text(
        "hosts:\n  ecoli:\n    note: the bottleneck here is the gc_window limit\n")
    (sandbox / "short.json").write_text(json.dumps(
        {"ecoli": "gc_window", "reason": "a sentence long enough to be prose, "
                                         "which is what this check reads"}))
    sentences = cold_sandbox.sandbox_sentences(sandbox)
    assert "the bottleneck here is the gc_window limit" in sentences
    assert "hosts:" not in sentences
    # The short mapping value is a label, not a sentence, so it is not prose.
    assert "gc_window\"" not in sentences
    assert "a sentence long enough to be prose" in sentences


# ------------------------------------------------------------ across campaigns


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_no_other_campaigns_answer_is_stated(campaign, sandboxes):
    """Check 23 over the whole built set, one reader at a time."""
    others = [c for c in _campaigns() if c.name != campaign.name]
    hits = cold_sandbox.cross_findings(
        campaign, sandboxes[campaign.name], others, _resolve_gold)
    assert not hits, (
        f"{campaign.name}: its cold sandbox states {len(hits)} graded answer(s) "
        "belonging to another campaign shipped alongside it:\n"
        + "\n".join(f"  {hit}" for hit in sorted(hits))
        + "\n\nRemove the value from this campaign's task.yaml/reference (the "
          "measurement belongs in docs/), or -- only if the match is a genuine "
          "collision -- declare it under cross_campaign_exemptions in "
          "audit/cold_sandbox_exemptions.json.")


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_every_declared_cross_collision_still_collides(campaign, sandboxes):
    others = [c for c in _campaigns() if c.name != campaign.name]
    stale = cold_sandbox.stale_cross_exemptions(
        campaign, sandboxes[campaign.name], others, _resolve_gold)
    assert not stale, f"{campaign.name}: {stale}"


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_declared_terms_equals_the_sandbox_vocabulary(campaign, sandboxes):
    """What makes check 23 affordable, asserted rather than assumed.

    The check needs 22 neighbours' vocabularies per run and reads them from
    `reference/` instead of building 22 sandboxes, which would copy every
    neighbour's inputs. That shortcut is only valid while the two agree.
    """
    assert (cold_sandbox.declared_terms(campaign)
            == cold_sandbox.vocabularies(sandboxes[campaign.name]))


def test_every_cross_exemption_names_a_pair_and_a_reason():
    declared = json.loads(cold_sandbox.EXEMPTIONS.read_text()).get(
        "cross_campaign_exemptions") or {}
    assert declared, "the cross-campaign waiver list is empty; the frameshift/" \
                     "genecall grouping should be in it"
    for reader, entries in sorted(declared.items()):
        assert (CAMPAIGNS / reader).is_dir(), f"{reader} is not a campaign"
        for entry in entries:
            assert (CAMPAIGNS / entry["states"]).is_dir(), entry
            assert entry["states"] != reader, entry
            assert entry.get("component"), entry
            assert len(entry.get("reason", "")) > 40, entry


def test_the_cross_scan_finds_a_planted_neighbours_answer(tmp_path):
    """The check has to be able to fail: one campaign's prose, another's gold."""
    reader = _one_component_campaign(tmp_path, {"n": 1}, "n", name="count")
    neighbour = tmp_path / "neighbour"
    (neighbour / "gold").mkdir(parents=True)
    (neighbour / "gold" / "gold.json").write_text(json.dumps({"loci_total": 4552}))
    (neighbour / "grading.json").write_text(json.dumps({"rungs": [
        {"rung_id": "r3", "components": [
            {"name": "loci_total", "scorer": "exact",
             "args": {"pred": "pred:loci_total", "gold": "gold:loci_total"}}]}]}))
    sandbox = _sandbox_with(tmp_path, (
        "notes_for_audit: >\n"
        "  the neighbouring corpus carries 4,552 loci, which is why this one\n"
        "  samples rather than enumerates.\n"))
    hits = cold_sandbox.cross_findings(reader, sandbox, [neighbour], _resolve_gold)
    assert hits == ["neighbour:loci_total = 4552"]


def test_the_cross_scan_leaves_small_scalars_alone(tmp_path):
    """The measured exclusion. A neighbour's graded 30 in this campaign's prose
    is a coincidence 77 times out of 77, so the shape is not run across
    campaigns -- stated here so removing the exclusion fails a test rather than
    merely looking tidier."""
    reader = _one_component_campaign(tmp_path, {"n": 1}, "n", name="count")
    neighbour = tmp_path / "neighbour"
    (neighbour / "gold").mkdir(parents=True)
    (neighbour / "gold" / "gold.json").write_text(json.dumps({"windows": 30}))
    (neighbour / "grading.json").write_text(json.dumps({"rungs": [
        {"rung_id": "r3", "components": [
            {"name": "windows", "scorer": "exact",
             "args": {"pred": "pred:windows", "gold": "gold:windows"}}]}]}))
    sandbox = _sandbox_with(tmp_path, "objective: read the 30 declared windows\n")
    assert cold_sandbox.cross_findings(
        reader, sandbox, [neighbour], _resolve_gold) == []


def test_a_number_inside_a_token_is_not_a_statement_of_it():
    """The false positive that tightened the boundary: the chemical-space
    campaign ships a compound named S2502 and the diff campaign grades 2502."""
    prose = '{"name":"S2502","accession":"BGC0002502"}'
    assert not cold_sandbox.states_scalar(prose, prose, "2502")
    assert cold_sandbox.states_scalar("source_entries: 2502", "source_entries: 2502",
                                      "2502")
    grouped = "the corpus holds 4,552 loci"
    assert cold_sandbox.states_scalar(grouped, cold_sandbox.plain_digits(grouped),
                                      "4552")


def test_distinctive_forms_is_check_21s_class():
    assert cold_sandbox.distinctive_forms(2502) == ["2502", "2,502"]
    assert cold_sandbox.distinctive_forms(30) == []
    assert cold_sandbox.distinctive_forms(0.721187) == ["0.721187"]
    assert cold_sandbox.distinctive_forms(0.72) == []
    assert cold_sandbox.distinctive_forms("gc_window") == []
    assert cold_sandbox.distinctive_forms([]) == []
    assert cold_sandbox.distinctive_forms({"a": 1}) == ['{"a": 1}']


def test_the_cross_scan_catches_the_leak_that_forced_it(tmp_path, sandboxes):
    """The historical case, restored verbatim.

    The malleobactin campaign's audit notes carried this sentence, which names
    the sevadicin instantiation's graded verdict while making a true point about
    the template. Checks 21 and 22 passed it: the verdict is one word of a
    vocabulary malleobactin itself declares, and the label is in sevadicin's
    gold path.
    """
    mall = CAMPAIGNS / "massbalance-nrps-malleobactin-01"
    sev = CAMPAIGNS / "massbalance-nrps-sevadicin-01"
    for campaign in (mall, sev):
        if not (campaign / "gold" / "gold.json").is_file():
            pytest.skip(f"{campaign.name} gold not generated")

    sandbox = tmp_path / "cold"
    shutil.copytree(sandboxes[mall.name], sandbox)
    task = sandbox / "task.yaml"
    task.write_text(task.read_text() + (
        "\n  - The verdict varies across instantiations of this template\n"
        "    (sevadicin BGC0000426 is exactly balanced, obafluorin BGC0001437\n"
        "    differs by one water).\n"))

    assert cold_sandbox.findings(mall, sandbox, _resolve_gold) == []
    assert cold_sandbox.cross_findings(mall, sandbox, [sev], _resolve_gold) == [
        "massbalance-nrps-sevadicin-01:verdict_sevadicin@sevadicin = 'balanced'"]

# ------------------------------------------------- check 26: one-word answers


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_no_single_term_answer_is_stated_in_the_cold_sandbox(campaign, sandboxes):
    hits = cold_sandbox.single_term_findings(
        campaign, sandboxes[campaign.name], _resolve_gold)
    assert not hits, (
        f"{campaign.name}: a component whose WHOLE answer is one declared term "
        f"has that term stated in the cold sandbox:\n"
        + "\n".join(f"  {hit}" for hit in sorted(hits))
        + "\n\nThree campaigns shipped this: the constraint control pair (one of "
          "them as a template parameter, `expect: feasible`, and both in their "
          "campaign NAMES) and a mass-balance instantiation whose audit notes "
          "stated its verdict. Fix the sandbox, or -- only for a genuine "
          "collision, where no sentence attaches the term to the case it "
          "answers -- declare it under single_term_exemptions."
    )


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_every_declared_single_term_collision_still_collides(campaign, sandboxes):
    stale = cold_sandbox.stale_single_term_exemptions(
        campaign, sandboxes[campaign.name], _resolve_gold)
    assert not stale, (
        f"{campaign.name}: {len(stale)} declared single-term collision(s) no "
        "longer describe anything:\n" + "\n".join(f"  {s}" for s in stale))


def test_a_vocabulary_declaration_is_not_a_disclosure():
    """The three shapes a campaign publishes its vocabulary in, none of which is
    the leak: a flow list, a block list, and a term keyed to its definition."""
    text = """
verdict_vocabulary: [balanced, cyclisation_only, additional_chemistry_required]
accepted:
  - confirmed
  - tentative
verdict_rules:
  substantive_change: the values differ after normalisation
"""
    left = cold_sandbox.declarations_removed(text)
    for term in ("balanced", "cyclisation_only", "additional_chemistry_required",
                 "confirmed", "tentative", "substantive_change"):
        assert term not in left, term


def test_a_term_stated_against_any_label_is_found(tmp_path):
    """The gap this check closes. Check 22's third shape pairs a term with the
    label in its gold PATH; `expect: feasible` states R3's whole answer against
    a label that appears in no path, and passed every scan for months."""
    campaign = tmp_path / "campaign"
    (campaign / "gold").mkdir(parents=True)
    (campaign / "gold" / "gold.json").write_text(json.dumps({"verdict": "feasible"}))
    (campaign / "grading.json").write_text(json.dumps({"rungs": [{
        "rung_id": "r3",
        "components": [{"name": "verdict", "scorer": "exact",
                        "args": {"pred": "pred:verdict", "gold": "gold:verdict"}}],
    }]}))
    sandbox = tmp_path / "sandbox"
    (sandbox / "reference").mkdir(parents=True)
    (sandbox / "reference" / "rules.json").write_text(
        json.dumps({"verdict_vocabulary": ["feasible", "infeasible"]}))

    # The vocabulary alone is not a hit.
    (sandbox / "task.yaml").write_text("objective: assess the constraint set\n")
    assert not cold_sandbox.single_term_findings(campaign, sandbox, _resolve_gold)

    # The answer against an unrelated label is.
    (sandbox / "task.yaml").write_text("template_params:\n  expect: feasible\n")
    hits = cold_sandbox.single_term_findings(campaign, sandbox, _resolve_gold)
    assert hits == ["verdict = 'feasible'"], hits

    # And so is the answer in prose.
    (sandbox / "task.yaml").write_text(
        "objective: >\n  the constraint set here is feasible, so supply a design\n")
    assert cold_sandbox.single_term_findings(campaign, sandbox, _resolve_gold)


def test_a_multi_term_answer_is_not_this_check_s_business(tmp_path):
    """Scope: this check is for a component whose ENTIRE answer is one word. A
    per-item table of the same terms is check 22's shape two, where the pairing
    with a label is what matters."""
    campaign = tmp_path / "campaign"
    (campaign / "gold").mkdir(parents=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps({"table": {"a": "feasible", "b": "infeasible"}}))
    (campaign / "grading.json").write_text(json.dumps({"rungs": [{
        "rung_id": "r3",
        "components": [{"name": "table", "scorer": "exact",
                        "args": {"pred": "pred:table", "gold": "gold:table"}}],
    }]}))
    sandbox = tmp_path / "sandbox"
    (sandbox / "reference").mkdir(parents=True)
    (sandbox / "reference" / "rules.json").write_text(
        json.dumps({"verdict_vocabulary": ["feasible", "infeasible"]}))
    (sandbox / "task.yaml").write_text("objective: both feasible and infeasible\n")
    assert not cold_sandbox.single_term_findings(campaign, sandbox, _resolve_gold)
