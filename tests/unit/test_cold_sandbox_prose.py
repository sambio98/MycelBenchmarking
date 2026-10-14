"""The cold sandbox must not name a graded ANSWER in prose, in any shape.

The readiness gate's twenty-first check, `cold_sandbox_withholds_answers`,
searches `task.yaml` and `reference/` for the string form of every graded value
and deliberately looks only at distinctive scalars: an integer needs four digits
and a float four decimals, because a graded count of 2 occurs in any prose by
coincidence. Its own docstring says so.

That leaves a shape it cannot see. Where a graded value is a MAP from a declared
label to a term out of a closed vocabulary -- a verdict per gene length, a policy
per threshold, a call per cell -- neither side of the pair is a distinctive
scalar. The label is a declared length and the term is published vocabulary;
both belong in the sandbox. What must not be in the sandbox is the two of them in
one sentence, which is the answer.

This caught a live leak in `frameshift-scoelicolor-01`. A note explaining why the
declared target ladder samples one length band densely went on to say what the
verdicts in that band are, which handed over twelve of R3's twenty-one answers in
a rules file the agent reads before doing anything. Same mechanism the project
named when it found the scalar version: a note explaining a design decision
drifts into quoting the measurement that motivated it. The explanation needs the
phenomenon, not the answer.

Kept as a test rather than a gate check for now, because the gate re-runs every
campaign's gold to certify it and the antiSMASH campaigns cost half an hour each.
It is the same guarantee either way, and open item 11 records the migration.
"""

from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
CAMPAIGNS = ROOT / "campaigns"

#: How close a label and a term have to be to read as one statement about that
#: label. Wide enough to cross a clause, narrow enough that two unrelated
#: sentences in one note do not collide.
WINDOW = 160

#: A term shorter than this is a word, not a vocabulary label, and collides with
#: ordinary prose ("one", "none", "mixed").
MIN_TERM = 6


def _vocabularies(campaign: pathlib.Path) -> set[str]:
    """Every term the campaign itself declares as a closed vocabulary.

    This is what bounds the check, and the bound is the point. A graded string
    that is NOT a declared term -- a molecular formula, a unit, an InChIKey, a
    ChEBI identifier -- is published in `reference/` because the campaign cannot
    be solved without it, and whether grading it back is a *free component* is a
    separate question with its own answer per campaign. What this check is for is
    narrower and has no legitimate case: a term out of the campaign's own closed
    enum, printed next to the label it is the answer for.
    """
    terms: set[str] = set()

    def walk(node) -> None:
        if isinstance(node, dict):
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            if node and all(isinstance(v, str) for v in node):
                terms.update(v for v in node if len(v) >= MIN_TERM)
            for value in node:
                walk(value)

    reference = campaign / "reference"
    if reference.is_dir():
        for path in sorted(reference.rglob("*.json")):
            try:
                walk(json.loads(path.read_text()))
            except json.JSONDecodeError:
                continue
    return terms


def _campaigns() -> list[pathlib.Path]:
    return [c for c in sorted(CAMPAIGNS.iterdir())
            if (c / "gold" / "gold.json").is_file() and (c / "grading.json").is_file()]


def _graded_keys(spec: dict) -> set[str]:
    """The top-level gold keys the grader reads, from the spec rather than a list
    maintained by hand: a component added to a rung is covered without an edit."""
    keys = set()
    for rung in spec["rungs"]:
        for component in rung["components"]:
            args = component.get("args")
            if not isinstance(args, dict):
                continue
            for value in args.values():
                if isinstance(value, str) and value.startswith("gold:"):
                    keys.add(value.split("gold:", 1)[1].split(".")[0])
    return keys


def _prose(campaign: pathlib.Path) -> str:
    """Everything the agent can read before it has computed anything.

    `inputs/` is excluded for the reason the gate gives: a value readable out of
    the inputs is a free component, which is a different defect.
    """
    parts = [(campaign / "task.yaml").read_text()]
    reference = campaign / "reference"
    if reference.is_dir():
        parts += [p.read_text() for p in sorted(reference.rglob("*")) if p.is_file()]
    return "\n".join(parts)


def _leaks(campaign: pathlib.Path) -> list[str]:
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    spec = json.loads((campaign / "grading.json").read_text())
    prose = _prose(campaign)
    vocabulary = _vocabularies(campaign)
    out = []
    for key in sorted(_graded_keys(spec)):
        block = gold.get(key)
        if not isinstance(block, dict):
            continue
        for label, value in sorted(block.items(), key=lambda kv: str(kv[0])):
            if not isinstance(value, str) or len(value) < MIN_TERM:
                continue
            if value not in vocabulary:
                continue
            for match in re.finditer(re.escape(str(label)), prose):
                window = prose[max(0, match.start() - WINDOW):
                               match.start() + WINDOW]
                if value in window:
                    out.append(f"{key}[{label}] = {value!r}")
                    break
    return out


#: The accepted pairs, each with its reason. Both are the tool-identity component,
#: graded as a single claim about what made the output -- which tool, which
#: version, which mode, which flags. The mode is a declared term and the rules
#: file has to publish it; the version cannot be read out of the sandbox, so the
#: component is not obtainable without running the tool.
ACCEPTED = {
    "frameshift-scoelicolor-01": {"environment[noise_mode] = 'single'"},
    "genecall-scoelicolor-01": {"environment[baseline_mode] = 'single'"},
}


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_no_graded_term_is_stated_against_its_own_label(campaign):
    found = set(_leaks(campaign)) - ACCEPTED.get(campaign.name, set())
    assert not found, (
        f"{campaign.name}: the cold sandbox states {len(found)} graded "
        "answer(s) in prose, which the gate's scalar check cannot see:\n"
        + "\n".join(f"  {hit}" for hit in sorted(found)))


def test_the_scan_reaches_every_built_campaign():
    """A scan that silently covered nothing would pass forever."""
    assert len(_campaigns()) >= 22


def test_the_scan_finds_a_planted_leak(tmp_path):
    """The other half of the contract: the check has to be able to fail.

    Rebuilds one campaign's sandbox with a note that states a verdict against its
    own declared label, which is exactly the sentence that was shipped.
    """
    source = CAMPAIGNS / "frameshift-scoelicolor-01"
    if not (source / "gold" / "gold.json").is_file():
        pytest.skip("frameshift campaign gold not generated")
    copy = tmp_path / source.name
    (copy / "gold").mkdir(parents=True)
    (copy / "reference").mkdir()
    # The real rules file, because it is what declares the verdict vocabulary the
    # check is bounded to.
    (copy / "reference" / "frameshift_rules.json").write_text(
        (source / "reference" / "frameshift_rules.json").read_text())
    gold = json.loads((source / "gold" / "gold.json").read_text())
    (copy / "gold" / "gold.json").write_text(json.dumps(gold))
    (copy / "grading.json").write_text((source / "grading.json").read_text())
    label, verdict = sorted(gold["verdict_by_length"].items(), key=lambda kv: int(kv[0]))[0]
    (copy / "task.yaml").write_text(
        f"objective: the gene selected for {label} bases comes back {verdict}\n")
    assert any(f"verdict_by_length[{label}]" in hit for hit in _leaks(copy))

# ===================================================================== scalars
#
# The second shape, and the one that cost the most. Check 21 searches for the
# string form of every graded value and requires an integer to carry four digits
# and a float four decimals, because a graded count of 2 occurs in any prose by
# coincidence. That threshold is what let `15`, `456`, `105`, `173`, `79`, `23`,
# `153`, `0.62`, `0.55` and `4.5` through -- ten graded answers in eight
# campaigns, each in a sentence that was true, useful and well-intentioned.
#
# The fix is to anchor on the component rather than on the digit count: a graded
# scalar is a leak if its exact value appears, on a word boundary, anywhere the
# agent can read before computing. That produces real coincidences -- a window
# size that equals a violation count, a histogram edge that equals a miss count
# -- so each one is listed below with the reason it is a coincidence, and the
# list is short enough to read.

#: Graded scalars whose value appears in the sandbox by coincidence, not as a
#: statement about that component. Each entry is "<campaign>:<component>" with
#: the collision named, because an exemption nobody can check is worse than no
#: check at all.
SCALAR_COINCIDENCES = {
    # The shipped chemistry table carries a heavy-atom count per record, and two
    # of those integers equal these corpus-level answers.
    "chemspace-mibig-4_0-01:connectivity_groups_with_multiple_compound_keys":
        "matches compound_heavy_atoms=111 in the shipped chemistry table",
    "chemspace-mibig-4_0-01:distinctions_lost_at_connectivity":
        "matches scaffold_heavy_atoms=119 in the shipped chemistry table",
    # The counting rules declare 50-nt GC windows and direct repeats as 11-mers;
    # RebH's two violation counts happen to equal both parameters.
    "construct-ecoli-rebh-01:uv_gc_window":
        "matches the declared 50-nt window length in the counting rules",
    "construct-ecoli-rebh-01:uv_direct_repeat":
        "matches the declared 11-mer repeat length in the counting rules",
    # The identity-band edges are declared as [70, 80], [80, 90] and so on.
    "transfer-ec1_14-uniprot-01:queries_without_a_hit":
        "matches a declared identity-band edge in the binning rule",
}

#: Graded scalars whose value is the campaign's own declared parameter, so the
#: sandbox must state it and the component would be a FREE COMPONENT rather than
#: a disclosure.
#:
#: The list is empty, and keeping it is the point: it held twenty entries across
#: nine campaigns, every one of them dropped from grading with the exclusion
#: recorded. An entry here is a component this project grades and should not, so
#: the empty dict is the claim that none remain. Adding one is a deliberate act
#: that has to survive `test_every_waiver_still_describes_a_real_collision`.
SCALAR_FREE_COMPONENTS: dict[str, str] = {}


def _dig(gold: dict, path: str):
    cursor = gold
    for part in path.split("."):
        if isinstance(cursor, dict) and part in cursor:
            cursor = cursor[part]
        else:
            return None
    return cursor


def _graded_scalars(campaign: pathlib.Path):
    """Every graded component whose gold value is a scalar, with that value."""
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    spec = json.loads((campaign / "grading.json").read_text())
    for rung in spec["rungs"]:
        for component in rung["components"]:
            args = component.get("args") or {}
            reference = args.get("gold")
            if not isinstance(reference, str) or not reference.startswith("gold:"):
                continue
            value = _dig(gold, reference.split("gold:", 1)[1])
            if isinstance(value, bool) or not isinstance(value, (int, float, str)):
                continue
            yield component["name"], value


def _forms(value) -> set[str]:
    if isinstance(value, float):
        return {repr(value), str(value), f"{value:g}"}
    return {str(value)}


def _scalar_leaks(campaign: pathlib.Path) -> list[str]:
    prose = _prose(campaign)
    vocabulary = _vocabularies(campaign)
    out = []
    for name, value in _graded_scalars(campaign):
        if isinstance(value, str):
            # A vocabulary term has to be published; which term applies is the
            # answer, and no label is paired with it here.
            if len(value) < 3 or value in vocabulary:
                continue
        elif isinstance(value, int) and abs(value) < 10:
            # A single digit matches anything. Check 21 skips these too, and for
            # the same reason.
            continue
        elif isinstance(value, float) and value in (0.0, 1.0):
            continue
        for form in _forms(value):
            if re.search(rf"(?<![\w.]){re.escape(form)}(?![\w.])", prose):
                out.append(f"{name} = {value!r}")
                break
    return out


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_no_graded_scalar_is_stated_in_the_cold_sandbox(campaign):
    waived = {key.split(":", 1)[1] for key in
              (*SCALAR_COINCIDENCES, *SCALAR_FREE_COMPONENTS)
              if key.startswith(campaign.name + ":")}
    found = [hit for hit in _scalar_leaks(campaign)
             if hit.split(" = ")[0] not in waived]
    assert not found, (
        f"{campaign.name}: the cold sandbox states {len(found)} graded value(s) "
        "that check 21's digit threshold does not reach:\n"
        + "\n".join(f"  {hit}" for hit in sorted(found))
        + "\n\nEither remove the value from task.yaml/reference (the "
          "measurement belongs in docs/), or, if the sandbox must state it, "
          "stop grading it and add it to SCALAR_FREE_COMPONENTS with its reason."
    )


def test_every_waiver_still_describes_a_real_collision():
    """A waiver list rots silently: an exempted component that no longer collides
    is a check nobody is getting, and one whose campaign was fixed is noise."""
    by_campaign = {c.name: c for c in _campaigns()}
    stale = []
    for key in (*SCALAR_COINCIDENCES, *SCALAR_FREE_COMPONENTS):
        name, component = key.split(":", 1)
        campaign = by_campaign.get(name)
        if campaign is None:
            continue
        if component not in {n for n, _ in _graded_scalars(campaign)}:
            stale.append(f"{key} (no such graded component)")
            continue
        if not any(hit.split(" = ")[0] == component
                   for hit in _scalar_leaks(campaign)):
            stale.append(f"{key} (no longer collides; drop the waiver)")
    assert not stale, "stale waivers:\n" + "\n".join(f"  {s}" for s in stale)


def test_no_free_component_is_left_in_grading():
    """What was open item 10, now closed and held closed. Twenty graded
    components across nine campaigns restated the campaign's own definition --
    an accession, an EC number, a snapshot release, a declared convention, an
    invocation parameter, a published cutoff -- and every one has been dropped
    from grading with its exclusion recorded. A new entry in the waiver list is
    a component this project grades and should not."""
    assert SCALAR_FREE_COMPONENTS == {}

