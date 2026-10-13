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
