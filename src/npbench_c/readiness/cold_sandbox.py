"""What the cold sandbox states that it should not, in the shapes check 21 misses.

Check 21, `cold_sandbox_withholds_answers`, searches the cold sandbox for the
string form of every graded value and deliberately looks only at DISTINCTIVE
scalars: four or more digits for an integer, four or more decimals for a float.
Its docstring says plainly that this leaves small graded values unchecked, and the
threshold is right for what it does -- a graded count of 2 matches any prose.

Two shapes live in that gap, and an audit found twelve answers in them across
nine campaigns:

**A small scalar.** `15` as a locus total, `456` as a corpus count, `79` as a
position, `0.62` as a GC floor. Each sat in a sentence that was true and useful,
and each handed over a graded component. Anchoring on the component rather than
on the digit count reaches them, at the cost of real collisions, which are
declared.

**A term stated against its own label.** Where a graded value maps a declared
label to a term out of a closed vocabulary -- a verdict per gene length, a policy
per threshold -- neither half of the pair is a distinctive scalar. The label is a
declared parameter and the term is published vocabulary; both belong in the
sandbox. The answer is the two of them in one sentence, which is how this
project's own frameshift campaign shipped twelve of R3's twenty-one answers in a
rules file before the leak was caught.

**Collisions are declared, in `audit/cold_sandbox_exemptions.json`.** Not in the
campaign: an exemption says which declared value a component's answer equals,
which states the answer. The file sits outside `campaigns/` so that no sandbox
can reach it. Every exemption must still collide -- a stale one is a check nobody
is getting.

**What this does not reach**, stated plainly rather than implied, because check
21's own docstring does the same and that honesty is what made this gap findable:

- A NUMBER inside a graded composite, keyed by a label. A census of
  `{"abolished": 10, ...}` is graded whole, so check 21 searches for the whole
  JSON, and a sentence saying "10" near "abolished" discloses one cell that
  neither check sees. The pair scan covers string terms only. Extending it to
  numeric cells was tried on paper and dropped for now: small integers near a
  published vocabulary word collide constantly, and no measured instance of the
  leak exists in this project to calibrate against. It is a real gap, not a
  solved one.
- A number spelled out in words. "ten of the twenty-one targets" is a
  disclosure and no regex finds it. Nothing mechanical will; this is what the
  audit notes' own rule is for -- say where to recompute a number, never what it
  comes to.
- A value an agent can DERIVE from the sandbox without doing the work. That is a
  free component rather than a disclosure, it has a different fix (stop grading
  it), and no scan over prose will find it.
"""

from __future__ import annotations

import json
import pathlib
import re

#: How close a label and its term have to be to read as one statement about that
#: label. Wide enough to cross a clause, narrow enough that two unrelated
#: sentences in one note do not collide.
WINDOW = 160

#: Below this length a vocabulary term is an ordinary word ("one", "none",
#: "mixed") and matches any prose.
MIN_TERM = 6

#: Check 21 covers integers at or above this, so this check covers below it.
DISTINCTIVE_INT = 1000

#: Check 21 covers floats with at least this many decimals.
DISTINCTIVE_DECIMALS = 4

EXEMPTIONS = (pathlib.Path(__file__).resolve().parents[3]
              / "audit" / "cold_sandbox_exemptions.json")

#: Suffixes worth reading as prose. A sandbox also carries the inputs, which are
#: excluded below for the reason check 21 gives.
TEXT_SUFFIXES = (".yaml", ".yml", ".json", ".md", ".txt")


def exemptions(campaign_name: str, path: pathlib.Path | None = None) -> dict[str, str]:
    """Declared collisions for one campaign, as {component: reason}."""
    source = path or EXEMPTIONS
    if not source.is_file():
        return {}
    declared = json.loads(source.read_text()).get("exemptions") or {}
    return {entry["component"]: entry["reason"]
            for entry in declared.get(campaign_name, [])}


def sandbox_prose(sandbox: pathlib.Path) -> str:
    """Everything an agent can read before computing, minus `inputs/`.

    `inputs/` is excluded for the reason check 21 gives: a number there is the
    data the campaign asked the agent to read, and a graded value readable out of
    it without doing the work is a free component, which is a different defect.
    """
    parts = []
    for path in sorted(sandbox.rglob("*")):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        if path.relative_to(sandbox).parts[0] == "inputs":
            continue
        parts.append(path.read_text(errors="replace"))
    return "\n".join(parts)


def vocabularies(sandbox: pathlib.Path) -> set[str]:
    """Every term the campaign declares as a closed vocabulary.

    This bounds the pair scan, and the bound is the point. A graded string that
    is NOT a declared term -- a molecular formula, a unit, an InChIKey -- is
    published because the campaign cannot be solved without it, and whether
    grading it back is a free component is a separate question with its own
    answer per campaign.
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

    for path in sorted(sandbox.rglob("*.json")):
        if path.relative_to(sandbox).parts[0] == "inputs":
            continue
        try:
            walk(json.loads(path.read_text(errors="replace")))
        except json.JSONDecodeError:
            continue
    return terms


def _graded(spec: dict, gold: dict, resolve):
    """Every graded component with its gold value, skipping what cannot resolve."""
    for rung in spec["rungs"]:
        for component in rung["components"]:
            args = component.get("args") or {}
            reference = args.get("gold")
            if not isinstance(reference, str) or not reference.startswith("gold:"):
                continue
            try:
                value = resolve(gold, reference.split(":", 1)[1])
            except (KeyError, IndexError, TypeError, ValueError):
                continue
            yield component["name"], value


def _scalar_forms(value) -> list[str]:
    """The string forms worth searching for, for a value check 21 skips."""
    if isinstance(value, bool) or value is None:
        return []
    if isinstance(value, int):
        # Single digits match anything; four or more are check 21's job.
        return [str(value)] if 10 <= abs(value) < DISTINCTIVE_INT else []
    if isinstance(value, float):
        text = repr(value)
        _, _, decimals = text.partition(".")
        if len(decimals) >= DISTINCTIVE_DECIMALS or value in (0.0, 1.0):
            return []
        return sorted({text, str(value), f"{value:g}"})
    if isinstance(value, str):
        return [value] if len(value) >= 3 else []
    return []


def findings(campaign: pathlib.Path, sandbox: pathlib.Path, resolve) -> list[str]:
    """Graded answers the cold sandbox states, as "<component> = <value>".

    `resolve` is the gate's gold-path resolver, passed in so there is one
    implementation of that and it lives with the grader.
    """
    spec = json.loads((campaign / "grading.json").read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    prose = sandbox_prose(sandbox)
    terms = vocabularies(sandbox)
    waived = exemptions(campaign.name)

    hits: list[str] = []
    for name, value in _graded(spec, gold, resolve):
        if name in waived:
            continue
        # --- shape one: a small scalar stated anywhere
        for form in _scalar_forms(value):
            if isinstance(value, str) and value in terms:
                # A vocabulary term has to be published, and no label is paired
                # with it here; which term applies is the answer.
                continue
            pattern = r"(?<![\w.])" + re.escape(form) + r"(?![\w.])"
            if re.search(pattern, prose):
                hits.append(f"{name} = {value!r}")
                break
        else:
            # --- shape two: a declared term stated against its own label
            if not isinstance(value, dict):
                continue
            for label, term in sorted(value.items(), key=lambda kv: str(kv[0])):
                if not isinstance(term, str) or term not in terms:
                    continue
                for match in re.finditer(re.escape(str(label)), prose):
                    window = prose[max(0, match.start() - WINDOW):
                                   match.start() + WINDOW]
                    if term in window:
                        hits.append(f"{name}[{label}] = {term!r}")
                        break
                else:
                    continue
                break
    return hits


def stale_exemptions(campaign: pathlib.Path, sandbox: pathlib.Path,
                     resolve) -> list[str]:
    """Declared collisions that no longer collide, or name no graded component.

    The reason this exists: the first waiver list written in this project went
    stale within a minute, because removing a leaked sentence also removed the
    collision a neighbouring waiver rested on.
    """
    spec = json.loads((campaign / "grading.json").read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    graded = {name for name, _ in _graded(spec, gold, resolve)}
    waived = exemptions(campaign.name)
    if not waived:
        return []

    # Re-run the scan with nothing waived, and see which exemptions earn their
    # place. `findings` skips waived components, so ask it again via a copy of
    # the campaign's exemption set emptied out.
    live = set()
    prose = sandbox_prose(sandbox)
    terms = vocabularies(sandbox)
    for name, value in _graded(spec, gold, resolve):
        if name not in waived:
            continue
        for form in _scalar_forms(value):
            if isinstance(value, str) and value in terms:
                continue
            if re.search(r"(?<![\w.])" + re.escape(form) + r"(?![\w.])", prose):
                live.add(name)
                break
        if name in live or not isinstance(value, dict):
            continue
        for label, term in sorted(value.items(), key=lambda kv: str(kv[0])):
            if not isinstance(term, str) or term not in terms:
                continue
            for match in re.finditer(re.escape(str(label)), prose):
                if term in prose[max(0, match.start() - WINDOW):
                                 match.start() + WINDOW]:
                    live.add(name)
                    break
            if name in live:
                break

    out = []
    for name in sorted(waived):
        if name not in graded:
            out.append(f"{name} (not a graded component)")
        elif name not in live:
            out.append(f"{name} (no longer collides; drop the exemption)")
    return out
