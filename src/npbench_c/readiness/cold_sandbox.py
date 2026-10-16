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

**A term whose label is in the gold PATH, spelled as prose spells it.** The scan
above looks inside the graded value, so it reaches a mapping and misses a bare
string: `hosts.ecoli.binding_constraint` grades one word of published
vocabulary, and the label it belongs to is in the path, not the value. That is
how two construct-design campaigns shipped both halves of their binding
constraint -- R3's answer and R4's -- in a note explaining what the two
instantiations share: "gc_window for E. coli, direct_repeat for Streptomyces".
A literal search finds neither pair, because the path writes the host the way a
key allows and the prose writes it the way a reader expects. So both sides are
normalised to lower-case alphanumerics before the windows are compared, which
also collapses `gc_window` onto "GC window".

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

#: Below this length a path segment is not a label worth pairing on: "id",
#: "r3", "aa". Four characters is the shortest host, gene or rung key this
#: project uses as a gold-path segment.
MIN_LABEL = 4

#: Check 21 covers integers at or above this, so this check covers below it.
DISTINCTIVE_INT = 1000

#: Check 21 covers floats with at least this many decimals.
DISTINCTIVE_DECIMALS = 4

EXEMPTIONS = (pathlib.Path(__file__).resolve().parents[3]
              / "audit" / "cold_sandbox_exemptions.json")

#: A JSON string shorter than this is a key, a label or an identifier rather
#: than a sentence. Shape three runs over sentences only: a nested structure
#: necessarily states a field name inside its container's key, so scanning the
#: structure pairs every label with every field under it.
MIN_PROSE = 40

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


def sandbox_sentences(sandbox: pathlib.Path) -> str:
    """The sandbox's prose with its structure dropped.

    Shape three pairs a gold-path label with a vocabulary term, and both halves
    of such a pair sit in any nested document by construction: a host block keyed
    `ecoli_bl21` that declares a field named `gc_window` states the pair without
    telling anyone anything. What matters is a SENTENCE that does. So a JSON
    string is read only when it is long enough to be one, and a YAML or Markdown
    line is read with its leading key dropped and its value kept.
    """
    parts: list[str] = []

    def walk(node) -> None:
        if isinstance(node, dict):
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)
        elif isinstance(node, str) and len(node) >= MIN_PROSE:
            parts.append(node)

    for path in sorted(sandbox.rglob("*")):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        if path.relative_to(sandbox).parts[0] == "inputs":
            continue
        text = path.read_text(errors="replace")
        if path.suffix == ".json":
            try:
                walk(json.loads(text))
            except json.JSONDecodeError:
                parts.append(text)
            continue
        # The key names the field; the value is what somebody wrote.
        parts.append(re.sub(r"(?m)^[ \t-]*[A-Za-z_][\w.-]*:[ \t]*", "", text))
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
    """Every graded component, its gold path, and its value.

    The path is yielded because it carries the label half of a pair: a component
    reached at `hosts.ecoli.binding_constraint` answers a question about `ecoli`,
    and nothing in the value itself says so.
    """
    for rung in spec["rungs"]:
        for component in rung["components"]:
            args = component.get("args") or {}
            reference = args.get("gold")
            if not isinstance(reference, str) or not reference.startswith("gold:"):
                continue
            path = reference.split(":", 1)[1]
            try:
                value = resolve(gold, path)
            except (KeyError, IndexError, TypeError, ValueError):
                continue
            yield component["name"], path, value


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


def _normalise(text: str) -> str:
    """Lower-case alphanumerics only, so a label matches however it is spelled.

    The leak that forced this: a graded binding constraint stated as "gc_window
    for E. coli", against a gold path of `hosts.ecoli.binding_constraint`. The
    label is the same word in both and no literal search finds it, because prose
    writes an organism the way a reader expects and a gold path writes it the way
    a key allows. Normalising both sides collapses "E. coli", "E.coli" and
    "ecoli" onto one string -- and collapses the term too, so `gc_window` is
    found in "GC window" as well.
    """
    return re.sub(r"[^a-z0-9]", "", text.lower())


def _path_labels(path: str) -> list[str]:
    """The one segment of the gold path that labels what the component answers.

    `hosts.ecoli.binding_constraint` labels `ecoli`. The last segment is the
    question rather than the thing asked about, and the segments outside the
    innermost one are containers -- `hosts`, `per_query` -- which name a shape:
    taking those as labels pairs a container with every field declared under it,
    which fires on a nested document for nothing. A list index labels nothing at
    all, and neither does a segment too short to be a name.
    """
    parts = path.split(".")[:-1]
    if not parts or parts[-1].isdigit() or len(parts[-1]) < MIN_LABEL:
        return []
    return [parts[-1]]


def _pair_in(normalised_prose: str, label: str, term: str,
             min_label: int = 1) -> bool:
    """Is `term` stated within a window of `label`, both normalised?

    `min_label` is 1 for a label the graded value itself carries -- a gene length
    of `170` is a label, and the value is what says so -- and MIN_LABEL for one
    taken from the gold path, where a short segment is a container abbreviation
    rather than a name.
    """
    needle, target = _normalise(label), _normalise(term)
    if len(needle) < min_label or len(target) < MIN_TERM:
        return False
    for match in re.finditer(re.escape(needle), normalised_prose):
        window = normalised_prose[max(0, match.start() - WINDOW):
                                  match.start() + WINDOW]
        if target in window:
            return True
    return False


def _hit(name: str, path: str, value, prose: str, normalised_prose: str,
         normalised_sentences: str, terms: set[str]) -> str | None:
    """What the cold sandbox states about this component, if anything.

    Three shapes, in the order they were found. The first two are check 22 as
    first written; the third is the one that shipped a leak past both.
    """
    # --- shape one: a small scalar stated anywhere. A published vocabulary term
    # is not a disclosure on its own -- which term applies is the answer -- so a
    # bare term falls through to shape three, where its label decides.
    if not (isinstance(value, str) and value in terms):
        for form in _scalar_forms(value):
            if re.search(r"(?<![\w.])" + re.escape(form) + r"(?![\w.])", prose):
                return f"{name} = {value!r}"

    # --- shape two: a mapping of declared label to declared term, both halves
    # inside the graded value.
    if isinstance(value, dict):
        for label, term in sorted(value.items(), key=lambda kv: str(kv[0])):
            if isinstance(term, str) and term in terms and _pair_in(
                    normalised_prose, str(label), term):
                return f"{name}[{label}] = {term!r}"

    # --- shape three: a declared term whose label is in the gold PATH. The
    # value is one word of published vocabulary and the sandbox may say it; what
    # it may not say is which label that word belongs to.
    if isinstance(value, str) and value in terms:
        for label in _path_labels(path):
            if _pair_in(normalised_sentences, label, value, MIN_LABEL):
                return f"{name}@{label} = {value!r}"
    return None


def findings(campaign: pathlib.Path, sandbox: pathlib.Path, resolve) -> list[str]:
    """Graded answers the cold sandbox states, as "<component> = <value>"."""
    spec = json.loads((campaign / "grading.json").read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    prose = sandbox_prose(sandbox)
    normalised = _normalise(prose)
    sentences = _normalise(sandbox_sentences(sandbox))
    terms = vocabularies(sandbox)
    waived = exemptions(campaign.name)

    hits = []
    for name, path, value in _graded(spec, gold, resolve):
        if name in waived:
            continue
        found = _hit(name, path, value, prose, normalised, sentences, terms)
        if found is not None:
            hits.append(found)
    return hits


def stale_exemptions(campaign: pathlib.Path, sandbox: pathlib.Path,
                     resolve) -> list[str]:
    """Declared collisions that no longer collide, or name no graded component.

    The reason this exists: the first waiver list written in this project went
    stale within a minute, because removing a leaked sentence also removed the
    collision a neighbouring waiver rested on. It has fired twice since, both
    times on a number that moved under a campaign rebuild rather than on a
    sentence anybody edited.

    The scan is the same one `findings` runs, with the waiver ignored: an
    exemption earns its place only if the hit it waives is still there.
    """
    spec = json.loads((campaign / "grading.json").read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    waived = exemptions(campaign.name)
    if not waived:
        return []

    prose = sandbox_prose(sandbox)
    normalised = _normalise(prose)
    sentences = _normalise(sandbox_sentences(sandbox))
    terms = vocabularies(sandbox)
    graded, live = set(), set()
    for name, path, value in _graded(spec, gold, resolve):
        graded.add(name)
        if name in waived and _hit(name, path, value, prose, normalised,
                                   sentences, terms) is not None:
            live.add(name)

    out = []
    for name in sorted(waived):
        if name not in graded:
            out.append(f"{name} (not a graded component)")
        elif name not in live:
            out.append(f"{name} (no longer collides; drop the exemption)")
    return out
