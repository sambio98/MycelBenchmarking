"""Molecular formula arithmetic. Pure stdlib, deterministic.

No chemistry toolkit at runtime: monomer formulas are resolved once at campaign
construction and pinned into reference/monomer_formulas.json, so the oracle and
the agent both read the same table and neither needs RDKit. That keeps this
campaign solvable at tool surface S0 and keeps a toolkit version bump out of the
scoring path.
"""

from __future__ import annotations

import json
import pathlib
import re
from typing import Mapping

_REFERENCE = pathlib.Path(__file__).resolve().parent.parent / "reference"

_RULES = json.loads((_REFERENCE / "assembly_rules.json").read_text())
_MASSES: dict[str, float] = json.loads((_REFERENCE / "atomic_masses.json").read_text())["monoisotopic"]

# Hill order: carbon, then hydrogen, then every other element alphabetically.
# Declared so that two correct answers cannot differ in their string form.
ELEMENT_ORDER: list[str] = list(_RULES["element_order"])

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


class FormulaError(ValueError):
    """A formula string that cannot be parsed. An agent's malformed formula is a
    wrong answer, not a grader error, so callers catch this and score zero."""


def parse(formula: str) -> dict[str, int]:
    """Parse a Hill-order formula string into element counts.

    Strict: the formula must consist only of element-count tokens, every element
    must be in the declared order list, and counts must be positive. Anything
    else raises, because silently accepting "C6H12O6!" would let a malformed
    submission score as if it parsed.
    """
    if not isinstance(formula, str) or not formula.strip():
        raise FormulaError("empty formula")
    text = formula.strip()
    out: dict[str, int] = {}
    pos = 0
    for m in _TOKEN.finditer(text):
        if m.start() != pos:
            raise FormulaError(f"unparsed characters in {formula!r} at {pos}")
        pos = m.end()
        el, count = m.group(1), m.group(2)
        if el not in ELEMENT_ORDER:
            raise FormulaError(f"element {el!r} not in the declared element order")
        n = int(count) if count else 1
        if n <= 0:
            raise FormulaError(f"non-positive count for {el!r}")
        out[el] = out.get(el, 0) + n
    if pos != len(text):
        raise FormulaError(f"trailing characters in {formula!r}")
    if not out:
        raise FormulaError(f"no elements in {formula!r}")
    return out


def render(counts: Mapping[str, int]) -> str:
    """Render element counts in Hill order, omitting zeros and count-1 suffixes."""
    parts = []
    for el in ELEMENT_ORDER:
        n = counts.get(el, 0)
        if n == 0:
            continue
        if n < 0:
            raise FormulaError(f"cannot render a negative count for {el!r}: {n}")
        parts.append(el if n == 1 else f"{el}{n}")
    if not parts:
        raise FormulaError("cannot render an empty formula")
    return "".join(parts)


def add(*operands: Mapping[str, int]) -> dict[str, int]:
    out: dict[str, int] = {}
    for op in operands:
        for el, n in op.items():
            out[el] = out.get(el, 0) + n
    return {el: out[el] for el in ELEMENT_ORDER if out.get(el)}


def subtract(a: Mapping[str, int], b: Mapping[str, int]) -> dict[str, int]:
    """Signed difference. Negative counts are kept: a residual is a signed
    vector, and dropping its sign would lose the whole point of R3."""
    out: dict[str, int] = {}
    for el in ELEMENT_ORDER:
        n = a.get(el, 0) - b.get(el, 0)
        if n:
            out[el] = n
    return out


def scale(counts: Mapping[str, int], k: int) -> dict[str, int]:
    return {el: counts[el] * k for el in ELEMENT_ORDER if counts.get(el)}


def monoisotopic_mass(counts: Mapping[str, int]) -> float:
    """Monoisotopic mass from the pinned atomic-mass table.

    Reduction is in declared element order so the sum is bit-stable, and the
    result is rounded to the precision the output contract states.
    """
    import math

    total = math.fsum(_MASSES[el] * counts[el] for el in ELEMENT_ORDER if counts.get(el))
    return round(total, _RULES["mass_decimals"])


WATER: dict[str, int] = parse(_RULES["water_formula"])
