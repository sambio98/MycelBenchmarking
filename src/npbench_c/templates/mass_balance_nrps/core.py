"""Shared oracle for the NRPS mass-balance template.

One oracle, many instantiations. A campaign of this template contains no code at
all -- only inputs/, reference/, task.yaml, grading.json and gold/ -- so adding
an instantiation is a data change plus one build command. Copying six oracle
files per campaign would be the template architecture failing, and drift between
copies would be invisible.

Everything here is pure stdlib and deterministic: no network, no clock, no RNG,
no chemistry toolkit at run time. Monomer formulas are resolved once at build
time (see build.py) and pinned into the campaign's reference/ directory, which
keeps a toolkit version bump out of the scoring path and the campaign solvable
at tool surface S0.
"""

from __future__ import annotations

import json
import math
import pathlib
import re
from dataclasses import dataclass
from typing import Mapping, Sequence

TEMPLATE_ID = "mass-balance-nrps"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


class FormulaError(ValueError):
    """A formula string that cannot be parsed. An agent's malformed formula is a
    wrong answer, not a grader error, so callers catch this and score zero."""


class AssemblyError(RuntimeError):
    """The entry does not support an assembly computation. A campaign bug if it
    fires on a shipped entry."""


@dataclass(frozen=True)
class MassBalance:
    """A campaign's pinned tables, loaded once.

    Parameterised by campaign directory rather than hardcoded relative paths, so
    the same oracle serves every instantiation.
    """

    element_order: tuple[str, ...]
    masses: Mapping[str, float]
    monomers: Mapping[str, Mapping]
    rules: Mapping
    mass_decimals: int

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "MassBalance":
        reference = pathlib.Path(campaign) / "reference"
        rules = json.loads((reference / "assembly_rules.json").read_text())
        masses = json.loads((reference / "atomic_masses.json").read_text())["monoisotopic"]
        monomers = json.loads((reference / "monomer_formulas.json").read_text())["monomers"]
        return cls(
            element_order=tuple(rules["element_order"]),
            masses=masses,
            monomers=monomers,
            rules=rules,
            mass_decimals=int(rules["mass_decimals"]),
        )

    # ---------------------------------------------------------------- formulas

    def parse(self, formula: str) -> dict[str, int]:
        """Parse a Hill-order formula into element counts.

        Strict: only element-count tokens, every element in the declared order,
        every count positive. Silently accepting "C6H12O6!" would let a
        malformed submission score as if it had parsed.
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
            if el not in self.element_order:
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

    def render(self, counts: Mapping[str, int]) -> str:
        """Render in Hill order, so two correct answers cannot differ in form."""
        parts = []
        for el in self.element_order:
            n = counts.get(el, 0)
            if n == 0:
                continue
            if n < 0:
                raise FormulaError(f"cannot render a negative count for {el!r}: {n}")
            parts.append(el if n == 1 else f"{el}{n}")
        if not parts:
            raise FormulaError("cannot render an empty formula")
        return "".join(parts)

    def add(self, *operands: Mapping[str, int]) -> dict[str, int]:
        out: dict[str, int] = {}
        for op in operands:
            for el, n in op.items():
                out[el] = out.get(el, 0) + n
        return {el: out[el] for el in self.element_order if out.get(el)}

    def subtract(self, a: Mapping[str, int], b: Mapping[str, int]) -> dict[str, int]:
        """Signed difference. Negative counts are kept: a residual is a signed
        vector, and dropping its sign would lose the point of R3, where -H2O
        means cyclisation and a positive residual means the product carries
        atoms the modules do not supply."""
        out: dict[str, int] = {}
        for el in self.element_order:
            n = a.get(el, 0) - b.get(el, 0)
            if n:
                out[el] = n
        return out

    def scale(self, counts: Mapping[str, int], k: int) -> dict[str, int]:
        return {el: counts[el] * k for el in self.element_order if counts.get(el)}

    def monoisotopic_mass(self, counts: Mapping[str, int]) -> float:
        """Mass from the pinned table, reduced in declared element order so the
        sum is bit-stable, rounded to the contract's precision."""
        total = math.fsum(self.masses[el] * counts[el]
                          for el in self.element_order if counts.get(el))
        return round(total, self.mass_decimals)

    @property
    def water(self) -> dict[str, int]:
        return self.parse(self.rules["water_formula"])

    # ---------------------------------------------------------------- assembly

    def module_monomers(self, entry: Mapping) -> list[str]:
        names = []
        for index, module in enumerate(entry["biosynthesis"]["modules"], start=1):
            subs = (module.get("a_domain") or {}).get("substrates") or []
            if len(subs) != 1:
                raise AssemblyError(
                    f"module {index} has {len(subs)} annotated substrates")
            names.append(subs[0]["name"])
        return names

    def monomer_formula(self, name: str) -> dict[str, int]:
        if name not in self.monomers:
            raise AssemblyError(f"monomer {name!r} is not in the pinned reference table")
        return self.parse(self.monomers[name]["formula"])

    def naive_assembly(self, monomers: Sequence[str]) -> dict:
        """Formula of the linear assembly the annotated modules would produce."""
        if len(monomers) < 2:
            raise AssemblyError("need at least two monomers for an assembly")
        total = self.add(*(self.monomer_formula(m) for m in monomers))
        bonds = len(monomers) - 1
        counts = self.subtract(total, self.scale(self.water, bonds))
        if any(v < 0 for v in counts.values()):
            raise AssemblyError(f"assembly has negative atom counts: {counts}")
        return {
            "formula": self.render(counts),
            "counts": counts,
            "peptide_bonds": bonds,
            "monoisotopic_mass": self.monoisotopic_mass(counts),
        }

    def classify(self, residual: Mapping[str, int]) -> str:
        """Classify a residual vector by the declared precedence.

        Precedence is published because a mixed-sign residual could otherwise be
        argued into either verdict, and an unstated tie-break is a grader defect
        rather than a hard task.
        """
        nonzero = {el: n for el, n in residual.items() if n}
        for verdict in self.rules["verdict_precedence"]:
            if verdict == "insufficient":
                continue
            if verdict == "balanced" and not nonzero:
                return verdict
            if verdict == "cyclisation_only" and nonzero == {"H": -2, "O": -1}:
                return verdict
            if verdict == "additional_chemistry_required" and any(n > 0 for n in nonzero.values()):
                return verdict
            if verdict == "deficit_unexplained" and nonzero and all(n < 0 for n in nonzero.values()):
                return verdict
        raise AssemblyError(f"no verdict rule matched residual {residual}")

    def analyse(self, entry: Mapping) -> dict:
        """Full R2/R3 analysis.

        Every formula-bearing compound the entry annotates is analysed, not just
        the first listed. 17 of the 43 admissible entries annotate more than one,
        so "the" product is not well defined and silently taking compounds[0]
        would bake an undeclared convention into gold and mark a correct answer
        about another congener wrong.
        """
        monomers = self.module_monomers(entry)
        assembly = self.naive_assembly(monomers)

        compounds = [c for c in entry.get("compounds", []) if c.get("formula")]
        if not compounds:
            raise AssemblyError("entry annotates no compound with a molecular formula")

        per_compound = {}
        for compound in sorted(compounds, key=lambda c: c["name"]):
            product = self.parse(compound["formula"])
            residual = self.subtract(product, assembly["counts"])
            per_compound[compound["name"]] = {
                "product_formula": self.render(product),
                "product_monoisotopic_mass": self.monoisotopic_mass(product),
                "residual": residual,
                "residual_full": {el: residual.get(el, 0) for el in self.element_order},
                "verdict": self.classify(residual),
            }

        return {
            "monomers": monomers,
            "monomer_formulas": {m: self.render(self.monomer_formula(m))
                                 for m in sorted(set(monomers))},
            "n_modules": len(monomers),
            "peptide_bonds": assembly["peptide_bonds"],
            "naive_assembly_formula": assembly["formula"],
            "naive_assembly_monoisotopic_mass": assembly["monoisotopic_mass"],
            "compounds": sorted(per_compound),
            "per_compound": per_compound,
        }

    def single_module_deletions(self, monomers: Sequence[str]) -> dict[str, dict]:
        """R4: for every single-module deletion, the resulting assembly.

        Keyed by 1-based module index as a string, so the answer is a map rather
        than an order-dependent list. Deletion is by index, not by monomer name,
        so a repeated monomer yields two distinct deletions.
        """
        base = self.naive_assembly(monomers)
        out = {}
        for i in range(len(monomers)):
            remaining = list(monomers[:i]) + list(monomers[i + 1:])
            if len(remaining) < 2:
                out[str(i + 1)] = {"feasible": False,
                                   "reason": "fewer than two monomers remain"}
                continue
            a = self.naive_assembly(remaining)
            out[str(i + 1)] = {
                "feasible": True,
                "deleted_monomer": monomers[i],
                "formula": a["formula"],
                "peptide_bonds": a["peptide_bonds"],
                "monoisotopic_mass": a["monoisotopic_mass"],
                "mass_shift": round(a["monoisotopic_mass"] - base["monoisotopic_mass"],
                                    self.mass_decimals),
            }
        return out
