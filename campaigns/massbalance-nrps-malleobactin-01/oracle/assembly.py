"""NRPS assembly arithmetic and residual analysis. Pure stdlib, deterministic.

What this computes is exactly what the annotated modules support, and no more:

  naive assembly  sum of the monomers the modules load, minus one water per
                  peptide bond
  residual        product formula minus naive assembly, signed, element by
                  element
  verdict         a mechanical classification of the residual vector, by the
                  precedence declared in reference/assembly_rules.json

It deliberately does NOT try to explain a residual biologically. MIBiG annotates
the enzymatic module count, not the number of monomer incorporations: iterative
module reuse, fatty-acid starters, prenyl groups and metal centres all break the
correspondence, and only 1 of 43 admissible entries is exactly balanced.
Quantifying the gap is well-posed; attributing it would be judgement.
"""

from __future__ import annotations

import json
import pathlib
from typing import Mapping, Sequence

import formula as F

_REFERENCE = pathlib.Path(__file__).resolve().parent.parent / "reference"
_RULES = json.loads((_REFERENCE / "assembly_rules.json").read_text())
_MONOMERS = json.loads((_REFERENCE / "monomer_formulas.json").read_text())["monomers"]

VERDICTS: list[str] = list(_RULES["verdict_vocabulary"])
PRECEDENCE: list[str] = list(_RULES["verdict_precedence"])


class AssemblyError(RuntimeError):
    """The entry does not support an assembly computation. A campaign bug if it
    fires on the shipped entry."""


def module_monomers(entry: Mapping) -> list[str]:
    """Monomer name per module, in module order.

    Raises if any module lacks exactly one annotated substrate: that is the
    'insufficient' case, and a campaign whose shipped entry hits it should not
    have been admitted.
    """
    names = []
    for index, module in enumerate(entry["biosynthesis"]["modules"], start=1):
        subs = (module.get("a_domain") or {}).get("substrates") or []
        if len(subs) != 1:
            raise AssemblyError(f"module {index} has {len(subs)} annotated substrates")
        names.append(subs[0]["name"])
    return names


def monomer_formula(name: str) -> dict[str, int]:
    if name not in _MONOMERS:
        raise AssemblyError(f"monomer {name!r} is not in the pinned reference table")
    return F.parse(_MONOMERS[name]["formula"])


def naive_assembly(monomers: Sequence[str]) -> dict:
    """Formula of the linear assembly the annotated modules would produce."""
    if len(monomers) < 2:
        raise AssemblyError("need at least two monomers for an assembly")
    total = F.add(*(monomer_formula(m) for m in monomers))
    bonds = len(monomers) - 1
    counts = F.subtract(total, F.scale(F.WATER, bonds))
    if any(v < 0 for v in counts.values()):
        raise AssemblyError(f"assembly has negative atom counts: {counts}")
    return {
        "formula": F.render(counts),
        "counts": counts,
        "peptide_bonds": bonds,
        "monoisotopic_mass": F.monoisotopic_mass(counts),
    }


def classify(residual: Mapping[str, int]) -> str:
    """Classify a residual vector by the declared precedence.

    Precedence is published because a mixed-sign residual could otherwise be
    argued into either 'additional_chemistry_required' or 'deficit_unexplained',
    and an unstated tie-break is a grader defect rather than a hard task.
    """
    nonzero = {el: n for el, n in residual.items() if n}
    for verdict in PRECEDENCE:
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


def analyse(entry: Mapping) -> dict:
    """Full R2/R3 analysis.

    Every formula-bearing compound the entry annotates is analysed, not just the
    first one listed. Many entries annotate a congener series -- this one
    annotates four, including a dimer -- and silently taking ``compounds[0]``
    would bake an undeclared convention into gold and mark a correct answer about
    a different congener wrong. The assembly itself is compound-independent; only
    the residual and the verdict vary.
    """
    monomers = module_monomers(entry)
    assembly = naive_assembly(monomers)

    compounds = [c for c in entry.get("compounds", []) if c.get("formula")]
    if not compounds:
        raise AssemblyError("entry annotates no compound with a molecular formula")

    per_compound = {}
    for compound in sorted(compounds, key=lambda c: c["name"]):
        product = F.parse(compound["formula"])
        residual = F.subtract(product, assembly["counts"])
        per_compound[compound["name"]] = {
            "product_formula": F.render(product),
            "product_monoisotopic_mass": F.monoisotopic_mass(product),
            "residual": residual,
            "residual_full": {el: residual.get(el, 0) for el in F.ELEMENT_ORDER},
            "verdict": classify(residual),
        }

    return {
        "monomers": monomers,
        "monomer_formulas": {m: F.render(monomer_formula(m)) for m in sorted(set(monomers))},
        "n_modules": len(monomers),
        "peptide_bonds": assembly["peptide_bonds"],
        "naive_assembly_formula": assembly["formula"],
        "naive_assembly_monoisotopic_mass": assembly["monoisotopic_mass"],
        "compounds": sorted(per_compound),
        "per_compound": per_compound,
    }


def single_module_deletions(monomers: Sequence[str]) -> dict[str, dict]:
    """R4: for every single-module deletion, the resulting assembly.

    Keyed by 1-based module index as a string, so the answer is a map rather
    than an order-dependent list and cannot be scored differently for being
    written in another order. Deletion is by index, not by monomer name, so a
    repeated monomer is handled correctly.
    """
    base = naive_assembly(monomers)
    out = {}
    for i in range(len(monomers)):
        remaining = list(monomers[:i]) + list(monomers[i + 1:])
        if len(remaining) < 2:
            out[str(i + 1)] = {"feasible": False, "reason": "fewer than two monomers remain"}
            continue
        a = naive_assembly(remaining)
        out[str(i + 1)] = {
            "feasible": True,
            "deleted_monomer": monomers[i],
            "formula": a["formula"],
            "peptide_bonds": a["peptide_bonds"],
            "monoisotopic_mass": a["monoisotopic_mass"],
            "mass_shift": round(a["monoisotopic_mass"] - base["monoisotopic_mass"],
                                _RULES["mass_decimals"]),
        }
    return out
