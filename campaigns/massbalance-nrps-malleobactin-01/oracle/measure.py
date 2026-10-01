"""Measure a submission. Deterministic, campaign-specific, agent-independent.

R1/R2 fields are recomputed from the agent's own answers where that is possible
(a formula string is re-parsed and re-rendered, so an unnormalised but correct
formula is not marked wrong on presentation). R3/R4 fields are the agent's
claims, normalised but not verified here: the grader compares them to gold.
"""

from __future__ import annotations

import json
import pathlib
import sys

import formula as F

CAMPAIGN = pathlib.Path(__file__).resolve().parent.parent
_VOCAB = json.loads((CAMPAIGN / "reference" / "assembly_rules.json").read_text())["verdict_vocabulary"]


def _canonical_formula(value: object) -> str | None:
    """Re-render a formula string in Hill order, or None if it does not parse.

    Normalising here means an agent writing H32C18N6O11 is credited: the output
    contract asks for Hill order, but grading presentation rather than chemistry
    would be exactly the brittle matching this benchmark exists to avoid.
    """
    if not isinstance(value, str):
        return None
    try:
        return F.render(F.parse(value))
    except F.FormulaError:
        return None


def _int_or_none(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float_or_none(value: object) -> float | None:
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def measure(submission_dir: pathlib.Path) -> dict:
    path = submission_dir / "route.json"
    out: dict = {"present": {"route": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["route"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["route"] = 0
        return out

    modules = sub.get("modules")
    modules = modules if isinstance(modules, list) else []
    out["n_modules"] = len(modules)
    out["monomers"] = [m.get("monomer") if isinstance(m, dict) else None for m in modules]
    out["monomer_formulas"] = {}
    parsed_all = bool(modules)
    for m in modules:
        if not isinstance(m, dict):
            parsed_all = False
            continue
        canon = _canonical_formula(m.get("formula"))
        if canon is None:
            parsed_all = False
        name = m.get("monomer")
        if isinstance(name, str):
            out["monomer_formulas"][name] = canon
    out["all_formulas_parse"] = int(parsed_all)

    out["naive_assembly_formula"] = _canonical_formula(sub.get("naive_assembly_formula"))
    out["peptide_bonds"] = _int_or_none(sub.get("peptide_bonds"))
    out["naive_assembly_monoisotopic_mass"] = _float_or_none(
        sub.get("naive_assembly_monoisotopic_mass"))

    # Per-compound answers. The entry annotates a congener series, so a single
    # residual would force an undeclared choice of which congener is "the"
    # product.
    per = sub.get("per_compound")
    per = per if isinstance(per, dict) else {}
    out["compounds"] = sorted(str(k) for k in per)
    out["per_compound"] = {}
    for name, block in sorted(per.items()):
        block = block if isinstance(block, dict) else {}
        residual = block.get("residual")
        residual = residual if isinstance(residual, dict) else {}
        verdict = block.get("verdict")
        out["per_compound"][str(name)] = {
            "residual_full": {
                el: (_int_or_none(residual.get(el)) if el in residual else 0)
                for el in F.ELEMENT_ORDER
            },
            "verdict": verdict if verdict in _VOCAB else None,
        }

    dels = sub.get("single_module_deletions")
    dels = dels if isinstance(dels, dict) else {}
    out["deletions"] = {}
    for key in (str(i + 1) for i in range(out["n_modules"] or 0)):
        block = dels.get(key) or {}
        block = block if isinstance(block, dict) else {}
        out["deletions"][key] = {
            "formula": _canonical_formula(block.get("formula")),
            "mass_shift": _float_or_none(block.get("mass_shift")),
        }
    return out


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: measure.py <submission_dir> <out_json>", file=sys.stderr)
        return 2
    result = measure(pathlib.Path(sys.argv[1]))
    pathlib.Path(sys.argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
