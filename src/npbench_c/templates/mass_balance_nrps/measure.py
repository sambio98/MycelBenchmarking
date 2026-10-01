"""Measure a submission against a campaign of this template. Deterministic.

Campaign-specific logic lives here rather than in the grader, which stays a
generic declarative comparator, and the agent cannot self-report properties of
its own answer.

R1/R2 fields are re-derived from the agent's own strings where possible: a
formula is re-parsed and re-rendered, so an unnormalised but chemically correct
answer is credited. Grading presentation rather than chemistry would be exactly
the brittle matching this benchmark exists to detect in others.

Usage: python -m npbench_c.templates.mass_balance_nrps.measure <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.mass_balance_nrps.core import FormulaError, MassBalance


def _int_or_none(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float_or_none(value: object) -> float | None:
    return (float(value) if isinstance(value, (int, float))
            and not isinstance(value, bool) else None)


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    mb = MassBalance.load(campaign)
    vocab = list(mb.rules["verdict_vocabulary"])

    def canonical(value: object) -> str | None:
        if not isinstance(value, str):
            return None
        try:
            return mb.render(mb.parse(value))
        except FormulaError:
            return None

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
        canon = canonical(m.get("formula"))
        if canon is None:
            parsed_all = False
        name = m.get("monomer")
        if isinstance(name, str):
            out["monomer_formulas"][name] = canon
    out["all_formulas_parse"] = int(parsed_all)

    out["naive_assembly_formula"] = canonical(sub.get("naive_assembly_formula"))
    out["peptide_bonds"] = _int_or_none(sub.get("peptide_bonds"))
    out["naive_assembly_monoisotopic_mass"] = _float_or_none(
        sub.get("naive_assembly_monoisotopic_mass"))

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
                for el in mb.element_order
            },
            "verdict": verdict if verdict in vocab else None,
        }

    dels = sub.get("single_module_deletions")
    dels = dels if isinstance(dels, dict) else {}
    out["deletions"] = {}
    for key in (str(i + 1) for i in range(out["n_modules"] or 0)):
        block = dels.get(key) or {}
        block = block if isinstance(block, dict) else {}
        out["deletions"][key] = {
            "formula": canonical(block.get("formula")),
            "mass_shift": _float_or_none(block.get("mass_shift")),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>", file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
