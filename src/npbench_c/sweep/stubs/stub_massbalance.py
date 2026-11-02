"""Deterministic stub systems for the mass-balance template.

Harness fixtures, not agents. They exist so the sweep's statistics and gate
logic can be exercised on this campaign; a stub-based sweep can never promote a
readiness check to PASS.

Unlike the construct-design stubs, these read ONLY the sandbox -- inputs/ and
reference/ -- with no privileged oracle path. That is deliberate: the
construct-design campaign shipped briefly unsolvable because its pinned tables
sat outside the sandbox, and stubs with oracle access hid it completely. A stub
restricted to the sandbox cannot hide that class of defect.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re

LEVELS = ("naive", "assembled", "analyst", "complete", "parametric")

TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse(f: str, order: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for el, n in TOKEN.findall(f):
        if el in order:
            out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def render(counts: dict[str, int], order: list[str]) -> str:
    return "".join(
        el if counts[el] == 1 else f"{el}{counts[el]}"
        for el in order if counts.get(el)
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    entry = json.loads((sandbox / "inputs" / "entry.json").read_text())
    rules = json.loads((sandbox / "reference" / "assembly_rules.json").read_text())
    masses = json.loads((sandbox / "reference" / "atomic_masses.json").read_text())["monoisotopic"]
    mono = json.loads((sandbox / "reference" / "monomer_formulas.json").read_text())["monomers"]
    order = rules["element_order"]
    dec = rules["mass_decimals"]

    monomers = [m["a_domain"]["substrates"][0]["name"]
                for m in entry["biosynthesis"]["modules"]]
    compounds = sorted(c for c in
                       (x["name"] for x in entry["compounds"] if x.get("formula")))
    formulas = {m: mono[m]["formula"] for m in monomers}

    sub = pathlib.Path(args.submission)
    sub.mkdir(parents=True, exist_ok=True)

    route: dict = {
        "modules": [{"module": i + 1, "monomer": m, "formula": formulas[m]}
                    for i, m in enumerate(monomers)],
        "per_compound": {},
        "single_module_deletions": {},
    }

    if args.level == "parametric":
        # Answers from plausible priors, campaign-agnostically: sums the monomers
        # without condensing any water (the error a non-computing system makes),
        # and assumes a curated assembly line accounts for its own product.
        total: dict[str, int] = {}
        for n in monomers:
            for el, v in parse(formulas[n], order).items():
                total[el] = total.get(el, 0) + v
        route["peptide_bonds"] = len(monomers)          # off by one
        route["naive_assembly_formula"] = render(total, order)
        route["naive_assembly_monoisotopic_mass"] = round(
            math.fsum(masses[el] * total[el] for el in order if total.get(el)), dec)
        for c in compounds:
            route["per_compound"][c] = {"residual": {}, "verdict": "balanced"}
        (sub / "route.json").write_text(json.dumps(route, indent=2, sort_keys=True))
        return 0

    if args.level == "naive":
        route["peptide_bonds"] = len(monomers)          # off by one
        (sub / "route.json").write_text(json.dumps(route, indent=2, sort_keys=True))
        return 0

    def assembly(names: list[str]) -> dict[str, int]:
        total: dict[str, int] = {}
        for n in names:
            for el, v in parse(formulas[n], order).items():
                total[el] = total.get(el, 0) + v
        bonds = len(names) - 1
        total["H"] = total.get("H", 0) - 2 * bonds
        total["O"] = total.get("O", 0) - bonds
        return {el: total[el] for el in order if total.get(el)}

    def mass(counts: dict[str, int]) -> float:
        return round(math.fsum(masses[el] * counts[el] for el in order if counts.get(el)), dec)

    base = assembly(monomers)
    route["peptide_bonds"] = len(monomers) - 1
    route["naive_assembly_formula"] = render(base, order)
    route["naive_assembly_monoisotopic_mass"] = mass(base)

    if args.level == "assembled":
        (sub / "route.json").write_text(json.dumps(route, indent=2, sort_keys=True))
        return 0

    products = {x["name"]: x["formula"] for x in entry["compounds"] if x.get("formula")}
    for c in compounds:
        prod = parse(products[c], order)
        resid = {el: prod.get(el, 0) - base.get(el, 0) for el in order}
        resid = {el: v for el, v in resid.items() if v}
        if not resid:
            verdict = "balanced"
        elif resid == {"H": -2, "O": -1}:
            verdict = "cyclisation_only"
        elif any(v > 0 for v in resid.values()):
            verdict = "additional_chemistry_required"
        else:
            verdict = "deficit_unexplained"
        route["per_compound"][c] = {"residual": resid, "verdict": verdict}

    if args.level == "analyst":
        (sub / "route.json").write_text(json.dumps(route, indent=2, sort_keys=True))
        return 0

    for i in range(len(monomers)):
        remaining = monomers[:i] + monomers[i + 1:]
        if len(remaining) < 2:
            continue
        a = assembly(remaining)
        route["single_module_deletions"][str(i + 1)] = {
            "formula": render(a, order),
            "mass_shift": round(mass(a) - mass(base), dec),
        }
    (sub / "route.json").write_text(json.dumps(route, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
