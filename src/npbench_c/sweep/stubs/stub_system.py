"""Deterministic stub systems: fixtures for testing the sweep harness itself.

These are NOT agents and are never part of a published result. They exist so
the harness's statistics and gate logic can be validated before any compute is
spent on real systems -- a sweep whose arithmetic is wrong is worse than no
sweep, because its numbers look authoritative.

A stub reaches the campaign's oracle module through an explicit --oracle path so
the privilege is visible in the command line. Real systems never get it.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

LEVELS = ("naive", "constrained", "analyst", "complete", "parametric")


def fasta(seq: str, name: str) -> str:
    return f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--oracle", required=True)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sys.path.insert(0, args.oracle)
    import construct as C
    from measure import read_fasta

    sub = pathlib.Path(args.submission)
    sub.mkdir(parents=True, exist_ok=True)
    protein = read_fasta(pathlib.Path("inputs/rebh.faa")) + "*"

    if args.level == "parametric":
        # The no-tool leakage ablation: answers from plausible priors without
        # computing anything. A high-GC actinomycete gene in E. coli "obviously"
        # strains the GC window, and round numbers look like measurements.
        (sub / "orf_ecoli.fasta").write_text(
            fasta(C.unconstrained_optimum(protein, "ecoli_bl21"), "guess"))
        (sub / "analysis.json").write_text(json.dumps({
            "r3": {"unconstrained_violations": {k: 0 for k in C.CONSTRAINT_ORDER},
                   "binding_constraint": "gc_window",
                   "non_binding_constraints": ["gc_global"]},
            "r4": {"unconstrained_violations": {k: 0 for k in C.CONSTRAINT_ORDER},
                   "binding_constraint": "gc_window",
                   "non_binding_constraints": ["gc_global"],
                   "flipped_constraints": []},
        }, sort_keys=True, indent=2))
        return 0

    if args.level == "naive":
        (sub / "orf_ecoli.fasta").write_text(
            fasta(C.unconstrained_optimum(protein, "ecoli_bl21"), "naive"))
        (sub / "analysis.json").write_text(json.dumps({"r3": {}, "r4": {}}))
        return 0

    ecoli = C.optimize(protein, "ecoli_bl21", "pet28a_ndei_xhoi")
    (sub / "orf_ecoli.fasta").write_text(fasta(ecoli, "designed"))

    if args.level == "constrained":
        (sub / "analysis.json").write_text(json.dumps({"r3": {}, "r4": {}}))
        return 0

    e = C.constraint_costs(protein, "ecoli_bl21", "pet28a_ndei_xhoi")
    e_uv = C.violation_counts(C.unconstrained_optimum(protein, "ecoli_bl21"),
                              "ecoli_bl21", "pet28a_ndei_xhoi")
    r3 = {"unconstrained_violations": e_uv,
          "binding_constraint": e["binding_constraint"],
          "non_binding_constraints": sorted(
              c for c in C.CONSTRAINT_ORDER if e["costs"][c] <= 0.0)}

    if args.level == "analyst":
        (sub / "analysis.json").write_text(json.dumps({"r3": r3, "r4": {}}, sort_keys=True, indent=2))
        return 0

    strep = C.optimize(protein, "streptomyces_coelicolor", "pset152_ndei_xhoi")
    (sub / "orf_streptomyces.fasta").write_text(fasta(strep, "designed"))
    s = C.constraint_costs(protein, "streptomyces_coelicolor", "pset152_ndei_xhoi")
    s_uv = C.violation_counts(
        C.unconstrained_optimum(protein, "streptomyces_coelicolor"),
        "streptomyces_coelicolor", "pset152_ndei_xhoi")
    s_nb = sorted(c for c in C.CONSTRAINT_ORDER if s["costs"][c] <= 0.0)
    r4 = {"unconstrained_violations": s_uv,
          "binding_constraint": s["binding_constraint"],
          "non_binding_constraints": s_nb,
          "flipped_constraints": sorted(set(r3["non_binding_constraints"]) ^ set(s_nb))}
    (sub / "analysis.json").write_text(json.dumps({"r3": r3, "r4": r4}, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
