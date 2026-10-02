"""Deterministic stub systems for the construct-design ladder.

Harness fixtures, not agents; a stub-based sweep can never promote a readiness
check to PASS.

These load the shared constraint engine against the SANDBOX, not the campaign:
the engine is library code, like any tool an agent might use, and every table
and parameter comes from inputs/, reference/ and task.yaml as shipped. The
earlier version took an --oracle path into the campaign, and that privilege is
exactly what hid a campaign whose pinned tables were outside the sandbox.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import yaml

from npbench_c.templates.construct_design.core import ConstructEngine

LEVELS = ("naive", "constrained", "analyst", "complete", "parametric")


def fasta(seq: str, name: str) -> str:
    return f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    engine = ConstructEngine.load(sandbox)
    params = yaml.safe_load((sandbox / "task.yaml").read_text())["template_params"]
    primary, secondary = sorted(params["hosts"])

    path = sandbox / params["protein_input"]
    protein = "".join(l.strip().upper()
                      for l in path.read_text().strip().splitlines()[1:]) + "*"

    sub = pathlib.Path(args.submission)
    sub.mkdir(parents=True, exist_ok=True)

    def host_of(key):
        return params["hosts"][key]["host"], params["hosts"][key]["strategy"]

    p_host, p_strategy = host_of(primary)
    s_host, s_strategy = host_of(secondary)

    if args.level == "parametric":
        # The no-computation ablation: reverse-translates with the most-adapted
        # codon per residue and assumes a high-GC gene must strain the GC window.
        (sub / f"orf_{primary}.fasta").write_text(
            fasta(engine.unconstrained_optimum(protein, p_host), "guess"))
        (sub / "analysis.json").write_text(json.dumps({
            "r3": {"unconstrained_violations": {c: 0 for c in engine.constraint_order},
                   "binding_constraint": "gc_window",
                   "non_binding_constraints": ["gc_global"]},
            "r4": {"unconstrained_violations": {c: 0 for c in engine.constraint_order},
                   "binding_constraint": "gc_window",
                   "non_binding_constraints": ["gc_global"],
                   "flipped_constraints": []},
        }, indent=2, sort_keys=True))
        return 0

    if args.level == "naive":
        (sub / f"orf_{primary}.fasta").write_text(
            fasta(engine.unconstrained_optimum(protein, p_host), "naive"))
        (sub / "analysis.json").write_text(json.dumps({"r3": {}, "r4": {}}))
        return 0

    (sub / f"orf_{primary}.fasta").write_text(
        fasta(engine.optimize(protein, p_host, p_strategy), "designed"))

    if args.level == "constrained":
        (sub / "analysis.json").write_text(json.dumps({"r3": {}, "r4": {}}))
        return 0

    def analysis_for(host, strategy):
        costs = engine.constraint_costs(protein, host, strategy)
        uv = engine.violation_counts(engine.unconstrained_optimum(protein, host),
                                     host, strategy)
        return {
            "unconstrained_violations": uv,
            "binding_constraint": costs["binding_constraint"],
            "non_binding_constraints": sorted(
                c for c in engine.constraint_order if costs["costs"][c] <= 0.0),
        }

    r3 = analysis_for(p_host, p_strategy)
    if args.level == "analyst":
        (sub / "analysis.json").write_text(
            json.dumps({"r3": r3, "r4": {}}, indent=2, sort_keys=True))
        return 0

    (sub / f"orf_{secondary}.fasta").write_text(
        fasta(engine.optimize(protein, s_host, s_strategy), "designed"))
    r4 = analysis_for(s_host, s_strategy)
    r4["flipped_constraints"] = sorted(
        set(r3["non_binding_constraints"]) ^ set(r4["non_binding_constraints"]))
    (sub / "analysis.json").write_text(
        json.dumps({"r3": r3, "r4": r4}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
