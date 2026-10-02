"""Deterministic stub systems for the constraint-conflict ladder.

Harness fixtures, not agents; a stub-based sweep can never promote a readiness
check to PASS. These read ONLY the sandbox -- inputs/ and reference/ -- with no
privileged oracle path, so they cannot hide an unsolvable campaign the way the
construct-design stubs did.
"""

from __future__ import annotations

import argparse
import json
import pathlib

LEVELS = ("naive", "bounds", "analyst", "complete", "optimist")

AA_CODONS = {
    "A": "GCT GCC GCA GCG", "R": "CGT CGC CGA CGG AGA AGG", "N": "AAT AAC",
    "D": "GAT GAC", "C": "TGT TGC", "Q": "CAA CAG", "E": "GAA GAG",
    "G": "GGT GGC GGA GGG", "H": "CAT CAC", "I": "ATT ATC ATA",
    "L": "TTA TTG CTT CTC CTA CTG", "K": "AAA AAG", "M": "ATG",
    "F": "TTT TTC", "P": "CCT CCC CCA CCG", "S": "TCT TCC TCA TCG AGT AGC",
    "T": "ACT ACC ACA ACG", "W": "TGG", "Y": "TAT TAC",
    "V": "GTT GTC GTA GTG", "*": "TAA TGA TAG",
}
AA_CODONS = {k: v.split() for k, v in AA_CODONS.items()}
GC = lambda c: sum(1 for x in c if x in "GC")  # noqa: E731


def bounds(protein: str, forbidden: set[str]) -> tuple[float, float, dict]:
    lo = hi = 0
    counts = {}
    for aa in sorted(set(protein)):
        counts[aa] = len([c for c in AA_CODONS[aa] if c not in forbidden])
    for aa in protein:
        allowed = [c for c in AA_CODONS[aa] if c not in forbidden]
        lo += min(GC(c) for c in allowed)
        hi += max(GC(c) for c in allowed)
    n = 3 * len(protein)
    return round(lo / n, 6), round(hi / n, 6), counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    spec = json.loads((sandbox / "reference" / "conflict_constraints.json").read_text())
    fasta = next(p for p in sorted((sandbox / "inputs").glob("*.faa")))
    protein = "".join(l.strip() for l in fasta.read_text().strip().splitlines()[1:]) + "*"

    forbidden = set(spec["forbidden_codons"])
    required = float(spec["gc_global_min"])
    vocab = list(spec["constraint_vocabulary"])
    _, _, counts = bounds(protein, forbidden)

    sub = pathlib.Path(args.submission)
    sub.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "orf_length_nt": 3 * len(protein),
        "constraints_assessed": vocab,
        "admissible_codon_counts": counts,
        "design_orf": None,
        "relaxations": {},
    }

    if args.level == "optimist":
        # The no-computation ablation: restates the constraints, assumes a
        # vendor-supplied set must be workable, and fabricates a design rather
        # than deriving the achievable range.
        report["verdict"] = "feasible"
        report["conflicting_constraints"] = []
        report["design_orf"] = "".join(
            sorted((c for c in AA_CODONS[aa] if c not in forbidden),
                   key=lambda c: (-GC(c), c))[0] for aa in protein)
        (sub / "feasibility_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if args.level == "naive":
        (sub / "feasibility_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    u_lo, u_hi, _ = bounds(protein, set())
    r_lo, r_hi, _ = bounds(protein, forbidden)
    report.update({"unrestricted_gc_min": u_lo, "unrestricted_gc_max": u_hi,
                   "restricted_gc_min": r_lo, "restricted_gc_max": r_hi})

    if args.level == "bounds":
        (sub / "feasibility_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    report["verdict"] = ("gc_requirement_unreachable" if required > r_hi else "feasible")
    if required > r_hi:
        report["verdict"] = "infeasible_gc_unreachable"
        report["conflicting_constraints"] = sorted(["forbidden_codons", "gc_global_min"])
    else:
        report["conflicting_constraints"] = []
    report["witness"] = {"required_gc_min": required,
                         "max_achievable_gc_under_forbidden_codons": r_hi}

    if args.level == "analyst":
        (sub / "feasibility_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    for c in vocab:
        if c == "forbidden_codons":
            gc_max, reachable = u_hi, required <= u_hi + 1e-06
        elif c == "gc_global_min":
            gc_max, reachable = r_hi, True
        else:
            gc_max, reachable = r_hi, required <= r_hi + 1e-06
        report["relaxations"][c] = {"max_achievable_gc": gc_max,
                                    "gc_requirement_reachable": bool(reachable)}
    (sub / "feasibility_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
