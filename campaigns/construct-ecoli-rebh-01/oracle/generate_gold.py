"""Generate gold artifacts and the tolerance record for construct-design.

Gold is produced by running the oracle inside the frozen container (tier G1).
Nothing here is lifted from a publication.
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import construct as C
from measure import read_fasta
from tables import CONSTRAINT_ORDER

HERE = pathlib.Path(__file__).resolve().parent
CAMPAIGN = HERE.parent

HOST_SPEC = {
    "ecoli": ("ecoli_bl21", "pet28a_ndei_xhoi"),
    "streptomyces": ("streptomyces_coelicolor", "pset152_ndei_xhoi"),
}


def analyse(protein: str, host: str, strategy: str) -> dict:
    opt = C.unconstrained_optimum(protein, host)
    uv = C.violation_counts(opt, host, strategy)
    per_policy = {
        pol: C.constraint_costs(protein, host, strategy, pol)
        for pol in C.OPTIMIZER_POLICIES
    }
    bindings = {a["binding_constraint"] for a in per_policy.values()}
    if len(bindings) != 1:
        raise C.OracleFailure(
            f"binding constraint unstable across optimiser policies: {bindings}. "
            "Gold is not defensible; redesign the rung (spec 8.5)."
        )
    cais = [a["cai_full"] for a in per_policy.values()]
    non_binding = sorted(
        c for c in CONSTRAINT_ORDER
        if all(a["costs"][c] <= 0.0 for a in per_policy.values())
    )
    ref = per_policy["first_improve"]
    return {
        "unconstrained_violations": {k: uv[k] for k in CONSTRAINT_ORDER},
        "binding_constraint": bindings.pop(),
        "non_binding_constraints": non_binding,
        "reference_sequence": ref["sequence"],
        "reference_cai": ref["cai_full"],
        "cai_min_across_policies": round(min(cais), 6),
        "cai_max_across_policies": round(max(cais), 6),
        "cai_spread_across_policies": round(max(cais) - min(cais), 6),
        "cost_spread_across_policies": {
            c: round(
                max(a["costs"][c] for a in per_policy.values())
                - min(a["costs"][c] for a in per_policy.values()), 6
            )
            for c in CONSTRAINT_ORDER
        },
    }


def write_warm_bundles(protein: str, gold: dict) -> None:
    """Materialise exactly what each rung's warm start may see.

    Warm inputs are rendered bundles, never whole gold files. Declaring
    "gold/gold.json#hosts.ecoli" and copying the file would hand an R4 warm run
    the Streptomyces answer it is being asked to compute -- a leak that nothing
    downstream would detect, because the run would simply look successful.
    """
    warm = CAMPAIGN / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    def fasta(seq: str, name: str) -> str:
        return f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"

    # R2 starts from R1's deliverable: a schema-valid ORF that does NOT yet
    # satisfy the constraints. The naive reverse translation is exactly that.
    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "orf_ecoli_start.fasta").write_text(
        fasta(C.unconstrained_optimum(protein, "ecoli_bl21"), "r1_schema_valid_start"))

    # R3 starts from the correct E. coli construct.
    (warm / "r3").mkdir(parents=True)
    shutil.copyfile(CAMPAIGN / "gold" / "orf_ecoli.fasta",
                    warm / "r3" / "orf_ecoli.fasta")

    # R4 adds the R3 answer for E. coli only -- never the Streptomyces block.
    (warm / "r4").mkdir(parents=True)
    shutil.copyfile(CAMPAIGN / "gold" / "orf_ecoli.fasta",
                    warm / "r4" / "orf_ecoli.fasta")
    e = gold["hosts"]["ecoli"]
    (warm / "r4" / "analysis_ecoli.json").write_text(json.dumps({
        "unconstrained_violations": e["unconstrained_violations"],
        "binding_constraint": e["binding_constraint"],
        "non_binding_constraints": e["non_binding_constraints"],
    }, sort_keys=True, indent=2) + "\n")


def main() -> int:
    protein = read_fasta(CAMPAIGN / "inputs" / "rebh.faa") + "*"
    gold: dict = {"protein_length_aa": len(protein) - 1, "hosts": {}}

    for key, (host, strategy) in sorted(HOST_SPEC.items()):
        a = analyse(protein, host, strategy)
        seq = a.pop("reference_sequence")
        counts = C.violation_counts(seq, host, strategy)
        if sum(counts.values()) != 0:
            raise C.OracleFailure(f"{key}: reference sequence violates constraints")
        if C.translate(seq) != protein:
            raise C.OracleFailure(f"{key}: reference sequence does not encode the protein")
        # CAI floor: below the worst legitimate optimiser policy by a margin,
        # so an agent is never failed for using a different valid optimiser.
        a["cai_floor"] = round(a["cai_min_across_policies"] - 0.025, 3)
        gold["hosts"][key] = a
        (CAMPAIGN / "gold" / f"orf_{key}.fasta").write_text(
            f">gold_{key}_rebh_orf\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
        )

    e, s = gold["hosts"]["ecoli"], gold["hosts"]["streptomyces"]
    gold["flipped_constraints"] = sorted(
        set(e["non_binding_constraints"]) ^ set(s["non_binding_constraints"])
    )
    gold["binding_constraint_flips"] = int(e["binding_constraint"] != s["binding_constraint"])

    (CAMPAIGN / "gold" / "gold.json").write_text(json.dumps(gold, sort_keys=True, indent=2) + "\n")
    write_warm_bundles(protein, gold)
    print(json.dumps(gold, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
