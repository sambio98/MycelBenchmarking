"""Build gold, warm bundles and the grading spec for the conflict ladder.

Usage: python -m npbench_c.templates.construct_design.conflict.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

from npbench_c.templates.construct_design.core import ConstructEngine, OracleFailure
from npbench_c.templates.construct_design.params import (
    CONFLICT_CONSTRAINTS,
    VERDICTS,
    ConflictSpec,
    copy_constants,
    load_task,
    read_protein,
)

GC_TOLERANCE = 1e-6

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.01, "r3": 0.01, "r4": 0.0005}

# Constraints that can change the achievable whole-ORF GC range. The others
# restrict which sequences are admissible without changing what GC any sequence
# can reach, which is why relaxing them cannot resolve an arithmetic GC conflict.
GC_RELEVANT = frozenset({"forbidden_codons", "gc_global_min"})


def analyse(campaign: pathlib.Path) -> dict:
    engine = ConstructEngine.load(campaign)
    spec = ConflictSpec.load(campaign)
    protein = read_protein(campaign, spec.protein_input)

    unrestricted = engine.achievable_gc_bounds(protein)
    restricted = engine.achievable_gc_bounds(protein, spec.forbidden_codons)

    if spec.expect not in ("infeasible", "feasible"):
        raise OracleFailure(f"conflict.expect must be infeasible or feasible, "
                            f"got {spec.expect!r}")
    if not restricted["feasible_codon_assignment"]:
        raise OracleFailure(
            "the forbidden-codon set starves a residue, so the conflict is a "
            "no-admissible-codon conflict rather than the intended GC conflict; "
            f"starved: {restricted['residues_without_codon']}")
    if spec.gc_global_min > unrestricted["gc_max"]:
        raise OracleFailure(
            f"gc_global_min {spec.gc_global_min} exceeds the unrestricted maximum "
            f"{unrestricted['gc_max']}, so the GC requirement is unsatisfiable on "
            "its own and the conflict is not a genuine PAIR conflict")

    # The declared expectation is checked against the arithmetic, so a campaign
    # cannot silently become the opposite of what it was authored to be.
    reachable = spec.gc_global_min <= restricted["gc_max"]
    if spec.expect == "infeasible" and reachable:
        raise OracleFailure(
            f"declared infeasible, but gc_global_min {spec.gc_global_min} is still "
            f"reachable under the forbidden-codon set (max {restricted['gc_max']}); "
            "the set is not provably infeasible and has no abstention gold")
    if spec.expect == "feasible" and not reachable:
        raise OracleFailure(
            f"declared feasible, but gc_global_min {spec.gc_global_min} exceeds the "
            f"restricted maximum {restricted['gc_max']}")

    # Feasibility is proved CONSTRUCTIVELY -- by exhibiting a design -- because
    # the achievable-GC arithmetic settles only whether the GC floor is
    # reachable, not whether the composition constraints leave any candidate.
    witness_design = None
    if spec.expect == "feasible":
        seq = engine.design_under_conflict_set(
            protein, spec.host, spec.strategy, spec.forbidden_codons,
            spec.gc_global_min)
        counts = engine.conflict_violations(seq, spec.host, spec.strategy,
                                             spec.forbidden_codons, spec.gc_global_min)
        if sum(counts.values()) != 0 or engine.translate(seq) != protein:
            raise OracleFailure(f"constructive witness is invalid: {counts}")
        witness_design = {
            "gc_fraction": round(engine.gc_fraction(seq), 6),
            "violations": counts,
            "encodes_declared_protein": True,
            "note": "Existence witness for the feasible verdict. Not graded: many "
                    "designs satisfy the set, and grading one would grade our "
                    "optimiser rather than the agent's reasoning.",
        }

    # Each single-constraint relaxation, by exact arithmetic. Only the two
    # GC-relevant constraints can move the bound.
    relaxations = {}
    for c in CONFLICT_CONSTRAINTS:
        if c == "forbidden_codons":
            gc_max = unrestricted["gc_max"]
            reachable = spec.gc_global_min <= gc_max + GC_TOLERANCE
            note = "all codons allowed again"
        elif c == "gc_global_min":
            gc_max = restricted["gc_max"]
            reachable = True
            note = "the GC requirement is removed, so it is vacuously met"
        else:
            gc_max = restricted["gc_max"]
            reachable = spec.gc_global_min <= gc_max + GC_TOLERANCE
            note = "does not change which codons are available, so the GC bound is unmoved"
        relaxations[c] = {"max_achievable_gc": gc_max,
                          "gc_requirement_reachable": bool(reachable),
                          "note": note}

    return {
        "protein_length_aa": len(protein) - 1,
        "orf_length_nt": 3 * len(protein),
        "host": spec.host,
        "cloning_strategy": spec.strategy,
        "forbidden_codons": list(spec.forbidden_codons),
        "gc_global_min": spec.gc_global_min,
        "unrestricted_gc_min": unrestricted["gc_min"],
        "unrestricted_gc_max": unrestricted["gc_max"],
        "restricted_gc_min": restricted["gc_min"],
        "restricted_gc_max": restricted["gc_max"],
        "admissible_codon_counts": restricted["admissible_codon_counts"],
        "verdict": ("infeasible_gc_unreachable" if spec.expect == "infeasible"
                    else "feasible"),
        "expect": spec.expect,
        "conflicting_constraints": (sorted(["forbidden_codons", "gc_global_min"])
                                    if spec.expect == "infeasible" else []),
        "constructive_witness": witness_design,
        "witness": {
            "required_gc_min": spec.gc_global_min,
            "max_achievable_gc_under_forbidden_codons": restricted["gc_max"],
            "margin": round(spec.gc_global_min - restricted["gc_max"], 6),
        },
        "each_constraint_satisfiable_alone": {
            "forbidden_codons": True,
            "gc_global_min": bool(spec.gc_global_min <= unrestricted["gc_max"]),
        },
        "single_constraint_relaxations": relaxations,
    }


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    # R2 starts from R1's deliverable: the constraint set restated, no bounds.
    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / "constraints_restated.json").write_text(json.dumps({
        "host": gold["host"], "cloning_strategy": gold["cloning_strategy"],
        "forbidden_codons": gold["forbidden_codons"],
        "gc_global_min": gold["gc_global_min"],
        "orf_length_nt": gold["orf_length_nt"],
    }, indent=2, sort_keys=True) + "\n")

    # R3 starts from R2's deliverable: the achievable bounds.
    (warm / "r3").mkdir(parents=True)
    (warm / "r3" / "bounds.json").write_text(json.dumps({
        "unrestricted_gc_min": gold["unrestricted_gc_min"],
        "unrestricted_gc_max": gold["unrestricted_gc_max"],
        "restricted_gc_min": gold["restricted_gc_min"],
        "restricted_gc_max": gold["restricted_gc_max"],
        "admissible_codon_counts": gold["admissible_codon_counts"],
    }, indent=2, sort_keys=True) + "\n")

    # R4 adds R3's verdict and proof -- and nothing about the relaxations.
    (warm / "r4").mkdir(parents=True)
    (warm / "r4" / "assessment.json").write_text(json.dumps({
        "unrestricted_gc_max": gold["unrestricted_gc_max"],
        "restricted_gc_max": gold["restricted_gc_max"],
        "verdict": gold["verdict"],
        "conflicting_constraints": gold["conflicting_constraints"],
        "witness": gold["witness"],
    }, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = load_task(campaign)
    overrides = ((task.get("grading") or {}).get("chance_levels") or {})
    chance = {**DEFAULT_CHANCE_LEVELS, **overrides}
    notes = ((task.get("grading") or {}).get("chance_level_notes") or {})

    def one(name, pred):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": 1}}

    def zero(name, pred):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": 0}}

    def ex(name, pred, g):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}"}}

    def tol(name, pred, g, abs_tol=1e-06):
        return {"name": name, "scorer": "abs_tol_match",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}",
                         "abs_tol": abs_tol}}

    r1 = [one("report_present", "present.report"),
          ex("orf_length_nt", "orf_length_nt", "orf_length_nt"),
          {"name": "constraints_restated", "scorer": "set_exact",
           "args": {"pred": "pred:constraints_assessed",
                    "gold": "gold:constraint_vocabulary"}},
          one("codon_counts_reported", "codon_counts_reported")]

    r2 = [tol("unrestricted_gc_max", "unrestricted_gc_max", "unrestricted_gc_max"),
          tol("unrestricted_gc_min", "unrestricted_gc_min", "unrestricted_gc_min"),
          tol("restricted_gc_max", "restricted_gc_max", "restricted_gc_max"),
          tol("restricted_gc_min", "restricted_gc_min", "restricted_gc_min"),
          ex("admissible_codon_counts", "admissible_codon_counts",
             "admissible_codon_counts")]

    r3 = [
        {"name": "verdict", "scorer": "enum_exact",
         "args": {"pred": "pred:verdict", "gold": "gold:verdict",
                  "vocabulary": list(VERDICTS)}},
        {"name": "conflicting_pair", "scorer": "set_exact",
         "args": {"pred": "pred:conflicting_constraints",
                  "gold": "gold:conflicting_constraints"}},
        tol("witness_required", "witness.required_gc_min",
            "witness.required_gc_min"),
        tol("witness_achievable", "witness.max_achievable_gc_under_forbidden_codons",
            "witness.max_achievable_gc_under_forbidden_codons"),
    ]
    if gold["expect"] == "infeasible":
        # Claiming a satisfying design when none can exist is fabrication, and
        # the behaviour the infeasible variant exists to catch. On the feasible
        # variant a design is admissible and ungraded, so the component is
        # omitted rather than inverted.
        r3.append(zero("no_design_claimed", "claims_satisfying_design"))

    r4 = []
    for c in CONFLICT_CONSTRAINTS:
        r4.append(tol(f"relax_{c}_gc_max", f"relaxations.{c}.max_achievable_gc",
                      f"single_constraint_relaxations.{c}.max_achievable_gc"))
        r4.append(ex(f"relax_{c}_reachable",
                     f"relaxations.{c}.gc_requirement_reachable",
                     f"single_constraint_relaxations.{c}.gc_requirement_reachable"))

    rungs = []
    for rid, ordinal, name, comps in (("r1", 1, "execute", r1), ("r2", 2, "correct", r2),
                                      ("r3", 3, "commit", r3),
                                      ("r4", 4, "counterfactual", r4)):
        rung = {"rung_id": rid, "ordinal": ordinal, "name": name, "gold_tier": "G1",
                "composition": "product", "pass_threshold": 1.0,
                "chance_level": chance[rid], "components": comps}
        if rid in notes:
            rung["chance_level_note"] = notes[rid]
        rungs.append(rung)

    spec = {"campaign_id": task["campaign_id"], "template_id": task["template_id"],
            "grading_spec_version": "1.0.0", "rungs": rungs}
    (pathlib.Path(campaign) / "grading.json").write_text(
        json.dumps(spec, indent=2) + "\n")
    return spec


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()

    copy_constants(campaign)
    ConflictSpec.load(campaign).write_reference(campaign)

    gold = analyse(campaign)
    gold["constraint_vocabulary"] = list(CONFLICT_CONSTRAINTS)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: {gold['protein_length_aa']} aa, host {gold['host']}")
    print(f"  achievable GC unrestricted [{gold['unrestricted_gc_min']}, "
          f"{gold['unrestricted_gc_max']}]")
    print(f"  achievable GC under forbidden codons [{gold['restricted_gc_min']}, "
          f"{gold['restricted_gc_max']}]")
    print(f"  required gc_global_min {gold['gc_global_min']} -> {gold['verdict']}")
    print(f"  conflicting pair {gold['conflicting_constraints']}, "
          f"margin {gold['witness']['margin']}")
    reachable = [c for c, v in gold["single_constraint_relaxations"].items()
                 if v["gc_requirement_reachable"]]
    print(f"  relaxations that resolve the conflict: {sorted(reachable)}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
