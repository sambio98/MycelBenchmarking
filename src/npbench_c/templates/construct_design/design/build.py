"""Build gold, warm bundles and the grading spec for the design ladder.

Usage: python -m npbench_c.templates.construct_design.design.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

from npbench_c.templates.construct_design.core import (
    OPTIMIZER_POLICIES,
    ConstructEngine,
    OracleFailure,
)
from npbench_c.templates.construct_design.params import (
    copy_constants,
    load_task,
    read_protein,
)

CAI_FLOOR_MARGIN = 0.025
DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.02, "r3": 0.001, "r4": 0.001}


def analyse(engine: ConstructEngine, protein: str, host: str, strategy: str) -> dict:
    """Per-host analysis, with the cross-optimiser spread that decides what is
    gradable.

    The determinism audit measures run-to-run spread, which is zero for a
    deterministic oracle. What matters for a quantity the agent computes with its
    own method is CROSS-METHOD spread, so several legitimate optimiser policies
    are run and their disagreement recorded. The binding constraint is stable
    across all of them and is graded; the CAI magnitudes are not and are graded
    as a floor instead.
    """
    opt = engine.unconstrained_optimum(protein, host)
    uv = engine.violation_counts(opt, host, strategy)
    per_policy = {p: engine.constraint_costs(protein, host, strategy, p)
                  for p in OPTIMIZER_POLICIES}
    bindings = {a["binding_constraint"] for a in per_policy.values()}
    if len(bindings) != 1:
        raise OracleFailure(
            f"binding constraint unstable across optimiser policies: {bindings}. "
            "Gold is not defensible; redesign the rung.")
    cais = [a["cai_full"] for a in per_policy.values()]
    non_binding = sorted(c for c in engine.constraint_order
                         if all(a["costs"][c] <= 0.0 for a in per_policy.values()))
    ref = per_policy["first_improve"]
    return {
        "unconstrained_violations": {k: uv[k] for k in engine.constraint_order},
        "binding_constraint": bindings.pop(),
        "non_binding_constraints": non_binding,
        "reference_sequence": ref["sequence"],
        "reference_cai": ref["cai_full"],
        "cai_min_across_policies": round(min(cais), 6),
        "cai_max_across_policies": round(max(cais), 6),
        "cai_spread_across_policies": round(max(cais) - min(cais), 6),
        "cost_spread_across_policies": {
            c: round(max(a["costs"][c] for a in per_policy.values())
                     - min(a["costs"][c] for a in per_policy.values()), 6)
            for c in engine.constraint_order
        },
    }


def build_gold(campaign: pathlib.Path) -> dict:
    engine = ConstructEngine.load(campaign)
    params = load_task(campaign)["template_params"]
    protein = read_protein(campaign, params["protein_input"])
    gold: dict = {"protein_length_aa": len(protein) - 1, "hosts": {}}

    for key, spec in sorted(params["hosts"].items()):
        a = analyse(engine, protein, spec["host"], spec["strategy"])
        seq = a.pop("reference_sequence")
        if sum(engine.violation_counts(seq, spec["host"], spec["strategy"]).values()) != 0:
            raise OracleFailure(f"{key}: reference sequence violates constraints")
        if engine.translate(seq) != protein:
            raise OracleFailure(f"{key}: reference sequence does not encode the protein")
        # Floor set below the worst legitimate optimiser policy, so an agent is
        # never failed for using a different valid optimiser.
        a["cai_floor"] = round(a["cai_min_across_policies"] - CAI_FLOOR_MARGIN, 3)
        gold["hosts"][key] = a
        (campaign / "gold").mkdir(exist_ok=True)
        (campaign / "gold" / f"orf_{key}.fasta").write_text(
            f">gold_{key}_orf\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")

    keys = sorted(gold["hosts"])
    if len(keys) != 2:
        raise OracleFailure("the design ladder's R4 is a host swap and needs exactly two hosts")
    a, b = gold["hosts"][keys[0]], gold["hosts"][keys[1]]
    gold["flipped_constraints"] = sorted(
        set(a["non_binding_constraints"]) ^ set(b["non_binding_constraints"]))
    gold["binding_constraint_flips"] = int(
        a["binding_constraint"] != b["binding_constraint"])
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    return gold


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    """Pre-rendered bundles, never whole gold files: copying gold.json to supply
    the previous answer would also supply the answer being asked for."""
    engine = ConstructEngine.load(campaign)
    params = load_task(campaign)["template_params"]
    protein = read_protein(campaign, params["protein_input"])
    primary, secondary = sorted(params["hosts"])

    warm = campaign / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)

    def fasta(seq: str, name: str) -> str:
        return f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"

    # R2 starts from R1's deliverable: a schema-valid ORF that does not yet
    # satisfy the constraints. The naive reverse translation is exactly that.
    (warm / "r2").mkdir(parents=True)
    (warm / "r2" / f"orf_{primary}_start.fasta").write_text(
        fasta(engine.unconstrained_optimum(protein, params["hosts"][primary]["host"]),
              "r1_schema_valid_start"))

    # R3 starts from the correct primary-host construct.
    (warm / "r3").mkdir(parents=True)
    shutil.copyfile(campaign / "gold" / f"orf_{primary}.fasta",
                    warm / "r3" / f"orf_{primary}.fasta")

    # R4 adds R3's answer for the primary host only -- never the secondary's.
    (warm / "r4").mkdir(parents=True)
    shutil.copyfile(campaign / "gold" / f"orf_{primary}.fasta",
                    warm / "r4" / f"orf_{primary}.fasta")
    p = gold["hosts"][primary]
    (warm / "r4" / f"analysis_{primary}.json").write_text(json.dumps({
        "unconstrained_violations": p["unconstrained_violations"],
        "binding_constraint": p["binding_constraint"],
        "non_binding_constraints": p["non_binding_constraints"],
    }, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    engine = ConstructEngine.load(campaign)
    task = load_task(campaign)
    params = task["template_params"]
    primary, secondary = sorted(params["hosts"])
    overrides = ((task.get("grading") or {}).get("chance_levels") or {})
    chance = {**DEFAULT_CHANCE_LEVELS, **overrides}
    notes = ((task.get("grading") or {}).get("chance_level_notes") or {})
    constraints = list(engine.constraint_order)
    schema_flags = ["valid_alphabet", "in_frame", "starts_atg", "ends_stop",
                    "no_internal_stop", "length_matches_protein"]

    def one(name, pred):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": 1}}

    def ex(name, pred, g):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{pred}", "gold": f"gold:{g}"}}

    r1 = [one("orf_present", f"present.orf_{primary}"),
          one("analysis_present", "present.analysis")]
    r1 += [one(f"schema_{k}", f"{primary}.schema.{k}") for k in schema_flags]

    r2 = [one("protein_identity", f"{primary}.protein_identity"),
          one("all_constraints_satisfied", f"{primary}.all_constraints_satisfied"),
          {"name": "cai_floor", "scorer": "floor_match",
           "args": {"pred": f"pred:{primary}.cai",
                    "floor": f"gold:hosts.{primary}.cai_floor"}}]
    r2 += [{"name": f"zero_{c}", "scorer": "exact",
            "args": {"pred": f"pred:{primary}.violations.{c}", "gold": 0}}
           for c in constraints]

    r3 = [ex(f"uv_{c}", f"claims.r3.unconstrained_violations.{c}",
             f"hosts.{primary}.unconstrained_violations.{c}") for c in constraints]
    r3 += [
        {"name": "binding_constraint", "scorer": "enum_exact",
         "args": {"pred": "pred:claims.r3.binding_constraint",
                  "gold": f"gold:hosts.{primary}.binding_constraint",
                  "vocabulary": constraints}},
        {"name": "non_binding_set", "scorer": "set_exact",
         "args": {"pred": "pred:claims.r3.non_binding_constraints",
                  "gold": f"gold:hosts.{primary}.non_binding_constraints"}},
    ]

    r4 = [one("swap_orf_present", f"present.orf_{secondary}"),
          one("swap_protein_identity", f"{secondary}.protein_identity"),
          one("swap_constraints_satisfied", f"{secondary}.all_constraints_satisfied"),
          {"name": "swap_cai_floor", "scorer": "floor_match",
           "args": {"pred": f"pred:{secondary}.cai",
                    "floor": f"gold:hosts.{secondary}.cai_floor"}}]
    r4 += [ex(f"swap_uv_{c}", f"claims.r4.unconstrained_violations.{c}",
              f"hosts.{secondary}.unconstrained_violations.{c}") for c in constraints]
    r4 += [
        {"name": "swap_binding_constraint", "scorer": "enum_exact",
         "args": {"pred": "pred:claims.r4.binding_constraint",
                  "gold": f"gold:hosts.{secondary}.binding_constraint",
                  "vocabulary": constraints}},
        {"name": "swap_non_binding_set", "scorer": "set_exact",
         "args": {"pred": "pred:claims.r4.non_binding_constraints",
                  "gold": f"gold:hosts.{secondary}.non_binding_constraints"}},
        {"name": "flipped_constraints", "scorer": "set_exact",
         "args": {"pred": "pred:claims.r4.flipped_constraints",
                  "gold": "gold:flipped_constraints"}},
    ]

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
            "grading_spec_version": "1.1.0", "rungs": rungs}
    (campaign / "grading.json").write_text(json.dumps(spec, indent=2) + "\n")
    return spec


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)
    gold = build_gold(campaign)
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)
    for key, block in sorted(gold["hosts"].items()):
        print(f"  {key:<16} binding={block['binding_constraint']:<16} "
              f"cai_floor={block['cai_floor']} "
              f"spread={block['cai_spread_across_policies']}")
    print(f"  flipped={gold['flipped_constraints']} "
          f"flips={gold['binding_constraint_flips']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
