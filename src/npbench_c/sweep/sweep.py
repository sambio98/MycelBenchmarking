"""The internal sweep: turn the readiness gate's five PENDING checks into numbers.

Run before any campaign reaches a collaborator. Expert attention is the
scarcest input in the project, so a campaign that will be dropped for
non-discrimination must be dropped by compute, not by an auditor who spent
three hours on it first.

Produces ``internal_sweep.json`` in the campaign directory, which the readiness
gate reads.
"""

from __future__ import annotations

import json
import pathlib
from collections import defaultdict
from typing import Sequence

import yaml

from npbench_c.grading.ladder import MIN_R1_CLEAR_RATE, check_monotonicity
from npbench_c.sweep.runner import (
    COMPLETED,
    INFRA_STATUSES,
    RunResult,
    SystemSpec,
    run_once,
)
from npbench_c.sweep.stats import (
    cluster_bootstrap_ci,
    leakage_z,
    mean,
    pass_at_1,
    pass_at_k,
)

REPEATS = 3
LEAKAGE_Z_MAX = 2.0
REFERENCE_LOWER_MIN = 0.50      # reference system must reach R1-R2
REFERENCE_UPPER_MAX = 0.50      # ...and must not reliably clear R4

__all__ = ["sweep_campaign", "REPEATS"]


def _run_summary(result: RunResult, rungs: Sequence[str]) -> dict:
    """Compact per-run record.

    Embedding each run's full grade would put ~190 KB per campaign into the
    committed artifact, regenerated on every sweep. The rung pass flags and the
    depth are what the gates and any later audit actually read; a full grade is
    reproducible from the sandbox at any time.
    """
    return {
        "system": result.system,
        "run_index": result.run_index,
        "rung_mode": result.rung_mode,
        "status": result.status,
        "depth": result.depth,
        "score": result.score,
        "rungs_passed": {r: result.cleared(r) for r in rungs},
        "monotonic": result.grade["monotonic"] if result.grade else None,
        "stderr_tail": result.stderr_tail,
    }


def _load_task(campaign: pathlib.Path) -> dict:
    return yaml.safe_load((campaign / "task.yaml").read_text())


def _rung_ids(task: dict) -> list[str]:
    return [r["rung_id"] for r in sorted(task["rungs"], key=lambda r: r["ordinal"])]


def _run_all(campaign: pathlib.Path, systems: Sequence[SystemSpec],
             sandbox_root: pathlib.Path, repeats: int) -> list[RunResult]:
    task = _load_task(campaign)
    measure_command = task["harness"]["measure_command"]
    results: list[RunResult] = []
    for system in systems:
        for i in range(repeats):
            results.append(run_once(campaign, system, i, sandbox_root,
                                    measure_command, "cold", None))
        # Warm runs give the per-rung difficulty profile. Cold runs cannot: a
        # failure at R2 censors every rung above it, so cold clear rates for
        # R3/R4 conflate "could not do it" with "never got there".
        for rung in task["rungs"]:
            bundle = rung.get("warm_start_bundle")
            if not bundle:
                continue
            for i in range(repeats):
                results.append(run_once(campaign, system, i, sandbox_root,
                                        measure_command, f"warm:{rung['rung_id']}",
                                        bundle))
    return results


def sweep_campaign(
    campaign_dir: pathlib.Path,
    systems: Sequence[SystemSpec],
    sandbox_root: pathlib.Path,
    repeats: int = REPEATS,
    write: bool = True,
) -> dict:
    campaign = pathlib.Path(campaign_dir).resolve()
    task = _load_task(campaign)
    rungs = _rung_ids(task)
    spec = json.loads((campaign / "grading.json").read_text())

    # The chance floor is recomputed from the grading spec rather than taken
    # from a grade, so the sweep's leakage comparison does not depend on any
    # particular run having succeeded.
    from npbench_c.grading.ladder import chance_floor
    chance = chance_floor([float(r.get("chance_level", 0.0)) for r in
                           sorted(spec["rungs"], key=lambda x: x["ordinal"])])
    results = _run_all(campaign, systems, sandbox_root, repeats)

    by_name = {s.name: s for s in systems}
    cold = [r for r in results if r.rung_mode == "cold"]
    completed_cold = [r for r in cold if r.status == COMPLETED]

    # Infra failures are reported, never silently scored as agent failures:
    # otherwise part of the headline is measuring our own container.
    infra = [r for r in results if r.status in INFRA_STATUSES]
    status_counts: dict[str, int] = defaultdict(int)
    for r in results:
        status_counts[r.status] += 1

    # --- per-rung clear rates (warm where available, cold for r1) -----------
    warm_rates: dict[str, float] = {}
    for rung in rungs:
        mode = "cold" if rung == rungs[0] else f"warm:{rung}"
        runs = [r for r in results
                if r.rung_mode == mode and r.status == COMPLETED
                and by_name[r.system].mode == "tooled"]
        if runs:
            warm_rates[rung] = pass_at_1([r.cleared(rung) for r in runs])

    # --- per-system profile -------------------------------------------------
    profile: dict[str, dict] = {}
    for system in systems:
        sys_cold = [r for r in cold if r.system == system.name and r.status == COMPLETED]
        per_rung = {}
        for rung in rungs:
            mode = "cold" if rung == rungs[0] else f"warm:{rung}"
            runs = [r for r in results if r.system == system.name
                    and r.rung_mode == mode and r.status == COMPLETED]
            if runs:
                per_rung[rung] = {
                    "pass_at_1": pass_at_1([r.cleared(rung) for r in runs]),
                    "pass_at_k": pass_at_k([[r.cleared(rung) for r in runs]]),
                }
        profile[system.name] = {
            "mode": system.mode,
            "is_reference": system.is_reference,
            "cold_mean_score": mean([r.score for r in sys_cold]) if sys_cold else None,
            "cold_depths": sorted(r.depth for r in sys_cold),
            "per_rung": per_rung,
        }

    # --- gates --------------------------------------------------------------
    checks: dict[str, dict] = {}

    rate_list = [warm_rates[r] for r in rungs if r in warm_rates]
    checks["monotonicity_measured"] = {
        "pass": bool(rate_list) and check_monotonicity(rate_list),
        "detail": f"clear rates by rung {dict(zip(rungs, rate_list))}; "
                  "non-increasing required (warm for r2+, cold for r1)",
    }

    generic = [r for r in cold if r.status == COMPLETED
               and by_name[r.system].mode == "tooled"
               and not by_name[r.system].is_reference]
    r1_rate = pass_at_1([r.cleared(rungs[0]) for r in generic]) if generic else None
    checks["r1_clear_rate"] = {
        "pass": r1_rate is not None and r1_rate >= MIN_R1_CLEAR_RATE,
        "detail": f"R1 pass@1 for generic tooled systems = {r1_rate}; "
                  f"minimum {MIN_R1_CLEAR_RATE} (below this the campaign measures "
                  "environment setup, not science)",
    }

    ablation = [r for r in cold if by_name[r.system].mode == "no_tool"
                and r.status == COMPLETED]
    if ablation:
        lk, se, z = leakage_z([r.score for r in ablation], chance)
        checks["leakage_no_tool_ablation"] = {
            "pass": z <= LEAKAGE_Z_MAX,
            "detail": f"no-tool score {lk:.4f} vs chance floor {chance:.4f}, "
                      f"SE {se:.4f}, z {z:.2f}; kill above z={LEAKAGE_Z_MAX}",
        }
    else:
        checks["leakage_no_tool_ablation"] = {
            "pass": False, "detail": "no no_tool ablation system was run",
        }

    ref = next((s for s in systems if s.is_reference), None)
    if ref and profile[ref.name]["per_rung"]:
        pr = profile[ref.name]["per_rung"]
        lower = [pr[r]["pass_at_1"] for r in rungs[:2] if r in pr]
        upper = [pr[r]["pass_at_1"] for r in rungs[2:] if r in pr]
        ok = (all(v >= REFERENCE_LOWER_MIN for v in lower)
              and bool(upper) and min(upper) <= REFERENCE_UPPER_MAX)
        checks["difficulty_gate"] = {
            "pass": ok,
            "detail": f"reference '{ref.name}' lower rungs {lower} "
                      f"(each >= {REFERENCE_LOWER_MIN}), upper rungs {upper} "
                      f"(min <= {REFERENCE_UPPER_MAX}); difficulty must sit at R3/R4",
        }
    else:
        checks["difficulty_gate"] = {
            "pass": False, "detail": "no reference system was run",
        }

    scored = [r for r in completed_cold if by_name[r.system].mode == "tooled"]
    depths = sorted({r.depth for r in scored})
    mean_score = mean([r.score for r in scored]) if scored else None
    dead_rungs = [r for r in rungs if warm_rates.get(r) == 0.0]
    checks["discrimination"] = {
        "pass": bool(scored) and len(depths) >= 2
                and mean_score is not None and 0.0 < mean_score < 1.0,
        "detail": f"distinct cold depths {depths}, mean score {mean_score}; "
                  f"needs >=2 depths and 0 < mean < 1"
                  + (f"; rungs never cleared by any system: {dead_rungs}" if dead_rungs else ""),
    }

    # --- headline -----------------------------------------------------------
    clusters = defaultdict(list)
    for r in scored:
        clusters[spec["campaign_id"]].append(r.score)
    point, lo, hi = cluster_bootstrap_ci(clusters) if clusters else (None, None, None)

    report = {
        **checks,
        "_meta": {
            "campaign_id": spec["campaign_id"],
            "repeats": repeats,
            "chance_floor": chance,
            "systems": [s.name for s in systems],
            "stub_systems": [s.name for s in systems if s.is_stub],
            "stub_based": any(s.is_stub for s in systems),
            "run_status_counts": dict(sorted(status_counts.items())),
            "infra_failures": len(infra),
            "note": "Scores exclude infra failures (timeout/crash); run_status_counts "
                    "reports them so the headline is not measuring the container.",
            "stub_warning": "A stub-based sweep validates the harness, not the "
                            "campaign: stub systems are built to reach preset "
                            "depths, so difficulty, discrimination and "
                            "monotonicity numbers from one are fixtures. The "
                            "readiness gate keeps these checks PENDING until a "
                            "sweep runs real systems.",
        },
        "_headline": {
            "cold_mean_score": point,
            "ci_low": lo,
            "ci_high": hi,
            "ci_note": "single-campaign cluster bootstrap; intervals tighten only "
                       "with more campaigns, not more repeats",
        },
        "_per_rung_clear_rates": warm_rates,
        "_per_system": profile,
        "_runs": [_run_summary(r, rungs) for r in results],
    }
    if write:
        (campaign / "internal_sweep.json").write_text(
            json.dumps(report, sort_keys=True, indent=2, default=str) + "\n")
    return report
