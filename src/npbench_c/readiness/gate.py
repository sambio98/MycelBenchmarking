"""The readiness gate: may this campaign leave the building?

A campaign is not ready because someone believes it is; it is ready when this
returns all-pass. Checks that require internal agent runs cannot be satisfied
mechanically and report PENDING rather than PASS -- the gate never launders an
unmeasured property into a green tick. Collaborator attention is the scarcest
input in the project, so nothing reaches an auditor until the mechanical checks
are green and the agent-run checks have real numbers.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys
from dataclasses import dataclass

from npbench_c.grading.grade import grade, grade_json
from npbench_c.sweep.runner import _absolute_pythonpath
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.grading.version import grader_version

PURITY_REPEATS = 100
FORBIDDEN_SOURCES = ("kegg", "metacyc", "biocyc", "drugbank")
VENDOR_TOKENS = ("mycel",)

PASS, FAIL, PENDING = "PASS", "FAIL", "PENDING"

AGENT_RUN_CHECKS = (
    "monotonicity_measured",
    "r1_clear_rate",
    "leakage_no_tool_ablation",
    "difficulty_gate",
    "discrimination",
)


@dataclass
class Check:
    name: str
    status: str
    detail: str


def _check(name: str, ok: bool, detail: str) -> Check:
    return Check(name, PASS if ok else FAIL, detail)


def _oracle_and_purity(campaign: pathlib.Path) -> list[Check]:
    import yaml

    checks = []
    sub = campaign / "oracle_submission"
    measurement = sub / "measurement.json"

    # The measure command comes from the campaign's own task.yaml, exactly as it
    # does for the sweep runner. The gate previously globbed for a *.faa input,
    # which is a construct-design assumption that broke the moment a campaign
    # had different inputs.
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    harness = task.get("harness") or {}
    measure_command = harness.get("measure_command")
    if not measure_command:
        return [Check("oracle_runs", FAIL,
                      "task.yaml declares no harness.measure_command")]

    # A campaign built from a shared template carries no oracle directory at
    # all, so the solve step is whatever task.yaml declares. The legacy form --
    # a campaign-local oracle/solve.py -- is still honoured.
    solve_command = harness.get("solve_command") or [sys.executable, "oracle/solve.py",
                                                     "{submission}"]
    env = _absolute_pythonpath(os.environ)
    try:
        subprocess.run(
            [tok.format(campaign=str(campaign), submission=str(sub)) for tok in solve_command],
            cwd=campaign, env=env, check=True, capture_output=True, timeout=1800,
        )
        subprocess.run(
            [tok.format(campaign=str(campaign), submission=str(sub),
                        measurement=str(measurement)) for tok in measure_command],
            cwd=campaign, env=env, check=True, capture_output=True, timeout=1800,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
        return [Check("oracle_runs", FAIL, f"oracle pipeline failed: {exc}")]

    measurement_data = json.loads(measurement.read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    spec = json.loads((campaign / "grading.json").read_text())
    result = grade(measurement_data, gold, spec)

    checks.append(_check("oracle_scores_1.0", result["score"] == 1.0,
                         f"score={result['score']} depth={result['depth']}"))
    checks.append(_check("oracle_ladder_monotonic", result["monotonic"],
                         f"rungs_cleared={result['rungs_cleared']} depth={result['depth']}"))
    # The cap applies to the discriminating rungs (R2 upward). R1 is an
    # execution precondition whose design target is an 85-95% clear rate, so
    # including it would make the cap unsatisfiable for every well-formed ladder.
    checks.append(_check("chance_floor_within_limit",
                         result["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR,
                         f"discriminating={result['discriminating_chance_floor']} "
                         f"(full={result['chance_floor']}) limit={MAX_CHANCE_FLOOR}"))

    reference = grade_json(measurement_data, gold, spec)
    identical = all(grade_json(measurement_data, gold, spec) == reference
                    for _ in range(PURITY_REPEATS))
    checks.append(_check("grader_purity", identical,
                         f"{PURITY_REPEATS} repeat gradings byte-identical={identical}"))
    checks.append(_check("grader_version_recorded",
                         spec.get("campaign_id") is not None,
                         f"grader_version={grader_version()}"))

    # The pinned version must match the grader that just ran. Without this the
    # campaign can silently carry a stale hash: a grader change is a benchmark
    # version bump, and a score whose grader_version does not identify the code
    # that produced it is not a comparable score.
    import yaml

    task = yaml.safe_load((campaign / "task.yaml").read_text())
    pinned, current = task.get("grader_version"), grader_version()
    checks.append(_check("grader_version_matches_pin", pinned == current,
                         f"pinned={pinned} current={current}"))
    return checks


def _oracle_determinism(campaign: pathlib.Path) -> Check:
    """Regenerate gold and require byte-identity with what is committed."""
    import yaml

    task = yaml.safe_load((campaign / "task.yaml").read_text())
    gold_command = ((task.get("harness") or {}).get("gold_command")
                    or [sys.executable, "oracle/generate_gold.py"])
    gold_path = campaign / "gold" / "gold.json"
    before = gold_path.read_bytes()
    try:
        subprocess.run([tok.format(campaign=str(campaign)) for tok in gold_command],
                       cwd=campaign, env=_absolute_pythonpath(os.environ),
                       check=True, capture_output=True, timeout=3600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
        return Check("oracle_determinism", FAIL, f"regeneration failed: {exc}")
    after = gold_path.read_bytes()
    if before != after:
        gold_path.write_bytes(before)
        return Check("oracle_determinism", FAIL, "regenerated gold differs from committed gold")
    return Check("oracle_determinism", PASS, "gold byte-identical on regeneration")


def _contract_and_provenance(campaign: pathlib.Path) -> list[Check]:
    import yaml  # declared dependency; not in the scoring path

    checks = []
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    spec = json.loads((campaign / "grading.json").read_text())

    tiers = [r.get("gold_tier") for r in spec["rungs"]]
    checks.append(_check("gold_tier_declared_per_rung", all(tiers),
                         f"tiers={tiers}"))
    checks.append(_check("at_most_one_g2_rung", tiers.count("G2") <= 1,
                         f"G2 rungs={tiers.count('G2')}"))
    checks.append(_check("no_g3_outside_repro_audit",
                         "G3" not in tiers or task.get("family") == "reproducibility_audit",
                         f"tiers={tiers} family={task.get('family')}"))
    checks.append(_check("license_class_declared",
                         task.get("license_class") in ("open", "nc"),
                         f"license_class={task.get('license_class')}"))

    blob = json.dumps(task).lower()
    bad = [s for s in FORBIDDEN_SOURCES if s in blob]
    checks.append(_check("no_license_landmine_sources", not bad,
                         f"forbidden sources referenced: {bad}" if bad else "none"))

    contract = task.get("output_contract") or {}
    checks.append(_check("output_contract_present", bool(contract),
                         f"{len(contract)} artifacts contracted"))

    # Every pinned table the task requires must be inside a path the sandbox
    # ships, or the campaign is unsolvable and no run can tell you so: the
    # agent simply fails and looks incapable.
    from npbench_c.sweep.runner import PUBLIC_CAMPAIGN_PATHS

    pinned = task.get("constraints", {}).get("pinned_tables") or []
    pinned = [pinned] if isinstance(pinned, str) else pinned
    reachable = [p for p in pinned
                 if pathlib.PurePath(p).parts[0] in PUBLIC_CAMPAIGN_PATHS
                 and (campaign / p).exists()]
    checks.append(_check("pinned_tables_reachable_by_agent",
                         bool(pinned) and len(reachable) == len(pinned),
                         f"{len(reachable)}/{len(pinned)} pinned tables are under a "
                         f"public sandbox path and present on disk"))

    # Every rung must carry the key so the declaration is explicit; R1 has no
    # earlier rung to start from, so null is the correct value there. Bundles
    # above R1 must exist on disk: a missing bundle silently turns a warm run
    # into a cold one and the per-rung difficulty profile becomes a fiction.
    rungs = sorted(task.get("rungs", []), key=lambda r: r["ordinal"])
    declared = [("warm_start_bundle" in r) for r in rungs]
    bundles_ok = all(
        (r.get("warm_start_bundle") is None) if i == 0
        else bool(r.get("warm_start_bundle"))
        and (campaign / r["warm_start_bundle"]).is_dir()
        for i, r in enumerate(rungs)
    )
    checks.append(_check("warm_start_declared_per_rung",
                         bool(rungs) and all(declared) and bundles_ok,
                         f"{sum(declared)}/{len(rungs)} rungs declare the key; "
                         f"bundles present and non-empty above R1: {bundles_ok}"))

    checks.append(_check("excluded_fields_justified",
                         all("reason" in e for e in task.get("excluded_from_grading", [])),
                         f"{len(task.get('excluded_from_grading', []))} exclusions, each with a measurement"))
    return checks


def _integrity(campaign: pathlib.Path) -> list[Check]:
    hits = []
    for path in sorted(campaign.rglob("*")):
        if not path.is_file() or path.suffix not in (".py", ".yaml", ".json", ".md"):
            continue
        text = path.read_text(errors="replace").lower()
        for token in VENDOR_TOKENS:
            if token in text:
                hits.append(f"{path.relative_to(campaign)}:{token}")
    return [_check("no_vendor_tool_names", not hits, f"hits={hits}" if hits else "none")]


def _pending_agent_runs(campaign: pathlib.Path) -> list[Check]:
    """Report measured values from an internal sweep if present, else PENDING."""
    sweep = campaign / "internal_sweep.json"
    if not sweep.is_file():
        return [Check(name, PENDING, "no internal_sweep.json; requires internal agent runs")
                for name in AGENT_RUN_CHECKS]
    data = json.loads(sweep.read_text())

    # A stub-based sweep is a harness self-test. Its systems are built to reach
    # preset depths, so promoting its numbers to PASS would launder a fixture
    # into evidence -- the exact thing this gate exists to prevent.
    stub_based = bool(data.get("_meta", {}).get("stub_based"))
    stubs = data.get("_meta", {}).get("stub_systems", [])

    out = []
    for name in AGENT_RUN_CHECKS:
        if name not in data:
            out.append(Check(name, PENDING, "absent from internal_sweep.json"))
            continue
        entry = data[name]
        detail = str(entry.get("detail", ""))
        if stub_based:
            out.append(Check(name, PENDING,
                             f"stub-based sweep (fixtures: {', '.join(stubs)}); "
                             f"real systems required. Harness result: {detail}"))
        else:
            out.append(Check(name, PASS if entry.get("pass") else FAIL, detail))
    return out


def run(campaign_dir: pathlib.Path) -> list[Check]:
    campaign = pathlib.Path(campaign_dir).resolve()
    checks: list[Check] = []
    checks += _oracle_and_purity(campaign)
    checks.append(_oracle_determinism(campaign))
    checks += _contract_and_provenance(campaign)
    checks += _integrity(campaign)
    checks += _pending_agent_runs(campaign)
    return checks


def report(campaign_dir: pathlib.Path) -> dict:
    checks = run(campaign_dir)
    counts = {s: sum(1 for c in checks if c.status == s) for s in (PASS, FAIL, PENDING)}
    return {
        "campaign_dir": str(campaign_dir),
        "grader_version": grader_version(),
        "ready_to_audit": counts[FAIL] == 0 and counts[PENDING] == 0,
        "mechanical_checks_green": counts[FAIL] == 0,
        "summary": counts,
        "checks": [{"name": c.name, "status": c.status, "detail": c.detail} for c in checks],
    }


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print("usage: python -m npbench_c.readiness.gate <campaign_dir>", file=sys.stderr)
        return 2
    rep = report(pathlib.Path(argv[0]))
    width = max(len(c["name"]) for c in rep["checks"])
    for c in rep["checks"]:
        mark = {"PASS": "  ok  ", "FAIL": " FAIL ", "PENDING": " pend "}[c["status"]]
        print(f"[{mark}] {c['name']:<{width}}  {c['detail']}")
    s = rep["summary"]
    print(f"\n{s['PASS']} pass / {s['FAIL']} fail / {s['PENDING']} pending")
    print(f"mechanical checks green: {rep['mechanical_checks_green']}")
    print(f"ready to send for audit: {rep['ready_to_audit']}")
    return 0 if rep["mechanical_checks_green"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
