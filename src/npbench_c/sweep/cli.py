"""Run an internal sweep. `python -m npbench_c.sweep.cli <campaign> [--stubs]`"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys

import yaml

from npbench_c.sweep.runner import SystemSpec
from npbench_c.sweep.sweep import REPEATS, sweep_campaign


def stub_systems(campaign: pathlib.Path) -> list[SystemSpec]:
    """Harness self-test systems, declared by the campaign in task.yaml.

    Never part of a published result: a stub-based sweep marks itself and the
    readiness gate keeps every agent-run check PENDING.
    """
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    harness = task.get("harness") or {}
    module = harness.get("stub_module")
    levels = harness.get("stub_levels") or []
    if not module or not levels:
        raise SystemExit(
            f"{campaign.name}: task.yaml declares no harness.stub_module/stub_levels"
        )

    subs = {"oracle": str((campaign / "oracle").resolve())}
    extra = [a.format(**subs) for a in (harness.get("stub_args") or [])]
    env = {"PYTHONPATH": os.pathsep.join(
        [str(pathlib.Path("src").resolve()), os.environ.get("PYTHONPATH", "")]
    )}

    systems = []
    for spec in levels:
        systems.append(SystemSpec(
            name=spec["name"],
            command=[sys.executable, "-m", module, "--level", spec["level"],
                     *extra, "--submission", "{submission}"],
            mode=spec.get("mode", "tooled"),
            is_reference=bool(spec.get("reference")),
            is_stub=True,
            timeout_s=int(spec.get("timeout_s", 600)),
            env=env,
        ))
    return systems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign")
    ap.add_argument("--stubs", action="store_true",
                    help="run the harness self-test systems instead of real agents")
    ap.add_argument("--repeats", type=int, default=REPEATS)
    ap.add_argument("--sandbox-root", default=None)
    args = ap.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    if not args.stubs:
        print("no real systems configured; pass --stubs to self-test the harness",
              file=sys.stderr)
        return 2

    sandbox = pathlib.Path(args.sandbox_root) if args.sandbox_root else \
        campaign.parent.parent / ".sweep_sandbox" / campaign.name
    report = sweep_campaign(campaign, stub_systems(campaign), sandbox, args.repeats)

    meta = report["_meta"]
    print(f"campaign: {meta['campaign_id']}  repeats={meta['repeats']}")
    print(f"chance floor: {meta['chance_floor']:.4f}")
    print(f"run status: {meta['run_status_counts']}")
    h = report["_headline"]
    if h["cold_mean_score"] is not None:
        print(f"cold mean score: {h['cold_mean_score']:.4f} "
              f"CI [{h['ci_low']:.4f}, {h['ci_high']:.4f}]")
    print("\nper-rung clear rates (tooled systems):")
    for rung, rate in report["_per_rung_clear_rates"].items():
        print(f"  {rung}: {rate:.3f}")
    print("\nper-system cold depths:")
    for name, p in report["_per_system"].items():
        print(f"  {name:<18} mode={p['mode']:<8} depths={p['cold_depths']}")
    print("\ngates:")
    gates = sorted((k, v) for k, v in report.items() if not k.startswith("_"))
    width = max(len(k) for k, _ in gates)
    for name, c in gates:
        mark = "  ok  " if c["pass"] else " FAIL "
        print(f"[{mark}] {name:<{width}}  {c['detail']}")
    failed = [k for k, v in gates if not v["pass"]]
    print(f"\n{len(failed)} gate(s) failing: {failed or 'none'}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
