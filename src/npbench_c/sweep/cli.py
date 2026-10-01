"""Run an internal sweep. `python -m npbench_c.sweep.cli <campaign> [--stubs]`"""

from __future__ import annotations

import argparse
import pathlib
import sys

from npbench_c.sweep.runner import SystemSpec
from npbench_c.sweep.sweep import REPEATS, sweep_campaign

STUB = [sys.executable, "-m", "npbench_c.sweep.stubs.stub_system"]


def stub_systems(campaign: pathlib.Path) -> list[SystemSpec]:
    """Harness self-test systems. Never part of a published result."""
    oracle = str((campaign / "oracle").resolve())
    import os

    env = {"PYTHONPATH": str(pathlib.Path("src").resolve()) + os.pathsep + os.environ.get("PYTHONPATH", "")}

    def spec(name, level, **kw):
        return SystemSpec(
            name=name,
            command=[*STUB, "--level", level, "--oracle", oracle, "--submission", "{submission}"],
            timeout_s=600, env=env, is_stub=True, **kw,
        )

    return [
        spec("stub-naive", "naive"),
        spec("stub-constrained", "constrained"),
        spec("stub-analyst", "analyst", is_reference=True),
        spec("stub-complete", "complete"),
        spec("stub-parametric", "parametric", mode="no_tool"),
    ]


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

    print(f"campaign: {report['_meta']['campaign_id']}  repeats={report['_meta']['repeats']}")
    print(f"chance floor: {report['_meta']['chance_floor']:.4f}")
    print(f"run status: {report['_meta']['run_status_counts']}")
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
    width = max(len(k) for k in report if not k.startswith("_"))
    for name, c in sorted((k, v) for k, v in report.items() if not k.startswith("_")):
        mark = "  ok  " if c["pass"] else " FAIL "
        print(f"[{mark}] {name:<{width}}  {c['detail']}")
    failed = [k for k, v in report.items() if not k.startswith("_") and not v["pass"]]
    print(f"\n{len(failed)} gate(s) failing: {failed or 'none'}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
