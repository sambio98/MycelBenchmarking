"""Reference solution for the annotation-transfer ladder.

Usage: python -m npbench_c.templates.annotation_transfer.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.annotation_transfer.build import (
    SEARCH_KEYS,
    SET_KEYS,
    VERDICT_KEYS,
)

REPORT = "transfer_report.json"


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in SET_KEYS + SEARCH_KEYS + VERDICT_KEYS}
    report["removal_census"] = gold["removal_census"]
    report["removal_changed"] = gold["removal_changed"]
    report["removal_transitions"] = gold["removal_transitions"]
    report["constant_depth_agreement"] = gold["constant_depth_agreement"]
    return report


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: solve <campaign_dir> <submission_dir>", file=sys.stderr)
        return 2
    campaign, out = pathlib.Path(argv[0]), pathlib.Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    (out / REPORT).write_text(json.dumps(render(gold), indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
