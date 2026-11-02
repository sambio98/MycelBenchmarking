"""Reference solution for the domain-architecture ladder.

Usage: python -m npbench_c.templates.domain_architecture.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.domain_architecture.build import (
    ARCHITECTURE_KEYS,
    COUNTERFACTUAL_KEYS,
    INVENTORY_KEYS,
    RECONCILIATION_KEYS,
)

REPORT = "architecture_report.json"


def render(gold: dict) -> dict:
    keys = (INVENTORY_KEYS + ARCHITECTURE_KEYS + RECONCILIATION_KEYS
            + COUNTERFACTUAL_KEYS)
    return {k: gold[k] for k in keys}


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
