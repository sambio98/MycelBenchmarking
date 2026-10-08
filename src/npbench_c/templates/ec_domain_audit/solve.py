"""Reference solution for the EC domain-audit ladder.

Usage: python -m npbench_c.templates.ec_domain_audit.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.ec_domain_audit.build import (
    ARCHITECTURE_KEYS,
    INVENTORY_KEYS,
    PRESENCE_KEYS,
)

REPORT = "domain_audit.json"


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in INVENTORY_KEYS + PRESENCE_KEYS
              + ARCHITECTURE_KEYS}
    report["misannotated_by_policy"] = gold["misannotated_by_policy"]
    report["threshold_absent"] = gold["threshold_absent"]
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
