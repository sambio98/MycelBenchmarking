"""Reference solution for the GCF-cutoff ladder.

Usage: python -m npbench_c.templates.gcf_cutoff.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.gcf_cutoff.build import (
    INVENTORY_KEYS,
    PARTITION_KEYS,
    RECONCILIATION_KEYS,
)

REPORT = "clustering_report.json"


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in INVENTORY_KEYS + PARTITION_KEYS
              + RECONCILIATION_KEYS}
    report["cutoff_sweep"] = gold["cutoff_sweep"]
    for key in ("best_cutoff", "best_pair_jaccard",
                "first_cutoff_with_a_false_join",
                "first_cutoff_absorbing_a_control",
                "first_cutoff_recovering_each_family"):
        report[key] = gold[key]
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
