"""Reference solution for the selectivity ladder.

Usage: python -m npbench_c.templates.chembl_selectivity.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.chembl_selectivity.build import (
    ADMISSION_KEYS,
    CENSUS_KEYS,
    VERDICT_KEYS,
)

REPORT = "selectivity_report.json"


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in CENSUS_KEYS + ADMISSION_KEYS + VERDICT_KEYS}
    report["admission_variants"] = gold["admission_variants"]
    report["ratio_extremes"] = gold["ratio_extremes"]
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
