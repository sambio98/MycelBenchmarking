"""Reference solution for the chemical-space ladder.

Renders gold into the submission shape the output contract describes, which makes
the contract testable: if the oracle cannot express gold in the contracted shape,
the contract is wrong.

Usage: python -m npbench_c.templates.chemical_space.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.chemical_space.build import (
    CENSUS_KEYS,
    COMMITMENT_KEYS,
    GRANULARITY_KEYS,
)

REPORT = "chemical_space_report.json"


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in
              CENSUS_KEYS + GRANULARITY_KEYS + COMMITMENT_KEYS}
    report["threshold_variants"] = gold["threshold_variants"]
    report["evidence_gate"] = gold["evidence_gate"]
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
