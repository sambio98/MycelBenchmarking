"""Reference solution for the design ladder.

Must score exactly 1.0 through the real grader or the campaign does not ship.

Usage:
  python -m npbench_c.templates.construct_design.design.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

from npbench_c.templates.construct_design.params import load_task


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: solve <campaign_dir> <submission_dir>", file=sys.stderr)
        return 2
    campaign, out = pathlib.Path(argv[0]), pathlib.Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    params = load_task(campaign)["template_params"]
    primary, secondary = sorted(params["hosts"])

    for key in sorted(params["hosts"]):
        shutil.copyfile(campaign / "gold" / f"orf_{key}.fasta", out / f"orf_{key}.fasta")

    p, s = gold["hosts"][primary], gold["hosts"][secondary]
    analysis = {
        "r3": {
            "unconstrained_violations": p["unconstrained_violations"],
            "binding_constraint": p["binding_constraint"],
            "non_binding_constraints": p["non_binding_constraints"],
        },
        "r4": {
            "unconstrained_violations": s["unconstrained_violations"],
            "binding_constraint": s["binding_constraint"],
            "non_binding_constraints": s["non_binding_constraints"],
            "flipped_constraints": gold["flipped_constraints"],
        },
    }
    (out / "analysis.json").write_text(json.dumps(analysis, indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
