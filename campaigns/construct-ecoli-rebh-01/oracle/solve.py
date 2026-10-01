"""Reference solution: produces a complete, correct submission.

Must score exactly 1.0 through the real grader or the campaign does not ship
(spec 2.7). Used by the readiness gate, never shipped to agents.
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

CAMPAIGN = pathlib.Path(__file__).resolve().parent.parent


def main() -> int:
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else CAMPAIGN / "oracle_submission"
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())

    for key in ("ecoli", "streptomyces"):
        shutil.copyfile(CAMPAIGN / "gold" / f"orf_{key}.fasta", out / f"orf_{key}.fasta")

    e, s = gold["hosts"]["ecoli"], gold["hosts"]["streptomyces"]
    analysis = {
        "r3": {
            "unconstrained_violations": e["unconstrained_violations"],
            "binding_constraint": e["binding_constraint"],
            "non_binding_constraints": e["non_binding_constraints"],
        },
        "r4": {
            "unconstrained_violations": s["unconstrained_violations"],
            "binding_constraint": s["binding_constraint"],
            "non_binding_constraints": s["non_binding_constraints"],
            "flipped_constraints": gold["flipped_constraints"],
        },
    }
    (out / "analysis.json").write_text(json.dumps(analysis, sort_keys=True, indent=2) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
