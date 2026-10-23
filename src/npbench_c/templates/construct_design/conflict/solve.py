"""Reference solution for the conflict ladder.

Must score exactly 1.0 through the real grader or the campaign does not ship.
Note what it does NOT do: it supplies no design, because none exists.

Usage:
  python -m npbench_c.templates.construct_design.conflict.solve <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: solve <campaign_dir> <submission_dir>", file=sys.stderr)
        return 2
    campaign, out = pathlib.Path(argv[0]), pathlib.Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((campaign / "gold" / "gold.json").read_text())

    report = {
        "orf_length_nt": gold["orf_length_nt"],
        "constraints_assessed": gold["constraint_vocabulary"],
        "admissible_codon_counts": gold["admissible_codon_counts"],
        "unrestricted_gc_min": gold["unrestricted_gc_min"],
        "unrestricted_gc_max": gold["unrestricted_gc_max"],
        "restricted_gc_min": gold["restricted_gc_min"],
        "restricted_gc_max": gold["restricted_gc_max"],
        "verdict": gold["verdict"],
        "conflicting_constraints": gold["conflicting_constraints"],
        "witness": {
            "required_gc_min": gold["witness"]["required_gc_min"],
            "max_achievable_gc_under_forbidden_codons":
                gold["witness"]["max_achievable_gc_under_forbidden_codons"],
        },
        "design_orf": None,
        "relaxations": {
            c: {"max_achievable_gc": v["max_achievable_gc"],
                "gc_requirement_reachable": v["gc_requirement_reachable"]}
            for c, v in gold["single_constraint_relaxations"].items()
        },
    }
    (out / "feasibility_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
