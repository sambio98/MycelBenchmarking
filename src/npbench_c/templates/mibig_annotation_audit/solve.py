"""Reference solution for the annotation audit.

Usage: python -m npbench_c.templates.mibig_annotation_audit.solve <campaign> <submission>
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

    report = {k: gold[k] for k in (
        "corpus_entries", "entries_with_annotations", "total_function_annotations",
        "function_vocabulary", "backed_by_function", "unbacked_by_function",
        "function_evidence_method_totals", "admissible_by_function",
        "backed_minus_admissible", "focal_admissible_entries")}
    report["policy_variants"] = {
        k: {"admissible_by_function": v["admissible_by_function"], "total": v["total"]}
        for k, v in gold["policy_variants"].items()}
    (out / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
