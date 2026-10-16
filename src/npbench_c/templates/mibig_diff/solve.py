"""Reference solution for the release-diff ladder.

Usage: python -m npbench_c.templates.mibig_diff.solve <campaign> <submission>
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
        "source_entries": gold["source_entries"],
        "target_entries": gold["target_entries"],
        "fields": gold["fields"],
        "only_in_source": gold["only_in_source"],
        "only_in_target": gold["only_in_target"],
        "shared_entries": gold["shared_entries"],
        "shared_status_in_target": gold["shared_status_in_target"],
        "raw_difference_counts": gold["raw_difference_counts"],
        "classified_verdict_counts": gold["classified_verdict_counts"],
        "retired_class_terms": gold["retired_class_terms"],
        "perturbations": {k: {"verdict": v["verdict"]}
                          for k, v in gold["perturbations"].items()},
    }
    (out / "diff_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
