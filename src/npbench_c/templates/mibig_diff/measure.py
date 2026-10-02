"""Measure a release-diff submission. Deterministic, agent-independent.

Usage:
  python -m npbench_c.templates.mibig_diff.measure <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.mibig_diff.core import ReleaseDiff


def _int_or_none(v: object) -> int | None:
    return int(v) if isinstance(v, int) and not isinstance(v, bool) else None


def _int_map(value: object) -> dict:
    if not isinstance(value, dict):
        return {}
    # Zero-valued classes are dropped so a submission that lists every verdict
    # with zeros compares equal to gold, which omits them. Presentation is not
    # what is being graded.
    out = {}
    for k, v in value.items():
        n = _int_or_none(v)
        if n:
            out[str(k)] = n
    return dict(sorted(out.items()))


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    diff = ReleaseDiff.load(campaign)
    gold = json.loads((pathlib.Path(campaign) / "gold" / "gold.json").read_text())
    vocab = set(diff.verdict_vocabulary)

    path = submission_dir / "diff_report.json"
    out: dict = {"present": {"report": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["report"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["report"] = 0
        return out

    for key in ("source_entries", "target_entries", "only_in_source",
                "only_in_target", "shared_entries"):
        out[key] = _int_or_none(sub.get(key))

    fields = sub.get("fields")
    out["fields"] = sorted(str(f) for f in fields) if isinstance(fields, list) else []
    out["shared_status_in_target"] = _int_map(sub.get("shared_status_in_target"))

    raw = sub.get("raw_difference_counts")
    raw = raw if isinstance(raw, dict) else {}
    out["raw_difference_counts"] = {
        f: _int_or_none(raw.get(f)) for f in sorted(diff.fields)
    }

    classified = sub.get("classified_verdict_counts")
    classified = classified if isinstance(classified, dict) else {}
    out["classified_verdict_counts"] = {
        f: _int_map(classified.get(f)) for f in sorted(diff.fields)
    }

    terms = sub.get("retired_class_terms")
    out["retired_class_terms"] = sorted(str(t) for t in terms) if isinstance(terms, list) else []

    perturbations = sub.get("perturbations")
    perturbations = perturbations if isinstance(perturbations, dict) else {}
    out["perturbations"] = {}
    for name in sorted(gold["perturbations"]):
        block = perturbations.get(name) or {}
        block = block if isinstance(block, dict) else {}
        verdict = block.get("verdict")
        out["perturbations"][name] = {
            "verdict": verdict if verdict in vocab else None}
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>", file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
