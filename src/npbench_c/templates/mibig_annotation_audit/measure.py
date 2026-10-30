"""Measure an annotation-audit submission. Deterministic, agent-independent.

Usage:
  python -m npbench_c.templates.mibig_annotation_audit.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.mibig_annotation_audit.core import AnnotationAudit


def _int_or_none(v: object) -> int | None:
    return int(v) if isinstance(v, int) and not isinstance(v, bool) else None


def _int_map(value: object, keys: tuple[str, ...] | None = None) -> dict:
    if not isinstance(value, dict):
        return {}
    if keys is not None:
        return {k: _int_or_none(value.get(k)) for k in keys}
    # Zero-valued entries are dropped so a submission listing every method with
    # zeros compares equal to gold, which omits none. Presentation is not graded.
    return {str(k): n for k, v in sorted(value.items())
            if (n := _int_or_none(v)) is not None}


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    audit = AnnotationAudit.load(campaign)
    functions = tuple(sorted(audit.function_vocabulary))
    gold = json.loads((pathlib.Path(campaign) / "gold" / "gold.json").read_text())

    path = submission_dir / "audit_report.json"
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

    for key in ("corpus_entries", "entries_with_annotations",
                "total_function_annotations"):
        out[key] = _int_or_none(sub.get(key))

    vocab = sub.get("function_vocabulary")
    out["function_vocabulary"] = sorted(str(f) for f in vocab) if isinstance(vocab, list) else []

    for key in ("backed_by_function", "unbacked_by_function",
                "admissible_by_function", "backed_minus_admissible"):
        out[key] = _int_map(sub.get(key), functions)
    out["function_evidence_method_totals"] = _int_map(
        sub.get("function_evidence_method_totals"))

    focal = sub.get("focal_admissible_entries")
    out["focal_admissible_entries"] = sorted(str(a) for a in focal) if isinstance(focal, list) else []

    variants = sub.get("policy_variants")
    variants = variants if isinstance(variants, dict) else {}
    out["policy_variants"] = {}
    for name in sorted(gold["policy_variants"]):
        block = variants.get(name) or {}
        block = block if isinstance(block, dict) else {}
        out["policy_variants"][name] = {
            "admissible_by_function": _int_map(block.get("admissible_by_function"),
                                               functions),
            "total": _int_or_none(block.get("total")),
        }
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
