"""Measure an annotation-transfer submission. Deterministic, agent-independent.

Every vocabulary comes from the campaign's own declared rules -- the verdict enum,
the identity bands, the graded depths, the removal variants -- so a submission is
never scored against something it was not told.

Usage:
  python -m npbench_c.templates.annotation_transfer.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.annotation_transfer.core import Rules

REPORT = "transfer_report.json"


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped: gold omits them."""
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def _agreement_block(value: object, keys) -> dict:
    """decidable / agree / rate per declared key. A rate is a float or None, and
    a submission that reports a rate where gold has None has said something about
    an empty band, which is not the same as saying nothing."""
    block = value if isinstance(value, dict) else {}
    out = {}
    for key in keys:
        inner = block.get(key)
        inner = inner if isinstance(inner, dict) else {}
        out[key] = {"decidable": _int(inner.get("decidable")),
                    "agree": _int(inner.get("agree")),
                    "rate": _float(inner.get("rate"))}
    return out


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    bands = tuple(rules.band_label(lo, hi) for lo, hi in rules.identity_bands)
    depths = tuple(str(d) for d in rules.graded_depths)
    verdicts = rules.verdict_order
    variants = tuple(sorted(rules.removal_variants))

    path = pathlib.Path(submission_dir) / REPORT
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

    # ---- R1: the set, and what searched it
    for key in ("release", "ec_prefix", "tie_break", "diamond_version"):
        out[key] = _str(sub.get(key))
    for key in ("entries_total", "distinct_ec_codes", "partial_ec_mentions",
                "entries_with_multiple_ec"):
        out[key] = _int(sub.get(key))
    out["ec_codes_per_entry"] = _free_int_map(sub.get("ec_codes_per_entry"))

    # ---- R2: the search
    out["queries_with_a_hit"] = _int(sub.get("queries_with_a_hit"))
    out["queries_without_a_hit"] = _int(sub.get("queries_without_a_hit"))
    out["tied_top_queries"] = _int(sub.get("tied_top_queries"))
    out["top_hit_identity_census"] = _int_map(sub.get("top_hit_identity_census"),
                                              bands)

    # ---- R3: the verdict
    out["agreement_by_depth"] = _agreement_block(sub.get("agreement_by_depth"),
                                                 depths)
    out["agreement_by_identity"] = _agreement_block(
        sub.get("agreement_by_identity"), bands)
    out["verdict_census"] = _int_map(sub.get("verdict_census"), verdicts)
    errors = sub.get("transfer_errors_at_subsubclass")
    out["transfer_errors_at_subsubclass"] = (
        sorted({str(a) for a in errors}) if isinstance(errors, list) else [])

    # ---- R4: the removal counterfactual
    declared = sub.get("removal_census")
    declared = declared if isinstance(declared, dict) else {}
    out["removal_census"] = {name: _int_map(declared.get(name), verdicts)
                             for name in variants}
    out["removal_changed"] = _int_map(sub.get("removal_changed"), variants)
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>",
              file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
