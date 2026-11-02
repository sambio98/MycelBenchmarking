"""Measure a BGC-detection submission. Deterministic, agent-independent.

Every vocabulary comes from the campaign's own declared rules -- the verdict enum,
the deletion variant names, the curated locus list -- so a submission is never
scored against something it was not told.

Usage:
  python -m npbench_c.templates.bgc_detection.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.bgc_detection.core import Rules

REPORT = "detection_report.json"


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _bool(value: object) -> bool | None:
    return value if isinstance(value, bool) else None


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped: gold omits them."""
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def _str_list(value: object) -> list:
    return sorted({str(v) for v in value}) if isinstance(value, list) else []


def _region_table(value: object) -> list:
    """Normalised region rows, sorted by coordinate.

    The row carries both conventions, because the 0-based half-open form and the
    1-based inclusive form are where an off-by-one hides: a submission that
    reports one correctly and the other shifted has not understood the
    convention, and grading only one would let that pass.
    """
    rows = value if isinstance(value, list) else []
    out = []
    for row in rows:
        row = row if isinstance(row, dict) else {}
        out.append({
            "number": _int(row.get("number")),
            "start": _int(row.get("start")),
            "end": _int(row.get("end")),
            "start_one_based": _int(row.get("start_one_based")),
            "end_one_based": _int(row.get("end_one_based")),
            "length": _int(row.get("length")),
            "products": _str_list(row.get("products")),
            "contig_edge": _bool(row.get("contig_edge")),
        })
    return sorted(out, key=lambda r: (r["start"] is None, r["start"],
                                      r["end"] is None, r["end"]))


def _per_locus(value: object, verdicts) -> list:
    rows = value if isinstance(value, list) else []
    out = []
    for row in rows:
        row = row if isinstance(row, dict) else {}
        verdict = _str(row.get("verdict"))
        out.append({
            "bgc": _str(row.get("bgc")),
            "verdict": verdict if verdict in verdicts else None,
            "region": _int(row.get("region")),
            "coverage": _float(row.get("coverage")),
            "jaccard": _float(row.get("jaccard")),
        })
    return sorted(out, key=lambda r: (r["bgc"] is None, r["bgc"]))


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    verdicts = rules.verdict_order
    variants = tuple(sorted(rules.deletion_variants))

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

    # ---- R1: the inventory, and which rules produced it
    out["accession"] = _str(sub.get("accession"))
    out["record_length"] = _int(sub.get("record_length"))
    out["coordinate_convention"] = _str(sub.get("coordinate_convention"))
    out["search_flags"] = (list(sub["search_flags"])
                           if isinstance(sub.get("search_flags"), list) else [])
    out["regions"] = _int(sub.get("regions"))
    out["distinct_products"] = _int(sub.get("distinct_products"))
    out["contig_edge_regions"] = _int(sub.get("contig_edge_regions"))
    out["environment"] = {
        "antismash_version": _str(sub.get("antismash_version")),
        "detection_rule_fingerprint": _str(sub.get("detection_rule_fingerprint")),
    }

    # ---- R2: the boundaries
    out["region_table"] = _region_table(sub.get("region_table"))
    span = sub.get("region_span")
    span = span if isinstance(span, dict) else {}
    out["region_span"] = {"min": _int(span.get("min")),
                          "median": _float(span.get("median")),
                          "max": _int(span.get("max"))}
    out["product_census"] = _free_int_map(sub.get("product_census"))
    out["region_bases"] = _int(sub.get("region_bases"))
    out["fraction_of_record"] = _float(sub.get("fraction_of_record"))

    # ---- R3: the reconciliation against the curated boundaries
    out["verdict_census"] = _int_map(sub.get("verdict_census"), verdicts)
    out["per_locus"] = _per_locus(sub.get("per_locus"), verdicts)
    out["loci_total"] = _int(sub.get("loci_total"))
    out["regions_without_a_locus"] = _int(sub.get("regions_without_a_locus"))
    shared = sub.get("regions_holding_several_loci")
    out["regions_holding_several_loci"] = (
        sorted({n for v in shared if (n := _int(v)) is not None})
        if isinstance(shared, list) else [])
    out["jaccard_summary"] = {
        "median": _float(sub.get("jaccard_median")),
        "min": _float(sub.get("jaccard_min")),
        "max": _float(sub.get("jaccard_max")),
    }

    # ---- R4: the deletion counterfactual
    declared = sub.get("deletion_outcomes")
    declared = declared if isinstance(declared, dict) else {}
    summary = {}
    for name in variants:
        block = declared.get(name)
        block = block if isinstance(block, dict) else {}
        summary[name] = {
            "regions": _int(block.get("regions")),
            "regions_lost": _int(block.get("regions_lost")),
            "regions_gained": _int(block.get("regions_gained")),
            "regions_unchanged": _int(block.get("regions_unchanged")),
            "lost_products": _str_list(block.get("lost_products")),
            "gained_products": _str_list(block.get("gained_products")),
        }
    out["deletion_summary"] = summary
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
