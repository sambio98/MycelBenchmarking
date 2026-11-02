"""Measure a kinetics submission. Deterministic, agent-independent.

Nothing here reads gold: the field names come from the template and the variant
names from the campaign's own task.yaml.

Usage:
  python -m npbench_c.templates.kinetics_consistency.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.kinetics_consistency.analyse import load_params
from npbench_c.templates.kinetics_consistency.core import FIELDS, Rules

REPORT = "kinetics_report.json"
BANDS = ("agree", "factor_1000", "factor_60", "near", "other")
COMPLETENESS = ("1", "2", "3")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _str_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: (str(block[k]) if isinstance(block.get(k), str) else None)
            for k in keys}


def _nested(value: object, *keys: str) -> object:
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    params = load_params(campaign)
    fields = tuple(sorted(FIELDS))
    variants = tuple(sorted(params.get("unit_error_variants") or {}))

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

    # ---- R1: the parse census
    for key in ("release", "ec_prefix"):
        value = sub.get(key)
        out[key] = str(value) if isinstance(value, str) else None
    out["distinct_ec_numbers"] = _int(sub.get("distinct_ec_numbers"))
    out["declared_units"] = _str_map(sub.get("declared_units"), fields)
    out["records_by_field"] = _int_map(sub.get("records_by_field"), fields)
    out["parsed_by_field"] = _int_map(sub.get("parsed_by_field"), fields)
    out["sentinel_and_unparsed"] = {
        "sentinel": _int_map(sub.get("sentinel_by_field"), fields),
        "unparsed": _int_map(sub.get("unparsed_by_field"), fields),
    }

    # ---- R2: the triple inventory
    out["triple_keys_total"] = _int(sub.get("triple_keys_total"))
    out["completeness_counts"] = _int_map(sub.get("completeness_counts"),
                                          COMPLETENESS)
    for key in ("complete_triples", "ec_with_complete_triples",
                "distinct_substrates_in_complete_triples"):
        out[key] = _int(sub.get(key))

    # ---- R3: the consistency verdict
    out["checked"] = _int(sub.get("checked"))
    out["excluded_non_positive"] = _int(sub.get("excluded_non_positive"))
    # Zero-valued bands are kept, unlike a free census: the band vocabulary is
    # closed and declared, so an omitted band is a missing answer rather than a
    # presentation choice.
    out["band_counts"] = _int_map(sub.get("band_counts"),
                                  tuple(sorted({b.name for b in rules.bands}
                                               | {"other"})))
    out["agreement_rate"] = _float(sub.get("agreement_rate"))
    out["unit_mixup_total"] = _int(sub.get("unit_mixup_total"))

    # ---- R4: the unit-error counterfactual
    declared = sub.get("unit_error_variants")
    declared = declared if isinstance(declared, dict) else {}
    rates, agree, factors = {}, {}, {}
    for name in variants:
        block = declared.get(name)
        block = block if isinstance(block, dict) else {}
        rates[name] = _float(block.get("agreement_rate"))
        agree[name] = _int(block.get("agree"))
        factors[name] = {"km_factor": _float(block.get("km_factor")),
                         "kcat_factor": _float(block.get("kcat_factor"))}
    out["unit_error_rates"] = rates
    out["unit_error_agree"] = agree
    out["unit_error_factors"] = factors
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
