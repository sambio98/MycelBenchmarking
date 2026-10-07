"""Measure a RiPP-precursor submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies -- the localisation verdicts, the cleavage conventions, the
reconciliation verdicts, the model panel, the boundary relations, the deletion
outcomes and the ingestion policies -- all of which are published in
`reference/precursor_rules.json`.

The row-shaped answers (the localisation table, the cleavage table, the domain
table, the deletion table) are projected to their declared columns and left in
the submission's own order, because gold's order is declared too and a
submission that sorts differently has got the order wrong.

Usage:
  python -m npbench_c.templates.ripp_precursor.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.ripp_precursor.core import Rules

REPORT = "precursor_report.json"

LOCALISATION_COLUMNS = ("bgc", "gene", "cores", "core_length", "occurrences",
                        "core_start", "leader_length", "follower_length",
                        "verdict")
CLEAVAGE_COLUMNS = ("bgc", "gene", "declared_from", "declared_to",
                    "implied_core_start", "derived_core_start", "verdict")
DOMAIN_COLUMNS = ("bgc", "gene", "model", "env_from", "env_to", "relation",
                  "boundary_offset")
DELETION_COLUMNS = ("bgc", "gene", "before", "after", "outcome")
SUMMARY_KEYS = ("count", "min", "median", "max")
OFFSET_KEYS = ("count", "min", "median", "max", "exactly_at_the_boundary")
SELECTION_KEYS = ("entries_in_class", "entries_with_a_named_precursor",
                  "entries_excluded_by_status",
                  "entries_absent_from_reference_set",
                  "entries_excluded_by_both",
                  "reference_records_shipped", "records_named",
                  "records_with_an_unresolved_gene", "corpus_entries",
                  "corpus_records")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _opt_int(value: object) -> int | None:
    """An integer column that is legitimately null in gold.

    `core_start` and `boundary_offset` are null wherever the campaign cannot
    derive them, so None has to survive the projection rather than being coerced
    to a number a submission never claimed.
    """
    return None if value is None else _int(value)


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _opt_float(value: object) -> float | None:
    return None if value is None else _float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _str_list(value: object) -> list[str]:
    return [str(v) for v in value] if isinstance(value, list) else []


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _rows(value: object, columns, casts) -> list[dict]:
    """Project a submission's rows to the declared columns, in its own order."""
    rows = value if isinstance(value, list) else []
    out = []
    for row in rows:
        block = row if isinstance(row, dict) else {}
        out.append({c: casts[c](block.get(c)) for c in columns})
    return out


_LOCALISATION_CASTS = {
    "bgc": _str, "gene": _str, "cores": _int, "core_length": _int,
    "occurrences": _opt_int, "core_start": _opt_int, "leader_length": _opt_int,
    "follower_length": _opt_int, "verdict": _str,
}
_CLEAVAGE_CASTS = {
    "bgc": _str, "gene": _str, "declared_from": _int, "declared_to": _int,
    "implied_core_start": _int, "derived_core_start": _opt_int, "verdict": _str,
}
_DOMAIN_CASTS = {
    "bgc": _str, "gene": _str, "model": _str, "env_from": _int, "env_to": _int,
    "relation": _str, "boundary_offset": _opt_int,
}
_DELETION_CASTS = {
    "bgc": _str, "gene": _str, "before": _str_list, "after": _str_list,
    "outcome": _str,
}


def _summary(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    out: dict = {}
    for key in keys:
        if key == "median":
            out[key] = _opt_float(block.get(key))
        elif key in ("min", "max"):
            out[key] = _opt_int(block.get(key))
        else:
            out[key] = _int(block.get(key))
    return out


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)

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

    # ---- R1: the selection, the shapes, and what produced the domain calls
    out["selection_census"] = _int_map(sub.get("selection_census"),
                                      SELECTION_KEYS)
    out["core_sequence_shape_census"] = _drop_zeros(
        _int_map(sub.get("core_sequence_shape_census"), rules.shape_order))
    out["core_sequence_case_census"] = _drop_zeros(
        _int_map(sub.get("core_sequence_case_census"), rules.case_order))
    out["records_declaring_a_cleavage_location"] = _int(
        sub.get("records_declaring_a_cleavage_location"))
    out["ripp_type_census"] = _free_int_map(sub.get("ripp_type_census"))
    out["status_census"] = _free_int_map(sub.get("status_census"))
    out["precursor_length_summary"] = _summary(
        sub.get("precursor_length_summary"), SUMMARY_KEYS)
    out["environment"] = {
        "hmmer_version": _str(sub.get("hmmer_version")),
        "pfam_version": _str(sub.get("pfam_version")),
        "pfam_sha256": _str(sub.get("pfam_sha256")),
    }

    # ---- R2: where each core sits
    out["localisation_table"] = _rows(sub.get("localisation_table"),
                                      LOCALISATION_COLUMNS, _LOCALISATION_CASTS)
    out["localisation_census"] = _int_map(sub.get("localisation_census"),
                                          rules.verdict_order)
    out["leader_length_summary"] = _summary(sub.get("leader_length_summary"),
                                            SUMMARY_KEYS)

    # ---- R3: the coordinate convention, and the domain evidence
    tally = sub.get("convention_tally")
    tally = tally if isinstance(tally, dict) else {}
    out["convention_tally"] = {
        convention: _int_map(tally.get(convention), rules.reconciliation_verdicts)
        for convention in rules.conventions}
    out["elected_convention"] = _str(sub.get("elected_convention"))
    out["cleavage_table"] = _rows(sub.get("cleavage_table"), CLEAVAGE_COLUMNS,
                                  _CLEAVAGE_CASTS)
    out["domain_table"] = _rows(sub.get("domain_table"), DOMAIN_COLUMNS,
                                _DOMAIN_CASTS)
    out["panel_census"] = _int_map(sub.get("panel_census"), rules.panel)
    out["records_with_a_panel_hit"] = _int(sub.get("records_with_a_panel_hit"))
    out["boundary_relation_census"] = _int_map(
        sub.get("boundary_relation_census"), rules.relation_order)
    out["boundary_offset_summary"] = _summary(
        sub.get("boundary_offset_summary"), OFFSET_KEYS)

    # ---- R4: the two counterfactual axes
    sweep = sub.get("policy_sweep")
    sweep = sweep if isinstance(sweep, dict) else {}
    out["policy_sweep"] = {}
    for policy in sorted(rules.policies):
        block = sweep.get(policy)
        block = block if isinstance(block, dict) else {}
        out["policy_sweep"][policy] = {
            "verdict_census": _int_map(block.get("verdict_census"),
                                       rules.verdict_order),
            "located": _int(block.get("located")),
            "cores_declared": _int(block.get("cores_declared")),
        }
    out["deletion_table"] = _rows(sub.get("deletion_table"), DELETION_COLUMNS,
                                  _DELETION_CASTS)
    out["deletion_census"] = _int_map(sub.get("deletion_census"),
                                      rules.deletion_outcomes)
    return out


def _drop_zeros(block: dict) -> dict:
    """Keys whose count is zero are omitted, as gold omits them.

    The shape and case censuses are over a closed set of possibilities but gold
    lists only what occurs, so a submission that spells out the absent shapes
    should still equal it.
    """
    return {k: v for k, v in block.items() if v}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped.

    The RiPP-type and status censuses are open-ended: their keys are whichever
    terms the selected entries happen to carry, so they cannot be read against a
    declared key list.
    """
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


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
