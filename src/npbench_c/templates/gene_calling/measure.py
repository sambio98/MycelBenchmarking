"""Measure a gene-calling submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies in `reference/calling_rules.json` -- the modes, the fragment grid,
the loss verdicts, the call columns and the partial flags -- or from this
template's own closed column lists.

Usage:
  python -m npbench_c.templates.gene_calling.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.gene_calling.core import Rules

REPORT = "calling_report.json"

CENSUS_KEYS = ("calls", "on_the_forward_strand", "on_the_reverse_strand",
               "shortest", "median_length", "longest", "coding_bases",
               "coding_fraction")
SPAN_KEYS = ("first_call_start", "last_call_end",
             "bases_before_the_first_call", "bases_after_the_last_call")
SUMMARY_KEYS = ("baseline_calls", "fragmented_calls", "retained", "lost",
                "gained", "retained_but_truncated")
GRID_CELL_KEYS = ("baseline_calls", "fragmented_calls", "retained", "lost",
                  "gained", "retained_but_truncated")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _str_list(value: object) -> list[str]:
    return [str(v) for v in value] if isinstance(value, list) else []


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped: gold omits them.

    The start-type census is keyed by whichever start codons the contig turns out
    to use, so it cannot be read against a declared key list.
    """
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def _census(value: object) -> dict:
    block = value if isinstance(value, dict) else {}
    out = {}
    for key in CENSUS_KEYS:
        if key in ("median_length", "coding_fraction"):
            out[key] = _float(block.get(key))
        else:
            out[key] = _int(block.get(key))
    return out


_CALL_CASTS = {"start": _int, "end": _int, "strand": _str, "partial": _str,
               "start_type": _str, "fragment": _str, "verdict": _str}


def _rows(value: object, columns) -> list[dict]:
    """Project a submission's rows to the declared columns, in its own order."""
    rows = value if isinstance(value, list) else []
    out = []
    for row in rows:
        block = row if isinstance(row, dict) else {}
        out.append({c: _CALL_CASTS[c](block.get(c)) for c in columns})
    return out


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    columns = tuple(rules.call_columns)

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

    # ---- R1: the contig and one run over it
    out["accession"] = _str(sub.get("accession"))
    out["contig_length"] = _int(sub.get("contig_length"))
    out["contig_sha256"] = _str(sub.get("contig_sha256"))
    out["gc_fraction"] = _float(sub.get("gc_fraction"))
    out["call_census"] = _census(sub.get("call_census"))
    out["partial_census"] = _int_map(sub.get("partial_census"),
                                     rules.partial_flags)
    out["start_type_census"] = _free_int_map(sub.get("start_type_census"))
    out["breakpoint_census"] = _int_map(
        sub.get("breakpoint_census"),
        tuple(str(n) for n in rules.fragment_counts))
    out["environment"] = {
        "prodigal_version": _str(sub.get("prodigal_version")),
        "baseline_mode": _str(sub.get("baseline_mode")),
        "fixed_flags": _str_list(sub.get("fixed_flags")),
    }

    # ---- R2: the calls
    out["call_table"] = _rows(sub.get("call_table"), columns)
    out["call_span"] = _int_map(sub.get("call_span"), SPAN_KEYS)

    # ---- R3: the fragmented contig against the whole one
    # Assembled from the six reported counts rather than read as a composite,
    # the same way `environment` is: what a submission reports are the counts,
    # and the grouping into one claim is the campaign's, not something the
    # output contract should make an agent restate.
    out["fragmentation_summary"] = {
        key: _int(sub.get(f"fragmented_{key}")) for key in SUMMARY_KEYS}
    out["fragmented_loss_verdicts"] = _int_map(
        sub.get("fragmented_loss_verdicts"), rules.loss_verdicts)
    out["fragmented_lost_table"] = _rows(sub.get("fragmented_lost_table"),
                                         (*columns, "verdict"))
    out["fragmented_partial_census"] = _int_map(
        sub.get("fragmented_partial_census"), rules.partial_flags)
    out["fragmented_gained_partial_census"] = _int_map(
        sub.get("fragmented_gained_partial_census"), rules.partial_flags)
    out["fragmented_call_census"] = _census(sub.get("fragmented_call_census"))

    # ---- R4: the declared grid and the mode comparison
    grid = sub.get("fragmentation_grid")
    grid = grid if isinstance(grid, dict) else {}
    out["fragmentation_grid"] = {}
    for mode in rules.modes:
        for count in rules.fragment_counts:
            label = f"{mode}|{count}"
            cell = grid.get(label)
            cell = cell if isinstance(cell, dict) else {}
            out["fragmentation_grid"][label] = {
                **_int_map(cell, GRID_CELL_KEYS),
                "loss_verdicts": _int_map(cell.get("loss_verdicts"),
                                          rules.loss_verdicts),
            }
    difference = sub.get("mode_difference")
    difference = difference if isinstance(difference, dict) else {}
    out["mode_difference"] = {}
    reference = rules.baseline_mode
    for mode in rules.modes:
        if mode == reference:
            continue
        label = f"{reference}_vs_{mode}"
        block = difference.get(label)
        block = block if isinstance(block, dict) else {}
        out["mode_difference"][label] = _int_map(
            block, ("shared", f"only_in_{reference}", f"only_in_{mode}"))
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
