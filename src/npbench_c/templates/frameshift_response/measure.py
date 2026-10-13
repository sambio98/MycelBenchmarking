"""Measure a frameshift-response submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies in `reference/frameshift_rules.json` -- the modes, the target
lengths, the offset grid, the noise subset and the verdicts -- or from this
template's own closed column lists.

Usage:
  python -m npbench_c.templates.frameshift_response.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.frameshift_response.core import Rules

REPORT = "frameshift_report.json"

TARGET_KEYS = ("start", "end", "strand", "length", "deleted_at")
SUMMARY_KEYS = ("targets", "shortest", "longest")
RESPONSE_KEYS = ("deleted_at", "target_start", "target_end", "target_strand",
                 "target_length", "verdict", "overlapping_new_calls",
                 "new_call_lengths", "bases_of_the_target_still_called",
                 "calls_before", "calls_after", "other_calls_lost",
                 "other_calls_gained")
SWEEP_KEYS = ("deleted_at", "bases_into_the_gene", "verdict",
              "new_call_lengths", "bases_of_the_target_still_called")
NOISE_KEYS = ("deleted_at", "verdict", "calls_explained_by_the_target",
              "other_calls_moved", "noise_to_signal")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _opt_float(value: object) -> float | None:
    """A ratio that is legitimately null in gold, where a target explains no call."""
    return None if value is None else _float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _str_list(value: object) -> list[str]:
    return [str(v) for v in value] if isinstance(value, list) else []


def _int_list(value: object) -> list[int] | None:
    """A list of lengths. None when the field is not a list at all, so a
    submission that omitted it is not silently equal to one reporting none."""
    if not isinstance(value, list):
        return None
    return [_int(v) for v in value]


def _bool(value: object) -> bool | None:
    return value if isinstance(value, bool) else None


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _block(value: object, keys, casts) -> dict:
    inner = value if isinstance(value, dict) else {}
    return {k: casts[k](inner.get(k)) for k in keys}


_RESPONSE_CASTS = {
    "deleted_at": _int, "target_start": _int, "target_end": _int,
    "target_strand": _str, "target_length": _int, "verdict": _str,
    "overlapping_new_calls": _int, "new_call_lengths": _int_list,
    "bases_of_the_target_still_called": _int, "calls_before": _int,
    "calls_after": _int, "other_calls_lost": _int, "other_calls_gained": _int,
}
_SWEEP_CASTS = {
    "deleted_at": _int, "bases_into_the_gene": _int, "verdict": _str,
    "new_call_lengths": _int_list,
    "bases_of_the_target_still_called": _int,
}
_NOISE_CASTS = {
    "deleted_at": _int, "verdict": _str,
    "calls_explained_by_the_target": _int, "other_calls_moved": _int,
    "noise_to_signal": _opt_float,
}
_TARGET_CASTS = {"start": _int, "end": _int, "strand": _str, "length": _int,
                 "deleted_at": _int}


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    labels = tuple(str(n) for n in rules.target_lengths)

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

    # ---- R1: the contig, the baseline run and the declared targets
    out["contig_length"] = _int(sub.get("contig_length"))
    out["contig_sha256"] = _str(sub.get("contig_sha256"))
    out["baseline_calls"] = _int(sub.get("baseline_calls"))
    table = sub.get("target_table")
    table = table if isinstance(table, dict) else {}
    out["target_table"] = {
        label: _block(table.get(label), TARGET_KEYS, _TARGET_CASTS)
        for label in labels}
    out["target_length_summary"] = _int_map(sub.get("target_length_summary"),
                                            SUMMARY_KEYS)
    out["environment"] = {
        "prodigal_version": _str(sub.get("prodigal_version")),
        "baseline_mode": _str(sub.get("baseline_mode")),
        "noise_mode": _str(sub.get("noise_mode")),
        "fixed_flags": _str_list(sub.get("fixed_flags")),
    }

    # ---- R2: the response to each planted frameshift
    responses = sub.get("response_table")
    responses = responses if isinstance(responses, dict) else {}
    out["response_table"] = {
        label: _block(responses.get(label), RESPONSE_KEYS, _RESPONSE_CASTS)
        for label in labels}

    # ---- R3: the verdicts
    out["verdict_census"] = _int_map(sub.get("verdict_census"), rules.verdicts)
    by_length = sub.get("verdict_by_length")
    by_length = by_length if isinstance(by_length, dict) else {}
    out["verdict_by_length"] = {label: _str(by_length.get(label))
                                for label in labels}
    out["verdict_is_monotone_in_length"] = _bool(
        sub.get("verdict_is_monotone_in_length"))

    # ---- R4: the offset sweep and the noise arm
    sweep = sub.get("offset_sweep")
    sweep = sweep if isinstance(sweep, dict) else {}
    out["offset_sweep"] = {
        f"{fraction:g}": _block(sweep.get(f"{fraction:g}"), SWEEP_KEYS,
                                _SWEEP_CASTS)
        for fraction in rules.offset_fractions}
    noise = sub.get("noise_arm")
    noise = noise if isinstance(noise, dict) else {}
    out["noise_arm"] = {
        str(length): _block(noise.get(str(length)), NOISE_KEYS, _NOISE_CASTS)
        for length in rules.noise_target_lengths}
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
