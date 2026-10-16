"""Measure a remote-homology submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies in `reference/homology_rules.json` -- the bands, the verdicts, the
shuffle windows, the declared query count -- or from this template's own closed
column lists. The query accessions are the one thing the measurer cannot declare
in advance, because they follow from the declared selection rule over the shipped
corpus; it recomputes them from the shipped label table rather than reading them
out of the submission, so a submission that invented its own query set scores
zero instead of being read on its own terms.

Usage:
  python -m npbench_c.templates.remote_homology.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.remote_homology.analyse import load_params
from npbench_c.templates.remote_homology.core import Corpus

REPORT = "homology_report.json"

QUERY_KEYS = ("accession", "length", "family_members", "family_name")
RECOVERY_KEYS = ("hits", "same_family_asserted", "same_family_detected",
                 "same_family_missed", "recovery_rate",
                 "superfamily_only_detected")
LABELS = ("same_family", "superfamily_only", "unclassified_pair")
LADDER_KEYS = ("hits", "lost", "lost_by_band", "loss_rate_by_band",
               "band_order_held")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _opt_float(value: object) -> float | None:
    """A rate that is legitimately null in gold, where the denominator is zero."""
    return None if value is None else _float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _str_list(value: object) -> list[str]:
    return [str(v) for v in value] if isinstance(value, list) else []


def _bool(value: object) -> bool | None:
    return value if isinstance(value, bool) else None


def _block(value: object, keys, casts) -> dict:
    inner = value if isinstance(value, dict) else {}
    return {k: casts[k](inner.get(k)) for k in keys}


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _opt_float_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _opt_float(block.get(k)) for k in keys}


_QUERY_CASTS = {"accession": _str, "length": _int, "family_members": _int,
                "family_name": _str}
_RECOVERY_CASTS = {"hits": _int, "same_family_asserted": _int,
                   "same_family_detected": _int, "same_family_missed": _int,
                   "recovery_rate": _opt_float,
                   "superfamily_only_detected": _int}


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["corpus_input"], params["labels_input"])
    rules = corpus.rules
    queries = corpus.queries()
    entries = sorted(queries)
    accessions = sorted(queries.values())
    bands = rules.band_labels
    windows = [str(w) for w in rules.shuffle_windows]

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

    # ---- R1: the corpus as parsed, the declared queries, and the search
    out["corpus_sequences"] = _int(sub.get("corpus_sequences"))
    out["corpus_sha256"] = _str(sub.get("corpus_sha256"))
    table = sub.get("query_table")
    table = table if isinstance(table, dict) else {}
    out["query_table"] = {
        entry: _block(table.get(entry), QUERY_KEYS, _QUERY_CASTS)
        for entry in entries}
    out["baseline_hits"] = _int_map(sub.get("baseline_hits"), accessions)
    sizes = sub.get("baseline_band_sizes")
    sizes = sizes if isinstance(sizes, dict) else {}
    out["baseline_band_sizes"] = {
        accession: _int_map(sizes.get(accession), bands)
        for accession in accessions}
    out["environment"] = {
        "mmseqs_version": _str(sub.get("mmseqs_version")),
        "sensitivity": _float(sub.get("sensitivity")),
        "evalue": _float(sub.get("evalue")),
        "max_seqs": _int(sub.get("max_seqs")),
        "threads": _int(sub.get("threads")),
        "format_fields": _str_list(sub.get("format_fields")),
    }

    # ---- R2: what the search found of what the resource asserts
    recovery = sub.get("recovery")
    recovery = recovery if isinstance(recovery, dict) else {}
    out["recovery"] = {
        accession: _block(recovery.get(accession), RECOVERY_KEYS, _RECOVERY_CASTS)
        for accession in accessions}
    census = sub.get("band_census")
    census = census if isinstance(census, dict) else {}
    out["band_census"] = {
        accession: {
            band: _int_map((census.get(accession) or {}).get(band)
                           if isinstance(census.get(accession), dict) else None,
                           LABELS)
            for band in bands}
        for accession in accessions}

    # ---- R3: the declared verdicts
    verdicts = sub.get("band_verdicts")
    verdicts = verdicts if isinstance(verdicts, dict) else {}
    out["band_verdicts"] = {
        accession: {
            band: _str((verdicts.get(accession) or {}).get(band)
                       if isinstance(verdicts.get(accession), dict) else None)
            for band in bands}
        for accession in accessions}
    out["verdict_census"] = _int_map(sub.get("verdict_census"), rules.verdicts)

    # ---- R4: the shuffle ladder
    ladder = sub.get("shuffle_ladder")
    ladder = ladder if isinstance(ladder, dict) else {}
    out["shuffle_ladder"] = {}
    for window in windows:
        rows = ladder.get(window)
        rows = rows if isinstance(rows, dict) else {}
        out["shuffle_ladder"][window] = {}
        for accession in accessions:
            row = rows.get(accession)
            row = row if isinstance(row, dict) else {}
            out["shuffle_ladder"][window][accession] = {
                "hits": _int(row.get("hits")),
                "lost": _int(row.get("lost")),
                "lost_by_band": _int_map(row.get("lost_by_band"), bands),
                "loss_rate_by_band": _opt_float_map(
                    row.get("loss_rate_by_band"), bands),
                "band_order_held": _bool(row.get("band_order_held")),
            }
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
