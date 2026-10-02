"""Measure a chemical-space submission. Deterministic, agent-independent.

Groups the submission's counts into the same views the grading spec compares, so
the agent's contract stays in the readable nested shape while the spec keeps one
component per claim. Nothing here reads gold: the class names, identity keys and
variant names all come from the campaign's own declared vocabularies.

Usage:
  python -m npbench_c.templates.chemical_space.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.chemical_space.analyse import load_params
from npbench_c.templates.chemical_space.core import KEYS, Rules

REPORT = "chemical_space_report.json"

SHARED_KEYS = ("all", "informative")
CROSS_KEYS = ("all", "informative", "naive_all", "naive_informative")
VARIANT_SHARING_KEYS = ("shared_informative", "cross_class_informative",
                        "naive_cross_class_informative")
GATE_COUNT_KEYS = ("entries", "compound_records", "distinct_compound_keys",
                   "distinct_scaffold_keys", "informative_records",
                   "informative_scaffold_keys")
CORPUS_COUNT_KEYS = ("corpus_entries", "active_entries", "compound_records",
                     "entries_with_compounds")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped.

    Gold omits an empty bucket, so a submission that lists every ring count
    including the zero-valued ones still compares equal: presentation is not
    graded.
    """
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


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
    classes = tuple(sorted(rules.class_vocabulary))
    variants = tuple(sorted(params.get("threshold_variants") or {}))

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

    # ---- R1: the corpus census
    release = sub.get("corpus_release")
    out["corpus_release"] = str(release) if release is not None else None
    out["corpus_counts"] = _int_map(sub, CORPUS_COUNT_KEYS)
    out["distinct_keys"] = {
        "compound": _int(sub.get("distinct_compound_keys")),
        "connectivity": _int(sub.get("distinct_connectivity_keys")),
        "scaffold": _int(sub.get("distinct_scaffold_keys")),
    }
    out["acyclic_compounds"] = _int(sub.get("acyclic_compounds"))
    out["duplicate_records"] = _int(sub.get("duplicate_records"))
    vocabulary = sub.get("class_vocabulary")
    out["class_vocabulary"] = (sorted({str(c) for c in vocabulary})
                               if isinstance(vocabulary, list) else [])
    out["entries_per_class"] = _int_map(sub.get("entries_per_class"), classes)
    out["multi_class_entries"] = _int(sub.get("multi_class_entries"))

    # ---- R2: what each key merges, and the per-class spread
    out["connectivity_groups_with_multiple_compound_keys"] = _int(
        sub.get("connectivity_groups_with_multiple_compound_keys"))
    out["distinctions_lost_at_connectivity"] = _int(
        sub.get("distinctions_lost_at_connectivity"))
    out["cross_entry_duplicates"] = _int_map(sub.get("cross_entry_duplicates"), KEYS)
    for key in ("records_per_class", "distinct_compounds_per_class",
                "distinct_scaffolds_per_class"):
        out[key] = _int_map(sub.get(key), classes)
    out["scaffold_ring_count_histogram"] = _free_int_map(
        sub.get("scaffold_ring_count_histogram"))

    # ---- R3: the sharing claim under the declared filter
    out["informative_rule"] = {
        "min_ring_count": _int(_nested(sub, "informative_rule", "min_ring_count")),
        "min_heavy_atom_coverage": _float(
            _nested(sub, "informative_rule", "min_heavy_atom_coverage")),
    }
    out["informative_records"] = _int(sub.get("informative_records"))
    out["informative_scaffold_keys"] = _int(sub.get("informative_scaffold_keys"))
    out["shared_scaffolds"] = _int_map(sub.get("shared_scaffolds"), SHARED_KEYS)
    out["cross_class_scaffolds"] = _int_map(sub.get("cross_class_scaffolds"),
                                           CROSS_KEYS)

    entries = sub.get("cross_class_informative_entries")
    entries = entries if isinstance(entries, dict) else {}
    out["cross_class_informative_entries"] = {
        str(k): sorted({str(a) for a in v}) for k, v in sorted(entries.items())
        if isinstance(v, list) and v}
    out["most_shared_informative_scaffold"] = {
        "scaffold_inchikey": (
            str(k) if isinstance(
                k := _nested(sub, "most_shared_informative_scaffold",
                             "scaffold_inchikey"), str) else None),
        "entries": _int(_nested(sub, "most_shared_informative_scaffold", "entries")),
    }

    # ---- R4: the threshold sweep and the evidence gate
    declared = sub.get("threshold_variants")
    declared = declared if isinstance(declared, dict) else {}
    records, scaffolds, sharing = {}, {}, {}
    for name in variants:
        block = declared.get(name)
        block = block if isinstance(block, dict) else {}
        records[name] = _int(block.get("informative_records"))
        scaffolds[name] = _int(block.get("informative_scaffold_keys"))
        sharing[name] = _int_map(block, VARIANT_SHARING_KEYS)
    out["threshold_variant_records"] = records
    out["threshold_variant_scaffolds"] = scaffolds
    out["threshold_variant_sharing"] = sharing

    gate = sub.get("evidence_gate")
    gate = gate if isinstance(gate, dict) else {}
    out["evidence_gate_counts"] = _int_map(gate, GATE_COUNT_KEYS)
    out["evidence_gate_sharing"] = {
        "shared_scaffolds": _int_map(gate.get("shared_scaffolds"), SHARED_KEYS),
        "cross_class_scaffolds": _int_map(gate.get("cross_class_scaffolds"),
                                         CROSS_KEYS),
    }
    keys = gate.get("cross_class_informative_keys")
    out["evidence_gate_cross_class_keys"] = (
        sorted({str(k) for k in keys}) if isinstance(keys, list) else [])
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
