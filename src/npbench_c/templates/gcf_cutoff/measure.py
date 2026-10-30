"""Measure a GCF-cutoff submission. Deterministic, agent-independent.

Every vocabulary comes from the campaign's own declared rules -- the cutoff grid,
the verdict enum, the family names, the control list -- so a submission is never
scored against something it was not told.

Usage:
  python -m npbench_c.templates.gcf_cutoff.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.gcf_cutoff.analyse import load_params
from npbench_c.templates.gcf_cutoff.core import Clustering, Rules

REPORT = "clustering_report.json"

SWEEP_FIELDS = ("groups", "singletons", "recovered_pairs", "missed_pairs",
                "false_join_pairs", "pair_jaccard", "families_exact")
PAIR_FIELDS = ("designed_pairs", "co_clustered_pairs", "recovered_pairs",
               "missed_pairs", "false_join_pairs", "pair_jaccard")
CENSUS_FIELDS = ("groups", "largest_group", "singletons", "clustered_records")


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


def _str_list(value: object) -> list:
    return sorted({str(v) for v in value}) if isinstance(value, list) else []


def _partition(value: object) -> list:
    """Groups as sorted member lists, sorted among themselves.

    A family label is a name rather than a result, so what is compared is the
    partition: two runs can agree on which clusters belong together and disagree
    on what the groups are called.
    """
    groups = value if isinstance(value, list) else []
    out = []
    for group in groups:
        if isinstance(group, list):
            out.append(sorted({str(m) for m in group}))
    return sorted(out)


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    params = load_params(campaign)
    clustering = Clustering.load(campaign, params["cluster_dir"],
                                 params["designed_families_input"])
    cutoffs = tuple(rules.label(c) for c in rules.cutoffs)
    families = tuple(sorted(clustering.designed))
    verdicts = rules.verdict_order

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

    # ---- R1: the corpus, and what clustered it
    out["corpus_size"] = _int(sub.get("corpus_size"))
    out["designed_families"] = _int_map(sub.get("designed_families"), families)
    out["controls"] = _str_list(sub.get("controls"))
    out["cutoffs"] = (list(sub["cutoffs"]) if isinstance(sub.get("cutoffs"), list)
                      else [])
    out["baseline_cutoff"] = _str(sub.get("baseline_cutoff"))
    out["include_gbk"] = _str(sub.get("include_gbk"))
    out["class_census"] = _free_int_map(sub.get("class_census"))
    out["environment"] = {
        "bigscape_version": _str(sub.get("bigscape_version")),
        "pfam_version": _str(sub.get("pfam_version")),
        "pfam_sha256": _str(sub.get("pfam_sha256")),
    }

    # ---- R2: the partition at the declared cutoff
    out["baseline_partition"] = _partition(sub.get("baseline_partition"))
    census = _int_map(sub, CENSUS_FIELDS)
    census["size_census"] = _free_int_map(sub.get("size_census"))
    out["group_census"] = census

    # ---- R3: the reconciliation with the chemistry-derived families
    declared = sub.get("per_family")
    declared = declared if isinstance(declared, dict) else {}
    per_family = {}
    for name in families:
        block = declared.get(name)
        block = block if isinstance(block, dict) else {}
        verdict = _str(block.get("verdict"))
        sizes = block.get("group_sizes")
        per_family[name] = {
            "verdict": verdict if verdict in verdicts else None,
            "groups_spanned": _int(block.get("groups_spanned")),
            "group_sizes": (sorted((n for s in sizes if (n := _int(s)) is not None),
                                   reverse=True)
                            if isinstance(sizes, list) else []),
            "outsiders": _str_list(block.get("outsiders")),
        }
    out["per_family"] = per_family
    out["verdict_census"] = _int_map(sub.get("verdict_census"), verdicts)
    out["families_exact"] = _int(sub.get("families_exact"))
    out["controls_absorbed"] = _str_list(sub.get("controls_absorbed"))
    pairs = _int_map(sub, PAIR_FIELDS)
    pairs["pair_jaccard"] = _float(sub.get("pair_jaccard"))
    out["pair_statistics"] = pairs

    # ---- R4: the sweep
    sweep = sub.get("cutoff_sweep")
    sweep = sweep if isinstance(sweep, dict) else {}
    summary = {}
    for label in cutoffs:
        row = sweep.get(label)
        row = row if isinstance(row, dict) else {}
        block = {k: _int(row.get(k)) for k in SWEEP_FIELDS if k != "pair_jaccard"}
        block["pair_jaccard"] = _float(row.get("pair_jaccard"))
        summary[label] = block
    out["cutoff_sweep_summary"] = summary
    out["sweep_boundaries"] = {
        "best_cutoff": _str(sub.get("best_cutoff")),
        "best_pair_jaccard": _float(sub.get("best_pair_jaccard")),
        "first_cutoff_with_a_false_join":
            _str(sub.get("first_cutoff_with_a_false_join")),
        "first_cutoff_absorbing_a_control":
            _str(sub.get("first_cutoff_absorbing_a_control")),
        # A family need not be recoverable at any declared cutoff, so null is a
        # real answer here and must not be normalised to a string.
        "first_cutoff_recovering_each_family": {
            name: _str((sub.get("first_cutoff_recovering_each_family") or {}).get(name))
            if isinstance(sub.get("first_cutoff_recovering_each_family"), dict)
            else None
            for name in families},
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
