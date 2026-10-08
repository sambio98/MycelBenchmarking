"""Measure a domain-architecture submission. Deterministic, agent-independent.

Nothing here reads gold. Every key list comes from the campaign's own declared
vocabularies in `reference/architecture_rules.json` -- the panel, the cutoffs,
the coordinate spans, the resolution policies, the domain mappings, the curated
domain keys and the reconciliation verdicts -- or from this template's own closed
column lists.

The row-shaped answers are projected to their declared columns and left in the
submission's own order, because gold's order is declared too and a submission
that sorts differently has got the order wrong.

Usage:
  python -m npbench_c.templates.domain_architecture.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.domain_architecture.core import Rules

REPORT = "architecture_report.json"

ARCHITECTURE_COLUMNS = ("bgc", "gene", "length", "hits", "domains",
                        "resolved_away", "architecture", "first_domain_start",
                        "last_domain_end")
RECONCILIATION_COLUMNS = ("bgc", "gene", "stratum", "curated_modules",
                          "curated_primary", "computed_primary",
                          "primary_delta", "primary_verdict",
                          "curated_secondary", "computed_secondary",
                          "secondary_verdict")
MODULE_COLUMNS = ("bgc", "gene", "complete_modules", "curated_modules",
                  "verdict")
LENGTH_KEYS = ("count", "min", "median", "max", "residues")
ARCH_CENSUS_KEYS = ("genes", "genes_with_no_domain", "distinct_architectures",
                    "domains_total")
OVERLAP_KEYS = ("hits_total", "overlapping_pairs", "genes_with_an_overlap")
SELECTION_KEYS = ("entries_with_modules", "modules",
                  "modules_spanning_several_genes",
                  "genes_named_by_a_single_gene_module",
                  "drop_excluded_by_status", "drop_other_module_type",
                  "drop_no_shipped_translation", "shipped_proteins",
                  "corpus_genes", "corpus_entries")
SWEEP_CELL_KEYS = ("domains_total", "distinct_architectures",
                   "a_domain_agreement", "genes_differing_from_baseline")


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _opt_int(value: object) -> int | None:
    return None if value is None else _int(value)


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _str_list(value: object) -> list[str]:
    return [str(v) for v in value] if isinstance(value, list) else []


def _drop_zeros(block: dict) -> dict:
    """Keys whose count is zero are omitted, as gold omits them."""
    return {k: v for k, v in block.items() if v}


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped: gold omits them.

    The curated module census is keyed by module count, the delta census by a
    signed offset, and the `active`-flag census by the literal values the field
    holds. None is a declared key list; they are whatever the data turns out to
    contain.
    """
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def _rows(value: object, columns, casts) -> list[dict]:
    rows = value if isinstance(value, list) else []
    out = []
    for row in rows:
        block = row if isinstance(row, dict) else {}
        out.append({c: casts[c](block.get(c)) for c in columns})
    return out


_ARCHITECTURE_CASTS = {
    "bgc": _str, "gene": _str, "length": _int, "hits": _int, "domains": _int,
    "resolved_away": _int, "architecture": _str,
    "first_domain_start": _opt_int, "last_domain_end": _opt_int,
}
_RECONCILIATION_CASTS = {
    "bgc": _str, "gene": _str, "stratum": _str, "curated_modules": _int,
    "curated_primary": _int, "computed_primary": _int, "primary_delta": _int,
    "primary_verdict": _str, "curated_secondary": _int,
    "computed_secondary": _int, "secondary_verdict": _str,
}
_MODULE_CASTS = {
    "bgc": _str, "gene": _str, "complete_modules": _int,
    "curated_modules": _int, "verdict": _str,
}


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

    # ---- R1: the selection, the curated side's own shape, and the environment
    out["selection_census"] = _int_map(sub.get("selection_census"),
                                       SELECTION_KEYS)
    out["curated_coordinate_usability"] = _int_map(
        sub.get("curated_coordinate_usability"), ("placeholder", "usable"))
    out["curated_active_flag"] = _free_int_map(sub.get("curated_active_flag"))
    out["curated_module_census"] = _free_int_map(sub.get("curated_module_census"))
    out["strata_census"] = _free_int_map(sub.get("strata_census"))
    per_type = sub.get("domains_by_module_type")
    per_type = per_type if isinstance(per_type, dict) else {}
    out["domains_by_module_type"] = {
        str(name): _drop_zeros(_int_map(block, rules.curated_domain_keys))
        for name, block in sorted(per_type.items()) if isinstance(block, dict)}
    out["curated_domain_totals"] = _int_map(sub.get("curated_domain_totals"),
                                            rules.curated_domain_keys)
    length = sub.get("protein_length_summary")
    length = length if isinstance(length, dict) else {}
    out["protein_length_summary"] = {
        k: (_float(length.get(k)) if k == "median" else _int(length.get(k)))
        for k in LENGTH_KEYS}
    out["panel_census"] = _int_map(sub.get("panel_census"), rules.panel)
    out["environment"] = {
        "hmmer_version": _str(sub.get("hmmer_version")),
        "pfam_version": _str(sub.get("pfam_version")),
        "pfam_sha256": _str(sub.get("pfam_sha256")),
    }

    # ---- R2: the architecture
    out["architecture_table"] = _rows(sub.get("architecture_table"),
                                      ARCHITECTURE_COLUMNS, _ARCHITECTURE_CASTS)
    census = sub.get("architecture_census")
    census = census if isinstance(census, dict) else {}
    per_gene = census.get("domains_per_gene")
    per_gene = per_gene if isinstance(per_gene, dict) else {}
    out["architecture_census"] = {
        **_int_map(census, ARCH_CENSUS_KEYS),
        "domains_per_gene": {"min": _int(per_gene.get("min")),
                             "median": _float(per_gene.get("median")),
                             "max": _int(per_gene.get("max"))},
    }
    out["overlap_census"] = _int_map(sub.get("overlap_census"), OVERLAP_KEYS)

    # ---- R3: against the curated modules
    out["reconciliation_table"] = _rows(sub.get("reconciliation_table"),
                                        RECONCILIATION_COLUMNS,
                                        _RECONCILIATION_CASTS)
    out["primary_verdicts"] = _int_map(sub.get("primary_verdicts"),
                                       rules.verdicts)
    out["secondary_verdicts"] = _int_map(sub.get("secondary_verdicts"),
                                         rules.verdicts)
    by_stratum = sub.get("secondary_verdicts_by_stratum")
    by_stratum = by_stratum if isinstance(by_stratum, dict) else {}
    out["secondary_verdicts_by_stratum"] = {
        name: _int_map(by_stratum.get(name), rules.verdicts)
        for name in sorted(set(rules.strata) | {"mixed"})}
    out["primary_delta_census"] = _free_int_map(sub.get("primary_delta_census"))
    out["module_table"] = _rows(sub.get("module_table"), MODULE_COLUMNS,
                                _MODULE_CASTS)
    out["module_verdicts"] = _int_map(sub.get("module_verdicts"), rules.verdicts)
    out["complete_modules_total"] = _int(sub.get("complete_modules_total"))

    # ---- R4: the declared readings
    sweep = sub.get("policy_sweep")
    sweep = sweep if isinstance(sweep, dict) else {}
    cells = [f"{c}|{p}|{s}" for c in rules.cutoffs
             for p in rules.resolution_policies for s in rules.coordinate_spans]
    out["policy_sweep"] = {cell: _int_map(sweep.get(cell), SWEEP_CELL_KEYS)
                           for cell in cells}
    mapping = sub.get("mapping_sweep")
    mapping = mapping if isinstance(mapping, dict) else {}
    out["mapping_sweep"] = {}
    for name in sorted(rules.domain_mappings):
        if name == rules.secondary_domain:
            continue
        block = mapping.get(name)
        block = block if isinstance(block, dict) else {}
        out["mapping_sweep"][name] = {
            "families": _str_list(block.get("families")),
            "verdicts": _int_map(block.get("verdicts"), rules.verdicts),
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
