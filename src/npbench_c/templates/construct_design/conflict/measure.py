"""Measure a conflict-ladder submission. Deterministic, agent-independent.

Usage:
  python -m npbench_c.templates.construct_design.conflict.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.construct_design.core import (
    CODON_TO_AA,
    ConstructEngine,
    gc_count,
)
from npbench_c.templates.construct_design.params import (
    CONFLICT_CONSTRAINTS,
    VERDICTS,
    ConflictSpec,
)


def _float_or_none(v: object) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _int_or_none(v: object) -> int | None:
    return int(v) if isinstance(v, int) and not isinstance(v, bool) else None


def _bool_or_none(v: object) -> bool | None:
    return bool(v) if isinstance(v, bool) else None


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    engine = ConstructEngine.load(campaign)
    spec = ConflictSpec.load(campaign)

    path = submission_dir / "feasibility_report.json"
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

    out["orf_length_nt"] = _int_or_none(sub.get("orf_length_nt"))
    assessed = sub.get("constraints_assessed")
    out["constraints_assessed"] = sorted(str(x) for x in assessed) if isinstance(assessed, list) else []

    counts = sub.get("admissible_codon_counts")
    counts = counts if isinstance(counts, dict) else {}
    out["admissible_codon_counts"] = {
        str(k): _int_or_none(v) for k, v in sorted(counts.items())
    }
    out["codon_counts_reported"] = int(bool(out["admissible_codon_counts"]))

    for key in ("unrestricted_gc_min", "unrestricted_gc_max",
                "restricted_gc_min", "restricted_gc_max"):
        out[key] = _float_or_none(sub.get(key))

    verdict = sub.get("verdict")
    out["verdict"] = verdict if verdict in VERDICTS else None

    pair = sub.get("conflicting_constraints")
    out["conflicting_constraints"] = (
        sorted(str(x) for x in pair) if isinstance(pair, list) else [])

    witness = sub.get("witness")
    witness = witness if isinstance(witness, dict) else {}
    out["witness"] = {
        "required_gc_min": _float_or_none(witness.get("required_gc_min")),
        "max_achievable_gc_under_forbidden_codons": _float_or_none(
            witness.get("max_achievable_gc_under_forbidden_codons")),
    }

    # Supplying a design is a fabricated feasibility claim: no sequence encoding
    # this protein can satisfy the set, so there is nothing to supply. The
    # measurer also checks any supplied sequence, so the audit packet records
    # what was wrong with it rather than only that it was offered.
    design = sub.get("design_orf")
    out["claims_satisfying_design"] = int(isinstance(design, str) and bool(design.strip()))
    out["supplied_design_check"] = _check_design(engine, spec, design)

    relax = sub.get("relaxations")
    relax = relax if isinstance(relax, dict) else {}
    out["relaxations"] = {}
    for c in CONFLICT_CONSTRAINTS:
        block = relax.get(c) or {}
        block = block if isinstance(block, dict) else {}
        out["relaxations"][c] = {
            "max_achievable_gc": _float_or_none(block.get("max_achievable_gc")),
            "gc_requirement_reachable": _bool_or_none(block.get("gc_requirement_reachable")),
        }
    return out


def _check_design(engine: ConstructEngine, spec: ConflictSpec,
                  design: object) -> dict:
    """Diagnostics on a supplied design. Never a pass condition."""
    if not isinstance(design, str) or not design.strip():
        return {"supplied": False}
    nt = design.strip().upper()
    info: dict = {"supplied": True, "length_nt": len(nt)}
    if len(nt) % 3 != 0 or not set(nt) <= set("ACGT"):
        info["parses"] = False
        return info
    info["parses"] = True
    codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    forbidden = set(spec.forbidden_codons)
    info["forbidden_codon_uses"] = sum(1 for c in codons if c in forbidden)
    info["gc_fraction"] = round(engine.gc_fraction(nt), 6)
    info["meets_gc_minimum"] = bool(info["gc_fraction"] >= spec.gc_global_min)
    info["encodes_declared_protein"] = all(c in CODON_TO_AA for c in codons)
    return info


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
