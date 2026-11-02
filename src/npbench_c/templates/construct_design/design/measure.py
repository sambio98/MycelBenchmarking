"""Measure a design-ladder submission. Deterministic, agent-independent.

R1/R2 fields are computed from the agent's own sequence, so they cannot be
faked. R3/R4 fields are the agent's analytical claims, normalised but not
verified here -- the grader compares them to oracle-computed gold.

Usage:
  python -m npbench_c.templates.construct_design.design.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.construct_design.core import ConstructEngine
from npbench_c.templates.construct_design.params import load_task, read_protein

STOPS = ("TAA", "TAG", "TGA")


def read_fasta(path: pathlib.Path) -> str:
    lines = path.read_text().strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"{path.name}: not a FASTA record")
    return "".join(l.strip().upper() for l in lines[1:] if not l.startswith(">"))


def _schema(nt: str, protein: str) -> dict:
    """Structural checks, each a 0/1 flag so the failing check stays visible
    component-wise and the grader stays declarative."""
    valid_alphabet = bool(nt) and set(nt) <= set("ACGT")
    in_frame = valid_alphabet and len(nt) % 3 == 0
    codons = [nt[i:i + 3] for i in range(0, len(nt), 3)] if in_frame else []
    internal_stops = sum(1 for c in codons[:-1] if c in STOPS) if codons else 1
    return {
        "valid_alphabet": int(valid_alphabet),
        "in_frame": int(in_frame),
        "starts_atg": int(bool(codons) and codons[0] == "ATG"),
        "ends_stop": int(bool(codons) and codons[-1] in STOPS),
        "no_internal_stop": int(internal_stops == 0),
        "length_matches_protein": int(bool(codons) and len(codons) == len(protein)),
    }


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    engine = ConstructEngine.load(campaign)
    params = load_task(campaign)["template_params"]
    protein = read_protein(campaign, params["protein_input"])
    constraints = list(engine.constraint_order)
    out: dict = {"present": {}, "claims": {}}

    for key, spec in sorted(params["hosts"].items()):
        path = submission_dir / f"orf_{key}.fasta"
        # Keys are dot-free: the grader resolves dotted paths, so a filename
        # with an extension would be read as two path segments.
        out["present"][f"orf_{key}"] = int(path.is_file())
        if not path.is_file():
            out[key] = {"schema": {k: 0 for k in _schema("", protein)},
                        "protein_identity": 0,
                        "violations": {c: -1 for c in constraints},
                        "total_violations": -1, "all_constraints_satisfied": 0,
                        "cai": 0.0, "gc": 0.0}
            continue
        try:
            nt = read_fasta(path)
        except ValueError:
            nt = ""
        schema = _schema(nt, protein)
        if all(schema.values()):
            counts = engine.violation_counts(nt, spec["host"], spec["strategy"])
            try:
                identity = int(engine.translate(nt) == protein)
            except (ValueError, KeyError):
                identity = 0
            body = {
                "protein_identity": identity,
                "violations": {c: counts[c] for c in constraints},
                "total_violations": sum(counts.values()),
                "all_constraints_satisfied": int(sum(counts.values()) == 0),
                "cai": round(engine.cai(nt, engine.hosts[spec["host"]]["weights"]), 6),
                "gc": round(engine.gc_fraction(nt), 6),
            }
        else:
            body = {"protein_identity": 0,
                    "violations": {c: -1 for c in constraints},
                    "total_violations": -1, "all_constraints_satisfied": 0,
                    "cai": 0.0,
                    "gc": round(engine.gc_fraction(nt), 6) if nt else 0.0}
        out[key] = {"schema": schema, **body}

    claims_path = submission_dir / "analysis.json"
    out["present"]["analysis"] = int(claims_path.is_file())
    if claims_path.is_file():
        try:
            raw = json.loads(claims_path.read_text())
        except json.JSONDecodeError:
            raw = {}
        raw = raw if isinstance(raw, dict) else {}
        for rung in ("r3", "r4"):
            block = raw.get(rung) or {}
            block = block if isinstance(block, dict) else {}
            uv = block.get("unconstrained_violations") or {}
            uv = uv if isinstance(uv, dict) else {}
            out["claims"][rung] = {
                "unconstrained_violations": {
                    c: (int(uv[c]) if isinstance(uv.get(c), int)
                        and not isinstance(uv.get(c), bool) else -1)
                    for c in constraints
                },
                "binding_constraint": block.get("binding_constraint"),
                "non_binding_constraints": sorted(
                    str(x) for x in (block.get("non_binding_constraints") or [])),
            }
        r4 = raw.get("r4") or {}
        r4 = r4 if isinstance(r4, dict) else {}
        out["claims"]["r4"]["flipped_constraints"] = sorted(
            str(x) for x in (r4.get("flipped_constraints") or []))
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>", file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
