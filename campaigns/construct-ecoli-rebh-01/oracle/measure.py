"""Measure a submission. Deterministic, campaign-specific, agent-independent.

The grader is a generic declarative comparator and must not contain campaign
logic; the agent cannot be trusted to self-report properties of its own
sequence. So measurement is a separate step that runs inside the container:

    submission/ --> measure.py --> measurement.json --> grade() --> grade.json

R1/R2 fields are *computed from the agent's sequence*, so they cannot be
faked. R3/R4 fields are the agent's analytical claims, normalised but not
verified here -- the grader compares them against oracle-computed gold.
"""

from __future__ import annotations

import json
import pathlib
import sys

import construct as C
from tables import CONSTRAINT_ORDER, HOSTS

REQUIRED_ORF_FILES = {"ecoli": "orf_ecoli.fasta", "streptomyces": "orf_streptomyces.fasta"}


def read_fasta(path: pathlib.Path) -> str:
    text = path.read_text().strip().splitlines()
    if not text or not text[0].startswith(">"):
        raise ValueError(f"{path.name}: not a FASTA record")
    return "".join(l.strip().upper() for l in text[1:] if not l.startswith(">"))


def _orf_schema(nt: str, protein: str) -> dict:
    """Structural checks. Every field is a 0/1 flag so the grader stays
    declarative and the failing check is visible component-wise."""
    valid_alphabet = bool(nt) and set(nt) <= set("ACGT")
    in_frame = valid_alphabet and len(nt) % 3 == 0
    codons = [nt[i:i + 3] for i in range(0, len(nt), 3)] if in_frame else []
    internal_stops = sum(
        1 for c in codons[:-1] if c in ("TAA", "TAG", "TGA")
    ) if codons else 1
    return {
        "valid_alphabet": int(valid_alphabet),
        "in_frame": int(in_frame),
        "starts_atg": int(bool(codons) and codons[0] == "ATG"),
        "ends_stop": int(bool(codons) and codons[-1] in ("TAA", "TAG", "TGA")),
        "no_internal_stop": int(internal_stops == 0),
        "length_matches_protein": int(bool(codons) and len(codons) == len(protein)),
    }


def _orf_constraints(nt: str, protein: str, host: str, strategy: str) -> dict:
    counts = C.violation_counts(nt, host, strategy)
    try:
        translated_ok = int(C.translate(nt) == protein)
    except (ValueError, KeyError):
        translated_ok = 0
    return {
        "protein_identity": translated_ok,
        "violations": {k: counts[k] for k in CONSTRAINT_ORDER},
        "total_violations": sum(counts.values()),
        "all_constraints_satisfied": int(sum(counts.values()) == 0),
        "cai": round(C.cai(nt, HOSTS[host]["weights"]), 6),
        "gc": round(C.gc_fraction(nt), 6),
    }


def measure(submission_dir: pathlib.Path, protein: str) -> dict:
    out: dict = {"present": {}, "ecoli": {}, "streptomyces": {}, "claims": {}}

    for key, fname in sorted(REQUIRED_ORF_FILES.items()):
        path = submission_dir / fname
        # Keys are dot-free: the grader resolves dotted paths, so a filename
        # with an extension would be read as two path segments.
        out["present"][f"orf_{key}"] = int(path.is_file())
        host = "ecoli_bl21" if key == "ecoli" else "streptomyces_coelicolor"
        strategy = "pet28a_ndei_xhoi" if key == "ecoli" else "pset152_ndei_xhoi"
        if not path.is_file():
            out[key] = {
                "schema": {k: 0 for k in _orf_schema("", protein)},
                "protein_identity": 0,
                "violations": {k: -1 for k in CONSTRAINT_ORDER},
                "total_violations": -1,
                "all_constraints_satisfied": 0,
                "cai": 0.0,
                "gc": 0.0,
            }
            continue
        try:
            nt = read_fasta(path)
        except ValueError:
            nt = ""
        schema = _orf_schema(nt, protein)
        if all(schema.values()):
            body = _orf_constraints(nt, protein, host, strategy)
        else:
            body = {
                "protein_identity": 0,
                "violations": {k: -1 for k in CONSTRAINT_ORDER},
                "total_violations": -1,
                "all_constraints_satisfied": 0,
                "cai": 0.0,
                "gc": round(C.gc_fraction(nt), 6) if nt else 0.0,
            }
        out[key] = {"schema": schema, **body}

    # Agent analytical claims, normalised only.
    claims_path = submission_dir / "analysis.json"
    out["present"]["analysis"] = int(claims_path.is_file())
    if claims_path.is_file():
        try:
            raw = json.loads(claims_path.read_text())
        except json.JSONDecodeError:
            raw = {}
        for rung in ("r3", "r4"):
            block = raw.get(rung) or {}
            uv = block.get("unconstrained_violations") or {}
            out["claims"][rung] = {
                "unconstrained_violations": {
                    k: (int(uv[k]) if isinstance(uv.get(k), (int, float)) else -1)
                    for k in CONSTRAINT_ORDER
                },
                "binding_constraint": block.get("binding_constraint"),
                "non_binding_constraints": sorted(
                    str(x) for x in (block.get("non_binding_constraints") or [])
                ),
            }
        fc = raw.get("r4", {}).get("flipped_constraints") or []
        out["claims"]["r4"]["flipped_constraints"] = sorted(str(x) for x in fc)
    return out


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: measure.py <submission_dir> <protein_fasta> <out_json>", file=sys.stderr)
        return 2
    sub = pathlib.Path(sys.argv[1])
    protein = read_fasta(pathlib.Path(sys.argv[2])) + "*"
    result = measure(sub, protein)
    pathlib.Path(sys.argv[3]).write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
