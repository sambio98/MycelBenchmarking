"""Reference solution for the residues ladder.

Renders gold into the submission shape the output contract describes, which is
what makes the contract itself testable: if the oracle cannot express gold in
the contracted shape, the contract is wrong.

Usage: python -m npbench_c.templates.structure_features.residues.solve \
    <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

FLAT_KEYS = (
    "accession", "sequence_length", "functional_record_counts",
    "secondary_structure_counts", "annotated_positions", "ligand_vocabulary",
    "pdb_cross_references", "mutagenesis_original_residue_agreements",
    "mutagenesis_original_residue_disagreements", "eco_code_citation_totals",
    "record_class_counts", "binding_positions_by_ligand", "pdb_evidence_index",
    "unevidenced_records", "inferred_only_records",
)


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in FLAT_KEYS}
    report["focal_structure"] = gold["focal_structure"]
    focal = gold["confirmation"]
    report["confirmation"] = {
        "structure": focal["structure"],
        "chain": focal["chain"],
        "numbering": focal["numbering"],
        "per_ligand": focal["per_ligand"],
        "totals": focal["totals"],
        "active_site_position": focal["active_site_position"],
        "active_site_residue": focal["active_site_residue"],
        "active_site_contacts": focal["active_site_contacts"],
    }
    swap = gold["structure_swap"]
    report["structure_swap"] = {
        "structure": swap["structure"],
        "numbering": swap["numbering"],
        "per_ligand": swap["per_ligand"],
        "confirmed_delta_by_ligand": swap["confirmed_delta_by_ligand"],
    }
    report["policy_variants"] = gold["policy_variants"]
    return report


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: solve <campaign_dir> <submission_dir>", file=sys.stderr)
        return 2
    campaign, out = pathlib.Path(argv[0]), pathlib.Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    (out / "residue_report.json").write_text(
        json.dumps(render(gold), indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
