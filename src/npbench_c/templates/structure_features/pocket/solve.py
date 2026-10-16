"""Reference solution for the pocket ladder.

Renders gold into the submission shape the output contract describes, which
makes the contract testable: if the oracle cannot express gold in the contracted
shape, the contract is wrong.

Usage: python -m npbench_c.templates.structure_features.pocket.solve \
    <campaign> <submission>
"""

from __future__ import annotations

import json
import pathlib
import sys

FLAT_KEYS = (
    "pdb_id", "accession", "atom_inventory", "polymer_chains", "ligand_copies",
    "numbering", "resolved_residues_per_chain", "resolved_span_per_chain",
    "sequence_agreement",
)


def render(gold: dict) -> dict:
    report = {k: gold[k] for k in FLAT_KEYS}
    report["shells"] = gold["shells"]
    report["descriptors"] = gold["descriptors"]
    report["counterfactuals"] = gold["counterfactuals"]
    return report


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: solve <campaign_dir> <submission_dir>", file=sys.stderr)
        return 2
    campaign, out = pathlib.Path(argv[0]), pathlib.Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    (out / "pocket_report.json").write_text(
        json.dumps(render(gold), indent=2, sort_keys=True) + "\n")
    print(f"oracle submission written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
