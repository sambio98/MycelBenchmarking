"""Assemble tasks/09_salinispora_paired_omics/data.manifest.json.

Three sources, because the task stages three kinds of input and each already
has a trustworthy record:

  * **Genomes** reuse task 02's manifest entries, filtered to this task's 111
    strains. The files are byte-identical copies, so copying the entries keeps
    one set of checksums rather than two that could drift.
  * **Spectra** come from work/spectra_fingerprints.json, written by
    scripts/stage_task09_spectra.py as it streamed each mzXML through.
  * **The consensus .mgf** reuses task 01's entry, for the same reason as the
    genomes.

    python scripts/build_task09_manifest.py --task tasks/09_salinispora_paired_omics

Fails loudly rather than writing a short manifest: every strain must have both
a genome entry and a spectra entry, or the bundle would stage a dataset the
solution cannot run on.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

GENOME_SOURCE = Path("tasks/02_salinispora_bgc/data.manifest.json")
MGF_SOURCE = Path("tasks/01_crusemann_networking/data.manifest.json")
MGF_PATH = "data/MSV000078836_specs_ms.mgf"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--genome-source", type=Path, default=GENOME_SOURCE)
    parser.add_argument("--mgf-source", type=Path, default=MGF_SOURCE)
    args = parser.parse_args()

    with (args.task / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))

    genome_entries = {e["path"]: e for e in json.loads(args.genome_source.read_text())}
    mgf_entries = {e["path"]: e for e in json.loads(args.mgf_source.read_text())}

    fingerprints_path = args.task / "work" / "spectra_fingerprints.json"
    if not fingerprints_path.is_file():
        raise SystemExit(f"no spectra fingerprints at {fingerprints_path}; "
                         f"run scripts/stage_task09_spectra.py first")
    spectra_entries = {e["path"]: e for e in json.loads(fingerprints_path.read_text())}

    entries, missing = [], []
    if MGF_PATH not in mgf_entries:
        raise SystemExit(f"{MGF_PATH} not in {args.mgf_source}")
    entries.append(mgf_entries[MGF_PATH])

    for strain in strains:
        for key, table, label in (
            (strain["genome_file"], genome_entries, "genome"),
            (strain["spectra_file"], spectra_entries, "spectra"),
        ):
            entry = table.get(key)
            if entry is None:
                missing.append(f"{strain['accession']}: no {label} entry for {key}")
            else:
                entries.append(entry)

    if missing:
        raise SystemExit("incomplete manifest, refusing to write:\n  "
                         + "\n  ".join(missing[:12])
                         + (f"\n  ... and {len(missing) - 12} more" if len(missing) > 12 else ""))

    # Sorted by path so the file is diffable, and with one key order throughout.
    fields = ("path", "url", "accession", "release", "bytes", "sha256")
    ordered = [
        {k: e[k] for k in fields if k in e}
        for e in sorted(entries, key=lambda e: e["path"])
    ]
    (args.task / "data.manifest.json").write_text(json.dumps(ordered, indent=2) + "\n")

    total = sum(e["bytes"] for e in ordered)
    spectra_total = sum(e["bytes"] for e in ordered if "/spectra/" in e["path"])
    print(f"{args.task / 'data.manifest.json'}: {len(ordered)} entries, "
          f"{total / 1e9:.2f} GB")
    print(f"  spectra: {spectra_total / 1e9:.2f} GB over "
          f"{sum(1 for e in ordered if '/spectra/' in e['path'])} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
