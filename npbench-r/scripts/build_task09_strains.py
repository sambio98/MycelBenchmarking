"""Write tasks/09_salinispora_paired_omics/strains.tsv from PoDP and task 02.

The Paired Omics Data Platform curates, for the Crusemann deposit, which mass
spectrometry file belongs to which sequenced strain. That mapping is the thing
that makes a paired-omics task possible at all: the clustered .mgf in task 01
records `FILENAME=specs_ms.pklbin` for every one of its 100,263 spectra, so it
cannot say which strain a consensus spectrum came from, and MassIVE's own
listing endpoints are unreachable from a current certificate store.

PoDP gives 360 genome-to-file links: 120 strains across three extraction
methods. 111 of those strains also have an NCBI assembly already staged for
task 02. This script keeps those 111 at one extraction method, which holds the
analytical batch constant -- same sample preparation, same instrumentation,
same deposit directory -- because comparing metabolomes across extraction
methods would confound every occurrence call downstream.

    python scripts/build_task09_strains.py --podp <project.json> \
        --strains tasks/02_salinispora_bgc/strains.tsv --out <strains.tsv>
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

EXTRACTION = "Ethyl acetate"

# RefSeq holds two independent assemblies of Salinispora tropica CNS416, both at
# scaffold level. The strain is the unit of every occurrence and co-occurrence
# count downstream, so one has to go or that strain is counted twice in the
# linkage matrix. The rule is deterministic and recorded rather than silently
# applied: prefer the better assembly level, then the lexicographically smallest
# accession.
ASSEMBLY_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
MASSIVE_PREFIX = "ftp://massive.ucsd.edu/MSV000078836/"
PROXY = "https://massiveproxy.gnps2.org/massiveproxy/MSV000078836/"


def normalise_strain(label: str) -> str:
    """Strain codes disagree between sources: CNT849 here, CNT-849 there."""
    return re.sub(r"[^A-Z0-9]", "", label.upper())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--podp", type=Path, required=True)
    parser.add_argument("--strains", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    project = json.loads(args.podp.read_text())
    links = [
        link for link in project["genome_metabolome_links"]
        if link["extraction_method_label"] == EXTRACTION
    ]
    by_strain = {
        normalise_strain(link["genome_label"].split()[-1]): link for link in links
    }

    with args.strains.open() as handle:
        assemblies = list(csv.DictReader(handle, delimiter="\t"))

    rows = []
    for assembly in assemblies:
        key = normalise_strain(assembly["strain"])
        link = by_strain.get(key)
        if link is None:
            continue
        path = link["metabolomics_file"]
        if not path.startswith(MASSIVE_PREFIX):
            raise SystemExit(f"unexpected file location: {path}")
        relative = path[len(MASSIVE_PREFIX):]
        rows.append({
            "accession": assembly["accession"],
            "organism": assembly["organism"],
            "strain": assembly["strain"],
            "assembly_level": assembly["assembly_level"],
            "genome_file": assembly["file"],
            "spectra_file": f"data/spectra/{Path(relative).name}",
            "spectra_url": PROXY + relative,
            "podp_genome_label": link["genome_label"],
        })

    # One assembly per strain, by the rule above.
    best: dict[str, dict] = {}
    for row in rows:
        key = normalise_strain(row["strain"])
        rank = (ASSEMBLY_RANK.get(row["assembly_level"], 99), row["accession"])
        if key not in best or rank < best[key]["_rank"]:
            best[key] = {**row, "_rank": rank}
    dropped = len(rows) - len(best)
    rows = [{k: v for k, v in row.items() if k != "_rank"}
            for row in best.values()]
    rows.sort(key=lambda r: r["accession"])
    if dropped:
        print(f"dropped {dropped} duplicate assembly/assemblies of an already-paired strain")

    if len({r["spectra_file"] for r in rows}) != len(rows):
        raise SystemExit("two strains still map to the same spectrum file")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{args.out}: {len(rows)} paired strains at '{EXTRACTION}'")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
