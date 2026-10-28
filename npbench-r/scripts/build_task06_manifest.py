"""Build tasks/06_streptomyces_ripp/data.manifest.json from NCBI RefSeq.

One-off operational script. The genome set is defined by the fixed query below
rather than a hand-curated list, so it can be re-derived; the resulting
manifest is the frozen record.

Selection rule: every NCBI RefSeq **reference genome** of a named *Streptomyces*
species whose assembly is complete and current. NCBI's own reference-genome
curation does the sampling, which avoids both an arbitrary cap and the
over-representation of heavily sequenced species -- there are 1,543 complete
Streptomyces assemblies but only 213 reference genomes, exactly one per species.

    python scripts/build_task06_manifest.py
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from npbench.build.manifest import Entry, describe, dump_manifest, sha256_file  # noqa: E402

TASK = Path(__file__).resolve().parents[1] / "tasks" / "06_streptomyces_ripp"
TAXON = "Streptomyces"
DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


def _get(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "npbench-r/0.1"})
    with urllib.request.urlopen(request, timeout=300) as response:
        return json.loads(response.read())


def select_assemblies() -> list[dict]:
    query = urllib.parse.urlencode(
        {
            "page_size": 1000,
            "filters.assembly_source": "refseq",
            "filters.assembly_level": "complete_genome",
            "filters.assembly_version": "current",
            "filters.reference_only": "true",
        }
    )
    report = _get(f"{DATASETS_API}/genome/taxon/{TAXON}/dataset_report?{query}")
    chosen = []
    for item in report.get("reports", []):
        organism = item["organism"]["organism_name"]
        if organism.split()[1] == "sp.":
            continue  # an unnamed species cannot carry a species-level claim
        chosen.append(
            {
                "accession": item["accession"],
                "organism": organism,
                "species": " ".join(organism.split()[:2]),
                "strain": item["organism"].get("infraspecific_names", {}).get("strain", ""),
                "name": item["assembly_info"]["assembly_name"],
            }
        )
    return sorted(chosen, key=lambda a: a["accession"])


def ftp_url(accession: str, assembly_name: str) -> str:
    digits = accession.split("_")[1].split(".")[0]
    triplets = "/".join(digits[i:i + 3] for i in range(0, 9, 3))
    stem = f"{accession}_{assembly_name.replace(' ', '_')}"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{triplets}/{stem}/"
        f"{stem}_genomic.fna.gz"
    )


def main() -> int:
    assemblies = select_assemblies()
    print(f"selected {len(assemblies)} reference genomes", flush=True)

    entries: list[Entry] = []
    rows = []
    for index, assembly in enumerate(assemblies, start=1):
        accession = assembly["accession"]
        target = f"data/genomes/{accession}.fna.gz"
        staged = TASK / target
        url = ftp_url(accession, assembly["name"])
        try:
            if staged.is_file():
                entry = Entry(path=target, bytes=staged.stat().st_size,
                              sha256=sha256_file(staged), url=url, accession=accession)
            else:
                entry = describe(url, target, accession=accession)
        except Exception as exc:
            print(f"  !! {accession}: {exc}", flush=True)
            return 1
        entries.append(entry)
        rows.append({**assembly, "file": target})
        if index % 25 == 0 or index == len(assemblies):
            print(f"  {index}/{len(assemblies)}", flush=True)

    dump_manifest(entries, TASK / "data.manifest.json")
    total = sum(e.bytes for e in entries)
    print(f"wrote manifest: {len(entries)} files, {total / 1e6:.1f} MB")

    (TASK / "strains.tsv").write_text(
        "accession\torganism\tspecies\tstrain\tfile\n"
        + "".join(
            f"{r['accession']}\t{r['organism']}\t{r['species']}\t{r['strain']}\t{r['file']}\n"
            for r in rows
        )
    )
    print(f"wrote {TASK / 'strains.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
