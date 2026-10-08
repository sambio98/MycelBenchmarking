"""Fetch the chromosome and extract the curated boundaries for it.

Usage:
  python -m npbench_c.templates.bgc_detection.prepare_input <campaign> \
      --mibig <mibig_json_4.0.tar.gz> [--accession AL645882.2]

Run once at instantiation. The campaign ships the complete chromosome -- 6.5 MB
gzipped -- rather than a slice, because a slice both truncates CDS features at its
edges (which antiSMASH refuses) and sets the boundary of any region that reaches
it. That is the one input in this benchmark whose size cannot be reduced without
changing the answer.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pathlib
import tarfile
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TIMEOUT_S = 900


def fetch_genbank(accession: str) -> str:
    params = {"db": "nuccore", "id": accession, "rettype": "gbwithparts",
              "retmode": "text"}
    url = f"{EFETCH}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=TIMEOUT_S) as response:
        text = response.read().decode()
    if not text.startswith("LOCUS"):
        raise RuntimeError(f"efetch did not return a GenBank record: {text[:120]!r}")
    if "ORIGIN" not in text:
        raise RuntimeError("the record carries no sequence; rettype=gbwithparts is "
                           "required for a CON record such as a complete genome")
    return text


def mibig_loci(archive: pathlib.Path, accession: str) -> list[dict]:
    """Every active MIBiG locus mapped to this accession, with a real span."""
    loci: list[dict] = []
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".json"):
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            entry = json.loads(handle.read().decode())
            if entry.get("status") != "active":
                continue
            for locus in entry.get("loci") or []:
                if locus.get("accession") != accession:
                    continue
                location = locus.get("location") or {}
                start, end = location.get("from"), location.get("to")
                if not (isinstance(start, int) and isinstance(end, int)
                        and end > start):
                    continue
                loci.append({
                    "bgc": entry["accession"],
                    "from": start,
                    "to": end,
                    "completeness": entry.get("completeness") or "unknown",
                    "classes": sorted({c.get("class") for c in
                                       (entry.get("biosynthesis") or {}).get("classes")
                                       or []}),
                })
    loci.sort(key=lambda item: (item["from"], item["to"], item["bgc"]))
    return loci


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--mibig", required=True)
    parser.add_argument("--accession", default="AL645882.2")
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)

    text = fetch_genbank(args.accession)
    name = f"{args.accession}.gb.gz"
    target = inputs / name
    # mtime=0 so the shipped bytes are a function of the content alone.
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(target, "wb"),
                       mtime=0, compresslevel=9) as handle:
        handle.write(text.encode())

    loci = mibig_loci(pathlib.Path(args.mibig), args.accession)
    if not loci:
        raise RuntimeError(f"MIBiG maps no active locus to {args.accession}")
    loci_name = f"mibig_loci_{args.accession}.json"
    (inputs / loci_name).write_text(json.dumps(
        {"accession": args.accession, "release": "MIBiG 4.0", "loci": loci},
        indent=2, sort_keys=True) + "\n")

    header = text.splitlines()[0]
    provenance = {
        "sequence_source": "NCBI Nucleotide efetch, rettype=gbwithparts",
        "accession": args.accession,
        "locus_line": header.strip(),
        "license": "NCBI GenBank records carry no copyright claim; "
                   "MIBiG 4.0 is CC BY 4.0",
        "attribution": ("Streptomyces coelicolor A3(2) chromosome "
                        f"{args.accession} from NCBI Nucleotide; curated cluster "
                        "boundaries from MIBiG 4.0 (CC BY 4.0, "
                        "https://mibig.secondarymetabolites.org/)"),
        "files": {"sequence": name, "curated_loci": loci_name},
        "curated_loci": len(loci),
        "sha256": {
            name: hashlib.sha256(target.read_bytes()).hexdigest(),
            loci_name: hashlib.sha256(
                (inputs / loci_name).read_bytes()).hexdigest(),
        },
        "note": ("The COMPLETE chromosome, unsliced. A slice truncates CDS "
                 "features at its edges, which antiSMASH refuses, and sets the "
                 "boundary of any region that reaches it. All 29 regions on the "
                 "whole record come back contig_edge=False."),
    }
    (inputs / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"{name}: {header.split()[2]} bp, "
          f"{len(text.encode())/1e6:.1f} MB raw, "
          f"{target.stat().st_size/1e6:.1f} MB gzipped")
    print(f"{loci_name}: {len(loci)} curated loci on {args.accession}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
