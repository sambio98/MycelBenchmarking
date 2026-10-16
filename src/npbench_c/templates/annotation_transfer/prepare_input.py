"""Fetch the reference set from UniProt and pin it as a closed FASTA plus annotations.

Usage: python -m npbench_c.templates.annotation_transfer.prepare_input <campaign> \
           [--ec-prefix 1.14]

Run once at instantiation. The campaign ships both files, so gold does not depend
on UniProt being reachable or unchanged at grading time, and the reference set the
transfer searches against is closed rather than live -- which is the difference
between a reproducible experiment and a snapshot of one afternoon.

The whole reviewed set for the prefix, no sampling: the queries and the reference
are the same file, so a sampled set would silently change both sides at once.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pathlib
import urllib.parse
import urllib.request

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
TIMEOUT_S = 900


def _get(url: str, params: dict) -> tuple[bytes, dict]:
    query = urllib.parse.urlencode(params)
    request = urllib.request.Request(f"{url}?{query}",
                                     headers={"Accept-Encoding": "identity"})
    with urllib.request.urlopen(request, timeout=TIMEOUT_S) as response:
        return response.read(), {k.lower(): v for k, v in response.headers.items()}


def _write_gzip(path: pathlib.Path, text: str) -> None:
    # mtime=0 so the shipped bytes are a function of the content alone.
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(path, "wb"),
                       mtime=0) as handle:
        handle.write(text.encode())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--ec-prefix", default="1.14")
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)

    query = f"ec:{args.ec_prefix}.* AND reviewed:true"
    head, headers = _get(SEARCH, {"query": query, "size": "1",
                                  "fields": "accession", "format": "tsv"})
    expected = int(headers["x-total-results"])
    release = headers["x-uniprot-release"]

    fasta_body, _ = _get(STREAM, {"query": query, "format": "fasta"})
    fasta = fasta_body.decode()
    got = sum(1 for line in fasta.splitlines() if line.startswith(">"))
    if got != expected:
        raise RuntimeError(
            f"the search reports {expected} entries and the stream returned "
            f"{got}; a partial download would change both the query set and the "
            "reference set at once")

    table_body, _ = _get(STREAM, {"query": query, "format": "tsv",
                                  "fields": "accession,ec"})

    # The search needs bare accessions: a hit table keyed on
    # 'sp|P00001|NAME_ORG' would have to be re-parsed everywhere downstream, and
    # the header's description field can itself contain the separator.
    rewritten = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            if len(parts) < 3:
                raise RuntimeError(f"unexpected FASTA header: {line[:60]!r}")
            rewritten.append(">" + parts[1])
        else:
            rewritten.append(line)

    annotations: dict[str, list[str]] = {}
    rows = table_body.decode().splitlines()
    for row in rows[1:]:
        accession, _, codes = row.partition("\t")
        annotations[accession.strip()] = [c.strip() for c in codes.split(";")
                                          if c.strip()]
    if len(annotations) != expected:
        raise RuntimeError(
            f"{len(annotations)} annotation rows for {expected} entries")

    prefix = args.ec_prefix.replace(".", "_")
    fasta_name = f"ec_{prefix}_reviewed_uniprot_{release}.fasta.gz"
    annotation_name = f"ec_{prefix}_reviewed_annotations_{release}.json.gz"
    _write_gzip(inputs / fasta_name, "\n".join(rewritten) + "\n")
    _write_gzip(inputs / annotation_name, json.dumps(
        {"release": release, "ec_numbers": dict(sorted(annotations.items()))},
        indent=2, sort_keys=True) + "\n")

    provenance = {
        "source": "UniProtKB REST",
        "release": release,
        "license": "CC BY 4.0",
        "attribution": ("UniProt Consortium, UniProtKB release "
                        f"{release} (https://www.uniprot.org/), CC BY 4.0"),
        "query": query,
        "entries": expected,
        "files": {
            "sequences": fasta_name,
            "annotations": annotation_name,
        },
        "sha256": {
            fasta_name: hashlib.sha256((inputs / fasta_name).read_bytes()).hexdigest(),
            annotation_name: hashlib.sha256(
                (inputs / annotation_name).read_bytes()).hexdigest(),
        },
        "note": ("The whole reviewed set for the prefix, no sampling. FASTA "
                 "headers are reduced to bare accessions so the hit table needs "
                 "no re-parsing; nothing else is altered. Queries and reference "
                 "are the same file, and the self hit is excluded by a declared "
                 "rule."),
    }
    (inputs / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"{fasta_name}: {expected} sequences, UniProt {release}")
    print(f"{annotation_name}: {len(annotations)} annotation records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
