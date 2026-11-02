"""Ship the contig the gene-calling campaign reads.

Usage:
  python -m npbench_c.templates.gene_calling.prepare_input <campaign> \
      --genbank <annotated record .gb or .gb.gz> --accession <id> \
      --source <description> --license <licence> --url <url>

Run once at instantiation.

The input is a FASTA of the sequence and nothing else. The source record is
annotated -- it is where the accession's curated CDS features live -- and those
features are what a gene caller predicts, so shipping it would hand over the
answer. This script reads the annotated record, extracts the bases, and writes
the sequence alone.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pathlib
import re

CONTIG = "contig.fna.gz"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--genbank", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--license", required=True)
    parser.add_argument("--url", required=True)
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)

    source = pathlib.Path(args.genbank)
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt") as handle:
        text = handle.read()
    if "\nORIGIN" not in text:
        raise RuntimeError(f"{source} carries no ORIGIN block to read bases from")
    body = text[text.index("\nORIGIN") + 7:].split("//")[0]
    sequence = "".join(re.findall(r"[acgtnACGTN]+", body)).upper()
    if not sequence:
        raise RuntimeError(f"{source}'s ORIGIN block holds no sequence")

    # How many CDS features the source carries, for the record. It is NOT shipped
    # and NOT graded: it is the answer to the question the campaign asks, and it
    # is written here only so the provenance says plainly what was discarded.
    curated = len(re.findall(r"\n     CDS {13}", text))

    fasta = (f">{args.accession}\n"
             + "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70))
             + "\n").encode()
    # mtime=0 so the archive is byte-reproducible across instantiations.
    with gzip.GzipFile(inputs / CONTIG, "wb", mtime=0) as handle:
        handle.write(fasta)

    (inputs / "provenance.json").write_text(json.dumps({
        "accession": args.accession,
        "source": args.source,
        "url": args.url,
        "license": args.license,
        "file": CONTIG,
        # The archive's own digest, for provenance. The SEQUENCE's length and
        # digest are deliberately absent: both are graded at R1, and a graded
        # value readable straight out of inputs/ without doing the work is a
        # free component. The digest is over the bases rather than the file, so
        # computing it is also where the declared normalisation gets checked.
        "sha256": _sha256_bytes((inputs / CONTIG).read_bytes()),
        "extraction_rule": "the bases of the source record's ORIGIN block, "
                           "upper-cased and wrapped at 70 columns, under the "
                           "accession as the single FASTA identifier",
        "annotation_discarded": (
            f"the source record carries {curated} CDS features and none of them "
            "is shipped. A gene caller predicts coding sequences, so an "
            "annotated record is the answer to this campaign's first rung. The "
            "count is recorded here for the provenance and is not graded."),
        "note": "Sequence only. The detection campaign on this same accession "
                "ships the annotated GenBank, because there the annotation is "
                "what the tool reads; here it is what the tool is asked to "
                "reproduce.",
    }, indent=2, sort_keys=True) + "\n")

    print(f"{campaign.name}: {args.accession}, {len(sequence)} bases, "
          f"{curated} CDS features discarded")  # printed, not written
    print(f"  {CONTIG}: {(inputs / CONTIG).stat().st_size / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
