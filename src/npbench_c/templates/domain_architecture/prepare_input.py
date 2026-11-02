"""Ship the annotation set and the translations the architecture audit reads.

Usage:
  python -m npbench_c.templates.domain_architecture.prepare_input <campaign> \
      --annotations <mibig json tar.gz> \
      --reference-source <dir of antiSMASH-processed .gbk> \
      --release <MIBiG release> --reference-archive-sha256 <sha256> \
      --reference-url <url>

Run once at instantiation.

The translations shipped are those for EVERY gene a curated module names, on
every entry the extraction can resolve -- not the subset one instantiation's
module type selects. That matters twice. It keeps the module-type selection in
the agent's hands, which is where the first rung's census comes from; and it
means the NRPS and PKS instantiations of this template ship the same file, so
neither reveals the other's corpus by its size.

The GenBank records themselves are not shipped. Reading several hundred processed
cluster records to pull out a thousand translations is forty megabytes to deliver
five, and the extraction is a mechanical step with no judgement in it -- so the
rule is published here, the source archive is fingerprinted, and the FASTA is
what travels.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))

from npbench_c.templates.ripp_precursor.core import cds_translations

ANNOTATIONS = "mibig_json_4.0.tar.gz"
PROTEINS = "module_gene_proteins.faa.gz"


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--reference-source", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--reference-archive-sha256", required=True)
    parser.add_argument("--reference-url", required=True)
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(pathlib.Path(args.annotations), inputs / ANNOTATIONS)

    import tarfile
    wanted: dict[str, set[str]] = {}
    with tarfile.open(inputs / ANNOTATIONS, "r:*") as archive:
        for member in sorted(archive.getmembers(), key=lambda m: m.name):
            if not member.isfile() or not member.name.endswith(".json"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            document = json.loads(handle.read().decode())
            accession = document.get("accession") or pathlib.Path(member.name).stem
            for module in (document.get("biosynthesis") or {}).get("modules") or []:
                for gene in module.get("genes") or []:
                    wanted.setdefault(accession, set()).add(gene)
    if not wanted:
        raise RuntimeError(f"{args.annotations} holds no curated modules")

    source = pathlib.Path(args.reference_source)
    records = []
    missing_record = missing_gene = 0
    for accession, genes in sorted(wanted.items()):
        gbk = source / f"{accession}.gbk"
        if not gbk.is_file():
            missing_record += len(genes)
            continue
        translations = cds_translations(gbk.read_text())
        for gene in sorted(genes):
            sequence = translations.get(gene)
            if sequence is None:
                missing_gene += 1
                continue
            records.append((accession, gene, sequence))
    if not records:
        raise RuntimeError(f"{source} resolved no module gene to a translation")

    text = "".join(f">{a}|{g}\n{s}\n" for a, g, s in sorted(records))
    # mtime=0 so the archive is byte-reproducible across instantiations.
    with gzip.GzipFile(inputs / PROTEINS, "wb", mtime=0) as handle:
        handle.write(text.encode())

    (inputs / "provenance.json").write_text(json.dumps({
        "release": args.release,
        "annotations": {
            "source": "MIBiG",
            "file": ANNOTATIONS,
            "license": "CC BY 4.0",
            "attribution": f"MIBiG {args.release} (https://mibig.secondarymetabolites.org/), licensed CC BY 4.0",
            "sha256": _sha256(inputs / ANNOTATIONS),
            "note": "Shipped whole and unfiltered. The module-type selection is "
                    "part of the task, so filtering the set here would answer "
                    "the first rung in the input.",
        },
        "proteins": {
            "source": "MIBiG antiSMASH-processed reference set",
            "file": PROTEINS,
            "url": args.reference_url,
            "license": "CC BY 4.0",
            "attribution": f"MIBiG antiSMASH-processed reference set ({args.reference_url}), licensed CC BY 4.0",
            "archive_sha256": args.reference_archive_sha256,
            "sha256": _sha256(inputs / PROTEINS),
            "extraction_rule": "the CDS translation of every gene named by any "
                               "curated module in the annotation set, taken "
                               "unchanged from that cluster's processed record "
                               "and named <cluster>|<gene>; no module type, "
                               "status or domain content enters the rule",
            # The COUNTS are deliberately not written here. The corpus size and
            # every drop the selection makes are graded at R1, and a graded value
            # readable straight out of inputs/ without doing the work is a free
            # component.
            "note": "Translations only. The processed GenBank records are forty "
                    "megabytes to deliver five, and pulling a CDS translation "
                    "out of one is mechanical, so the rule is published and the "
                    "FASTA travels. The records are NOT re-run through the "
                    "image's antiSMASH: this campaign never compares output "
                    "from two antiSMASH versions.",
        },
    }, indent=2, sort_keys=True) + "\n")

    print(f"{campaign.name}: {sum(len(g) for g in wanted.values())} genes named "
          f"by a curated module over {len(wanted)} entries")
    print(f"  resolved {len(records)}; {missing_record} had no reference record, "
          f"{missing_gene} no matching CDS")
    print(f"  {ANNOTATIONS}: {(inputs / ANNOTATIONS).stat().st_size / 1e6:.2f} MB")
    print(f"  {PROTEINS}: {(inputs / PROTEINS).stat().st_size / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
