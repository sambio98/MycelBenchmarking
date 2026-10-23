"""Ship the annotation set and the reference records the audit reads.

Usage:
  python -m npbench_c.templates.ripp_precursor.prepare_input <campaign> \
      --annotations <mibig json tar.gz> \
      --reference-source <dir of antiSMASH-processed .gbk> \
      --release <MIBiG release> --reference-archive-sha256 <sha256> \
      --reference-url <url>

Run once at instantiation.

The reference records shipped are those for EVERY entry of the declared
biosynthetic class that the source set carries -- not the subset the audit ends
up using. That matters: the selection rule's last clause is that an entry has a
reference record at all, and shipping only the entries that pass it would make
the count of those that do not readable off a directory listing instead of
computed from the two resources. The campaign's own readiness gate calls that a
free component.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import shutil
import tarfile
import tempfile

ANNOTATIONS = "mibig_json_4.0.tar.gz"
REFERENCE = "mibig_antismash_ribosomal.tar.gz"
CLASS = "ribosomal"


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _accessions_in_class(archive: pathlib.Path) -> list[str]:
    out = []
    with tarfile.open(archive, "r:*") as tar:
        for member in tar.getmembers():
            if not member.isfile() or not member.name.endswith(".json"):
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            document = json.loads(handle.read().decode())
            classes = (document.get("biosynthesis") or {}).get("classes") or []
            if any(block.get("class") == CLASS for block in classes):
                out.append(document.get("accession")
                           or pathlib.Path(member.name).stem)
    return sorted(out)


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

    annotations = pathlib.Path(args.annotations)
    shutil.copyfile(annotations, inputs / ANNOTATIONS)

    in_class = _accessions_in_class(inputs / ANNOTATIONS)
    if not in_class:
        raise RuntimeError(f"{annotations} holds no {CLASS}-class entries")

    source = pathlib.Path(args.reference_source)
    shipped, absent = [], []
    for accession in in_class:
        (shipped if (source / f"{accession}.gbk").is_file() else absent).append(
            accession)
    if not shipped:
        raise RuntimeError(f"{source} holds no reference record for any entry")
    if not absent:
        raise RuntimeError(
            f"{source} carries a reference record for all {len(in_class)} "
            f"{CLASS} entries, so the selection rule's reference clause cannot "
            "remove anything and the campaign would not exercise it")

    # A reproducible archive: sorted members, and the per-file metadata that
    # would otherwise record this machine's mtimes and ownership stripped out.
    with tempfile.TemporaryDirectory() as tmp:
        staged = pathlib.Path(tmp) / "records"
        staged.mkdir()
        for accession in shipped:
            shutil.copyfile(source / f"{accession}.gbk",
                            staged / f"{accession}.gbk")
        with tarfile.open(inputs / REFERENCE, "w:gz") as tar:
            for accession in shipped:
                member = tar.gettarinfo(str(staged / f"{accession}.gbk"),
                                        arcname=f"{accession}.gbk")
                member.mtime = 0
                member.uid = member.gid = 0
                member.uname = member.gname = ""
                with open(staged / f"{accession}.gbk", "rb") as handle:
                    tar.addfile(member, handle)

    (inputs / "provenance.json").write_text(json.dumps({
        "release": args.release,
        "annotations": {
            "source": "MIBiG",
            "file": ANNOTATIONS,
            "license": "CC BY 4.0",
            "attribution": f"MIBiG {args.release} (https://mibig.secondarymetabolites.org/), licensed CC BY 4.0",
            "sha256": _sha256(inputs / ANNOTATIONS),
            "note": "Shipped whole and unfiltered. The selection rule is part of "
                    "the task, so filtering the set here would answer the first "
                    "rung in the input.",
        },
        "reference_records": {
            "source": "MIBiG antiSMASH-processed reference set",
            "file": REFERENCE,
            "url": args.reference_url,
            "license": "CC BY 4.0",
            "attribution": f"MIBiG antiSMASH-processed reference set ({args.reference_url}), licensed CC BY 4.0",
            "archive_sha256": args.reference_archive_sha256,
            "sha256": _sha256(inputs / REFERENCE),
            "selection_rule": f"every MIBiG entry of biosynthetic class "
                              f"{CLASS!r} for which the source set carries a "
                              f"record, unchanged",
            # The COUNTS are deliberately not written here. The number of
            # entries in the class, the number of records shipped and the gap
            # between them are all graded at R1, and a graded value readable
            # straight out of inputs/ without doing the work is a free
            # component. The rule is published; the arithmetic is the task.
            "note": "The records are the processed GenBank files as published, "
                    "including their CDS translations, which is what the audit "
                    "reads. They are NOT re-run through the image's antiSMASH: "
                    "this campaign never compares output from two antiSMASH "
                    "versions, which antiSMASH's own rules forbid.",
        },
    }, indent=2, sort_keys=True) + "\n")

    print(f"{campaign.name}: {len(in_class)} {CLASS} entries, "
          f"{len(shipped)} reference records shipped, "
          f"{len(absent)} entries without one")
    print(f"  {ANNOTATIONS}: {(inputs / ANNOTATIONS).stat().st_size / 1e6:.2f} MB")
    print(f"  {REFERENCE}: {(inputs / REFERENCE).stat().st_size / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
