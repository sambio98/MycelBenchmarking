"""Stage task 09's spectra in batches, fingerprint them, and extract features.

17.8 GB of mzXML will not sit on this container's disk alongside everything
else, and it does not need to: the only things the bundle keeps are the
manifest entry for each file and the compact feature set extracted from it. So
this downloads a batch, records bytes and sha256, runs the solution's own
stage-2 extraction over it, and deletes the mzXML before fetching the next.

Deleting is safe because the manifest is the record: `npbench fetch` restages
any file from its URL and verifies the checksum. An agent attempting the task
stages all 111 at once, which is within the 20 GB target in PROJECT.md
section 10.

    python scripts/stage_task09_spectra.py --task tasks/09_salinispora_paired_omics \
        --batch 8

Re-running resumes: a strain whose feature file already exists is skipped
without being downloaded again.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

CHUNK = 1 << 20


def sha256_and_size(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    total = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
            total += len(block)
    return digest.hexdigest(), total


# The GNPS2 proxy rate-limits: a third consecutive request returns 429. These
# are the delays that got 111 files through without a failure.
PAUSE_BETWEEN = 8.0
RETRY_DELAYS = (15, 30, 60, 120, 240)


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    last = ""
    for attempt, delay in enumerate((0, *RETRY_DELAYS)):
        if delay:
            print(f"    rate limited, waiting {delay}s", flush=True)
            time.sleep(delay)
        result = subprocess.run(
            ["curl", "-sS", "--fail", "--max-time", "1800", "-o", str(partial), url],
            capture_output=True, text=True,
        )
        if result.returncode == 0:
            partial.rename(target)
            return
        partial.unlink(missing_ok=True)
        last = result.stderr.strip()
        # Anything other than a rate limit is not worth retrying.
        if "429" not in last:
            break
    raise SystemExit(f"download failed for {url}: {last}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--batch", type=int, default=8)
    parser.add_argument("--fingerprints", type=Path, default=None,
                        help="default: <task>/work/spectra_fingerprints.json")
    args = parser.parse_args()

    work = args.task / "work"
    work.mkdir(parents=True, exist_ok=True)
    record_path = args.fingerprints or work / "spectra_fingerprints.json"
    records: dict[str, dict] = {}
    if record_path.is_file():
        records = {r["path"]: r for r in json.loads(record_path.read_text())}

    with (args.task / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))

    todo = [
        s for s in strains
        if s["spectra_file"] not in records
        or not (work / "features" / f"{s['accession']}.npz").is_file()
    ]
    print(f"{len(strains)} strains; {len(todo)} still to stage", flush=True)

    for start in range(0, len(todo), args.batch):
        batch = todo[start:start + args.batch]
        staged = []
        for strain in batch:
            target = args.task / strain["spectra_file"]
            if not target.is_file():
                download(strain["spectra_url"], target)
            digest, size = sha256_and_size(target)
            records[strain["spectra_file"]] = {
                "path": strain["spectra_file"], "url": strain["spectra_url"],
                "accession": "MSV000078836", "bytes": size, "sha256": digest,
            }
            staged.append(target)
            print(f"  staged {target.name}: {size / 1e6:.0f} MB", flush=True)
            time.sleep(PAUSE_BETWEEN)

        record_path.write_text(json.dumps(sorted(records.values(),
                                                 key=lambda r: r["path"]), indent=2) + "\n")
        subprocess.run(
            [sys.executable, str(args.task / "solution" / "extract_features.py"),
             "--task", str(args.task)],
            check=True,
        )
        for path in staged:
            path.unlink()
        print(f"  batch done, {start + len(batch)}/{len(todo)} staged and cleared",
              flush=True)

    print(f"fingerprints for {len(records)} file(s) in {record_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
