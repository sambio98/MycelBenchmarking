"""Stage 1 of the reference solution: run GECCO over every staged genome.

Writes one GECCO output directory per assembly under <work>/gecco/<accession>/.
Re-running skips assemblies whose clusters.tsv already exists, so an
interrupted run resumes.

    python solution/run_gecco.py --task tasks/02_salinispora_bgc --jobs 4
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

GECCO_VERSION = "0.11.0"


def check_version() -> None:
    out = subprocess.run(["gecco", "--version"], capture_output=True, text=True).stdout
    if GECCO_VERSION not in out:
        raise SystemExit(
            f"this task's truth was produced by GECCO {GECCO_VERSION}, found: {out.strip()!r}. "
            f"A cluster inventory is only valid against the exact version that produced it."
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=None, help="default: <task>/work")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()

    check_version()
    work = args.work or args.task / "work"
    genomes = sorted((args.task / "data" / "genomes").glob("*.fna.gz"))
    if not genomes:
        raise SystemExit(f"no genomes under {args.task / 'data' / 'genomes'}; run `make fetch` first")

    print(f"{len(genomes)} genomes -> {work / 'gecco'}", flush=True)
    started = time.time()
    for index, archive in enumerate(genomes, start=1):
        accession = archive.name.removesuffix(".fna.gz")
        out = work / "gecco" / accession
        if (out / f"{accession}.clusters.tsv").is_file():
            continue
        out.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as scratch:
            plain = Path(scratch) / f"{accession}.fna"
            with gzip.open(archive, "rb") as src, plain.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            result = subprocess.run(
                ["gecco", "--quiet", "--no-color", "--no-markup", "--jobs", str(args.jobs),
                 "run", "--genome", str(plain), "--output-dir", str(out)],
                capture_output=True, text=True,
            )
        if result.returncode != 0:
            print(f"!! {accession} failed:\n{result.stderr[-2000:]}", flush=True)
            return 1
        elapsed = time.time() - started
        print(f"  {index}/{len(genomes)} {accession} ({elapsed/index:.0f}s/genome avg)", flush=True)
    print(f"done in {(time.time() - started)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
