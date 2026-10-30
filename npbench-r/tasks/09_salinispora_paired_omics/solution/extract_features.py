"""Stage 2 of the reference solution: one filtered MS2 spectrum set per strain.

Reads each strain's LC-MS/MS file from data/spectra/ and writes the processed
MS2 spectra to <work>/features/<accession>.npz. Re-running skips strains
already written, so an interrupted run resumes and a 17.8 GB staging set can
be processed a few files at a time.

    python solution/extract_features.py --task tasks/09_salinispora_paired_omics

The processing is prescribed rather than left open, because every occurrence
call downstream depends on it:

  * intensities normalised to the base peak;
  * peaks within 17 Da of the precursor removed, as in task 01, so the
    unfragmented precursor does not dominate a cosine;
  * the 50 most intense surviving peaks kept -- these are profile-mode spectra
    with a median of 2,686 peaks each, nearly all of them noise, and an
    unfiltered cosine is both slower and worse;
  * at least 6 peaks required, matching the floor task 01 applies to the
    consensus spectra it will be compared against.

Ties in the top-50 cut are broken by m/z ascending, and the kept peaks are
written back in m/z order, so the output does not depend on the sort's
stability.
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata as metadata
from pathlib import Path

import numpy as np
from matchms import Spectrum
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz
from pyteomics import mzxml

PYTEOMICS_VERSION = "5.0.2"

PRECURSOR_WINDOW_DA = 17.0
TOP_PEAKS = 50
MIN_PEAKS = 6
DECIMALS = 4


def check_version() -> None:
    # Read from the distribution metadata: `pyteomics.version` is a submodule
    # that `from pyteomics import mzxml` does not pull in, so the attribute
    # access people reach for first raises AttributeError.
    found = metadata.version("pyteomics")
    if found != PYTEOMICS_VERSION:
        raise SystemExit(
            f"this task's truth was produced by pyteomics {PYTEOMICS_VERSION}, found "
            f"{found}. Scan parsing and precursor reporting are version-dependent."
        )


def process(mz: np.ndarray, intensities: np.ndarray, precursor: float):
    """One MS2 scan, filtered. Returns (mz, intensities) or None if discarded."""
    if not len(mz):
        return None
    order = np.argsort(mz, kind="stable")
    spectrum = Spectrum(
        mz=np.asarray(mz, dtype=float)[order],
        intensities=np.asarray(intensities, dtype=float)[order],
        metadata={"precursor_mz": float(precursor)},
        metadata_harmonization=False,
    )
    spectrum = normalize_intensities(spectrum)
    if spectrum is None:
        return None
    spectrum = remove_peaks_around_precursor_mz(
        spectrum, mz_tolerance=PRECURSOR_WINDOW_DA)
    if spectrum is None or len(spectrum.peaks.mz) < MIN_PEAKS:
        return None

    peak_mz, peak_intensity = spectrum.peaks.mz, spectrum.peaks.intensities
    if len(peak_mz) > TOP_PEAKS:
        # Rank on descending intensity, then ascending m/z, so the cut is
        # defined even when intensities tie after normalisation.
        order = sorted(range(len(peak_mz)),
                       key=lambda i: (-peak_intensity[i], peak_mz[i]))[:TOP_PEAKS]
        keep = np.sort(np.asarray(order))
        peak_mz, peak_intensity = peak_mz[keep], peak_intensity[keep]
    return (np.round(peak_mz, DECIMALS), np.round(peak_intensity, DECIMALS))


def extract(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Filtered MS2 spectra from one mzXML, flattened for compact storage.

    Returns precursor masses, a flat peak array, the offsets that cut it back
    into spectra, the scan numbers, and how many MS2 scans the file held.
    """
    precursors, flat_mz, flat_intensity, offsets, scans = [], [], [], [0], []
    n_ms2 = 0
    for scan in mzxml.read(str(path)):
        if int(scan.get("msLevel", 0)) != 2:
            continue
        n_ms2 += 1
        reported = scan.get("precursorMz")
        if not reported:
            continue
        kept = process(scan.get("m/z array", []), scan.get("intensity array", []),
                       reported[0]["precursorMz"])
        if kept is None:
            continue
        peak_mz, peak_intensity = kept
        precursors.append(round(float(reported[0]["precursorMz"]), DECIMALS))
        scans.append(int(scan.get("num", len(scans) + 1)))
        flat_mz.append(peak_mz)
        flat_intensity.append(peak_intensity)
        offsets.append(offsets[-1] + len(peak_mz))
    return (
        np.asarray(precursors, dtype=float),
        np.concatenate(flat_mz) if flat_mz else np.zeros(0),
        np.concatenate(flat_intensity) if flat_intensity else np.zeros(0),
        np.asarray(offsets, dtype=np.int64),
        np.asarray(scans, dtype=np.int64),
    ), n_ms2


def load_strains(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=None, help="default: <task>/work")
    args = parser.parse_args()
    check_version()

    work = args.work or args.task / "work"
    out_dir = work / "features"
    out_dir.mkdir(parents=True, exist_ok=True)

    strains = load_strains(args.task / "strains.tsv")
    done = skipped = 0
    for index, strain in enumerate(strains, start=1):
        accession = strain["accession"]
        target = out_dir / f"{accession}.npz"
        if target.is_file():
            skipped += 1
            continue
        source = args.task / strain["spectra_file"]
        if not source.is_file():
            print(f"[{index}/{len(strains)}] {accession}: {source} not staged, skipping",
                  flush=True)
            continue
        (precursors, flat_mz, flat_intensity, offsets, scans), n_ms2 = extract(source)
        np.savez_compressed(
            target, precursors=precursors, peak_mz=flat_mz,
            peak_intensity=flat_intensity, offsets=offsets, scans=scans,
            n_ms2_scans=np.asarray([n_ms2]),
        )
        done += 1
        print(f"[{index}/{len(strains)}] {accession}: {n_ms2} MS2 scans, "
              f"{len(precursors)} kept", flush=True)
    print(f"features written for {done} strain(s); {skipped} already present", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
