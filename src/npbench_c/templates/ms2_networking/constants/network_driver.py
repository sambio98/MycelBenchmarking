"""Score MS2 spectra with matchms, as a subprocess of the prefix interpreter.

    python3.11 network_driver.py <job.json> <out.json>

Two jobs in one driver, because a networking campaign needs both and they are
the same call with different arguments:

  "pairs"   every pair within one set of spectra, for the network itself
  "windows" one set of library candidates per node, for the annotation layer

Never imported into the campaign runtime. matchms brings numpy, scipy and numba,
and the benchmark's grading core is stdlib-only by design.
"""
import json
import pathlib
import sys

from matchms import Spectrum, calculate_scores
from matchms.similarity import CosineGreedy
import numpy as np


def spectrum(row: dict) -> Spectrum:
    peaks = row["peaks"]
    return Spectrum(
        mz=np.array([p[0] for p in peaks], dtype="float64"),
        intensities=np.array([p[1] for p in peaks], dtype="float64"),
        metadata={"compound_name": row["accession"],
                  "precursor_mz": row["precursor_mz"]},
        # Declared: matchms harmonises metadata by default, and a default that
        # runs unrecorded is a parameter of gold nobody declared.
        metadata_harmonization=False)


def pairs(job: dict) -> list:
    """Every pair in one set, upper triangle only, as [a, b, score, matches].

    `is_symmetric=True` is the whole reason this is one call: the matrix is
    computed once and read twice.
    """
    rows = job["spectra"]
    specs = [spectrum(row) for row in rows]
    names = [row["accession"] for row in rows]
    scores = calculate_scores(
        references=specs, queries=specs,
        similarity_function=CosineGreedy(tolerance=float(job["peak_tolerance"])),
        is_symmetric=True)
    matrix = scores.to_array()
    floor = float(job.get("score_floor", 0.0))
    out = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            score = float(matrix[i][j]["CosineGreedy_score"])
            if score < floor:
                continue
            out.append([names[i], names[j], score,
                        int(matrix[i][j]["CosineGreedy_matches"])])
    return out


def windows(job: dict) -> dict:
    """Per node, its candidate library spectra scored. Same shape as the
    dereplication driver's output, and for the same reason: the filter and the
    ranking stay in the campaign's own code."""
    out: dict[str, list] = {}
    tolerance = float(job["peak_tolerance"])
    for task in job["tasks"]:
        query = spectrum(task["query"])
        candidates = task["candidates"]
        if not candidates:
            out[task["query"]["accession"]] = []
            continue
        scores = calculate_scores(
            references=[spectrum(row) for row in candidates], queries=[query],
            similarity_function=CosineGreedy(tolerance=tolerance))
        ranked = []
        # Two values per pair, and the first is the REFERENCE SPECTRUM rather
        # than an index: scoring by position gave 0.007 accuracy on the
        # dereplication template where the right answer was 0.80.
        for reference, score in scores.scores_by_query(query, sort=False):
            ranked.append([reference.get("compound_name"),
                           float(score[0]), int(score[1])])
        out[task["query"]["accession"]] = ranked
    return out


def main() -> int:
    job = json.loads(pathlib.Path(sys.argv[1]).read_text())
    kind = job["kind"]
    if kind == "pairs":
        result = pairs(job)
    elif kind == "windows":
        result = windows(job)
    else:
        print(f"unknown job kind {kind!r}", file=sys.stderr)
        return 2
    pathlib.Path(sys.argv[2]).write_text(
        json.dumps(result, indent=1, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
