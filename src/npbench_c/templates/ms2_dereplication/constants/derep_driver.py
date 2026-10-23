"""Score MS2 spectra with matchms, as a subprocess of the prefix interpreter.

    python3.11 derep_driver.py <job.json> <out.json>

The job names a peak tolerance and, per query, the library spectra that passed
the campaign's precursor filter. The driver scores those pairs and returns the
ranked candidates. It does the comparison and nothing else: the filter, the
ranking rule and every verdict stay in the campaign's own stdlib code, so what
the tool decides is exactly the similarity and not the answer.

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


def main() -> int:
    job = json.loads(pathlib.Path(sys.argv[1]).read_text())
    tolerance = float(job["peak_tolerance"])
    out: dict[str, list] = {}
    for task in job["tasks"]:
        query = spectrum(task["query"])
        candidates = task["candidates"]
        if not candidates:
            out[task["query"]["accession"]] = []
            continue
        references = [spectrum(row) for row in candidates]
        scores = calculate_scores(
            references=references, queries=[query],
            similarity_function=CosineGreedy(tolerance=tolerance))
        ranked = []
        # Two values per pair, and the first is the REFERENCE SPECTRUM rather
        # than an index: scoring by position here gave 0.007 accuracy where the
        # right answer is 0.80.
        for reference, score in scores.scores_by_query(query, sort=False):
            ranked.append([reference.get("compound_name"),
                           float(score[0]), int(score[1])])
        out[task["query"]["accession"]] = ranked
    pathlib.Path(sys.argv[2]).write_text(
        json.dumps(out, indent=1, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
