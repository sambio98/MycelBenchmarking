"""Score an MS2 fixture with matchms CosineGreedy, deterministically.

Run by the PREFIX interpreter as a subprocess -- never imported into the
benchmark's runtime, which is stdlib-only by design.

    python3.11 matchms_cosine.py --fixture <json[.gz]> --tolerance 0.01 \
        --out <tsv> [--shuffle-library SEED]

The artifact is a TSV of every query's top match: the query accession, the
library accession, the score at full float precision, and the number of matched
peaks. Full precision on purpose -- a rounded score hides exactly the bit-level
difference an invariance suite exists to catch.
"""
import argparse, gzip, json, pathlib, random, sys

from matchms import Spectrum, calculate_scores
from matchms.similarity import CosineGreedy
import numpy as np


def load(path: pathlib.Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def spectra(rows: list[dict]) -> list[Spectrum]:
    out = []
    for row in rows:
        peaks = sorted(row["peaks"], key=lambda p: p[0])
        out.append(Spectrum(
            mz=np.array([p[0] for p in peaks], dtype="float64"),
            intensities=np.array([p[1] for p in peaks], dtype="float64"),
            metadata={"compound_name": row["accession"],
                      "precursor_mz": row["precursor_mz"],
                      "inchikey_block": row["block"]},
            metadata_harmonization=False)) 
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixture", required=True)
    ap.add_argument("--tolerance", type=float, default=0.01)
    ap.add_argument("--out", required=True)
    ap.add_argument("--shuffle-library", type=int, default=None)
    args = ap.parse_args()

    payload = load(pathlib.Path(args.fixture))
    queries = spectra(payload["queries"])
    library_rows = list(payload["library"])
    if args.shuffle_library is not None:
        random.Random(args.shuffle_library).shuffle(library_rows)
    library = spectra(library_rows)

    scores = calculate_scores(references=library, queries=queries,
                              similarity_function=CosineGreedy(tolerance=args.tolerance))
    lines = []
    for query in queries:
        best = None
        # Two values, not three, and the reference is the SPECTRUM OBJECT rather
        # than an index: an earlier pass indexed the library positionally here
        # and scored 0.007 where the correct answer is 0.80.
        for reference, score in scores.scores_by_query(query, sort=False):
            value, matched = float(score[0]), int(score[1])
            key = (value, matched, reference.get("compound_name"))
            if best is None or (value, matched) > (best[0], best[1]) or (
                    (value, matched) == (best[0], best[1])
                    and reference.get("compound_name") < best[2]):
                best = key
        lines.append("\t".join([
            query.get("compound_name"), best[2], repr(best[0]), str(best[1])]))
    pathlib.Path(args.out).write_text("\n".join(sorted(lines)) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
