"""Ship the spectra an MS2 dereplication campaign reads.

Usage:
  python -m npbench_c.templates.ms2_dereplication.prepare_input <campaign> \
      --massbank <MassBank-data clone at a release TAG> [--queries 24]

Run once at instantiation. Everything here is a declared rule over a pinned
release, so the inputs are reproducible from the release and this file alone.

The rules, and why each one is a rule rather than a preference:

**[M+H]+ positive-mode MS2 only.** Polarity cannot be mixed in a reference
library, and holding the adduct fixed is what makes the R4 misassignment a
perturbation of one thing rather than a comparison across adducts.

**A compound is its InChIKey's first block.** The first block is the skeleton:
two spectra of the same compound from different labs agree on it, while the full
key disagrees on stereochemistry and protonation state that MS2 cannot see.

**Eligible compounds span two instrument types AND two contributors.**
Dereplication is gradeable only when the query compound is in the library
through a DIFFERENT spectrum, and non-trivial only when that spectrum comes from
another instrument in another lab. 1,693 of MassBank 2026.03's compounds qualify
under the clauses above.

**The query is the lexicographically first contributor's spectrum; the library
is everything else.** A declared, order-free split: no sampling, no seed, and no
way for the choice of query to depend on how well it happens to score.

**The release is named by its TAG.** MassBank's dev branch carries a VERSION
string that reads like a release (`2026.03.1-SNAPSHOT`), and an earlier pass
recorded it as one -- a corpus no clone can reproduce, because dev moves. This
script refuses a clone whose VERSION is not a plain release.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import pathlib
import re

QUERIES = "queries.json.gz"
LIBRARY = "library.json.gz"
PROVENANCE = "provenance.json"

#: Matched on the FULL tag, never as a substring: `ION_MODE` inside
#: `FRAGMENTATION_MODE` reported CID and HCD as ion modes the first time this
#: corpus was read.
FIELDS = {
    "accession": re.compile(r"^ACCESSION:\s*(\S+)"),
    "formula": re.compile(r"^CH\$FORMULA:\s*(\S+)"),
    "exact_mass": re.compile(r"^CH\$EXACT_MASS:\s*(\S+)"),
    "inchikey": re.compile(r"^CH\$LINK: INCHIKEY\s+(\S+)"),
    "instrument_type": re.compile(r"^AC\$INSTRUMENT_TYPE:\s*(\S+)"),
    "ms_type": re.compile(r"^AC\$MASS_SPECTROMETRY: MS_TYPE\s+(\S+)"),
    "ion_mode": re.compile(r"^AC\$MASS_SPECTROMETRY: ION_MODE\s+(\S+)"),
    "precursor_mz": re.compile(r"^MS\$FOCUSED_ION: PRECURSOR_M/Z\s+(\S+)"),
    "precursor_type": re.compile(r"^MS\$FOCUSED_ION: PRECURSOR_TYPE\s+(\S+)"),
}
PEAK = re.compile(r"^\s{2}(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*$")

ADDUCT = "[M+H]+"
MIN_PEAKS = 5
VERSION_TAG = re.compile(r"^version\s*=\s*(\d{4}\.\d{2}(?:\.\d+)?)\s*$", re.M)


class InputError(RuntimeError):
    """The release cannot be read as the campaign declares it."""


def read_record(path: pathlib.Path) -> dict:
    record: dict = {}
    peaks: list[tuple[float, float]] = []
    in_peaks = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("PK$PEAK:"):
            in_peaks = True
            continue
        if in_peaks:
            match = PEAK.match(line)
            if match:
                peaks.append((float(match.group(1)), float(match.group(2))))
            continue
        for name, pattern in FIELDS.items():
            if name not in record:
                match = pattern.match(line)
                if match:
                    record[name] = match.group(1)
    record["peaks"] = peaks
    record["contributor"] = path.parent.name
    return record


def release_of(clone: pathlib.Path) -> str:
    version = clone / "VERSION"
    if not version.is_file():
        raise InputError(
            f"{clone} carries no VERSION file; a MassBank clone has one and the "
            "campaign's provenance names the release it was cut from")
    match = VERSION_TAG.search(version.read_text())
    if match is None:
        raise InputError(
            f"{clone}'s VERSION does not name a release: "
            f"{version.read_text().splitlines()[0]!r}")
    release = match.group(1)
    if "SNAPSHOT" in version.read_text():
        raise InputError(
            f"{clone} is a development snapshot ({release}-SNAPSHOT). The dev "
            "branch moves, so a corpus cut from it cannot be reproduced; clone "
            "a release tag instead")
    return release


def corpus(clone: pathlib.Path) -> tuple[dict[str, list[dict]], dict[str, int]]:
    """Every eligible spectrum, keyed by compound, plus the census of the cuts."""
    census: collections.Counter = collections.Counter()
    by_block: dict[str, list[dict]] = {}
    for path in sorted(clone.rglob("*.txt")):
        if path.parent == clone:
            continue
        census["records"] += 1
        record = read_record(path)
        if record.get("ms_type") != "MS2" or record.get("ion_mode") != "POSITIVE":
            continue
        census["ms2_positive"] += 1
        if record.get("precursor_type") != ADDUCT:
            continue
        census["protonated"] += 1
        if not (record.get("inchikey") and record.get("precursor_mz")
                and len(record["peaks"]) >= MIN_PEAKS):
            continue
        try:
            precursor = float(record["precursor_mz"])
        except ValueError:
            continue
        census["usable"] += 1
        block = record["inchikey"].split("-")[0]
        by_block.setdefault(block, []).append({
            "accession": record["accession"],
            "compound_block": block,
            "formula": record.get("formula", ""),
            "contributor": record["contributor"],
            "instrument_type": record.get("instrument_type", ""),
            "precursor_mz": precursor,
            "precursor_type": record["precursor_type"],
            "peaks": [[round(mz, 6), round(intensity, 4)]
                      for mz, intensity in sorted(record["peaks"])],
        })
    census["compounds"] = len(by_block)
    census["compounds_with_2plus_spectra"] = sum(
        1 for rows in by_block.values() if len(rows) >= 2)
    eligible = {
        block: rows for block, rows in by_block.items()
        if len({r["contributor"] for r in rows}) >= 2
        and len({r["instrument_type"] for r in rows}) >= 2}
    census["eligible_compounds"] = len(eligible)
    census["eligible_spectra"] = sum(len(rows) for rows in eligible.values())
    if census["eligible_compounds"] < 500:
        raise InputError(
            f"only {census['eligible_compounds']} compounds are eligible; this "
            "release is not the one the campaign was grounded on")
    return eligible, dict(census)


def split(eligible: dict[str, list[dict]], queries: int) -> tuple[list[dict], list[dict]]:
    chosen_queries, library = [], []
    for block in sorted(eligible):
        rows = sorted(eligible[block], key=lambda r: (r["contributor"], r["accession"]))
        first = rows[0]["contributor"]
        query = next(r for r in rows if r["contributor"] == first)
        chosen_queries.append(query)
        library.extend(r for r in rows if r["accession"] != query["accession"])
    graded = sorted(chosen_queries, key=lambda r: r["accession"])[:queries]
    return graded, sorted(library, key=lambda r: r["accession"])


def write_gz(path: pathlib.Path, payload: object) -> str:
    body = (json.dumps(payload, indent=1, sort_keys=True) + "\n").encode()
    with gzip.GzipFile(path, "wb", compresslevel=9, mtime=0) as handle:
        handle.write(body)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--massbank", required=True)
    parser.add_argument("--queries", type=int, default=24)
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    clone = pathlib.Path(args.massbank).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)

    release = release_of(clone)
    eligible, census = corpus(clone)
    queries, library = split(eligible, args.queries)
    if len(queries) != args.queries:
        raise InputError(
            f"asked for {args.queries} graded queries and the rule yields "
            f"{len(queries)}")

    query_digest = write_gz(inputs / QUERIES, queries)
    library_digest = write_gz(inputs / LIBRARY, library)

    (inputs / PROVENANCE).write_text(json.dumps({
        "source": f"MassBank-data, release tag {release}, cloned shallow from "
                  "https://github.com/MassBank/MassBank-data",
        "url": "https://github.com/MassBank/MassBank-data",
        "license": "CC BY 4.0",
        "attribution": f"MassBank consortium, MassBank-data release {release} "
                       "(https://github.com/MassBank/MassBank-data), CC BY 4.0",
        "release": release,
        "release_is_a_tag": (
            "The tag, not the branch. MassBank's dev branch carries a VERSION "
            "string that reads like a release -- 2026.03.1-SNAPSHOT -- and an "
            "earlier grounding pass recorded that as the corpus, which no clone "
            "can reproduce because dev moves. prepare_input refuses a snapshot."),
        "selection_rule": (
            f"Every record of the release that is MS2, POSITIVE, {ADDUCT}, "
            f"carries an InChIKey and a precursor m/z, and has at least "
            f"{MIN_PEAKS} peaks. Grouped by the InChIKey's first block, a "
            "compound is eligible when its spectra span at least two instrument "
            "types AND at least two contributors. For each eligible compound the "
            "QUERY is the spectrum from the lexicographically first contributor "
            "and every other spectrum of it joins the LIBRARY. The graded "
            f"queries are the first {args.queries} of those by accession; the "
            "library is every eligible compound's non-query spectra, unfiltered, "
            "so the right answer is in it from another lab and so are thousands "
            "of distractors."),
        "census": census,
        "graded_queries": len(queries),
        "library_spectra": len(library),
        "files": {QUERIES: query_digest, LIBRARY: library_digest},
        "peaks": "sorted by m/z, m/z rounded to 6 decimals, intensity to 4",
        "what_is_not_shipped": (
            "No annotation of the queries: the compound each query IS is the "
            "answer. The library carries its own compound blocks, because a "
            "reference library without identities is not a reference library."),
    }, indent=2, sort_keys=True) + "\n")

    print(f"  release {release}")
    print(f"  {census['records']} records -> {census['usable']} usable -> "
          f"{census['eligible_compounds']} eligible compounds")
    print(f"  {len(queries)} graded queries, {len(library)} library spectra")
    print(f"  {QUERIES} {query_digest[:16]}  {LIBRARY} {library_digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
