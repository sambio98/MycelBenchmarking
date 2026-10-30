"""Ship the spectra an MS2 networking campaign reads.

Usage:
  python -m npbench_c.templates.ms2_networking.prepare_input <campaign> \
      --massbank <MassBank-data clone at a release TAG> \
      --contributor AAFC --instrument LC-ESI-ITFT

Run once at instantiation. Everything here is a declared rule over a pinned
release, so the inputs are reproducible from the release and this file alone.

The record cuts are the dereplication template's, imported rather than copied:
MS2, POSITIVE, [M+H]+, an InChIKey, a precursor m/z, at least five peaks. The
parser behind them has been wrong twice -- a `/note` line starting with a slash,
and `ION_MODE` matched inside `FRAGMENTATION_MODE` -- and one copy of it is one
place for that to be fixed.

What this template chooses differently:

**The node set is ONE laboratory's spectra on ONE instrument.** That is what an
experiment is, and it is what a molecular network is built from: a set of spectra
acquired together, most of them unidentified. Mixing instruments inside a network
mixes fragmentation behaviour into the edge weights, so the instrument is part of
the declared node-set rule rather than left open.

**The library is every OTHER contributor's spectra, which is what makes an
annotation cross-laboratory.** Same invariant the dereplication template
enforces, reached more simply: the node contributor is excluded wholesale, so no
node can be annotated from its own lab's replicate.

**The library is filtered to the declared precursor window.** A library spectrum
outside every node's window cannot be a candidate for any node, by the definition
of the filter, so shipping it would add megabytes and change no answer. The rule
is a declared function of the shipped nodes and the declared window -- the agent
can recompute the filter it applies -- and the whole-library census is recorded
in provenance so an auditor can see what was filtered from what.

**A node carries no identity.** The compound a node is, is the answer, so the
node records carry what an experimental spectrum carries and the identities go to
`gold/nodes_truth.json`, which no sandbox copies and no release ships.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import json
import pathlib

from npbench_c.templates.ms2_dereplication.prepare_input import (
    ADDUCT,
    LIBRARY_FIELDS,
    MIN_PEAKS,
    QUERY_FIELDS,
    project,
    read_record,
    release_of,
    write_gz,
)

NODES = "nodes.json.gz"
LIBRARY = "library.json.gz"
PROVENANCE = "provenance.json"
TRUTH = "nodes_truth.json"

#: The precursor window the annotation layer declares. Read from the campaign's
#: pinned rules rather than hard-coded, so one file defines it.
RULES = "networking_rules.json"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"

#: Below this a network is a handful of edges and a census of nothing.
MIN_NODES = 200
#: And a node set that is one compound measured many times is not a network.
MIN_COMPOUNDS = 40
#: The library has to be a reference library rather than one lab's opinion.
MIN_LIBRARY_CONTRIBUTORS = 5


class InputError(RuntimeError):
    """The release cannot be read as the campaign declares it."""


def usable(clone: pathlib.Path) -> tuple[list[dict], dict[str, int]]:
    """Every record that passes the declared cuts, with the census of the cuts."""
    census: collections.Counter = collections.Counter()
    rows: list[dict] = []
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
        rows.append({
            "accession": record["accession"],
            "compound_block": record["inchikey"].split("-")[0],
            "formula": record.get("formula", ""),
            "contributor": record["contributor"],
            "instrument_type": record.get("instrument_type", ""),
            "precursor_mz": precursor,
            "precursor_type": record["precursor_type"],
            "peaks": [[round(mz, 6), round(intensity, 4)]
                      for mz, intensity in sorted(record["peaks"])],
        })
    return rows, dict(census)


def split(rows: list[dict], contributor: str, instrument: str,
          window: float) -> tuple[list[dict], list[dict], dict]:
    """The node set, the library it is annotated from, and their census."""
    nodes = sorted((row for row in rows
                    if row["contributor"] == contributor
                    and row["instrument_type"] == instrument),
                   key=lambda row: row["accession"])
    if len(nodes) < MIN_NODES:
        raise InputError(
            f"{contributor} on {instrument} has {len(nodes)} usable spectra, "
            f"fewer than the {MIN_NODES} a network needs to have a shape")
    blocks = {row["compound_block"] for row in nodes}
    if len(blocks) < MIN_COMPOUNDS:
        raise InputError(
            f"the node set covers {len(blocks)} compounds, fewer than the "
            f"{MIN_COMPOUNDS} this template requires; a set of replicates of a "
            "few compounds is a cluster, not a network")

    elsewhere = sorted((row for row in rows if row["contributor"] != contributor),
                       key=lambda row: row["precursor_mz"])
    masses = [row["precursor_mz"] for row in elsewhere]
    keep: dict[str, dict] = {}
    for node in nodes:
        lo = bisect.bisect_left(masses, node["precursor_mz"] - window)
        hi = bisect.bisect_right(masses, node["precursor_mz"] + window)
        for row in elsewhere[lo:hi]:
            keep[row["accession"]] = row
    library = sorted(keep.values(), key=lambda row: row["accession"])
    if not library:
        raise InputError(
            "no library spectrum falls inside any node's precursor window, so "
            "no node can be annotated and every rung above R1 is empty")
    labs = {row["contributor"] for row in library}
    if len(labs) < MIN_LIBRARY_CONTRIBUTORS:
        raise InputError(
            f"the library draws on {len(labs)} contributors, fewer than the "
            f"{MIN_LIBRARY_CONTRIBUTORS} this template requires")
    if contributor in labs:
        raise InputError(
            f"the library carries {contributor}'s own spectra, so a node could "
            "be annotated from its own laboratory's replicate")

    census = {
        "nodes": len(nodes),
        "node_compounds": len(blocks),
        "node_instrument": instrument,
        "library_spectra": len(library),
        "library_compounds": len({row["compound_block"] for row in library}),
        "library_contributors": len(labs),
        "cross_contributor_spectra_before_the_window_filter": len(elsewhere),
        "nodes_whose_compound_is_in_the_library": sum(
            1 for node in nodes
            if any(row["compound_block"] == node["compound_block"]
                   for row in library)),
    }
    return nodes, library, census


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--massbank", required=True)
    parser.add_argument("--contributor", required=True)
    parser.add_argument("--instrument", required=True)
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    clone = pathlib.Path(args.massbank).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    reference = campaign / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    rules_path = reference / RULES
    if not rules_path.is_file():
        rules_path.write_bytes((CONSTANTS / RULES).read_bytes())
    window = float(json.loads(rules_path.read_text())
                   ["annotation"]["precursor_window"])

    release = release_of(clone)
    rows, cuts = usable(clone)
    nodes, library, census = split(rows, args.contributor, args.instrument, window)

    truth = {row["accession"]: row["compound_block"] for row in nodes}
    gold = campaign / "gold"
    gold.mkdir(parents=True, exist_ok=True)
    (gold / TRUTH).write_text(json.dumps(truth, indent=1, sort_keys=True) + "\n")

    node_digest = write_gz(inputs / NODES, project(nodes, QUERY_FIELDS))
    library_digest = write_gz(inputs / LIBRARY, project(library, LIBRARY_FIELDS))

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
            "string that reads like a release -- 2026.03.1-SNAPSHOT -- which no "
            "clone can reproduce because dev moves. prepare_input refuses one."),
        "node_set_rule": (
            f"Every usable record contributed by {args.contributor} on "
            f"{args.instrument}: one laboratory, one instrument, which is what a "
            "molecular network is built from. Usable means MS2, POSITIVE, "
            f"{ADDUCT}, carrying an InChIKey and a precursor m/z, with at least "
            f"{MIN_PEAKS} peaks."),
        "library_rule": (
            "Every usable record from a DIFFERENT contributor whose precursor "
            f"m/z is within {window} Da -- the declared annotation window -- of "
            "at least one node's. Excluding the node contributor wholesale is "
            "what makes an annotation cross-laboratory. The window filter "
            "changes no answer: a spectrum outside every node's window cannot be "
            "a candidate for any node, by the definition of the filter, so the "
            "filter is a declared function of the shipped nodes and the declared "
            "window rather than a choice about what to hide. The unfiltered "
            "cross-contributor count is in the census below."),
        "cuts": cuts,
        "census": census,
        "files": {NODES: node_digest, LIBRARY: library_digest},
        "peaks": "sorted by m/z, m/z rounded to 6 decimals, intensity to 4",
        "what_is_not_shipped": (
            "No identity for the nodes: the compound a node is, is the answer. A "
            f"node record carries {', '.join(QUERY_FIELDS)} and no compound "
            "block, InChIKey or formula -- a formula would name the compound "
            f"almost uniquely inside a precursor window. The identities go to "
            f"gold/{TRUTH}, which no sandbox copies and no release ships, and "
            "the grading core refuses a node record that carries one. The "
            "library does carry compound blocks, because a reference library "
            "without identities is not a reference library. The network itself "
            "is not shipped either: it follows from the declared edge rule over "
            "the shipped nodes."),
    }, indent=2, sort_keys=True) + "\n")

    print(f"  release {release}")
    print(f"  {cuts['records']} records -> {cuts['usable']} usable")
    print(f"  {census['nodes']} nodes ({census['node_compounds']} compounds) "
          f"from {args.contributor} on {args.instrument}")
    print(f"  {census['library_spectra']} library spectra "
          f"({census['library_compounds']} compounds, "
          f"{census['library_contributors']} contributors) out of "
          f"{census['cross_contributor_spectra_before_the_window_filter']} "
          "cross-contributor spectra")
    print(f"  {census['nodes_whose_compound_is_in_the_library']} nodes have "
          "their own compound somewhere in the library")
    print(f"  identities withheld to gold/{TRUTH} ({len(truth)} entries)")
    print(f"  {NODES} {node_digest[:16]}  {LIBRARY} {library_digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
