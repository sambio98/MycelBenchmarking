"""Deterministic stub systems for the MS2 networking ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each one scores spectra through the tool prefix.

  reader      clears R1: the graph the declared edge rule builds
  annotate    clears R2: what the reference library says about the nodes
  commit      clears R3: what the network claims about the rest
  complete    clears R4
  seedall     COMPETENCE probe, tooled: seeds the network from EVERY node the
              library matched, ignoring the declared confidence tiers. Clears R1
              and fails R2.
  propagate   COMPETENCE probe, tooled: propagates identity across every edge,
              ignoring the precursor difference that separates a same-compound
              claim from an analogue claim. Clears R1 and R2, fails R3.
  noncompute  LEAKAGE ablation: never scores a spectrum and never opens an input.

Two probes, because this campaign has two independent places to be plausibly
wrong, and the measurements put them on different rungs.

**`seedall`.** The declared rule seeds the network only from the confidence tiers
it accepts. Seeding from every match is the obvious shortcut and it is how a
network gets contaminated in practice. Measured on the shipped node set it
doubles the seed count, and the per-component annotation census moves in every
bucket -- the dark components nearly halve and the disagreeing ones multiply --
so R2 fails where the census is. R1 is untouched: the graph does not depend on
the annotation layer at all.

**`propagate`.** This is the mistake molecular networking invites. An edge means
two spectra fragment alike, which at the same precursor supports "the same
compound" and at a different precursor supports "an analogue of it". Treating
every edge as the first is a one-line simplification that reads as diligence.
Measured: every analogue claim on this node set becomes a same-compound claim, so
the claim table, the census and the mass-shift table all fail at R3, and nothing
at R1 or R2 moves. The grounding also measured what the truth says about the two
claims, which is why this is the probe: the analogue claims are borne out and the
same-compound claims mostly are not.

A third plausible slip is recorded rather than fixtured: dropping the declared
matched-peak floor from the edge rule takes the network from 676 edges to 822. It
fails R1, and a fixture that fails the first rung drags the harness's own R1
clear-rate gate below its floor and reports a fixture artefact as a campaign
defect.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.ms2_networking import analyse as net
from npbench_c.templates.ms2_networking.core import Network

LEVELS = ("reader", "annotate", "commit", "complete", "seedall", "propagate",
          "noncompute")
REPORT = "network_report.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    out = pathlib.Path(args.submission)
    out.mkdir(parents=True, exist_ok=True)
    target = out / REPORT

    def write(payload: dict) -> int:
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return 0

    if args.level == "noncompute":
        return write({
            "environment": {}, "network_census": {}, "edges": [],
            "annotation_census": {}, "top_matches": {},
            "component_annotation_census": {}, "propagation_claims": {},
            "propagation_census": {}, "mass_shifts": {}, "threshold_grid": {},
            "dropped_tier": {},
        })

    params = net.load_params(sandbox)
    network = Network.load(sandbox, params["nodes_input"], params["library_input"])

    # R1. The graph, which needs the node spectra and nothing else: this level
    # never opens the library.
    report = {
        "environment": net.environment(network),
        "network_census": net.network_census(network),
        "edges": [[a, b, score, matches]
                  for a, b, score, matches in network.edges()],
    }
    if args.level == "reader":
        return write(report)

    # R2.
    rows = network.annotations()
    report.update(net.annotation_layer(network, rows))
    if args.level == "seedall":
        # Every node the library matched becomes a seed, whatever the declared
        # confidence rule says. The census below is the only thing that changes,
        # and it changes in every bucket.
        seeds = {row["node"]: row["block"] for row in rows
                 if row["verdict"] == "annotated"}
        buckets = {bucket: 0 for bucket in net.COMPONENT_BUCKETS}
        for members in network.components():
            blocks = {seeds[node] for node in members if node in seeds}
            buckets[net.COMPONENT_BUCKETS[min(len(blocks), 3)]] += 1
        report["component_annotation_census"] = buckets
        return write(report)

    report["component_annotation_census"] = net.component_annotation_census(
        network, rows)
    if args.level == "annotate":
        return write(report)

    # R3.
    if args.level == "propagate":
        claims = network.claims(rows)
        census = {label: 0 for label in network.rules.propagation_vocabulary}
        table = {}
        for node, record in sorted(claims.items()):
            claim = ("same_compound"
                     if record["claim"] == "analogue_with_mass_shift"
                     else record["claim"])
            table[node] = claim
            census[claim] += 1
        report["propagation_claims"] = table
        report["propagation_census"] = {label: census[label]
                                        for label in sorted(census)}
        # And no mass shift is reported, because on this reading there are no
        # analogues: the table is empty rather than wrong, which is the same
        # failure stated honestly.
        report["mass_shifts"] = {}
        return write(report)

    report.update(net.propagation(network, rows))
    if args.level == "commit":
        return write(report)

    # R4.
    report["threshold_grid"] = net.threshold_grid(network, rows)
    report["dropped_tier"] = net.dropped_tier(network, rows)
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
