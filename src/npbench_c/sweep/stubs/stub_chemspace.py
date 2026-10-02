"""Deterministic stub systems for the chemical-space ladder.

Harness fixtures, never a published result: a stub-based sweep marks itself and
the readiness gate keeps every agent-run check PENDING.

  census       clears R1: the corpus and the three distinct-key counts
  granularity  clears R2: what each key merges, and the per-class spread
  informative  clears R3: the sharing claim under the declared filter
  complete     clears R4
  unfiltered   COMPETENCE probe, tooled: does all the bookkeeping but reports
               the UNFILTERED sharing statistic as the informative one, which is
               the mistake the rung exists to catch -- benzene in 76 entries read
               as a shared scaffold. Clears R2, fails R3.
  noncompute   LEAKAGE ablation: never opens the table or the archive.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.chemical_space import analyse as chem
from npbench_c.templates.chemical_space.core import ChemicalSpace

LEVELS = ("census", "granularity", "informative", "complete", "unfiltered",
          "noncompute")
REPORT = "chemical_space_report.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    out = pathlib.Path(args.submission)
    out.mkdir(parents=True, exist_ok=True)
    target = out / REPORT

    if args.level == "noncompute":
        # Shaped like a report, informed by nothing. The record and distinct-key
        # counts cannot be guessed, so this fails R1 and sets the leakage floor.
        target.write_text(json.dumps({
            "corpus_release": "4.0", "corpus_entries": 3000,
            "active_entries": 2400, "compound_records": 3000,
            "entries_with_compounds": 1900, "distinct_compound_keys": 3000,
            "distinct_connectivity_keys": 3000, "distinct_scaffold_keys": 1500,
            "acyclic_compounds": 0, "duplicate_records": 0,
            "class_vocabulary": ["NRPS", "PKS", "other", "ribosomal",
                                 "saccharide", "terpene"],
            "entries_per_class": {}, "multi_class_entries": 0,
            "cross_entry_duplicates": {}, "informative_rule": {},
            "shared_scaffolds": {}, "cross_class_scaffolds": {},
            "cross_class_informative_entries": {},
            "most_shared_informative_scaffold": {},
            "threshold_variants": {}, "evidence_gate": {},
        }, indent=2, sort_keys=True) + "\n")
        return 0

    params = chem.load_params(sandbox)
    space = ChemicalSpace.load(sandbox, params["corpus_archive"])

    report = chem.census(space)
    if args.level == "census":
        target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        return 0

    report.update(chem.granularity(space))
    if args.level == "granularity":
        target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        return 0

    commitment = chem.commitment(space)
    if args.level == "unfiltered":
        # Reports the unfiltered statistic in the informative slots: no filter
        # applied, so benzene and the single heterocycles count as shared
        # scaffolds and as crossing classes.
        unfiltered = chem.sharing(space, space.records, min_rings=1, min_coverage=0.0)
        entries = unfiltered.pop("_cross_class_informative_entries")
        unfiltered.pop("_informative")
        report.update({
            "informative_rule": commitment["informative_rule"],
            "informative_records": unfiltered["informative_records"],
            "informative_scaffold_keys": unfiltered["informative_scaffold_keys"],
            "shared_scaffolds": unfiltered["shared_scaffolds"],
            "cross_class_scaffolds": unfiltered["cross_class_scaffolds"],
            "cross_class_informative_entries": entries,
            "most_shared_informative_scaffold":
                commitment["most_shared_informative_scaffold"],
        })
        target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        return 0

    report.update(commitment)
    if args.level == "informative":
        target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        return 0

    report.update(chem.counterfactuals(space, params))
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
