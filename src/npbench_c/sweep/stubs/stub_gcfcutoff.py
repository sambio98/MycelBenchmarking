"""Deterministic stub systems for the GCF-cutoff ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each runs BiG-SCAPE over the declared cutoff grid.

  reader      clears R1: the corpus inventory, the class census and the tool
              identity
  partition   clears R2: the reconstructed partition at the baseline cutoff
  reconcile   clears R3: the comparison with the chemistry-derived families
  complete    clears R4
  tableonly   COMPETENCE probe, tooled: runs the tool correctly and reads the
              clustering table AS the partition, omitting the singleton
              reconstruction. Clears R1 and fails R2 -- every group it reports is
              right, and a third of the corpus is missing.
  noncompute  LEAKAGE ablation: never runs the tool and never opens the corpus.

`tableonly` is the probe this campaign needs, because its symptom is a report that
looks tidier than gold rather than broken: the clustering table lists only
clusters that joined a family, so a system reading it directly reports fewer
groups, no singletons, and a family structure that looks cleaner than the tool
actually produced.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.gcf_cutoff import analyse as gcf
from npbench_c.templates.gcf_cutoff.core import Clustering

LEVELS = ("reader", "partition", "reconcile", "complete", "tableonly",
          "noncompute")
REPORT = "clustering_report.json"


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
            "corpus_size": 30, "designed_families": {}, "controls": [],
            "bigscape_version": "BiG-SCAPE 2.0.3", "pfam_version": "Pfam 35.0",
            "pfam_sha256": "0" * 64,
            "cutoffs": [], "baseline_cutoff": "0.3", "include_gbk": "BGC",
            "class_census": {},
            "baseline_partition": [], "groups": 15, "size_census": {},
            "largest_group": 5, "singletons": 10, "clustered_records": 20,
            "per_family": {}, "verdict_census": {}, "families_exact": 1,
            "controls_absorbed": [], "designed_pairs": 80,
            "co_clustered_pairs": 40, "recovered_pairs": 35, "missed_pairs": 45,
            "false_join_pairs": 5, "pair_jaccard": 0.4,
            "cutoff_sweep": {}, "best_cutoff": "0.5", "best_pair_jaccard": 0.8,
            "first_cutoff_with_a_false_join": "0.6",
            "first_cutoff_absorbing_a_control": "0.6",
            "first_cutoff_recovering_each_family": {},
        })

    params = gcf.load_params(sandbox)
    clustering = Clustering.load(sandbox, params["cluster_dir"],
                                 params["designed_families_input"])
    rules = clustering.rules
    pfam = rules.resources["pfam_a_hmm"]
    from npbench_c.tools.cache import file_digest

    assignments = clustering.assignments()
    baseline_label = rules.label(rules.baseline_cutoff)
    assignment = assignments[baseline_label]

    report = {
        "corpus_size": len(clustering.corpus),
        "designed_families": {name: len(members) for name, members
                              in sorted(clustering.designed.items())},
        "controls": list(clustering.controls),
        "bigscape_version": clustering.bigscape_version(),
        "pfam_version": pfam.version,
        "pfam_sha256": file_digest(pfam.resolve()),
        "cutoffs": [rules.label(c) for c in rules.cutoffs],
        "baseline_cutoff": baseline_label,
        "include_gbk": rules.include_gbk,
        "class_census": clustering.class_census(assignment),
    }
    if args.level == "reader":
        return write(report)

    if args.level == "tableonly":
        # The tool ran and every group is right; the records the clustering table
        # omits are dropped instead of being reconstructed as singletons.
        groups: dict[tuple[str, str], list[str]] = {}
        for stem, key in sorted(assignment.items()):
            groups.setdefault(key, []).append(stem)
        listed = sorted(sorted(members) for members in groups.values())
        report["baseline_partition"] = listed
        report.update(clustering.group_census(listed))
        return write(report)

    partition = clustering.partition(assignment)
    report["baseline_partition"] = partition
    report.update(clustering.group_census(partition))
    if args.level == "partition":
        return write(report)

    report.update(clustering.reconciliation(partition))
    report.update(clustering.pair_statistics(partition))
    if args.level == "reconcile":
        return write(report)

    sweep = clustering.sweep(assignments)
    report["cutoff_sweep"] = sweep
    report.update(clustering.sweep_summary(sweep))
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
