"""Deterministic stub systems for the MIBiG release-diff ladder.

Harness fixtures, not agents. They read ONLY the sandbox -- the two release
archives in inputs/ and the mapping tables in reference/ -- using the template
core as library code, the way an agent would use any tool.

Two of these are worth distinguishing, because conflating them makes the
leakage gate ask the wrong question:

  stub-nonormalise  a COMPETENCE probe. It loads both corpora and computes the
                    raw diff correctly, then reports it as the classified
                    result, so it clears R2 and fails R3. It is a tooled system
                    missing one skill, not an ablation -- it did the work.
  stub-noncompute   the LEAKAGE ablation. It never opens the archives and
                    answers with plausible round numbers, which is what
                    "answerable without doing the work" actually means for an S0
                    campaign where no tool can be withheld.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import yaml

from npbench_c.templates.mibig_diff.core import ReleaseDiff, load_release

LEVELS = ("naive", "rawdiff", "analyst", "complete", "nonormalise", "noncompute")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    diff = ReleaseDiff.load(sandbox)
    params = yaml.safe_load((sandbox / "task.yaml").read_text())["template_params"]
    sub = pathlib.Path(args.submission)
    sub.mkdir(parents=True, exist_ok=True)

    if args.level == "noncompute":
        # Never opens the archives. Reports the sort of round numbers a system
        # answering from priors would offer, so it fails at R1.
        (sub / "diff_report.json").write_text(json.dumps({
            "source_entries": 2500, "target_entries": 3000,
            "fields": sorted(diff.fields),
            "only_in_source": 0, "only_in_target": 500, "shared_entries": 2500,
            "shared_status_in_target": {"active": 2500},
            "raw_difference_counts": {f: 0 for f in sorted(diff.fields)},
            "classified_verdict_counts": {f: {"unchanged": 2500}
                                          for f in sorted(diff.fields)},
            "retired_class_terms": [], "perturbations": {},
        }, indent=2, sort_keys=True))
        return 0

    source = load_release(sandbox / params["source_archive"])
    target = load_release(sandbox / params["target_archive"])
    fields = sorted(diff.fields)
    report: dict = {
        "source_entries": len(source),
        "target_entries": len(target),
        "fields": fields,
        "perturbations": {},
    }

    if args.level == "naive":
        (sub / "diff_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    corpus = diff.corpus_diff(source, target)
    for key in ("only_in_source", "only_in_target", "shared_entries",
                "shared_status_in_target", "raw_difference_counts"):
        report[key] = corpus[key]

    if args.level == "nonormalise":
        # The no-normalisation ablation: the raw counts are right, and are then
        # reported as though they were the classified result.
        report["classified_verdict_counts"] = {
            f: {"substantive_change": corpus["raw_difference_counts"][f]}
            for f in fields if corpus["raw_difference_counts"][f]
        }
        report["retired_class_terms"] = []
        (sub / "diff_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if args.level == "rawdiff":
        (sub / "diff_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    report["classified_verdict_counts"] = corpus["classified_verdict_counts"]
    report["retired_class_terms"] = sorted(diff.retired_terms)

    if args.level == "analyst":
        (sub / "diff_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    from npbench_c.templates.mibig_diff.build import _apply

    base = params["perturbation_base"]
    for name, operation in sorted(params["perturbations"].items()):
        mutated = _apply(source[base], operation["operation"])
        report["perturbations"][name] = {
            "verdict": diff.classify_field(operation["field"], mutated, target[base])}
    (sub / "diff_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
