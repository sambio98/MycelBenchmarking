"""Deterministic stub systems for the annotation-audit ladder.

Harness fixtures. Sandbox-only. Following the distinction the MIBiG diff
campaign forced:

  stub-census      clears R1-R2: the backed/unbacked census is right
  stub-analyst     clears R3: applies the three-condition admissibility policy
  stub-complete    clears R4
  stub-conflated   a COMPETENCE probe, tooled: reports the merely-backed counts
                   as admissible, conflating the function evidence axis with the
                   locus and status conditions. Clears R2, fails R3.
  stub-noncompute  the LEAKAGE ablation: never opens the archive.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import yaml

from npbench_c.templates.mibig_annotation_audit.core import AnnotationAudit, Policy
from npbench_c.templates.mibig_diff.core import load_release

LEVELS = ("naive", "census", "analyst", "complete", "conflated", "noncompute")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    audit = AnnotationAudit.load(sandbox)
    params = yaml.safe_load((sandbox / "task.yaml").read_text())["template_params"]
    functions = sorted(audit.function_vocabulary)

    sub = pathlib.Path(args.submission)
    sub.mkdir(parents=True, exist_ok=True)
    report: dict = {"function_vocabulary": functions, "policy_variants": {}}

    if args.level == "noncompute":
        report.update({"corpus_entries": 3000, "entries_with_annotations": 600,
                       "total_function_annotations": 2500,
                       "backed_by_function": {f: 100 for f in functions},
                       "unbacked_by_function": {f: 100 for f in functions},
                       "function_evidence_method_totals": {},
                       "admissible_by_function": {f: 100 for f in functions},
                       "backed_minus_admissible": {f: 0 for f in functions},
                       "focal_admissible_entries": []})
        (sub / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    corpus = load_release(sandbox / params["corpus_archive"])
    report["corpus_entries"] = len(corpus)
    census = audit.census(corpus)

    # R1 asks for the corpus and annotation counts, which is counting rather
    # than classification, so the naive level supplies them and stops short of
    # the backed/unbacked split. One stub per rung boundary.
    report.update({k: census[k] for k in (
        "entries_with_annotations", "total_function_annotations")})

    if args.level == "naive":
        (sub / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    report.update({k: census[k] for k in (
        "backed_by_function", "unbacked_by_function",
        "function_evidence_method_totals")})

    if args.level == "conflated":
        # Conflates the three evidence axes: reports merely-backed as admissible.
        report["admissible_by_function"] = census["backed_by_function"]
        report["backed_minus_admissible"] = {f: 0 for f in functions}
        report["focal_admissible_entries"] = []
        (sub / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if args.level == "census":
        (sub / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    baseline = audit.admissible(corpus, audit.baseline)
    report["admissible_by_function"] = baseline["admissible_by_function"]
    report["backed_minus_admissible"] = {
        f: census["backed_by_function"][f] - baseline["admissible_by_function"][f]
        for f in functions}
    report["focal_admissible_entries"] = baseline[
        "admissible_entries_by_function"][params["focal_function"]]

    if args.level == "analyst":
        (sub / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        return 0

    for name, spec in sorted(params["policy_variants"].items()):
        policy = Policy(name=name,
                        function_evidence=frozenset(spec["function_evidence_accepted"]),
                        locus_evidence=frozenset(spec["locus_evidence_accepted"]),
                        require_active=bool(spec.get("require_active_status", True)))
        result = audit.admissible(corpus, policy)
        report["policy_variants"][name] = {
            "admissible_by_function": result["admissible_by_function"],
            "total": sum(result["admissible_by_function"].values())}
    (sub / "audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
