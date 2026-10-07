"""Deterministic stub systems for the EC domain-audit ladder.

Harness fixtures. Sandbox-only -- plus the pinned image, which is what makes this
the first S1 campaign: every level above `noncompute` runs HMMER.

  reader        clears R1: the sequence inventory and what the audit ran against
  presence      clears R2: canonical-domain presence per section, and the length
                profile that tests the fragment explanation
  architecture  clears R3: what the canonical-absent sequences carry instead
  complete      clears R4
  inverted      COMPETENCE probe, tooled: runs the tool correctly and then reads
                the hit table the other way round, calling a sequence absent when
                the model matched. Clears R1 and fails R2 -- the tool ran, the
                panel is right, and every verdict is backwards.
  noncompute    LEAKAGE ablation: never runs the tool and never opens the FASTA.

The `inverted` level matters here more than in the stdlib campaigns. A system
with the tool can produce a plausible-looking audit while mishandling the
domtblout columns -- column 1 is the sequence and column 4 the model, and reading
them the other way round reports model names as sequences without failing. This
level is the graded form of that mistake.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.ec_domain_audit import analyse as eca
from npbench_c.templates.ec_domain_audit.core import Audit

LEVELS = ("reader", "presence", "architecture", "complete", "inverted",
          "noncompute")
REPORT = "domain_audit.json"


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
            "release": "2026_03", "ec_number": "1.1.3.15",
            "canonical_model": "FMN_dh", "threshold": "--cut_ga",
            "hmmer_version": "HMMER 3.4", "pfam_version": "Pfam 35.0",
            "pfam_sha256": "0" * 64,
            "sequences_total": 5000, "sequences_by_section": {},
            "model_panel_size": 13,
            "canonical_present_by_section": {},
            "canonical_absent_by_section": {},
            "canonical_absent_total": 1000, "length_profile": {},
            "family_census": {}, "architecture_census": {},
            "absent_family_census": {}, "absent_architecture_census": {},
            "dominant_absent_architecture": "FMN_dh",
            "misannotated_by_policy": {}, "threshold_absent": {},
        })

    params = eca.load_params(sandbox)
    audit = Audit.load(sandbox, params["sequence_input"])
    rules = audit.rules
    domains = audit.domains()
    sections = audit.sequences.by_section()
    resource = rules.resources["pfam_a_hmm"].verify()

    report = {
        "release": audit.sequences.release,
        "ec_number": rules.ec_number,
        "canonical_model": rules.canonical_model,
        "threshold": rules.threshold,
        "hmmer_version": audit.hmmer_version(),
        "pfam_version": resource["version"],
        "pfam_sha256": resource["found_sha256"],
        "sequences_total": len(audit.sequences.records),
        "sequences_by_section": {k: len(v) for k, v in sections.items()},
        "model_panel_size": len(rules.panel),
    }
    if args.level == "reader":
        return write(report)

    if args.level == "inverted":
        # The tool ran and the panel is right; the hit table is read backwards,
        # so presence and absence swap.
        flipped = {a: (frozenset() if rules.canonical_model in m
                       else frozenset({rules.canonical_model}))
                   for a, m in domains.items()}
        report.update(audit.presence_census(flipped))
        report["canonical_absent_total"] = len(audit.absent(flipped))
        report["length_profile"] = audit.length_profile(flipped)
        return write(report)

    absent = audit.absent(domains)
    report.update(audit.presence_census(domains))
    report["canonical_absent_total"] = len(absent)
    report["length_profile"] = audit.length_profile(domains)
    if args.level == "presence":
        return write(report)

    report["family_census"] = audit.family_census(domains)
    report["architecture_census"] = audit.architecture_census(domains)
    report["absent_family_census"] = audit.family_census(domains, absent)
    report["absent_architecture_census"] = audit.architecture_census(domains,
                                                                     absent)
    report["dominant_absent_architecture"] = next(
        iter(report["absent_architecture_census"]))
    if args.level == "architecture":
        return write(report)

    report["misannotated_by_policy"] = audit.policy_census(domains)
    report["threshold_absent"] = {
        variant: len(audit.absent(audit.domains(variant)))
        for variant in sorted(rules.threshold_variants)}
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
