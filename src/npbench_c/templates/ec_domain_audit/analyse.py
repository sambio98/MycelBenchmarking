"""Everything the EC domain-audit ladder grades, from the sandbox plus the image.

The one call that is not pure computation is the HMMER run, and it happens here
once per threshold: the baseline plus each declared variant. The declared panel
keeps that to seconds rather than the minutes a full Pfam-A pass costs.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.ec_domain_audit.core import Audit, AuditError


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    audit = Audit.load(campaign, params["sequence_input"])
    rules = audit.rules

    sections = audit.sequences.by_section()
    domains = audit.domains()
    absent = audit.absent(domains)

    # The counterfactual: the same audit under every other declared cutoff. The
    # measured answer here is that it hardly moves, which is worth grading
    # precisely because the opposite was the plausible expectation.
    threshold_absent = {}
    for variant in rules.threshold_variants:
        variant_domains = audit.domains(variant)
        threshold_absent[variant] = len(audit.absent(variant_domains))

    resource = rules.resources["pfam_a_hmm"].verify()
    if not resource["matches"]:
        raise AuditError(
            f"{resource['name']} at {resource['path']} hashes to "
            f"{resource['found_sha256'][:16]} and the campaign declares "
            f"{resource['declared_sha256'][:16]}. A domain verdict against a "
            "different Pfam release is not this campaign's gold.")

    return {
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
        **audit.presence_census(domains),
        "canonical_absent_total": len(absent),
        "family_census": audit.family_census(domains),
        "architecture_census": audit.architecture_census(domains),
        "absent_architecture_census": audit.architecture_census(domains, absent),
        "absent_family_census": audit.family_census(domains, absent),
        "dominant_absent_architecture": next(
            iter(audit.architecture_census(domains, absent))),
        "length_profile": audit.length_profile(domains),
        "misannotated_by_policy": audit.policy_census(domains),
        "threshold_absent": dict(sorted(threshold_absent.items())),
    }
