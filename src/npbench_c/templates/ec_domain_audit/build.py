"""Build gold, warm bundles and the grading spec for the EC domain-audit ladder.

Usage: python -m npbench_c.templates.ec_domain_audit.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.ec_domain_audit.analyse import analyse, load_params
from npbench_c.templates.ec_domain_audit.core import (
    CONSTANTS,
    RULES,
    AuditError,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.002, "r3": 0.001, "r4": 0.001}

INVENTORY_KEYS = ("release", "ec_number", "canonical_model", "threshold",
                  "hmmer_version", "pfam_version", "pfam_sha256",
                  "sequences_total", "sequences_by_section", "model_panel_size")
PRESENCE_KEYS = ("canonical_present_by_section", "canonical_absent_by_section",
                 "canonical_absent_total", "length_profile")
ARCHITECTURE_KEYS = ("family_census", "architecture_census",
                     "absent_family_census", "absent_architecture_census",
                     "dominant_absent_architecture")


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    absent = gold["canonical_absent_total"]
    if absent < 100:
        raise AuditError(
            f"only {absent} sequences lack the canonical domain; the audit would "
            "be an anecdote rather than a measurement")
    if absent == gold["sequences_total"]:
        raise AuditError(
            "no sequence carries the canonical domain, so either the panel or "
            "the threshold is wrong -- not an annotation finding")

    present = gold["canonical_present_by_section"]
    if not present.get("reviewed"):
        raise AuditError(
            "no reviewed entry carries the canonical domain, so the control "
            "side of the comparison is empty and the unreviewed rate means "
            "nothing on its own")

    architectures = gold["absent_architecture_census"]
    if len(architectures) < 2:
        raise AuditError(
            "the canonical-absent set has one architecture, so R3's census is a "
            "single number dressed up as a distribution")

    policies = gold["misannotated_by_policy"]
    if len(set(policies.values())) < 2:
        raise AuditError(
            f"every declared policy gives the same count {policies}; R4's "
            "policy axis does not move and the rung would show nothing")

    # The threshold axis is allowed to be inert -- that is a measured robustness
    # result -- but only because the policy axis is not. If both go flat the rung
    # is vacuous, which is the defect that re-scoped the selectivity campaign.
    spread = max(gold["threshold_absent"].values()) - min(
        gold["threshold_absent"].values())
    strict = policies.get("p0_canonical_only")
    widest = max(policies.values()) - min(policies.values())
    if spread == 0 and widest == 0:
        raise AuditError(
            "neither the threshold nor the policy axis moves, so R4 is not a "
            "counterfactual")
    # The baseline cutoff appears in the threshold variants, which makes it the
    # identity control of that axis: it must reproduce the headline absence count
    # or the two are not measuring the same thing. The kinetics ladder carries the
    # same control for the same reason.
    baseline = gold["threshold_absent"].get(gold["threshold"])
    if baseline != absent:
        raise AuditError(
            f"the baseline cutoff {gold['threshold']} re-runs to {baseline} and "
            f"the headline absence count is {absent}; R4's identity control is "
            "not controlling anything")

    if strict is not None and strict != absent:
        raise AuditError(
            f"the strict policy counts {strict} and the canonical-absent set is "
            f"{absent}; the policy that accepts only the canonical model must "
            "reproduce the absence count exactly or the two are measuring "
            "different things")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", PRESENCE_KEYS, "presence.json"),
                             ("r4", ARCHITECTURE_KEYS, "architecture.json")):
        (warm / rung).mkdir(parents=True)
        (warm / rung / name).write_text(
            json.dumps({k: gold[k] for k in keys}, indent=2, sort_keys=True) + "\n")


def build_grading(campaign: pathlib.Path, gold: dict) -> dict:
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    grading = task.get("grading") or {}
    chance = {**DEFAULT_CHANCE_LEVELS, **(grading.get("chance_levels") or {})}
    notes = grading.get("chance_level_notes") or {}

    def ex(name, path=None):
        return {"name": name, "scorer": "exact",
                "args": {"pred": f"pred:{name}", "gold": f"gold:{path or name}"}}

    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("release"), ex("ec_number"), ex("canonical_model"), ex("threshold"),
          ex("sequences_total"), ex("sequences_by_section"),
          ex("model_panel_size"), ex("environment")]

    r2 = [ex("canonical_present_by_section"), ex("canonical_absent_by_section"),
          ex("canonical_absent_total"), ex("length_profile")]

    r3 = [ex("family_census"), ex("architecture_census"),
          ex("absent_family_census"), ex("absent_architecture_census"),
          ex("dominant_absent_architecture")]

    r4 = [ex("misannotated_by_policy"), ex("threshold_absent")]

    rungs = []
    for rid, ordinal, label, components in (("r1", 1, "execute", r1),
                                            ("r2", 2, "correct", r2),
                                            ("r3", 3, "commit", r3),
                                            ("r4", 4, "counterfactual", r4)):
        rung = {"rung_id": rid, "ordinal": ordinal, "name": label, "gold_tier": "G1",
                "composition": "product", "pass_threshold": 1.0,
                "chance_level": chance[rid], "components": components}
        if rid in notes:
            rung["chance_level_note"] = notes[rid]
        rungs.append(rung)

    spec = {"campaign_id": task["campaign_id"], "template_id": task["template_id"],
            "grading_spec_version": "1.0.0", "rungs": rungs}
    (pathlib.Path(campaign) / "grading.json").write_text(
        json.dumps(spec, indent=2) + "\n")
    return spec


def _flatten_for_grading(gold: dict) -> dict:
    """Group the tool and library identity into one graded composite.

    The three fields are one claim -- which tool and which library produced this
    verdict -- and grading them separately would let a system that pinned nothing
    lose one component instead of the claim.
    """
    gold["environment"] = {"hmmer_version": gold["hmmer_version"],
                           "pfam_version": gold["pfam_version"],
                           "pfam_sha256": gold["pfam_sha256"]}
    return gold


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)

    gold = _flatten_for_grading(analyse(campaign))
    _assert_has_content(gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    print(f"{campaign.name}: UniProt {gold['release']}, EC {gold['ec_number']}, "
          f"{gold['sequences_total']} sequences {gold['sequences_by_section']}")
    print(f"  {gold['hmmer_version']} against Pfam {gold['pfam_version']}, "
          f"panel of {gold['model_panel_size']} models at {gold['threshold']}")
    for section in sorted(gold["canonical_present_by_section"]):
        present = gold["canonical_present_by_section"][section]
        absent = gold["canonical_absent_by_section"][section]
        total = present + absent
        print(f"    {section:<11} {gold['canonical_model']} present {present:>5}"
              f"/{total:<5} absent {absent:>5} ({absent / total:.4f})")
    profile = gold["length_profile"]
    print(f"  length median: present {profile['canonical_present']['median']:.0f}, "
          f"absent {profile['canonical_absent']['median']:.0f} -- the absent side "
          "is not a set of fragments")
    print(f"  dominant absent architecture: "
          f"{gold['dominant_absent_architecture']} "
          f"({next(iter(gold['absent_architecture_census'].values()))})")
    print(f"  misannotated by policy: {gold['misannotated_by_policy']}")
    print(f"  absent by threshold: {gold['threshold_absent']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
