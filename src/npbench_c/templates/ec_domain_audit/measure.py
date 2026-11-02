"""Measure an EC domain-audit submission. Deterministic, agent-independent.

Nothing here reads gold: the sections, the model panel, the policy names, the
threshold variants and the length bands all come from the campaign's own declared
vocabularies.

Usage:
  python -m npbench_c.templates.ec_domain_audit.measure \
      <campaign> <submission> <out>
"""

from __future__ import annotations

import json
import pathlib
import sys

from npbench_c.templates.ec_domain_audit.core import SECTIONS, Rules

REPORT = "domain_audit.json"


def _int(value: object) -> int | None:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def _float(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _int_map(value: object, keys) -> dict:
    block = value if isinstance(value, dict) else {}
    return {k: _int(block.get(k)) for k in keys}


def _free_int_map(value: object) -> dict:
    """Int map over the submission's own keys, zeros dropped: gold omits them.

    The architecture census is open-ended -- its keys are whichever domain
    contents the set turns out to have -- so it cannot be read against a declared
    key list. Dropping zeros keeps a submission that spells out absent
    architectures equal to gold, which omits them.
    """
    block = value if isinstance(value, dict) else {}
    return {str(k): n for k, v in sorted(block.items())
            if (n := _int(v)) is not None and n != 0}


def measure(campaign: pathlib.Path, submission_dir: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    rules = Rules.load(campaign)
    sections = tuple(sorted(set(SECTIONS.values())))
    bands = tuple(str(b) for b in rules.length_bands)
    policies = tuple(sorted(rules.misannotation_policies))
    thresholds = tuple(sorted(rules.threshold_variants))
    sides = ("canonical_absent", "canonical_present")

    path = pathlib.Path(submission_dir) / REPORT
    out: dict = {"present": {"report": int(path.is_file())}}
    if not path.is_file():
        return out
    try:
        sub = json.loads(path.read_text())
    except json.JSONDecodeError:
        out["present"]["report"] = 0
        return out
    if not isinstance(sub, dict):
        out["present"]["report"] = 0
        return out

    # ---- R1: the inventory, and what produced it
    for key in ("release", "ec_number", "canonical_model", "threshold"):
        out[key] = _str(sub.get(key))
    out["sequences_total"] = _int(sub.get("sequences_total"))
    out["sequences_by_section"] = _int_map(sub.get("sequences_by_section"), sections)
    out["model_panel_size"] = _int(sub.get("model_panel_size"))
    out["environment"] = {
        "hmmer_version": _str(sub.get("hmmer_version")),
        "pfam_version": _str(sub.get("pfam_version")),
        "pfam_sha256": _str(sub.get("pfam_sha256")),
    }

    # ---- R2: presence, and the fragment explanation
    out["canonical_present_by_section"] = _int_map(
        sub.get("canonical_present_by_section"), sections)
    out["canonical_absent_by_section"] = _int_map(
        sub.get("canonical_absent_by_section"), sections)
    out["canonical_absent_total"] = _int(sub.get("canonical_absent_total"))
    profile = sub.get("length_profile")
    profile = profile if isinstance(profile, dict) else {}
    measured_profile = {}
    for side in sides:
        block = profile.get(side)
        block = block if isinstance(block, dict) else {}
        measured_profile[side] = {
            "count": _int(block.get("count")),
            "min": _int(block.get("min")),
            "median": _float(block.get("median")),
            "max": _int(block.get("max")),
            "shorter_than": _int_map(block.get("shorter_than"), bands),
        }
    out["length_profile"] = measured_profile

    # ---- R3: what the canonical-absent sequences actually are
    out["family_census"] = _int_map(sub.get("family_census"), rules.panel)
    out["absent_family_census"] = _int_map(sub.get("absent_family_census"),
                                           rules.panel)
    out["architecture_census"] = _free_int_map(sub.get("architecture_census"))
    out["absent_architecture_census"] = _free_int_map(
        sub.get("absent_architecture_census"))
    out["dominant_absent_architecture"] = _str(
        sub.get("dominant_absent_architecture"))

    # ---- R4: the two counterfactual axes
    out["misannotated_by_policy"] = _int_map(sub.get("misannotated_by_policy"),
                                             policies)
    out["threshold_absent"] = _int_map(sub.get("threshold_absent"), thresholds)
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print("usage: measure <campaign_dir> <submission_dir> <out_json>",
              file=sys.stderr)
        return 2
    result = measure(pathlib.Path(argv[0]), pathlib.Path(argv[1]))
    pathlib.Path(argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
