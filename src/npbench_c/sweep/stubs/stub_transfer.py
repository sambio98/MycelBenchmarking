"""Deterministic stub systems for the annotation-transfer ladder.

Harness fixtures. Sandbox-only plus the pinned DIAMOND, which every level above
`noncompute` runs.

  reader      clears R1: the set inventory, including how many EC mentions stop
              early and how many entries carry several codes
  search      clears R2: the search census, the tied-top count and the identity
              distribution of the closest hits
  verdicts    clears R3: agreement by depth and by identity band, the verdict
              census and the committed error set
  complete    clears R4
  dashvalue   COMPETENCE probe, tooled: runs the search correctly and compares EC
              numbers as STRINGS, so `1.14.-.-` is treated as a value rather than
              as a missing claim. Clears R2 and fails R3 -- every search number
              right, every agreement number wrong, and the undecidable cell empty.
  noncompute  LEAKAGE ablation: never runs the tool and never opens the input.

`dashvalue` is the probe worth having here. The search is the visible work and it
is correct; what is wrong is one rule in the comparison, and the symptom is a
report that looks more complete than gold rather than less -- 902 queries move out
of `undecidable` and into verdicts nobody can support.
"""

from __future__ import annotations

import argparse
import collections
import json
import pathlib

from npbench_c.templates.annotation_transfer import analyse as atr
from npbench_c.templates.annotation_transfer.core import Transfer

LEVELS = ("reader", "search", "verdicts", "complete", "dashvalue", "noncompute")
REPORT = "transfer_report.json"


def _dash_agreement(query, subject, depth):
    """The mistake: compare the first `depth` fields as strings, dash and all."""
    def level(code):
        parts = code.split(".")
        return ".".join(parts[:depth]) if len(parts) == 4 else None
    left = {level(c) for c in query} - {None}
    right = {level(c) for c in subject} - {None}
    if not left or not right:
        return None
    return bool(left & right)


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
            "release": "2026_03", "ec_prefix": "1.14",
            "tie_break": "bitscore", "diamond_version": "diamond version 2.2.8",
            "entries_total": 4000, "distinct_ec_codes": 500,
            "ec_codes_per_entry": {}, "partial_ec_mentions": 1000,
            "entries_with_multiple_ec": 400,
            "queries_with_a_hit": 4000, "queries_without_a_hit": 0,
            "tied_top_queries": 0, "top_hit_identity_census": {},
            "agreement_by_depth": {}, "agreement_by_identity": {},
            "verdict_census": {}, "transfer_errors_at_subsubclass": [],
            "removal_census": {}, "removal_changed": {},
        })

    params = atr.load_params(sandbox)
    transfer = Transfer.load(sandbox, params["sequence_input"],
                             params["annotation_input"])
    rules = transfer.rules
    hits = transfer.hits()

    report = {
        "release": transfer.release,
        "ec_prefix": rules.ec_prefix,
        "tie_break": rules.tie_break,
        "diamond_version": transfer.diamond_version(),
        "entries_total": len(transfer.sequences),
        **transfer.annotation_census(),
    }
    if args.level == "reader":
        return write(report)

    report.update({
        "queries_with_a_hit": len(hits),
        "queries_without_a_hit": len(transfer.sequences) - len(hits),
        "tied_top_queries": transfer.tied_top_queries(hits),
        "top_hit_identity_census": transfer.identity_census(hits),
    })
    if args.level == "search":
        return write(report)

    deepest = max(rules.graded_depths)
    if args.level == "dashvalue":
        # The search is right; the comparison treats a dash as a value.
        depth_block = {}
        for depth in rules.graded_depths:
            decidable = agree = 0
            for query in sorted(hits):
                top = transfer.ranked(hits[query])[0]
                verdict = _dash_agreement(transfer.annotations[query],
                                          transfer.annotations[top.subject], depth)
                if verdict is None:
                    continue
                decidable += 1
                agree += verdict
            depth_block[str(depth)] = {
                "decidable": decidable, "agree": agree,
                "rate": round(agree / decidable, 6) if decidable else None}
        band_block = {}
        for lo, hi in rules.identity_bands:
            decidable = agree = 0
            for query in sorted(hits):
                top = transfer.ranked(hits[query])[0]
                if not (lo <= top.identity < hi):
                    continue
                verdict = _dash_agreement(transfer.annotations[query],
                                          transfer.annotations[top.subject], deepest)
                if verdict is None:
                    continue
                decidable += 1
                agree += verdict
            band_block[rules.band_label(lo, hi)] = {
                "decidable": decidable, "agree": agree,
                "rate": round(agree / decidable, 6) if decidable else None}
        census = collections.Counter()
        errors = []
        for query in sorted(hits):
            top = transfer.ranked(hits[query])[0]
            mine = transfer.annotations[query]
            theirs = transfer.annotations[top.subject]
            at_three = _dash_agreement(mine, theirs, 3)
            at_four = _dash_agreement(mine, theirs, 4)
            if at_three is False:
                census["wrong_at_subsubclass"] += 1
                errors.append(query)
            elif at_four is False:
                census["wrong_at_serial"] += 1
            else:
                census["supported"] += 1
        report["agreement_by_depth"] = depth_block
        report["agreement_by_identity"] = band_block
        report["verdict_census"] = {v: census.get(v, 0) for v in rules.verdict_order}
        report["transfer_errors_at_subsubclass"] = sorted(errors)
        return write(report)

    verdicts = transfer.verdicts(hits)
    report["agreement_by_depth"] = transfer.agreement_by_depth(
        hits, rules.graded_depths)
    report["agreement_by_identity"] = transfer.agreement_by_identity(hits, deepest)
    report["verdict_census"] = transfer.verdict_census(verdicts)
    report["transfer_errors_at_subsubclass"] = sorted(
        q for q, v in verdicts.items() if v == "wrong_at_subsubclass")
    if args.level == "verdicts":
        return write(report)

    report["removal_census"] = {
        name: transfer.census_after_removal(hits, remove)
        for name, remove in sorted(rules.removal_variants.items())}
    report["removal_changed"] = {
        name: sum(transfer.verdict_transitions(hits, remove).values()) if remove else 0
        for name, remove in sorted(rules.removal_variants.items())}
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
