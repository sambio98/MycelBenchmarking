"""Build gold, warm bundles and the grading spec for the remote-homology ladder.

Usage: python -m npbench_c.templates.remote_homology.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.remote_homology.analyse import analyse
from npbench_c.templates.remote_homology.core import CONSTANTS, RULES, HomologyError

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's answers, and therefore exactly what R2's warm bundle may hand over.
INVENTORY_KEYS = ("label_resource", "corpus_superfamily", "sensitivity", "evalue",
                  "max_seqs", "threads", "format_fields", "mmseqs_version",
                  "corpus_sequences", "corpus_sha256", "query_table",
                  "baseline_hits", "baseline_band_sizes")
#: R2's answers: what the search found of what the resource asserts.
RECOVERY_KEYS = ("recovery", "band_census")
#: R3's answers: the declared verdicts.
VERDICT_KEYS = ("band_verdicts", "verdict_census")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("shuffle_ladder",)


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    if gold["corpus_sequences"] < 500:
        raise HomologyError(
            f"a corpus of {gold['corpus_sequences']} sequences cannot carry a "
            "per-band census of a twilight zone")

    queries = gold["query_table"]
    if len(queries) < 3:
        raise HomologyError(
            f"{len(queries)} declared queries cannot carry a per-query census")
    accessions = [row["accession"] for row in queries.values()]
    if len(set(accessions)) != len(accessions):
        raise HomologyError(
            "two declared families selected the same query, so one family is "
            "audited twice and a family stands for nothing")
    thin = sorted(entry for entry, row in queries.items()
                  if row["family_members"] < 10)
    if thin:
        raise HomologyError(
            f"the declared families {thin} cover fewer than ten corpus members "
            "each, so their recovery rate is an anecdote rather than a census")

    hits = gold["baseline_hits"]
    if set(hits) != set(accessions):
        raise HomologyError(
            "the baseline hit counts and the query table describe different "
            "query sets, so a count belongs to a query nobody declared")
    silent = sorted(accession for accession, count in hits.items() if not count)
    if silent:
        raise HomologyError(
            f"the search reports no hit for {silent}, which is a broken run "
            "rather than a finding about remote homology")

    recovery = gold["recovery"]
    missed = {accession: row["same_family_missed"] for accession, row in recovery.items()}
    if not any(missed.values()):
        raise HomologyError(
            "the search recovers every relationship the resource asserts, so the "
            "gap this campaign exists to measure does not occur on this corpus "
            "at this sensitivity")
    if not all(row["same_family_asserted"] for row in recovery.values()):
        raise HomologyError(
            "some declared query has no asserted family relative at all, so its "
            "recovery rate is undefined and its row measures nothing")

    census = gold["band_census"]
    overlap = sum(1 for cells in census.values() for cell in cells.values()
                  if cell["same_family"] and cell["superfamily_only"])
    if not overlap:
        raise HomologyError(
            "no band holds both same-family and superfamily-only hits, so the "
            "overlap that makes identity alone insufficient does not occur and "
            "the campaign is measuring a clean separation")

    verdicts = gold["verdict_census"]
    populated = sorted(v for v, n in verdicts.items() if n)
    if len(populated) < 2:
        raise HomologyError(
            f"only {populated} verdicts occur across every query and band, so the "
            "vocabulary is not separating anything")
    if not verdicts["overlapping"]:
        raise HomologyError(
            "no band is classified as overlapping, which is the twilight-zone "
            "finding the campaign is named for")
    if sum(verdicts.values()) != sum(len(cells) for cells in census.values()):
        raise HomologyError(
            f"the verdict census counts {sum(verdicts.values())} cells and the "
            f"band census has {sum(len(cells) for cells in census.values())}")

    ladder = gold["shuffle_ladder"]
    if "0" not in ladder:
        raise HomologyError(
            "the shuffle ladder has no zero-window rung, so there is no identity "
            "control and every loss is uncalibrated")
    control = ladder["0"]
    noisy = sorted(accession for accession, row in control.items() if row["lost"])
    if noisy:
        raise HomologyError(
            f"the zero-window control loses hits for {noisy}, which means the "
            "unperturbed query does not reproduce its own baseline and the "
            "search is not deterministic under this invocation")
    top = ladder[str(max(int(w) for w in ladder))]
    if not any(row["lost"] for row in top.values()):
        raise HomologyError(
            "the widest declared window costs no query a single hit, so the "
            "ladder never reaches a perturbation the search notices")
    if len(ladder) < 3:
        raise HomologyError(
            f"a ladder of {len(ladder)} windows cannot show a response to window "
            "length")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", INVENTORY_KEYS, "inventory.json"),
                             ("r3", RECOVERY_KEYS, "recovery.json"),
                             ("r4", VERDICT_KEYS, "verdicts.json")):
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

    # `label_resource`, `corpus_superfamily` and every search flag are declared in
    # reference/homology_rules.json, so a component for one measures a copy out of
    # the sandbox. The tool identity is graded as ONE claim that also carries the
    # search binary's version, which the sandbox does not state.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("corpus_sequences"), ex("corpus_sha256"), ex("query_table"),
          ex("baseline_hits"), ex("baseline_band_sizes"), ex("environment")]

    r2 = [ex("recovery"), ex("band_census")]
    r3 = [ex("band_verdicts"), ex("verdict_census")]
    r4 = [ex("shuffle_ladder")]

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
    """The search identity is one claim: which binary, at which sensitivity,
    under which reporting threshold, cap, thread count and output projection."""
    gold["environment"] = {
        "mmseqs_version": gold["mmseqs_version"],
        "sensitivity": gold["sensitivity"],
        "evalue": gold["evalue"],
        "max_seqs": gold["max_seqs"],
        "threads": gold["threads"],
        "format_fields": gold["format_fields"],
    }
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

    print(f"{campaign.name}: {gold['corpus_sequences']} corpus sequences in "
          f"{gold['corpus_superfamily']}, {gold['mmseqs_version']} at "
          f"s={gold['sensitivity']:g}")
    print("  declared queries and what the search found:")
    for entry, row in sorted(gold["query_table"].items()):
        acc = row["accession"]
        rec = gold["recovery"][acc]
        print(f"    {entry} {acc:<11} fam={row['family_members']:<4} "
              f"hits={rec['hits']:<5} asserted={rec['same_family_asserted']:<4} "
              f"found={rec['same_family_detected']:<4} "
              f"missed={rec['same_family_missed']:<4} "
              f"rate={rec['recovery_rate']}")
    print("  band census (same-family / superfamily-only / unclassified):")
    for acc, cells in sorted(gold["band_census"].items()):
        parts = " ".join(
            f"{band}:{c['same_family']}/{c['superfamily_only']}/{c['unclassified_pair']}"
            for band, c in cells.items())
        print(f"    {acc:<11} {parts}")
    print(f"  verdicts: {gold['verdict_census']}")
    print("  shuffle ladder (lost, and loss rate by band):")
    for window, rows in sorted(gold["shuffle_ladder"].items(), key=lambda kv: int(kv[0])):
        for acc, row in sorted(rows.items()):
            rates = " ".join(f"{b}:{row['loss_rate_by_band'][b]}"
                             for b in row["loss_rate_by_band"])
            print(f"    w={window:<4} {acc:<11} hits={row['hits']:<5} "
                  f"lost={row['lost']:<4} order={row['band_order_held']} {rates}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
