"""Build gold, warm bundles and the grading spec for the dereplication ladder.

Usage: python -m npbench_c.templates.ms2_dereplication.build <campaign>
"""

from __future__ import annotations

import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.templates.ms2_dereplication.analyse import analyse
from npbench_c.templates.ms2_dereplication.core import (
    CONSTANTS,
    RULES,
    DereplicationError,
    Rules,
)

DEFAULT_CHANCE_LEVELS = {"r1": 1.0, "r2": 0.001, "r3": 0.001, "r4": 0.001}

#: R1's report, and therefore exactly what R2's warm bundle may hand over.
EXECUTE_KEYS = ("environment", "corpus_census", "candidate_census")
#: R2's answers: what the search retrieved for each query. NOT which of them are
#: right: that needs the compound each query is, which is the answer and is not
#: in the sandbox. Gold carries `correct_queries` and `misannotated_queries` for
#: the audit and no component reads them.
RETRIEVAL_KEYS = ("retrieval_census", "top_matches")
#: R3's answers: the declared confidence class, and the ambiguity admitted to.
COMMITMENT_KEYS = ("confidence_by_query", "confidence_census", "margin_summary",
                   "ambiguous_queries")
#: R4's answers.
COUNTERFACTUAL_KEYS = ("tolerance_movement", "adduct_misassignment_census")

#: Below this the retrieval is broken rather than imperfect, and R2's census
#: would be measuring a bad split instead of a hard search.
MIN_CORRECT_FRACTION = 0.6
#: A competing compound has to exist often enough for the margin to be a
#: distribution rather than an anecdote.
MIN_COMPETING = 10


def copy_constants(campaign: pathlib.Path) -> None:
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSTANTS / RULES, reference / RULES)


def _assert_has_content(campaign: pathlib.Path, gold: dict) -> None:
    """Refuse to emit gold for a campaign whose rungs have nothing to say."""
    rules = Rules.load(campaign)

    if gold["environment"]["matchms_version"] != rules.tool_version:
        raise DereplicationError(
            f"the pin names matchms {rules.tool_version} and the image answered "
            f"{gold['environment']['matchms_version']}; every score in gold was "
            "produced by a library the campaign does not declare")

    census = gold["corpus_census"]
    if census["library_compounds"] < 100:
        raise DereplicationError(
            f"a library of {census['library_compounds']} compounds cannot carry "
            "a dereplication question: the right answer is nearly the only answer")
    retrieval = gold["retrieval_census"]
    if sum(retrieval.values()) != census["graded_queries"]:
        raise DereplicationError(
            f"the retrieval census covers {sum(retrieval.values())} queries and "
            f"{census['graded_queries']} were graded")

    candidates = gold["candidate_census"]
    if candidates["max"] < 2:
        raise DereplicationError(
            "no query's precursor window holds two candidates, so nothing is "
            "ever chosen between and the search is a lookup")

    correct = gold["correct_queries"]
    floor = MIN_CORRECT_FRACTION * census["graded_queries"]
    if len(correct) < floor:
        raise DereplicationError(
            f"only {len(correct)} of {census['graded_queries']} queries retrieve "
            f"their own compound, below the {floor:g} this template requires; "
            "that is a broken split rather than a hard search")
    if not gold["misannotated_queries"]:
        raise DereplicationError(
            "every annotated query retrieves its own compound, so the "
            "misannotation R2's census exists to count does not occur here")
    if set(correct) & set(gold["misannotated_queries"]):
        raise DereplicationError(
            "a query is counted both correct and misannotated")
    if len(gold["top_matches"]) != census["graded_queries"]:
        raise DereplicationError(
            f"{len(gold['top_matches'])} queries carry a top match and "
            f"{census['graded_queries']} were graded")

    confidence = gold["confidence_census"]
    substantive = [label for label in confidence
                   if label not in ("no_candidate", "no_match")]
    empty = sorted(label for label in substantive if not confidence[label])
    if empty:
        raise DereplicationError(
            f"the confidence classes {empty} are empty on this corpus, so the "
            "declared rule is not separating the annotations it classifies")
    if sum(confidence.values()) != census["graded_queries"]:
        raise DereplicationError(
            f"the confidence census counts {sum(confidence.values())} queries "
            f"and {census['graded_queries']} were graded")
    if len(gold["confidence_by_query"]) != census["graded_queries"]:
        raise DereplicationError(
            f"{len(gold['confidence_by_query'])} queries carry a confidence "
            f"class and {census['graded_queries']} were graded")

    margin = gold["margin_summary"]
    if margin["queries_with_a_competing_compound"] < MIN_COMPETING:
        raise DereplicationError(
            f"a competing compound exists for only "
            f"{margin['queries_with_a_competing_compound']} queries, fewer than "
            f"the {MIN_COMPETING} this template requires for a margin summary")
    if not gold["ambiguous_queries"]:
        raise DereplicationError(
            "no annotation sits inside the declared ambiguity margin, so R3's "
            "ambiguous set is empty and the rung grades a constant")
    if len(gold["ambiguous_queries"]) == margin["queries_with_a_competing_compound"]:
        raise DereplicationError(
            "every query with a competitor is inside the margin, so the margin "
            "separates nothing and the declared threshold does no work")

    movement = gold["tolerance_movement"]
    if len(movement) < 3:
        raise DereplicationError(
            f"a grid of {len(movement)} tolerances cannot show a response to the "
            "tolerance")
    pinned = f"{rules.peak_tolerance}"
    if movement.get(pinned):
        raise DereplicationError(
            f"the pinned tolerance {pinned} moves "
            f"{len(movement[pinned])} answers against itself, which means the "
            "search is not deterministic under this invocation and every "
            "movement above it is uncalibrated")
    if not any(queries for key, queries in movement.items() if key != pinned):
        raise DereplicationError(
            "no declared tolerance moves a single answer, so the grid is a "
            "parameter the campaign could have left undeclared")
    adduct = gold["adduct_misassignment_census"]
    if sum(adduct.values()) != census["graded_queries"]:
        raise DereplicationError(
            f"the adduct census counts {sum(adduct.values())} queries and "
            f"{census['graded_queries']} were graded")
    if not adduct["no_candidate"]:
        raise DereplicationError(
            "the misassigned adduct leaves every query a candidate window, so "
            "the counterfactual this campaign measures does not occur")
    # The reachable degeneracy: a perturbation that leaves every annotation
    # naming the same compound it named before, which would make R4's census a
    # restatement of R2's answer.
    if adduct["unchanged"] >= retrieval["annotated"]:
        raise DereplicationError(
            f"{adduct['unchanged']} annotations name the same compound after the "
            f"misassigned adduct as before it and {retrieval['annotated']} "
            "queries were annotated, so the perturbation costs the search "
            "nothing and R4's census restates R2's")


def write_warm_bundles(campaign: pathlib.Path, gold: dict) -> None:
    warm = pathlib.Path(campaign) / "gold" / "warm"
    if warm.exists():
        shutil.rmtree(warm)
    for rung, keys, name in (("r2", EXECUTE_KEYS, "execution.json"),
                             ("r3", RETRIEVAL_KEYS, "retrieval.json"),
                             ("r4", COMMITMENT_KEYS, "commitments.json")):
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

    # `corpus_census` is reported and NOT graded: all three of its numbers are in
    # inputs/provenance.json, which the agent is given, so a component for it
    # would measure a copy out of the sandbox. The tool identity is graded as ONE
    # claim that also carries the library's version, which the sandbox does not
    # state, and `candidate_census` is the filter actually applied.
    r1 = [{"name": "report_present", "scorer": "exact",
           "args": {"pred": "pred:present.report", "gold": 1}},
          ex("environment"), ex("candidate_census")]

    # Which retrievals are RIGHT is not graded and not reported: it needs the
    # compound each query is, the queries ship without it, and a component no
    # system can compute is not a component. Gold keeps both sets for the audit.
    r2 = [ex("retrieval_census"), ex("top_matches")]
    r3 = [ex("confidence_by_query"), ex("confidence_census"), ex("margin_summary"),
          ex("ambiguous_queries")]
    r4 = [ex("tolerance_movement"), ex("adduct_misassignment_census")]

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


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: build <campaign_dir>", file=sys.stderr)
        return 2
    campaign = pathlib.Path(argv[0]).resolve()
    copy_constants(campaign)

    gold = analyse(campaign)
    _assert_has_content(campaign, gold)
    (campaign / "gold").mkdir(exist_ok=True)
    (campaign / "gold" / "gold.json").write_text(
        json.dumps(gold, indent=2, sort_keys=True) + "\n")
    write_warm_bundles(campaign, gold)
    spec = build_grading(campaign, gold)

    env, census = gold["environment"], gold["corpus_census"]
    print(f"{campaign.name}: {census['graded_queries']} queries against "
          f"{census['library_spectra']} library spectra of "
          f"{census['library_compounds']} compounds, matchms "
          f"{env['matchms_version']} {env['similarity']} at "
          f"{env['peak_tolerance']} Da")
    print(f"  candidate sets: {gold['candidate_census']}")
    print(f"  retrieval: {gold['retrieval_census']}; of those, not graded: "
          f"{len(gold['correct_queries'])} correct, "
          f"{len(gold['misannotated_queries'])} misannotated")
    print(f"  confidence: {gold['confidence_census']}")
    print(f"  margins: {gold['margin_summary']}, ambiguous "
          f"{len(gold['ambiguous_queries'])}")
    print("  tolerance movement: "
          f"{ {k: len(v) for k, v in gold['tolerance_movement'].items()} }")
    print(f"  adduct misassignment: {gold['adduct_misassignment_census']}")
    print(f"  adduct against the truth (not graded): "
          f"{gold['adduct_against_the_truth']}")
    print(f"  grading.json: {sum(len(r['components']) for r in spec['rungs'])} "
          "components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
