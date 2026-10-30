"""Grader purity (NPBench-C spec 3.3, 9.1).

The grader must be a pure function of (submission, gold, grader_version): same
inputs, byte-identical output, on any machine at any time. Purity does not make
a score correct -- that is the separate construct-validity study -- but it is
the part that is fully under our control.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.grading.grade import grade, grade_json
from npbench_c.grading.version import content_hash, grader_version

REPEATS = 100
CAMPAIGN = pathlib.Path(__file__).resolve().parents[2] / "campaigns" / "construct-ecoli-rebh-01"


def _fixtures():
    return (
        json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text()),
        json.loads((CAMPAIGN / "gold" / "gold.json").read_text()),
        json.loads((CAMPAIGN / "grading.json").read_text()),
    )


@pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                    reason="campaign gold not generated")
def test_repeat_grading_is_byte_identical():
    m, g, s = _fixtures()
    reference = grade_json(m, g, s)
    assert all(grade_json(m, g, s) == reference for _ in range(REPEATS))


@pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                    reason="campaign gold not generated")
def test_grading_does_not_mutate_its_inputs():
    m, g, s = _fixtures()
    before = (json.dumps(m, sort_keys=True), json.dumps(g, sort_keys=True),
              json.dumps(s, sort_keys=True))
    grade(m, g, s)
    after = (json.dumps(m, sort_keys=True), json.dumps(g, sort_keys=True),
             json.dumps(s, sort_keys=True))
    assert before == after


@pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                    reason="campaign gold not generated")
def test_key_insertion_order_does_not_change_the_score():
    """A submission serialised with keys in a different order is the same
    submission. Hash-ordering dependence would show up here."""
    m, g, s = _fixtures()
    reference = grade_json(m, g, s)
    reordered = json.loads(json.dumps(m, sort_keys=True))
    shuffled = {k: reordered[k] for k in reversed(list(reordered))}
    assert grade_json(shuffled, g, s) == reference


def test_grader_version_is_stable_within_a_process():
    assert grader_version() == grader_version()
    assert content_hash() == content_hash()
    assert "+" in grader_version()


def test_scoring_path_imports_no_third_party_modules():
    """Third-party numeric libraries are excluded from the scoring path so that
    a numpy or pandas version bump cannot move a published score."""
    import ast

    allowed = {
        "__future__", "math", "re", "json", "hashlib", "pathlib", "typing",
        "dataclasses", "npbench_c",
    }
    grading_dir = pathlib.Path(__file__).resolve().parents[2] / "src" / "npbench_c" / "grading"
    for path in sorted(grading_dir.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert alias.name.split(".")[0] in allowed, f"{path.name}: {alias.name}"
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                assert node.module.split(".")[0] in allowed, f"{path.name}: {node.module}"


def test_null_prediction_scores_zero_instead_of_crashing():
    """An omitted field is a wrong answer, not a grader error.

    Regression: a measurer records an omitted numeric field as JSON null, which
    resolves successfully to None rather than missing. Passing that to a numeric
    scorer raised TypeError and aborted the entire grade -- one omitted field
    would have taken down a whole sweep.
    """
    spec = {
        "campaign_id": "t",
        "rungs": [{
            "rung_id": "r1", "ordinal": 1, "chance_level": 0.1,
            "components": [
                {"name": "mass", "scorer": "abs_tol_match",
                 "args": {"pred": "pred:mass", "gold": "gold:mass", "abs_tol": 0.001}},
                {"name": "label", "scorer": "exact",
                 "args": {"pred": "pred:label", "gold": "gold:label"}},
            ],
        }],
    }
    gold = {"mass": 508.2129, "label": "x"}
    result = grade({"mass": None, "label": None}, gold, spec)
    assert result["rungs"][0]["components"] == {"label": 0.0, "mass": 0.0}
    assert result["score"] == 0.0


def test_null_gold_raises_as_a_campaign_bug():
    spec = {
        "campaign_id": "t",
        "rungs": [{
            "rung_id": "r1", "ordinal": 1,
            "components": [{"name": "m", "scorer": "exact",
                            "args": {"pred": "pred:m", "gold": "gold:m"}}],
        }],
    }
    with pytest.raises(Exception):
        grade({"m": 1}, {"m": None}, spec)
