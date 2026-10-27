"""Task-level aggregation: pass fraction, soft artifacts, format compliance."""

import pytest
import yaml

from npbench.grading.config import load_grading_config
from npbench.grading.score import score_task

ROWS = "bgc_id,product_class\nBGC0001,nrps\nBGC0002,t1pks\nBGC0003,terpene\nBGC0004,nrps\n"


def _artifact(name, soft=False, **overrides):
    return {
        "path": f"{name}.csv",
        "truth": f"{name}.csv",
        "format": "csv",
        "soft": soft,
        "row_key": ["bgc_id"],
        "normalisers": ["strip", "casefold"],
        "column_aliases": {},
        "value_aliases": {},
        "metrics": [{"type": "identifier_f1"}],
        **overrides,
    }


@pytest.fixture
def task(tmp_path, write):
    """Build a task with three artifacts and return a runner over submissions."""

    def _task(artifacts, submissions):
        for spec in artifacts:
            write(f"truth/{spec['truth']}", ROWS)
        for name, content in submissions.items():
            write(f"outputs/{name}", content)
        path = tmp_path / "grading.yaml"
        path.write_text(yaml.safe_dump({"version": 1, "artifacts": artifacts}))
        config = load_grading_config(path)
        return score_task(config, tmp_path / "truth", tmp_path / "outputs")

    return _task


def test_task_score_is_the_fraction_of_hard_artifacts_passing(task):
    # a: perfect. b: one of four rows missing -> F1 = 6/7 = 0.857, passes.
    # c: not written at all -> 0.00, fails.
    result = task(
        [_artifact("a"), _artifact("b"), _artifact("c")],
        {
            "a.csv": ROWS,
            "b.csv": "bgc_id,product_class\nBGC0001,nrps\nBGC0002,t1pks\nBGC0003,terpene\n",
        },
    )
    assert result.task_score == pytest.approx(2 / 3)
    assert (result.n_hard, result.n_hard_passed) == (3, 2)
    assert [a.score for a in result.artifacts] == [
        pytest.approx(1.0), pytest.approx(6 / 7), 0.0
    ]


def test_soft_artifacts_are_scored_but_kept_out_of_the_headline(task):
    # PROJECT.md section 9: a library match table depends on a live service, so
    # it is reported and excluded rather than dropped.
    result = task(
        [_artifact("a"), _artifact("library_hits", soft=True)],
        {"a.csv": ROWS},
    )
    assert result.task_score == 1.0
    assert (result.n_hard, result.n_hard_passed) == (1, 1)
    assert len(result.soft_artifacts) == 1
    assert result.soft_artifacts[0].score == 0.0


def test_format_compliance_is_a_separate_axis_from_score(task):
    # Both artifacts are scientifically perfect. One uses a renamed column.
    renamed = ROWS.replace("bgc_id", "BGC_ID", 1)
    result = task(
        [_artifact("a"), _artifact("b", column_aliases={"bgc_id": ["BGC_ID"]})],
        {"a.csv": ROWS, "b.csv": renamed},
    )
    assert result.task_score == 1.0
    assert result.format_compliance == pytest.approx(0.5)


def test_format_compliance_ignores_artifacts_that_were_never_written(task):
    # A missing file is neither compliant nor non-compliant; counting it as
    # non-compliant would conflate two different failures.
    result = task([_artifact("a"), _artifact("c")], {"a.csv": ROWS})
    assert result.format_compliance == 1.0
    assert result.task_score == 0.5


def test_a_failing_run_scores_zero_without_crashing(task):
    result = task([_artifact("a"), _artifact("b")], {})
    assert result.task_score == 0.0
    assert all(a.reason == "file_missing" for a in result.artifacts)
    assert "file_missing" in result.summary()
