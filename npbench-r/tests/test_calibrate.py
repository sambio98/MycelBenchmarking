"""Difficulty calibration arithmetic, with hand-computed expectations."""

import json

import pytest

from npbench.calibrate import (
    BAND_HIGH,
    BAND_LOW,
    INFORMATIVE,
    MIN_RUNS_FOR_A_VERDICT,
    TOO_EASY,
    TOO_HARD,
    UNDECIDED,
    CalibrationError,
    calibrate,
    classify,
    load_run,
    load_runs,
    report,
)

TASK = "02_salinispora_bgc"


def _record(model, epoch, task_score, passes, soft=()):
    return {
        "task": TASK, "model": model, "epoch": epoch, "task_score": task_score,
        "artifacts": [
            {"path": path, "passed": passed, "soft": path in soft}
            for path, passed in passes.items()
        ],
    }


@pytest.fixture
def runs(tmp_path):
    """Four runs: two models, two epochs each.

    A passes 4/4, B 2/4, C 0/4, D 3/4.
    """
    records = [
        _record("m1", 1, 1.00, {"A": True, "B": True, "C": False, "D": True}),
        _record("m1", 2, 0.50, {"A": True, "B": False, "C": False, "D": True}),
        _record("m2", 1, 0.75, {"A": True, "B": True, "C": False, "D": True}),
        _record("m2", 2, 0.25, {"A": True, "B": False, "C": False, "D": False}),
    ]
    for index, record in enumerate(records):
        (tmp_path / f"run{index}.json").write_text(json.dumps(record))
    return load_runs(tmp_path)


class TestClassify:
    def test_too_few_runs_is_undecided(self):
        # A rate from 2 runs can only be 0.00, 0.50 or 1.00.
        for passed in range(MIN_RUNS_FOR_A_VERDICT):
            assert classify(MIN_RUNS_FOR_A_VERDICT - 1, passed) == UNDECIDED

    def test_everything_passes(self):
        assert classify(10, 10) == TOO_EASY

    def test_nothing_passes(self):
        assert classify(10, 0) == TOO_HARD

    def test_middle(self):
        assert classify(10, 5) == INFORMATIVE

    def test_the_band_edges_are_inclusive(self):
        # Exactly 0.80 and exactly 0.10 are kept: the band is "roughly 0.1 to
        # 0.8" in PROJECT.md section 8, so the boundary should not drop one.
        assert classify(10, 8) == INFORMATIVE
        assert classify(10, 1) == INFORMATIVE
        assert classify(10, 9) == TOO_EASY
        assert BAND_LOW == 0.10 and BAND_HIGH == 0.80


class TestPassRates:
    def test_rates_are_counted_per_artifact(self, runs):
        result = calibrate(runs)
        rates = {a.path: (a.n_runs, a.n_passed, round(a.pass_rate, 3))
                 for a in result.artifacts}
        assert rates == {"A": (4, 4, 1.0), "B": (4, 2, 0.5),
                         "C": (4, 0, 0.0), "D": (4, 3, 0.75)}

    def test_verdicts(self, runs):
        verdicts = {a.path: a.verdict for a in calibrate(runs).artifacts}
        assert verdicts == {"A": TOO_EASY, "B": INFORMATIVE,
                            "C": TOO_HARD, "D": INFORMATIVE}

    def test_soft_artifacts_are_not_tallied(self, tmp_path):
        # PROJECT.md section 9: a soft artifact is reported, never in the headline.
        (tmp_path / "r.json").write_text(json.dumps(
            _record("m1", 1, 1.0, {"A": True, "library_hits": False},
                    soft=("library_hits",))))
        result = calibrate(load_runs(tmp_path))
        assert [a.path for a in result.artifacts] == ["A"]


class TestModelScores:
    def test_mean_per_model_and_task(self, runs):
        scores = calibrate(runs).model_scores
        assert scores[("m1", TASK)] == pytest.approx(0.75)   # (1.00 + 0.50) / 2
        assert scores[("m2", TASK)] == pytest.approx(0.50)   # (0.75 + 0.25) / 2


class TestDrops:
    def test_only_uninformative_artifacts_are_dropped(self, runs):
        dropped = {a.path for a in calibrate(runs).drops()}
        assert dropped == {"A", "C"}

    def test_undecided_artifacts_are_never_dropped(self, tmp_path):
        for index in range(MIN_RUNS_FOR_A_VERDICT - 1):
            (tmp_path / f"r{index}.json").write_text(json.dumps(
                _record("m1", index, 1.0, {"A": True})))
        result = calibrate(load_runs(tmp_path))
        assert result.artifacts[0].verdict == UNDECIDED
        assert result.drops() == []

    def test_task_score_is_recomputed_over_the_kept_artifacts(self, runs):
        # A and C go; each run is rescored over {B, D} only:
        #   m1 epoch1 B+ D+ -> 1.00 ; m1 epoch2 B- D+ -> 0.50 ; mean 0.75
        #   m2 epoch1 B+ D+ -> 1.00 ; m2 epoch2 B- D- -> 0.00 ; mean 0.50
        after = calibrate(runs).task_scores_after_drops()
        assert after[("m1", TASK)] == pytest.approx(0.75)
        assert after[("m2", TASK)] == pytest.approx(0.50)

    def test_dropping_can_move_the_headline_score(self, tmp_path):
        # One run, four artifacts, three of which everything passes. As built
        # the score is 0.75; over the one informative artifact it is 0.00.
        records = [
            _record("m1", e, 0.75, {"easy1": True, "easy2": True, "easy3": True,
                                    "real": False})
            for e in range(MIN_RUNS_FOR_A_VERDICT)
        ]
        for index, record in enumerate(records):
            (tmp_path / f"r{index}.json").write_text(json.dumps(record))
        result = calibrate(load_runs(tmp_path))
        assert result.model_scores[("m1", TASK)] == pytest.approx(0.75)
        # every artifact is uninformative here, so nothing survives to rescore
        assert result.task_scores_after_drops() == {}


class TestLoading:
    def test_missing_field(self, tmp_path):
        (tmp_path / "r.json").write_text(json.dumps({"task": TASK, "model": "m"}))
        with pytest.raises(CalibrationError, match="missing required field 'artifacts'"):
            load_run(tmp_path / "r.json")

    def test_artifact_without_passed(self, tmp_path):
        (tmp_path / "r.json").write_text(json.dumps(
            {"task": TASK, "model": "m", "artifacts": [{"path": "A"}]}))
        with pytest.raises(CalibrationError, match="lacks 'path' or 'passed'"):
            load_run(tmp_path / "r.json")

    def test_duplicate_artifact(self, tmp_path):
        (tmp_path / "r.json").write_text(json.dumps(
            {"task": TASK, "model": "m",
             "artifacts": [{"path": "A", "passed": True}, {"path": "A", "passed": False}]}))
        with pytest.raises(CalibrationError, match="appears twice"):
            load_run(tmp_path / "r.json")

    def test_malformed_json(self, tmp_path):
        (tmp_path / "r.json").write_text("{")
        with pytest.raises(CalibrationError, match="not valid JSON"):
            load_run(tmp_path / "r.json")

    def test_empty_directory(self, tmp_path):
        with pytest.raises(CalibrationError, match="holds no run JSON files"):
            load_runs(tmp_path)


def test_report_names_the_drops_and_both_scores(runs):
    text = report(calibrate(runs))
    assert "too_easy" in text and "too_hard" in text
    assert "recommended drops" in text
    assert "after drops" in text


def test_report_warns_about_thin_evidence(tmp_path):
    (tmp_path / "r.json").write_text(json.dumps(_record("m1", 1, 1.0, {"A": True})))
    text = report(calibrate(load_runs(tmp_path)))
    assert "fewer than" in text and "not acted on" in text
