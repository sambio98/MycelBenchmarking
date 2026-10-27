"""Artifact-level scoring: loading, alignment, degenerate submissions."""

import pytest

from npbench.grading.artifact import score_artifact
from npbench.grading.metrics import TruthError

TRUTH = """
    bgc_id,product_class,length_bp
    BGC0001,nrps,45000
    BGC0002,t1pks,32000
    BGC0003,nrps,51000
    BGC0004,terpene,21000
"""

BASE = {
    "path": "outputs/bgc.csv",
    "truth": "bgc.csv",
    "format": "csv",
    "soft": False,
    "row_key": ["bgc_id"],
    "normalisers": ["strip", "collapse_whitespace", "casefold"],
    "column_aliases": {},
    "value_aliases": {},
    "metrics": [{"type": "identifier_f1"}],
}


def _spec(spec_from, **overrides):
    return spec_from({**BASE, **overrides})


@pytest.fixture
def truth(write):
    return write("truth/bgc.csv", TRUTH)


def test_perfect_submission_scores_one(write, truth, spec_from):
    pred = write("outputs/bgc.csv", TRUTH)
    result = score_artifact(_spec(spec_from), truth, pred)
    assert result.score == 1.0
    assert result.passed
    assert result.format_compliant


def test_missing_file_scores_zero(write, truth, spec_from, tmp_path):
    result = score_artifact(_spec(spec_from), truth, tmp_path / "outputs/absent.csv")
    assert (result.score, result.passed, result.reason) == (0.0, False, "file_missing")


def test_empty_file_scores_zero(write, truth, spec_from):
    pred = write("outputs/bgc.csv", "")
    assert score_artifact(_spec(spec_from), truth, pred).reason == "file_empty"


def test_header_only_file_scores_zero(write, truth, spec_from):
    pred = write("outputs/bgc.csv", "bgc_id,product_class,length_bp\n")
    assert score_artifact(_spec(spec_from), truth, pred).reason == "no_data_rows"


def test_unreadable_file_scores_zero(write, truth, spec_from):
    # Ragged rows: a real symptom of an agent writing a table by hand.
    pred = write("outputs/bgc.csv", "bgc_id,product_class\nBGC0001,nrps,45000,extra\n")
    assert score_artifact(_spec(spec_from), truth, pred).score == 0.0


def test_all_rows_missing_scores_zero(write, truth, spec_from):
    pred = write(
        "outputs/bgc.csv",
        """
        bgc_id,product_class,length_bp
        BGC9001,nrps,45000
        """,
    )
    assert score_artifact(_spec(spec_from), truth, pred).score == 0.0


def test_duplicate_row_keys_fail_the_artifact(write, truth, spec_from):
    # An inflated inner join: every truth row present, one of them twice.
    pred = write(
        "outputs/bgc.csv",
        TRUTH + "    BGC0001,nrps,45000\n",
    )
    result = score_artifact(_spec(spec_from), truth, pred)
    assert (result.score, result.reason) == (0.0, "duplicate_row_keys")


def test_duplicate_row_keys_in_truth_is_a_bundle_bug(write, spec_from):
    truth = write("truth/bgc.csv", TRUTH + "    BGC0001,nrps,45000\n")
    pred = write("outputs/bgc.csv", TRUTH)
    with pytest.raises(TruthError, match="duplicate row keys in truth"):
        score_artifact(_spec(spec_from), truth, pred)


def test_missing_declared_column_scores_zero(write, truth, spec_from):
    pred = write("outputs/bgc.csv", "cluster,product_class\nBGC0001,nrps\n")
    result = score_artifact(_spec(spec_from), truth, pred)
    assert result.score == 0.0
    assert "missing_column" in result.reason


class TestNormalisationAndAliases:
    def test_case_and_whitespace_differences_do_not_cost_score(self, write, truth, spec_from):
        pred = write(
            "outputs/bgc.csv",
            """
            bgc_id,product_class,length_bp
            bgc0001 ,nrps,45000
            BGC0002,t1pks,32000
            BGC0003 ,nrps,51000
            BGC0004,terpene,21000
            """,
        )
        assert score_artifact(_spec(spec_from), truth, pred).score == 1.0

    def test_renamed_column_scores_but_is_recorded_as_non_compliant(
        self, write, truth, spec_from
    ):
        # The BixBench3 failure mode of PROJECT.md section 3: correct science,
        # wrong column name. It must not cost artifact score.
        pred = write(
            "outputs/bgc.csv",
            """
            BGC_ID,product_class,length_bp
            BGC0001,nrps,45000
            BGC0002,t1pks,32000
            BGC0003,nrps,51000
            BGC0004,terpene,21000
            """,
        )
        spec = _spec(spec_from, column_aliases={"bgc_id": ["BGC_ID", "cluster_id"]})
        result = score_artifact(spec, truth, pred)
        assert result.score == 1.0
        assert not result.format_compliant
        assert result.renamed_columns == {"bgc_id": "BGC_ID"}

    def test_value_alias_resolves_a_label_synonym(self, write, truth, spec_from):
        pred = write(
            "outputs/bgc.csv",
            """
            bgc_id,product_class,length_bp
            BGC0001,NRPS,45000
            BGC0002,Type I PKS,32000
            BGC0003,NRPS,51000
            BGC0004,terpene,21000
            """,
        )
        spec = _spec(
            spec_from,
            metrics=[{"type": "identifier_f1"}, {"type": "macro_f1", "column": "product_class"}],
            value_aliases={"product_class": {"Type I PKS": "t1pks"}},
        )
        assert score_artifact(spec, truth, pred).score == pytest.approx(1.0)


class TestMetricComposition:
    CCC_METRICS = [
        {"type": "identifier_f1"},
        {"type": "lins_ccc", "column": "length_bp", "log_transform": False, "pseudocount": None},
    ]

    def test_product_of_two_metrics(self, write, truth, spec_from):
        # Three of four rows present and correct, one missing.
        # identifier_f1: TP=3 FP=0 FN=1 -> P=1 R=0.75 F1=6/7
        # macro_f1 over the 3 aligned rows, all labels correct:
        #   classes seen are nrps and t1pks, both F1=1 -> 1.0
        # product = 6/7 = 0.857 -> passes at 0.80
        pred = write(
            "outputs/bgc.csv",
            """
            bgc_id,product_class,length_bp
            BGC0001,nrps,45000
            BGC0002,t1pks,32000
            BGC0003,nrps,51000
            """,
        )
        spec = _spec(
            spec_from,
            metrics=[{"type": "identifier_f1"}, {"type": "macro_f1", "column": "product_class"}],
        )
        result = score_artifact(spec, truth, pred)
        assert result.score == pytest.approx(6 / 7)
        assert result.passed

    def test_identifier_f1_carries_the_penalty_for_dropped_rows(self, write, truth, spec_from):
        # Only one row survives, but its value is right. CCC sees one aligned
        # row and returns 0, so the artifact fails. This is why the validator
        # insists identifier_f1 rides alongside the aligned metrics.
        pred = write("outputs/bgc.csv", "bgc_id,product_class,length_bp\nBGC0001,nrps,45000\n")
        assert score_artifact(_spec(spec_from, metrics=self.CCC_METRICS), truth, pred).score == 0.0

    def test_small_offset_is_tolerated_because_the_spread_dominates(
        self, write, truth, spec_from
    ):
        # Every length 1000 bp too large. CCC = 0.996, so the artifact passes.
        # That is the metric behaving correctly: against a 21-51 kb spread a
        # 1 kb offset is a 2% error. CCC only bites when the offset is
        # comparable to the spread -- see test_metrics.py, where an offset of 1
        # on [1,2,3,4] drops CCC to 0.71 while Pearson stays at 1.0.
        # Consequence for task design: a numeric column whose between-row
        # variance swamps the error an agent would plausibly make is a weak
        # artifact. Pick quantities where a method error moves the rank order.
        pred = write(
            "outputs/bgc.csv",
            """
            bgc_id,product_class,length_bp
            BGC0001,nrps,46000
            BGC0002,t1pks,33000
            BGC0003,nrps,52000
            BGC0004,terpene,22000
            """,
        )
        result = score_artifact(_spec(spec_from, metrics=self.CCC_METRICS), truth, pred)
        assert result.score == pytest.approx(0.99632, abs=1e-5)
        assert result.passed

    def test_scale_error_fails_the_artifact(self, write, truth, spec_from):
        # Every length doubled: the kind of unit or normalisation slip CCC exists
        # to catch. Pearson would read 1.0. Worked by hand:
        #   t = [45000, 32000, 51000, 21000], p = 2t
        #   var_t = 135,187,500; var_p = 540,750,000; cov = 270,375,000
        #   mean difference = 37,250, squared = 1,387,562,500
        #   CCC = 2*270375000 / (135187500 + 540750000 + 1387562500) = 0.26206
        pred = write(
            "outputs/bgc.csv",
            """
            bgc_id,product_class,length_bp
            BGC0001,nrps,90000
            BGC0002,t1pks,64000
            BGC0003,nrps,102000
            BGC0004,terpene,42000
            """,
        )
        result = score_artifact(_spec(spec_from, metrics=self.CCC_METRICS), truth, pred)
        assert result.score == pytest.approx(0.26206, abs=1e-5)
        assert not result.passed

    def test_rows_are_aligned_by_key_not_by_file_order(self, write, truth, spec_from):
        pred = write(
            "outputs/bgc.csv",
            """
            bgc_id,product_class,length_bp
            BGC0004,terpene,21000
            BGC0002,t1pks,32000
            BGC0003,nrps,51000
            BGC0001,nrps,45000
            """,
        )
        spec = _spec(spec_from, metrics=self.CCC_METRICS)
        assert score_artifact(spec, truth, pred).score == pytest.approx(1.0)


class TestIntervals:
    TRUTH_INTERVALS = """
        bgc_id,contig,start,end
        BGC0001,contig_1,1000,46000
        BGC0002,contig_2,5000,37000
    """
    METRICS = [
        {"type": "identifier_f1"},
        {
            "type": "interval_overlap_f1",
            "group_column": "contig",
            "start_column": "start",
            "end_column": "end",
            "min_reciprocal_overlap": 0.5,
        },
    ]

    def test_perfect_coordinates(self, write, spec_from):
        truth = write("truth/i.csv", self.TRUTH_INTERVALS)
        pred = write("outputs/i.csv", self.TRUTH_INTERVALS)
        spec = _spec(spec_from, truth="i.csv", path="outputs/i.csv", metrics=self.METRICS)
        assert score_artifact(spec, truth, pred).score == 1.0

    def test_non_integer_coordinate_in_submission_scores_zero(self, write, spec_from):
        truth = write("truth/i.csv", self.TRUTH_INTERVALS)
        pred = write(
            "outputs/i.csv",
            """
            bgc_id,contig,start,end
            BGC0001,contig_1,1000,46000
            BGC0002,contig_2,5000,~37000
            """,
        )
        spec = _spec(spec_from, truth="i.csv", path="outputs/i.csv", metrics=self.METRICS)
        result = score_artifact(spec, truth, pred)
        assert result.score == 0.0
        assert "non_integer_coordinate" in result.reason


def test_composite_row_key(write, spec_from):
    truth = write(
        "truth/link.csv",
        """
        strain,feature_id,score
        CNB-440,FT001,0.91
        CNB-440,FT002,0.44
        CNS-205,FT001,0.78
        """,
    )
    # Same strains, same features, but only the pairs make a unique key.
    pred = write(
        "outputs/link.csv",
        """
        strain,feature_id,score
        CNB-440,FT001,0.91
        CNS-205,FT001,0.78
        CNB-440,FT002,0.44
        """,
    )
    spec = _spec(
        spec_from,
        path="outputs/link.csv",
        truth="link.csv",
        row_key=["strain", "feature_id"],
        metrics=[
            {"type": "identifier_f1"},
            {"type": "lins_ccc", "column": "score", "log_transform": False, "pseudocount": None},
        ],
    )
    assert score_artifact(spec, truth, pred).score == pytest.approx(1.0)
