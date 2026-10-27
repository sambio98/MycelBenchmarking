"""Hand-computed fixtures for the four metrics.

Every expected value in this file is worked out by hand in the comment above
it. A grader that is plausible but subtly wrong is worse than no grader, so
nothing here is asserted against the implementation's own output.
"""

import math

import pytest

from npbench.grading.metrics import (
    TruthError,
    combine,
    identifier_f1,
    interval_overlap_f1,
    lins_ccc,
    macro_f1,
)
from npbench.grading.metrics import MetricResult


class TestIdentifierF1:
    def test_perfect(self):
        assert identifier_f1(["a", "b"], ["b", "a"]).score == 1.0

    def test_one_missing_one_extra(self):
        # TP=2 {a,b}, FP=1 {d}, FN=1 {c}. P=R=2/3, F1=2/3.
        assert identifier_f1(["a", "b", "c"], ["a", "b", "d"]).score == pytest.approx(2 / 3)

    def test_extra_rows_are_penalised_not_ignored(self):
        # An agent that emits every truth row plus five invented ones:
        # TP=2, FP=5, FN=0. P=2/7, R=1, F1=2*(2/7)/(2/7+1)=4/9.
        result = identifier_f1(["a", "b"], ["a", "b", "x", "y", "z", "w", "v"])
        assert result.score == pytest.approx(4 / 9)
        assert result.detail["n_extra"] == 5

    def test_empty_prediction(self):
        assert identifier_f1(["a", "b"], []).score == 0.0

    def test_both_empty(self):
        assert identifier_f1([], []).score == 1.0

    def test_disjoint(self):
        assert identifier_f1(["a"], ["b"]).score == 0.0


class TestLinsCCC:
    def test_perfect(self):
        assert lins_ccc([1, 2, 3], [1, 2, 3], log_transform=False, pseudocount=None).score == 1.0

    def test_constant_offset_scores_well_on_pearson_and_poorly_on_ccc(self):
        # The case CLAUDE.md names. t=[1,2,3,4], p=t+1.
        # mean_t=2.5, mean_p=3.5; var_t=var_p=cov=1.25 (population).
        # CCC = 2*1.25 / (1.25 + 1.25 + 1.0) = 2.5/3.5 = 0.7142857...
        # Pearson on the same input is exactly 1.0, which is why CCC is the metric.
        result = lins_ccc([1, 2, 3, 4], [2, 3, 4, 5], log_transform=False, pseudocount=None)
        assert result.score == pytest.approx(2.5 / 3.5)
        import numpy as np
        assert np.corrcoef([1, 2, 3, 4], [2, 3, 4, 5])[0, 1] == pytest.approx(1.0)

    def test_scale_error_is_penalised(self):
        # t=[1,2,3,4], p=2t=[2,4,6,8]. mean_t=2.5, mean_p=5.
        # var_t=1.25, var_p=5.0, cov=2.5.
        # CCC = 5 / (1.25 + 5 + 6.25) = 5/12.5 = 0.4
        result = lins_ccc([1, 2, 3, 4], [2, 4, 6, 8], log_transform=False, pseudocount=None)
        assert result.score == pytest.approx(0.4)

    def test_anticorrelation_clamps_to_zero(self):
        result = lins_ccc([1, 2, 3, 4], [4, 3, 2, 1], log_transform=False, pseudocount=None)
        assert result.score == 0.0
        assert result.detail["ccc"] < 0

    def test_log_transform(self):
        # log10 of [10,100,1000] vs [10,100,1000] is exact agreement.
        assert lins_ccc(
            [10, 100, 1000], [10, 100, 1000], log_transform=True, pseudocount=None
        ).score == pytest.approx(1.0)

    def test_fewer_than_two_aligned_rows(self):
        result = lins_ccc([1.0], [1.0], log_transform=False, pseudocount=None)
        assert result.score == 0.0
        assert result.detail["reason"] == "fewer_than_two_aligned_rows"

    def test_no_aligned_rows(self):
        assert lins_ccc([], [], log_transform=False, pseudocount=None).score == 0.0

    def test_non_finite_prediction_scores_zero(self):
        result = lins_ccc([1, 2, 3], [1, 2, float("nan")], log_transform=False, pseudocount=None)
        assert result.score == 0.0
        assert result.detail["reason"] == "non_finite_prediction"

    def test_non_positive_prediction_under_log_scores_zero(self):
        result = lins_ccc([1, 2, 3], [1, 2, 0], log_transform=True, pseudocount=None)
        assert result.score == 0.0
        assert result.detail["reason"] == "non_positive_prediction_under_log"

    def test_pseudocount_rescues_zeros(self):
        result = lins_ccc([0, 9, 99], [0, 9, 99], log_transform=True, pseudocount=1)
        assert result.score == pytest.approx(1.0)

    def test_constant_truth_column_is_a_bundle_bug(self):
        with pytest.raises(TruthError, match="constant"):
            lins_ccc([5, 5, 5], [5, 5, 5], log_transform=False, pseudocount=None)

    def test_non_finite_truth_is_a_bundle_bug(self):
        with pytest.raises(TruthError, match="non-finite"):
            lins_ccc([1, 2, float("nan")], [1, 2, 3], log_transform=False, pseudocount=None)


class TestMacroF1:
    def test_perfect(self):
        assert macro_f1(["nrps", "pks"], ["nrps", "pks"]).score == 1.0

    def test_minority_class_collapse(self):
        # truth [nrps,pks,nrps], pred all nrps.
        # nrps: TP=2 FP=1 FN=0 -> P=2/3 R=1 F1=0.8
        # pks:  TP=0 -> F1=0
        # macro = (0.8 + 0)/2 = 0.4
        # A micro-averaged score would read 2/3 here; macro is the metric because
        # a rare BGC class that the agent never predicts should be visible.
        result = macro_f1(["nrps", "pks", "nrps"], ["nrps", "nrps", "nrps"])
        assert result.score == pytest.approx(0.4)
        assert result.detail["per_class"] == {"nrps": pytest.approx(0.8), "pks": 0.0}

    def test_invented_class_is_penalised(self):
        # truth [a,a], pred [a,zzz].
        # a:   TP=1 FP=0 FN=1 -> P=1 R=0.5 F1=2/3
        # zzz: TP=0 -> F1=0
        # macro = (2/3)/2 = 1/3
        assert macro_f1(["a", "a"], ["a", "zzz"]).score == pytest.approx(1 / 3)

    def test_no_aligned_rows(self):
        result = macro_f1([], [])
        assert result.score == 0.0
        assert result.detail["reason"] == "no_aligned_rows"

    def test_length_mismatch_is_a_bundle_bug(self):
        with pytest.raises(TruthError):
            macro_f1(["a", "b"], ["a"])


class TestIntervalOverlapF1:
    def test_perfect(self):
        assert interval_overlap_f1(
            [("c1", 1, 100)], [("c1", 1, 100)], min_reciprocal_overlap=0.5
        ).score == 1.0

    def test_shifted_interval_above_threshold(self):
        # 1-based inclusive. truth 100-200 (len 101), pred 110-210 (len 101).
        # overlap = 200-110+1 = 91. reciprocal = 91/101 = 0.90099 >= 0.5 -> match.
        assert interval_overlap_f1(
            [("c1", 100, 200)], [("c1", 110, 210)], min_reciprocal_overlap=0.5
        ).score == 1.0

    def test_same_pair_below_a_stricter_threshold(self):
        # 0.90099 < 0.95, so no match: TP=0, FP=1, FN=1 -> F1=0.
        assert interval_overlap_f1(
            [("c1", 100, 200)], [("c1", 110, 210)], min_reciprocal_overlap=0.95
        ).score == 0.0

    def test_one_sprawling_prediction_cannot_claim_two_truth_intervals(self):
        # truth c1:1-100 and c1:101-200; a single pred c1:1-200.
        # reciprocal overlap with each is min(100/100, 100/200) = 0.5, so both
        # qualify at threshold 0.5 -- but matching is one-to-one.
        # TP=1, FP=0, FN=1. P=1, R=0.5, F1=2*1*0.5/1.5 = 2/3.
        result = interval_overlap_f1(
            [("c1", 1, 100), ("c1", 101, 200)], [("c1", 1, 200)], min_reciprocal_overlap=0.5
        )
        assert result.score == pytest.approx(2 / 3)
        assert result.detail["n_matched"] == 1

    def test_different_contigs_never_match(self):
        assert interval_overlap_f1(
            [("c1", 1, 100)], [("c2", 1, 100)], min_reciprocal_overlap=0.5
        ).score == 0.0

    def test_both_empty(self):
        assert interval_overlap_f1([], [], min_reciprocal_overlap=0.5).score == 1.0

    def test_empty_prediction(self):
        assert interval_overlap_f1(
            [("c1", 1, 100)], [], min_reciprocal_overlap=0.5
        ).score == 0.0

    def test_matching_is_deterministic_under_ties(self):
        truth = [("c1", 1, 100), ("c1", 1, 100)]
        preds = [("c1", 1, 100), ("c1", 1, 100)]
        first = interval_overlap_f1(truth, preds, min_reciprocal_overlap=0.5)
        for _ in range(5):
            assert interval_overlap_f1(truth, preds, min_reciprocal_overlap=0.5).score == first.score

    def test_inverted_truth_interval_is_a_bundle_bug(self):
        with pytest.raises(TruthError, match="start > end"):
            interval_overlap_f1([("c1", 200, 100)], [], min_reciprocal_overlap=0.5)

    def test_threshold_out_of_range_is_a_bundle_bug(self):
        with pytest.raises(TruthError, match="min_reciprocal_overlap"):
            interval_overlap_f1([("c1", 1, 2)], [], min_reciprocal_overlap=0.0)


class TestCombine:
    def test_product(self):
        results = [MetricResult("a", 0.9), MetricResult("b", 0.5)]
        assert combine(results) == pytest.approx(0.45)

    def test_one_zero_sinks_the_artifact(self):
        assert combine([MetricResult("a", 1.0), MetricResult("b", 0.0)]) == 0.0

    def test_no_metrics_is_a_bundle_bug(self):
        with pytest.raises(TruthError):
            combine([])
