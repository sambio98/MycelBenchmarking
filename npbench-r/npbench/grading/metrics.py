"""The four artifact metrics of PROJECT.md §7.

Every metric returns a score in [0, 1] plus a detail dict. Where a metric is
undefined on its input it returns 0.0 and names the reason, never NaN: a NaN
would propagate through the product combiner and make an artifact unscoreable
rather than failed.

Truth-side degeneracy is a different matter. If our own truth file is empty of
variance, holds a non-finite number, or is otherwise unusable, that is a bug in
the task bundle, not an agent failure, so it raises.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class MetricResult:
    name: str
    score: float
    detail: dict = field(default_factory=dict)


class TruthError(ValueError):
    """The truth side of a comparison is unusable. A bundle bug, not a failure."""


def _f1(true_positive: int, false_positive: int, false_negative: int) -> float:
    if true_positive == 0:
        return 0.0
    precision = true_positive / (true_positive + false_positive)
    recall = true_positive / (true_positive + false_negative)
    return 2 * precision * recall / (precision + recall)


def identifier_f1(truth_keys: Sequence[str], pred_keys: Sequence[str]) -> MetricResult:
    """F1 over the row-key universe, penalising missing and extra entries alike.

    Called with already-normalised keys. Duplicates are rejected upstream in
    `artifact.py`, so these are treated as sets.
    """
    truth, pred = set(truth_keys), set(pred_keys)
    if not truth and not pred:
        return MetricResult("identifier_f1", 1.0, {"reason": "both_empty"})
    shared = truth & pred
    score = _f1(len(shared), len(pred - truth), len(truth - pred))
    return MetricResult(
        "identifier_f1",
        score,
        {
            "n_truth": len(truth),
            "n_pred": len(pred),
            "n_shared": len(shared),
            "n_missing": len(truth - pred),
            "n_extra": len(pred - truth),
        },
    )


def lins_ccc(
    truth: Sequence[float],
    pred: Sequence[float],
    *,
    log_transform: bool,
    pseudocount: float | None,
) -> MetricResult:
    """Lin's concordance correlation coefficient over aligned rows.

    CCC = 2*cov / (var_t + var_p + (mean_t - mean_p)^2), population moments.

    Chosen over Pearson because it penalises a changed slope or offset: a
    normalisation error that preserves rank order is a real error (PROJECT.md
    §7). The result is clamped to [0, 1] for the product combiner, since
    negative concordance is not better than none. The raw coefficient is kept
    in the detail dict.
    """
    name = "lins_ccc"
    t = np.asarray(truth, dtype=float)
    p = np.asarray(pred, dtype=float)
    if t.shape != p.shape:
        raise TruthError(f"aligned arrays differ in length: {t.shape} vs {p.shape}")

    if not np.all(np.isfinite(t)):
        raise TruthError("truth column holds a non-finite value")
    if not np.all(np.isfinite(p)):
        return MetricResult(name, 0.0, {"reason": "non_finite_prediction"})
    if t.size < 2:
        return MetricResult(name, 0.0, {"reason": "fewer_than_two_aligned_rows", "n": int(t.size)})

    if log_transform:
        offset = 0.0 if pseudocount is None else pseudocount
        if np.any(t + offset <= 0):
            raise TruthError("truth column holds a non-positive value under log transform")
        if np.any(p + offset <= 0):
            return MetricResult(name, 0.0, {"reason": "non_positive_prediction_under_log"})
        t, p = np.log10(t + offset), np.log10(p + offset)

    var_t, var_p = t.var(ddof=0), p.var(ddof=0)
    if var_t == 0:
        raise TruthError(
            "truth column is constant, so CCC is undefined; this artifact needs a "
            "different metric"
        )
    covariance = float(((t - t.mean()) * (p - p.mean())).mean())
    denominator = var_t + var_p + (t.mean() - p.mean()) ** 2
    ccc = 2 * covariance / float(denominator)
    return MetricResult(
        name,
        float(max(0.0, min(1.0, ccc))),
        {"ccc": float(ccc), "n": int(t.size), "log_transform": log_transform},
    )


def macro_f1(truth_labels: Sequence[str], pred_labels: Sequence[str]) -> MetricResult:
    """Macro F1 across label classes over aligned rows.

    The class universe is the union of labels seen in truth and in the
    prediction, so an invented class scores 0 on its own row and drags the
    average down rather than being ignored.
    """
    name = "macro_f1"
    if len(truth_labels) != len(pred_labels):
        raise TruthError(
            f"aligned label sequences differ in length: {len(truth_labels)} vs {len(pred_labels)}"
        )
    if not truth_labels:
        return MetricResult(name, 0.0, {"reason": "no_aligned_rows"})

    classes = sorted(set(truth_labels) | set(pred_labels))
    per_class = {}
    for label in classes:
        tp = sum(1 for t, p in zip(truth_labels, pred_labels) if t == label and p == label)
        fp = sum(1 for t, p in zip(truth_labels, pred_labels) if t != label and p == label)
        fn = sum(1 for t, p in zip(truth_labels, pred_labels) if t == label and p != label)
        per_class[label] = _f1(tp, fp, fn)
    score = sum(per_class.values()) / len(per_class)
    return MetricResult(
        name, score, {"n": len(truth_labels), "n_classes": len(classes), "per_class": per_class}
    )


def _reciprocal_overlap(a: tuple[int, int], b: tuple[int, int]) -> float:
    """Smaller of the two overlap fractions. Coordinates are 1-based inclusive."""
    overlap = min(a[1], b[1]) - max(a[0], b[0]) + 1
    if overlap <= 0:
        return 0.0
    return min(overlap / (a[1] - a[0] + 1), overlap / (b[1] - b[0] + 1))


def interval_overlap_f1(
    truth_intervals: Sequence[tuple[str, int, int]],
    pred_intervals: Sequence[tuple[str, int, int]],
    *,
    min_reciprocal_overlap: float,
) -> MetricResult:
    """Overlap F1 over (group, start, end) intervals under a reciprocal overlap rule.

    Used for BGC coordinates on contigs, where `group` is the contig. Matching
    is one-to-one and greedy on descending reciprocal overlap, so one sprawling
    predicted interval cannot claim credit for several truth intervals. Ties
    break on the interval tuples, which makes the result deterministic.
    """
    name = "interval_overlap_f1"
    if not 0 < min_reciprocal_overlap <= 1:
        raise TruthError(
            f"min_reciprocal_overlap must be in (0, 1], got {min_reciprocal_overlap}"
        )
    for group, start, end in truth_intervals:
        if start > end:
            raise TruthError(f"truth interval {group}:{start}-{end} has start > end")
    if not truth_intervals and not pred_intervals:
        return MetricResult(name, 1.0, {"reason": "both_empty"})

    candidates = []
    for i, t in enumerate(truth_intervals):
        for j, p in enumerate(pred_intervals):
            if t[0] != p[0] or p[1] > p[2]:
                continue
            overlap = _reciprocal_overlap((t[1], t[2]), (p[1], p[2]))
            if overlap >= min_reciprocal_overlap:
                candidates.append((-overlap, t, p, i, j))
    candidates.sort()

    matched_truth: set[int] = set()
    matched_pred: set[int] = set()
    for _, _, _, i, j in candidates:
        if i not in matched_truth and j not in matched_pred:
            matched_truth.add(i)
            matched_pred.add(j)

    tp = len(matched_truth)
    score = _f1(tp, len(pred_intervals) - tp, len(truth_intervals) - tp)
    return MetricResult(
        name,
        score,
        {
            "n_truth": len(truth_intervals),
            "n_pred": len(pred_intervals),
            "n_matched": tp,
            "min_reciprocal_overlap": min_reciprocal_overlap,
        },
    )


def combine(results: Sequence[MetricResult]) -> float:
    """Multiply the metrics applicable to one artifact into a single score."""
    if not results:
        raise TruthError("an artifact must declare at least one metric")
    return math.prod(r.score for r in results)
