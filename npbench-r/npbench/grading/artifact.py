"""Score one artifact: load the submission, align it to truth, apply the metrics.

Two rules here are judgement calls worth stating.

Duplicate row keys in a submission fail the artifact outright. Row alignment is
undefined when a key appears twice, and silently keeping the first row would
reward exactly the failure PROJECT.md §8 wants to catch: an agent that
inner-joins badly and inflates its row count. A duplicate key in *truth* raises
instead, because that is our bug.

Aligned metrics (CCC, macro F1) see only the truth/prediction intersection.
They cannot notice dropped rows on their own, so the validator requires
identifier_f1 alongside them and the product combiner carries the penalty.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import metrics as M
from .config import PASS_THRESHOLD, ArtifactSpec
from .normalise import build_alias_map, normalise_all

KEY_SEPARATOR = "\x1f"


@dataclass
class ArtifactScore:
    path: str
    soft: bool
    score: float
    passed: bool
    metrics: list[M.MetricResult] = field(default_factory=list)
    renamed_columns: dict[str, str] = field(default_factory=dict)
    reason: str | None = None

    @property
    def format_compliant(self) -> bool:
        """True when every declared column was found under its declared name."""
        return not self.renamed_columns


class SubmissionError(Exception):
    """The submission cannot be scored. Always becomes a 0.00, never a crash."""


def read_table(path: Path, separator: str) -> pd.DataFrame:
    if not path.exists():
        raise SubmissionError("file_missing")
    if path.stat().st_size == 0:
        raise SubmissionError("file_empty")
    try:
        frame = pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False)
    except Exception as exc:  # unparseable: wrong delimiter, ragged rows, binary
        raise SubmissionError(f"file_unreadable: {type(exc).__name__}") from exc
    if frame.empty:
        raise SubmissionError("no_data_rows")
    return frame


def resolve_columns(
    frame: pd.DataFrame, wanted: list[str], spec: ArtifactSpec
) -> tuple[dict[str, str], dict[str, str]]:
    """Map each declared column name onto a header actually present.

    Headers are matched case- and whitespace-insensitively, then through the
    declared alias list. A column found under anything but its declared name is
    recorded as a rename, which feeds the format-compliance axis rather than
    the artifact score (PROJECT.md §7).
    """
    present = {" ".join(str(c).split()).casefold(): str(c) for c in frame.columns}
    resolved: dict[str, str] = {}
    renamed: dict[str, str] = {}
    missing: list[str] = []
    for column in wanted:
        for candidate in [column, *spec.column_aliases.get(column, [])]:
            actual = present.get(" ".join(candidate.split()).casefold())
            if actual is not None:
                resolved[column] = actual
                if actual != column:
                    renamed[column] = actual
                break
        else:
            missing.append(column)
    if missing:
        raise SubmissionError(f"missing_column(s): {', '.join(missing)}")
    return resolved, renamed


def _row_keys(frame: pd.DataFrame, spec: ArtifactSpec, resolved: dict[str, str]) -> list[str]:
    columns = [
        normalise_all(frame[resolved[c]], spec.rules_for(c)) for c in spec.row_key
    ]
    return [KEY_SEPARATOR.join(parts) for parts in zip(*columns)]


def _values(frame: pd.DataFrame, header: str, column: str, spec: ArtifactSpec) -> list[str]:
    """Values of `header`, with the alias map declared for canonical name `column`."""
    rules = spec.rules_for(column)
    aliases = build_alias_map(spec.value_aliases.get(column, {}), rules)
    return normalise_all(frame[header], rules, aliases)


def _intervals(
    frame: pd.DataFrame,
    spec: ArtifactSpec,
    resolved: dict[str, str],
    params: dict,
    side: str,
) -> list[tuple[str, int, int]]:
    group_column = params["group_column"]
    groups = normalise_all(frame[resolved[group_column]], spec.rules_for(group_column))
    out = []
    for group, start, end in zip(
        groups, frame[resolved[params["start_column"]]], frame[resolved[params["end_column"]]]
    ):
        try:
            out.append((group, int(str(start).strip()), int(str(end).strip())))
        except ValueError as exc:
            message = f"non_integer_coordinate: {start!r}-{end!r}"
            if side == "truth":
                raise M.TruthError(message) from exc
            raise SubmissionError(message) from exc
    return out


def declared_columns(spec: ArtifactSpec) -> list[str]:
    """Every column the spec needs, in a stable order, without duplicates."""
    wanted = list(spec.row_key)
    for metric in spec.metrics:
        for key in ("column", "group_column", "start_column", "end_column"):
            if key in metric.params and metric.params[key] not in wanted:
                wanted.append(metric.params[key])
    return wanted


def score_artifact(spec: ArtifactSpec, truth_path: Path, pred_path: Path) -> ArtifactScore:
    """Grade one submitted artifact against its published counterpart."""
    wanted = declared_columns(spec)

    # Truth-side problems raise. They mean the bundle is broken, not the agent.
    try:
        truth = read_table(truth_path, spec.separator)
        truth_columns, _ = resolve_columns(truth, wanted, spec)
    except SubmissionError as exc:
        raise M.TruthError(f"{truth_path}: {exc}") from exc
    truth_keys = _row_keys(truth, spec, truth_columns)
    if len(set(truth_keys)) != len(truth_keys):
        raise M.TruthError(f"{truth_path}: duplicate row keys in truth")

    try:
        pred = read_table(pred_path, spec.separator)
        pred_columns, renamed = resolve_columns(pred, wanted, spec)
        pred_keys = _row_keys(pred, spec, pred_columns)
        if len(set(pred_keys)) != len(pred_keys):
            raise SubmissionError("duplicate_row_keys")
        results = _apply_metrics(spec, truth, pred, truth_columns, pred_columns,
                                 truth_keys, pred_keys)
    except SubmissionError as exc:
        return ArtifactScore(spec.path, spec.soft, 0.0, False, reason=str(exc))

    score = M.combine(results)
    return ArtifactScore(
        spec.path, spec.soft, score, score >= PASS_THRESHOLD, results, renamed
    )


def _apply_metrics(
    spec: ArtifactSpec,
    truth: pd.DataFrame,
    pred: pd.DataFrame,
    truth_columns: dict[str, str],
    pred_columns: dict[str, str],
    truth_keys: list[str],
    pred_keys: list[str],
) -> list[M.MetricResult]:
    # Rows shared by both sides, in truth order, so aligned metrics are deterministic.
    pred_row = {key: i for i, key in enumerate(pred_keys)}
    shared = [(i, pred_row[key]) for i, key in enumerate(truth_keys) if key in pred_row]

    results = []
    for metric in spec.metrics:
        params = metric.params
        if metric.type == "identifier_f1":
            results.append(M.identifier_f1(truth_keys, pred_keys))
        elif metric.type == "interval_overlap_f1":
            results.append(
                M.interval_overlap_f1(
                    _intervals(truth, spec, truth_columns, params, "truth"),
                    _intervals(pred, spec, pred_columns, params, "pred"),
                    min_reciprocal_overlap=params["min_reciprocal_overlap"],
                )
            )
        elif metric.type == "macro_f1":
            column = params["column"]
            t = _values(truth, truth_columns[column], column, spec)
            p = _values(pred, pred_columns[column], column, spec)
            results.append(M.macro_f1([t[i] for i, _ in shared], [p[j] for _, j in shared]))
        elif metric.type == "lins_ccc":
            column = params["column"]
            t = pd.to_numeric(truth[truth_columns[column]], errors="coerce")
            p = pd.to_numeric(pred[pred_columns[column]], errors="coerce")
            results.append(
                M.lins_ccc(
                    [t.iloc[i] for i, _ in shared],
                    [p.iloc[j] for _, j in shared],
                    log_transform=params["log_transform"],
                    pseudocount=params["pseudocount"],
                )
            )
        else:  # pragma: no cover - config.py rejects unknown types
            raise M.TruthError(f"unhandled metric type {metric.type!r}")
    return results
