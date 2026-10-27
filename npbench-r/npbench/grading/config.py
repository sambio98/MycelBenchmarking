"""grading.yaml parsing.

Strict by design. Every per-artifact parameter the spec does not fix — the
reciprocal overlap threshold, the log transform, the normalisation rules, the
alias maps — must be stated in the file. There are no defaults, so a score can
never be computed under a threshold nobody chose.

The 0.80 pass threshold is the one constant that is *not* per-task. It is a
property of the benchmark (PROJECT.md §7), to be re-derived once by expert
calibration, and tasks must not quietly differ on it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from .normalise import RULES

PASS_THRESHOLD = 0.80

TABLE_FORMATS = {"csv": ",", "tsv": "\t"}


class ConfigError(ValueError):
    """grading.yaml does not satisfy the contract."""


def _has_type(value: object, types: tuple[type, ...]) -> bool:
    """isinstance, except that a bool never satisfies int or float.

    `soft: 1` and `min_reciprocal_overlap: true` are both mistakes worth catching,
    and bare isinstance accepts them because bool subclasses int.
    """
    if isinstance(value, bool) and bool not in types:
        return False
    return isinstance(value, types)


def _type_names(types: tuple[type, ...]) -> str:
    return " or ".join(t.__name__ for t in types)


def _require(mapping: dict, key: str, where: str, types: tuple[type, ...]) -> object:
    if key not in mapping:
        raise ConfigError(f"{where}: missing required field {key!r}")
    value = mapping[key]
    if not _has_type(value, types):
        raise ConfigError(
            f"{where}: field {key!r} must be {_type_names(types)}, got {type(value).__name__}"
        )
    return value


@dataclass(frozen=True)
class MetricSpec:
    type: str
    params: dict


@dataclass(frozen=True)
class ArtifactSpec:
    path: str                       # where the agent writes it, relative to outputs/
    truth: str                      # the published artifact, relative to truth/
    format: str                     # csv | tsv
    soft: bool                      # excluded from the headline score (PROJECT.md §9)
    row_key: tuple[str, ...]
    normalisers: tuple[str, ...]
    column_normalisers: dict[str, tuple[str, ...]]
    column_aliases: dict[str, list[str]]
    value_aliases: dict[str, dict[str, str]]
    metrics: tuple[MetricSpec, ...]

    @property
    def separator(self) -> str:
        return TABLE_FORMATS[self.format]

    def rules_for(self, column: str) -> tuple[str, ...]:
        """Normalisation rules for one column, falling back to the artifact's.

        Needed because case is meaningful in some columns and noise in others:
        in a SMILES string `c` is an aromatic carbon and `C` an aliphatic one,
        so casefolding a structure changes the molecule, while casefolding a
        compound name only forgives a formatting slip.
        """
        return self.column_normalisers.get(column, self.normalisers)


@dataclass(frozen=True)
class GradingConfig:
    version: int
    artifacts: tuple[ArtifactSpec, ...]


# Required parameters per metric type. Anything else in the metric block is an error,
# which catches a misspelled parameter rather than silently ignoring it.
METRIC_PARAMS: dict[str, dict[str, tuple[type, ...]]] = {
    "identifier_f1": {},
    "lins_ccc": {"column": (str,), "log_transform": (bool,), "pseudocount": (int, float, type(None))},
    "macro_f1": {"column": (str,)},
    "interval_overlap_f1": {
        "group_column": (str,),
        "start_column": (str,),
        "end_column": (str,),
        "min_reciprocal_overlap": (int, float),
    },
}

# Metrics that compare rows pairwise, so they only see the truth/prediction
# intersection. On their own they would not notice a prediction that dropped
# half the rows, so the validator insists identifier_f1 rides alongside and the
# product combiner carries the penalty.
ALIGNED_METRICS = frozenset({"lins_ccc", "macro_f1"})


def _parse_metric(block: object, where: str) -> MetricSpec:
    if not isinstance(block, dict):
        raise ConfigError(f"{where}: each metric must be a mapping")
    kind = _require(block, "type", where, (str,))
    if kind not in METRIC_PARAMS:
        raise ConfigError(
            f"{where}: unknown metric type {kind!r}; known types: "
            f"{', '.join(sorted(METRIC_PARAMS))}"
        )
    where = f"{where} metric {kind!r}"
    required = METRIC_PARAMS[kind]
    unexpected = set(block) - set(required) - {"type"}
    if unexpected:
        raise ConfigError(f"{where}: unexpected parameter(s) {', '.join(sorted(unexpected))}")
    params = {}
    for param, types in required.items():
        if param not in block:
            raise ConfigError(
                f"{where}: missing required parameter {param!r}; this metric has no default"
            )
        value = block[param]
        if not _has_type(value, types):
            raise ConfigError(
                f"{where}: parameter {param!r} must be {_type_names(types)}, "
                f"got {type(value).__name__}"
            )
        params[param] = value
    return MetricSpec(kind, params)


def _parse_artifact(block: object, index: int) -> ArtifactSpec:
    where = f"artifact {index}"
    if not isinstance(block, dict):
        raise ConfigError(f"{where}: each artifact must be a mapping")
    path = _require(block, "path", where, (str,))
    where = f"artifact {path!r}"

    fmt = _require(block, "format", where, (str,))
    if fmt not in TABLE_FORMATS:
        raise ConfigError(
            f"{where}: format must be one of {', '.join(sorted(TABLE_FORMATS))}, got {fmt!r}"
        )

    row_key = _require(block, "row_key", where, (list,))
    if not row_key or not all(isinstance(c, str) for c in row_key):
        raise ConfigError(f"{where}: row_key must be a non-empty list of column names")

    normalisers = _require(block, "normalisers", where, (list,))
    unknown = [r for r in normalisers if r not in RULES]
    if unknown:
        raise ConfigError(
            f"{where}: unknown normalisation rule(s) {', '.join(map(repr, unknown))}; "
            f"known rules: {', '.join(RULES)}"
        )

    overrides = block.get("column_normalisers") or {}
    if not isinstance(overrides, dict):
        raise ConfigError(f"{where}: column_normalisers must be a mapping")
    column_normalisers = {}
    for column, rules in overrides.items():
        if not isinstance(rules, list):
            raise ConfigError(f"{where}: column_normalisers[{column!r}] must be a list")
        unknown = [r for r in rules if r not in RULES]
        if unknown:
            raise ConfigError(
                f"{where}: column_normalisers[{column!r}] has unknown rule(s) "
                f"{', '.join(map(repr, unknown))}; known rules: {', '.join(RULES)}"
            )
        column_normalisers[column] = tuple(rules)

    column_aliases = _require(block, "column_aliases", where, (dict,))
    value_aliases = _require(block, "value_aliases", where, (dict,))
    for column, mapping in value_aliases.items():
        if not isinstance(mapping, dict):
            raise ConfigError(f"{where}: value_aliases[{column!r}] must be a mapping")

    metrics = _require(block, "metrics", where, (list,))
    if not metrics:
        raise ConfigError(f"{where}: must declare at least one metric")

    return ArtifactSpec(
        path=path,
        truth=_require(block, "truth", where, (str,)),
        format=fmt,
        soft=_require(block, "soft", where, (bool,)),
        row_key=tuple(row_key),
        normalisers=tuple(normalisers),
        column_normalisers=column_normalisers,
        column_aliases={k: list(v) for k, v in column_aliases.items()},
        value_aliases={k: dict(v) for k, v in value_aliases.items()},
        metrics=tuple(_parse_metric(m, where) for m in metrics),
    )


def load_grading_config(path: Path) -> GradingConfig:
    try:
        raw = yaml.safe_load(path.read_text())
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path}: not valid YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: top level must be a mapping")

    version = _require(raw, "version", str(path), (int,))
    if version != 1:
        raise ConfigError(f"{path}: unsupported grading.yaml version {version}")

    artifacts = _require(raw, "artifacts", str(path), (list,))
    if not artifacts:
        raise ConfigError(f"{path}: declares no artifacts")

    parsed = tuple(_parse_artifact(a, i) for i, a in enumerate(artifacts))
    seen = [a.path for a in parsed]
    duplicates = {p for p in seen if seen.count(p) > 1}
    if duplicates:
        raise ConfigError(f"{path}: duplicate artifact path(s) {', '.join(sorted(duplicates))}")
    return GradingConfig(version=version, artifacts=parsed)
