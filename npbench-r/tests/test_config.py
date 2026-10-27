"""grading.yaml is strict: every parameter the spec does not fix must be stated."""

import pytest
import yaml

from npbench.grading.config import ConfigError, load_grading_config

VALID = {
    "version": 1,
    "artifacts": [
        {
            "path": "outputs/bgc.csv",
            "truth": "bgc.csv",
            "format": "csv",
            "soft": False,
            "row_key": ["bgc_id"],
            "normalisers": ["strip", "casefold"],
            "column_aliases": {"bgc_id": ["BGC_ID"]},
            "value_aliases": {"product_class": {"Type I PKS": "t1pks"}},
            "metrics": [
                {"type": "identifier_f1"},
                {"type": "macro_f1", "column": "product_class"},
            ],
        }
    ],
}


@pytest.fixture
def load(tmp_path):
    def _load(document):
        path = tmp_path / "grading.yaml"
        path.write_text(yaml.safe_dump(document))
        return load_grading_config(path)

    return _load


def test_valid_config_round_trips(load):
    config = load(VALID)
    artifact = config.artifacts[0]
    assert artifact.path == "outputs/bgc.csv"
    assert artifact.row_key == ("bgc_id",)
    assert artifact.separator == ","
    assert [m.type for m in artifact.metrics] == ["identifier_f1", "macro_f1"]


def _mutate(**changes):
    artifact = {**VALID["artifacts"][0], **changes}
    return {"version": 1, "artifacts": [artifact]}


class TestRequiredParameters:
    def test_interval_metric_must_state_its_overlap_threshold(self, load):
        # The parameter PROJECT.md leaves open. There is deliberately no default.
        metric = {
            "type": "interval_overlap_f1",
            "group_column": "contig",
            "start_column": "start",
            "end_column": "end",
        }
        with pytest.raises(ConfigError, match="min_reciprocal_overlap.*no default"):
            load(_mutate(metrics=[metric]))

    def test_ccc_must_state_its_log_transform(self, load):
        with pytest.raises(ConfigError, match="log_transform.*no default"):
            load(_mutate(metrics=[{"type": "lins_ccc", "column": "x", "pseudocount": None}]))

    def test_ccc_must_state_its_pseudocount_even_when_unused(self, load):
        with pytest.raises(ConfigError, match="pseudocount.*no default"):
            load(_mutate(metrics=[{"type": "lins_ccc", "column": "x", "log_transform": False}]))

    @pytest.mark.parametrize(
        "field", ["path", "truth", "format", "soft", "row_key", "normalisers",
                  "column_aliases", "value_aliases", "metrics"]
    )
    def test_every_artifact_field_is_required(self, load, field):
        artifact = {k: v for k, v in VALID["artifacts"][0].items() if k != field}
        with pytest.raises(ConfigError, match="missing required field"):
            load({"version": 1, "artifacts": [artifact]})


class TestRejections:
    def test_misspelled_parameter_is_not_silently_ignored(self, load):
        metric = {"type": "lins_ccc", "column": "x", "log_transform": False,
                  "pseudocount": None, "log_base": 2}
        with pytest.raises(ConfigError, match="unexpected parameter"):
            load(_mutate(metrics=[metric]))

    def test_unknown_metric_type(self, load):
        with pytest.raises(ConfigError, match="unknown metric type"):
            load(_mutate(metrics=[{"type": "spearman"}]))

    def test_unknown_normalisation_rule(self, load):
        with pytest.raises(ConfigError, match="unknown normalisation rule"):
            load(_mutate(normalisers=["lowercase"]))

    def test_unknown_format(self, load):
        with pytest.raises(ConfigError, match="format must be one of"):
            load(_mutate(format="xlsx"))

    def test_empty_row_key(self, load):
        with pytest.raises(ConfigError, match="non-empty list"):
            load(_mutate(row_key=[]))

    def test_no_metrics(self, load):
        with pytest.raises(ConfigError, match="at least one metric"):
            load(_mutate(metrics=[]))

    def test_no_artifacts(self, load):
        with pytest.raises(ConfigError, match="declares no artifacts"):
            load({"version": 1, "artifacts": []})

    def test_duplicate_artifact_paths(self, load):
        document = {"version": 1, "artifacts": [VALID["artifacts"][0], VALID["artifacts"][0]]}
        with pytest.raises(ConfigError, match="duplicate artifact path"):
            load(document)

    def test_unsupported_version(self, load):
        with pytest.raises(ConfigError, match="unsupported grading.yaml version"):
            load({**VALID, "version": 2})

    def test_soft_must_be_a_boolean_not_an_int(self, load):
        # bool subclasses int, so this needs a deliberate check.
        with pytest.raises(ConfigError, match="must be bool"):
            load(_mutate(soft=1))

    def test_threshold_must_be_a_number_not_a_boolean(self, load):
        metric = {"type": "interval_overlap_f1", "group_column": "c", "start_column": "s",
                  "end_column": "e", "min_reciprocal_overlap": True}
        with pytest.raises(ConfigError, match="must be int or float"):
            load(_mutate(metrics=[metric]))

    def test_malformed_yaml(self, tmp_path):
        path = tmp_path / "grading.yaml"
        path.write_text("version: 1\nartifacts: [\n")
        with pytest.raises(ConfigError, match="not valid YAML"):
            load_grading_config(path)
