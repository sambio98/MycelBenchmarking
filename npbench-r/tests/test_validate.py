"""Bundle validator: passes the fixture task, fails each deliberate break.

Each case mutates a copy of tests/fixtures/minimal_task and asserts the
validator says something a task author could act on, not just that it failed.
"""

import json
import shutil
from pathlib import Path

import pytest
import yaml

from npbench.build.validate import Finding, report, validate_task

FIXTURE = Path(__file__).parent / "fixtures" / "minimal_task"


@pytest.fixture
def bundle(tmp_path):
    target = tmp_path / "minimal_task"
    shutil.copytree(FIXTURE, target)
    return target


def _errors(bundle: Path) -> list[Finding]:
    return [f for f in validate_task(bundle) if f.level == "error"]


def _edit_yaml(path: Path, mutate):
    document = yaml.safe_load(path.read_text())
    mutate(document)
    path.write_text(yaml.safe_dump(document, sort_keys=False))


def test_the_fixture_task_passes_clean(bundle):
    assert validate_task(bundle) == []


def test_the_fixture_truth_grades_itself_at_one(bundle):
    # If truth cannot score 1.00 against itself, the grader or the truth is wrong
    # and nothing downstream is trustworthy (CLAUDE.md build order, step 3).
    from npbench.grading import load_grading_config, score_task

    config = load_grading_config(bundle / "grading.yaml")
    result = score_task(config, bundle / "truth", bundle / "truth")
    assert result.task_score == 1.0
    assert result.format_compliance == 1.0


# --------------------------------------------------- the five deliberate breaks


def test_break_1_truth_file_declared_but_absent(bundle):
    (bundle / "truth" / "bgc_inventory.csv").unlink()
    assert "not present in truth/" in _errors(bundle)[0].message


def test_break_2_manifest_entry_without_a_checksum(bundle):
    path = bundle / "data.manifest.json"
    entries = json.loads(path.read_text())
    del entries[0]["sha256"]
    path.write_text(json.dumps(entries))
    assert any("sha256 must be 64 lowercase hex" in e.message for e in _errors(bundle))


def test_break_3_artifact_with_no_metrics(bundle):
    _edit_yaml(bundle / "grading.yaml", lambda d: d["artifacts"][0].update(metrics=[]))
    assert any("at least one metric" in e.message for e in _errors(bundle))


def test_break_4_metric_names_a_column_absent_from_truth(bundle):
    def retarget_ccc(document):
        for metric in document["artifacts"][0]["metrics"]:
            if metric["type"] == "lins_ccc":
                metric["column"] = "gc_content"

    _edit_yaml(bundle / "grading.yaml", retarget_ccc)
    assert any("missing_column(s): gc_content" in e.message for e in _errors(bundle))


def test_break_5_aligned_metric_without_identifier_f1(bundle):
    # The rule that keeps dropped rows from going unpunished: CCC and macro F1
    # only ever see the intersection.
    def drop_identifier_f1(document):
        metrics = document["artifacts"][0]["metrics"]
        document["artifacts"][0]["metrics"] = [
            m for m in metrics if m["type"] != "identifier_f1"
        ]

    _edit_yaml(bundle / "grading.yaml", drop_identifier_f1)
    message = next(e.message for e in _errors(bundle) if "identifier_f1" in e.message)
    assert "lins_ccc, macro_f1" in message
    assert "dropped half the rows" in message


# ------------------------------------------------------- the rest of the contract


class TestLayout:
    @pytest.mark.parametrize("name", ["task.yaml", "prompt.md", "grading.yaml",
                                      "data.manifest.json"])
    def test_missing_required_file(self, bundle, name):
        (bundle / name).unlink()
        assert any("required file is missing" in e.message for e in _errors(bundle))

    def test_solution_directory_without_a_script(self, bundle):
        (bundle / "solution" / "make_truth.py").unlink()
        assert any("not usable as ground truth" in e.message for e in _errors(bundle))

    def test_environment_without_a_dockerfile_or_lock(self, bundle):
        (bundle / "environment" / "Dockerfile").unlink()
        assert any("pinned environment" in e.message for e in _errors(bundle))


class TestTaskYaml:
    def test_unpinned_tool_version(self, bundle):
        _edit_yaml(bundle / "task.yaml",
                   lambda d: d.update(tool_versions={"antismash": "latest"}))
        assert any("not an exact version" in e.message for e in _errors(bundle))

    def test_tool_version_without_a_number(self, bundle):
        _edit_yaml(bundle / "task.yaml",
                   lambda d: d.update(tool_versions={"antismash": "stable"}))
        assert any("not an exact version" in e.message for e in _errors(bundle))

    def test_task_id_must_match_the_directory(self, bundle):
        _edit_yaml(bundle / "task.yaml", lambda d: d.update(task_id="something_else"))
        assert any("does not match directory name" in e.message for e in _errors(bundle))

    def test_runtime_over_the_cap(self, bundle):
        _edit_yaml(bundle / "task.yaml",
                   lambda d: d["runtime"].update(wall_clock_hours=48))
        assert any("over the 24 hour cap" in e.message for e in _errors(bundle))

    def test_missing_doi(self, bundle):
        _edit_yaml(bundle / "task.yaml", lambda d: d["source"].pop("doi"))
        assert any("source is missing 'doi'" in e.message for e in _errors(bundle))


class TestNetworkPolicy:
    def test_missing_network_block(self, bundle):
        _edit_yaml(bundle / "task.yaml", lambda d: d.pop("network"))
        assert any("declaring deny_terms" in e.message for e in _errors(bundle))

    def test_empty_deny_terms(self, bundle):
        _edit_yaml(bundle / "task.yaml", lambda d: d.update(network={"deny_terms": []}))
        assert any("must be a non-empty list" in e.message for e in _errors(bundle))

    def test_deny_terms_that_miss_the_source_doi(self, bundle):
        # The leak check CLAUDE.md asks for: the paper must not be one search away.
        _edit_yaml(
            bundle / "task.yaml",
            lambda d: d.update(network={"deny_terms": ["something-unrelated"]}),
        )
        assert any("does not cover the source DOI" in e.message for e in _errors(bundle))

    def test_missing_solution_steps(self, bundle):
        _edit_yaml(bundle / "task.yaml", lambda d: d.pop("solution"))
        assert any("solution.steps" in e.message for e in _errors(bundle))

    def test_missing_truth_provenance(self, bundle):
        _edit_yaml(bundle / "task.yaml", lambda d: d.pop("truth_provenance"))
        assert any("truth_provenance" in e.message for e in _errors(bundle))


class TestPrompt:
    def test_missing_a_required_section(self, bundle):
        path = bundle / "prompt.md"
        path.write_text(path.read_text().replace("Method guidance", "Notes"))
        assert any("no method guidance section" in e.message for e in _errors(bundle))

    def test_graded_artifact_the_prompt_never_asks_for(self, bundle):
        path = bundle / "prompt.md"
        path.write_text(path.read_text().replace("`bgc_inventory.csv`", "`inventory.csv`"))
        assert any("not a fair test" in e.message for e in _errors(bundle))


class TestManifest:
    def test_entry_without_a_url_or_accession(self, bundle):
        path = bundle / "data.manifest.json"
        entries = json.loads(path.read_text())
        del entries[0]["url"], entries[0]["accession"]
        path.write_text(json.dumps(entries))
        assert any("needs a 'url' or an 'accession'" in e.message for e in _errors(bundle))

    def test_reference_database_without_a_release(self, bundle):
        path = bundle / "data.manifest.json"
        entries = json.loads(path.read_text())
        del entries[1]["release"]
        path.write_text(json.dumps(entries))
        assert any("must record its release" in e.message for e in _errors(bundle))

    def test_path_escaping_the_bundle(self, bundle):
        path = bundle / "data.manifest.json"
        entries = json.loads(path.read_text())
        entries[0]["path"] = "../../etc/passwd"
        path.write_text(json.dumps(entries))
        assert any("stay inside the bundle" in e.message for e in _errors(bundle))

    def test_zero_byte_entry(self, bundle):
        path = bundle / "data.manifest.json"
        entries = json.loads(path.read_text())
        entries[0]["bytes"] = 0
        path.write_text(json.dumps(entries))
        assert any("bytes must be a positive integer" in e.message for e in _errors(bundle))


class TestTruthContents:
    def test_duplicate_row_keys_in_truth(self, bundle):
        path = bundle / "truth" / "bgc_inventory.csv"
        path.write_text(path.read_text() + "BGC0001,nrps,contig_1,1000,46000,45001\n")
        assert any("duplicate row keys" in e.message for e in _errors(bundle))

    def test_constant_column_graded_with_ccc(self, bundle):
        path = bundle / "truth" / "bgc_inventory.csv"
        rows = path.read_text().splitlines()
        path.write_text(
            "\n".join([rows[0]] + [r.rsplit(",", 1)[0] + ",40000" for r in rows[1:]]) + "\n"
        )
        assert any("is constant, so CCC is undefined" in e.message for e in _errors(bundle))

    def test_inverted_interval_in_truth(self, bundle):
        path = bundle / "truth" / "bgc_inventory.csv"
        path.write_text(
            path.read_text().replace("contig_1,1000,46000", "contig_1,46000,1000")
        )
        assert any("start > end" in e.message for e in _errors(bundle))

    def test_non_numeric_value_in_a_ccc_column(self, bundle):
        path = bundle / "truth" / "bgc_inventory.csv"
        path.write_text(path.read_text().replace(",45001", ",~45001"))
        assert any("non-numeric values" in e.message for e in _errors(bundle))


class TestGradingConsistency:
    def test_two_ccc_metrics_on_different_columns_are_allowed(self, bundle):
        # An artifact can legitimately have more than one numeric column worth
        # grading; the product combiner handles it.
        def add_second_ccc(document):
            document["artifacts"][0]["metrics"].append(
                {"type": "lins_ccc", "column": "start",
                 "log_transform": False, "pseudocount": None}
            )

        _edit_yaml(bundle / "grading.yaml", add_second_ccc)
        assert _errors(bundle) == []

    def test_two_ccc_metrics_on_the_same_column_are_rejected(self, bundle):
        # A copy-paste slip that would silently square that column's score.
        def duplicate_ccc(document):
            metrics = document["artifacts"][0]["metrics"]
            existing = next(m for m in metrics if m["type"] == "lins_ccc")
            metrics.append(dict(existing))

        _edit_yaml(bundle / "grading.yaml", duplicate_ccc)
        assert any("lins_ccc twice on 'length_bp'" in e.message for e in _errors(bundle))

    def test_two_identifier_f1_metrics_are_rejected(self, bundle):
        _edit_yaml(
            bundle / "grading.yaml",
            lambda d: d["artifacts"][0]["metrics"].append({"type": "identifier_f1"}),
        )
        assert any("identifier_f1 twice" in e.message for e in _errors(bundle))


    def test_artifact_path_must_not_repeat_the_outputs_prefix(self, bundle):
        _edit_yaml(
            bundle / "grading.yaml",
            lambda d: d["artifacts"][0].update(path="outputs/bgc_inventory.csv"),
        )
        assert any("must not start with 'outputs/'" in e.message for e in _errors(bundle))

    def test_two_artifacts_sharing_one_truth_file(self, bundle):
        def duplicate(document):
            second = dict(document["artifacts"][0], path="copy.csv")
            document["artifacts"].append(second)

        _edit_yaml(bundle / "grading.yaml", duplicate)
        assert any("share truth file" in e.message for e in _errors(bundle))


def test_missing_truth_directory_is_a_note_not_an_error(bundle):
    # `make validate` runs before `make truth`, so an unbuilt truth/ must not
    # fail the contract check.
    shutil.rmtree(bundle / "truth")
    findings = validate_task(bundle)
    assert [f.level for f in findings] == ["note"]
    assert "run `make truth`" in findings[0].message


def test_report_counts_errors_and_notes(bundle):
    (bundle / "truth" / "bgc_inventory.csv").unlink()
    text = report(validate_task(bundle))
    assert "1 error(s), 0 note(s)" in text
