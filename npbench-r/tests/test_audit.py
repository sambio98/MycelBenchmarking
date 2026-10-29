"""Audit checks: each one catches the defect it was written for, and stays quiet otherwise.

The leak and column checks exist because both defects were found by hand during
the build -- a prompt that stated a graded row count, and artifacts whose
columns were required of the agent but read by no metric. These tests pin the
behaviour that found them, including the three false-positive classes that made
the first version of the column check useless.
"""

import shutil
from pathlib import Path

import pytest
import yaml

from npbench.audit import (
    check_columns,
    check_leakage,
    column_list,
    dossier,
    render,
    required_outputs,
)
from npbench.grading import load_grading_config

FIXTURE = Path(__file__).parent / "fixtures" / "minimal_task"


@pytest.fixture
def bundle(tmp_path):
    target = tmp_path / "minimal_task"
    shutil.copytree(FIXTURE, target)
    return target


def _config(bundle: Path):
    return load_grading_config(bundle / "grading.yaml")


def test_the_fixture_task_has_nothing_to_report(bundle):
    assert dossier(bundle).findings == []


# --------------------------------------------------------------------- leakage


def _grow_truth(bundle: Path, rows: int) -> int:
    """Give the fixture artifact enough rows to clear the row-count floor."""
    path = bundle / "truth" / "bgc_inventory.csv"
    lines = path.read_text().splitlines()
    header, body = lines[0], lines[1:]
    template = body[0].split(",")
    while len(body) < rows:
        template[0] = f"BGC{len(body):04d}"
        body.append(",".join(template))
    path.write_text("\n".join([header, *body]) + "\n")
    return len(body)


def test_a_stated_row_count_is_reported(bundle):
    rows = _grow_truth(bundle, 12)
    prompt = bundle / "prompt.md"
    prompt.write_text(
        prompt.read_text() + f"\n\nYou should find {rows} clusters in total.\n"
    )
    findings = check_leakage(bundle, _config(bundle))
    assert [f.level for f in findings] == ["leak"]
    assert "row count of bgc_inventory.csv" in findings[0].message


def test_a_method_parameter_is_not_reported_merely_for_being_a_number(bundle):
    prompt = bundle / "prompt.md"
    prompt.write_text(prompt.read_text() + "\n\nUse a 30000 bp window.\n")
    assert check_leakage(bundle, _config(bundle)) == []


def test_small_numbers_are_below_the_cell_floor(bundle):
    """A prompt is full of 1-based, 2-decimal, 6-peak instructions.

    Every one of them matches a cell in some truth file, so without a floor the
    check would report dozens of findings and be ignored.
    """
    prompt = bundle / "prompt.md"
    prompt.write_text(prompt.read_text() + "\n\nCoordinates are 1-based; keep 2 decimals.\n")
    assert check_leakage(bundle, _config(bundle)) == []


def test_an_acknowledged_finding_leaves_the_open_list(bundle):
    rows = _grow_truth(bundle, 12)
    prompt = bundle / "prompt.md"
    prompt.write_text(prompt.read_text() + f"\n\nYou should find {rows} clusters.\n")
    assert len(dossier(bundle).findings) == 1

    meta = yaml.safe_load((bundle / "task.yaml").read_text())
    meta["audit_acknowledged"] = [
        {"contains": "row count of bgc_inventory.csv", "why": "countable from data/"}
    ]
    (bundle / "task.yaml").write_text(yaml.safe_dump(meta, sort_keys=False))

    report = dossier(bundle)
    assert report.findings == []
    assert len(report.settled) == 1
    assert report.settled[0][1] == "countable from data/"
    assert "already judged" in render(report)


def test_an_acknowledgement_matches_text_not_a_line_number(bundle):
    """Acknowledgements are matched on the message so an edit above them holds."""
    rows = _grow_truth(bundle, 12)
    prompt = bundle / "prompt.md"
    prompt.write_text(
        "<!-- a line inserted above everything -->\n"
        + prompt.read_text()
        + f"\n\nYou should find {rows} clusters.\n"
    )
    meta = yaml.safe_load((bundle / "task.yaml").read_text())
    meta["audit_acknowledged"] = [
        {"contains": "row count of bgc_inventory.csv", "why": "structural"}
    ]
    (bundle / "task.yaml").write_text(yaml.safe_dump(meta, sort_keys=False))
    assert dossier(bundle).findings == []


# --------------------------------------------------------------------- columns


def test_a_graded_column_the_prompt_never_names_is_unfair(bundle):
    prompt = bundle / "prompt.md"
    prompt.write_text(prompt.read_text().replace("`length_bp`", "`length`"))
    findings = check_columns(bundle, _config(bundle))
    assert [f.level for f in findings] == ["unfair"]
    assert "length_bp" in findings[0].message


def test_a_prompt_without_per_artifact_blocks_is_still_checked(bundle):
    """The fixture states its columns in a paragraph, not under a heading.

    PROJECT.md §6 asks only that every column be stated in the required-outputs
    section, so a prompt that does not use one heading per artifact must still
    be checked rather than skipped.
    """
    text = (bundle / "prompt.md").read_text()
    assert "### `bgc_inventory.csv`" not in text
    assert required_outputs(text).strip()
    assert check_columns(bundle, _config(bundle)) == []


def test_values_a_column_may_take_are_not_mistaken_for_columns():
    """`start`, `end` and `both` are the values of a `boundary` column."""
    block = (
        "One row per cluster on a contig boundary.\n\n"
        "`cluster_id`, `accession`, `boundary`\n\n"
        "`boundary` is `start`, `end`, or `both`.\n"
    )
    assert column_list(block) == {"cluster_id", "accession", "boundary"}


def test_input_field_names_in_prose_are_not_mistaken_for_columns():
    block = (
        "One row per substrate whose `proteinogenic` flag is false.\n\n"
        "`substrate_name`, `formula`\n\n"
        "`formula` is the Hill-order formula.\n"
    )
    assert column_list(block) == {"substrate_name", "formula"}


def test_a_schema_line_may_wrap_across_two_lines():
    block = (
        "One row per compound.\n\n"
        "`npaid`, `mol_formula`, `exact_mass`,\n"
        "`origin_type`, `genus`\n"
    )
    assert column_list(block) == {
        "npaid", "mol_formula", "exact_mass", "origin_type", "genus",
    }


# --------------------------------------------------------------------- dossier


def test_the_dossier_counts_hard_artifacts_and_metrics(bundle):
    report = dossier(bundle)
    config = _config(bundle)
    assert report.n_hard == sum(1 for a in config.artifacts if not a.soft)
    assert report.n_metrics == sum(len(a.metrics) for a in config.artifacts)
    assert report.truth_rows["bgc_inventory.csv"] > 0


def test_the_dossier_renders_without_a_task_yaml(bundle):
    (bundle / "task.yaml").unlink()
    text = render(dossier(bundle))
    assert "minimal_task" in text
