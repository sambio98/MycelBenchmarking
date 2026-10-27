"""CLI exit codes. A broken bundle, a missing truth and a failed submission must
be distinguishable, because `make selfgrade` is a gate in the build order."""

import shutil
from pathlib import Path

import pytest

from npbench.build.validate import solution_steps
from npbench.cli import main

FIXTURE = Path(__file__).parent / "fixtures" / "minimal_task"


@pytest.fixture
def bundle(tmp_path):
    target = tmp_path / "minimal_task"
    shutil.copytree(FIXTURE, target)
    return target


def test_validate_passes_the_fixture(bundle, capsys):
    assert main(["validate", str(bundle)]) == 0
    assert "bundle contract OK" in capsys.readouterr().out


def test_validate_fails_a_broken_bundle(bundle, capsys):
    (bundle / "task.yaml").unlink()
    assert main(["validate", str(bundle)]) == 1
    assert "error(s)" in capsys.readouterr().out


def test_grade_a_perfect_submission(bundle, capsys):
    assert main(["grade", str(bundle), str(bundle / "truth"), "--expect-perfect"]) == 0
    assert "task score      1.00" in capsys.readouterr().out


def test_grade_an_empty_submission_scores_zero(bundle, tmp_path, capsys):
    empty = tmp_path / "nothing"
    empty.mkdir()
    assert main(["grade", str(bundle), str(empty)]) == 0
    assert "0.00" in capsys.readouterr().out


def test_expect_perfect_fails_an_imperfect_submission(bundle, tmp_path, capsys):
    out = tmp_path / "outputs"
    out.mkdir()
    rows = (bundle / "truth" / "bgc_inventory.csv").read_text().splitlines()
    (out / "bgc_inventory.csv").write_text("\n".join(rows[:3]) + "\n")
    assert main(["grade", str(bundle), str(out), "--expect-perfect"]) == 1
    assert "Expected 1.00" in capsys.readouterr().err


def test_missing_truth_is_reported_not_crashed(bundle, capsys):
    shutil.rmtree(bundle / "truth")
    assert main(["grade", str(bundle), str(bundle)]) == 2
    assert "Run `make truth" in capsys.readouterr().err


def test_a_broken_truth_file_is_reported_as_a_bundle_bug(bundle, capsys):
    # Constant numeric column: undefined for CCC, and our fault, not the agent's.
    path = bundle / "truth" / "bgc_inventory.csv"
    rows = path.read_text().splitlines()
    path.write_text("\n".join([rows[0]] + [r.rsplit(",", 1)[0] + ",40000" for r in rows[1:]]) + "\n")
    assert main(["grade", str(bundle), str(bundle / "truth")]) == 3
    assert "bundle is broken, not the submission" in capsys.readouterr().err


def test_solution_steps_are_read_from_task_yaml(bundle):
    assert solution_steps(bundle) == ["python {task}/solution/make_truth.py --out {out}"]


def test_solution_steps_rejects_a_bundle_without_them(bundle):
    path = bundle / "task.yaml"
    path.write_text(path.read_text().replace("solution:", "not_solution:"))
    with pytest.raises(ValueError, match="solution.steps"):
        solution_steps(bundle)
