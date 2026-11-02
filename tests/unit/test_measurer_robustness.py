"""No measurer may crash on an answer a real system could actually write.

A measurer that raises does not score a bad answer zero -- it turns the run into
a crash. The sweep records `measurement failed` rather than a score, and the
campaign's clear rates are then computed over whatever happened not to crash,
which is a silent selection on the systems being compared.

Every measurer in the suite already survives these shapes; the test is here so
the next one does too. The shapes are not invented: an empty object is what a
system writes when it decides it has nothing to report, wrong types are what a
system writes when it misreads the contract, and a bare list is what a system
writes when it answers with rows instead of a report.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.sweep.runner import _absolute_pythonpath

ROOT = pathlib.Path(__file__).resolve().parents[2] / "campaigns"
CAMPAIGNS = sorted(p for p in ROOT.iterdir()
                   if (p / "gold" / "gold.json").is_file())

pytestmark = pytest.mark.skipif(not CAMPAIGNS, reason="no campaign gold generated")

#: Each one a shape of answer, not a fuzz case: the point is that the measurer
#: reaches a verdict about it rather than that it survives arbitrary bytes.
SHAPES = {
    "empty_object": {},
    "wrong_types": {"present": "yes", "r3": 7, "modules": "none",
                    "accession": [], "selection_census": 3.5},
    "a_list": [1, 2, 3],
    "null": None,
}


@pytest.mark.parametrize("campaign", CAMPAIGNS, ids=lambda p: p.name)
def test_a_malformed_submission_scores_zero_rather_than_crashing(campaign, tmp_path):
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    measure_command = (task.get("harness") or {}).get("measure_command")
    assert measure_command, "campaign declares no harness.measure_command"
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    spec = json.loads((campaign / "grading.json").read_text())

    for shape, body in SHAPES.items():
        submission = tmp_path / shape / "submission"
        submission.mkdir(parents=True)
        for artefact in task["output_contract"]:
            if artefact.endswith(".json"):
                (submission / artefact).write_text(json.dumps(body))
            else:
                (submission / artefact).write_text(">x\nNOTASEQUENCE\n")
        measurement = tmp_path / shape / "measurement.json"
        proc = subprocess.run(
            [tok.format(campaign=str(campaign), submission=str(submission),
                        measurement=str(measurement)) for tok in measure_command],
            cwd=campaign, env=_absolute_pythonpath(os.environ),
            capture_output=True, text=True, timeout=1800,
        )
        assert proc.returncode == 0, f"{shape}: {(proc.stderr or '')[-400:]}"
        result = grade(json.loads(measurement.read_text()), gold, spec)
        assert result["depth"] == 0, f"{shape}: depth={result['depth']}"
        assert result["score"] == 0.0, f"{shape}: score={result['score']}"
