"""Harness pieces that can be tested without a sandbox or a model.

CLAUDE.md build order, step 4, asks for two outcomes: a deliberately failing
run scores 0.00 with no crash, and the reference solution scores 1.00. Both are
exercised here with a stand-in reader in place of the sandbox.
"""

import asyncio
from pathlib import Path

import pytest

from npbench.grading import load_grading_config, score_task
from npbench.harness.adjudicator import Adjudicator, Policy
from npbench.harness.scorer import METHODS, collect, summarise
from npbench.harness.task import load_policy, sandbox_spec
from npbench.harness.tools import decide_and_describe

FIXTURE = Path(__file__).parent / "fixtures" / "minimal_task"


def _reader(files: dict[str, str]):
    async def read(path: str) -> str:
        if path not in files:
            raise FileNotFoundError(path)
        return files[path]

    return read


def _run(coro):
    return asyncio.run(coro)


@pytest.fixture
def config():
    return load_grading_config(FIXTURE / "grading.yaml")


class TestCollectAndScore:
    def test_a_perfect_run_scores_one(self, config, tmp_path):
        truth = (FIXTURE / "truth" / "bgc_inventory.csv").read_text()
        read = _reader({"outputs/bgc_inventory.csv": truth, f"outputs/{METHODS}": "ran GECCO"})
        collected = _run(collect(config, tmp_path, read))
        result = score_task(config, FIXTURE / "truth", tmp_path)
        assert result.task_score == 1.0
        assert sorted(collected) == [METHODS, "bgc_inventory.csv"]
        assert summarise(result, collected)["methods_written"] is True

    def test_an_agent_that_wrote_nothing_scores_zero_without_crashing(self, config, tmp_path):
        collected = _run(collect(config, tmp_path, _reader({})))
        result = score_task(config, FIXTURE / "truth", tmp_path)
        assert collected == []
        assert result.task_score == 0.0
        summary = summarise(result, collected)
        assert summary["artifacts"][0]["reason"] == "file_missing"
        assert summary["methods_written"] is False

    def test_a_reader_that_raises_is_not_fatal(self, config, tmp_path):
        async def read(path: str) -> str:
            raise RuntimeError("sandbox went away")

        assert _run(collect(config, tmp_path, read)) == []
        assert score_task(config, FIXTURE / "truth", tmp_path).task_score == 0.0

    def test_a_partial_run_is_scored_not_rejected(self, config, tmp_path):
        rows = (FIXTURE / "truth" / "bgc_inventory.csv").read_text().splitlines()
        half = "\n".join(rows[:3]) + "\n"
        read = _reader({"outputs/bgc_inventory.csv": half})
        _run(collect(config, tmp_path, read))
        result = score_task(config, FIXTURE / "truth", tmp_path)
        assert 0.0 < result.artifacts[0].score < 1.0

    def test_summarise_records_the_format_axis(self, config, tmp_path):
        truth = (FIXTURE / "truth" / "bgc_inventory.csv").read_text()
        renamed = truth.replace("bgc_id", "BGC_ID", 1)
        _run(collect(config, tmp_path, _reader({"outputs/bgc_inventory.csv": renamed})))
        result = score_task(config, FIXTURE / "truth", tmp_path)
        summary = summarise(result, ["bgc_inventory.csv"])
        assert summary["task_score"] == 1.0
        assert summary["format_compliance"] == 0.0
        assert summary["artifacts"][0]["renamed_columns"] == {"bgc_id": "BGC_ID"}


class TestWebAccessTool:
    @pytest.fixture
    def adj(self):
        return Adjudicator(Policy(deny_terms=("10.0000/fixture",)))

    def test_allow_message(self, adj):
        allowed, message = decide_and_describe(adj, "https://pypi.org/simple/gecco/")
        assert allowed and message.startswith("ALLOWED:")

    def test_denial_explains_the_boundary(self, adj):
        allowed, message = decide_and_describe(adj, "https://doi.org/10.1111/x")
        assert not allowed
        assert "has to be run, not" in message

    def test_a_task_deny_term_is_refused(self, adj):
        allowed, _ = decide_and_describe(adj, "https://github.com/q?=10.0000/fixture")
        assert not allowed


class TestTaskWiring:
    def test_policy_comes_from_task_yaml(self):
        policy = load_policy(FIXTURE)
        assert "10.0000/fixture" in policy.deny_terms

    def test_sandbox_prefers_the_tasks_own_dockerfile(self):
        assert sandbox_spec(FIXTURE) == ("docker", str(FIXTURE / "environment" / "Dockerfile"))

    def test_sandbox_falls_back_to_local_without_one(self, tmp_path):
        assert sandbox_spec(tmp_path) == "local"
