"""Harness pieces that can be tested without a sandbox or a model.

CLAUDE.md build order, step 4, asks for two outcomes: a deliberately failing
run scores 0.00 with no crash, and the reference solution scores 1.00. Both are
exercised here with a stand-in reader in place of the sandbox.
"""

import asyncio
from pathlib import Path

import importlib.util

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


class TestTheDocumentedEntrypointLoads:
    """`inspect eval npbench/harness/task.py@npbench_r` must actually start.

    This is the invocation `make eval` uses and the one Inspect's docs show.
    Both of the bugs it had were invisible to every other test, because the
    harness was constructed in-process by its own package and never loaded the
    way Inspect loads it:

      * Inspect loads the task file **by path**, as a module with no package
        context, so `from .adjudicator import ...` raised ModuleNotFoundError
        before the task was built.
      * Inspect then changes the working directory to the task file's own
        directory, so a relative `-T task_dir=tasks/03_...` passed from the
        repository root resolved against `npbench/harness/` and the bundle was
        not found.

    Neither is reachable by importing the module normally, which is why they
    survived until the entrypoint was run end to end.
    """

    TASK_FILE = Path(__file__).resolve().parents[1] / "npbench" / "harness" / "task.py"

    def test_the_task_file_has_no_relative_imports(self):
        source = self.TASK_FILE.read_text()
        offenders = [
            line for line in source.splitlines()
            if line.startswith("from .") or line.startswith("from ..")
        ]
        assert offenders == [], (
            "a relative import in this file breaks loading it by path: " + "; ".join(offenders)
        )

    def test_it_loads_as_a_standalone_module_by_path(self):
        """Exactly how Inspect loads it: a path, as its own top-level module."""
        spec = importlib.util.spec_from_file_location(str(self.TASK_FILE), self.TASK_FILE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        assert hasattr(module, "npbench_r")

    def test_a_relative_task_dir_resolves_from_another_working_directory(self, tmp_path, monkeypatch):
        """Inspect cds into npbench/harness/ before calling the task function."""
        from npbench.harness.task import resolve_task_dir

        monkeypatch.chdir(self.TASK_FILE.parent)
        resolved = resolve_task_dir("tasks/03_mibig_cheminformatics")
        assert (resolved / "task.yaml").is_file()

    def test_the_task_itself_builds_from_another_working_directory(self, monkeypatch):
        """The resolution has to be wired into npbench_r, not merely available.

        Testing resolve_task_dir on its own passes even when npbench_r still
        calls Path(task_dir) directly, which was the actual bug. This builds
        the task the way Inspect does and is the test that catches it.
        """
        from npbench.harness.task import npbench_r

        monkeypatch.chdir(self.TASK_FILE.parent)
        built = npbench_r(task_dir="tasks/03_mibig_cheminformatics")
        assert built.name == "npbench-r-03_mibig_cheminformatics"
        assert len(built.dataset) == 1
        assert built.scorer is not None

    def test_an_absolute_task_dir_is_used_as_given(self):
        from npbench.harness.task import resolve_task_dir

        bundle = Path(__file__).resolve().parents[1] / "tasks" / "03_mibig_cheminformatics"
        assert resolve_task_dir(str(bundle)) == bundle

    def test_a_missing_bundle_says_where_it_looked(self):
        from npbench.harness.task import resolve_task_dir

        with pytest.raises(SystemExit) as caught:
            resolve_task_dir("tasks/99_does_not_exist")
        message = str(caught.value)
        assert "99_does_not_exist" in message and "repository root" in message
