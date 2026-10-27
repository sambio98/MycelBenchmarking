"""Inspect AI task assembly.

Builds one Inspect `Task` from an NPBench-R bundle: the prompt as the sample,
a ReAct loop with the five tools, the caps from PROJECT.md §10, and the
artifact scorer.

    inspect eval npbench/harness/task.py -T task_dir=tasks/02_salinispora_bgc \
        --model anthropic/claude-opus-5

Caveat, recorded rather than hidden: the docker sandbox wiring below has not
been executed. The session this was written in had no docker daemon, so only
the scorer and the tools are covered by tests. Treat `sandbox_spec` as the
first thing to verify when a provisioned machine is available.
"""

from __future__ import annotations

from pathlib import Path

from inspect_ai import Task, task
from inspect_ai.agent import react
from inspect_ai.dataset import Sample
from inspect_ai.tool import bash_session, python, text_editor

from .adjudicator import Adjudicator, Policy
from .scorer import artifact_scorer
from .tools import request_web_access

# PROJECT.md §10: 32 vCPUs, 128 GB, no GPU, 24 hour cap, matching BixBench3 so
# results stay comparable.
TIME_LIMIT_SECONDS = 24 * 60 * 60
MESSAGE_LIMIT = 500

WORKSPACE_PROMPT = """\
You are working in a container. `data/` is read-only input. Write every \
deliverable to `outputs/`; nothing outside `outputs/` is collected or graded.

Direct network access is limited to package registries and OS mirrors. For \
anything else use the web access request tool. Requests that would reveal \
published results for this dataset are denied: the analysis has to be run, not \
looked up.

Work by writing and running scripts. The inputs are too large to read by eye.

When you are finished, call submit().
"""


def load_policy(task_dir: Path) -> Policy:
    """Read the task's network policy out of task.yaml."""
    import yaml

    meta = yaml.safe_load((task_dir / "task.yaml").read_text()) or {}
    network = meta.get("network") or {}
    return Policy(
        deny_terms=tuple(network.get("deny_terms") or ()),
        deny_hosts=tuple(network.get("deny_hosts") or ()),
        allow_hosts=tuple(network.get("allow_hosts") or ()),
    )


def sandbox_spec(task_dir: Path) -> tuple[str, str] | str:
    """Prefer a compose file if the task ships one, else its Dockerfile."""
    compose = task_dir / "environment" / "compose.yaml"
    if compose.is_file():
        return ("docker", str(compose))
    dockerfile = task_dir / "environment" / "Dockerfile"
    if dockerfile.is_file():
        return ("docker", str(dockerfile))
    return "local"


@task
def npbench_r(
    task_dir: str,
    message_limit: int = MESSAGE_LIMIT,
    time_limit: int = TIME_LIMIT_SECONDS,
) -> Task:
    """One NPBench-R task bundle as an Inspect task."""
    directory = Path(task_dir)
    adjudicator = Adjudicator(load_policy(directory))

    sample = Sample(
        id=directory.name,
        input=(directory / "prompt.md").read_text(),
        # Small committed inputs travel with the sample. The staged payload
        # under data/ is mounted by the sandbox, not copied per sample.
        files={
            name: str(directory / name)
            for name in ("strains.tsv", "data.manifest.json")
            if (directory / name).is_file()
        },
        metadata={"task_dir": str(directory)},
    )

    return Task(
        dataset=[sample],
        solver=react(
            prompt=WORKSPACE_PROMPT,
            tools=[
                bash_session(),
                python(),
                text_editor(),
                request_web_access(adjudicator),
            ],
        ),
        scorer=artifact_scorer(directory),
        sandbox=sandbox_spec(directory),
        message_limit=message_limit,
        time_limit=time_limit,
        name=f"npbench-r-{directory.name}",
    )
