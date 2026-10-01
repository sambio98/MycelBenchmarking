"""Execute one system against one campaign, in an isolated sandbox.

Agent-agnostic: a system is a command. Real agents (a bash-only coding agent, a
general biomedical agent, a tool-equipped system) plug in through SystemSpec
without the harness knowing anything about them.

The load-bearing property here is **gold isolation**. If gold or oracle sources
leak into the sandbox, every number the sweep produces is silently invalid and
nothing downstream will notice. So the sandbox is built from an explicit
allowlist of public material and then *audited* before the system runs; a
violation raises rather than warns.
"""

from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Mapping, Sequence

from npbench_c.grading.grade import grade

# Only these campaign paths are ever visible to a system. "reference" carries
# the pinned tables the task requires an agent to use: without them this
# campaign is literally unsolvable, since R3's gold is a function of the codon
# weights. Stub systems hid that, because they were handed the oracle directory.
PUBLIC_CAMPAIGN_PATHS = ("inputs", "reference", "task.yaml")

# These must never appear in a sandbox, except as an explicitly declared
# warm-start input.
SECRET_CAMPAIGN_PATHS = ("gold", "oracle", "oracle_submission", "grading.json",
                         "internal_sweep.json")

COMPLETED, TIMEOUT, CRASHED, NO_SUBMISSION = "completed", "timeout", "crashed", "no_submission"
INFRA_STATUSES = (TIMEOUT, CRASHED)


class GoldLeak(RuntimeError):
    """Gold or oracle material reached a system sandbox. Invalidates the sweep."""


@dataclass(frozen=True)
class SystemSpec:
    name: str
    command: Sequence[str]
    mode: str = "tooled"              # "tooled" | "no_tool" (leakage ablation)
    timeout_s: int = 1800
    env: Mapping[str, str] = field(default_factory=dict)
    is_reference: bool = False        # the system the difficulty gate is judged on
    is_stub: bool = False             # harness self-test fixture, not an agent


@dataclass
class RunResult:
    system: str
    run_index: int
    rung_mode: str                    # "cold" or "warm:<rung_id>"
    status: str
    grade: dict | None
    stderr_tail: str = ""

    @property
    def is_infra_failure(self) -> bool:
        return self.status in INFRA_STATUSES

    @property
    def depth(self) -> int:
        return self.grade["depth"] if self.grade else 0

    @property
    def score(self) -> float:
        return self.grade["score"] if self.grade else 0.0

    def cleared(self, rung_id: str) -> bool:
        if not self.grade:
            return False
        return any(r["rung_id"] == rung_id and r["passed"] for r in self.grade["rungs"])


# Never admissible in a sandbox, not even as a warm-start input: these carry
# answers for rungs the run has not reached.
NEVER_IN_SANDBOX = ("gold.json", "grading.json", "measurement.json")


def _audit_sandbox(sandbox: pathlib.Path) -> None:
    """Fail loudly if anything secret is visible.

    Warm-start material arrives as a pre-rendered bundle whose contents the
    campaign author vetted, so the audit needs no per-file allowlist: the
    directory names and the never-admissible filenames are enough.
    """
    for path in sorted(sandbox.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(sandbox)
        if rel.parts[0] in SECRET_CAMPAIGN_PATHS or rel.name in NEVER_IN_SANDBOX:
            raise GoldLeak(f"{rel} is visible in {sandbox}")


def build_sandbox(
    campaign: pathlib.Path,
    sandbox: pathlib.Path,
    warm_start_bundle: str | None = None,
) -> pathlib.Path:
    """Create a sandbox holding public campaign material plus a warm bundle."""
    if sandbox.exists():
        shutil.rmtree(sandbox)
    sandbox.mkdir(parents=True)
    for name in PUBLIC_CAMPAIGN_PATHS:
        src = campaign / name
        if src.is_dir():
            shutil.copytree(src, sandbox / name)
        elif src.is_file():
            shutil.copyfile(src, sandbox / name)

    # Warm start: the rendered bundle's contents, flattened into the sandbox
    # root so a system never sees the gold directory structure.
    if warm_start_bundle:
        bundle = campaign / warm_start_bundle
        if not bundle.is_dir():
            raise FileNotFoundError(f"warm-start bundle missing: {bundle}")
        for src in sorted(bundle.iterdir()):
            if src.is_file():
                shutil.copyfile(src, sandbox / src.name)

    (sandbox / "submission").mkdir()
    _audit_sandbox(sandbox)
    return sandbox


def _format(tokens: Sequence[str], **subs: str) -> list[str]:
    return [t.format(**subs) for t in tokens]


def run_once(
    campaign: pathlib.Path,
    system: SystemSpec,
    run_index: int,
    sandbox_root: pathlib.Path,
    measure_command: Sequence[str],
    rung_mode: str = "cold",
    warm_start_bundle: str | None = None,
) -> RunResult:
    """One system, one repeat. Builds the sandbox, runs, measures, grades."""
    import os

    sandbox = build_sandbox(
        campaign, sandbox_root / system.name / rung_mode.replace(":", "_") / str(run_index),
        warm_start_bundle,
    )
    submission = sandbox / "submission"

    env = {**os.environ, **system.env, "NPBENCH_MODE": system.mode}
    try:
        proc = subprocess.run(
            _format(system.command, submission=str(submission), sandbox=str(sandbox)),
            cwd=sandbox, env=env, capture_output=True, timeout=system.timeout_s, text=True,
        )
    except subprocess.TimeoutExpired:
        return RunResult(system.name, run_index, rung_mode, TIMEOUT, None,
                         "timed out")
    if proc.returncode != 0:
        return RunResult(system.name, run_index, rung_mode, CRASHED, None,
                         (proc.stderr or "")[-400:])
    if not any(submission.iterdir()):
        return RunResult(system.name, run_index, rung_mode, NO_SUBMISSION, None,
                         "submission directory empty")

    measurement = sandbox / "measurement.json"
    try:
        subprocess.run(
            _format(measure_command, submission=str(submission),
                    measurement=str(measurement)),
            cwd=campaign, check=True, capture_output=True, timeout=1800, text=True,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        return RunResult(system.name, run_index, rung_mode, CRASHED, None,
                         f"measurement failed: {exc}")

    g = grade(
        json.loads(measurement.read_text()),
        json.loads((campaign / "gold" / "gold.json").read_text()),
        json.loads((campaign / "grading.json").read_text()),
    )
    return RunResult(system.name, run_index, rung_mode, COMPLETED, g)
