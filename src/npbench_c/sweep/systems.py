"""Real systems for the sweep, declared in a config file rather than in code.

The harness has always been agent-agnostic -- a system is a command -- and the
readiness gate has always kept its five agent-run checks PENDING until one runs.
What was missing was never the runner. It was two things: the instructions a real
agent needs (`npbench_c.sweep.prompt`) and a declared, validated set of systems
to run. This is the second.

A config is deliberately small, because everything else is derived:

    {
      "systems": [
        {"name": "tool-equipped", "reference": true,
         "command": ["bash", "-lc", "my-agent --prompt {prompt} --out {submission}"]},
        {"name": "bash-only",
         "command": ["bash", "-lc", "my-agent --no-tools --prompt {prompt} --out {submission}"]},
        {"name": "noncompute", "mode": "no_tool",
         "command": ["bash", "-lc", "my-agent --prompt {prompt} --out {submission}"]}
      ]
    }

What this module refuses, and why each refusal exists:

**A command that never mentions `{submission}`.** The run would be scored on an
empty directory and look like a capability result.

**No `no_tool` ablation, or more than one.** The leakage check compares that
system's score against the chance floor; without it the sweep cannot say whether
a campaign is answerable without doing the work, and the gate's own
`leakage_no_tool_ablation` fails rather than passing quietly.

**No reference system, or more than one.** The difficulty gate is judged on the
reference system's per-rung profile. Two of them is an ambiguous gate.

**A `no_tool` system that can still see the tool prefix.** The ablation's whole
claim is that it did not compute. The loader strips the prefix from its PATH and
refuses a config that puts it back.

**A stub.** `is_stub` marks a harness fixture, and a fixture in a real sweep
would be laundered into evidence -- which is the one thing the gate exists to
prevent.

The per-system TIMEOUT is not configured: it comes from the campaign's own
declared `budgets.wall_clock_minutes`, so no system is silently given more or
less wall clock than the campaign budgeted for. A config may lower it, never
raise it.
"""

from __future__ import annotations

import json
import os
import pathlib
from typing import Sequence

import yaml

from npbench_c.sweep.runner import SystemSpec
from npbench_c.tools import resources

#: Added to the campaign's declared budget before it becomes a hard timeout, so
#: a system that uses its whole budget is killed by the grader's clock rather
#: than by ours.
TIMEOUT_MARGIN = 1.25

#: Modes a config may declare.
MODES = ("tooled", "no_tool")


class SystemConfigError(RuntimeError):
    """The declared systems cannot produce a sweep worth reading."""


def _budget_seconds(campaign: pathlib.Path) -> int:
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    minutes = ((task.get("budgets") or {}).get("wall_clock_minutes"))
    if not minutes:
        raise SystemConfigError(
            f"{campaign.name} declares no budgets.wall_clock_minutes, so there "
            "is no budget to derive a timeout from")
    return int(float(minutes) * 60 * TIMEOUT_MARGIN)


def _without_tool_prefix(path: str) -> str:
    prefix = str(resources.tool_prefix())
    kept = [part for part in path.split(os.pathsep)
            if part and not part.startswith(prefix)]
    return os.pathsep.join(kept)


def load(config: pathlib.Path) -> list[dict]:
    """The declared systems, validated as a set rather than one at a time."""
    document = json.loads(pathlib.Path(config).read_text())
    declared = document.get("systems")
    if not isinstance(declared, list) or not declared:
        raise SystemConfigError(
            f"{config} declares no systems; a sweep needs at least a tooled "
            "system and a no-tool ablation")

    names = [entry.get("name") for entry in declared]
    if len(set(names)) != len(names) or not all(names):
        raise SystemConfigError(f"system names must be present and unique: {names}")
    for entry in declared:
        if entry.get("is_stub"):
            raise SystemConfigError(
                f"{entry['name']} declares is_stub: a fixture in a real sweep "
                "would be laundered into evidence")
        mode = entry.get("mode", "tooled")
        if mode not in MODES:
            raise SystemConfigError(
                f"{entry['name']} declares mode {mode!r}; one of {MODES}")
        command = entry.get("command")
        if not isinstance(command, list) or not command:
            raise SystemConfigError(f"{entry['name']} declares no command")
        if not any("{submission}" in str(part) for part in command):
            raise SystemConfigError(
                f"{entry['name']}'s command never mentions {{submission}}, so "
                "the run would be scored on an empty directory")
        env = entry.get("env") or {}
        if mode == "no_tool" and "PATH" in env:
            if str(resources.tool_prefix()) in env["PATH"]:
                raise SystemConfigError(
                    f"{entry['name']} is the no-tool ablation and its PATH "
                    "carries the tool prefix; the ablation's claim is that it "
                    "did not compute")

    ablations = [e["name"] for e in declared if e.get("mode") == "no_tool"]
    if len(ablations) != 1:
        raise SystemConfigError(
            f"exactly one no_tool ablation is required, got {ablations}; the "
            "leakage check is computed against it")
    references = [e["name"] for e in declared if e.get("reference")]
    if len(references) != 1:
        raise SystemConfigError(
            f"exactly one reference system is required, got {references}; the "
            "difficulty gate is judged on its per-rung profile")
    tooled = [e["name"] for e in declared if e.get("mode", "tooled") == "tooled"]
    if not tooled:
        raise SystemConfigError("no tooled system declared")
    # The R1 clear rate is computed over GENERIC tooled systems -- tooled and not
    # the reference -- because the reference is chosen to be strong and would
    # flatter the rate it is meant to police. A config without one makes that
    # gate uncomputable, and the sweep reports it as a failure with a rate of
    # None, which reads as a campaign defect rather than a missing system.
    generic = [name for name in tooled if name not in references]
    if not generic:
        raise SystemConfigError(
            "every tooled system is the reference; at least one more is needed, "
            "because the R1 clear-rate gate is measured over the tooled systems "
            "that are NOT the reference")
    if references[0] in ablations:
        raise SystemConfigError(
            f"{references[0]} is both the reference system and the ablation")
    return declared


def real_systems(campaign: pathlib.Path, config: pathlib.Path) -> list[SystemSpec]:
    """The declared systems as SystemSpecs for one campaign.

    The timeout comes from the campaign, not the config, so a system cannot be
    given more wall clock than the campaign budgeted. A config may ask for less.
    """
    declared = load(config)
    budget = _budget_seconds(pathlib.Path(campaign))
    specs = []
    for entry in declared:
        mode = entry.get("mode", "tooled")
        env = dict(entry.get("env") or {})
        if mode == "no_tool":
            env["PATH"] = _without_tool_prefix(
                env.get("PATH") or os.environ.get("PATH", ""))
        asked = int(entry.get("timeout_s") or budget)
        specs.append(SystemSpec(
            name=entry["name"],
            command=[str(part) for part in entry["command"]],
            mode=mode,
            timeout_s=min(asked, budget),
            env=env,
            is_reference=bool(entry.get("reference")),
            is_stub=False,
        ))
    return specs
