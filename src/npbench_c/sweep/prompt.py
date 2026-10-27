"""The instructions a real system is given, built from the cold sandbox alone.

A stub fixture reads `task.yaml` itself, because it is code in this repository.
A real agent cannot: it needs to be told what to do, where to write it, and what
it is allowed to use. Nothing in the harness produced that, which is why every
sweep on record is stub-based -- the missing piece was never the runner (a system
is a command) but the prompt.

Three properties make this generator worth having rather than a prompt per
campaign:

**It is built from the SANDBOX, not from the campaign.** The only inputs are the
files a system already has -- `task.yaml`, the names under `reference/` and
`inputs/`, and any warm-start bundle handed over. It cannot leak gold because it
never opens it, and a test asserts the rendered prompt against the gate's own
cold-sandbox scans.

**It is the same for every system.** One campaign, one rung mode, one prompt,
byte-identical for every system that runs it. A difficulty comparison between two
systems given different instructions measures the instructions.

**It states the contract and declines to state the method.** `method_guidance:
partial` across this suite means the campaign says WHAT to report and leaves HOW
to the system. So the prompt renders the objective, the scope note, the output
contract field by field, the pinned tables, the budget and the tool surface --
and never a worked procedure.
"""

from __future__ import annotations

import json
import pathlib
from typing import Mapping

import yaml

PROMPT = "PROMPT.md"

#: Keys rendered as prose sections, in this order, with the heading to use.
#: Everything else in `task.yaml` is either machinery (harness, grading, rungs)
#: or an auditor's note that would read as a hint.
SECTIONS = (
    ("objective", "What to do"),
    ("scope_note", "What is being graded"),
    ("catalog_deviation", "What this campaign does not ask"),
    ("abstention_note", "On abstention"),
    ("control_pair_note", "On abstention"),
    ("gold_construction", "How the answer is established"),
    ("method_guidance", None),          # rendered inside "How you are judged"
)

#: Never rendered. `notes_for_audit` is written for someone holding gold, and
#: the rest is harness wiring.
WITHHELD = ("notes_for_audit", "harness", "grading", "rungs", "gold_tier",
            "grader_version", "grading_spec_version", "template_params",
            "required_sibling", "excluded_from_grading", "license_notice")


class PromptError(RuntimeError):
    """The sandbox does not carry what a prompt needs."""


def _wrap(text: str, width: int = 88) -> str:
    """Re-flow a YAML folded block into paragraphs, keeping blank lines."""
    import textwrap

    out = []
    for block in str(text).strip().split("\n\n"):
        collapsed = " ".join(block.split())
        out.append("\n".join(textwrap.wrap(collapsed, width)) if collapsed else "")
    return "\n\n".join(out)


def _tree(sandbox: pathlib.Path, name: str) -> list[str]:
    directory = sandbox / name
    if not directory.is_dir():
        return []
    return sorted(f"{name}/{p.relative_to(directory)}"
                  for p in directory.rglob("*") if p.is_file())


#: What the sandbox always has at its root, bundle or no bundle.
FURNITURE = ("task.yaml", PROMPT)


def _warm_files(sandbox: pathlib.Path) -> list[str]:
    """A warm bundle's files, which arrive at the sandbox root beside task.yaml.

    Any file type: the first draft looked for `*.json` and raised on the two
    construct campaigns, whose R2 bundle hands over a FASTA. A bundle is whatever
    the campaign pre-rendered, and the prompt's job is to name the files rather
    than to have an opinion about their format.
    """
    return sorted(p.name for p in sandbox.iterdir()
                  if p.is_file() and p.name not in FURNITURE)


def _fields(fields: Mapping, depth: int = 0) -> list[str]:
    """The contract, as a nested bullet list.

    A contract entry's description is sometimes itself a mapping -- a group of
    sub-fields under one heading, which is how the larger reports are declared.
    Rendering that with `str()` put a Python dict repr in front of the system,
    which is how the first draft of this generator read.
    """
    out = []
    pad = "  " * depth
    for field, description in fields.items():
        if isinstance(description, Mapping):
            out.append(f"{pad}- **`{field}`**")
            out += _fields(description, depth + 1)
        else:
            text = " ".join(str(description).split())
            out.append(f"{pad}- **`{field}`** — {text}")
    return out


def render(sandbox: pathlib.Path, rung_mode: str = "cold") -> str:
    """The prompt for one sandbox. Deterministic; reads nothing but the sandbox."""
    sandbox = pathlib.Path(sandbox)
    task_path = sandbox / "task.yaml"
    if not task_path.is_file():
        raise PromptError(f"{sandbox} carries no task.yaml, so there is no task")
    task = yaml.safe_load(task_path.read_text())

    contract = task.get("output_contract") or {}
    if not contract:
        raise PromptError(
            f"{task.get('campaign_id')} declares no output_contract, so a system "
            "cannot know what to produce")

    lines = [f"# {task.get('campaign_id', sandbox.name)}", ""]
    lines.append(_wrap(
        "You are given a working directory and asked to produce one or more JSON "
        "reports in it. Everything you need is in this directory; there is no "
        "network. Read `task.yaml` as well -- this file is a rendering of it, and "
        "where the two differ, `task.yaml` is the contract."))
    lines.append("")

    for key, heading in SECTIONS:
        if key == "method_guidance" or key not in task:
            continue
        lines += [f"## {heading}", "", _wrap(task[key]), ""]

    # ---- what to write, and where
    lines += ["## What to produce", ""]
    lines.append(_wrap(
        "Write each artefact below into the `submission/` directory of this "
        "working directory, as JSON. An absent artefact scores zero for the whole "
        "run, and so does an empty `submission/`. Report every field: a field you "
        "cannot compute should be present and null rather than missing, so that "
        "the difference between 'no answer' and 'this answer' is yours to state "
        "rather than ours to infer."))
    lines.append("")
    for artefact, fields in contract.items():
        lines += [f"### `submission/{artefact}`", ""]
        lines += _fields(fields or {})
        lines.append("")

    # ---- what is in the directory
    inputs, reference = _tree(sandbox, "inputs"), _tree(sandbox, "reference")
    lines += ["## What is in this directory", ""]
    for name, listing in (("inputs", inputs), ("reference", reference)):
        if listing:
            lines.append(f"- `{name}/` — " + ", ".join(f"`{p}`" for p in listing))
    pinned = (task.get("constraints") or {}).get("pinned_tables") or []
    if pinned:
        lines.append("")
        lines.append(_wrap(
            "The pinned tables below are part of the contract rather than "
            "reference material: the graded answers are defined in terms of them, "
            "so a different table or a different convention is a different "
            "answer. Use them as given."))
        lines.append("")
        for table in ([pinned] if isinstance(pinned, str) else pinned):
            lines.append(f"- `{table}`")
    external = (task.get("constraints") or {}).get("external_resources_note")
    if external:
        lines += ["", _wrap(external)]
    lines.append("")

    # ---- warm start
    handed = _warm_files(sandbox)
    if rung_mode != "cold":
        if not handed:
            raise PromptError(
                f"{rung_mode} is a warm run and no bundle reached {sandbox}; a "
                "warm prompt that hands over nothing is a cold run mislabelled")
        lines += ["## What you have been given", ""]
        lines.append(_wrap(
            "The JSON files at the root of this directory carry verified answers "
            "from the earlier stages of this task: " +
            ", ".join(f"`{name}`" for name in handed) + ". Take them as correct. "
            "You are not being asked to redo that work -- produce the full report, "
            "copying those answers where they cover a field, and compute the rest."))
        lines.append("")

    # ---- how it is judged
    lines += ["## How you are judged", ""]
    judged = [
        "Your report is read by a deterministic measurer and compared against a "
        "reference answer computed in advance. No model scores you and no "
        "partial credit is given for a field: a value either matches or it does "
        "not.",
        "The task is graded in four stages, each a conjunction of its own "
        "claims, and a stage counts only if every stage below it also holds. So "
        "the reporting basics matter as much as the hardest number: a malformed "
        "field costs the stages above it.",
    ]
    if task.get("method_guidance") == "partial":
        judged.append(
            "How you compute any of this is yours to choose. The contract fixes "
            "what to report and the conventions to report it under, not the "
            "method.")
    for paragraph in judged:
        lines += [_wrap(paragraph), ""]

    budgets = task.get("budgets") or {}
    if budgets:
        parts = []
        if budgets.get("wall_clock_minutes"):
            parts.append(f"{budgets['wall_clock_minutes']} minutes of wall clock")
        if budgets.get("cpu"):
            parts.append(f"{budgets['cpu']} CPU")
        if budgets.get("memory_gb"):
            parts.append(f"{budgets['memory_gb']} GB of memory")
        parts.append("no network")
        lines += ["## Budget", "", _wrap(
            "This task is budgeted for " + ", ".join(parts) + ". The budget is "
            "declared so you can plan: a run that spends it all on one stage and "
            "reports nothing scores zero."), ""]

    surface = task.get("tool_surface_note")
    if surface:
        lines += ["## Tools", "", _wrap(surface), ""]

    return "\n".join(lines).rstrip() + "\n"


def write(sandbox: pathlib.Path, rung_mode: str = "cold") -> pathlib.Path:
    """Render the prompt into the sandbox and return its path."""
    path = pathlib.Path(sandbox) / PROMPT
    path.write_text(render(sandbox, rung_mode))
    return path


def main(argv: list[str] | None = None) -> int:
    """Render one campaign's cold prompt to stdout, for review."""
    import argparse
    import sys
    import tempfile

    from npbench_c.sweep.runner import build_sandbox

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("campaign")
    ap.add_argument("--rung", default=None,
                    help="render the warm prompt for this rung id, e.g. r3")
    args = ap.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    bundle = None
    if args.rung:
        task = yaml.safe_load((campaign / "task.yaml").read_text())
        bundle = next((r.get("warm_start_bundle") for r in task.get("rungs", [])
                       if r["rung_id"] == args.rung), None)
        if not bundle:
            print(f"{campaign.name} declares no warm bundle for {args.rung}",
                  file=sys.stderr)
            return 2
    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(campaign, pathlib.Path(tmp) / "cold", bundle)
        print(render(sandbox, "cold" if not args.rung else f"warm:{args.rung}"),
              end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
