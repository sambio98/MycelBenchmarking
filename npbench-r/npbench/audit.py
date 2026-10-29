"""Checks for the human curation pass, and the per-task dossier it reads from.

`npbench validate` answers "does this bundle satisfy the contract". This module
answers a different question: "would a reviewer, reading this bundle, find
anything that makes the task unfair, unreproducible or already answered". The
three checks here are the ones that were found by hand during the build and
would otherwise have to be found by hand again:

- **Answer leakage.** CLAUDE.md: do not leak answers. A prompt that names a row
  count, a family count or any other result hands over part of the grade. The
  prompts describe row universes qualitatively for exactly this reason, so a
  number in a prompt that equals a number in truth/ is worth a look.
- **Graded but unstated columns.** `validate` already rejects an artifact that
  grading.yaml grades and prompt.md never requests. The same unfairness applies
  one level down, per column: a column that is graded but never named in the
  prompt is a column the agent was not asked for.
- **Stated but ungraded columns.** The mirror image, and not unfair, only waste:
  a column the prompt demands that no metric reads. It was how a degenerate
  artifact survived in task 01 until it was measured.

Nothing here is a contract violation by itself, so the checks report findings
for a person to judge rather than failing a build.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .grading.artifact import declared_columns, read_table
from .grading.config import GradingConfig, load_grading_config

# A number as a prompt would write one, thousands separators included.
NUMBER = re.compile(r"\b\d{1,3}(?:,\d{3})+\b|\b\d+(?:\.\d+)?\b")

# A number matching some cell in some truth file is weak evidence: a prompt is
# full of tolerances, thresholds and decimal places, and every one of them
# matches a cell somewhere. So that check only considers numbers large enough
# to be unlikely coincidences.
CELL_FLOOR = 100

# A number matching an artifact's row count is strong evidence whatever its
# magnitude, so it gets a much lower floor -- just high enough to not fire on
# the "1-based", "top 10", "at least 6" of ordinary method guidance.
ROW_COUNT_FLOOR = 10


@dataclass(frozen=True)
class Finding:
    level: str  # "leak" | "unfair" | "waste"
    where: str
    message: str

    def __str__(self) -> str:
        return f"{self.level.upper():7} {self.where}: {self.message}"


def acknowledged(meta: dict) -> list[dict]:
    """Findings a reviewer has already judged, from task.yaml.

    A leak check cannot tell a result from a number the agent could have worked
    out for itself: that 120 genomes are staged is countable, and that they make
    7,140 pairs is arithmetic. Both look exactly like a leaked row count. So the
    reviewer's judgement is recorded in the bundle once, against the message
    text rather than a line number, and the finding stops competing for
    attention with the ones nobody has looked at yet.
    """
    entries = meta.get("audit_acknowledged") or []
    out = []
    for entry in entries:
        if isinstance(entry, dict) and "contains" in entry:
            out.append({"contains": str(entry["contains"]),
                        "why": str(entry.get("why", ""))})
    return out


@dataclass
class Dossier:
    task_id: str
    title: str = ""
    depth: object = None
    categories: list[str] = field(default_factory=list)
    artifacts: list[tuple[str, bool, int]] = field(default_factory=list)
    staged_bytes: int = 0
    truth_rows: dict[str, int] = field(default_factory=dict)
    findings: list[Finding] = field(default_factory=list)
    settled: list[tuple[Finding, str]] = field(default_factory=list)
    scope_note: str = ""

    @property
    def n_hard(self) -> int:
        return sum(1 for _, soft, _ in self.artifacts if not soft)

    @property
    def n_metrics(self) -> int:
        return sum(n for _, _, n in self.artifacts)


def prompt_sections(text: str) -> dict[str, str]:
    """The prompt's `### \\`file.csv\\`` blocks, keyed by artifact path.

    Every task's prompt specifies its outputs this way, one heading per
    artifact, so the per-column checks can be scoped to the right block rather
    than searching the whole prompt and matching a column name that happens to
    appear in the method guidance.
    """
    blocks: dict[str, str] = {}
    current = None
    for line in text.splitlines():
        heading = re.match(r"^###\s+`([^`]+)`\s*$", line)
        if heading:
            current = heading.group(1)
            blocks[current] = ""
        elif line.startswith("## "):
            current = None
        elif current is not None:
            blocks[current] += line + "\n"
    return blocks


def truth_values(task_dir: Path, config: GradingConfig) -> tuple[set[str], dict[str, int]]:
    """Every cell in truth/, plus each artifact's row count.

    The row counts are kept separately because they are the leak that matters
    most: a prompt that states how many rows an artifact has told the agent
    whether its filter was right.
    """
    cells: set[str] = set()
    rows: dict[str, int] = {}
    for artifact in config.artifacts:
        path = task_dir / "truth" / artifact.truth
        if not path.is_file():
            continue
        frame = read_table(path, artifact.separator)
        rows[artifact.path] = len(frame)
        for column in frame.columns:
            cells.update(frame[column].astype(str).str.strip())
    return cells, rows


def check_leakage(task_dir: Path, config: GradingConfig) -> list[Finding]:
    prompt = task_dir / "prompt.md"
    if not prompt.is_file():
        return []
    text = prompt.read_text()
    cells, rows = truth_values(task_dir, config)
    by_count = {count: name for name, count in rows.items()}

    findings = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        for token in NUMBER.findall(line):
            bare = token.replace(",", "")
            try:
                value = float(bare)
            except ValueError:
                continue
            if (value.is_integer() and value >= ROW_COUNT_FLOOR
                    and int(value) in by_count):
                findings.append(Finding(
                    "leak", f"prompt.md:{line_number}",
                    f"{token} is the row count of {by_count[int(value)]}; a prompt "
                    f"that states an artifact's size hands over its filter",
                ))
            elif value >= CELL_FLOOR and bare in cells:
                findings.append(Finding(
                    "leak", f"prompt.md:{line_number}",
                    f"{token} also appears as a value in truth/; check it is a "
                    f"method parameter and not a result",
                ))
    return findings


def required_outputs(text: str) -> str:
    """The prompt's required-outputs section, whatever its heading level."""
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if re.match(r"^#+\s", line) and "required output" in line.casefold():
            rest = []
            for following in lines[index + 1:]:
                if re.match(r"^#{1,2}\s", following):
                    break
                rest.append(following)
            return "\n".join(rest)
    return ""


def column_list(block: str) -> set[str]:
    """The columns an output block declares, from its column-list lines only.

    Each block states its schema as one or two lines holding nothing but
    backticked names separated by commas. Reading every backtick in the block
    instead was the first attempt and it was wrong three ways: it picked up the
    values a column may take (`start`, `end`, `both`), the input fields the
    prose refers to (`proteinogenic`, `CHARGE`) and the columns of *other*
    artifacts mentioned in passing (`cluster_id`). Restricting to the schema
    lines removes all three without a heuristic about case or content.
    """
    names: set[str] = set()
    for line in block.splitlines():
        stripped = line.strip()
        if not stripped or "`" not in stripped:
            continue
        # A schema line is backticked names and separators, nothing else.
        if re.fullmatch(r"(?:`[A-Za-z_][A-Za-z0-9_]*`[,\s]*)+", stripped):
            names.update(re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`", stripped))
    return names


def check_columns(task_dir: Path, config: GradingConfig) -> list[Finding]:
    prompt = task_dir / "prompt.md"
    if not prompt.is_file():
        return []
    text = prompt.read_text()
    blocks = prompt_sections(text)
    outputs_names = set(re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`",
                                   required_outputs(text)))

    findings = []
    for artifact in config.artifacts:
        block = blocks.get(artifact.path)
        graded = set(declared_columns(artifact))
        if block is None:
            # A prompt need not use one heading per artifact; the contract in
            # PROJECT.md §6 only asks that every column be stated somewhere in
            # the required-outputs section. Fall back to that whole section and
            # check the fairness direction only -- without per-artifact scoping
            # there is no way to tell which artifact an unused name belongs to.
            for column in sorted(graded - outputs_names):
                findings.append(Finding(
                    "unfair", artifact.path,
                    f"column {column!r} is graded but the prompt's required-outputs "
                    f"section never names it",
                ))
            continue
        named = column_list(block)
        for column in sorted(graded - named):
            findings.append(Finding(
                "unfair", f"{artifact.path}",
                f"column {column!r} is graded but the prompt's own output block "
                f"never names it",
            ))
        # A column the prompt names that nothing reads. Aliases are excluded:
        # the prompt may legitimately mention one.
        aliases = {a for spec in config.artifacts
                   for group in spec.column_aliases.values() for a in group}
        for column in sorted(named - graded - aliases):
            findings.append(Finding(
                "waste", f"{artifact.path}",
                f"the prompt names {column!r} in this artifact's block but no "
                f"metric reads it",
            ))
    return findings


def dossier(task_dir: Path) -> Dossier:
    """Everything a reviewer needs about one task, gathered in one place."""
    import json

    import yaml

    config = load_grading_config(task_dir / "grading.yaml")
    report = Dossier(task_id=task_dir.name)

    meta_path = task_dir / "task.yaml"
    meta: dict = {}
    if meta_path.is_file():
        meta = yaml.safe_load(meta_path.read_text()) or {}
        report.title = str(meta.get("title", ""))
        report.depth = meta.get("depth")
        report.categories = list(meta.get("categories") or [])
        report.scope_note = str(meta.get("data_scope") or "").strip()

    manifest = task_dir / "data.manifest.json"
    if manifest.is_file():
        report.staged_bytes = sum(
            int(e.get("bytes", 0)) for e in json.loads(manifest.read_text()))

    report.artifacts = [
        (a.path, a.soft, len(a.metrics)) for a in config.artifacts
    ]
    _, report.truth_rows = truth_values(task_dir, config)

    found = check_leakage(task_dir, config) + check_columns(task_dir, config)
    rules = acknowledged(meta if meta_path.is_file() else {})
    for finding in found:
        match = next((r for r in rules if r["contains"] in finding.message), None)
        if match is None:
            report.findings.append(finding)
        else:
            report.settled.append((finding, match["why"]))
    return report


def render(report: Dossier) -> str:
    lines = [f"# {report.task_id}"]
    if report.title:
        lines.append(f"  {report.title}")
    lines.append("")
    lines.append(f"  depth {report.depth}   categories: {', '.join(report.categories) or '-'}")
    lines.append(
        f"  {len(report.artifacts)} artifacts ({report.n_hard} hard), "
        f"{report.n_metrics} metric instances, {report.staged_bytes / 1e6:.1f} MB staged"
    )
    lines.append("")
    lines.append("  artifact                            rows  metrics")
    for path, soft, n_metrics in report.artifacts:
        rows = report.truth_rows.get(path)
        tag = "  [soft]" if soft else ""
        lines.append(
            f"    {path:32} {rows if rows is not None else '-':>7} {n_metrics:>8}{tag}"
        )
    if report.scope_note:
        lines.append("")
        lines.append("  scope note (task.yaml data_scope):")
        for chunk in _wrap(report.scope_note, 72):
            lines.append(f"    {chunk}")
    lines.append("")
    if report.findings:
        lines.append(f"  {len(report.findings)} finding(s) for review:")
        for finding in report.findings:
            lines.append(f"    {finding}")
    else:
        lines.append("  no findings for review")
    if report.settled:
        lines.append("")
        lines.append(f"  {len(report.settled)} already judged (task.yaml audit_acknowledged):")
        for finding, why in report.settled:
            lines.append(f"    {finding.level}: {finding.message}")
            lines.append(f"      -> {why}")
    return "\n".join(lines)


def _wrap(text: str, width: int) -> list[str]:
    words = text.split()
    out, line = [], ""
    for word in words:
        if len(line) + len(word) + 1 > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        out.append(line)
    return out
