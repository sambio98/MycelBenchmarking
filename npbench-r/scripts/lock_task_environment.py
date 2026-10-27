"""Write a task's environment/requirements.lock from the installed runtime closure.

CLAUDE.md: environments are pinned or they are wrong. The lock must be the
closure actually needed at runtime, so extras-gated requirements
(`extra == "test"`, `extra == "dev"`) are skipped. Following them is how pytest,
setuptools and wheel ended up in the first version of this file.

    python scripts/lock_task_environment.py tasks/02_salinispora_bgc gecco gecco-tool
"""

from __future__ import annotations

import argparse
import importlib.metadata as md
import re
import sys
from pathlib import Path

HEADER = """\
# Pinned runtime environment for {task}.
#
# CLAUDE.md: environments are pinned or they are wrong. These are the exact
# versions truth/ was produced under. GECCO's cluster counts and type
# assignments are version-dependent, so truth is only valid against these.
#
# Regenerate with:
#   python scripts/lock_task_environment.py {task} {roots}
# then regenerate truth/ and re-run `make selfgrade`, because a version bump
# invalidates the old truth.
#
# Install with:
#   pip install --no-deps -r environment/requirements.lock
"""


def requirement_name(spec: str) -> str | None:
    """The distribution name of a requirement, or None if it is extras-gated."""
    requirement, _, marker = spec.partition(";")
    if re.search(r"\bextra\s*==", marker):
        return None
    name = re.split(r"[\[<>=!~ (]", requirement.strip(), maxsplit=1)[0].strip()
    return name or None


def closure(roots: list[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    pending = list(roots)
    while pending:
        name = pending.pop()
        key = name.lower().replace("_", "-")
        if key in found:
            continue
        try:
            dist = md.distribution(name)
        except md.PackageNotFoundError:
            print(f"  (not installed, skipped: {name})", file=sys.stderr)
            continue
        found[key] = dist.version
        for spec in dist.requires or []:
            dependency = requirement_name(spec)
            if dependency:
                pending.append(dependency)
    return found


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("task_dir", type=Path)
    parser.add_argument("roots", nargs="+", help="top-level distributions to resolve")
    args = parser.parse_args()

    resolved = closure(args.roots)
    target = args.task_dir / "environment" / "requirements.lock"
    target.write_text(
        HEADER.format(task=args.task_dir, roots=" ".join(args.roots))
        + "".join(f"{name}=={version}\n" for name, version in sorted(resolved.items()))
    )
    print(f"wrote {target}: {len(resolved)} packages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
