"""Regenerate image/environment.lock.json from the prefix that was built.

Usage:
  python -m npbench_c.tools.lockfile [--prefix /opt/npbench-tools] [--out <path>]
  python -m npbench_c.tools.lockfile --check     # regenerate and diff, write nothing

The lock records the closure the pins resolved to: every package with its
artifact URL, sha256 and size, plus how the prefix was asked for. It is a
record AND a build input: `--explicit` renders it as a conda explicit file, and
the Dockerfile creates the prefix from that rather than from the eleven
requested specs. It used to build from the specs and the difference was
measured -- a fresh solve of the same eleven pins produced 232 packages against
the lock's 231, with brotli, curl, cpython, fonttools, ca-certificates and a
dozen more transitives at newer versions than the lock records. The eleven pins
held; everything underneath them drifted, which is the whole of what a lock is
for. So the one thing it must never be is a hand-maintained approximation of a
prefix nobody re-read. It was hand-assembled once and immediately went stale: the version
after it said, in prose, "THE NEXT BUILD MUST RE-SOLVE AND REGENERATE THIS FILE".
This module is that regeneration, so the next person runs a command instead of
reading an instruction.

Everything mechanical is read from `$PREFIX/conda-meta`: one JSON per installed
package carrying exactly the fields the lock wants, and a `history` file
recording each create and install command with its timestamp. Everything
editorial -- the description, the rationale for each pin, the declared absences
-- is carried over from the existing lock, because it is writing rather than
measurement and regenerating it would silently delete it.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys

DEFAULT_PREFIX = pathlib.Path("/opt/npbench-tools")
DEFAULT_OUT = (pathlib.Path(__file__).resolve().parents[3]
               / "image" / "environment.lock.json")

#: Fields the lock carries that are prose rather than measurement. Preserved
#: from the file being replaced; a regeneration that dropped them would delete
#: the reasoning behind every pin.
EDITORIAL = ("description", "requested_rationale", "not_in_this_image",
             "banned_binaries")

#: Fields this module derives from the prefix, and therefore overwrites.
DERIVED = ("artifact_bytes", "build_history", "channels", "package_count",
           "packages", "platform", "prefix", "requested", "resolved_on",
           "solver")

SPEC = re.compile(r"^([A-Za-z0-9_.-]+)=([^=\s]+)=([^=\s]+)$")


class LockError(RuntimeError):
    """The prefix cannot be described, so the lock would be a guess."""


def _packages(prefix: pathlib.Path) -> list[dict]:
    meta = prefix / "conda-meta"
    if not meta.is_dir():
        raise LockError(f"{prefix} has no conda-meta directory to read")
    out = []
    for path in sorted(meta.glob("*.json")):
        document = json.loads(path.read_text())
        missing = [k for k in ("name", "version", "build", "url", "sha256", "size")
                   if k not in document]
        if missing:
            raise LockError(f"{path.name} records no {missing}")
        out.append({key: document[key] for key in
                    ("build", "name", "sha256", "size", "url")}
                   | {"version": document["version"]})
    if not out:
        raise LockError(f"{meta} holds no package records")
    return sorted(out, key=lambda p: (p["name"], p["version"], p["build"]))


def _history(prefix: pathlib.Path) -> list[dict]:
    """Each create/install command recorded for this prefix, with its date.

    The Dockerfile presents one `micromamba create` with the whole spec list.
    The prefix was actually built in stages, and saying so here is the point: a
    closure assembled over five solves is not guaranteed to equal the closure one
    solve would produce, and the lock is where that is visible rather than
    inferred.
    """
    path = prefix / "conda-meta" / "history"
    if not path.is_file():
        return []
    out, when = [], None
    for line in path.read_text(errors="replace").splitlines():
        stamp = re.match(r"==>\s*(.+?)\s*<==", line)
        if stamp:
            when = stamp.group(1)
            continue
        if line.startswith("# cmd:"):
            out.append({"at": when, "command": " ".join(line[6:].split())})
    return out


def _requested(history: list[dict]) -> dict[str, str]:
    """Every `name=version=build` spec asked of this prefix, last one winning."""
    requested: dict[str, str] = {}
    for entry in history:
        for token in entry["command"].split():
            match = SPEC.match(token)
            if match:
                requested[match.group(1)] = f"{match.group(2)}={match.group(3)}"
    return dict(sorted(requested.items()))


def _channels(packages: list[dict]) -> list[str]:
    seen = []
    for package in packages:
        match = re.search(r"https://conda\.anaconda\.org/([^/]+)/", package["url"])
        if match and match.group(1) not in seen:
            seen.append(match.group(1))
    return seen


def _platform(packages: list[dict]) -> str:
    subdirs = {package["url"].rstrip("/").split("/")[-2] for package in packages}
    concrete = sorted(subdirs - {"noarch"})
    if len(concrete) != 1:
        raise LockError(
            f"the prefix mixes platforms {concrete}; a lock that named one of "
            "them would describe a prefix nobody can rebuild")
    return concrete[0]


def _solver() -> str:
    try:
        done = subprocess.run(["micromamba", "--version"],
                              capture_output=True, text=True, timeout=60)
        version = (done.stdout or "").strip().splitlines()
        return f"micromamba {version[0]}" if version else "micromamba"
    except (OSError, subprocess.TimeoutExpired):
        return "micromamba"


def build(prefix: pathlib.Path, existing: dict) -> dict:
    packages = _packages(prefix)
    history = _history(prefix)
    document = {key: existing[key] for key in EDITORIAL if key in existing}
    document.update({
        "artifact_bytes": sum(int(p["size"]) for p in packages),
        "build_history": history,
        "channels": _channels(packages),
        "package_count": len(packages),
        "packages": packages,
        "platform": _platform(packages),
        "prefix": str(prefix),
        "requested": _requested(history),
        "resolved_on": (history[-1]["at"].split()[0] if history
                        else existing.get("resolved_on")),
        "solver": _solver(),
    })
    return document


EXPLICIT_HEADER = """# Generated from environment.lock.json by
# `python -m npbench_c.tools.lockfile --explicit`. Do not edit: every line is a
# recorded artifact URL with the sha256 the lock recorded for it, so a rebuild
# links the closure the lock describes instead of whatever the solver resolves
# today. Regenerate it when the lock changes; a test fails if the two disagree.
@EXPLICIT
"""


def explicit(document: dict) -> str:
    """The lock as a conda explicit file: no solver, one artifact per line."""
    packages = document.get("packages") or []
    if not packages:
        raise LockError("the lock carries no packages, so there is nothing to pin")
    missing = [p["name"] for p in packages if not p.get("url") or not p.get("sha256")]
    if missing:
        raise LockError(
            f"{len(missing)} package(s) carry no url or no sha256, so an explicit "
            f"file would be unverifiable: {sorted(missing)[:5]}")
    lines = [f"{p['url']}#{p['sha256']}"
             for p in sorted(packages, key=lambda p: p["name"])]
    return EXPLICIT_HEADER + "\n".join(lines) + "\n"


DEFAULT_EXPLICIT = (pathlib.Path(__file__).resolve().parents[3]
                    / "image" / "environment.explicit.txt")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default=str(DEFAULT_PREFIX))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--check", action="store_true",
                        help="regenerate in memory and report what would change")
    parser.add_argument("--explicit", nargs="?", const=str(DEFAULT_EXPLICIT),
                        default=None,
                        help="render the committed lock as the conda explicit "
                             "file the image builds from")
    args = parser.parse_args(argv)

    if args.explicit:
        source = pathlib.Path(args.out)
        target = pathlib.Path(args.explicit)
        rendered = explicit(json.loads(source.read_text()))
        if args.check:
            same = target.is_file() and target.read_text() == rendered
            print(f"{target}: {'up to date' if same else 'STALE'}")
            return 0 if same else 1
        target.write_text(rendered)
        print(f"wrote {target} ({len(rendered.splitlines()) - 6} artifacts)")
        return 0

    out = pathlib.Path(args.out)
    existing = json.loads(out.read_text()) if out.is_file() else {}
    document = build(pathlib.Path(args.prefix), existing)
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"

    if args.check:
        same = out.is_file() and out.read_text() == rendered
        print(f"{out}: {'up to date' if same else 'STALE'} "
              f"({document['package_count']} packages, "
              f"{document['artifact_bytes']} bytes)")
        return 0 if same else 1

    out.write_text(rendered)
    print(f"wrote {out}")
    print(f"  {document['package_count']} packages, "
          f"{document['artifact_bytes']} artifact bytes, "
          f"{len(document['build_history'])} build step(s)")
    print(f"  requested: {', '.join(sorted(document['requested']))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
