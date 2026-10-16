"""Grader identity (NPBench-C spec 3.3).

Every reported score carries a ``grader_version``: a semver string plus a
content hash over the grading package. A grader change is a benchmark version
bump, and a score without this string is not a comparable score.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

GRADER_SEMVER = "2.0.0"

_SCORING_PATH_FILES = (
    "primitives.py",
    "compose.py",
    "ladder.py",
    "grade.py",
    "version.py",
)


def content_hash() -> str:
    """SHA-256 over the scoring-path sources, in fixed filename order.

    Only files that can influence a score are hashed. Directory listing order
    is never consulted -- the file tuple above is explicit so the hash is stable
    across filesystems.
    """
    h = hashlib.sha256()
    here = Path(__file__).resolve().parent
    for name in _SCORING_PATH_FILES:
        h.update(name.encode("utf-8"))
        h.update(b"\x00")
        h.update((here / name).read_bytes())
        h.update(b"\x00")
    return h.hexdigest()[:8]


def grader_version() -> str:
    return f"{GRADER_SEMVER}+{content_hash()}"
