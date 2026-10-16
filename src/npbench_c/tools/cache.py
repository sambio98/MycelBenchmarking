"""A memo for expensive tool output, enabled only where it cannot mislead.

An S1 or S2 campaign's stub sweep re-runs the same tool call for every stub
level, every rung and every repeat, over an input that never changes. Measured on
the first such campaign that is 29m32s of wall clock and 47m of CPU for a sweep
whose purpose is to check the grader and the ladder, not the tool: the same
HMMER output recomputed seventy-two times.

So the memo exists, and the whole design is about where it is allowed to apply.

**It is off unless `NPBENCH_TOOL_CACHE` names a directory.** No environment
variable, no caching: a build, a test, a real agent run and a human invocation
all compute from scratch. The sweep runner sets the variable for **stub systems
in tooled mode only** -- never for a real system, because a sweep over real
systems is partly a measurement of whether they can drive the tool at all, and
never for a `no_tool` ablation, because that run's whole claim is that it did not
compute.

**The key is the full identity of the computation**, including the tool version
and the fingerprint of every data file behind it. A memo keyed on less would
survive a change it should not survive -- a new Pfam release, a different
threshold -- and hand back an answer from the wrong world. Caller passes those
parts explicitly; this module will not guess them.

**The value is tool output derived from public inputs.** Nothing from `gold/` or
`oracle/` may be keyed or stored here, and the sweep's gold-isolation audit still
runs on every sandbox: the cache lives outside them.
"""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
from typing import Callable, Sequence

#: Names the cache directory. Unset means no caching anywhere.
ENV_VAR = "NPBENCH_TOOL_CACHE"

#: Bumped when the on-disk layout or the meaning of a key changes, so stale
#: entries from an older scheme are never read back as current.
CACHE_FORMAT = "1"


class CacheError(RuntimeError):
    """The cache was asked to do something that would make it a lie."""


def cache_dir() -> pathlib.Path | None:
    raw = os.environ.get(ENV_VAR, "").strip()
    return pathlib.Path(raw) if raw else None


def enabled() -> bool:
    return cache_dir() is not None


def file_digest(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _key(parts: Sequence[str]) -> str:
    if not parts or any(not isinstance(p, str) or not p for p in parts):
        raise CacheError(
            "every key part must be a non-empty string; a memo with a hole in "
            "its key returns an answer from the wrong world")
    payload = json.dumps([CACHE_FORMAT, *parts], sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def memoize_json(parts: Sequence[str], compute: Callable[[], object]) -> object:
    """Return `compute()`, reading and writing a JSON memo when one is enabled.

    `parts` must name every input the result depends on -- the tool version, the
    fingerprint of each data file, the flags. The value must be JSON-serialisable;
    anything richer should be converted by the caller, so what lands on disk is
    inspectable by whoever is wondering what a sweep actually ran.
    """
    # Validated before the early return, so a malformed key is a failure
    # everywhere -- in a build, in a test, on a developer's machine -- and not
    # only on the one run where caching happens to be switched on. A key defect
    # that surfaces only inside a sweep is a key defect nobody sees.
    key = _key(parts)

    directory = cache_dir()
    if directory is None:
        return compute()

    directory.mkdir(parents=True, exist_ok=True)
    entry = directory / f"{key}.json"
    if entry.is_file():
        try:
            stored = json.loads(entry.read_text())
        except json.JSONDecodeError:
            stored = None
        if isinstance(stored, dict) and stored.get("format") == CACHE_FORMAT:
            return stored["value"]

    value = compute()
    # Written through a temporary name so a concurrent reader never sees half an
    # entry, and so a crash mid-write leaves no entry rather than a truncated one.
    staging = entry.with_suffix(".json.partial")
    staging.write_text(json.dumps(
        {"format": CACHE_FORMAT, "key_parts": list(parts), "value": value},
        indent=2, sort_keys=True) + "\n")
    staging.replace(entry)
    return value
