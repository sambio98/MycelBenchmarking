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

#: Version strings already queried in this process, keyed as the memo keys them.
#: Not a cache of anything that can change: a tool does not change version while
#: one process runs.
_VERSIONS: dict[str, str] = {}


def memoize_version(parts: Sequence[str], compute: Callable[[], str]) -> str:
    """A tool's self-reported version, read once per process and once per sweep.

    Measured, because only two of the nine tools this benchmark drives are worth
    memoising: `antismash --version` costs **2.12s** and `bigscape --version`
    **2.38s**, both of them Python packages whose import dominates the answer,
    while hmmscan, hmmsearch, prodigal, mmseqs and diamond all answer in 10ms or
    less and are left to answer. One warm bgcdetect sweep run spent 4.4 of its
    4.5 seconds asking antiSMASH its version twice.

    Two scopes with two different justifications:

    - **Per process, always.** A tool cannot change version while one process
      runs, so this needs no opt-in and applies to a build, a test and a real
      agent's own invocation alike.
    - **On disk, only where the tool memo is enabled** -- a stub in tooled mode,
      per `memoize_json`. The key is a PROXY for what is installed: a command
      path and whatever the caller pins, not a digest of the package, because a
      Python tool reports a version from its installed library rather than from
      the wrapper script a digest would cover. Staleness therefore surfaces as a
      failed run rather than a wrong pass: gold carries the string the oracle's
      own invocation produced, so a memo that returned an older one makes the
      run fail the environment component. What keeps the pin honest is the
      registry's verify, which the readiness gate runs.
    """
    key = _key(["tool_version", *parts])
    if key in _VERSIONS:
        return _VERSIONS[key]
    value = memoize_json(["tool_version", *parts], compute)
    if not isinstance(value, str) or not value.strip():
        raise CacheError(
            f"a version query returned {value!r}; a memo entry that is not a "
            "non-empty string would be reported as a tool identity")
    _VERSIONS[key] = value
    return value
