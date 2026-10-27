"""Manifest-driven staging of raw data.

CLAUDE.md: raw data never enters git. The manifest holds the URL or accession,
the byte count and the sha256; a fetch verifies both before the file counts as
staged. If an accession stops resolving, the manifest is the record of what we
had, so a failure names the entry rather than quietly staging less.

A partial download is deleted rather than left on disk. Half a genome that
passes no checksum is worse than no genome, because the next run would skip it
if we kept it.
"""

from __future__ import annotations

import hashlib
import json
import os
import ssl
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CHUNK = 1 << 20


@dataclass(frozen=True)
class Entry:
    path: str                   # relative to the task directory
    bytes: int
    sha256: str
    url: str | None = None
    accession: str | None = None
    release: str | None = None   # required for anything under reference/

    @property
    def source(self) -> str:
        return self.url or self.accession or "<no source>"


class FetchError(RuntimeError):
    """A manifest entry could not be staged. Never silently tolerated."""


def load_manifest(path: Path) -> list[Entry]:
    raw = json.loads(path.read_text())
    if not isinstance(raw, list):
        raise FetchError(f"{path}: manifest must be a JSON list")
    entries = []
    for item in raw:
        known = {f for f in Entry.__dataclass_fields__}
        unexpected = set(item) - known
        if unexpected:
            raise FetchError(
                f"{path}: entry {item.get('path')!r} has unknown field(s) "
                f"{', '.join(sorted(unexpected))}"
            )
        entries.append(Entry(**item))
    return entries


def dump_manifest(entries: list[Entry], path: Path) -> None:
    """Write a manifest, dropping unset optional fields so the file stays readable."""
    records = []
    for entry in entries:
        record = {"path": entry.path}
        for field in ("url", "accession", "release"):
            if getattr(entry, field) is not None:
                record[field] = getattr(entry, field)
        record["bytes"] = entry.bytes
        record["sha256"] = entry.sha256
        records.append(record)
    path.write_text(json.dumps(records, indent=2) + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def _ssl_context() -> ssl.SSLContext:
    """Honour SSL_CERT_FILE so a session behind a TLS-terminating proxy works."""
    bundle = os.environ.get("SSL_CERT_FILE") or os.environ.get("REQUESTS_CA_BUNDLE")
    return ssl.create_default_context(cafile=bundle) if bundle else ssl.create_default_context()


def download(url: str, target: Path, timeout: int = 300) -> None:
    """Download to a sibling .part file, then move it into place."""
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "npbench-r/0.1"})
        with urllib.request.urlopen(request, timeout=timeout,
                                    context=_ssl_context()) as response:
            with partial.open("wb") as handle:
                while chunk := response.read(CHUNK):
                    handle.write(chunk)
    except Exception as exc:
        partial.unlink(missing_ok=True)
        raise FetchError(f"{url}: download failed: {type(exc).__name__}: {exc}") from exc
    partial.replace(target)


def check(task_dir: Path, entry: Entry) -> str | None:
    """Return None if the staged file matches the manifest, else why it does not."""
    path = task_dir / entry.path
    if not path.is_file():
        return "not staged"
    size = path.stat().st_size
    if size != entry.bytes:
        return f"size is {size}, manifest says {entry.bytes}"
    actual = sha256_file(path)
    if actual != entry.sha256:
        return f"sha256 is {actual}, manifest says {entry.sha256}"
    return None


def fetch(task_dir: Path, entries: list[Entry], force: bool = False) -> list[str]:
    """Stage every entry, verifying each. Returns a line per entry for the log.

    Raises on the first entry that cannot be staged. A task is either fully
    staged or it is not staged; there is no useful middle state.
    """
    log = []
    for entry in entries:
        target = task_dir / entry.path
        if not force and check(task_dir, entry) is None:
            log.append(f"ok       {entry.path} (already staged)")
            continue
        if entry.url is None:
            raise FetchError(
                f"{entry.path}: no URL, so it cannot be fetched automatically. "
                f"Accession {entry.accession!r} must be staged by hand"
            )
        download(entry.url, target)
        problem = check(task_dir, entry)
        if problem is not None:
            target.unlink(missing_ok=True)
            raise FetchError(
                f"{entry.path}: {problem}. Fetched from {entry.url}. The file at "
                f"that URL is not the file this task was built against"
            )
        log.append(f"fetched  {entry.path} ({entry.bytes} bytes)")
    return log


def describe(url: str, path: str, timeout: int = 300, **fields) -> Entry:
    """Download a file to a temporary location to record its size and checksum.

    Used when building a manifest for a new task, not at runtime.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / "payload"
        download(url, target, timeout=timeout)
        return Entry(
            path=path, bytes=target.stat().st_size, sha256=sha256_file(target),
            url=url, **fields,
        )
