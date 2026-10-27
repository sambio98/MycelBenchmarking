"""Manifest staging, served from a local HTTP server so no network is needed."""

import functools
import hashlib
import http.server
import json
import threading

import pytest

from npbench.build.manifest import (
    Entry,
    FetchError,
    check,
    describe,
    dump_manifest,
    fetch,
    load_manifest,
    sha256_file,
)

PAYLOAD = b">contig_1\n" + b"ACGT" * 4096


@pytest.fixture
def origin(tmp_path):
    """Serve tmp_path/'origin' over HTTP and yield its base URL."""
    root = tmp_path / "origin"
    root.mkdir()
    (root / "genome.fna").write_bytes(PAYLOAD)
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_port}", root
    server.shutdown()


@pytest.fixture
def entry(origin):
    base, _ = origin
    return Entry(
        path="data/genome.fna",
        bytes=len(PAYLOAD),
        sha256=hashlib.sha256(PAYLOAD).hexdigest(),
        url=f"{base}/genome.fna",
    )


def test_fetch_stages_and_verifies(tmp_path, entry):
    task = tmp_path / "task"
    log = fetch(task, [entry])
    assert (task / "data/genome.fna").read_bytes() == PAYLOAD
    assert "fetched" in log[0]
    assert check(task, entry) is None


def test_second_fetch_skips_a_verified_file(tmp_path, entry):
    task = tmp_path / "task"
    fetch(task, [entry])
    assert "already staged" in fetch(task, [entry])[0]


def test_force_refetches(tmp_path, entry):
    task = tmp_path / "task"
    fetch(task, [entry])
    assert "fetched" in fetch(task, [entry], force=True)[0]


def test_checksum_mismatch_raises_and_removes_the_file(tmp_path, entry):
    # The manifest says one thing and the URL now serves another: the exact
    # record-drift case frozen bundles exist to catch.
    wrong = Entry(**{**entry.__dict__, "sha256": "0" * 64})
    task = tmp_path / "task"
    with pytest.raises(FetchError, match="not the file this task was built against"):
        fetch(task, [wrong])
    assert not (task / "data/genome.fna").exists()


def test_size_mismatch_raises(tmp_path, entry):
    wrong = Entry(**{**entry.__dict__, "bytes": 11})
    with pytest.raises(FetchError, match="size is"):
        fetch(tmp_path / "task", [wrong])


def test_missing_url_names_the_accession(tmp_path):
    entry = Entry(path="data/x.fna", bytes=1, sha256="a" * 64, accession="IMG-2512875010")
    with pytest.raises(FetchError, match="IMG-2512875010.*must be staged by hand"):
        fetch(tmp_path / "task", [entry])


def test_unreachable_url_raises_without_leaving_a_partial(tmp_path):
    entry = Entry(
        path="data/x.fna", bytes=1, sha256="a" * 64,
        url="http://127.0.0.1:1/nothing",
    )
    task = tmp_path / "task"
    with pytest.raises(FetchError, match="download failed"):
        fetch(task, [entry])
    assert not list(task.rglob("*.part"))


def test_check_reports_each_failure_mode(tmp_path, entry):
    task = tmp_path / "task"
    assert check(task, entry) == "not staged"
    fetch(task, [entry])
    (task / "data/genome.fna").write_bytes(PAYLOAD + b"X")
    assert "size is" in check(task, entry)
    (task / "data/genome.fna").write_bytes(b"Y" * len(PAYLOAD))
    assert "sha256 is" in check(task, entry)


def test_describe_records_size_and_checksum(origin):
    base, _ = origin
    entry = describe(f"{base}/genome.fna", "data/genome.fna")
    assert entry.bytes == len(PAYLOAD)
    assert entry.sha256 == hashlib.sha256(PAYLOAD).hexdigest()


def test_manifest_round_trip(tmp_path, entry):
    path = tmp_path / "data.manifest.json"
    dump_manifest([entry], path)
    assert load_manifest(path) == [entry]
    # Optional fields that were never set stay out of the file.
    assert "release" not in json.loads(path.read_text())[0]


def test_unknown_manifest_field_is_rejected(tmp_path):
    path = tmp_path / "data.manifest.json"
    path.write_text(json.dumps([{"path": "a", "bytes": 1, "sha256": "a" * 64, "size_mb": 1}]))
    with pytest.raises(FetchError, match="unknown field"):
        load_manifest(path)


def test_sha256_file_matches_hashlib(tmp_path):
    path = tmp_path / "x"
    path.write_bytes(PAYLOAD)
    assert sha256_file(path) == hashlib.sha256(PAYLOAD).hexdigest()


# The environment lock generator lives in scripts/, but the marker logic decides
# what ends up pinned, so it is worth a test: following extras-gated
# requirements once pulled pytest, setuptools and wheel into a task lock.
def test_requirement_name_skips_extras_gated_requirements():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from lock_task_environment import requirement_name

    assert requirement_name("pyhmmer ~=0.10") == "pyhmmer"
    assert requirement_name("scikit-learn<2,>=1.0") == "scikit-learn"
    assert requirement_name("rich-argparse") == "rich-argparse"
    assert requirement_name("numpy; python_version >= '3.9'") == "numpy"
    assert requirement_name('pytest; extra == "test"') is None
    assert requirement_name("wheel ; extra == 'dev'") is None
    assert requirement_name(
        'wheel ; (os_name == "nt" and implementation_name != "pypy") and extra == \'test\''
    ) is None
    assert requirement_name("statsmodels<0.15,>=0.13; extra == \"train\"") is None
