"""The two costs a stub sweep was paying over and over, and the fixes.

Item 8 of the open list was "sweep cost: the tool memo covers the tool, not the
work around it". Profiled rather than guessed, the work around it was three
things, and two of them were defects rather than prices:

1. **A .tar.gz read out of order.** The RiPP corpus opened its reference archive
   once per record and the architecture corpus extracted members in name order.
   A gzip stream cannot seek backwards, so each out-of-order read re-inflates
   from the start: 5.3s of a 7.6s corpus load, on every one of 72 runs.
2. **A 1.57 GB Pfam file hashed per run** to report a digest the campaign's own
   `reference/` already pins -- 4.5s a run, six minutes a sweep.
3. **A tool asked its version twice a run**, where that tool is a Python package
   that imports itself to answer: `antismash --version` 2.12s, `bigscape
   --version` 2.38s, against 10ms or less for the other six tools.

These tests hold each fix in place. The sweep times are in docs; what is
asserted here is the behaviour that produced them.
"""

from __future__ import annotations

import gzip
import io
import json
import pathlib
import tarfile
import tempfile

import pytest

from npbench_c.tools import cache
from npbench_c.tools.archive import ArchiveError, read_members

ROOT = pathlib.Path(__file__).resolve().parents[2]
CAMPAIGNS = ROOT / "campaigns"


def _archive(path: pathlib.Path, members: dict[str, bytes]) -> pathlib.Path:
    with tarfile.open(path, "w:gz") as tar:
        for name, blob in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    return path


# ------------------------------------------------------------------- the reader


def test_every_matching_member_is_returned_in_archive_order(tmp_path):
    path = _archive(tmp_path / "a.tar.gz", {
        "z.json": b"3", "a.json": b"1", "m.txt": b"skip", "k.json": b"2"})
    out = read_members(path, ".json")
    assert list(out) == ["z.json", "a.json", "k.json"], (
        "the reading order is the archive's; sorting belongs to the caller")
    assert out["a.json"] == b"1"


def test_a_missing_archive_refuses(tmp_path):
    with pytest.raises(ArchiveError, match="not present"):
        read_members(tmp_path / "absent.tar.gz", ".json")


def test_an_archive_without_matching_members_refuses(tmp_path):
    path = _archive(tmp_path / "b.tar.gz", {"only.txt": b"x"})
    with pytest.raises(ArchiveError, match="no .json members"):
        read_members(path, ".json")


def test_a_duplicate_member_name_refuses(tmp_path):
    path = tmp_path / "dupe.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        for blob in (b"first", b"second"):
            info = tarfile.TarInfo("same.json")
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    with pytest.raises(ArchiveError, match="appears twice"):
        read_members(path, ".json")


# ------------------------------------- the corpora that were re-reading archives


@pytest.mark.parametrize("campaign,module,loader", [
    ("ripp-precursor-mibig-4_0-01", "npbench_c.templates.ripp_precursor.core",
     "ripp"),
    ("arch-nrps-mibig-4_0-01", "npbench_c.templates.domain_architecture.core",
     "arch"),
])
def test_a_corpus_opens_each_archive_once(campaign, module, loader, monkeypatch):
    """The regression guard for the defect itself, not for its symptom.

    Counting opens is what distinguishes "reads the archive" from "reads the
    archive once per record": the second is what cost 5.3 seconds a run, and a
    future refactor that reintroduces a per-member reader passes every
    correctness test in the suite while failing this one.
    """
    directory = CAMPAIGNS / campaign
    if not (directory / "task.yaml").is_file():
        pytest.skip(f"{campaign} is not built")

    opens: list[str] = []
    real_open = tarfile.open

    def counting_open(name=None, *args, **kwargs):
        opens.append(str(name))
        return real_open(name, *args, **kwargs)

    monkeypatch.setattr(tarfile, "open", counting_open)

    import importlib
    import yaml
    params = yaml.safe_load((directory / "task.yaml").read_text())["template_params"]
    core = importlib.import_module(module)
    if loader == "ripp":
        core.Corpus.load(directory, params["annotations_input"],
                         params["reference_records_input"])
    else:
        rules = core.Rules.load(directory)
        core.read_curated(directory / params["annotations_input"], rules)

    counts: dict[str, int] = {}
    for name in opens:
        counts[name] = counts.get(name, 0) + 1
    assert counts and all(n == 1 for n in counts.values()), counts


# --------------------------------------------------------- the declared digest


@pytest.mark.parametrize("campaign,stub_module", [
    ("ripp-precursor-mibig-4_0-01", "npbench_c.sweep.stubs.stub_ripp"),
    ("arch-nrps-mibig-4_0-01", "npbench_c.sweep.stubs.stub_arch"),
    ("gcfcutoff-mibig-32-01", "npbench_c.sweep.stubs.stub_gcfcutoff"),
])
def test_the_stubs_declared_digest_is_the_one_gold_measured(campaign, stub_module):
    """The shortcut cannot drift from gold, because gold hashed the file.

    A stub reports the digest the campaign declares instead of hashing 1.57 GB
    per run. That is only sound while the declared digest IS what the oracle
    measured -- the oracle calls Resource.verify(), which hashes the file and
    refuses to build gold on a mismatch -- so the two are asserted equal here.
    """
    directory = CAMPAIGNS / campaign
    gold_path = directory / "gold" / "gold.json"
    if not gold_path.is_file():
        pytest.skip(f"{campaign} gold not generated")
    source = pathlib.Path(ROOT / "src" / pathlib.Path(
        stub_module.replace(".", "/") + ".py"))
    assert "pfam.sha256" in source.read_text(), (
        f"{stub_module} no longer reports the declared digest; if it hashes the "
        "file again, say why in the stub and delete this test")

    declared = None
    for rules in sorted((directory / "reference").glob("*.json")):
        document = json.loads(rules.read_text())
        for block in (document.get("external_resources") or {}).values():
            if "pfam" in str(block.get("relative_path", "")).lower():
                declared = block["sha256"]
    assert declared, f"{campaign} declares no Pfam resource"
    assert json.loads(gold_path.read_text())["environment"]["pfam_sha256"] == declared


# ----------------------------------------------------------- the version memo


def test_a_version_is_queried_once_per_process(monkeypatch):
    monkeypatch.setattr(cache, "_VERSIONS", {})
    monkeypatch.delenv(cache.ENV_VAR, raising=False)
    calls = []

    def query() -> str:
        calls.append(1)
        return "TOOL 1.2.3"

    assert cache.memoize_version(["tool", "--version"], query) == "TOOL 1.2.3"
    assert cache.memoize_version(["tool", "--version"], query) == "TOOL 1.2.3"
    assert len(calls) == 1, "a tool does not change version inside one process"


def test_a_different_key_is_a_different_version(monkeypatch):
    monkeypatch.setattr(cache, "_VERSIONS", {})
    monkeypatch.delenv(cache.ENV_VAR, raising=False)
    assert cache.memoize_version(["a"], lambda: "A") == "A"
    assert cache.memoize_version(["b"], lambda: "B") == "B"


def test_the_version_memo_survives_a_new_process_only_on_disk(monkeypatch, tmp_path):
    """The on-disk half, which is what a sweep's 72 runs actually use."""
    monkeypatch.setattr(cache, "_VERSIONS", {})
    monkeypatch.setenv(cache.ENV_VAR, str(tmp_path))
    assert cache.memoize_version(["tool"], lambda: "TOOL 9") == "TOOL 9"

    # A fresh process: the per-process dict is empty, the directory is not.
    monkeypatch.setattr(cache, "_VERSIONS", {})

    def refuse() -> str:
        raise AssertionError("the disk memo should have answered")

    assert cache.memoize_version(["tool"], refuse) == "TOOL 9"


def test_without_a_memo_directory_a_new_process_queries_again(monkeypatch):
    monkeypatch.setattr(cache, "_VERSIONS", {})
    monkeypatch.delenv(cache.ENV_VAR, raising=False)
    assert cache.memoize_version(["tool"], lambda: "TOOL 9") == "TOOL 9"
    monkeypatch.setattr(cache, "_VERSIONS", {})
    assert cache.memoize_version(["tool"], lambda: "TOOL 10") == "TOOL 10"


def test_an_empty_version_refuses(monkeypatch):
    monkeypatch.setattr(cache, "_VERSIONS", {})
    monkeypatch.delenv(cache.ENV_VAR, raising=False)
    with pytest.raises(cache.CacheError, match="version query returned"):
        cache.memoize_version(["tool"], lambda: "  ")


def test_only_the_two_slow_tools_are_memoised():
    """Measured, and recorded as a boundary: antiSMASH 2.12s and BiG-SCAPE 2.38s
    import a Python package to answer, the other six answer in 10ms. A memo
    around a 10ms call is a mechanism nobody needs."""
    templates = ROOT / "src" / "npbench_c" / "templates"
    memoised = sorted(
        path.parent.name for path in templates.glob("*/core.py")
        if "memoize_version" in path.read_text())
    assert memoised == ["bgc_detection", "gcf_cutoff"], memoised
