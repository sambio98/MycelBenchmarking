"""Phase 1: the pinned-tool image and its thread-invariance suite.

The rule the provisioning notes set is that every tool runs at 1 and 8 threads on
fixed input and must produce identical output, and that anything failing which
cannot be pinned into determinism is made single-threaded or excluded from
gold-producing paths. These tests enforce the rule and pin what measuring it
found -- including the two cases where the originally declared control was wrong.

The registry and fixture tests need no tools and always run. The invariance tests
need the image; without it they SKIP with a reason, because a determinism suite
that quietly measures nothing is worse than one that fails.
"""

from __future__ import annotations

import json
import pathlib
import shutil

import pytest

from npbench_c.tools import build_fixtures
from npbench_c.tools.invariance import (
    FAIL,
    FIXTURES,
    PASS,
    PROJECTIONS,
    SKIP,
    XFAIL,
    ProjectionError,
    run_tool,
)
from npbench_c.tools.registry import (
    BANNED_BINARIES,
    NOT_IN_IMAGE,
    REPEATS,
    THREAD_COUNTS,
    TOOLS,
    banned_present,
    fingerprint,
    verify,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
LOCK = ROOT / "image" / "environment.lock.json"
REPORT = ROOT / "image" / "thread_invariance.json"
#: Where the image puts the binaries. An override keeps the suite usable from a
#: different prefix without editing code.
import os

PREFIX = os.environ.get("NPBENCH_TOOL_PREFIX", "/opt/npbench-tools/bin")

#: {pfam} resolves to the Pfam library inside the antiSMASH databases, which
#: BiG-SCAPE reuses rather than fetching a second 1.5 GB copy.
ALLOWED_PLACEHOLDERS = {"bin", "fixtures", "setup", "work", "threads", "pfam"}

image_present = pytest.mark.skipif(
    not pathlib.Path(PREFIX).is_dir(),
    reason=f"pinned-tool image not present at {PREFIX}; nothing measured")

#: A tool marked slow takes minutes per pass -- BiG-SCAPE runs hmmscan over Pfam
#: for every cluster in the fixture. The suite still measures it and the verdict
#: is committed in image/thread_invariance.json; what this controls is whether the
#: default test run re-measures it. A skip with a reason is not a green tick.
RUN_SLOW = os.environ.get("NPBENCH_SLOW_TOOLS") == "1"


def skip_if_slow(tool) -> None:
    if tool.slow and not RUN_SLOW:
        pytest.skip(
            f"{tool.name} takes minutes per pass; its measured verdict is in "
            "image/thread_invariance.json. Set NPBENCH_SLOW_TOOLS=1 to re-measure")

#: Each tool is measured once per session. Several tests read the same result and
#: re-running the suite per test would multiply the wall clock for nothing -- the
#: runs inside run_tool are what establish determinism, not repeating run_tool.
_MEASURED: dict[str, dict] = {}


def measure(tool) -> dict:
    if tool.name not in _MEASURED:
        _MEASURED[tool.name] = run_tool(tool, PREFIX)
    return _MEASURED[tool.name]


# ------------------------------------------------------------- the registry

def test_the_suite_runs_at_the_declared_thread_counts():
    """1 and 8, twice each. The rule is specific and the code must not drift."""
    assert THREAD_COUNTS == (1, 8)
    assert REPEATS >= 2


@pytest.mark.parametrize("tool", TOOLS, ids=lambda t: t.name)
def test_every_tool_declares_what_it_is_for_and_how_to_run_it(tool):
    assert tool.pinned_version and tool.needed_by and tool.controls
    assert tool.invocations, "a tool with no invocation is a tool not measured"
    for invocation in tool.invocations:
        assert invocation.artifact
        placeholders = {
            token.split("{", 1)[1].split("}", 1)[0]
            for command in (*invocation.setup, invocation.argv)
            for token in command if "{" in token}
        assert placeholders <= ALLOWED_PLACEHOLDERS, placeholders


@pytest.mark.parametrize("tool", TOOLS, ids=lambda t: t.name)
def test_every_declared_licence_carries_a_reason(tool):
    """A normalisation or a canonicalisation without a reason is a loophole."""
    for rule in tool.normalisations:
        assert rule.pattern and len(rule.reason) > 40, rule.pattern
    if tool.canonicalisation:
        assert tool.canonicalisation.mode == "sort_lines"
        assert len(tool.canonicalisation.reason) > 40


@pytest.mark.parametrize("tool", TOOLS, ids=lambda t: t.name)
def test_a_tool_with_a_thread_flag_is_tested_across_thread_counts(tool):
    if tool.thread_flag:
        assert tool.thread_counts == THREAD_COUNTS
        # The flag has to appear in something the suite actually runs, or the
        # thread count under test never reaches the tool.
        assert any(tool.thread_flag in invocation.argv
                   for invocation in tool.invocations)
    else:
        assert tool.thread_counts == (1,)


def test_the_registry_pins_match_the_lockfile():
    """Two files could drift apart; the image build would then install one thing
    and the suite would verify another."""
    lock = json.loads(LOCK.read_text())
    requested = lock["requested"]
    assert {tool.name for tool in TOOLS} <= set(requested)
    for tool in TOOLS:
        assert requested[tool.name].split("=")[0] == tool.pinned_version
    # Anything pinned that is not a registry tool -- the interpreter -- still has
    # to say why it is pinned.
    for name in requested:
        assert name in lock["requested_rationale"], name


def test_the_lockfile_can_rebuild_the_image_byte_for_byte():
    lock = json.loads(LOCK.read_text())
    assert lock["packages"], "an empty closure is not a lock"
    assert len(lock["packages"]) == lock["package_count"]
    for package in lock["packages"]:
        assert package["url"] and package["sha256"], package["name"]
        assert len(package["sha256"]) == 64


def test_the_absences_are_declared_rather_than_implied():
    """What is missing is a decision with a reason, not a gap for an auditor to
    notice by counting rows."""
    assert {"matchms", "iqtree"} <= set(NOT_IN_IMAGE)
    for reason in NOT_IN_IMAGE.values():
        assert len(reason) > 60
    installed = {tool.name for tool in TOOLS}
    assert not (installed & set(NOT_IN_IMAGE))


def test_a_ban_is_enforced_rather_than_written_down():
    """FastTreeMP is not a deferral, it is a ban -- and it arrives whether the ban
    is remembered or not, because antiSMASH depends on the fasttree package that
    ships it. So it lives in BANNED_BINARIES, which is checked."""
    assert "FastTreeMP" in BANNED_BINARIES
    assert "FastTreeMP" not in NOT_IN_IMAGE, (
        "a ban and a deferral are different things; keeping it in both gives two "
        "sources of truth")
    assert len(BANNED_BINARIES["FastTreeMP"]) > 100


# ------------------------------------------------------------- projections

def test_the_antismash_projection_keeps_regions_and_drops_run_metadata():
    document = {
        "version": "8.0.4", "input_file": "/some/run/path.embl",
        "timestamp": "whenever",
        "records": [{
            "id": "AY312585.1",
            "features": [
                {"type": "region", "location": "[100:2000](+)",
                 "qualifiers": {"product": ["NRPS", "siderophore"]}},
                {"type": "CDS", "location": "[1:50](+)", "qualifiers": {}},
                {"type": "region", "location": "[10:90](+)",
                 "qualifiers": {"product": ["terpene"]}},
            ],
        }],
    }
    projected = PROJECTIONS["antismash_regions"](json.dumps(document).encode())
    assert projected == ("AY312585.1\t[10:90](+)\tterpene\n"
                         "AY312585.1\t[100:2000](+)\tNRPS,siderophore\n"), (
        "regions must come out in coordinate order, not lexical order: a string "
        "sort would put [100:2000] before [10:90], and the projection is what a "
        "campaign reads")
    # The run metadata is gone and the regions are sorted, so a reordered
    # features list is not mistaken for a different detection.
    assert "8.0.4" not in projected and "/some/run/path" not in projected


def test_the_domain_projection_keeps_architecture_and_drops_scores():
    """E-values and bitscores are floats whose last digits are a reduction-order
    artefact. The architecture -- which domains, in which order, over which
    residues -- is the claim, and it is integers and strings."""
    document = {"records": [{"id": "R1", "modules": {
        "antismash.detection.nrps_pks_domains": {"cds_results": {
            "cdsB": {"domain_hmms": [
                {"hit_id": "PKS_KS", "query_start": 10, "query_end": 40,
                 "evalue": 1e-30, "bitscore": 99.5}]},
            "cdsA": {"domain_hmms": [
                {"hit_id": "PP-binding", "query_start": 300, "query_end": 350,
                 "evalue": 2e-10, "bitscore": 30.1},
                {"hit_id": "AMP-binding", "query_start": 44, "query_end": 449,
                 "evalue": 7.3e-72, "bitscore": 234.1}]},
            "cdsC": {"domain_hmms": []},
        }}}}]}
    projected = PROJECTIONS["antismash_nrps_pks_domains"](
        json.dumps(document).encode())
    assert projected == (
        "R1\tcdsA\tAMP-binding:44-449;PP-binding:300-350\n"
        "R1\tcdsB\tPKS_KS:10-40\n")
    assert "234.1" not in projected and "7.3e-72" not in projected
    # A CDS with no domain hit contributes no line rather than an empty one.
    assert "cdsC" not in projected


def test_the_bigscape_projection_keeps_the_partition_and_drops_the_labels():
    """FAM_00001 is a name, not a result. Two runs can agree completely on which
    clusters belong together and disagree on what the groups are called."""
    table = (
        "# other_clustering_c0.3.tsv\n"
        "Record\tGBK\tRecord_Type\tRecord_Number\tCC\tFamily\n"
        "B.gbk_region_1\tB\tregion\t1\t26\tFAM_00002\n"
        "A.gbk_region_1\tA\tregion\t1\t26\tFAM_00002\n"
        "C.gbk_region_1\tC\tregion\t1\t31\tFAM_00009\n"
        "# PKS_clustering_c0.3.tsv\n"
        "Record\tGBK\tRecord_Type\tRecord_Number\tCC\tFamily\n"
        "D.gbk_region_1\tD\tregion\t1\t30\tFAM_00005\n")
    projected = PROJECTIONS["bigscape_gcf_partition"](table.encode())
    assert projected == ("PKS\tD.gbk_region_1\n"
                         "other\tA.gbk_region_1,B.gbk_region_1\n"
                         "other\tC.gbk_region_1\n")

    # The same partition under different labels and a different row order must
    # project identically -- that is the whole point of the projection.
    relabelled = (table.replace("FAM_00002", "FAM_01111")
                       .replace("FAM_00009", "FAM_02222")
                       .replace("\t26\t", "\t99\t"))
    assert PROJECTIONS["bigscape_gcf_partition"](relabelled.encode()) == projected


def test_the_bigscape_projection_refuses_to_find_nothing():
    for table in ("", "# other_clustering_c0.3.tsv\nRecord\tGBK\n"):
        with pytest.raises(ProjectionError):
            PROJECTIONS["bigscape_gcf_partition"](table.encode())


def test_the_sqlalchemy_bound_is_recorded_where_it_will_be_read():
    """BiG-SCAPE 2.0.3 crashes on sqlalchemy 2.1 before reading an input file, and
    its own recipe does not say so. If that fact lives only in a commit message it
    will be lost the next time someone re-solves the environment."""
    from npbench_c.tools.registry import BIGSCAPE

    controls = " ".join(BIGSCAPE.controls).lower()
    assert "sqlalchemy" in controls and "2.1" in controls
    lock = json.loads(LOCK.read_text())
    assert lock["requested"]["sqlalchemy"].startswith("2.0.")
    assert "sqlalchemy" in lock["requested_rationale"]


def test_the_antismash_projection_refuses_to_find_nothing():
    """A projection that silently yields nothing turns the determinism test into a
    tautology: every run would agree on the same emptiness."""
    for document in ({"records": []},
                     {"records": [{"id": "x", "features": [
                         {"type": "CDS", "location": "[1:2]", "qualifiers": {}}]}]}):
        with pytest.raises(ProjectionError):
            PROJECTIONS["antismash_regions"](json.dumps(document).encode())


# ------------------------------------------------------------- the fixtures

def test_the_fixtures_are_committed_and_hashed():
    provenance = json.loads((FIXTURES / "provenance.json").read_text())
    import hashlib

    for name in ("halogenases.faa", "halogenases.afa", "synthetic.fna",
                 "bgc_triplet.embl"):
        path = FIXTURES / name
        assert path.is_file(), name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert provenance[name]["sha256"] == digest, name


def test_the_nucleotide_fixture_regenerates_byte_identically():
    """Pure Python, no download and no toolkit: the declared codon cycle and
    spacer make the contig a function of the protein fixture alone."""
    records = build_fixtures.read_fasta(FIXTURES / "halogenases.faa")
    rebuilt = build_fixtures.build_contig(records)
    committed = "".join(
        line.strip() for line in (FIXTURES / "synthetic.fna").read_text().splitlines()
        if not line.startswith(">"))
    assert rebuilt == committed


def test_the_bigscape_fixture_has_designed_family_structure():
    """32 clusters chosen so the set has real families rather than 32 singletons:
    a clustering fixture made of unrelated BGCs tests almost nothing."""
    provenance = json.loads((FIXTURES / "provenance.json").read_text())
    block = provenance["bigscape_input"]
    files = sorted(p.stem for p in (FIXTURES / "bigscape_input").glob("*.gbk"))
    listed = sorted({a for accs in block["families"].values() for a in accs}
                    | set(block["controls"]))
    assert files == listed
    assert len(block["families"]) == 4 and len(files) == 32
    # The placeholder exclusion is the one judgement in the selection, so it is
    # named rather than left implicit.
    assert "capsular polysaccharide" in block["excluded_placeholder_names"]
    assert "as8b1" in block["url"]


def test_the_protein_fixture_is_a_real_homolog_set():
    records = build_fixtures.read_fasta(FIXTURES / "halogenases.faa")
    assert len(records) >= 20
    headers = " ".join(h for h, _ in records)
    # RebH and PrnA are already in the benchmark, so the fixture introduces no
    # provenance the project does not already carry.
    assert "Q8KHZ8" in headers and "P95480" in headers
    # One fragment, deliberately: a short record is where tie-breaks show.
    lengths = sorted(len(s) for _, s in records)
    assert lengths[0] < lengths[-1] / 2


# --------------------------------------------------------- measured results

@image_present
def test_every_pinned_version_is_installed():
    rows = verify(PREFIX)
    wrong = [(name, pinned, found) for name, pinned, found, ok in rows if not ok]
    assert not wrong, f"pin drift: {wrong}"


@image_present
def test_no_banned_binary_is_present():
    assert banned_present(PREFIX) == []


@image_present
def test_the_antismash_rule_set_is_fingerprinted():
    """A region set is a function of the detection rules that produced it, and the
    rule count has moved every major version: 58 at v5, 88 at v7.1, 103 here. The
    standing instruction never to compare across versions needs a mechanical form,
    and a hash over the rule files is it -- a changed rule set is then detectable
    even when the version string is not what moved."""
    from npbench_c.tools.registry import ANTISMASH

    package_root = pathlib.Path(PREFIX).parent / "lib"
    candidates = sorted(package_root.glob("python3.*/site-packages/antismash"))
    if not candidates:
        pytest.skip("antiSMASH package root not found")
    digest = fingerprint(ANTISMASH, str(candidates[-1]))
    assert digest and len(digest) == 16
    rules = sum(
        1 for relative in ANTISMASH.fingerprint.paths
        for line in (candidates[-1] / relative).read_text().splitlines()
        if line.startswith("RULE "))
    assert rules == 103, (
        f"the detection rule set has {rules} rules, not the 103 this pin was "
        "fingerprinted against; any campaign built on it needs regenerating")


@image_present
@pytest.mark.parametrize("tool", TOOLS, ids=lambda t: t.name)
def test_tool_is_thread_invariant(tool):
    skip_if_slow(tool)
    result = measure(tool)
    if result["status"] == SKIP:
        pytest.skip(result["detail"])
    assert result["status"] == PASS, result["detail"]
    for invocation in result["invocations"]:
        assert invocation["status"] in (PASS, XFAIL), invocation["detail"]
        if invocation["status"] == PASS:
            assert invocation["canonical_identical"], invocation["detail"]


@image_present
def test_diamond_no_reorder_is_still_the_witness_it_is_declared_to_be():
    """The originally declared control had it backwards.

    The provisioning notes listed --no-reorder as a determinism control for
    DIAMOND. Measured, that flag is what BREAKS determinism: with it, 2.2.8 emits
    the same hit set in a different query order on each multithreaded run, and
    without it every run at 1, 4 and 8 threads is byte-identical. The witness
    invocation keeps that fact under test, so a future release that fixes it is
    noticed rather than assumed.
    """
    from npbench_c.tools.registry import DIAMOND

    result = measure(DIAMOND)
    if result["status"] == SKIP:
        pytest.skip(result["detail"])
    by_name = {i["invocation"]: i for i in result["invocations"]}
    assert by_name["blastp_tabular"]["raw_identical"]
    witness = by_name["blastp_tabular_no_reorder"]
    assert witness["status"] == XFAIL, (
        "the --no-reorder witness passed; re-read the pin before trusting it: "
        + witness["detail"])
    assert not witness["repeatable_at_fixed_threads"]


@image_present
def test_mmseqs2_needs_a_canonical_order_and_loses_nothing_to_it():
    """No flag fixes MMseqs2's multithreaded output order, so the order is fixed
    by the invocation contract. The content is untouched: at eight threads the hit
    multiset is identical and only its order moves, which is exactly why sorting
    is sound here and would not be on an artifact with a header."""
    from npbench_c.tools.registry import MMSEQS2

    assert MMSEQS2.canonicalisation is not None
    result = measure(MMSEQS2)
    if result["status"] == SKIP:
        pytest.skip(result["detail"])
    invocation = result["invocations"][0]
    assert invocation["canonical_identical"]
    assert not invocation["normalised_identical"], (
        "MMseqs2 came out order-stable without canonicalisation; if that is real "
        "the contract can be relaxed, but check before relaxing it")


@image_present
def test_a_projected_invocation_never_claims_a_raw_match():
    """antiSMASH's JSON carries the input path and a timestamp, so its raw bytes
    cannot be identical across runs. Reporting 'raw-identical' for a projected
    invocation would overstate the result exactly the way the incidental-match
    rule exists to prevent."""
    from npbench_c.tools.registry import ANTISMASH

    result = measure(ANTISMASH)
    if result["status"] == SKIP:
        pytest.skip(result["detail"])
    for invocation in result["invocations"]:
        assert invocation["projection"], invocation["invocation"]
        assert invocation["compared"].startswith("projection:")
        assert not invocation["raw_identical"]
        assert invocation["projection_identical"], invocation["detail"]


@image_present
def test_antismash_detects_the_three_fixture_classes():
    """The fixture earns its place only if the detection rules it exercises are
    not all of a kind. Three records, three classes, four product names."""
    from npbench_c.tools.registry import ANTISMASH

    result = measure(ANTISMASH)
    if result["status"] == SKIP:
        pytest.skip(result["detail"])
    assert result["status"] == PASS
    assert result["fingerprint"], "the rule set must be fingerprinted"


@image_present
def test_the_committed_report_is_not_stale():
    """image/thread_invariance.json records the measured verdicts. Digests move
    with wall-clock stamps, so only the verdicts are compared."""
    if not REPORT.is_file():
        pytest.skip("no committed invariance report")
    committed = json.loads(REPORT.read_text())
    assert set(committed["tools"]) == {tool.name for tool in TOOLS}
    for tool in TOOLS:
        if tool.slow and not RUN_SLOW:
            continue
        recorded = committed["tools"][tool.name]
        fresh = measure(tool)
        if fresh["status"] == SKIP:
            pytest.skip(fresh["detail"])
        assert fresh["status"] == recorded["status"], tool.name
        fresh_by_name = {i["invocation"]: i["status"] for i in fresh["invocations"]}
        recorded_by_name = {i["invocation"]: i["status"]
                            for i in recorded["invocations"]}
        assert fresh_by_name == recorded_by_name, tool.name
