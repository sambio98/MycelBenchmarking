"""The release must not carry the answers, and the exclusions must be earning it.

`build_sandbox` assembles what one system sees for one campaign; this assembles
what leaves the building. Same defect at a different scale: a release that ships
`gold/` hands over every answer at once, and unlike a sandbox leak it cannot be
taken back.

Three of the exclusions are less obvious than `gold/`, and two tests below exist
to prove they are load-bearing rather than tidy. `docs/` holds this project's
build log, and the rule for a measurement that motivated a design decision is
that it belongs there, "which the sandbox never sees" -- true of a sandbox and
false of a release that shipped it. `src/npbench_c/templates/` holds each
template's `analyse.py` and `build.py`, which ARE the oracle. A test asserts that
both really do contain graded values, so that if someone ever argues the
exclusion is paranoid, the measurement answers.
"""

from __future__ import annotations

import json
import pathlib

import pytest
import yaml

from npbench_c.readiness import cold_sandbox
from npbench_c.readiness.gate import _resolve_gold
from npbench_c.release.package import (
    ATTRIBUTION,
    EXCLUDED_TREES,
    MANIFEST,
    README,
    LICENCE_KEYS,
    ReleaseError,
    _attribution_for,
    _licence_subtrees,
    build_release,
    campaigns,
)
from npbench_c.sweep.runner import (
    NEVER_IN_SANDBOX,
    PUBLIC_CAMPAIGN_PATHS,
    SECRET_CAMPAIGN_PATHS,
    GoldLeak,
)

ROOT = pathlib.Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def release(tmp_path_factory) -> pathlib.Path:
    """One release of the whole repo, built once for the module."""
    out = tmp_path_factory.mktemp("release") / "npbench-c"
    build_release(ROOT, out)
    return out


# ----------------------------------------------------------- what is withheld


def test_no_withheld_path_survives(release):
    for path in sorted(release.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(release)
        assert rel.name not in NEVER_IN_SANDBOX, rel
        assert rel.parts[0] not in EXCLUDED_TREES, rel
        if rel.parts[0] == "campaigns" and len(rel.parts) > 2:
            assert rel.parts[2] not in SECRET_CAMPAIGN_PATHS, rel


def test_no_graded_answer_is_stated_in_the_release(release):
    """The same scan the gate runs as check 22, pointed at the release: gold and
    the grading spec come from the build repo, the prose from what would ship."""
    for campaign in campaigns(ROOT):
        if not (campaign / "gold" / "gold.json").is_file():
            continue
        packaged = release / "campaigns" / campaign.name
        hits = cold_sandbox.findings(campaign, packaged, _resolve_gold)
        assert not hits, (campaign.name, hits)


def test_the_docs_exclusion_is_load_bearing():
    """Not tidiness. `docs/` states graded values on purpose, because the rule
    for a measurement behind a design decision is that it lives there rather
    than in the sandbox. A release that shipped it would publish the answers."""
    assert "docs" in EXCLUDED_TREES
    blob = "\n".join(p.read_text(errors="replace")
                     for p in sorted((ROOT / "docs").glob("*.md")))
    stated = []
    for campaign in campaigns(ROOT):
        gold_path = campaign / "gold" / "gold.json"
        if not gold_path.is_file():
            continue
        gold = json.loads(gold_path.read_text())
        for key, value in gold.items():
            if isinstance(value, int) and not isinstance(value, bool) \
                    and value >= 1000 and str(value) in blob:
                stated.append(f"{campaign.name}:{key}={value}")
    assert len(stated) >= 5, (
        "docs/ was expected to state graded values, which is why it is "
        f"excluded; found {stated}")


def test_the_oracle_is_not_shipped(release):
    """Each template's analyse.py and build.py compute gold from the inputs, so
    shipping src/ would let a participant regenerate every answer."""
    assert "src" in EXCLUDED_TREES
    assert not list(release.rglob("*.py"))
    assert not list(release.rglob("analyse.py"))
    assert not list(release.rglob("build.py"))


def test_the_grading_spec_is_withheld(release):
    """Not an answer, but a map of where the marks are: which fields are graded,
    under which scorer, at what tolerance."""
    assert "grading.json" in SECRET_CAMPAIGN_PATHS
    assert not list(release.rglob("grading.json"))


# ------------------------------------------------------------ what is shipped


def test_every_campaign_is_still_solvable_from_the_release(release):
    for campaign in campaigns(ROOT):
        packaged = release / "campaigns" / campaign.name
        assert packaged.is_dir(), campaign.name
        for name in PUBLIC_CAMPAIGN_PATHS:
            if (campaign / name).exists():
                assert (packaged / name).exists(), f"{campaign.name}/{name}"
        assert (packaged / "submission").is_dir()
        assert not any((packaged / "submission").iterdir())


def test_the_public_material_is_byte_identical_to_the_repo(release):
    """A release is a copy, not a rendering. Anything that transformed the task
    description or a pinned table in transit would make the shipped campaign a
    different campaign from the one the gate passed."""
    for campaign in campaigns(ROOT):
        packaged = release / "campaigns" / campaign.name
        for source in sorted(campaign.rglob("*")):
            if not source.is_file():
                continue
            rel = source.relative_to(campaign)
            if rel.parts[0] not in PUBLIC_CAMPAIGN_PATHS:
                continue
            assert (packaged / rel).read_bytes() == source.read_bytes(), rel


def test_the_manifest_covers_every_file_with_a_digest(release):
    manifest = json.loads((release / MANIFEST).read_text())
    on_disk = {str(p.relative_to(release)) for p in release.rglob("*")
               if p.is_file() and p.name != MANIFEST}
    assert set(manifest["files"]) == on_disk
    assert manifest["file_count"] == len(on_disk)
    assert len(manifest["campaigns"]) == len(campaigns(ROOT))
    assert manifest["grader_version"]
    for digest in manifest["files"].values():
        assert len(digest) == 64


def test_the_manifest_does_not_carry_a_gold_digest(release):
    """A digest of gold is not gold, but it lets a system confirm a guess, which
    is the same thing one bit at a time."""
    manifest = json.loads((release / MANIFEST).read_text())
    flat = json.dumps(manifest)
    for campaign in campaigns(ROOT):
        gold = campaign / "gold" / "gold.json"
        if gold.is_file():
            import hashlib
            assert hashlib.sha256(gold.read_bytes()).hexdigest() not in flat


def test_the_readme_says_what_is_absent(release):
    text = (release / README).read_text()
    for word in ("gold", "oracle", "solvable from its own directory"):
        assert word in text


# ------------------------------------------------------------- the attribution


def test_every_campaign_carries_a_licence_statement():
    """A release that cannot attribute its data cannot ship. Several inputs here
    are CC BY, where attribution is a condition rather than a courtesy."""
    for campaign in campaigns(ROOT):
        entry = _attribution_for(campaign)
        assert entry["from_task_yaml"] or entry["from_provenance"], campaign.name


def test_a_campaign_with_no_licence_statement_is_refused(tmp_path):
    """The check that found the one real gap: construct-ecoli-rebh-01 shipped
    UniProt data with no provenance file at all, while its sibling instantiation
    had one."""
    campaign = tmp_path / "campaigns" / "nameless-01"
    (campaign / "inputs").mkdir(parents=True)
    (campaign / "task.yaml").write_text(
        yaml.safe_dump({"campaign_id": "nameless-01", "license_class": "open"}))
    with pytest.raises(ReleaseError, match="states no licence or attribution"):
        _attribution_for(campaign)


def test_a_licence_nested_in_a_provenance_record_counts(tmp_path):
    """The second thing this check got wrong. `pocket-rebh-01` and
    `residues-prna-01` state their licences inside `annotation_source` and
    `structure_sources`, so a top-level scan called two correctly-attributed
    campaigns unattributed."""
    document = {"note": "x", "annotation_source": {"database": "UniProtKB",
                                                   "license": "CC BY 4.0"},
                "structure_sources": [{"database": "RCSB PDB",
                                       "license": "CC0 1.0"}]}
    found = _licence_subtrees(document)
    assert set(found) == {"annotation_source", "structure_sources[0]"}
    assert found["annotation_source"]["license"] == "CC BY 4.0"


def test_a_provenance_source_field_is_not_a_licence():
    """The first false pass: `"source": "PubMed"` inside an evidence record says
    where a claim came from, not what may be redistributed. A rule that accepts
    it accepts almost any JSON."""
    assert "source" not in LICENCE_KEYS
    assert not _licence_subtrees({"evidences": [{"source": "PubMed",
                                                 "id": "17876823"}]})


def test_the_attribution_file_names_every_campaign(release):
    text = (release / ATTRIBUTION).read_text()
    for campaign in campaigns(ROOT):
        assert campaign.name in text, campaign.name


# ---------------------------------------------------------------- the refusals


def test_a_planted_secret_file_is_caught(tmp_path, release):
    """The audit runs over the assembled tree rather than trusting the copy
    loop, so something added afterwards is still caught."""
    from npbench_c.release.package import audit_release
    planted = release / "campaigns" / campaigns(ROOT)[0].name / "gold.json"
    planted.write_text("{}")
    try:
        with pytest.raises(GoldLeak, match="gold.json"):
            audit_release(ROOT, release, campaigns(ROOT))
    finally:
        planted.unlink()


def test_writing_the_release_inside_the_repo_is_refused(tmp_path):
    """The next packaging run would find it and package it again."""
    with pytest.raises(ReleaseError, match="inside the build repo"):
        build_release(ROOT, ROOT / "release")


def test_an_unknown_campaign_is_refused(tmp_path):
    with pytest.raises(ReleaseError, match="no such campaign"):
        build_release(ROOT, tmp_path / "out", ("not-a-campaign",))
