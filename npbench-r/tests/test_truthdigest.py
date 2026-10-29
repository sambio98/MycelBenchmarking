"""The committed fingerprint of an uncommitted truth/.

`tasks/*/truth/` is gitignored, so this digest is the only committed record of
what the answer key holds. These tests pin the two behaviours that give it any
value: it must notice a changed answer key, and it must stay silent for a
checkout whose data has never been fetched.
"""

import json
import shutil
from pathlib import Path

import pytest

from npbench.build import truthdigest
from npbench.build.validate import validate_task

FIXTURE = Path(__file__).parent / "fixtures" / "minimal_task"


@pytest.fixture
def bundle(tmp_path):
    target = tmp_path / "minimal_task"
    shutil.copytree(FIXTURE, target)
    return target


def test_a_fresh_digest_agrees_with_the_truth_it_measured(bundle):
    truthdigest.write(bundle)
    assert truthdigest.compare(bundle) == []


def test_it_records_a_row_count_and_a_hash_per_artifact(bundle):
    entries = truthdigest.write(bundle)
    assert entries
    for entry in entries:
        assert len(entry.sha256) == 64
        assert entry.rows >= 1          # data rows, header excluded
    stored = json.loads((bundle / truthdigest.DIGEST_FILE).read_text())
    assert {e["path"] for e in stored} == {e.path for e in entries}


def test_a_changed_answer_key_is_caught(bundle):
    truthdigest.write(bundle)
    target = next((bundle / "truth").iterdir())
    with target.open("a") as handle:
        handle.write("BGC9999,nrps,contig_1,1,100,100\n")

    problems = truthdigest.compare(bundle)
    assert len(problems) == 1
    assert target.name in problems[0]
    assert "invalidates every score" in problems[0]


def test_a_changed_answer_key_fails_the_validator(bundle):
    truthdigest.write(bundle)
    target = next((bundle / "truth").iterdir())
    with target.open("a") as handle:
        handle.write("BGC9999,nrps,contig_1,1,100,100\n")

    errors = [f for f in validate_task(bundle) if f.level == "error"]
    assert any(truthdigest.DIGEST_FILE in f.where for f in errors)


def test_a_missing_artifact_is_caught(bundle):
    truthdigest.write(bundle)
    next((bundle / "truth").iterdir()).unlink()
    problems = truthdigest.compare(bundle)
    assert problems and "missing from truth/" in problems[0]


def test_an_extra_artifact_is_caught(bundle):
    truthdigest.write(bundle)
    (bundle / "truth" / "surprise.csv").write_text("a\n1\n")
    problems = truthdigest.compare(bundle)
    assert problems and "not in the digest" in problems[0]


def test_an_unfetched_checkout_is_not_a_failure(bundle):
    """A clone that has not run `make truth` has no truth to check.

    truth/ is generated, so its absence is the normal state of a fresh
    checkout. Reporting it here would make `npbench validate` fail for
    everyone who has not spent an hour staging data first.
    """
    truthdigest.write(bundle)
    shutil.rmtree(bundle / "truth")
    assert truthdigest.compare(bundle) == []
    assert not [f for f in validate_task(bundle)
                if f.level == "error" and truthdigest.DIGEST_FILE in f.where]


def test_a_bundle_with_no_digest_is_not_a_failure(bundle):
    """The digest is additive: a task that has not got one still validates."""
    assert not (bundle / truthdigest.DIGEST_FILE).exists()
    assert truthdigest.compare(bundle) == []
