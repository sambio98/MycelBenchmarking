"""EC domain audit: the first campaign that runs a tool, and the first that reads
a resource it is not allowed to ship.

Three things are specific to this campaign and each has a test that would fail if
it drifted. The declared panel must partition the set, or the architecture census
is ambiguous. The canonical-absent sequences must not be fragments, or the whole
audit is an artefact of sequencing completeness rather than a statement about
annotation. And the domtblout columns must be read the right way round, which is
the one parsing mistake here that produces a plausible-looking wrong answer
instead of an error.
"""

from __future__ import annotations

import json
import os
import pathlib

import pytest
import yaml

from npbench_c.grading.grade import grade
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.sweep.runner import build_sandbox
from npbench_c.templates.ec_domain_audit import analyse as eca
from npbench_c.templates.ec_domain_audit.build import _assert_has_content
from npbench_c.templates.ec_domain_audit.core import (
    SECTIONS,
    Audit,
    AuditError,
    Resource,
    Rules,
    tool_prefix,
)
from npbench_c.templates.ec_domain_audit.measure import measure
from npbench_c.tools.resources import ResourceError

CAMPAIGN = (pathlib.Path(__file__).resolve().parents[2]
            / "campaigns" / "ecaudit-fmn-dh-1_1_3_15-01")
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")

image_present = pytest.mark.skipif(
    not tool_prefix().is_dir(),
    reason=f"pinned-tool image not present at {tool_prefix()}; this campaign "
           "cannot be solved without it, which is what min_tool_surface says")


def _gold() -> dict:
    return json.loads((CAMPAIGN / "gold" / "gold.json").read_text())


def _task() -> dict:
    return yaml.safe_load((CAMPAIGN / "task.yaml").read_text())


@pytest.fixture(scope="module")
def audit():
    params = eca.load_params(CAMPAIGN)
    return Audit.load(CAMPAIGN, params["sequence_input"])


@pytest.fixture(scope="module")
def domains(audit):
    """One HMMER pass for the whole module: 13 models over 7,675 sequences."""
    return audit.domains()


# ---------------------------------------------------------------- the gates


def test_oracle_scores_exactly_one():
    m = json.loads((CAMPAIGN / "oracle_submission" / "measurement.json").read_text())
    s = json.loads((CAMPAIGN / "grading.json").read_text())
    r = grade(m, _gold(), s)
    assert r["score"] == 1.0 and r["depth"] == 4 and r["monotonic"] is True
    assert r["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR


def test_all_rungs_are_g1():
    """Gold is the tool's output computed in-container. Nothing here is curated,
    and nothing here is the literature's number."""
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    assert [r["gold_tier"] for r in spec["rungs"]] == ["G1"] * 4


def test_the_tool_surface_is_declared_and_is_not_s0():
    """The first campaign no S0 system can answer. Declaring it is what makes the
    S0 column a finding about tooling rather than a confound."""
    task = _task()
    assert task["min_tool_surface"] == "S1"
    assert "S0 cannot answer" in task["tool_surface_note"]
    assert task["snapshots_required"]["image"] == ["hmmer@3.4", "pfam@35.0"]


# ------------------------------------------------------- the declared panel


@image_present
def test_the_panel_partitions_the_set(audit, domains):
    """Every canonical-present sequence carries the canonical model and nothing
    else from the panel.

    If a panel model overlapped FMN_dh the architecture census would stop being a
    partition and the dominant-architecture claim would depend on which of two
    overlapping families scored higher. Measured, the 6,329 canonical-present
    sequences are all exactly 'FMN_dh'.
    """
    canonical = audit.rules.canonical_model
    present = [a for a, m in domains.items() if canonical in m]
    assert present
    assert all(domains[a] == frozenset({canonical}) for a in present)
    census = audit.architecture_census(domains)
    assert census[canonical] == len(present)


@image_present
def test_every_declared_model_exists_in_the_pinned_release(audit, domains):
    """The build refuses to run if a panel model is missing from the library, so
    the panel and the pinned Pfam release cannot drift apart silently."""
    found = {m for models in domains.values() for m in models}
    assert found <= set(audit.rules.panel)
    # Not every model has to hit, but a panel where most never fire would be a
    # vocabulary chosen for looks.
    assert len(found) >= len(audit.rules.panel) // 2


def test_the_panel_is_published_not_hinted():
    """A declared vocabulary in reference/, reachable by the agent. Making the
    agent discover the families instead would hide the answer behind a full
    Pfam-A scan per run."""
    rules = json.loads((CAMPAIGN / "reference" / "audit_rules.json").read_text())
    assert len(rules["audit"]["model_panel"]) == _gold()["model_panel_size"]
    assert "Declared, not discovered" in rules["audit"]["model_panel_note"]


# --------------------------------------------------------------- the control


def test_the_reviewed_section_is_a_clean_control():
    """29 of 29 reviewed entries carry the canonical domain. Without that the
    unreviewed rate would be uninterpretable: a low rate could mean the model is
    wrong rather than the annotations are."""
    gold = _gold()
    present = gold["canonical_present_by_section"]
    absent = gold["canonical_absent_by_section"]
    assert absent["reviewed"] == 0
    assert present["reviewed"] == gold["sequences_by_section"]["reviewed"]
    assert absent["unreviewed"] > 0
    assert set(present) == set(SECTIONS.values())


def test_the_measured_rate_is_not_the_published_one():
    """Rembeza & Engqvist found 78% in BRENDA 2017.1. This campaign measures
    UniProt 2026_03 and gets 17.6%, and grading the paper's figure would be a G3
    recall test rather than a computation."""
    gold = _gold()
    rate = (gold["canonical_absent_by_section"]["unreviewed"]
            / gold["sequences_by_section"]["unreviewed"])
    assert 0.15 < rate < 0.20
    assert "78%" in _task()["scope_note"]
    excluded = [e["field"] for e in _task()["excluded_from_grading"]]
    assert any("literature" in f for f in excluded)


# ------------------------------------------------- the fragment explanation


@image_present
def test_the_absent_sequences_are_not_fragments(audit, domains):
    """The explanation that would make the campaign an artefact, tested.

    A sequence can miss a 348-column model by being a fragment. If that drove the
    result, 'lacks the domain' would be a statement about sequencing. Measured,
    the canonical-absent side is LONGER than the canonical-present side and has a
    smaller short-sequence fraction, so the absence is not about length.
    """
    profile = audit.length_profile(domains)
    absent, present = profile["canonical_absent"], profile["canonical_present"]
    assert absent["median"] > present["median"]
    for band in ("300", "348"):
        absent_fraction = absent["shorter_than"][band] / absent["count"]
        present_fraction = present["shorter_than"][band] / present["count"]
        assert absent_fraction < present_fraction
    assert profile == _gold()["length_profile"]


def test_the_profile_is_graded_so_an_auditor_sees_it():
    """A refutation nobody can check is an assertion. The profile is a graded
    component, and it is the first thing in notes_for_audit."""
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    graded = {c["name"] for r in spec["rungs"] for c in r["components"]}
    assert "length_profile" in graded
    assert "fragment explanation" in _task()["notes_for_audit"][0]


# ------------------------------------------------------ what they really are


@image_present
def test_the_absent_set_is_a_coherent_family_not_noise(audit, domains):
    """The finding that decides whether 'absent' may be read as 'misannotated'.

    967 of the 1,346 canonical-absent sequences carry exactly
    FAD-oxidase_C + FAD_binding_4 -- the FAD-linked glycolate oxidase subunit
    architecture. They are not junk: they are a different enzyme for the same
    substrate, which is why the campaign grades counts under declared policies
    instead of asserting a misannotation rate.
    """
    absent = audit.absent(domains)
    census = audit.architecture_census(domains, absent)
    dominant, count = next(iter(census.items()))
    assert dominant == "FAD-oxidase_C + FAD_binding_4"
    assert count / len(absent) > 0.5
    assert len(census) >= 2
    assert dominant == _gold()["dominant_absent_architecture"]


@image_present
def test_the_domtblout_columns_are_read_the_right_way_round(audit, domains):
    """The one parsing mistake here that does not fail, it just lies.

    In hmmsearch --domtblout column 1 is the target (a sequence) and column 4 is
    the query (a model). Reading them the other way round reports model names as
    sequence accessions and sequence accessions as models -- no error, a
    plausible-looking table, every number wrong. The check is that what came back
    is keyed by the sequences we shipped and valued by the models we declared.
    """
    assert set(domains) == set(audit.sequences.records)
    seen = {m for models in domains.values() for m in models}
    assert seen <= set(audit.rules.panel)
    assert not (seen & set(audit.sequences.records))


# ------------------------------------------------------------ the two axes


def test_the_policy_axis_moves_and_the_strict_policy_is_the_absence_count():
    """R4's load-bearing axis. The three policies give 1,346 / 319 / 22, and the
    strict one must reproduce the canonical-absent count exactly or the two are
    measuring different things."""
    gold = _gold()
    policies = gold["misannotated_by_policy"]
    assert policies["p0_canonical_only"] == gold["canonical_absent_total"]
    assert (policies["p0_canonical_only"] > policies["p1_accept_fad_oxidase"]
            > policies["p2_accept_any_flavin"])
    assert len(set(policies.values())) == 3


def test_the_threshold_axis_is_inert_and_that_is_the_result():
    """--cut_ga, --cut_nc and --cut_tc give an identical answer and the E-value
    cutoffs move it by under 3%. A counterfactual that does not move is normally a
    defect; here it is a robustness measurement, and it is admissible only because
    the policy axis in the same rung moves by a factor of sixty."""
    gold = _gold()
    absent = gold["threshold_absent"]
    curated = {absent["--cut_ga"], absent["--cut_nc"], absent["--cut_tc"]}
    assert len(curated) == 1
    spread = max(absent.values()) - min(absent.values())
    assert spread / gold["canonical_absent_total"] < 0.03
    policies = gold["misannotated_by_policy"]
    assert (max(policies.values()) - min(policies.values())) > spread * 10
    assert "robustness result, not a defect" in _task()["notes_for_audit"][2]


def test_the_build_refuses_gold_whose_rungs_say_nothing():
    """Every refusal is a defect this campaign could actually have."""
    gold = _gold()

    ok = dict(gold)
    _assert_has_content(ok)

    thin = dict(gold, canonical_absent_total=4,
                misannotated_by_policy={**gold["misannotated_by_policy"],
                                        "p0_canonical_only": 4})
    with pytest.raises(AuditError, match="anecdote"):
        _assert_has_content(thin)

    no_control = dict(gold, canonical_present_by_section={
        **gold["canonical_present_by_section"], "reviewed": 0})
    with pytest.raises(AuditError, match="control"):
        _assert_has_content(no_control)

    flat = dict(gold, misannotated_by_policy={
        k: gold["canonical_absent_total"] for k in gold["misannotated_by_policy"]})
    with pytest.raises(AuditError, match="does not move"):
        _assert_has_content(flat)

    one_architecture = dict(gold, absent_architecture_census={"FMN_dh": 1346})
    with pytest.raises(AuditError, match="single number"):
        _assert_has_content(one_architecture)

    mismatched = dict(gold, misannotated_by_policy={
        **gold["misannotated_by_policy"], "p0_canonical_only": 7})
    with pytest.raises(AuditError, match="reproduce the absence count"):
        _assert_has_content(mismatched)

    # The baseline cutoff is R4's identity control: it re-runs the audit under the
    # same flag that produced the headline and must land on the same number.
    drifted = dict(gold, threshold_absent={**gold["threshold_absent"],
                                           gold["threshold"]: 999})
    with pytest.raises(AuditError, match="identity control"):
        _assert_has_content(drifted)


# -------------------------------------------------- the external resource


def test_an_external_resource_is_declared_with_a_fingerprint():
    """The campaign reads a 1.5 GB library it does not ship. What replaces the
    bytes is a path, a provider, a version and a hash -- otherwise the sandbox
    looks complete while the answer depends on something nobody pinned."""
    rules = Rules.load(CAMPAIGN)
    resource = rules.resources["pfam_a_hmm"]
    assert resource.sha256 == _gold()["pfam_sha256"]
    assert resource.version == "Pfam 35.0"
    declared = _task()["constraints"]["external_resources"]
    assert [e["name"] for e in declared] == ["pfam_a_hmm"]
    assert all(e.get("reason") for e in declared)


@image_present
def test_a_changed_resource_is_caught_rather_than_silently_used():
    """A domain verdict against a different Pfam release is not this campaign's
    gold, so a fingerprint mismatch has to be an error and not a warning."""
    rules = Rules.load(CAMPAIGN)
    real = rules.resources["pfam_a_hmm"]
    assert real.verify()["matches"] is True

    wrong = Resource(name=real.name, relative_path=real.relative_path,
                     provider=real.provider, version=real.version,
                     sha256="f" * 64, note=real.note)
    assert wrong.verify()["matches"] is False

    missing = Resource(name="absent", relative_path="no/such/file.hmm",
                       provider=real.provider, version=real.version,
                       sha256=real.sha256, note="test")
    # The declaration itself now lives in npbench_c.tools.resources, shared by
    # every campaign that reads a resource from the image, so the refusal is a
    # ResourceError rather than this template's own error type.
    with pytest.raises(ResourceError, match="reads it from the pinned image"):
        missing.resolve()


def test_the_library_is_not_redistributed():
    """Not shipped, and the notice says why. The reason here is size, not licence
    -- Pfam is CC0 -- and recording which it is matters for the next S2 campaign,
    where it will be licence."""
    shipped = {p.name for p in (CAMPAIGN / "inputs").iterdir()}
    assert not any(n.endswith(".hmm") for n in shipped)
    rules = json.loads((CAMPAIGN / "reference" / "audit_rules.json").read_text())
    note = rules["external_resources"]["pfam_a_hmm"]["note"]
    assert "size, not licence" in note
    assert _task()["license_notice"]["redistribution"].count("NOT redistributed")


# -------------------------------------------------------- the environment


def test_a_tool_call_puts_the_prefix_on_path_last():
    """The lesson the invariance suite paid for twice. The prefix must be ON the
    subprocess PATH, because a tool resolves its helpers there; and it must be
    LAST, because the image carries its own Python and the grading core is
    stdlib-only on the interpreter this code runs under."""
    from npbench_c.templates.ec_domain_audit.core import _run

    done = _run(["python3", "-c", "import os; print(os.environ['PATH'])"])
    path = done.stdout.strip().split(os.pathsep)
    assert path[-1] == str(tool_prefix())
    assert path.count(str(tool_prefix())) == 1


@image_present
def test_gold_records_what_produced_it(audit):
    """A tool verdict is only as pinned as the tool and the library behind it, so
    the three identity fields are graded as one composite: a system that pinned
    nothing loses the claim rather than one component of it."""
    gold = _gold()
    assert gold["hmmer_version"].startswith("HMMER 3.4")
    assert gold["environment"] == {"hmmer_version": gold["hmmer_version"],
                                   "pfam_version": gold["pfam_version"],
                                   "pfam_sha256": gold["pfam_sha256"]}
    spec = json.loads((CAMPAIGN / "grading.json").read_text())
    r1 = next(r for r in spec["rungs"] if r["rung_id"] == "r1")
    assert "environment" in {c["name"] for c in r1["components"]}


# ------------------------------------------------------------- the measurer


def test_the_measurer_never_reads_gold(tmp_path):
    """Measured rather than asserted: the measurer is run against a campaign copy
    with no gold/ at all, and must produce the same measurement.

    Every vocabulary it needs -- the sections, the panel, the policies, the
    threshold variants, the length bands -- comes from the campaign's own declared
    rules, so a submission cannot be scored against something it was not told.
    """
    import shutil

    stripped = tmp_path / "campaign"
    shutil.copytree(CAMPAIGN, stripped,
                    ignore=shutil.ignore_patterns("gold", "oracle",
                                                  "oracle_submission"))
    assert not (stripped / "gold").exists()
    submission = CAMPAIGN / "oracle_submission"
    assert measure(stripped, submission) == measure(CAMPAIGN, submission)


def test_the_measurer_survives_a_missing_or_broken_report(tmp_path):
    out = measure(CAMPAIGN, tmp_path)
    assert out == {"present": {"report": 0}}
    (tmp_path / "domain_audit.json").write_text("{not json")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}
    (tmp_path / "domain_audit.json").write_text("[]")
    assert measure(CAMPAIGN, tmp_path) == {"present": {"report": 0}}


def test_the_measurer_drops_zeros_from_the_open_census(tmp_path):
    """The architecture census has open-ended keys, so a submission that spells
    out an architecture with no members must compare equal to gold, which omits
    it."""
    gold = _gold()
    report = json.loads(
        (CAMPAIGN / "oracle_submission" / "domain_audit.json").read_text())
    report["architecture_census"] = {**report["architecture_census"],
                                     "Fer4_17 + Pyr_redox": 0}
    (tmp_path / "domain_audit.json").write_text(json.dumps(report))
    out = measure(CAMPAIGN, tmp_path)
    assert out["architecture_census"] == gold["architecture_census"]


# ------------------------------------------------------------- gold isolation


def test_the_sandbox_hands_over_no_answers(tmp_path):
    sandbox = build_sandbox(CAMPAIGN, tmp_path / "s")
    listed = {p.name for p in sandbox.iterdir()}
    assert listed == {"inputs", "reference", "task.yaml", "submission"}
    blob = json.dumps({p.name: p.read_text() for p in
                       (sandbox / "reference").iterdir()})
    gold = _gold()
    for key in ("canonical_absent_total", "sequences_total"):
        assert str(gold[key]) not in blob
    assert gold["dominant_absent_architecture"] not in blob


def test_the_warm_bundles_withhold_their_own_rung(tmp_path):
    """R3's bundle may hand over R2's presence census; it must not hand over the
    architecture census R3 grades."""
    gold = _gold()
    for rung, forbidden in (("r2", "canonical_absent_total"),
                            ("r3", "dominant_absent_architecture"),
                            ("r4", "misannotated_by_policy")):
        bundle = CAMPAIGN / "gold" / "warm" / rung
        blob = json.dumps({p.name: json.loads(p.read_text())
                           for p in bundle.iterdir()})
        assert json.dumps(gold[forbidden]) not in blob


# --------------------------------------------------------- the tool memo


def test_the_tool_memo_is_off_unless_asked_for(audit, monkeypatch):
    """No environment variable, no caching. A build, a test, a real agent run and
    a human invocation all compute from scratch."""
    from npbench_c.tools import cache

    monkeypatch.delenv(cache.ENV_VAR, raising=False)
    assert cache.enabled() is False
    assert cache.memoize_json(["x"], lambda: 1) == 1


@image_present
def test_the_memo_returns_what_the_tool_returned(audit, domains, tmp_path,
                                                 monkeypatch):
    """A memo that changes the answer is worse than no memo."""
    from npbench_c.tools import cache

    monkeypatch.setenv(cache.ENV_VAR, str(tmp_path / "memo"))
    first = audit.domains()
    assert first == domains
    assert list((tmp_path / "memo").glob("*.json"))
    # Second call must not re-run the tool, and must agree.
    second = audit.domains()
    assert second == domains


@image_present
def test_the_memo_key_separates_thresholds(audit, tmp_path, monkeypatch):
    """A key that dropped the flags would hand back an answer from another run.
    The two cutoffs here differ by 25 sequences, so a collision would show."""
    from npbench_c.tools import cache

    monkeypatch.setenv(cache.ENV_VAR, str(tmp_path / "memo"))
    strict = audit.domains("-E 1e-10")
    loose = audit.domains("--cut_ga")
    assert len(audit.absent(strict)) != len(audit.absent(loose))
    assert len(list((tmp_path / "memo").glob("*.json"))) == 2


def test_a_memo_with_a_hole_in_its_key_is_refused():
    """Every key part must be a non-empty string: an unresolved tool version or a
    missing fingerprint would silently widen what the memo matches."""
    from npbench_c.tools.cache import CacheError, memoize_json

    for bad in ([], [""], ["ok", None], ["ok", 3]):
        with pytest.raises(CacheError, match="key part"):
            memoize_json(bad, lambda: 1)


def test_the_runner_refuses_a_memo_for_a_real_system_or_an_ablation(tmp_path):
    """The restriction is the whole point of the mechanism: a sweep over real
    systems is partly a measurement of whether they can drive the tool, and a
    no-tool ablation's claim is that it did not compute."""
    from npbench_c.sweep.runner import SystemSpec, ToolCacheMisuse, run_once
    from npbench_c.sweep.sweep import _tool_cache_for

    real = SystemSpec(name="agent", command=["true"], is_stub=False)
    ablation = SystemSpec(name="stub-noncompute", command=["true"],
                          mode="no_tool", is_stub=True)
    stub = SystemSpec(name="stub-reader", command=["true"], is_stub=True)

    assert _tool_cache_for(real, tmp_path) is None
    assert _tool_cache_for(ablation, tmp_path) is None
    assert _tool_cache_for(stub, tmp_path) == tmp_path / ".toolcache"

    for system in (real, ablation):
        with pytest.raises(ToolCacheMisuse):
            run_once(CAMPAIGN, system, 0, tmp_path / "sb", ["true"],
                     tool_cache=tmp_path / ".toolcache")


def test_the_sweep_report_says_the_wall_clock_is_not_a_timing():
    """A cached sweep's wall clock would otherwise read as a cost measurement of
    the campaign, which it is not."""
    report = json.loads((CAMPAIGN / "internal_sweep.json").read_text())
    meta = report["_meta"]
    assert meta["tool_cache_systems"]
    assert "stub-noncompute" not in meta["tool_cache_systems"]
    assert "NOT a timing measurement" in meta["tool_cache_note"]
