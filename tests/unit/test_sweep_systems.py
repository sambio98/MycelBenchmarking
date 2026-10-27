"""The two pieces a real sweep needs: the prompt, and a declared system set.

Every sweep on record is stub-based, and the readiness gate keeps its five
agent-run checks PENDING because of it. What was missing was never the runner --
a system has always been a command -- but the instructions a real agent needs and
a validated set of systems to run. These tests hold both, plus the harness defect
that only a stub ladder could hide: `str.format` over a command line.
"""

from __future__ import annotations

import json
import pathlib
import tempfile

import pytest
import yaml

from npbench_c.readiness import cold_sandbox
from npbench_c.readiness.gate import _resolve_gold
from npbench_c.readiness import gate
from npbench_c.sweep import prompt as prompt_module
from npbench_c.sweep import systems as systems_module
from npbench_c.sweep.runner import _format, build_sandbox

ROOT = pathlib.Path(__file__).resolve().parents[2]
CAMPAIGNS = ROOT / "campaigns"


def _campaigns() -> list[pathlib.Path]:
    return [c for c in sorted(CAMPAIGNS.iterdir())
            if (c / "gold" / "gold.json").is_file()
            and (c / "grading.json").is_file()]


@pytest.fixture(scope="module")
def cold(tmp_path_factory) -> dict[str, pathlib.Path]:
    root = tmp_path_factory.mktemp("prompt-cold")
    out = {}
    for campaign in _campaigns():
        sandbox = build_sandbox(campaign, root / campaign.name, None)
        prompt_module.write(sandbox, "cold")
        out[campaign.name] = sandbox
    return out


# ------------------------------------------------------------------ the prompt


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_the_cold_prompt_discloses_nothing(campaign, cold):
    """The prompt is built from the sandbox, so it cannot carry anything the
    system was not already given -- asserted with the gate's own scans rather
    than by reading it."""
    sandbox = cold[campaign.name]
    hits = (cold_sandbox.findings(campaign, sandbox, _resolve_gold)
            + cold_sandbox.single_term_findings(campaign, sandbox, _resolve_gold))
    assert not hits, f"{campaign.name}: the rendered prompt states {sorted(hits)}"


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_the_prompt_names_every_contracted_artefact(campaign, cold):
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    text = (cold[campaign.name] / prompt_module.PROMPT).read_text()
    for artefact in task["output_contract"]:
        assert f"submission/{artefact}" in text, artefact
    assert "## What to produce" in text
    assert "## How you are judged" in text


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_the_prompt_withholds_the_auditor_s_notes(campaign, cold):
    """`notes_for_audit` is written for someone holding gold. It says what to
    recompute and from where, which in front of a system reads as a method."""
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    text = (cold[campaign.name] / prompt_module.PROMPT).read_text()
    for note in task.get("notes_for_audit", []):
        sentence = " ".join(str(note).split())[:60]
        if len(sentence) > 30:
            assert sentence not in " ".join(text.split()), sentence


def test_every_declared_warm_bundle_renders_a_warm_prompt(tmp_path):
    """A warm prompt must say what was handed over, whatever the bundle's file
    type: the first draft looked for `*.json` and raised on the two construct
    campaigns, whose R2 bundle is a FASTA."""
    rendered = 0
    for campaign in _campaigns():
        task = yaml.safe_load((campaign / "task.yaml").read_text())
        for rung in task.get("rungs", []):
            bundle = rung.get("warm_start_bundle")
            if not bundle:
                continue
            sandbox = build_sandbox(
                campaign, tmp_path / campaign.name / rung["rung_id"], bundle)
            text = prompt_module.render(sandbox, f"warm:{rung['rung_id']}")
            assert "## What you have been given" in text
            handed = [p.name for p in sandbox.iterdir()
                      if p.is_file() and p.name not in prompt_module.FURNITURE]
            assert handed, (campaign.name, rung["rung_id"])
            for name in handed:
                assert f"`{name}`" in text, (campaign.name, name)
            rendered += 1
    assert rendered >= 60, rendered


def test_a_warm_prompt_that_hands_over_nothing_is_refused(tmp_path):
    campaign = _campaigns()[0]
    sandbox = build_sandbox(campaign, tmp_path / "cold", None)
    with pytest.raises(prompt_module.PromptError, match="cold run mislabelled"):
        prompt_module.render(sandbox, "warm:r2")


def test_a_sandbox_without_a_task_is_refused(tmp_path):
    (tmp_path / "empty").mkdir()
    with pytest.raises(prompt_module.PromptError, match="no task.yaml"):
        prompt_module.render(tmp_path / "empty")


def test_the_prompt_is_byte_identical_for_every_system(tmp_path):
    """A difficulty comparison between two systems given different instructions
    measures the instructions."""
    campaign = _campaigns()[0]
    first = build_sandbox(campaign, tmp_path / "a", None)
    second = build_sandbox(campaign, tmp_path / "b", None)
    assert (prompt_module.render(first, "cold")
            == prompt_module.render(second, "cold"))


def test_a_nested_contract_group_is_rendered_as_a_list():
    """It was rendered with `str()` at first, which put a Python dict repr in
    front of the system."""
    rows = prompt_module._fields({
        "shells": {"cutoff": "the pinned cutoff", "contacts": "per copy"},
        "flat": "a plain field",
    })
    assert "- **`shells`**" in rows
    assert "  - **`cutoff`** — the pinned cutoff" in rows
    assert "- **`flat`** — a plain field" in rows
    assert not any("{" in row for row in rows)


# --------------------------------------------------------- the command format


def test_a_brace_in_a_command_is_not_a_format_placeholder():
    """The defect a stub ladder could not have found. A stub's argv is written in
    this repository and never contains a brace; a real agent's command plausibly
    contains several, and `printf '{}' > {submission}/report.json` raised before
    it reached a shell."""
    out = _format(["bash", "-lc", "printf '{}' > {submission}/r.json; echo ${HOME}"],
                  submission="/tmp/sub", campaign="/c", sandbox="/s", prompt="/p")
    assert out[2] == "printf '{}' > /tmp/sub/r.json; echo ${HOME}"


def test_every_placeholder_is_substituted():
    out = _format(["{campaign}", "{submission}", "{sandbox}", "{prompt}"],
                  campaign="C", submission="S", sandbox="B", prompt="P")
    assert out == ["C", "S", "B", "P"]


# ------------------------------------------------------------- the system set


def _config(tmp_path, systems) -> pathlib.Path:
    path = tmp_path / "systems.json"
    path.write_text(json.dumps({"systems": systems}))
    return path


SOUND = [
    {"name": "tooled-reference", "reference": True,
     "command": ["agent", "--prompt", "{prompt}", "--out", "{submission}"]},
    {"name": "tooled-generic",
     "command": ["agent", "--plain", "--out", "{submission}"]},
    {"name": "ablation", "mode": "no_tool",
     "command": ["agent", "--out", "{submission}"]},
]


def test_a_sound_config_loads_with_the_campaign_s_own_budget(tmp_path):
    campaign = _campaigns()[0]
    budget = yaml.safe_load((campaign / "task.yaml").read_text())["budgets"]
    specs = systems_module.real_systems(campaign, _config(tmp_path, SOUND))
    assert [s.name for s in specs] == [s["name"] for s in SOUND]
    assert all(not s.is_stub for s in specs)
    expected = int(budget["wall_clock_minutes"] * 60 * systems_module.TIMEOUT_MARGIN)
    assert all(s.timeout_s == expected for s in specs)
    ablation = next(s for s in specs if s.mode == "no_tool")
    assert str(ablation.env["PATH"]).find("npbench-tools") == -1


def test_a_config_cannot_raise_the_campaign_s_budget(tmp_path):
    campaign = _campaigns()[0]
    greedy = [dict(entry, timeout_s=10 ** 6) for entry in SOUND]
    specs = systems_module.real_systems(campaign, _config(tmp_path, greedy))
    budget = yaml.safe_load((campaign / "task.yaml").read_text())["budgets"]
    ceiling = int(budget["wall_clock_minutes"] * 60 * systems_module.TIMEOUT_MARGIN)
    assert all(s.timeout_s == ceiling for s in specs)
    modest = [dict(entry, timeout_s=60) for entry in SOUND]
    assert all(s.timeout_s == 60 for s in
               systems_module.real_systems(campaign, _config(tmp_path, modest)))


@pytest.mark.parametrize("mutate,pattern", [
    (lambda s: [], "declares no systems"),
    (lambda s: [dict(s[0], name=None)] + s[1:], "present and unique"),
    (lambda s: [s[0], dict(s[1], name=s[0]["name"]), s[2]], "present and unique"),
    (lambda s: [dict(s[0], is_stub=True)] + s[1:], "laundered into evidence"),
    (lambda s: [dict(s[0], mode="offline")] + s[1:], "declares mode"),
    (lambda s: [dict(s[0], command=[])] + s[1:], "declares no command"),
    (lambda s: [dict(s[0], command=["agent", "--plain"])] + s[1:],
     "never mentions"),
    (lambda s: s[:2], "exactly one no_tool ablation"),
    (lambda s: s + [dict(s[2], name="second-ablation")],
     "exactly one no_tool ablation"),
    (lambda s: [dict(s[0], reference=False)] + s[1:],
     "exactly one reference system"),
    (lambda s: [s[0], dict(s[1], reference=True), s[2]],
     "exactly one reference system"),
    (lambda s: [dict(s[0], reference=False), dict(s[1], reference=False),
                dict(s[2], reference=True)], "both the reference system"),
    (lambda s: [s[0], s[2]], "every tooled system is the reference"),
    (lambda s: s[:2] + [dict(s[2], env={"PATH": "/opt/npbench-tools/bin:/usr/bin"})],
     "carries the tool prefix"),
])
def test_the_loader_refuses_a_config_that_cannot_measure(tmp_path, mutate, pattern):
    import copy

    config = _config(tmp_path, mutate(copy.deepcopy(SOUND)))
    with pytest.raises(systems_module.SystemConfigError, match=pattern):
        systems_module.load(config)


def test_a_campaign_without_a_budget_cannot_derive_a_timeout(tmp_path):
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "task.yaml").write_text("campaign_id: x\n")
    with pytest.raises(systems_module.SystemConfigError, match="no budgets"):
        systems_module.real_systems(campaign, _config(tmp_path, SOUND))


def test_the_cli_refuses_to_mix_fixtures_with_real_systems(tmp_path, capsys):
    from npbench_c.sweep import cli

    code = cli.main([str(_campaigns()[0]), "--stubs",
                     "--systems", str(_config(tmp_path, SOUND))])
    assert code == 2
    assert "not both" in capsys.readouterr().err


def test_the_cli_says_what_a_sweep_needs_when_given_neither(capsys):
    from npbench_c.sweep import cli

    code = cli.main([str(_campaigns()[0])])
    assert code == 2
    assert "--systems" in capsys.readouterr().err


def test_a_trial_run_cannot_become_the_campaign_s_evidence(tmp_path, monkeypatch):
    """The near-miss this flag exists for.

    `sweep_campaign` records into the campaign by default, which is right for a
    real sweep. The first smoke test of the real-system path used two throwaway
    systems that write `{}` -- and wrote the result into
    `constraint-gc62-ecoli-rebh-01/internal_sweep.json`, marked `stub_based:
    false`. The readiness gate READS that file and promotes a non-stub sweep to a
    measurement, so a smoke test had quietly become a campaign's evidence, with
    three gates failing. It was committed before it was noticed.
    """
    from npbench_c.sweep import cli

    seen = {}

    def fake_sweep(campaign, systems, sandbox_root, repeats, write=True):
        seen["write"] = write
        return {"a_gate": {"pass": True, "detail": "fixture"},
                "_meta": {"campaign_id": "x", "repeats": repeats,
                          "chance_floor": 0.0, "systems": [], "stub_systems": [],
                          "stub_based": False, "run_status_counts": {}},
                "_headline": {"cold_mean_score": None, "ci_low": None,
                              "ci_high": None},
                "_per_rung_clear_rates": {}, "_per_system": {}}

    monkeypatch.setattr(cli, "sweep_campaign", fake_sweep)
    config = _config(tmp_path, SOUND)
    campaign = str(_campaigns()[0])

    cli.main([campaign, "--systems", str(config), "--no-record"])
    assert seen["write"] is False
    cli.main([campaign, "--systems", str(config)])
    assert seen["write"] is True


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_the_prompt_names_every_field_the_contract_declares(campaign, cold):
    """Every graded field, at every depth the contract nests to.

    The contract is what a system is told, so a field missing from the rendering
    is a field no system will write -- and the gate's contract check would not
    see it, because it reads `task.yaml` rather than the prompt.
    """
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    text = (cold[campaign.name] / prompt_module.PROMPT).read_text()
    for artefact, declared in task["output_contract"].items():
        for field in gate._contract_fields(declared):
            assert f"`{field}`" in text, f"{artefact}: {field}"


@pytest.mark.parametrize("campaign", _campaigns(), ids=lambda c: c.name)
def test_the_prompt_separates_conventions_from_fields(campaign, cold):
    """A convention is rendered, and rendered as something other than a field.

    Four campaigns declare a format, a closed vocabulary or a property of a
    FASTA record. Each is part of the contract and none of them is a key to
    write, so the prompt says which it is -- the alternative is a system
    reporting `formula_format` as an answer.
    """
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    text = (cold[campaign.name] / prompt_module.PROMPT).read_text()
    declared = {artefact: (spec or {}).get(prompt_module.CONVENTIONS_KEY)
                for artefact, spec in task["output_contract"].items()}
    if not any(declared.values()):
        pytest.skip("campaign declares no conventions")
    assert "not fields to report" in text
    for conventions in declared.values():
        for name in conventions or {}:
            assert f"`{name}`" in text, name
