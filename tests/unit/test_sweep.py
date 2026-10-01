"""Internal sweep: gold isolation, statistics, and gate semantics."""

from __future__ import annotations

import json
import pathlib
import shutil

import pytest

from npbench_c.sweep.runner import (
    NEVER_IN_SANDBOX,
    GoldLeak,
    SystemSpec,
    _audit_sandbox,
    build_sandbox,
)
from npbench_c.sweep.stats import (
    cluster_bootstrap_ci,
    leakage_z,
    pass_at_1,
    pass_at_k,
    standard_error,
)

CAMPAIGN = pathlib.Path(__file__).resolve().parents[2] / "campaigns" / "construct-ecoli-rebh-01"
pytestmark = pytest.mark.skipif(not (CAMPAIGN / "gold" / "gold.json").is_file(),
                                reason="campaign gold not generated")


# ------------------------------------------------------------------ isolation

def test_cold_sandbox_contains_no_gold(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "cold")
    names = {p.name for p in s.rglob("*") if p.is_file()}
    assert "gold.json" not in names
    assert "grading.json" not in names
    assert not (s / "oracle").exists()
    assert (s / "inputs" / "rebh.faa").is_file()


def test_warm_bundle_for_r4_withholds_the_r4_answer(tmp_path):
    """The leak that a whole-gold-file warm start would cause.

    Declaring gold/gold.json as an R4 warm input would hand the run the
    Streptomyces block it is being asked to compute, and nothing downstream
    would notice: the run would simply look successful.
    """
    s = build_sandbox(CAMPAIGN, tmp_path / "warm_r4", "gold/warm/r4")
    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())
    strep = gold["hosts"]["streptomyces"]

    # task.yaml names Streptomyces on purpose -- the agent must know what it is
    # being asked -- so the test is on gold *values*, not on the word.
    blob = "\n".join(p.read_text(errors="replace") for p in s.rglob("*") if p.is_file())
    assert json.dumps(strep["unconstrained_violations"], sort_keys=True) not in blob
    assert json.dumps(gold["flipped_constraints"], sort_keys=True) not in blob
    assert not (s / "orf_streptomyces.fasta").exists()

    handed_over = json.loads((s / "analysis_ecoli.json").read_text())
    assert set(handed_over) == {
        "unconstrained_violations", "binding_constraint", "non_binding_constraints"}
    assert handed_over["binding_constraint"] == gold["hosts"]["ecoli"]["binding_constraint"]
    assert handed_over["binding_constraint"] != strep["binding_constraint"]

    # ...while the R3 answer it is entitled to IS present.
    assert (s / "orf_ecoli.fasta").is_file()


def test_audit_catches_a_planted_leak(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "leaky")
    shutil.copytree(CAMPAIGN / "gold", s / "gold")
    with pytest.raises(GoldLeak):
        _audit_sandbox(s)


@pytest.mark.parametrize("name", NEVER_IN_SANDBOX)
def test_never_admissible_filenames_are_rejected_anywhere(tmp_path, name):
    s = build_sandbox(CAMPAIGN, tmp_path / f"leak_{name}")
    (s / name).write_text("{}")
    with pytest.raises(GoldLeak):
        _audit_sandbox(s)


def test_missing_warm_bundle_raises_rather_than_degrading(tmp_path):
    """A missing bundle would silently turn a warm run into a cold one and make
    the per-rung difficulty profile a fiction."""
    with pytest.raises(FileNotFoundError):
        build_sandbox(CAMPAIGN, tmp_path / "nobundle", "gold/warm/r99")


# ------------------------------------------------------------------ statistics

def test_pass_at_1_and_pass_at_k_differ_on_flaky_cells():
    cell = [True, True, False]
    assert pass_at_1(cell) == pytest.approx(2 / 3)
    assert pass_at_k([cell]) == 0.0
    assert pass_at_k([[True, True, True]]) == 1.0


def test_standard_error_of_single_observation_is_zero():
    assert standard_error([0.5]) == 0.0


def test_bootstrap_resamples_clusters_not_observations():
    """One cluster with many observations must not masquerade as many clusters:
    four instantiations of a template share a failure mode."""
    one_big = {"a": [1.0] * 40}
    point, lo, hi = cluster_bootstrap_ci(one_big)
    assert (point, lo, hi) == (1.0, 1.0, 1.0)

    spread = {"a": [1.0], "b": [0.0], "c": [1.0], "d": [0.0]}
    _, lo2, hi2 = cluster_bootstrap_ci(spread)
    assert hi2 - lo2 > 0.4, "four disagreeing clusters must give a wide interval"


def test_bootstrap_is_deterministic():
    c = {"a": [1.0, 0.75], "b": [0.5], "c": [0.25, 0.0]}
    assert cluster_bootstrap_ci(c) == cluster_bootstrap_ci(c)


def test_bootstrap_weights_clusters_equally():
    """A campaign contributing more runs must not dominate the mean."""
    lopsided = {"a": [1.0] * 10, "b": [0.0]}
    point, _, _ = cluster_bootstrap_ci(lopsided)
    assert point == 0.5


def test_leakage_z_flags_a_leaking_campaign_and_clears_a_clean_one():
    _, _, leaking = leakage_z([0.60, 0.58, 0.62], 0.0765)
    _, _, clean = leakage_z([0.08, 0.05, 0.10], 0.0765)
    assert leaking > 2.0
    assert clean <= 2.0


def test_leakage_z_handles_zero_spread():
    """Degenerate SE must not raise: three identical runs are common at n=3."""
    _, se, above = leakage_z([0.5, 0.5, 0.5], 0.1)
    assert se == 0.0 and above == float("inf")
    _, _, below = leakage_z([0.05, 0.05, 0.05], 0.1)
    assert below == float("-inf")


# ------------------------------------------------------------------ provenance

def test_stub_based_sweep_is_not_promoted_to_pass():
    """A fixture must never be laundered into evidence."""
    from npbench_c.readiness.gate import PENDING, _pending_agent_runs

    sweep = CAMPAIGN / "internal_sweep.json"
    if not sweep.is_file():
        pytest.skip("no internal_sweep.json")
    data = json.loads(sweep.read_text())
    if not data.get("_meta", {}).get("stub_based"):
        pytest.skip("sweep on record is not stub-based")
    for check in _pending_agent_runs(CAMPAIGN):
        assert check.status == PENDING
        assert "real systems required" in check.detail


def test_system_spec_defaults_are_conservative():
    s = SystemSpec(name="x", command=["true"])
    assert s.mode == "tooled" and s.is_reference is False and s.is_stub is False


# ------------------------------------------------------------------ solvability

def test_campaign_is_solvable_from_its_sandbox(tmp_path):
    """R3's gold must be recomputable using only files the agent can see.

    This campaign was briefly unsolvable: the codon weight table lived in
    oracle/, which the sandbox excludes, so R3's unconstrained_violations could
    not be derived from anything an agent was given. The stub systems hid it
    because they were handed the oracle directory. The check here reimplements
    the rule from the reference files alone and compares to gold.
    """
    s = build_sandbox(CAMPAIGN, tmp_path / "solvable")
    codon = json.loads((s / "reference" / "codon_tables.json").read_text())
    cloning = json.loads((s / "reference" / "cloning_strategies.json").read_text())
    limits = json.loads((s / "reference" / "composition_limits.json").read_text())

    protein = "".join(
        l.strip() for l in (s / "inputs" / "rebh.faa").read_text().splitlines()[1:]
    ) + "*"

    # Rebuild the genetic code from the weight table's own keys, then apply the
    # documented rule: highest weight, ties to the lexicographically smallest.
    weights = codon["hosts"]["ecoli_bl21"]["weights"]
    from npbench_c.grading import primitives  # noqa: F401  (import sanity only)
    import collections

    stops = {"TAA", "TAG", "TGA"}
    table = {
        "A": "GCT GCC GCA GCG", "R": "CGT CGC CGA CGG AGA AGG", "N": "AAT AAC",
        "D": "GAT GAC", "C": "TGT TGC", "Q": "CAA CAG", "E": "GAA GAG",
        "G": "GGT GGC GGA GGG", "H": "CAT CAC", "I": "ATT ATC ATA",
        "L": "TTA TTG CTT CTC CTA CTG", "K": "AAA AAG", "M": "ATG",
        "F": "TTT TTC", "P": "CCT CCC CCA CCG", "S": "TCT TCC TCA TCG AGT AGC",
        "T": "ACT ACC ACA ACG", "W": "TGG", "Y": "TAT TAC",
        "V": "GTT GTC GTA GTG", "*": " ".join(sorted(stops)),
    }
    nt = "".join(
        min(table[aa].split(), key=lambda c: (-weights[c], c)) for aa in protein
    )

    sites = cloning["strategies"]["pet28a_ndei_xhoi"]["forbidden_sites"]
    n_sites = sum(
        sum(1 for i in range(len(nt)) if nt.startswith(motif, i))
        for motif in sites.values()
    )

    k = limits["max_direct_repeat"] + 1
    counts = collections.Counter(nt[i:i + k] for i in range(len(nt) - k + 1))
    n_repeats = sum(1 for c in counts.values() if c > 1)

    gold = json.loads((CAMPAIGN / "gold" / "gold.json").read_text())
    uv = gold["hosts"]["ecoli"]["unconstrained_violations"]
    assert n_sites == uv["forbidden_sites"]
    assert n_repeats == uv["direct_repeat"]


def test_reference_dir_is_public_and_oracle_dir_is_not(tmp_path):
    s = build_sandbox(CAMPAIGN, tmp_path / "paths")
    assert (s / "reference").is_dir()
    assert not (s / "oracle").exists()
