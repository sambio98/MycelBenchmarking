"""Pinned codon and constraint tables for the construct-design template.

The values live in ``../reference/`` as JSON and are loaded here, so the tables
the oracle computes gold from and the tables the agent is given are the *same
files*. Duplicating them would let the two drift, and a drift in the codon
weights silently changes every R3 answer.

These tables are campaign data, not claims about biological truth. Gold does
not depend on them being the single correct measurement of codon usage, only on
their being frozen, hashed and published alongside the campaign.

PROVENANCE GAP (must close before the audit packet ships): the relative
adaptiveness values need a cited source release -- a named Kazusa / HIVE-CUT
release, or a table computed over a declared ribosomal-protein gene set --
recorded with a hash. They are currently representative values for
highly-expressed genes and are flagged as such in the JSON.
"""

from __future__ import annotations

import json
import pathlib

_REFERENCE = pathlib.Path(__file__).resolve().parent.parent / "reference"


def _load(name: str) -> dict:
    return json.loads((_REFERENCE / name).read_text())


_CODON = _load("codon_tables.json")
_CLONING = _load("cloning_strategies.json")
_COMPOSITION = _load("composition_limits.json")

AA_CODONS: dict[str, tuple[str, ...]] = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "N": ("AAT", "AAC"),
    "D": ("GAT", "GAC"),
    "C": ("TGT", "TGC"),
    "Q": ("CAA", "CAG"),
    "E": ("GAA", "GAG"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "K": ("AAA", "AAG"),
    "M": ("ATG",),
    "F": ("TTT", "TTC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "*": ("TAA", "TGA", "TAG"),
}

CODON_TO_AA: dict[str, str] = {
    codon: aa for aa, codons in AA_CODONS.items() for codon in codons
}

ECOLI_W: dict[str, float] = _CODON["hosts"]["ecoli_bl21"]["weights"]
STREPTOMYCES_W: dict[str, float] = _CODON["hosts"]["streptomyces_coelicolor"]["weights"]

HOSTS: dict[str, dict] = {
    name: {
        "weights": spec["weights"],
        "gc_global": tuple(spec["gc_global"]),
        "gc_window": tuple(spec["gc_window"]),
        "label": spec["label"],
    }
    for name, spec in _CODON["hosts"].items()
}

CLONING_STRATEGIES: dict[str, dict] = _CLONING["strategies"]

MAX_HOMOPOLYMER = _COMPOSITION["max_homopolymer"]
MAX_DIRECT_REPEAT = _COMPOSITION["max_direct_repeat"]
GC_WINDOW_SIZE = _COMPOSITION["gc_window_size"]
INTERNAL_RBS_MOTIF = _COMPOSITION["internal_rbs_motif"]
CONSTRAINT_ORDER = tuple(_COMPOSITION["constraint_order"])
