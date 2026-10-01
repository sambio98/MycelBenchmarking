"""Pinned codon and constraint tables for the construct-design template.

These tables are *campaign data*, not claims about biological truth. The gold
artifacts do not depend on the tables being the single correct measurement of
codon usage -- only on their being frozen, hashed and published alongside the
campaign, so that any party recomputes the same gold. The environment supplies
the numeric truth (spec 6.3, G1).

PROVENANCE GAP (must close before the audit packet ships): the relative
adaptiveness values below need a cited source release -- a named Kazusa /
HIVE-CUT release or a computed table over a declared ribosomal-protein gene
set -- recorded in provenance.md with a hash. They are currently
representative values for highly-expressed genes and are flagged as such.
"""

from __future__ import annotations

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

# Relative adaptiveness w in (0, 1]; 1.0 is the most-adapted synonymous codon.
ECOLI_W: dict[str, float] = {
    "GCT": 1.00, "GCC": 0.72, "GCA": 0.59, "GCG": 0.41,
    "CGT": 1.00, "CGC": 0.72, "CGA": 0.07, "CGG": 0.10, "AGA": 0.06, "AGG": 0.03,
    "AAT": 0.42, "AAC": 1.00,
    "GAT": 0.66, "GAC": 1.00,
    "TGT": 0.56, "TGC": 1.00,
    "CAA": 0.31, "CAG": 1.00,
    "GAA": 1.00, "GAG": 0.44,
    "GGT": 1.00, "GGC": 0.92, "GGA": 0.09, "GGG": 0.14,
    "CAT": 0.56, "CAC": 1.00,
    "ATT": 0.56, "ATC": 1.00, "ATA": 0.04,
    "TTA": 0.07, "TTG": 0.09, "CTT": 0.08, "CTC": 0.11, "CTA": 0.02, "CTG": 1.00,
    "AAA": 1.00, "AAG": 0.30,
    "ATG": 1.00,
    "TTT": 0.42, "TTC": 1.00,
    "CCT": 0.30, "CCC": 0.14, "CCA": 0.28, "CCG": 1.00,
    "TCT": 1.00, "TCC": 0.85, "TCA": 0.18, "TCG": 0.20, "AGT": 0.17, "AGC": 0.72,
    "ACT": 0.80, "ACC": 1.00, "ACA": 0.15, "ACG": 0.25,
    "TGG": 1.00,
    "TAT": 0.44, "TAC": 1.00,
    "GTT": 1.00, "GTC": 0.56, "GTA": 0.42, "GTG": 0.56,
    "TAA": 1.00, "TGA": 0.22, "TAG": 0.05,
}

# Actinomycete high-GC preference: strong third-position G/C bias.
STREPTOMYCES_W: dict[str, float] = {
    "GCT": 0.18, "GCC": 1.00, "GCA": 0.14, "GCG": 0.74,
    "CGT": 0.21, "CGC": 1.00, "CGA": 0.12, "CGG": 0.63, "AGA": 0.03, "AGG": 0.08,
    "AAT": 0.17, "AAC": 1.00,
    "GAT": 0.26, "GAC": 1.00,
    "TGT": 0.16, "TGC": 1.00,
    "CAA": 0.14, "CAG": 1.00,
    "GAA": 0.25, "GAG": 1.00,
    "GGT": 0.23, "GGC": 1.00, "GGA": 0.19, "GGG": 0.40,
    "CAT": 0.21, "CAC": 1.00,
    "ATT": 0.14, "ATC": 1.00, "ATA": 0.03,
    "TTA": 0.02, "TTG": 0.13, "CTT": 0.09, "CTC": 0.60, "CTA": 0.03, "CTG": 1.00,
    "AAA": 0.18, "AAG": 1.00,
    "ATG": 1.00,
    "TTT": 0.16, "TTC": 1.00,
    "CCT": 0.11, "CCC": 0.56, "CCA": 0.10, "CCG": 1.00,
    "TCT": 0.08, "TCC": 0.56, "TCA": 0.05, "TCG": 1.00, "AGT": 0.06, "AGC": 0.74,
    "ACT": 0.09, "ACC": 1.00, "ACA": 0.07, "ACG": 0.65,
    "TGG": 1.00,
    "TAT": 0.18, "TAC": 1.00,
    "GTT": 0.12, "GTC": 1.00, "GTA": 0.08, "GTG": 0.84,
    "TAA": 0.25, "TGA": 1.00, "TAG": 0.40,
}

HOSTS: dict[str, dict] = {
    "ecoli_bl21": {
        "weights": ECOLI_W,
        "gc_global": (0.48, 0.555),
        "gc_window": (0.35, 0.65),
        "label": "Escherichia coli BL21(DE3)",
    },
    "streptomyces_coelicolor": {
        "weights": STREPTOMYCES_W,
        "gc_global": (0.62, 0.75),
        "gc_window": (0.45, 0.88),
        "label": "Streptomyces coelicolor A3(2)",
    },
}

# Recognition sites forbidden inside the ORF for the declared cloning strategy.
# All six are palindromic, so one strand suffices.
CLONING_STRATEGIES: dict[str, dict] = {
    "pet28a_ndei_xhoi": {
        "vector": "pET-28a(+)",
        "forbidden_sites": {
            "NdeI": "CATATG",
            "XhoI": "CTCGAG",
            "NcoI": "CCATGG",
            "BamHI": "GGATCC",
            "EcoRI": "GAATTC",
            "HindIII": "AAGCTT",
        },
    },
    "pset152_ndei_xhoi": {
        "vector": "pSET152",
        "forbidden_sites": {
            "NdeI": "CATATG",
            "XhoI": "CTCGAG",
            "EcoRI": "GAATTC",
            "HindIII": "AAGCTT",
        },
    },
}

# Sequence-composition limits, host-independent.
MAX_HOMOPOLYMER = 5
MAX_DIRECT_REPEAT = 10
GC_WINDOW_SIZE = 50
INTERNAL_RBS_MOTIF = "AGGAGG"

# Fixed tie-break order for reporting a binding constraint.
CONSTRAINT_ORDER = (
    "forbidden_sites",
    "gc_window",
    "gc_global",
    "homopolymer",
    "direct_repeat",
    "internal_rbs",
)
