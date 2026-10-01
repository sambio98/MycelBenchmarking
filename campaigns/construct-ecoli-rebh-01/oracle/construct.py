"""Deterministic construct design and constraint analysis.

Used both to generate gold (oracle) and to measure a submission (measure.py).
Pure: no RNG, no clock, no network. Every tie is broken by an explicit,
documented rule so that two independent runs agree byte for byte.
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence

from tables import (
    AA_CODONS,
    CLONING_STRATEGIES,
    CODON_TO_AA,
    CONSTRAINT_ORDER,
    GC_WINDOW_SIZE,
    HOSTS,
    INTERNAL_RBS_MOTIF,
    MAX_DIRECT_REPEAT,
    MAX_HOMOPOLYMER,
)

MAX_REPAIR_ITERATIONS = 20000


class OracleFailure(RuntimeError):
    """The oracle could not satisfy the declared constraints. Always a campaign
    bug: a campaign whose oracle cannot reach 1.0 does not ship (spec 2.7)."""


# --------------------------------------------------------------------------
# sequence properties
# --------------------------------------------------------------------------

def translate(nt: str) -> str:
    if len(nt) % 3 != 0:
        raise ValueError("sequence length is not a multiple of three")
    return "".join(CODON_TO_AA[nt[i:i + 3]] for i in range(0, len(nt), 3))


def gc_fraction(nt: str) -> float:
    if not nt:
        raise ValueError("empty sequence")
    return sum(1 for c in nt if c in "GC") / len(nt)


def cai(nt: str, weights: Mapping[str, float]) -> float:
    """Codon adaptation index: geometric mean of relative adaptiveness.

    Single-codon families (Met, Trp) and stop codons are excluded, per the
    standard definition -- they carry no choice and would otherwise inflate the
    index. math.fsum fixes the reduction order.
    """
    logs = []
    for i in range(0, len(nt), 3):
        codon = nt[i:i + 3]
        aa = CODON_TO_AA[codon]
        if aa == "*" or len(AA_CODONS[aa]) == 1:
            continue
        w = weights[codon]
        if w <= 0:
            raise ValueError(f"non-positive weight for {codon}")
        logs.append(math.log(w))
    if not logs:
        raise ValueError("no scorable codons")
    return math.exp(math.fsum(logs) / len(logs))


# --------------------------------------------------------------------------
# constraint violations
# --------------------------------------------------------------------------

def _forbidden_site_hits(nt: str, sites: Mapping[str, str]) -> list[tuple[str, int]]:
    hits = []
    for name, motif in sorted(sites.items()):
        start = 0
        while True:
            idx = nt.find(motif, start)
            if idx < 0:
                break
            hits.append((name, idx))
            start = idx + 1
    return sorted(hits, key=lambda t: (t[1], t[0]))


def _homopolymer_runs(nt: str, limit: int) -> list[tuple[int, int]]:
    out, i, n = [], 0, len(nt)
    while i < n:
        j = i
        while j + 1 < n and nt[j + 1] == nt[i]:
            j += 1
        if j - i + 1 > limit:
            out.append((i, j - i + 1))
        i = j + 1
    return out


def _repeat_kmers(nt: str, limit: int) -> list[tuple[str, int]]:
    """Every (limit+1)-mer occurring more than once, with its occurrence count.

    A direct repeat longer than ``limit`` necessarily contains a repeated
    (limit+1)-mer, so this detects all of them in linear time.
    """
    k = limit + 1
    seen: dict[str, int] = {}
    for i in range(0, max(0, len(nt) - k + 1)):
        seen[nt[i:i + k]] = seen.get(nt[i:i + k], 0) + 1
    return sorted((kmer, c) for kmer, c in seen.items() if c > 1)


def _gc_window_violations(nt: str, lo: float, hi: float) -> list[int]:
    w = GC_WINDOW_SIZE
    if len(nt) < w:
        return []
    out = []
    gc = sum(1 for c in nt[:w] if c in "GC")
    for i in range(0, len(nt) - w + 1):
        if i > 0:
            gc += (1 if nt[i + w - 1] in "GC" else 0) - (1 if nt[i - 1] in "GC" else 0)
        frac = gc / w
        if frac < lo or frac > hi:
            out.append(i)
    return out


def _rbs_hits(nt: str) -> list[int]:
    out, start = [], 0
    while True:
        idx = nt.find(INTERNAL_RBS_MOTIF, start)
        if idx < 0:
            return out
        out.append(idx)
        start = idx + 1


def violations(nt: str, host: str, strategy: str) -> dict[str, list]:
    """All constraint violations, keyed by constraint name in CONSTRAINT_ORDER."""
    h = HOSTS[host]
    s = CLONING_STRATEGIES[strategy]
    g_lo, g_hi = h["gc_global"]
    w_lo, w_hi = h["gc_window"]
    gc = gc_fraction(nt)
    return {
        "forbidden_sites": _forbidden_site_hits(nt, s["forbidden_sites"]),
        "gc_window": _gc_window_violations(nt, w_lo, w_hi),
        "gc_global": [] if g_lo <= gc <= g_hi else [round(gc, 6)],
        "homopolymer": _homopolymer_runs(nt, MAX_HOMOPOLYMER),
        "direct_repeat": _repeat_kmers(nt, MAX_DIRECT_REPEAT),
        "internal_rbs": _rbs_hits(nt),
    }


def violation_counts(nt: str, host: str, strategy: str) -> dict[str, int]:
    return {k: len(v) for k, v in violations(nt, host, strategy).items()}


def _violating_positions(nt: str, host: str, strategy: str,
                         active: Iterable[str]) -> list[int]:
    """Nucleotide indices participating in any active violation, ascending."""
    v = violations(nt, host, strategy)
    act = set(active)
    idx: set[int] = set()
    if "forbidden_sites" in act:
        for name, pos in v["forbidden_sites"]:
            motif_len = len(CLONING_STRATEGIES[strategy]["forbidden_sites"][name])
            idx.update(range(pos, pos + motif_len))
    if "internal_rbs" in act:
        for pos in v["internal_rbs"]:
            idx.update(range(pos, pos + len(INTERNAL_RBS_MOTIF)))
    if "homopolymer" in act:
        for pos, length in v["homopolymer"]:
            idx.update(range(pos, pos + length))
    if "direct_repeat" in act:
        k = MAX_DIRECT_REPEAT + 1
        for kmer, _ in v["direct_repeat"]:
            start = 0
            while True:
                p = nt.find(kmer, start)
                if p < 0:
                    break
                idx.update(range(p, p + k))
                start = p + 1
    if "gc_window" in act:
        for pos in v["gc_window"]:
            idx.update(range(pos, pos + GC_WINDOW_SIZE))
    if "gc_global" in act and v["gc_global"]:
        idx.update(range(len(nt)))
    return sorted(i for i in idx if 0 <= i < len(nt))


def _total_violations(nt: str, host: str, strategy: str,
                      active: Iterable[str]) -> int:
    c = violation_counts(nt, host, strategy)
    return sum(c[k] for k in active)


# --------------------------------------------------------------------------
# optimiser
# --------------------------------------------------------------------------

def _ranked_codons(aa: str, weights: Mapping[str, float]) -> tuple[str, ...]:
    """Synonymous codons by descending adaptiveness, ties broken
    lexicographically so the ordering is total and reproducible."""
    return tuple(sorted(AA_CODONS[aa], key=lambda c: (-weights[c], c)))


def unconstrained_optimum(protein: str, host: str) -> str:
    """Each residue assigned its most-adapted codon. The reference point
    against which every constraint's CAI cost is measured."""
    w = HOSTS[host]["weights"]
    return "".join(_ranked_codons(aa, w)[0] for aa in protein)


def optimize(protein: str, host: str, strategy: str,
             active: Sequence[str] | None = None,
             policy: str = "first_improve") -> str:
    """Design a construct satisfying the active constraint subset.

    Greedy initialisation at the unconstrained optimum, then repair. Each
    iteration locates the first violating nucleotide, and tries synonymous
    substitutions at the codons spanning it in descending adaptiveness order,
    accepting the first that strictly reduces the total violation count. On a
    plateau it accepts an equal-count move to an unvisited sequence, which with
    the visited set guarantees termination.
    """
    act = list(CONSTRAINT_ORDER if active is None else active)
    if policy == "reverse_order":
        act = list(reversed(act))
    w = HOSTS[host]["weights"]
    nt = unconstrained_optimum(protein, host)
    if not act:
        return nt

    visited = {nt}
    for _ in range(MAX_REPAIR_ITERATIONS):
        current = _total_violations(nt, host, strategy, act)
        if current == 0:
            return nt
        positions = _violating_positions(nt, host, strategy, act)
        if not positions:
            return nt

        best_improve = None
        best_plateau = None
        # Codon indices touched by violations, nearest-first.
        codon_idxs = sorted({p // 3 for p in positions})
        for ci in codon_idxs:
            aa = protein[ci]
            for codon in _ranked_codons(aa, w):
                if codon == nt[ci * 3:ci * 3 + 3]:
                    continue
                cand = nt[:ci * 3] + codon + nt[ci * 3 + 3:]
                score = _total_violations(cand, host, strategy, act)
                if score < current:
                    cand_cai = 0.0 if policy == "lex_tiebreak" else cai(cand, w)
                    key = (score, -cand_cai, cand)
                    if best_improve is None or key < best_improve[0]:
                        best_improve = (key, cand)
                elif score == current and cand not in visited and best_plateau is None:
                    best_plateau = cand
            if best_improve is not None and policy != "best_improve":
                break

        nxt = best_improve[1] if best_improve else best_plateau
        if nxt is None:
            raise OracleFailure(
                f"stuck with {current} violations under {act}; "
                "constraint set may be unsatisfiable for this protein"
            )
        nt = nxt
        visited.add(nt)

    raise OracleFailure(f"exceeded {MAX_REPAIR_ITERATIONS} repair iterations")


# --------------------------------------------------------------------------
# R3/R4 analysis: which constraint actually costs adaptation
# --------------------------------------------------------------------------

OPTIMIZER_POLICIES = ("first_improve", "best_improve", "reverse_order", "lex_tiebreak")


def constraint_costs(protein: str, host: str, strategy: str,
                     policy: str = "first_improve") -> dict:
    """CAI cost of each constraint (spec: R3 gold, tier G1).

    cost(c) = CAI(optimised without c) - CAI(optimised with everything)

    A common unit across constraints, and a property of (protein, host,
    constraint set) rather than of any particular optimiser. The argmax is
    stable even where the absolute values wobble, which is why the binding
    constraint is graded exactly while the costs are interval-scored against a
    measured width.
    """
    w = HOSTS[host]["weights"]
    full_seq = optimize(protein, host, strategy, list(CONSTRAINT_ORDER), policy)
    cai_full = cai(full_seq, w)

    costs = {}
    for c in CONSTRAINT_ORDER:
        subset = [x for x in CONSTRAINT_ORDER if x != c]
        relaxed = optimize(protein, host, strategy, subset, policy)
        costs[c] = cai(relaxed, w) - cai_full

    binding = min(CONSTRAINT_ORDER, key=lambda c: (-costs[c], CONSTRAINT_ORDER.index(c)))
    return {
        "cai_full": round(cai_full, 6),
        "cai_unconstrained": round(cai(unconstrained_optimum(protein, host), w), 6),
        "costs": {c: round(costs[c], 6) for c in CONSTRAINT_ORDER},
        "binding_constraint": binding,
        "binding_cost": round(costs[binding], 6),
        "sequence": full_seq,
    }
