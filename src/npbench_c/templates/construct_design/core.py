"""Shared constraint engine for the construct-design templates.

Two ladders are built on this engine and must not diverge:

  design/    T-L6-1  design an ORF satisfying the constraint set, then say which
                     constraint actually binds and how that changes with host
  conflict/  T-L6-2  the constraint set is provably unsatisfiable; detect that,
                     name the conflicting pair, and prove it

Pure stdlib and deterministic: no network, no clock, no RNG. Every tie is broken
by an explicit, documented rule so two independent runs agree byte for byte.

Tables are loaded from the campaign's reference/ directory -- the same files the
agent is given -- so the tables gold is computed from and the tables the agent
reads cannot drift.
"""

from __future__ import annotations

import json
import math
import pathlib
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

TEMPLATE_ID = "construct-design"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
MAX_REPAIR_ITERATIONS = 20000
OPTIMIZER_POLICIES = ("first_improve", "best_improve", "reverse_order", "lex_tiebreak")

# Violating codon positions examined per repair iteration in the constructive
# search. Bounds the per-iteration cost while leaving every position reachable,
# because the window rotates with the iteration counter.
CANDIDATE_POSITION_CAP = 48

# The genetic code. Not a tunable table: it belongs in code, unlike codon usage
# weights, which are campaign data.
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
CODON_TO_AA: dict[str, str] = {c: aa for aa, cs in AA_CODONS.items() for c in cs}
STOP_CODONS = frozenset(AA_CODONS["*"])


class OracleFailure(RuntimeError):
    """The oracle could not do what the campaign requires. Always a campaign
    bug: a campaign whose oracle cannot reach 1.0 does not ship."""


def gc_count(codon: str) -> int:
    return sum(1 for c in codon if c in "GC")


@dataclass(frozen=True)
class ConstructEngine:
    """A campaign's pinned tables plus the constraint arithmetic over them."""

    hosts: Mapping[str, Mapping]
    strategies: Mapping[str, Mapping]
    max_homopolymer: int
    max_direct_repeat: int
    gc_window_size: int
    internal_rbs_motif: str
    constraint_order: tuple[str, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "ConstructEngine":
        reference = pathlib.Path(campaign) / "reference"
        codon = json.loads((reference / "codon_tables.json").read_text())
        cloning = json.loads((reference / "cloning_strategies.json").read_text())
        comp = json.loads((reference / "composition_limits.json").read_text())
        hosts = {
            name: {"weights": spec["weights"],
                   "gc_global": tuple(spec["gc_global"]),
                   "gc_window": tuple(spec["gc_window"]),
                   "label": spec["label"]}
            for name, spec in codon["hosts"].items()
        }
        return cls(
            hosts=hosts,
            strategies=cloning["strategies"],
            max_homopolymer=int(comp["max_homopolymer"]),
            max_direct_repeat=int(comp["max_direct_repeat"]),
            gc_window_size=int(comp["gc_window_size"]),
            internal_rbs_motif=comp["internal_rbs_motif"],
            constraint_order=tuple(comp["constraint_order"]),
        )

    # ------------------------------------------------------------- sequences

    def translate(self, nt: str) -> str:
        if len(nt) % 3 != 0:
            raise ValueError("sequence length is not a multiple of three")
        return "".join(CODON_TO_AA[nt[i:i + 3]] for i in range(0, len(nt), 3))

    def gc_fraction(self, nt: str) -> float:
        if not nt:
            raise ValueError("empty sequence")
        return sum(1 for c in nt if c in "GC") / len(nt)

    def cai(self, nt: str, weights: Mapping[str, float]) -> float:
        """Codon adaptation index: geometric mean of relative adaptiveness.

        Single-codon families (Met, Trp) and stops are excluded per the standard
        definition -- they carry no choice and would inflate the index.
        math.fsum fixes the reduction order.
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

    # ------------------------------------------------------------ violations

    def _forbidden_site_hits(self, nt: str, sites: Mapping[str, str]) -> list[tuple[str, int]]:
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

    def _homopolymer_runs(self, nt: str) -> list[tuple[int, int]]:
        out, i, n = [], 0, len(nt)
        while i < n:
            j = i
            while j + 1 < n and nt[j + 1] == nt[i]:
                j += 1
            if j - i + 1 > self.max_homopolymer:
                out.append((i, j - i + 1))
            i = j + 1
        return out

    def _repeat_kmers(self, nt: str) -> list[tuple[str, int]]:
        """Every (limit+1)-mer occurring more than once.

        A direct repeat longer than the limit necessarily contains a repeated
        (limit+1)-mer, so this finds all of them in linear time.
        """
        k = self.max_direct_repeat + 1
        seen: dict[str, int] = {}
        for i in range(0, max(0, len(nt) - k + 1)):
            seen[nt[i:i + k]] = seen.get(nt[i:i + k], 0) + 1
        return sorted((kmer, c) for kmer, c in seen.items() if c > 1)

    def _gc_window_violations(self, nt: str, lo: float, hi: float) -> list[int]:
        w = self.gc_window_size
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

    def _rbs_hits(self, nt: str) -> list[int]:
        out, start = [], 0
        while True:
            idx = nt.find(self.internal_rbs_motif, start)
            if idx < 0:
                return out
            out.append(idx)
            start = idx + 1

    def violations(self, nt: str, host: str, strategy: str) -> dict[str, list]:
        h, s = self.hosts[host], self.strategies[strategy]
        g_lo, g_hi = h["gc_global"]
        w_lo, w_hi = h["gc_window"]
        gc = self.gc_fraction(nt)
        return {
            "forbidden_sites": self._forbidden_site_hits(nt, s["forbidden_sites"]),
            "gc_window": self._gc_window_violations(nt, w_lo, w_hi),
            "gc_global": [] if g_lo <= gc <= g_hi else [round(gc, 6)],
            "homopolymer": self._homopolymer_runs(nt),
            "direct_repeat": self._repeat_kmers(nt),
            "internal_rbs": self._rbs_hits(nt),
        }

    def violation_counts(self, nt: str, host: str, strategy: str) -> dict[str, int]:
        return {k: len(v) for k, v in self.violations(nt, host, strategy).items()}

    def _violating_positions(self, nt: str, host: str, strategy: str,
                             active: Iterable[str]) -> list[int]:
        v = self.violations(nt, host, strategy)
        act = set(active)
        idx: set[int] = set()
        if "forbidden_sites" in act:
            sites = self.strategies[strategy]["forbidden_sites"]
            for name, pos in v["forbidden_sites"]:
                idx.update(range(pos, pos + len(sites[name])))
        if "internal_rbs" in act:
            for pos in v["internal_rbs"]:
                idx.update(range(pos, pos + len(self.internal_rbs_motif)))
        if "homopolymer" in act:
            for pos, length in v["homopolymer"]:
                idx.update(range(pos, pos + length))
        if "direct_repeat" in act:
            k = self.max_direct_repeat + 1
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
                idx.update(range(pos, pos + self.gc_window_size))
        if "gc_global" in act and v["gc_global"]:
            idx.update(range(len(nt)))
        return sorted(i for i in idx if 0 <= i < len(nt))

    def _total_violations(self, nt: str, host: str, strategy: str,
                          active: Iterable[str]) -> int:
        c = self.violation_counts(nt, host, strategy)
        return sum(c[k] for k in active)

    # ------------------------------------------------------------- optimiser

    def ranked_codons(self, aa: str, weights: Mapping[str, float]) -> tuple[str, ...]:
        """Synonymous codons by descending adaptiveness, ties broken
        lexicographically so the ordering is total and reproducible."""
        return tuple(sorted(AA_CODONS[aa], key=lambda c: (-weights[c], c)))

    def unconstrained_optimum(self, protein: str, host: str) -> str:
        w = self.hosts[host]["weights"]
        return "".join(self.ranked_codons(aa, w)[0] for aa in protein)

    def optimize(self, protein: str, host: str, strategy: str,
                 active: Sequence[str] | None = None,
                 policy: str = "first_improve") -> str:
        """Design a construct satisfying the active constraint subset.

        Greedy initialisation at the unconstrained optimum, then repair. Each
        iteration locates the first violating nucleotide and tries synonymous
        substitutions at the codons spanning it in descending adaptiveness
        order, accepting the first that strictly reduces the total violation
        count. On a plateau it accepts an equal-count move to an unvisited
        sequence, which with the visited set guarantees termination.
        """
        act = list(self.constraint_order if active is None else active)
        if policy == "reverse_order":
            act = list(reversed(act))
        w = self.hosts[host]["weights"]
        nt = self.unconstrained_optimum(protein, host)
        if not act:
            return nt

        visited = {nt}
        for _ in range(MAX_REPAIR_ITERATIONS):
            current = self._total_violations(nt, host, strategy, act)
            if current == 0:
                return nt
            positions = self._violating_positions(nt, host, strategy, act)
            if not positions:
                return nt

            best_improve = None
            best_plateau = None
            for ci in sorted({p // 3 for p in positions}):
                aa = protein[ci]
                for codon in self.ranked_codons(aa, w):
                    if codon == nt[ci * 3:ci * 3 + 3]:
                        continue
                    cand = nt[:ci * 3] + codon + nt[ci * 3 + 3:]
                    score = self._total_violations(cand, host, strategy, act)
                    if score < current:
                        cand_cai = 0.0 if policy == "lex_tiebreak" else self.cai(cand, w)
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
                    "constraint set may be unsatisfiable for this protein")
            nt = nxt
            visited.add(nt)

        raise OracleFailure(f"exceeded {MAX_REPAIR_ITERATIONS} repair iterations")

    def constraint_costs(self, protein: str, host: str, strategy: str,
                         policy: str = "first_improve") -> dict:
        """CAI cost of each constraint: the design ladder's R3 gold.

        cost(c) = CAI(optimised without c) - CAI(optimised with everything).
        A common unit across constraints, and a property of (protein, host,
        constraint set) rather than of any particular optimiser. The argmax is
        stable across policies even where absolute values wobble, which is why
        the binding constraint is graded and the magnitudes are not.
        """
        w = self.hosts[host]["weights"]
        full_seq = self.optimize(protein, host, strategy, list(self.constraint_order), policy)
        cai_full = self.cai(full_seq, w)
        costs = {}
        for c in self.constraint_order:
            subset = [x for x in self.constraint_order if x != c]
            costs[c] = self.cai(self.optimize(protein, host, strategy, subset, policy), w) - cai_full
        binding = min(self.constraint_order,
                      key=lambda c: (-costs[c], self.constraint_order.index(c)))
        return {
            "cai_full": round(cai_full, 6),
            "cai_unconstrained": round(
                self.cai(self.unconstrained_optimum(protein, host), w), 6),
            "costs": {c: round(costs[c], 6) for c in self.constraint_order},
            "binding_constraint": binding,
            "binding_cost": round(costs[binding], 6),
            "sequence": full_seq,
        }

    # ------------------------------------------- achievable GC (feasibility)

    def achievable_gc_bounds(self, protein: str,
                             forbidden_codons: Iterable[str] = ()) -> dict:
        """Exact bounds on whole-ORF GC, given a set of disallowed codons.

        This is the arithmetic the conflict ladder's proof rests on. For each
        residue independently, the lowest and highest GC content reachable with
        an allowed codon; summed, divided by the ORF length. Because codon
        choices are independent, these bounds are *attained*, so they are exact
        rather than merely valid.

        A GC requirement above the upper bound cannot be met by ANY sequence
        encoding this protein, whatever other constraints apply -- adding
        constraints only shrinks the feasible set. That makes infeasibility
        provable without search, which is what abstention gold requires: a
        greedy optimiser failing to find a design proves nothing.
        """
        forbidden = set(forbidden_codons)
        lo = hi = 0
        per_aa: dict[str, int] = {}
        starved: list[str] = []
        for aa in sorted(set(protein)):
            allowed = [c for c in AA_CODONS[aa] if c not in forbidden]
            per_aa[aa] = len(allowed)
            if not allowed:
                starved.append(aa)
        if starved:
            return {"feasible_codon_assignment": False,
                    "residues_without_codon": starved,
                    "admissible_codon_counts": per_aa,
                    "gc_min": None, "gc_max": None}
        for aa in protein:
            allowed = [c for c in AA_CODONS[aa] if c not in forbidden]
            lo += min(gc_count(c) for c in allowed)
            hi += max(gc_count(c) for c in allowed)
        n = 3 * len(protein)
        return {
            "feasible_codon_assignment": True,
            "residues_without_codon": [],
            "admissible_codon_counts": per_aa,
            "gc_min": round(lo / n, 6),
            "gc_max": round(hi / n, 6),
        }

    # ------------------------------------- designing under a conflict set

    COMPOSITION_CONSTRAINTS = ("forbidden_sites", "gc_window", "homopolymer",
                               "direct_repeat", "internal_rbs")

    def conflict_violations(self, nt: str, host: str, strategy: str,
                            forbidden_codons: Iterable[str],
                            gc_global_min: float) -> dict[str, int]:
        """Violation counts under a conflict set.

        The host's gc_global window is replaced by the explicit gc_global_min
        requirement, and two constraint types the design ladder does not use are
        added: forbidden_codons and gc_global_min.
        """
        base = self.violation_counts(nt, host, strategy)
        forbidden = set(forbidden_codons)
        codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
        out = {c: base[c] for c in self.COMPOSITION_CONSTRAINTS}
        out["forbidden_codons"] = sum(1 for c in codons if c in forbidden)
        out["gc_global_min"] = int(self.gc_fraction(nt) < gc_global_min)
        return out

    def design_under_conflict_set(self, protein: str, host: str, strategy: str,
                                  forbidden_codons: Iterable[str],
                                  gc_global_min: float) -> str:
        """Construct a sequence satisfying a conflict set, or raise.

        Feasibility of a full constraint set is NOT decidable by the achievable-
        GC arithmetic: that arithmetic settles whether the GC floor is reachable,
        but the composition constraints could still exclude every candidate. So
        feasibility is established CONSTRUCTIVELY -- by exhibiting a sequence --
        which is a proof, where a failure to find one is not.

        The greedy repair used by ``optimize`` cannot climb a global GC floor: a
        single synonymous substitution rarely crosses the threshold, so the
        violation count never strictly improves and the search stalls. This
        starts from the GC-maximising admissible assignment, which satisfies the
        floor whenever it is reachable at all, and then repairs the composition
        constraints while refusing any substitution that would drop GC back
        below the floor.
        """
        forbidden = set(forbidden_codons)
        allowed: dict[str, tuple[str, ...]] = {}
        for aa in sorted(set(protein)):
            options = tuple(c for c in AA_CODONS[aa] if c not in forbidden)
            if not options:
                raise OracleFailure(f"residue {aa!r} has no admissible codon")
            allowed[aa] = options

        # GC-maximising start: highest GC per residue, ties lexicographic.
        nt = "".join(max(allowed[aa], key=lambda c: (gc_count(c), c)) for aa in protein)
        if self.gc_fraction(nt) < gc_global_min:
            raise OracleFailure(
                f"GC floor {gc_global_min} exceeds the achievable maximum "
                f"{self.gc_fraction(nt):.6f}; the set is infeasible, not feasible")

        active = list(self.COMPOSITION_CONSTRAINTS)
        gc_total = sum(1 for c in nt if c in "GC")
        floor_count = math.ceil(gc_global_min * len(nt))
        visited = {nt}

        for step in range(MAX_REPAIR_ITERATIONS):
            counts = self.conflict_violations(nt, host, strategy, forbidden, gc_global_min)
            total = sum(counts[c] for c in active)
            if total == 0:
                if counts["forbidden_codons"] or counts["gc_global_min"]:
                    raise OracleFailure(f"repair broke an invariant: {counts}")
                return nt

            # Candidate positions are the violating codons, scanned in order
            # and capped per iteration. Scanning all of them with a full
            # violation recount per candidate was quadratic in the ORF length
            # and far too slow for the gate's determinism re-run; scanning only
            # the earliest violation's span was too narrow, since a 50-nt GC
            # window cannot be repaired by two codon substitutions. The cap
            # rotates so every position is reachable across iterations.
            positions = self._violating_positions(nt, host, strategy, active)
            if not positions:
                return nt
            codon_idxs = sorted({p // 3 for p in positions})
            if len(codon_idxs) > CANDIDATE_POSITION_CAP:
                offset = step % len(codon_idxs)
                codon_idxs = (codon_idxs[offset:] + codon_idxs[:offset])[:CANDIDATE_POSITION_CAP]

            improve = plateau = None
            for ci in codon_idxs:
                aa = protein[ci]
                current = nt[ci * 3:ci * 3 + 3]
                # Prefer substitutions that keep GC high, so repeated repairs do
                # not erode the floor; ties lexicographic for determinism.
                for codon in sorted(allowed[aa], key=lambda c: (-gc_count(c), c)):
                    if codon == current:
                        continue
                    # GC in O(1) from the codon delta rather than rescanning.
                    if gc_total - gc_count(current) + gc_count(codon) < floor_count:
                        continue
                    cand = nt[:ci * 3] + codon + nt[ci * 3 + 3:]
                    cc = self.conflict_violations(cand, host, strategy, forbidden,
                                                  gc_global_min)
                    score = sum(cc[c] for c in active)
                    if score < total:
                        if improve is None or cand < improve:
                            improve = cand
                    elif score == total and cand not in visited and plateau is None:
                        plateau = cand
                if improve is not None:
                    break

            nxt = improve or plateau
            if nxt is None:
                raise OracleFailure(
                    f"stuck with {total} composition violations while holding GC "
                    f">= {gc_global_min}; cannot prove this set feasible")
            nt = nxt
            gc_total = sum(1 for c in nt if c in "GC")
            visited.add(nt)
        raise OracleFailure(f"exceeded {MAX_REPAIR_ITERATIONS} repair iterations")
