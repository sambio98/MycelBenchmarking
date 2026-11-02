"""Gene calls under assembly fragmentation, and the flag that matters more.

A gene caller reads a contig and emits coding-sequence predictions. Cut the same
contig into pieces and run it again and the answer changes -- genes crossing a
cut cannot be called whole. That much is geometry. What this campaign measures is
how much of the change is geometry and how much is the tool reacting to its
input, and the answer on this genome is that a single command-line flag moves
more calls than cutting the chromosome into fifty pieces.

Five things are load-bearing.

**Nothing annotated is shipped.** A gene-calling campaign cannot hand the agent a
record carrying CDS features: that is the answer. The input is the sequence and
nothing else, which is why this campaign ships a FASTA where the detection
campaign on the same chromosome ships a GenBank.

**The fragmentation is a declared RULE, not a shipped file.** The breakpoints are
arithmetic over the sequence length, so the agent cuts the contig itself and the
campaign ships one file. Shipping pre-cut fragments would hand over the offsets,
and the offsets are the coordinate contract.

**The coordinate contract is the trap.** The tool reports 1-based inclusive
positions ON THE FRAGMENT, and comparing those with whole-contig calls requires
adding each fragment's offset. A reader who skips that compares fragment
coordinates against chromosome coordinates and reports almost everything as lost.

**The mode is pinned and varied.** Prodigal's single mode trains on the input it
is given; its meta mode uses pre-trained models. So fragmentation changes the
training set in one mode and not the other, and the campaign separates the local
effect -- calls at a cut -- from the global one, which only retraining produces.

**Losses are classified by geometry, not by guesswork.** A lost call either spans
a cut, lies within a declared window of one, or is distant from every cut. The
third term is the retraining effect and it is empty in the mode that does not
retrain, which is the campaign's central claim and a measurement rather than an
assumption.

Everything outside the tool calls is pure stdlib and sorted order throughout.
"""

from __future__ import annotations

import collections
import gzip
import hashlib
import json
import os
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.tools import cache
from npbench_c.tools.resources import DEFAULT_PREFIX, tool_prefix

TEMPLATE_ID = "gene-calling"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "calling_rules.json"
PRECISION = 6

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 3600

#: The nucleotides a contig may contain. Anything else is a different input than
#: the campaign declares, and the tool's handling of it would be the subject.
ALPHABET = frozenset("ACGTN")


class CallingError(RuntimeError):
    """The contig or the environment does not support the experiment."""


def _run(argv: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a pinned tool with the prefix appended to PATH, last."""
    prefix = str(tool_prefix())
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if prefix in parts:
        parts.remove(prefix)
    env = {**os.environ, **FIXED_ENV, "PATH": os.pathsep.join([*parts, prefix])}
    kwargs.setdefault("stdout", subprocess.PIPE)
    kwargs.setdefault("stderr", subprocess.PIPE)
    return subprocess.run(list(argv), text=True, env=env, timeout=TIMEOUT_S,
                          **kwargs)


# ----------------------------------------------------------------------- rules


@dataclass(frozen=True)
class Rules:
    modes: tuple[str, ...]
    baseline_mode: str
    mode_flags: Mapping[str, tuple[str, ...]]
    fragment_counts: tuple[int, ...]
    baseline_fragments: int
    locality_window: int
    loss_verdicts: tuple[str, ...]
    call_columns: tuple[str, ...]
    partial_flags: tuple[str, ...]
    fixed_flags: tuple[str, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        a = spec["calling"]
        modes = tuple(a["modes"])
        flags = {name: tuple(f) for name, f in sorted(a["mode_flags"].items())}
        if set(flags) != set(modes):
            raise CallingError(
                f"the declared modes {sorted(modes)} and the mode flags "
                f"{sorted(flags)} do not describe the same set, so a mode would "
                "run with flags nobody declared")
        if a["baseline_mode"] not in modes:
            raise CallingError(
                f"the baseline mode {a['baseline_mode']!r} is not declared")
        counts = tuple(int(n) for n in a["fragment_counts"])
        if 1 not in counts:
            raise CallingError(
                "the fragment grid does not include 1, so there is no unfragmented "
                "arm and every comparison would be between two cut assemblies")
        if a["baseline_fragments"] not in counts:
            raise CallingError(
                f"the baseline fragmentation {a['baseline_fragments']} is not on "
                f"the declared grid {list(counts)}")
        if int(a["locality_window"]) <= 0:
            raise CallingError("a non-positive locality window classifies nothing")
        return cls(modes=modes, baseline_mode=a["baseline_mode"], mode_flags=flags,
                   fragment_counts=tuple(sorted(counts)),
                   baseline_fragments=int(a["baseline_fragments"]),
                   locality_window=int(a["locality_window"]),
                   loss_verdicts=tuple(a["loss_verdicts"]),
                   call_columns=tuple(a["call_columns"]),
                   partial_flags=tuple(a["partial_flags"]),
                   fixed_flags=tuple(a["fixed_flags"]))

    # ------------------------------------------------------ the declared rule

    def breakpoints(self, length: int, fragments: int) -> tuple[int, ...]:
        """The interior cut positions for a declared fragment count.

        0-based offsets of the first base of every fragment after the first, so
        fragment k runs from `breakpoints[k-1]` to `breakpoints[k]`. Evenly
        spaced by count rather than by a fixed size, so the rule reads the same
        for any contig and the grid is a count rather than a length nobody can
        check against this sequence.
        """
        if fragments < 1:
            raise CallingError(f"{fragments} fragments is not an assembly")
        if fragments > length:
            raise CallingError(
                f"{fragments} fragments of a {length} base contig would be "
                "shorter than a codon")
        return tuple(round(k * length / fragments) for k in range(1, fragments))

    def fragment_name(self, index: int) -> str:
        return f"frag_{index:04d}"


# ------------------------------------------------------------------- the contig


def read_contig(path: pathlib.Path) -> tuple[str, str]:
    """The one sequence in a FASTA, as (identifier, bases).

    One record, because the campaign's whole subject is what happens when one
    contig becomes several: an input that is already fragmented has no
    unfragmented arm to compare against.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        text = handle.read()
    identifier: str | None = None
    buffer: list[str] = []
    records = 0
    for line in text.splitlines():
        if line.startswith(">"):
            records += 1
            if records > 1:
                raise CallingError(
                    f"{path.name} holds more than one record; this campaign's "
                    "input is a single unfragmented contig")
            identifier = line[1:].split()[0]
        else:
            buffer.append(line.strip().upper())
    sequence = "".join(buffer)
    if not identifier or not sequence:
        raise CallingError(f"{path.name} holds no sequence")
    unknown = set(sequence) - ALPHABET
    if unknown:
        raise CallingError(
            f"{path.name} contains {sorted(unknown)}, which is not in the "
            f"declared alphabet {''.join(sorted(ALPHABET))}")
    return identifier, sequence


@dataclass(frozen=True)
class Call:
    """One coding-sequence prediction, in PARENT coordinates."""

    start: int
    end: int
    strand: str
    partial: str
    start_type: str
    fragment: str

    @property
    def key(self) -> tuple[int, int, str]:
        """What makes two calls the same call.

        Coordinates and strand, and deliberately not the partial flag or the
        start type: those are properties the fragmentation CHANGES about a call
        that is otherwise in the same place, and folding them into identity would
        report every truncated call as a different gene rather than the same gene
        truncated.
        """
        return (self.start, self.end, self.strand)

    @property
    def length(self) -> int:
        return self.end - self.start + 1


# ------------------------------------------------------------------ the engine


@dataclass(frozen=True)
class Calling:
    rules: Rules
    accession: str
    sequence: str
    input_path: pathlib.Path

    @classmethod
    def load(cls, campaign: pathlib.Path, relative: str) -> "Calling":
        campaign = pathlib.Path(campaign)
        accession, sequence = read_contig(campaign / relative)
        return cls(rules=Rules.load(campaign), accession=accession,
                   sequence=sequence, input_path=campaign / relative)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def sequence_sha256(self) -> str:
        """The bases, not the file: a FASTA rewrapped at a different width is the
        same contig and must fingerprint the same."""
        return hashlib.sha256(self.sequence.encode()).hexdigest()

    def gc_fraction(self) -> float:
        gc = self.sequence.count("G") + self.sequence.count("C")
        return round(gc / self.length, PRECISION)

    # ------------------------------------------------------------ the tool

    def prodigal_version(self) -> str:
        done = _run([str(tool_prefix() / "prodigal"), "-v"])
        text = (done.stderr or "") + (done.stdout or "")
        for line in text.splitlines():
            if line.strip().startswith("Prodigal"):
                return line.strip().rstrip(":").split(":")[0].strip()
        raise CallingError("prodigal did not report a version")

    def fasta(self, fragments: int) -> str:
        """The contig cut into a declared number of pieces, as FASTA.

        The header carries the fragment's name and nothing else. The offset is
        NOT written into it: it follows from the declared breakpoint rule, and
        putting it in the file would hand over the coordinate contract the
        campaign exists to test.
        """
        cuts = self.rules.breakpoints(self.length, fragments)
        bounds = list(zip([0, *cuts], [*cuts, self.length]))
        pieces = []
        for index, (lo, hi) in enumerate(bounds, 1):
            piece = self.sequence[lo:hi]
            wrapped = "\n".join(piece[i:i + 70] for i in range(0, len(piece), 70))
            pieces.append(f">{self.rules.fragment_name(index)}\n{wrapped}\n")
        return "".join(pieces)

    def offsets(self, fragments: int) -> dict[str, int]:
        cuts = self.rules.breakpoints(self.length, fragments)
        return {self.rules.fragment_name(i): lo
                for i, lo in enumerate([0, *cuts], 1)}

    def calls(self, mode: str | None = None,
              fragments: int | None = None) -> tuple[Call, ...]:
        """One gene-calling run, projected to the declared columns.

        Memoised when and only when a tool memo is enabled, which the sweep does
        for stub systems alone. The key carries the sequence, the tool version,
        the mode's flags and the fragment count -- everything the answer depends
        on.
        """
        chosen = mode or self.rules.baseline_mode
        count = self.rules.baseline_fragments if fragments is None else fragments
        if chosen not in self.rules.modes:
            raise CallingError(f"{chosen!r} is not a declared mode")
        if count not in self.rules.fragment_counts:
            raise CallingError(
                f"{count} is not on the declared fragment grid "
                f"{list(self.rules.fragment_counts)}")

        def compute() -> list[dict]:
            return [c.__dict__ for c in self._call(chosen, count)]

        if not cache.enabled():
            return self._call(chosen, count)
        memo = cache.memoize_json(
            ["gene_calling.calls", self.sequence_sha256(),
             self.prodigal_version(), chosen,
             " ".join(self.rules.mode_flags[chosen]),
             " ".join(self.rules.fixed_flags), str(count)],
            compute)
        return tuple(Call(**row) for row in memo)

    def _call(self, mode: str, fragments: int) -> tuple[Call, ...]:
        """The tool call itself. Never memoised."""
        offsets = self.offsets(fragments)
        with tempfile.TemporaryDirectory(prefix="npbench-calling-") as tmp:
            work = pathlib.Path(tmp)
            contig = work / "contig.fna"
            contig.write_text(self.fasta(fragments))
            gff = work / "calls.gff"
            done = _run([str(tool_prefix() / "prodigal"),
                         "-i", str(contig), "-o", str(gff),
                         *self.rules.mode_flags[mode], *self.rules.fixed_flags])
            if done.returncode != 0:
                raise CallingError(
                    f"prodigal failed in {mode} mode on {fragments} fragments: "
                    f"{(done.stderr or '')[-600:]}")
            calls = _parse_gff(gff.read_text(), offsets, self.rules)
        if not calls:
            raise CallingError(
                f"no call at all in {mode} mode on {fragments} fragments, which "
                "is a broken run rather than a finding about fragmentation")
        return calls

    # ------------------------------------------------------- the inventory

    def call_table(self, calls: Sequence[Call]) -> list[dict]:
        return [{c: getattr(call, c) for c in self.rules.call_columns}
                for call in calls]

    def call_census(self, calls: Sequence[Call]) -> dict:
        lengths = sorted(call.length for call in calls)
        strands = collections.Counter(call.strand for call in calls)
        coding = sum(lengths)
        return {
            "calls": len(calls),
            "on_the_forward_strand": strands.get("+", 0),
            "on_the_reverse_strand": strands.get("-", 0),
            "shortest": lengths[0],
            "median_length": _median(lengths),
            "longest": lengths[-1],
            "coding_bases": coding,
            "coding_fraction": round(coding / self.length, PRECISION),
        }

    def partial_census(self, calls: Sequence[Call]) -> dict[str, int]:
        counts = collections.Counter(call.partial for call in calls)
        unknown = set(counts) - set(self.rules.partial_flags)
        if unknown:
            raise CallingError(
                f"the tool reported partial flags outside the declared "
                f"vocabulary: {sorted(unknown)}")
        return {flag: counts.get(flag, 0) for flag in self.rules.partial_flags}

    def start_type_census(self, calls: Sequence[Call]) -> dict[str, int]:
        counts = collections.Counter(call.start_type for call in calls)
        return dict(sorted(counts.items()))

    def span(self, calls: Sequence[Call]) -> dict:
        starts = sorted(call.start for call in calls)
        ends = sorted(call.end for call in calls)
        return {"first_call_start": starts[0], "last_call_end": ends[-1],
                "bases_before_the_first_call": starts[0] - 1,
                "bases_after_the_last_call": self.length - ends[-1]}

    # -------------------------------------------------------- the comparison

    def loss_verdict(self, call: Call, cuts: Sequence[int]) -> str:
        """Where a lost call sits relative to the cuts, by arithmetic alone.

        `cuts` are 0-based offsets of a fragment's first base and the call is
        1-based inclusive, so the cut falls strictly inside the call when
        `start <= cut < end` in those two conventions -- which is the one place
        the coordinate contract has to be got right inside the engine rather
        than by the agent.
        """
        if not cuts:
            raise CallingError(
                "an unfragmented assembly has no cut to classify a loss against")
        if any(call.start <= cut < call.end for cut in cuts):
            return "crosses_a_cut"
        nearest = min(min(abs(call.start - cut), abs(call.end - cut))
                      for cut in cuts)
        if nearest <= self.rules.locality_window:
            return "near_a_cut"
        return "distant_from_every_cut"

    def reconcile(self, baseline: Sequence[Call], fragmented: Sequence[Call],
                  cuts: Sequence[int]) -> dict:
        """The whole-contig call set against a fragmented one.

        Both sides are in parent coordinates, so a call that survives
        fragmentation is the same key on both. `retained_but_truncated` is the
        interesting middle: the same coordinates, and the tool now flags the call
        as running off the end of its fragment.
        """
        base = {call.key: call for call in baseline}
        frag = {call.key: call for call in fragmented}
        lost = [base[k] for k in sorted(base.keys() - frag.keys())]
        gained = [frag[k] for k in sorted(frag.keys() - base.keys())]
        shared = sorted(base.keys() & frag.keys())
        truncated = [frag[k] for k in shared if frag[k].partial != base[k].partial]
        verdicts = collections.Counter(
            self.loss_verdict(call, cuts) for call in lost)
        unknown = set(verdicts) - set(self.rules.loss_verdicts)
        if unknown:
            raise CallingError(
                f"loss verdicts outside the declared vocabulary: {sorted(unknown)}")
        return {
            "baseline_calls": len(baseline),
            "fragmented_calls": len(fragmented),
            "retained": len(shared),
            "lost": len(lost),
            "gained": len(gained),
            "retained_but_truncated": len(truncated),
            "loss_verdicts": {v: verdicts.get(v, 0) for v in self.rules.loss_verdicts},
            "lost_table": [{**{c: getattr(call, c) for c in self.rules.call_columns},
                            "verdict": self.loss_verdict(call, cuts)}
                           for call in lost],
            "gained_partial_census": self.partial_census(gained) if gained else {
                flag: 0 for flag in self.rules.partial_flags},
        }

    # ------------------------------------------------------- the counterfactual

    def grid(self) -> dict:
        """Every declared (mode, fragment count) against that mode's own whole
        contig.

        Against its OWN whole contig and not the baseline mode's, because the two
        effects are independent: comparing meta-on-fragments with
        single-on-whole would add the mode difference to the fragmentation
        difference and report their sum as fragmentation.
        """
        out: dict[str, dict] = {}
        for mode in self.rules.modes:
            whole = self.calls(mode, 1)
            for count in self.rules.fragment_counts:
                cuts = self.rules.breakpoints(self.length, count)
                calls = self.calls(mode, count)
                if count == 1:
                    block = {
                        "baseline_calls": len(whole), "fragmented_calls": len(calls),
                        "retained": len(calls), "lost": 0, "gained": 0,
                        "retained_but_truncated": 0,
                        "loss_verdicts": {v: 0 for v in self.rules.loss_verdicts},
                    }
                else:
                    block = {k: v for k, v in
                             self.reconcile(whole, calls, cuts).items()
                             if k not in ("lost_table", "gained_partial_census")}
                out[f"{mode}|{count}"] = block
        return out

    def mode_difference(self) -> dict:
        """The two modes on the SAME unfragmented contig.

        The campaign's headline: this is the comparison that isolates the flag
        from the assembly, because the input bytes are identical on both sides.
        """
        out = {}
        reference = self.rules.baseline_mode
        base = {call.key for call in self.calls(reference, 1)}
        for mode in self.rules.modes:
            if mode == reference:
                continue
            other = {call.key for call in self.calls(mode, 1)}
            out[f"{reference}_vs_{mode}"] = {
                "shared": len(base & other),
                "only_in_" + reference: len(base - other),
                "only_in_" + mode: len(other - base),
            }
        return out


# --------------------------------------------------------------------- parsing


def _parse_gff(text: str, offsets: Mapping[str, int], rules: Rules) -> tuple[Call, ...]:
    """Prodigal's GFF, projected and lifted to parent coordinates.

    The tool writes 1-based inclusive positions on the fragment; the parent
    position is that plus the fragment's 0-based offset. Columns 4 and 5 are the
    start and the end, column 7 the strand, and the attributes carry `partial`
    and `start_type`. The scores are dropped: a float's last digits are not what
    a call is, and nothing downstream needs them to say where a gene is.
    """
    calls: list[Call] = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        name = fields[0]
        if name not in offsets:
            raise CallingError(
                f"the tool reported a call on {name!r}, which is not one of the "
                f"declared fragments; the fragmentation and the run disagree")
        offset = offsets[name]
        attributes = dict(
            kv.split("=", 1) for kv in fields[8].rstrip(";").split(";") if "=" in kv)
        calls.append(Call(
            start=int(fields[3]) + offset,
            end=int(fields[4]) + offset,
            strand=fields[6],
            partial=attributes.get("partial", "??"),
            start_type=attributes.get("start_type", "?"),
            fragment=name,
        ))
    return tuple(sorted(calls, key=lambda c: (c.start, c.end, c.strand)))


def _median(ordered: Sequence[int]) -> float:
    n = len(ordered)
    if not n:
        raise CallingError("median of an empty set")
    if n % 2:
        return float(ordered[n // 2])
    return round((ordered[n // 2 - 1] + ordered[n // 2]) / 2, PRECISION)
