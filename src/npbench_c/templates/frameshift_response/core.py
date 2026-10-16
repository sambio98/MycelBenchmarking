"""What a gene caller does to a gene with one base taken out of it.

A frameshift is the cheapest way to break a gene: delete one base and everything
downstream reads in the wrong frame until a stop codon turns up. A caller that
saw the intact gene as one long ORF now sees something else, and what it sees is
the measurement.

This campaign PLANTS the frameshift. That matters, because "which locus is a
pseudogene" is a judgement nobody can grade: the literature disagrees, the
annotation pipelines disagree, and there is no oracle. A frameshift the campaign
introduced itself is known by construction, so what is graded is not whether a
gene is broken but what the caller reports when it demonstrably is.

Five things are load-bearing.

**The mode is meta, and the reason is measured rather than conventional.** The
companion gene-calling campaign pins Prodigal's single mode because its
retraining is the finding there. Here retraining is NOISE: single mode trains on
the input, so deleting ONE base out of 8.7 megabases changes eleven to twenty-
seven calls across the chromosome, and the planted frameshift's own signature --
one call lost, two gained -- is unrecoverable inside that. Meta mode does not
train, and the same deletion changes exactly the calls it should. Two campaigns
on one tool, opposite modes, each for a stated measurement.

**The coordinate contract is a lift past the deletion.** Delete a base at
position p and every mutated coordinate at or after p corresponds to a parent
coordinate one higher. Skip that and every call downstream of the cut looks moved
by one, which turns a surgical result into thousands of spurious differences.

**The verdict is classified by overlap, not by a window.** A frameshifted gene can
produce an ORF running past where the original ended, so any fixed window around
the target misclassifies it. The classification asks which new calls overlap the
target's own span, which needs no window and no tolerance.

**The target selection is a declared rule over the baseline calls.** The campaign
picks the call closest to each declared length, so the target set is reproducible
from the caller's own output and needs no annotation.

**And the verdict is not a function of length.** Measured on this chromosome, a
189 bp gene survives as one truncated ORF while 219 and 261 bp genes are
abolished outright. Whether each half of a broken gene is callable depends on a
start codon and a ribosome-binding site falling in the right place, not only on
how much room there is, so the campaign reports the verdicts rather than a
threshold.

That last point decides the ladder's shape, and it is worth knowing why. An
earlier draft declared the usual round lengths -- 100, 150, 200, 250, 300, 450
and up -- and the verdict list came out in three contiguous blocks, which is two
cut points in disguise: R3 would have graded one threshold restated once per
target. The round numbers step straight over the band where the caller actually
turns over. The declared ladder samples that band at a uniform ten-base step
instead, and the build refuses gold if the list is contiguous anyway.

Everything outside the caller runs is pure stdlib and sorted order throughout.
"""

from __future__ import annotations

import collections
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
from npbench_c.templates.gene_calling.core import read_contig

TEMPLATE_ID = "frameshift-response"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "frameshift_rules.json"
PRECISION = 6

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 3600


class FrameshiftError(RuntimeError):
    """The contig or the environment does not support the experiment."""


def _run(argv: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
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
    noise_mode: str
    mode_flags: Mapping[str, tuple[str, ...]]
    fixed_flags: tuple[str, ...]
    target_lengths: tuple[int, ...]
    indel_fraction: float
    offset_fractions: tuple[float, ...]
    offset_target_length: int
    noise_target_lengths: tuple[int, ...]
    verdicts: tuple[str, ...]
    call_columns: tuple[str, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        a = spec["frameshift"]
        modes = tuple(a["modes"])
        flags = {name: tuple(f) for name, f in sorted(a["mode_flags"].items())}
        if set(flags) != set(modes):
            raise FrameshiftError(
                f"the declared modes {sorted(modes)} and the mode flags "
                f"{sorted(flags)} do not describe the same set")
        for field in ("baseline_mode", "noise_mode"):
            if a[field] not in modes:
                raise FrameshiftError(f"{field} {a[field]!r} is not a declared mode")
        if a["baseline_mode"] == a["noise_mode"]:
            raise FrameshiftError(
                "the baseline and noise modes are the same mode, so the "
                "counterfactual that shows the signature being swamped compares "
                "a run with itself")
        lengths = tuple(int(n) for n in a["target_lengths"])
        if len(set(lengths)) != len(lengths):
            raise FrameshiftError(
                f"the declared target lengths repeat: {sorted(lengths)}")
        if len(lengths) < 5:
            raise FrameshiftError(
                f"{len(lengths)} declared target lengths cannot span a range")
        fraction = float(a["indel_fraction"])
        if not 0.0 < fraction < 1.0:
            raise FrameshiftError(
                f"an indel fraction of {fraction} is not inside the gene")
        offsets = tuple(float(f) for f in a["offset_fractions"])
        if not all(0.0 < f < 1.0 for f in offsets):
            raise FrameshiftError(
                "every declared offset fraction must fall inside the gene")
        if fraction not in offsets:
            raise FrameshiftError(
                f"the baseline indel fraction {fraction} is not on the declared "
                f"offset grid {list(offsets)}, so the sweep has no shared point "
                "with the rest of the campaign")
        if int(a["offset_target_length"]) not in lengths:
            raise FrameshiftError(
                "the offset sweep's target length is not among the declared "
                "target lengths, so the sweep runs on a gene nothing else uses")
        noise = tuple(int(n) for n in a["noise_target_lengths"])
        if not set(noise) <= set(lengths):
            raise FrameshiftError(
                "the noise arm's target lengths must be a subset of the declared "
                "target lengths")
        if not noise:
            raise FrameshiftError("the noise arm has no target")
        return cls(modes=modes, baseline_mode=a["baseline_mode"],
                   noise_mode=a["noise_mode"], mode_flags=flags,
                   fixed_flags=tuple(a["fixed_flags"]),
                   target_lengths=tuple(sorted(lengths)),
                   indel_fraction=fraction,
                   offset_fractions=tuple(sorted(offsets)),
                   offset_target_length=int(a["offset_target_length"]),
                   noise_target_lengths=tuple(sorted(noise)),
                   verdicts=tuple(a["verdicts"]),
                   call_columns=tuple(a["call_columns"]))


# ------------------------------------------------------------------- the calls


@dataclass(frozen=True)
class Call:
    start: int
    end: int
    strand: str
    partial: str

    @property
    def key(self) -> tuple[int, int, str]:
        return (self.start, self.end, self.strand)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def overlaps(self, start: int, end: int) -> bool:
        return self.start <= end and self.end >= start


def parse_gff(text: str) -> tuple[Call, ...]:
    """Prodigal's GFF over a single contig, projected to the declared columns."""
    calls: list[Call] = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        attributes = dict(
            kv.split("=", 1) for kv in fields[8].rstrip(";").split(";") if "=" in kv)
        calls.append(Call(start=int(fields[3]), end=int(fields[4]),
                          strand=fields[6],
                          partial=attributes.get("partial", "??")))
    return tuple(sorted(calls, key=lambda c: c.key))


# ------------------------------------------------------------------ the engine


@dataclass(frozen=True)
class Frameshift:
    rules: Rules
    accession: str
    sequence: str

    @classmethod
    def load(cls, campaign: pathlib.Path, relative: str) -> "Frameshift":
        campaign = pathlib.Path(campaign)
        accession, sequence = read_contig(campaign / relative)
        return cls(rules=Rules.load(campaign), accession=accession,
                   sequence=sequence)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def sequence_sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode()).hexdigest()

    # --------------------------------------------------------------- the tool

    def prodigal_version(self) -> str:
        done = _run([str(tool_prefix() / "prodigal"), "-v"])
        text = (done.stderr or "") + (done.stdout or "")
        for line in text.splitlines():
            if line.strip().startswith("Prodigal"):
                return line.strip().rstrip(":").split(":")[0].strip()
        raise FrameshiftError("prodigal did not report a version")

    def _call(self, sequence: str, mode: str) -> tuple[Call, ...]:
        with tempfile.TemporaryDirectory(prefix="npbench-frameshift-") as tmp:
            work = pathlib.Path(tmp)
            contig = work / "contig.fna"
            contig.write_text(
                f">{self.accession}\n"
                + "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70))
                + "\n")
            gff = work / "calls.gff"
            done = _run([str(tool_prefix() / "prodigal"),
                         "-i", str(contig), "-o", str(gff),
                         *self.rules.mode_flags[mode], *self.rules.fixed_flags])
            if done.returncode != 0:
                raise FrameshiftError(
                    f"prodigal failed in {mode} mode: {(done.stderr or '')[-600:]}")
            calls = parse_gff(gff.read_text())
        if not calls:
            raise FrameshiftError(
                f"no call at all in {mode} mode, which is a broken run rather "
                "than a finding about frameshifts")
        return calls

    def calls(self, mode: str | None = None,
              deleted_at: int | None = None) -> tuple[Call, ...]:
        """The caller's output, in PARENT coordinates.

        With `deleted_at`, the contig is mutated by removing that one 1-based
        position and the result is lifted back: a mutated coordinate at or after
        the deleted position corresponds to a parent coordinate one higher. That
        lift is this campaign's coordinate contract, and skipping it turns a
        surgical result into thousands of spurious differences.

        Memoised when and only when a tool memo is enabled, which the sweep does
        for stub systems alone.
        """
        chosen = mode or self.rules.baseline_mode
        if chosen not in self.rules.modes:
            raise FrameshiftError(f"{chosen!r} is not a declared mode")
        if deleted_at is not None and not 1 <= deleted_at <= self.length:
            raise FrameshiftError(
                f"position {deleted_at} is outside the contig")

        def compute() -> list[dict]:
            return [c.__dict__ for c in self._lifted(chosen, deleted_at)]

        if not cache.enabled():
            return self._lifted(chosen, deleted_at)
        memo = cache.memoize_json(
            ["frameshift_response.calls", self.sequence_sha256(),
             self.prodigal_version(), chosen,
             " ".join(self.rules.mode_flags[chosen]),
             " ".join(self.rules.fixed_flags),
             "whole" if deleted_at is None else str(deleted_at)],
            compute)
        return tuple(Call(**row) for row in memo)

    def _lifted(self, mode: str, deleted_at: int | None) -> tuple[Call, ...]:
        if deleted_at is None:
            return self._call(self.sequence, mode)
        mutated = self.sequence[:deleted_at - 1] + self.sequence[deleted_at:]
        raw = self._call(mutated, mode)
        lifted = tuple(
            Call(start=c.start + (1 if c.start >= deleted_at else 0),
                 end=c.end + (1 if c.end >= deleted_at else 0),
                 strand=c.strand, partial=c.partial)
            for c in raw)
        return tuple(sorted(lifted, key=lambda c: c.key))

    # ------------------------------------------------------- the target rule

    def targets(self, baseline: Sequence[Call]) -> dict[str, Call]:
        """The declared target set, keyed by the declared length it stands for.

        For each declared length the campaign takes the baseline call whose own
        length is closest to it, ties broken by the lower start so the choice
        cannot depend on iteration order. A length that selects a call another
        length already took is a collision, and the rules are wrong rather than
        the data: the campaign refuses rather than silently auditing one gene
        twice.
        """
        out: dict[str, Call] = {}
        taken: set[tuple[int, int, str]] = set()
        for want in self.rules.target_lengths:
            ranked = sorted(baseline, key=lambda c: (abs(c.length - want), c.start))
            chosen = next((c for c in ranked if c.key not in taken), None)
            if chosen is None:
                raise FrameshiftError(
                    f"no call is left to stand for the declared length {want}")
            # Both directions. A one-sided test on the difference cannot catch
            # a declared length far LARGER than anything on the contig -- the
            # gap to a short call is bounded by the length itself -- so a ladder
            # asking for 20 kb on a contig whose longest gene is 300 bp would
            # silently audit that gene as its 20 kb target.
            if not want / 2 <= chosen.length <= 2 * want:
                raise FrameshiftError(
                    f"the closest call to the declared length {want} is "
                    f"{chosen.length} bp, more than a factor of two away in "
                    "either direction; the declared lengths do not describe "
                    "this contig")
            out[str(want)] = chosen
            taken.add(chosen.key)
        return out

    def indel_position(self, target: Call, fraction: float | None = None) -> int:
        """The 1-based position this campaign deletes inside a target.

        Taken as a FRACTION of the gene rather than a fixed distance, so the rule
        reads the same for a 100 bp gene and a 22 kb one. Rounded down, and
        clamped inside the gene so a fraction at either extreme cannot land on
        the flanking base.
        """
        chosen = self.rules.indel_fraction if fraction is None else fraction
        position = target.start + int((target.length - 1) * chosen)
        return min(max(position, target.start), target.end)

    # ---------------------------------------------------------- the response

    def response(self, baseline: Sequence[Call], target: Call,
                 mutated: Sequence[Call]) -> dict:
        """What the caller did to one target, and to everything else.

        The verdict counts the NEW calls overlapping the target's own span. By
        overlap rather than inside a window, because a frameshifted gene can
        produce an ORF running past where the original ended and any fixed window
        would misclassify it.

        `other_calls_lost` is the collateral: a single deleted base can cost a
        neighbouring call too, which is measured rather than assumed away.
        """
        before = {c.key for c in baseline}
        after = {c.key for c in mutated}
        survived = target.key in after
        new = [c for c in sorted(mutated, key=lambda c: c.key)
               if c.key not in before and c.overlaps(target.start, target.end)]
        if survived:
            verdict = "unchanged"
        elif len(new) >= 2:
            verdict = "split_into_two"
        elif len(new) == 1:
            verdict = "truncated_to_one"
        else:
            verdict = "abolished"
        if verdict not in self.rules.verdicts:
            raise FrameshiftError(
                f"{verdict!r} is not a declared verdict; the vocabulary and the "
                "classification disagree")
        covered = sum(min(target.end, c.end) - max(target.start, c.start) + 1
                      for c in new)
        return {
            "target_start": target.start, "target_end": target.end,
            "target_strand": target.strand, "target_length": target.length,
            "verdict": verdict,
            "overlapping_new_calls": len(new),
            "new_call_lengths": [c.length for c in new],
            "bases_of_the_target_still_called": covered,
            "calls_before": len(baseline), "calls_after": len(mutated),
            "other_calls_lost": len(before - after) - (0 if survived else 1),
            "other_calls_gained": len(after - before) - len(new),
        }

    def responses(self, baseline: Sequence[Call], targets: Mapping[str, Call],
                  mode: str | None = None) -> dict[str, dict]:
        chosen = mode or self.rules.baseline_mode
        whole = baseline if chosen == self.rules.baseline_mode else self.calls(chosen)
        out = {}
        for label, target in sorted(targets.items(), key=lambda kv: int(kv[0])):
            position = self.indel_position(target)
            mutated = self.calls(chosen, position)
            out[label] = {"deleted_at": position,
                          **self.response(whole, target, mutated)}
        return out

    def verdict_census(self, responses: Mapping[str, Mapping]) -> dict[str, int]:
        counts = collections.Counter(row["verdict"] for row in responses.values())
        return {v: counts.get(v, 0) for v in self.rules.verdicts}

    def verdict_by_length(self, responses: Mapping[str, Mapping]) -> dict[str, str]:
        return {label: row["verdict"]
                for label, row in sorted(responses.items(), key=lambda kv: int(kv[0]))}

    def verdict_is_monotone_in_length(self, responses: Mapping[str, Mapping]) -> bool:
        """Whether each verdict occupies one contiguous block as the gene gets
        longer, so that the per-length list is recoverable from its cut points.

        Not "changes at most once": with three verdicts on the ladder a monotone
        list changes twice, at two cut points. The test is one block per verdict,
        which is what makes the list compressible and R3 a threshold.

        Measured rather than assumed, and on this chromosome it is FALSE: a 189 bp
        gene survives as one truncated ORF while 219 and 261 bp genes are
        abolished, because whether a half is callable depends on a start codon
        falling in the right place and not only on how much room there is.
        """
        ordered = [row["verdict"] for _, row in
                   sorted(responses.items(), key=lambda kv: int(kv[0]))]
        changes = sum(1 for a, b in zip(ordered, ordered[1:]) if a != b)
        return changes <= len(set(ordered)) - 1

    # ------------------------------------------------------- the counterfactual

    def offset_sweep(self, baseline: Sequence[Call],
                     targets: Mapping[str, Call]) -> dict[str, dict]:
        """One target, the indel walked across it.

        The catalogued R4 was "an indel is introduced elsewhere -- predict the new
        call". Measured, the VERDICT does not change as the indel moves: what
        changes is where the gene splits, and it tracks the indel. So the sweep
        grades the split point rather than the verdict, which is the claim the
        data supports.
        """
        label = str(self.rules.offset_target_length)
        target = targets[label]
        out = {}
        for fraction in self.rules.offset_fractions:
            position = self.indel_position(target, fraction)
            mutated = self.calls(self.rules.baseline_mode, position)
            block = self.response(baseline, target, mutated)
            out[f"{fraction:g}"] = {
                "deleted_at": position,
                "bases_into_the_gene": position - target.start + 1,
                "verdict": block["verdict"],
                "new_call_lengths": block["new_call_lengths"],
                "bases_of_the_target_still_called":
                    block["bases_of_the_target_still_called"],
            }
        return out

    def noise_arm(self, targets: Mapping[str, Call]) -> dict:
        """The same deletions under the mode that retrains.

        The campaign's reason for pinning the mode it does, as a measurement: the
        signal is one call lost and two gained, and in the retraining mode a
        single deleted base moves an order of magnitude more calls than that
        across the rest of the chromosome. The signature is still there and it is
        no longer separable.
        """
        mode = self.rules.noise_mode
        whole = self.calls(mode)
        out = {}
        for want in self.rules.noise_target_lengths:
            label = str(want)
            target = targets[label]
            position = self.indel_position(target)
            mutated = self.calls(mode, position)
            block = self.response(whole, target, mutated)
            signal = block["overlapping_new_calls"] + (0 if block["verdict"] == "unchanged" else 1)
            noise = block["other_calls_lost"] + block["other_calls_gained"]
            out[label] = {
                "deleted_at": position,
                "verdict": block["verdict"],
                "calls_explained_by_the_target": signal,
                "other_calls_moved": noise,
                "noise_to_signal": (round(noise / signal, PRECISION)
                                    if signal else None),
            }
        return out
