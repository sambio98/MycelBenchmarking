"""Annotation transfer by top hit, and the rate at which it is wrong.

The oldest shortcut in functional genomics is to take a protein's closest
database hit and copy its annotation across. This engine does exactly that over a
shipped, closed reference set and then checks the result against the annotation
the query already carries. Gold is the arithmetic over a pinned search; no rung
asks what a protein "really" does.

Four things are load-bearing, and each one changes the answer.

**The reference set is closed and shipped.** Transfer against a live database
would make gold a function of the day it was computed. Here the queries and the
reference are the same shipped file, the self hit is excluded by a declared rule,
and the search is a pinned invocation of a pinned tool.

**The top hit needs a declared tie-break.** On this set 602 of 4,552 queries have
two or more hits sharing the top bitscore, so "the closest hit" is not a function
of the data until the tie-break is published. DIAMOND's own output order happens
to be stable for this version, but resting gold on an undeclared property of a
binary is how a campaign stops being reproducible across a tool upgrade.

**An EC number that stops early makes no claim.** 1,850 of this set's EC mentions
end in a dash, and an entry annotated `1.14.-.-` cannot agree or disagree at the
third level. Comparison at level k is therefore DECIDABLE only when both sides
specify level k, and the undecidable cases are counted rather than scored either
way. Treating a dash as a value is the quiet way to a wrong headline.

**Agreement is a function of identity, not a constant.** That is the finding the
ladder is built to expose, so identity bands are declared and graded.

Pure stdlib outside the search call, sorted order throughout.
"""

from __future__ import annotations

import collections
import gzip
import json
import os
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.tools import cache

TEMPLATE_ID = "annotation-transfer"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "transfer_rules.json"
PRECISION = 6

DEFAULT_PREFIX = "/opt/npbench-tools/bin"

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 1800

#: The six columns the declared invocation asks for, in order.
OUTFMT = ("qseqid", "sseqid", "pident", "length", "evalue", "bitscore")


class TransferError(RuntimeError):
    """The set or the environment does not support the transfer experiment."""


def tool_prefix() -> pathlib.Path:
    return pathlib.Path(os.environ.get("NPBENCH_TOOL_PREFIX", DEFAULT_PREFIX))


def _run(argv: Sequence[str]) -> subprocess.CompletedProcess:
    """Run a pinned tool with the prefix appended to PATH.

    Last, so a tool environment carrying its own interpreter cannot shadow the
    one this code runs under; present at all, because a tool that shells out
    resolves its helpers there.
    """
    prefix = str(tool_prefix())
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if prefix in parts:
        parts.remove(prefix)
    env = {**os.environ, **FIXED_ENV, "PATH": os.pathsep.join([*parts, prefix])}
    return subprocess.run(list(argv), capture_output=True, text=True, env=env,
                          timeout=TIMEOUT_S)


# ----------------------------------------------------------------- EC levels


def ec_level(code: str, depth: int) -> str | None:
    """The first `depth` levels of an EC number, or None if it stops earlier.

    `1.14.-.-` returns a value at depth 2 and None at depth 3. A code that is not
    four dot-separated fields returns None at every depth: this refuses to guess
    what a malformed annotation meant.
    """
    parts = code.split(".")
    if len(parts) != 4 or any(p == "-" or not p for p in parts[:depth]):
        return None
    return ".".join(parts[:depth])


def ec_agreement(query: Sequence[str], subject: Sequence[str],
                 depth: int) -> bool | None:
    """Does the subject's EC match the query's at `depth`?

    None means undecidable: at least one side does not specify that level, so
    there is nothing to agree or disagree about. An entry may carry several EC
    numbers, and the declared rule is any-vs-any -- a transfer is counted correct
    if any of the subject's codes matches any of the query's, because a
    multifunctional enzyme's annotation is a set and not a ranking.
    """
    left = {ec_level(c, depth) for c in query} - {None}
    right = {ec_level(c, depth) for c in subject} - {None}
    if not left or not right:
        return None
    return bool(left & right)


# ----------------------------------------------------------------- the rules


@dataclass(frozen=True)
class Rules:
    ec_prefix: str
    max_target_seqs: int
    search_flags: tuple[str, ...]
    graded_depths: tuple[int, ...]
    constant_depths: tuple[int, ...]
    verdict_order: tuple[str, ...]
    identity_bands: tuple[tuple[float, float], ...]
    tie_break: str
    removal_variants: Mapping[str, int]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        transfer = spec["transfer"]
        return cls(
            ec_prefix=transfer["ec_prefix"],
            max_target_seqs=int(transfer["max_target_seqs"]),
            search_flags=tuple(transfer["search_flags"]),
            graded_depths=tuple(int(d) for d in transfer["graded_depths"]),
            constant_depths=tuple(int(d) for d in transfer["constant_depths"]),
            verdict_order=tuple(transfer["verdict_order"]),
            identity_bands=tuple((float(lo), float(hi))
                                 for lo, hi in transfer["identity_bands"]),
            tie_break=transfer["tie_break"],
            removal_variants={
                name: int(spec["remove_top_hits"])
                for name, spec in sorted(transfer["removal_variants"].items())},
        )

    def band_label(self, lo: float, hi: float) -> str:
        return f"{lo:g}-{hi:g}"


# ------------------------------------------------------------------ the set


@dataclass(frozen=True)
class Hit:
    subject: str
    identity: float
    bitscore: float


@dataclass(frozen=True)
class Transfer:
    release: str
    sequences: Mapping[str, int]              # accession -> length
    annotations: Mapping[str, tuple[str, ...]]  # accession -> EC codes
    rules: Rules
    fasta_path: pathlib.Path

    @classmethod
    def load(cls, campaign: pathlib.Path, fasta: str, annotations: str) -> "Transfer":
        campaign = pathlib.Path(campaign)
        provenance = json.loads(
            (campaign / "inputs" / "provenance.json").read_text())
        path = campaign / fasta
        lengths: dict[str, int] = {}
        accession, buffer = None, []
        with (gzip.open(path, "rt") if path.suffix == ".gz" else open(path)) as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if accession:
                        lengths[accession] = len("".join(buffer))
                    accession, buffer = line[1:].split()[0], []
                else:
                    buffer.append(line.strip())
        if accession:
            lengths[accession] = len("".join(buffer))
        if not lengths:
            raise TransferError("the shipped FASTA holds no sequences")

        apath = campaign / annotations
        with (gzip.open(apath, "rt") if apath.suffix == ".gz" else open(apath)) as h:
            raw = json.load(h)
        annotated = {k: tuple(v) for k, v in sorted(raw["ec_numbers"].items())}
        missing = set(lengths) - set(annotated)
        if missing:
            raise TransferError(
                f"{len(missing)} sequences carry no annotation record, so a "
                f"transfer verdict for them is undefined: {sorted(missing)[:5]}")
        return cls(release=provenance["release"], sequences=dict(sorted(lengths.items())),
                   annotations=annotated, rules=Rules.load(campaign),
                   fasta_path=path)

    # ------------------------------------------------------------ the search

    def diamond_version(self) -> str:
        done = _run([str(tool_prefix() / "diamond"), "version"])
        if done.returncode != 0:
            raise TransferError(f"diamond version failed: {done.stderr[-400:]}")
        return done.stdout.strip()

    def _search(self) -> dict[str, tuple[Hit, ...]]:
        """One all-against-all search of the shipped set against itself."""
        with tempfile.TemporaryDirectory(prefix="npbench-transfer-") as tmp:
            work = pathlib.Path(tmp)
            plain = work / "set.fasta"
            with (gzip.open(self.fasta_path, "rt") if self.fasta_path.suffix == ".gz"
                  else open(self.fasta_path)) as src, open(plain, "w") as dst:
                for line in src:
                    dst.write(line)
            database = work / "set"
            made = _run([str(tool_prefix() / "diamond"), "makedb", "--in", str(plain),
                         "-d", str(database), "--threads", "2", "--quiet"])
            if made.returncode != 0:
                raise TransferError(f"diamond makedb failed: {made.stderr[-600:]}")
            table = work / "hits.tsv"
            done = _run([str(tool_prefix() / "diamond"), "blastp",
                         "-q", str(plain), "-d", str(database), "-o", str(table),
                         "--threads", "2",
                         "--max-target-seqs", str(self.rules.max_target_seqs),
                         *self.rules.search_flags,
                         "--outfmt", "6", *OUTFMT])
            if done.returncode != 0:
                raise TransferError(f"diamond blastp failed: {done.stderr[-600:]}")
            raw: dict[str, list[Hit]] = {}
            for line in table.read_text().splitlines():
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) != len(OUTFMT):
                    raise TransferError(
                        f"expected {len(OUTFMT)} columns from the declared "
                        f"outfmt, got {len(fields)}: {line[:80]!r}")
                query, subject = fields[0], fields[1]
                # The self hit is excluded here rather than filtered later, so no
                # downstream rule can accidentally transfer an annotation from the
                # query to itself and call it support.
                if query == subject:
                    continue
                raw.setdefault(query, []).append(
                    Hit(subject=subject, identity=float(fields[2]),
                        bitscore=float(fields[5])))
            unknown = ({q for q in raw} | {h.subject for v in raw.values() for h in v}
                       ) - set(self.sequences)
            if unknown:
                raise TransferError(
                    f"the search reported identifiers outside the shipped set: "
                    f"{sorted(unknown)[:5]}")
            return {q: tuple(v) for q, v in sorted(raw.items())}

    def hits(self) -> dict[str, tuple[Hit, ...]]:
        """Ranked hits per query, under the declared tie-break.

        Memoised only where a tool memo is enabled, which is stub sweeps alone.
        """
        def compute() -> dict[str, list[list]]:
            return {q: [[h.subject, h.identity, h.bitscore] for h in self.ranked(v)]
                    for q, v in self._search().items()}

        if not cache.enabled():
            return {q: self.ranked(v) for q, v in self._search().items()}
        memo = cache.memoize_json(
            ["annotation_transfer.hits",
             cache.file_digest(self.fasta_path),
             self.diamond_version(),
             str(self.rules.max_target_seqs),
             ",".join(self.rules.search_flags) or "none",
             self.rules.tie_break],
            compute)
        return {q: tuple(Hit(s, float(i), float(b)) for s, i, b in v)
                for q, v in memo.items()}

    def ranked(self, hits: Sequence[Hit]) -> tuple[Hit, ...]:
        """Declared order: bitscore descending, then subject accession ascending.

        The second key is not decoration. 602 of this set's queries share their
        top bitscore with another hit, so without it the top hit is whatever the
        binary happened to emit first.
        """
        return tuple(sorted(hits, key=lambda h: (-h.bitscore, h.subject)))

    # ----------------------------------------------------------- the verdict

    def verdict(self, query: str, hits: Sequence[Hit],
                excluded: frozenset[str] = frozenset()) -> tuple[str, str | None]:
        """The declared verdict for transferring from this query's closest hit.

        First match wins over the declared order, which is why the order is
        published: a query that is undecidable at the third level must not also be
        counted as an error at the fourth.
        """
        pool = [h for h in hits if h.subject not in excluded]
        if not pool:
            return "no_hit", None
        subject = pool[0].subject
        mine, theirs = self.annotations[query], self.annotations[subject]
        at_three = ec_agreement(mine, theirs, 3)
        if at_three is None:
            return "undecidable", subject
        if at_three is False:
            return "wrong_at_subsubclass", subject
        at_four = ec_agreement(mine, theirs, 4)
        if at_four is None:
            return "supported_to_subsubclass", subject
        return ("supported" if at_four else "wrong_at_serial"), subject

    def verdicts(self, hits: Mapping[str, Sequence[Hit]],
                 excluded: Mapping[str, frozenset[str]] | None = None) -> dict[str, str]:
        out = {}
        for query in sorted(hits):
            skip = (excluded or {}).get(query, frozenset())
            out[query] = self.verdict(query, hits[query], skip)[0]
        return out

    def verdict_census(self, verdicts: Mapping[str, str]) -> dict[str, int]:
        counts = collections.Counter(verdicts.values())
        return {v: counts.get(v, 0) for v in self.rules.verdict_order}

    def census_after_removal(self, hits: Mapping[str, Sequence[Hit]],
                             remove: int) -> dict[str, int]:
        """The verdict census when the top `remove` hits are taken out.

        The catalogued counterfactual, and it needs no second search: removing a
        hit from the reference set means reading further down a ranking that is
        already computed. remove=0 is the identity control and must reproduce the
        headline census exactly.
        """
        excluded = {q: frozenset(h.subject for h in self.ranked(v)[:remove])
                    for q, v in hits.items()}
        return self.verdict_census(self.verdicts(hits, excluded))

    def verdict_transitions(self, hits: Mapping[str, Sequence[Hit]],
                            remove: int) -> dict[str, int]:
        """How many queries move from one verdict to another under removal.

        The count of CHANGED verdicts is what makes the rung informative: a
        removal that leaves every verdict alone has demonstrated nothing.
        """
        before = self.verdicts(hits)
        excluded = {q: frozenset(h.subject for h in self.ranked(v)[:remove])
                    for q, v in hits.items()}
        after = self.verdicts(hits, excluded)
        moves = collections.Counter(
            f"{before[q]}->{after[q]}" for q in sorted(before)
            if before[q] != after[q])
        return dict(sorted(moves.items(), key=lambda kv: (-kv[1], kv[0])))

    # ------------------------------------------------------------ agreement

    def agreement_by_depth(self, hits: Mapping[str, Sequence[Hit]],
                           depths: Sequence[int]) -> dict[str, dict]:
        """Decidable count and agreement rate at each requested EC depth."""
        out = {}
        for depth in depths:
            decidable = agree = 0
            for query in sorted(hits):
                ranked = self.ranked(hits[query])
                if not ranked:
                    continue
                verdict = ec_agreement(self.annotations[query],
                                       self.annotations[ranked[0].subject], depth)
                if verdict is None:
                    continue
                decidable += 1
                agree += verdict
            out[str(depth)] = {
                "decidable": decidable,
                "agree": agree,
                "rate": round(agree / decidable, PRECISION) if decidable else None,
            }
        return out

    def agreement_by_identity(self, hits: Mapping[str, Sequence[Hit]],
                              depth: int) -> dict[str, dict]:
        """The curve the campaign exists to expose: agreement against identity."""
        out = {}
        for lo, hi in self.rules.identity_bands:
            decidable = agree = 0
            for query in sorted(hits):
                ranked = self.ranked(hits[query])
                if not ranked or not (lo <= ranked[0].identity < hi):
                    continue
                verdict = ec_agreement(self.annotations[query],
                                       self.annotations[ranked[0].subject], depth)
                if verdict is None:
                    continue
                decidable += 1
                agree += verdict
            out[self.rules.band_label(lo, hi)] = {
                "decidable": decidable,
                "agree": agree,
                "rate": round(agree / decidable, PRECISION) if decidable else None,
            }
        return out

    # -------------------------------------------------------------- the set

    def annotation_census(self) -> dict:
        codes = collections.Counter()
        per_entry = collections.Counter()
        partial = 0
        for accession in sorted(self.annotations):
            mine = self.annotations[accession]
            per_entry[len(mine)] += 1
            for code in mine:
                codes[code] += 1
                if ec_level(code, 4) is None:
                    partial += 1
        return {
            "distinct_ec_codes": len(codes),
            "ec_codes_per_entry": {str(k): v for k, v in sorted(per_entry.items())},
            "partial_ec_mentions": partial,
            "entries_with_multiple_ec": sum(v for k, v in per_entry.items() if k > 1),
        }

    def tied_top_queries(self, hits: Mapping[str, Sequence[Hit]]) -> int:
        """Queries whose top bitscore is shared, so the tie-break decides."""
        tied = 0
        for query in sorted(hits):
            ranked = self.ranked(hits[query])
            if len(ranked) > 1 and ranked[0].bitscore == ranked[1].bitscore:
                tied += 1
        return tied

    def identity_census(self, hits: Mapping[str, Sequence[Hit]]) -> dict[str, int]:
        counts: dict[str, int] = {}
        for lo, hi in self.rules.identity_bands:
            label = self.rules.band_label(lo, hi)
            counts[label] = sum(
                1 for q in hits
                if (r := self.ranked(hits[q])) and lo <= r[0].identity < hi)
        return counts
