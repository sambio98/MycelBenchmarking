"""What a pairwise search can see of a relationship a curated resource asserts.

Remote homology is the standard example of a question with no oracle. At thirty
per cent identity two proteins may share an ancestor or a fold, and asking a
system "are these homologous" is asking it to agree with whoever curated the
answer. This campaign does not ask that. It takes InterPro's own published
assertion -- two proteins carry the same `family` entry, or only the same
`homologous_superfamily` entry -- and grades what a pinned pairwise search
REPORTS against it, by identity band.

So the claim is never about homology. It is about the search: where does a
pairwise searcher stop recovering a relationship a profile-based resource
asserts, and what does it report instead.

Four things are load-bearing.

**The label is a published curation with a release, not a judgement.** InterPro
110.0 says which entries a protein carries and what type each entry is. The
campaign ships that table because it cannot be posed without it, and the table is
not the answer: the answer is what the search finds, which the table does not
contain. And the vocabulary refuses the obvious overreach -- a pair sharing no
entry is NOT labelled non-homologous, because InterPro records what has been
classified rather than what is related.

**Every pair in the corpus is already a pair InterPro places in one superfamily.**
The corpus is every reviewed entry cross-referencing IPR036188, so the question
is never whether two members are related at all; it is whether the search can see
a relationship that is asserted at family level, and what it reports for the rest.

**The twilight zone is measured, not assumed.** On this corpus at the declared
sensitivity, hits between 25 and 35 per cent identity hold same-family pairs and
superfamily-only pairs in comparable numbers, which is what makes that band the
one where identity alone cannot separate them. Above 35 per cent the same-family
pairs dominate; the campaign reports the census per band rather than asserting a
cutoff.

**And the counterfactual is a window ladder, because the catalogued one is
degenerate.** The catalogue asked for the region carrying the shared motif to be
shuffled. Measured, an InterPro family match on these proteins covers almost the
whole sequence, so that perturbation destroys the protein and grades one answer.
A declared ladder of window lengths centred on the query's midpoint asks how much
local disorder a hit survives instead -- and the measured answer is ordered by
identity: the distant hits die first, at a higher rate in every row where
anything is lost.

Everything outside the search runs is stdlib and sorted order throughout.
"""

from __future__ import annotations

import collections
import gzip
import hashlib
import json
import os
import pathlib
import random
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.tools import cache
from npbench_c.tools.resources import tool_prefix

TEMPLATE_ID = "remote-homology"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "homology_rules.json"
PRECISION = 6

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 3600

#: The declared sequence alphabet. Anything else in a shipped corpus means the
#: input is not what the campaign thinks it is.
ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZUO")


class HomologyError(RuntimeError):
    """The corpus, the labels or the environment does not support the campaign."""


def _run(argv: Sequence[str]) -> subprocess.CompletedProcess:
    env = {**os.environ, **FIXED_ENV}
    return subprocess.run(list(argv), capture_output=True, text=True, env=env,
                          timeout=TIMEOUT_S)


# ------------------------------------------------------------------- the rules


@dataclass(frozen=True)
class Rules:
    sensitivity: float
    evalue: float
    max_seqs: int
    threads: int
    format_fields: tuple[str, ...]
    family_types: tuple[str, ...]
    superfamily_types: tuple[str, ...]
    corpus_superfamily: str
    identity_bands: tuple[tuple[float, float], ...]
    band_labels: tuple[str, ...]
    queries: int
    shuffle_windows: tuple[int, ...]
    verdicts: tuple[str, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        a = spec["homology"]
        bands = tuple(tuple(float(x) for x in band) for band in a["identity_bands"])
        labels = tuple(a["band_labels"])
        if len(bands) != len(labels):
            raise HomologyError(
                f"{len(bands)} identity bands and {len(labels)} band labels do not "
                "describe the same partition")
        windows = tuple(int(w) for w in a["shuffle_windows"])
        if 0 not in windows:
            raise HomologyError(
                "the shuffle ladder has no zero-length window, so the "
                "counterfactual has no identity control to compare against")
        if len(set(windows)) != len(windows):
            raise HomologyError(f"the declared shuffle windows repeat: {windows}")
        if sorted(windows) != list(windows):
            raise HomologyError(
                f"the declared shuffle windows are not in increasing order "
                f"({windows}), so the ladder does not read as a ladder")
        if int(a["queries"]) < 3:
            raise HomologyError(
                f"{a['queries']} queries cannot carry a per-query census")
        if int(a["max_seqs"]) < 1:
            raise HomologyError("max_seqs must be positive")
        if not a["family_types"] or not a["superfamily_types"]:
            raise HomologyError(
                "the family and superfamily entry types are what the label IS; "
                "neither may be empty")
        if set(a["family_types"]) & set(a["superfamily_types"]):
            raise HomologyError(
                f"the entry types {sorted(set(a['family_types']) & set(a['superfamily_types']))} "
                "are declared as both family and superfamily, so one entry would "
                "make a pair same-family and superfamily-only at once")
        return cls(
            sensitivity=float(a["sensitivity"]), evalue=float(a["evalue"]),
            max_seqs=int(a["max_seqs"]), threads=int(a["threads"]),
            format_fields=tuple(a["format_fields"]),
            family_types=tuple(a["family_types"]),
            superfamily_types=tuple(a["superfamily_types"]),
            corpus_superfamily=str(a["corpus_superfamily"]),
            identity_bands=bands, band_labels=labels,
            queries=int(a["queries"]), shuffle_windows=windows,
            verdicts=tuple(a["verdicts"]))

    def band(self, identity: float) -> str:
        """The declared band a hit's identity falls in.

        Half-open on the low side so a hit at exactly a boundary lands in the
        upper band, and the last band is the sentinel for everything above.
        """
        for (lo, hi), label in zip(self.identity_bands[:-1], self.band_labels[:-1]):
            if lo <= identity < hi:
                return label
        return self.band_labels[-1]


# ------------------------------------------------------------------ the corpus


@dataclass(frozen=True)
class Hit:
    query: str
    target: str
    identity: float
    alnlen: int
    evalue: float
    bits: float

    @property
    def key(self) -> tuple[str, str]:
        return (self.query, self.target)


def read_corpus(path: pathlib.Path) -> dict[str, str]:
    """A gzipped FASTA as {accession: sequence}, under the declared alphabet."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        text = fh.read()
    out: dict[str, str] = {}
    accession = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if accession is not None:
                out[accession] = "".join(chunks)
            accession = line[1:].split()[0]
            if accession in out:
                raise HomologyError(
                    f"the accession {accession} appears twice in the corpus, so a "
                    "pair could be a protein against itself under two names")
            chunks = []
        elif accession is not None:
            chunks.append(line.strip().upper())
    if accession is not None:
        out[accession] = "".join(chunks)
    if not out:
        raise HomologyError("the corpus holds no sequence")
    for accession, sequence in sorted(out.items()):
        stray = set(sequence) - ALPHABET
        if stray:
            raise HomologyError(
                f"{accession} carries {sorted(stray)}, which is outside the "
                "declared amino-acid alphabet")
    return out


@dataclass(frozen=True)
class Corpus:
    rules: Rules
    sequences: Mapping[str, str]
    families: Mapping[str, frozenset]
    superfamilies: Mapping[str, frozenset]
    entry_names: Mapping[str, str]

    @classmethod
    def load(cls, campaign: pathlib.Path, fasta: str, labels: str) -> "Corpus":
        campaign = pathlib.Path(campaign)
        rules = Rules.load(campaign)
        sequences = read_corpus(campaign / fasta)
        with gzip.open(campaign / labels, "rt") as fh:
            table = json.load(fh)
        types = table["entry_types"]
        assignments = table["assignments"]

        missing = sorted(set(sequences) - set(assignments))
        if missing:
            raise HomologyError(
                f"{len(missing)} corpus accessions have no InterPro assignment "
                f"(first: {missing[:3]}), so the label does not cover the corpus")

        fam, sup = {}, {}
        for accession in sorted(sequences):
            entries = assignments[accession]
            unknown = [e for e in entries if e not in types]
            if unknown:
                raise HomologyError(
                    f"{accession} carries entries with no declared type "
                    f"({unknown[:3]}), so their label class is undefined")
            fam[accession] = frozenset(
                e for e in entries if types[e]["type"] in rules.family_types)
            sup[accession] = frozenset(
                e for e in entries if types[e]["type"] in rules.superfamily_types)

        off = sorted(a for a, s in sup.items()
                     if rules.corpus_superfamily not in s)
        if off:
            raise HomologyError(
                f"{len(off)} corpus members do not carry the declared superfamily "
                f"{rules.corpus_superfamily} (first: {off[:3]}); the corpus is not "
                "the set the selection rule describes")
        return cls(rules=rules, sequences=sequences, families=fam,
                   superfamilies=sup,
                   entry_names={e: meta["name"] for e, meta in sorted(types.items())})

    # ------------------------------------------------------------ the inventory

    @property
    def size(self) -> int:
        return len(self.sequences)

    def sequences_sha256(self) -> str:
        """A digest over the SEQUENCES, in accession order, not over the file.

        The same corpus rewrapped at another line width, or recompressed, must
        fingerprint the same, which is what makes this a parse check rather than
        a checksum of a download.
        """
        digest = hashlib.sha256()
        for accession in sorted(self.sequences):
            digest.update(accession.encode())
            digest.update(b"\n")
            digest.update(self.sequences[accession].encode())
            digest.update(b"\n")
        return digest.hexdigest()

    def family_census(self) -> dict[str, int]:
        """How many corpus members each declared family entry covers."""
        counts: collections.Counter = collections.Counter()
        for entries in self.families.values():
            for entry in entries:
                counts[entry] += 1
        return dict(sorted(counts.items()))

    # --------------------------------------------------------------- the tool

    def mmseqs_version(self) -> str:
        done = _run([str(tool_prefix() / "mmseqs"), "version"])
        text = (done.stdout or "") + (done.stderr or "")
        for line in text.splitlines():
            if line.strip():
                return line.strip()
        raise HomologyError("mmseqs did not report a version")

    def _flags(self) -> tuple[str, ...]:
        return ("--threads", str(self.rules.threads),
                "--max-seqs", str(self.rules.max_seqs),
                "-s", f"{self.rules.sensitivity:g}",
                "-e", f"{self.rules.evalue:g}",
                "--format-output", ",".join(self.rules.format_fields),
                "-v", "0")

    def _write_fasta(self, path: pathlib.Path, batch: Mapping[str, str]) -> None:
        path.write_text("".join(
            f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
            + "\n" for name, seq in sorted(batch.items())))

    def _search_uncached(self, batch: Mapping[str, str]) -> tuple[Hit, ...]:
        prefix = tool_prefix() / "mmseqs"
        with tempfile.TemporaryDirectory(prefix="npbench-homology-") as tmp:
            work = pathlib.Path(tmp)
            target = work / "corpus.faa"
            self._write_fasta(target, self.sequences)
            db = work / "db"
            done = _run([str(prefix), "createdb", str(target), str(db), "-v", "0"])
            if done.returncode != 0:
                raise HomologyError(
                    f"mmseqs createdb failed: {(done.stderr or '')[-500:]}")
            done = _run([str(prefix), "createindex", str(db), str(work / "idx"),
                         "--threads", str(self.rules.threads), "-v", "0"])
            if done.returncode != 0:
                raise HomologyError(
                    f"mmseqs createindex failed: {(done.stderr or '')[-500:]}")
            queries = work / "queries.faa"
            self._write_fasta(queries, batch)
            out = work / "hits.tsv"
            done = _run([str(prefix), "easy-search", str(queries), str(db),
                         str(out), str(work / "tmp"), *self._flags()])
            if done.returncode != 0:
                raise HomologyError(
                    f"mmseqs easy-search failed: {(done.stderr or '')[-500:]}")
            # The declared canonicalisation: sort by line before anything reads
            # it. The output has no header and each line is one record.
            lines = sorted(out.read_text().splitlines())
        hits = []
        for line in lines:
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) != len(self.rules.format_fields):
                raise HomologyError(
                    f"a hit line carries {len(fields)} fields and the declared "
                    f"projection has {len(self.rules.format_fields)}")
            query, target, fident, alnlen, evalue, bits = fields
            hits.append(Hit(query=query, target=target, identity=float(fident),
                            alnlen=int(alnlen), evalue=float(evalue),
                            bits=float(bits)))
        if not hits:
            raise HomologyError(
                "the search reported no hit at all, not even a self-hit, which is "
                "a broken run rather than a finding about remote homology")
        return tuple(sorted(hits, key=lambda h: (h.query, h.target)))

    def search(self, batch: Mapping[str, str]) -> tuple[Hit, ...]:
        """The pinned search for one batch of queries, in declared order.

        Memoised when and only when a tool memo is enabled, which the sweep does
        for stub systems alone.
        """
        if not cache.enabled():
            return self._search_uncached(batch)
        payload = cache.memoize_json(
            ["remote_homology.search", self.sequences_sha256(),
             self.mmseqs_version(), " ".join(self._flags()),
             hashlib.sha256(
                 json.dumps(dict(sorted(batch.items())), sort_keys=True).encode()
             ).hexdigest()],
            lambda: [h.__dict__ for h in self._search_uncached(batch)])
        return tuple(Hit(**row) for row in payload)

    # -------------------------------------------------------- the query rule

    def queries(self) -> dict[str, str]:
        """The declared query set, as {family entry: accession}.

        For each of the largest declared families, in descending size with ties
        broken by entry accession, the lexicographically first member not already
        taken. The taken set is load-bearing: one accession can be the first
        member of two families, and auditing it twice would make a family stand
        for nothing.
        """
        members: dict[str, list[str]] = collections.defaultdict(list)
        for accession, entries in sorted(self.families.items()):
            for entry in sorted(entries):
                members[entry].append(accession)
        out: dict[str, str] = {}
        taken: set[str] = set()
        for entry, group in sorted(members.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            pick = next((a for a in sorted(group) if a not in taken), None)
            if pick is None:
                continue
            out[entry] = pick
            taken.add(pick)
            if len(out) == self.rules.queries:
                break
        if len(out) < self.rules.queries:
            raise HomologyError(
                f"only {len(out)} families could supply a query that no other "
                f"family already took, and {self.rules.queries} are declared")
        return dict(sorted(out.items()))

    # ------------------------------------------------------- the label join

    def pair_label(self, left: str, right: str) -> str:
        """What the curated resource asserts about one pair.

        Three outcomes, and the third is deliberately not "not homologous":
        InterPro records what has been classified, which is not the same as what
        is related.
        """
        if self.families.get(left, frozenset()) & self.families.get(right, frozenset()):
            return "same_family"
        if self.superfamilies.get(left, frozenset()) & self.superfamilies.get(
                right, frozenset()):
            return "superfamily_only"
        return "unclassified_pair"

    def same_family_members(self, query: str) -> frozenset:
        """Every corpus member InterPro places in a family with this query."""
        wanted = self.families.get(query, frozenset())
        if not wanted:
            return frozenset()
        return frozenset(
            accession for accession, entries in self.families.items()
            if accession != query and entries & wanted)

    # ---------------------------------------------------------- the reconcile

    def hit_table(self, hits: Sequence[Hit], query: str) -> list[dict]:
        """The declared projection of one query's hits, in declared order."""
        return [
            {"target": h.target, "identity": round(h.identity, PRECISION),
             "alnlen": h.alnlen, "band": self.rules.band(h.identity),
             "label": self.pair_label(h.query, h.target)}
            for h in hits if h.query == query and h.target != query]

    def band_census(self, hits: Sequence[Hit], query: str) -> dict[str, dict]:
        """Per declared band, how many of this query's hits carry each label.

        Every band and every label is always present, so a census whose keys
        depend on the data cannot make the measurer's shape depend on the answer.
        """
        labels = ("same_family", "superfamily_only", "unclassified_pair")
        out = {band: {label: 0 for label in labels}
               for band in self.rules.band_labels}
        for hit in hits:
            if hit.query != query or hit.target == query:
                continue
            out[self.rules.band(hit.identity)][
                self.pair_label(hit.query, hit.target)] += 1
        return out

    def recovery(self, hits: Sequence[Hit], query: str) -> dict:
        """What the search found of what the resource asserts, and what it missed.

        `same_family_missed` is the campaign's point: a relationship InterPro
        asserts at family level that this pairwise search, at this sensitivity,
        does not report at all.
        """
        detected = {h.target for h in hits
                    if h.query == query and h.target != query}
        asserted = self.same_family_members(query)
        found = asserted & detected
        return {
            "hits": len(detected),
            "same_family_asserted": len(asserted),
            "same_family_detected": len(found),
            "same_family_missed": len(asserted - detected),
            "recovery_rate": (round(len(found) / len(asserted), PRECISION)
                              if asserted else None),
            "superfamily_only_detected": len(detected - asserted),
        }

    def verdict(self, cell: Mapping[str, int]) -> str:
        """The declared verdict for one band's census, computed and not judged."""
        same = cell["same_family"]
        only = cell["superfamily_only"]
        if same and only:
            chosen = "overlapping"
        elif same:
            chosen = "separable"
        elif only:
            chosen = "search_only"
        else:
            chosen = "unresolved"
        if chosen not in self.rules.verdicts:
            raise HomologyError(
                f"{chosen!r} is not a declared verdict; the vocabulary and the "
                "classification disagree")
        return chosen

    # ------------------------------------------------------- the counterfactual

    def shuffled(self, accession: str, window: int) -> str:
        """The query with a declared window of residues permuted in place.

        Seeded from the accession and the window, so the same window on two
        queries is not the same permutation and every run reproduces both.
        """
        sequence = self.sequences[accession]
        if window == 0:
            return sequence
        if window > len(sequence):
            raise HomologyError(
                f"the declared window {window} is longer than {accession} "
                f"({len(sequence)} residues), so the ladder's top rung is a whole "
                "sequence shuffle on this query and measures nothing local")
        seed = int(hashlib.sha256(f"{accession}:{window}".encode()).hexdigest(), 16)
        rnd = random.Random(seed)
        middle = len(sequence) // 2
        low = max(0, middle - window // 2)
        high = min(len(sequence), low + window)
        block = list(sequence[low:high])
        rnd.shuffle(block)
        return sequence[:low] + "".join(block) + sequence[high:]

    def ladder(self, baseline: Sequence[Hit]) -> dict[str, dict]:
        """The shuffle ladder: what each query keeps as the window grows.

        The rates are per identity band OF THE BASELINE HIT, which is the only
        normalisation that makes the bands comparable: a query with seven hundred
        hits above 35 per cent and five below 25 loses more of the former in
        absolute terms while losing all of the latter.
        """
        queries = self.queries()
        base = {
            accession: {h.target: h.identity for h in baseline
                        if h.query == accession and h.target != accession}
            for accession in queries.values()}
        sizes = {
            accession: collections.Counter(
                self.rules.band(identity) for identity in row.values())
            for accession, row in base.items()}

        out: dict[str, dict] = {}
        for window in self.rules.shuffle_windows:
            batch = {accession: self.shuffled(accession, window)
                     for accession in sorted(queries.values())}
            hits = self.search(batch)
            per_query = {}
            for accession in sorted(queries.values()):
                now = {h.target for h in hits
                       if h.query == accession and h.target != accession}
                lost = set(base[accession]) - now
                lost_by_band = collections.Counter(
                    self.rules.band(base[accession][t]) for t in lost)
                rates = {
                    band: (round(lost_by_band[band] / sizes[accession][band], PRECISION)
                           if sizes[accession][band] else None)
                    for band in self.rules.band_labels}
                per_query[accession] = {
                    "hits": len(now),
                    "lost": len(lost),
                    "lost_by_band": {band: lost_by_band[band]
                                     for band in self.rules.band_labels},
                    "loss_rate_by_band": rates,
                    "band_order_held": _ordered(
                        [rates[band] for band in self.rules.band_labels]),
                }
            out[str(window)] = dict(sorted(per_query.items()))
        return out

    def baseline_band_sizes(self, baseline: Sequence[Hit]) -> dict[str, dict]:
        """How many of each query's baseline hits fall in each band.

        Reported because every loss rate in the ladder is a fraction of one of
        these, and a rate without its denominator is not a measurement.
        """
        out = {}
        for accession in sorted(self.queries().values()):
            counts = collections.Counter(
                self.rules.band(h.identity) for h in baseline
                if h.query == accession and h.target != accession)
            out[accession] = {band: counts[band] for band in self.rules.band_labels}
        return dict(sorted(out.items()))


def _ordered(rates: Sequence[float | None]) -> bool:
    """Whether loss rates fall from the distant band to the close one.

    The campaign's R4 claim, as a boolean per cell: a local shuffle should cost
    the weakest hits first. Bands with no baseline hit carry no rate and are
    skipped rather than counted as zero, which would invent an ordering.
    """
    present = [r for r in rates if r is not None]
    return all(a >= b for a, b in zip(present, present[1:]))
