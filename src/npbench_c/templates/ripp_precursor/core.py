"""RiPP precursor annotation audit: where the core peptide sits, and who says so.

A ribosomally synthesised and post-translationally modified peptide is made from
a precursor protein that carries a leader, a core, and sometimes a follower. The
core is the part that becomes the product. MIBiG records, per precursor, the gene
that encodes it, the core peptide's sequence, and -- for some entries -- a
`leader_cleavage_location`. The precursor's own translation is in the
antiSMASH-processed reference record for the same cluster.

Three things can therefore be checked against each other, and none of the checks
needs a judgement: the declared core must be findable in the translation, the
declared cleavage coordinate must land where the core begins, and a domain model
that claims to recognise a precursor must have an envelope somewhere in it.

Four things are load-bearing.

**The field has four JSON shapes and the campaign does not repair them silently.**
`core_sequence` arrives as a bare string, as a list of one string, as a list of
several strings, and as a string holding the printed form of a list -- brackets
and quotes included. A reader that assumes any one shape mis-parses the others,
and a reader that repairs the printed form without saying so has made a curation
decision inside its parser. So the ingestion policy is declared, the baseline is
the conservative one, and the repair is an R4 counterfactual rather than a
default.

**Case is a normalisation, not a validity question.** Some cores are recorded in
lower case. Folding case is a declared transform; a character that is not an
amino-acid letter in any case is a malformed core. Conflating the two reports a
case difference as a data defect.

**The coordinate convention is established, not assumed.** Nothing in MIBiG says
whether `leader_cleavage_location` is 0-based or 1-based, or whether `from` or
`to` marks the core's first residue. All four readings are scored against the
core's measured position in the translation and the convention is ELECTED by
agreement, with the tie-break order declared. This is the same move the detection
campaign had to make for MIBiG's `from`/`to` locus bounds.

**A domain model is evidence about the precursor, not about the boundary.** The
declared Pfam panel is scanned at a pinned threshold and each hit's envelope is
compared with the measured leader/core boundary. Whether any family marks that
boundary is the measurement, and the leader-deletion counterfactual asks the tool
directly which side of the boundary its own signal lives on.

Everything outside the two HMMER calls is pure stdlib and sorted order throughout.
"""

from __future__ import annotations

import collections
import hashlib
import json
import os
import pathlib
import re
import subprocess
import tarfile
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.tools import cache
from npbench_c.tools.resources import (
    DEFAULT_PREFIX,
    Resource,
    load_declared,
    sha256,
    tool_prefix,
)

TEMPLATE_ID = "ripp-precursor-audit"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "precursor_rules.json"
PRECISION = 6

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 1800


class PrecursorError(RuntimeError):
    """The annotation set or the environment does not support the audit."""


def _run(argv: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a pinned tool with the prefix appended to PATH, last.

    First-or-not it must be on PATH at all: a tool that shells out to helpers
    resolves them there, and leaving it off measures the ambient shell instead of
    the image.
    """
    prefix = str(tool_prefix())
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if prefix in parts:
        parts.remove(prefix)
    env = {**os.environ, **FIXED_ENV, "PATH": os.pathsep.join([*parts, prefix])}
    kwargs.setdefault("stdout", subprocess.PIPE)
    kwargs.setdefault("stderr", subprocess.PIPE)
    return subprocess.run(list(argv), text=True, env=env, timeout=TIMEOUT_S,
                          **kwargs)


def _digest_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# ----------------------------------------------------------------------- rules


@dataclass(frozen=True)
class Policy:
    """A declared way of reading `core_sequence` into a list of cores."""

    name: str
    case_fold: bool
    repair_printed_list: bool


#: A string that is the printed form of a Python list of quoted strings. Matched
#: whole rather than searched for, so a core that merely CONTAINS a bracket is
#: still malformed rather than silently reinterpreted.
_PRINTED_LIST = re.compile(
    r"^\[\s*(?:'[^']*'|\"[^\"]*\")(?:\s*,\s*(?:'[^']*'|\"[^\"]*\"))*\s*\]$")


@dataclass(frozen=True)
class Rules:
    biosynthetic_class: str
    excluded_statuses: tuple[str, ...]
    alphabet: frozenset[str]
    policies: Mapping[str, Policy]
    baseline_policy: str
    shape_order: tuple[str, ...]
    case_order: tuple[str, ...]
    verdict_order: tuple[str, ...]
    conventions: tuple[str, ...]
    reconciliation_verdicts: tuple[str, ...]
    panel: tuple[str, ...]
    threshold: str
    relation_order: tuple[str, ...]
    deletion_outcomes: tuple[str, ...]
    resources: Mapping[str, Resource]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        audit = spec["audit"]
        policies = {
            name: Policy(name=name, case_fold=bool(block["case_fold"]),
                         repair_printed_list=bool(block["repair_printed_list"]))
            for name, block in sorted(audit["ingestion_policies"].items())}
        baseline = audit["baseline_policy"]
        if baseline not in policies:
            raise PrecursorError(
                f"the baseline policy {baseline!r} is not among the declared "
                f"policies {sorted(policies)}")
        alphabet = frozenset(audit["amino_acid_alphabet"])
        if len(alphabet) != 20:
            raise PrecursorError(
                f"the declared alphabet has {len(alphabet)} letters; a protein "
                "alphabet that is not the twenty is a different campaign")
        return cls(
            biosynthetic_class=audit["biosynthetic_class"],
            excluded_statuses=tuple(audit["excluded_statuses"]),
            alphabet=alphabet,
            policies=policies,
            baseline_policy=baseline,
            shape_order=tuple(audit["core_sequence_shapes"]),
            case_order=tuple(audit["core_sequence_cases"]),
            verdict_order=tuple(audit["localisation_verdicts"]),
            conventions=tuple(audit["cleavage_conventions"]),
            reconciliation_verdicts=tuple(audit["reconciliation_verdicts"]),
            panel=tuple(audit["model_panel"]),
            threshold=audit["threshold"],
            relation_order=tuple(audit["boundary_relations"]),
            deletion_outcomes=tuple(audit["deletion_outcomes"]),
            resources=load_declared(spec["external_resources"]),
        )

    # ---------------------------------------------------------- the transforms

    def ingest(self, raw: object, policy: str | None = None) -> list[str]:
        """Read `core_sequence` into a list of cores under a declared policy.

        The two legitimate JSON shapes are a string (one core) and a list of
        strings (several). Everything else is passed through untouched so the
        alphabet check can call it malformed, which is the honest verdict for a
        field holding the printed form of a list.
        """
        chosen = self.policies.get(policy or self.baseline_policy)
        if chosen is None:
            raise PrecursorError(f"{policy!r} is not a declared ingestion policy")
        if isinstance(raw, str):
            cores = [raw]
        elif isinstance(raw, list) and all(isinstance(x, str) for x in raw):
            cores = list(raw)
        else:
            raise PrecursorError(
                f"core_sequence is a {type(raw).__name__}; this campaign reads "
                "only a string or a list of strings")
        if chosen.repair_printed_list:
            repaired: list[str] = []
            for core in cores:
                text = core.strip()
                if _PRINTED_LIST.match(text):
                    repaired.extend(part.strip().strip("'\"")
                                    for part in text[1:-1].split(","))
                else:
                    repaired.append(core)
            cores = repaired
        if chosen.case_fold:
            cores = [core.upper() for core in cores]
        return cores

    def well_formed(self, core: str) -> bool:
        """Every character is an amino-acid letter, in either case.

        Case is a normalisation and is handled by the policy. A core that is
        malformed is malformed under every policy that does not repair it.
        """
        return bool(core) and all(ch.upper() in self.alphabet for ch in core)

    def convention_start(self, convention: str, location: Mapping) -> int:
        """The 0-based index a declared convention says the core begins at."""
        try:
            lower, upper = int(location["from"]), int(location["to"])
        except (KeyError, TypeError, ValueError) as exc:
            raise PrecursorError(
                f"a cleavage location without integer from/to bounds: "
                f"{location!r}") from exc
        table = {"to_as_0based_core_start": upper,
                 "to_as_1based_core_start": upper - 1,
                 "from_as_0based_core_start": lower,
                 "from_as_1based_core_start": lower - 1}
        if convention not in table:
            raise PrecursorError(f"{convention!r} is not a declared convention")
        return table[convention]


# ------------------------------------------------------------------- GenBank


#: GenBank feature table: a feature key starts at column 5, qualifiers at 21.
_FEATURE = re.compile(r"^ {5}(\S+)\s+(\S.*)$")
_QUALIFIER = re.compile(r"^ {21}/(\w+)(?:=(.*))?$")

#: The qualifiers a precursor gene can be named by, in the order they are tried.
#: MIBiG names a precursor either by protein accession or by gene name, and the
#: reference records carry both -- but not always on the same feature, so the
#: first match wins and the order is declared rather than incidental.
NAME_QUALIFIERS = ("protein_id", "locus_tag", "gene", "label")


def cds_translations(text: str) -> dict[str, str]:
    """Every CDS translation in a GenBank record, keyed by each of its names.

    A translation is reached by whichever of the declared name qualifiers the
    feature carries, and the FIRST feature to claim a name keeps it: a record
    where two CDS features share a gene name is ambiguous, and resolving the
    ambiguity by taking the last one would depend on file order.
    """
    out: dict[str, str] = {}
    in_features = False
    key: str | None = None
    qualifiers: dict[str, list[str]] = {}
    pending: tuple[str, list[str]] | None = None

    def flush_qualifier() -> None:
        nonlocal pending
        if pending is None:
            return
        name, parts = pending
        value = "".join(parts)
        if value.startswith('"'):
            value = value[1:]
        if value.endswith('"'):
            value = value[:-1]
        qualifiers.setdefault(name, []).append(value)
        pending = None

    def flush_feature() -> None:
        flush_qualifier()
        if key != "CDS":
            return
        translations = qualifiers.get("translation") or []
        if not translations:
            return
        translation = re.sub(r"\s+", "", translations[0])
        for qualifier in NAME_QUALIFIERS:
            for name in qualifiers.get(qualifier, []):
                out.setdefault(name, translation)

    for line in text.splitlines():
        if not in_features:
            in_features = line.startswith("FEATURES")
            continue
        # The feature table ends at the first line that starts in column 0 --
        # ORIGIN, CONTIG, a second record's LOCUS, or the // terminator.
        if line[:1] not in ("", " "):
            flush_feature()
            break
        feature = _FEATURE.match(line)
        if feature:
            flush_feature()
            key, qualifiers = feature.group(1), {}
            continue
        qualifier = _QUALIFIER.match(line)
        if qualifier:
            flush_qualifier()
            pending = (qualifier.group(1), [qualifier.group(2) or ""])
            continue
        if pending is not None:
            # A continuation line. Translations wrap mid-sequence with leading
            # whitespace, so the parts are joined and the whitespace stripped
            # once, at the end -- stripping per line would also join two words
            # of a free-text qualifier without a space.
            pending[1].append(line.strip())
    else:
        flush_feature()
    return out


# ------------------------------------------------------------------ the corpus


@dataclass(frozen=True)
class Record:
    """One precursor record that resolved to a translation."""

    bgc: str
    gene: str
    status: str
    ripp_type: str | None
    raw_core_sequence: object
    cleavage_location: Mapping | None
    translation: str

    @property
    def key(self) -> tuple[str, str]:
        return (self.bgc, self.gene)


@dataclass(frozen=True)
class Corpus:
    rules: Rules
    release: str
    entries_in_class: int
    entries_with_a_named_precursor: int
    entries_excluded_by_status: int
    entries_absent_from_reference_set: int
    entries_excluded_by_both: int
    records_named: int
    records_unresolved_gene: int
    records: tuple[Record, ...]
    reference_records_shipped: int

    @classmethod
    def load(cls, campaign: pathlib.Path, annotations: str,
             reference_records: str) -> "Corpus":
        campaign = pathlib.Path(campaign)
        rules = Rules.load(campaign)
        provenance = json.loads(
            (campaign / "inputs" / "provenance.json").read_text())

        entries = _read_mibig(campaign / annotations)
        shipped = _reference_index(campaign / reference_records)

        in_class, named, excluded, absent, both = 0, 0, 0, 0, 0
        records_named, unresolved = 0, 0
        records: list[Record] = []
        for accession, document in sorted(entries.items()):
            block = _ribosomal_block(document)
            if block is None:
                continue
            in_class += 1
            precursors = [p for p in (block.get("precursors") or [])
                          if p.get("gene") and p.get("core_sequence") is not None]
            if not precursors:
                continue
            named += 1
            records_named += len(precursors)
            # The two exclusion clauses are counted INDEPENDENTLY over the named
            # entries rather than in sequence. Applied one after the other, the
            # second clause reports only what the first left behind, so two
            # clauses that happen to remove the same entries would report one of
            # them as removing nothing -- which is a fact about the evaluation
            # order and not about the resources. Counted apart, their overlap is
            # itself a measurement.
            retired = document.get("status") in rules.excluded_statuses
            unreferenced = accession not in shipped
            excluded += int(retired)
            absent += int(unreferenced)
            both += int(retired and unreferenced)
            if retired or unreferenced:
                continue
            translations = cds_translations(shipped[accession]())
            for precursor in precursors:
                translation = translations.get(precursor["gene"])
                if translation is None:
                    unresolved += 1
                    continue
                records.append(Record(
                    bgc=accession, gene=precursor["gene"],
                    status=document.get("status"),
                    ripp_type=block.get("ripp_type"),
                    raw_core_sequence=precursor["core_sequence"],
                    cleavage_location=precursor.get("leader_cleavage_location"),
                    translation=translation))
        if not records:
            raise PrecursorError(
                "no precursor record resolved to a translation; the annotation "
                "set and the reference records do not describe the same clusters")
        duplicates = [k for k, n in collections.Counter(
            r.key for r in records).items() if n > 1]
        if duplicates:
            raise PrecursorError(
                f"{len(duplicates)} precursor records share a (cluster, gene) "
                f"key, e.g. {duplicates[0]}; the table cannot be keyed by it")
        return cls(rules=rules, release=provenance["release"],
                   entries_in_class=in_class,
                   entries_with_a_named_precursor=named,
                   entries_excluded_by_status=excluded,
                   entries_absent_from_reference_set=absent,
                   entries_excluded_by_both=both,
                   records_named=records_named,
                   records_unresolved_gene=unresolved,
                   records=tuple(sorted(records, key=lambda r: r.key)),
                   reference_records_shipped=len(shipped))

    # ---------------------------------------------------------- the inventory

    def shape_census(self) -> dict[str, int]:
        """The JSON shapes `core_sequence` arrives in, over the corpus.

        A schema fact: it depends on the field's serialisation and on nothing the
        later rungs compute, which is why it can sit in R1 while the localisation
        it affects sits in R2.
        """
        counts: collections.Counter = collections.Counter()
        for record in self.records:
            raw = record.raw_core_sequence
            if isinstance(raw, str):
                counts["printed_list_in_a_string" if _PRINTED_LIST.match(raw.strip())
                       else "bare_string"] += 1
            elif isinstance(raw, list) and len(raw) == 1:
                counts["list_of_one"] += 1
            elif isinstance(raw, list):
                counts["list_of_several"] += 1
            else:
                counts["other"] += 1
        unknown = set(counts) - set(self.rules.shape_order)
        if unknown:
            raise PrecursorError(
                f"shapes outside the declared vocabulary: {sorted(unknown)}")
        return {shape: counts[shape] for shape in self.rules.shape_order
                if counts.get(shape)}

    def case_census(self) -> dict[str, int]:
        counts: collections.Counter = collections.Counter()
        for record in self.records:
            text = "".join(
                record.raw_core_sequence if isinstance(record.raw_core_sequence, list)
                else [record.raw_core_sequence])
            letters = [ch for ch in text if ch.isalpha()]
            if letters and all(ch.isupper() for ch in letters):
                counts["upper"] += 1
            elif letters and all(ch.islower() for ch in letters):
                counts["lower"] += 1
            else:
                counts["mixed"] += 1
        return {case: counts[case] for case in self.rules.case_order
                if counts.get(case)}

    def ripp_type_census(self) -> dict[str, int]:
        counts = collections.Counter(
            r.ripp_type if r.ripp_type else "(unstated)" for r in self.records)
        return dict(sorted(counts.items()))

    def status_census(self) -> dict[str, int]:
        return dict(sorted(collections.Counter(
            r.status for r in self.records).items()))

    def declaring_a_cleavage_location(self) -> tuple[Record, ...]:
        return tuple(r for r in self.records if r.cleavage_location)

    def precursor_length_summary(self) -> dict:
        lengths = sorted(len(r.translation) for r in self.records)
        return {"count": len(lengths), "min": lengths[0],
                "median": _median(lengths), "max": lengths[-1]}

    # -------------------------------------------------------- the localisation

    def localise(self, record: Record, policy: str | None = None) -> dict:
        """Where the declared core sits in the declared precursor's translation.

        The verdict vocabulary is closed and every term is a property of two
        strings, so nothing here is a judgement. `core_start` is null wherever the
        core has no single unambiguous position, and the three downstream
        quantities are null with it rather than guessed.
        """
        cores = self.rules.ingest(record.raw_core_sequence, policy)
        head = cores[0] if cores else ""
        translation = record.translation
        block = {"bgc": record.bgc, "gene": record.gene, "cores": len(cores),
                 "core_length": len(head) if head else 0,
                 "occurrences": None, "core_start": None,
                 "leader_length": None, "follower_length": None}
        if not self.rules.well_formed(head):
            return {**block, "verdict": "malformed_core"}
        # The core is compared AS THE POLICY PRODUCED IT. Folding case here
        # instead would make every policy case-insensitive and the declared
        # normalisation unobservable: the baseline would be its own
        # counterfactual.
        occurrences = translation.count(head)
        block["occurrences"] = occurrences
        if occurrences == 0:
            return {**block, "verdict": "core_absent"}
        if occurrences > 1:
            return {**block, "verdict": "core_ambiguous"}
        start = translation.find(head)
        if head == translation:
            return {**block, "verdict": "core_is_whole_precursor"}
        block["core_start"] = start
        block["leader_length"] = start
        block["follower_length"] = len(translation) - start - len(head)
        if start == 0:
            return {**block, "verdict": "core_at_N_terminus"}
        if block["follower_length"] == 0:
            return {**block, "verdict": "core_at_C_terminus"}
        return {**block, "verdict": "core_internal"}

    def localisation_table(self, policy: str | None = None) -> list[dict]:
        return [self.localise(r, policy) for r in self.records]

    def localisation_census(self, table: Sequence[Mapping]) -> dict[str, int]:
        counts = collections.Counter(row["verdict"] for row in table)
        unknown = set(counts) - set(self.rules.verdict_order)
        if unknown:
            raise PrecursorError(
                f"localisation produced verdicts outside the declared "
                f"vocabulary: {sorted(unknown)}")
        return {verdict: counts.get(verdict, 0)
                for verdict in self.rules.verdict_order}

    def leader_length_summary(self, table: Sequence[Mapping]) -> dict:
        lengths = sorted(row["leader_length"] for row in table
                         if row["leader_length"] is not None)
        if not lengths:
            raise PrecursorError(
                "no core was located, so there is no leader to measure and the "
                "corpus carries no R2")
        return {"count": len(lengths), "min": lengths[0],
                "median": _median(lengths), "max": lengths[-1]}

    # ------------------------------------------------------ the reconciliation

    def reconcile(self, table: Sequence[Mapping]) -> dict:
        """Score every declared reading of `leader_cleavage_location`.

        The convention is ELECTED by agreement rather than assumed, and the
        tie-break is the declared order, so the answer does not depend on how a
        dict happened to iterate. `underivable` is a real outcome and not a gap:
        it says the core has no single position for a coordinate to agree with.
        """
        located = {(row["bgc"], row["gene"]): row["core_start"] for row in table}
        tally: dict[str, dict[str, int]] = {}
        per_record: list[dict] = []
        for convention in self.rules.conventions:
            counts: collections.Counter = collections.Counter()
            for record in self.declaring_a_cleavage_location():
                derived = located.get(record.key)
                declared = self.rules.convention_start(
                    convention, record.cleavage_location)
                if derived is None:
                    counts["underivable"] += 1
                else:
                    counts["agrees" if declared == derived else "disagrees"] += 1
            tally[convention] = {v: counts.get(v, 0)
                                 for v in self.rules.reconciliation_verdicts}
        elected = max(self.rules.conventions,
                      key=lambda c: (tally[c]["agrees"],
                                     -self.rules.conventions.index(c)))
        for record in self.declaring_a_cleavage_location():
            derived = located.get(record.key)
            declared = self.rules.convention_start(elected, record.cleavage_location)
            per_record.append({
                "bgc": record.bgc, "gene": record.gene,
                "declared_from": int(record.cleavage_location["from"]),
                "declared_to": int(record.cleavage_location["to"]),
                "implied_core_start": declared,
                "derived_core_start": derived,
                "verdict": ("underivable" if derived is None
                            else "agrees" if declared == derived else "disagrees"),
            })
        return {"convention_tally": {c: tally[c] for c in self.rules.conventions},
                "elected_convention": elected,
                "cleavage_table": sorted(per_record,
                                         key=lambda r: (r["bgc"], r["gene"]))}

    # ------------------------------------------------------------- the tool

    def hmmer_version(self) -> str:
        done = _run([str(tool_prefix() / "hmmsearch"), "-h"])
        for line in done.stdout.splitlines():
            if line.startswith("# HMMER"):
                return line[2:].strip()
        raise PrecursorError("hmmsearch did not report a version")

    def fasta(self, table: Sequence[Mapping] | None = None) -> str:
        """The corpus as FASTA, optionally with each declared leader removed.

        With a localisation table, every record whose core was located is written
        from the core's first residue onward and every other record is left out --
        a record with no located core has no declared leader to remove, so
        including it whole would quietly compare a full precursor against the
        others' cores.
        """
        if table is None:
            rows = [(r.bgc, r.gene, r.translation) for r in self.records]
        else:
            starts = {(row["bgc"], row["gene"]): row["core_start"] for row in table}
            rows = [(r.bgc, r.gene, r.translation[starts[r.key]:])
                    for r in self.records if starts.get(r.key) is not None]
        return "".join(f">{bgc}|{gene}\n{seq}\n" for bgc, gene, seq in sorted(rows))

    def _panel(self, work: pathlib.Path) -> pathlib.Path:
        """Extract the declared models from the pinned Pfam library."""
        library = self.rules.resources["pfam_a_hmm"].resolve()
        names = work / "panel.txt"
        names.write_text("\n".join(self.rules.panel) + "\n")
        panel = work / "panel.hmm"
        with open(panel, "w") as handle:
            done = _run([str(tool_prefix() / "hmmfetch"), "-f", str(library),
                         str(names)], stdout=handle)
        if done.returncode != 0:
            raise PrecursorError(f"hmmfetch failed: {done.stderr[-600:]}")
        found = sum(1 for line in panel.read_text().splitlines()
                    if line.startswith("NAME"))
        if found != len(self.rules.panel):
            raise PrecursorError(
                f"the declared panel names {len(self.rules.panel)} models and "
                f"{found} came out of {library}; a campaign cannot score a model "
                "the pinned release does not carry")
        return panel

    def domains(self, text: str) -> dict[tuple[str, str], list[dict]]:
        """Panel hits on a FASTA, keyed by (cluster, gene).

        Memoised when and only when a tool memo is enabled, which the sweep does
        for stub systems alone. The key carries the sequences themselves, so the
        leader-removed scan cannot be answered from the full scan's entry.
        """
        library = self.rules.resources["pfam_a_hmm"]

        def compute() -> list[dict]:
            return [{"bgc": bgc, "gene": gene, **hit}
                    for (bgc, gene), hits in sorted(self._domains(text).items())
                    for hit in hits]

        if not cache.enabled():
            return self._domains(text)
        memo = cache.memoize_json(
            ["ripp_precursor.domains", _digest_text(text), library.sha256,
             self.hmmer_version(), ",".join(self.rules.panel),
             self.rules.threshold],
            compute)
        out: dict[tuple[str, str], list[dict]] = {}
        for row in memo:
            key = (row["bgc"], row["gene"])
            out.setdefault(key, []).append(
                {k: v for k, v in row.items() if k not in ("bgc", "gene")})
        return out

    def _domains(self, text: str) -> dict[tuple[str, str], list[dict]]:
        """The tool call itself. Never memoised, so a caller that wants the real
        thing can always have it."""
        flags = self.rules.threshold.split()
        if not flags:
            raise PrecursorError("an empty threshold is not a declared threshold")
        if not text.strip():
            raise PrecursorError("an empty FASTA is not a scan")
        with tempfile.TemporaryDirectory(prefix="npbench-ripp-") as tmp:
            work = pathlib.Path(tmp)
            panel = self._panel(work)
            fasta = work / "precursors.faa"
            fasta.write_text(text)
            table = work / "hits.domtbl"
            done = _run([str(tool_prefix() / "hmmsearch"), *flags,
                         "--cpu", "2", "--noali", "--domtblout", str(table),
                         str(panel), str(fasta)], stdout=subprocess.DEVNULL)
            if done.returncode != 0:
                raise PrecursorError(f"hmmsearch failed: {done.stderr[-600:]}")
            hits: dict[tuple[str, str], list[dict]] = {}
            for line in table.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                # --domtblout from hmmsearch: column 1 is the target (sequence)
                # and column 4 the query (model). Reading them the other way
                # round reports model names as sequences and is silent about it.
                target, model = fields[0], fields[3]
                if "|" not in target:
                    raise PrecursorError(
                        f"unexpected target name {target!r}; the scan's FASTA is "
                        "written by this campaign and names every record "
                        "cluster|gene")
                bgc, gene = target.split("|", 1)
                hits.setdefault((bgc, gene), []).append({
                    "model": model,
                    "env_from": int(fields[19]),
                    "env_to": int(fields[20]),
                })
        unknown = {h["model"] for v in hits.values() for h in v} - set(self.rules.panel)
        if unknown:
            raise PrecursorError(
                f"hmmsearch reported models outside the declared panel: "
                f"{sorted(unknown)}")
        return {key: sorted(value, key=lambda h: (h["model"], h["env_from"]))
                for key, value in sorted(hits.items())}

    # ----------------------------------------------------- the domain evidence

    def domain_audit(self, table: Sequence[Mapping],
                     hits: Mapping[tuple[str, str], Sequence[Mapping]]) -> dict:
        """Each panel hit against the measured leader/core boundary.

        The envelope is 1-based inclusive as HMMER writes it; the boundary is the
        0-based `core_start`, so the leader occupies 1-based positions 1 through
        `core_start` and the core begins at `core_start + 1`. `boundary_offset` is
        the envelope's last position minus the leader's length: zero means the
        model stops exactly where the leader does.
        """
        starts = {(row["bgc"], row["gene"]): row["core_start"] for row in table}
        rows: list[dict] = []
        for (bgc, gene), entries in sorted(hits.items()):
            start = starts.get((bgc, gene))
            for hit in entries:
                if start is None:
                    relation, offset = "core_not_located", None
                elif hit["env_to"] <= start:
                    relation, offset = "leader_only", hit["env_to"] - start
                elif hit["env_from"] >= start + 1:
                    relation, offset = "core_only", hit["env_to"] - start
                else:
                    relation, offset = "spans_boundary", hit["env_to"] - start
                rows.append({"bgc": bgc, "gene": gene, "model": hit["model"],
                             "env_from": hit["env_from"], "env_to": hit["env_to"],
                             "relation": relation, "boundary_offset": offset})
        counts = collections.Counter(row["relation"] for row in rows)
        unknown = set(counts) - set(self.rules.relation_order)
        if unknown:
            raise PrecursorError(
                f"boundary relations outside the declared vocabulary: "
                f"{sorted(unknown)}")
        panel_counts = collections.Counter(row["model"] for row in rows)
        offsets = sorted(row["boundary_offset"] for row in rows
                         if row["boundary_offset"] is not None)
        return {
            "domain_table": sorted(rows, key=lambda r: (r["bgc"], r["gene"],
                                                        r["model"])),
            "panel_census": {m: panel_counts.get(m, 0) for m in self.rules.panel},
            "records_with_a_panel_hit": len(hits),
            "boundary_relation_census": {r: counts.get(r, 0)
                                         for r in self.rules.relation_order},
            "boundary_offset_summary": (
                {"count": len(offsets), "min": offsets[0],
                 "median": _median(offsets), "max": offsets[-1],
                 "exactly_at_the_boundary": sum(1 for o in offsets if o == 0)}
                if offsets else
                {"count": 0, "min": None, "median": None, "max": None,
                 "exactly_at_the_boundary": 0}),
        }

    # ------------------------------------------------------- the counterfactual

    def policy_sweep(self) -> dict:
        """The localisation census under each declared ingestion policy."""
        out = {}
        for name in sorted(self.rules.policies):
            table = self.localisation_table(name)
            census = self.localisation_census(table)
            out[name] = {
                "verdict_census": census,
                "located": sum(1 for row in table if row["core_start"] is not None),
                "cores_declared": sum(row["cores"] for row in table),
            }
        return out

    def deletion_outcome(self, table: Sequence[Mapping],
                         before: Mapping[tuple[str, str], Sequence[Mapping]],
                         after: Mapping[tuple[str, str], Sequence[Mapping]]) -> dict:
        """Which panel calls survive once the declared leader is removed.

        The scan after deletion covers exactly the records whose core was located,
        so a record with no located core is `not_applicable` rather than counted
        as a loss. `unchanged_absent` is the negative control: a record the panel
        never recognised and still does not.
        """
        rows: list[dict] = []
        starts = {(row["bgc"], row["gene"]): row["core_start"] for row in table}
        for record in self.records:
            was = sorted({h["model"] for h in before.get(record.key, ())})
            now = sorted({h["model"] for h in after.get(record.key, ())})
            if starts.get(record.key) is None:
                outcome = "not_applicable"
            elif not was and not now:
                outcome = "unchanged_absent"
            elif not was:
                outcome = "gained"
            elif not now:
                outcome = "lost"
            elif was == now:
                outcome = "retained"
            else:
                outcome = "changed"
            rows.append({"bgc": record.bgc, "gene": record.gene,
                         "before": was, "after": now, "outcome": outcome})
        counts = collections.Counter(row["outcome"] for row in rows)
        unknown = set(counts) - set(self.rules.deletion_outcomes)
        if unknown:
            raise PrecursorError(
                f"deletion outcomes outside the declared vocabulary: "
                f"{sorted(unknown)}")
        return {
            "deletion_table": [r for r in sorted(
                rows, key=lambda r: (r["bgc"], r["gene"]))
                if r["outcome"] not in ("unchanged_absent", "not_applicable")],
            "deletion_census": {o: counts.get(o, 0)
                                for o in self.rules.deletion_outcomes},
        }


# --------------------------------------------------------------------- loading


def _ribosomal_block(document: Mapping) -> Mapping | None:
    for block in (document.get("biosynthesis") or {}).get("classes") or []:
        if block.get("class") == "ribosomal":
            return block
    return None


def _read_mibig(path: pathlib.Path) -> dict[str, dict]:
    """Every MIBiG entry in the shipped archive, keyed by accession."""
    if not path.is_file():
        raise PrecursorError(f"the annotation archive {path} is not present")
    entries: dict[str, dict] = {}
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".json"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            document = json.loads(handle.read().decode())
            accession = document.get("accession") or pathlib.Path(member.name).stem
            if accession in entries:
                raise PrecursorError(
                    f"{accession} appears twice in {path.name}")
            entries[accession] = document
    if not entries:
        raise PrecursorError(f"{path.name} holds no annotation records")
    return entries


def _reference_index(path: pathlib.Path):
    """Accession -> a callable returning that reference record's text.

    The records are read on demand rather than all at once: the archive holds
    every entry of the class and the corpus uses a fraction of them, and holding
    twelve megabytes of GenBank in memory to read fifty of them is a cost the
    sweep pays once per run.
    """
    if not path.is_file():
        raise PrecursorError(f"the reference archive {path} is not present")
    index: dict[str, str] = {}
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".gbk"):
                continue
            accession = pathlib.Path(member.name).stem
            if accession in index:
                raise PrecursorError(
                    f"{accession} appears twice in {path.name}")
            index[accession] = member.name
    if not index:
        raise PrecursorError(f"{path.name} holds no reference records")

    def reader(name: str):
        def read() -> str:
            with tarfile.open(path, "r:*") as archive:
                handle = archive.extractfile(name)
                if handle is None:
                    raise PrecursorError(f"{name} is not readable in {path.name}")
                return handle.read().decode()
        return read

    return {accession: reader(name) for accession, name in sorted(index.items())}


def _median(ordered: Sequence[int]) -> float:
    n = len(ordered)
    if not n:
        raise PrecursorError("median of an empty set")
    if n % 2:
        return float(ordered[n // 2])
    return round((ordered[n // 2 - 1] + ordered[n // 2]) / 2, PRECISION)
