"""Domain-content audit of an EC-annotated sequence set, computed in-container.

This is the benchmark's first campaign that runs a tool. The claim it grades is
not "what does the literature say about this EC number" -- that would be a G3
gold tier and a different kind of benchmark. It is "which of these sequences
carry the domain the EC number's canonical enzyme requires", answered by running
a pinned HMMER against a pinned Pfam release at a declared threshold. Gold is the
tool's output. Nothing here is a curated label.

Three things are load-bearing.

**The model panel is declared, not discovered.** The campaign publishes the Pfam
models it scores and the agent scans with exactly those. Discovering which
families a mis-annotated set belongs to is a build-time finding, recorded in the
design notes; making it the task would hide the answer in an open-ended search
whose cost is a full Pfam-A scan per run.

**The threshold is declared.** HMMER's verdict is a function of its cutoff, so
`--cut_ga` is pinned and the R4 counterfactual re-runs the audit under other
cutoffs. Measured on this set the choice barely matters, which is a result rather
than an assumption.

**Absence is not the same as misannotation.** A sequence can lack the canonical
domain and still be a real enzyme of a neighbouring family, which is exactly what
this set contains. So "lacks the domain" is computed, and "is misannotated" is a
declared policy over the computed architecture -- published, and varied at R4.

Everything outside the tool call is pure stdlib and sorted order throughout.
"""

from __future__ import annotations

import collections
import gzip
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

TEMPLATE_ID = "ec-domain-audit"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "audit_rules.json"
PRECISION = 6

#: Where the pinned image puts its binaries. The suite and the templates share
#: the override so a different prefix needs no code change.
DEFAULT_PREFIX = "/opt/npbench-tools/bin"

#: Thread-count and locale pinning for every tool call, matching the invariance
#: suite. HMMER's tblout is thread-invariant after one declared normalisation
#: rule, but the environment is fixed anyway so the two measurements are of the
#: same thing.
FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 900


class AuditError(RuntimeError):
    """The sequence set or the environment does not support the audit."""


def tool_prefix() -> pathlib.Path:
    return pathlib.Path(os.environ.get("NPBENCH_TOOL_PREFIX", DEFAULT_PREFIX))


def _run(argv: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a pinned tool with the prefix on PATH.

    The prefix goes last, so a tool environment that carries its own interpreter
    cannot shadow the one this code runs under, and first-or-not it must be on
    PATH at all: a tool that shells out to helpers resolves them there, and
    leaving it off measures the ambient shell instead of the image.
    """
    prefix = str(tool_prefix())
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if prefix in parts:
        parts.remove(prefix)
    env = {**os.environ, **FIXED_ENV, "PATH": os.pathsep.join([*parts, prefix])}
    # capture_output is a shorthand for both pipes, so it cannot be combined with
    # an explicit stdout; set the two separately and let a caller redirect one.
    kwargs.setdefault("stdout", subprocess.PIPE)
    kwargs.setdefault("stderr", subprocess.PIPE)
    return subprocess.run(list(argv), text=True, env=env, timeout=TIMEOUT_S,
                          **kwargs)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# ------------------------------------------------------------------- resources


@dataclass(frozen=True)
class Resource:
    """A pinned file the campaign reads from the image and does not redistribute.

    Some reference data cannot be shipped inside a campaign: it is too large, or
    it carries a licence the benchmark is not in a position to pass on. Such a
    resource is declared here with the path, the provider, the version and a
    fingerprint, so a campaign records exactly what it was built against even
    though the bytes live outside it.
    """

    name: str
    relative_path: str
    provider: str
    version: str
    sha256: str
    note: str

    @property
    def path(self) -> pathlib.Path:
        return tool_prefix().parent / self.relative_path

    def resolve(self) -> pathlib.Path:
        if not self.path.is_file():
            raise AuditError(
                f"{self.name} is not present at {self.path}. This campaign reads "
                f"it from the pinned image ({self.provider}) rather than shipping "
                f"it: {self.note}")
        return self.path

    def verify(self) -> dict:
        """Fingerprint the resource as found, for the record."""
        path = self.resolve()
        found = sha256(path)
        return {"name": self.name, "path": str(path), "provider": self.provider,
                "version": self.version, "declared_sha256": self.sha256,
                "found_sha256": found, "matches": found == self.sha256}


# ----------------------------------------------------------------------- rules


@dataclass(frozen=True)
class Rules:
    ec_number: str
    canonical_model: str
    panel: tuple[str, ...]
    threshold: str
    threshold_variants: tuple[str, ...]
    misannotation_policies: Mapping[str, tuple[str, ...]]
    length_bands: tuple[int, ...]
    resources: Mapping[str, Resource]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        audit = spec["audit"]
        resources = {
            name: Resource(name=name, relative_path=r["relative_path"],
                           provider=r["provider"], version=r["version"],
                           sha256=r["sha256"], note=r["note"])
            for name, r in sorted(spec["external_resources"].items())}
        return cls(
            ec_number=audit["ec_number"],
            canonical_model=audit["canonical_model"],
            panel=tuple(audit["model_panel"]),
            threshold=audit["threshold"],
            threshold_variants=tuple(audit["threshold_variants"]),
            misannotation_policies={
                name: tuple(models)
                for name, models in sorted(audit["misannotation_policies"].items())},
            length_bands=tuple(int(b) for b in audit["length_bands"]),
            resources=resources,
        )


# ------------------------------------------------------------------- sequences


#: UniProt FASTA header: db|accession|entry_name description.
_HEADER = re.compile(r"^(?P<db>[a-z]{2})\|(?P<accession>[^|]+)\|(?P<entry>\S+)")

#: Which UniProt section an entry came from. 'sp' is reviewed (Swiss-Prot) and
#: 'tr' is unreviewed (TrEMBL); the distinction is the campaign's control, so it
#: is read from the header rather than taken on trust from a separate file.
SECTIONS = {"sp": "reviewed", "tr": "unreviewed"}


@dataclass(frozen=True)
class Sequences:
    release: str
    records: Mapping[str, dict]

    @classmethod
    def load(cls, campaign: pathlib.Path, relative: str) -> "Sequences":
        path = pathlib.Path(campaign) / relative
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            text = handle.read()
        provenance = json.loads(
            (pathlib.Path(campaign) / "inputs" / "provenance.json").read_text())
        records: dict[str, dict] = {}
        accession, buffer = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if accession:
                    records[accession]["length"] = len("".join(buffer))
                match = _HEADER.match(line[1:])
                if match is None:
                    raise AuditError(f"unparseable FASTA header: {line[:60]!r}")
                accession = match.group("accession")
                section = SECTIONS.get(match.group("db"))
                if section is None:
                    raise AuditError(
                        f"{accession}: unknown UniProt section "
                        f"{match.group('db')!r}; the reviewed/unreviewed control "
                        "depends on this field")
                records[accession] = {"accession": accession, "section": section,
                                      "entry": match.group("entry")}
                buffer = []
            else:
                buffer.append(line.strip())
        if accession:
            records[accession]["length"] = len("".join(buffer))
        if not records:
            raise AuditError("the shipped FASTA holds no sequences")
        return cls(release=provenance["release"],
                   records=dict(sorted(records.items())))

    def by_section(self) -> dict[str, tuple[str, ...]]:
        out: dict[str, list[str]] = {s: [] for s in sorted(set(SECTIONS.values()))}
        for accession, record in self.records.items():
            out[record["section"]].append(accession)
        return {k: tuple(v) for k, v in sorted(out.items())}


# ----------------------------------------------------------------------- audit


@dataclass(frozen=True)
class Audit:
    sequences: Sequences
    rules: Rules
    input_path: pathlib.Path

    @classmethod
    def load(cls, campaign: pathlib.Path, relative: str) -> "Audit":
        campaign = pathlib.Path(campaign)
        return cls(sequences=Sequences.load(campaign, relative),
                   rules=Rules.load(campaign),
                   input_path=campaign / relative)

    # ------------------------------------------------------------- tool calls

    def hmmer_version(self) -> str:
        done = _run([str(tool_prefix() / "hmmsearch"), "-h"])
        for line in done.stdout.splitlines():
            if line.startswith("# HMMER"):
                return line[2:].strip()
        raise AuditError("hmmsearch did not report a version")

    def _fasta(self, work: pathlib.Path) -> pathlib.Path:
        """A plain FASTA for HMMER, which does not read gzip."""
        if self.input_path.suffix != ".gz":
            return self.input_path
        plain = work / "sequences.fasta"
        with gzip.open(self.input_path, "rt") as src, open(plain, "w") as dst:
            shutil.copyfileobj(src, dst)
        return plain

    def _panel(self, work: pathlib.Path) -> pathlib.Path:
        """Extract the declared models from the pinned Pfam library."""
        library = self.rules.resources["pfam_a_hmm"].resolve()
        names = work / "panel.txt"
        names.write_text("\n".join(self.rules.panel) + "\n")
        panel = work / "panel.hmm"
        done = _run([str(tool_prefix() / "hmmfetch"), "-f", str(library),
                     str(names)], stdout=open(panel, "w"))
        if done.returncode != 0:
            raise AuditError(f"hmmfetch failed: {done.stderr[-600:]}")
        found = sum(1 for line in panel.read_text().splitlines()
                    if line.startswith("NAME"))
        if found != len(self.rules.panel):
            raise AuditError(
                f"the declared panel names {len(self.rules.panel)} models and "
                f"{found} came out of {library}; a campaign cannot score a model "
                "the pinned release does not carry")
        return panel

    def domains(self, threshold: str | None = None) -> dict[str, frozenset[str]]:
        """Accession -> the declared models it carries, at a declared threshold.

        One hmmsearch of the declared panel over the whole set. The models are
        fetched from the pinned library rather than shipped, so what the campaign
        records is the fingerprint and not the bytes.
        """
        flags = (threshold or self.rules.threshold).split()
        if not flags:
            raise AuditError("an empty threshold is not a declared threshold")
        with tempfile.TemporaryDirectory(prefix="npbench-ecaudit-") as tmp:
            work = pathlib.Path(tmp)
            panel = self._panel(work)
            fasta = self._fasta(work)
            table = work / "hits.domtbl"
            done = _run([str(tool_prefix() / "hmmsearch"), *flags,
                         "--cpu", "2", "--noali", "--domtblout", str(table),
                         str(panel), str(fasta)], stdout=subprocess.DEVNULL)
            if done.returncode != 0:
                raise AuditError(f"hmmsearch failed: {done.stderr[-600:]}")
            hits: dict[str, set[str]] = {}
            for line in table.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                # --domtblout: column 1 is the target (sequence), column 4 the
                # query (model). Reading them the other way round reports model
                # names as sequences and is silent about it.
                target, model = fields[0], fields[3]
                accession = _HEADER.match(target)
                key = accession.group("accession") if accession else target
                hits.setdefault(key, set()).add(model)
        unknown = {m for models in hits.values() for m in models} - set(self.rules.panel)
        if unknown:
            raise AuditError(
                f"hmmsearch reported models outside the declared panel: "
                f"{sorted(unknown)}")
        return {a: frozenset(hits.get(a, ())) for a in sorted(self.sequences.records)}

    # -------------------------------------------------------------- the census

    def architecture(self, models: frozenset[str]) -> str:
        """The declared label for a domain content: sorted names, or '(none)'."""
        return " + ".join(sorted(models)) if models else "(none)"

    def presence_census(self, domains: Mapping[str, frozenset[str]]) -> dict:
        """Per UniProt section: canonical model present, absent."""
        canonical = self.rules.canonical_model
        present: dict[str, int] = {}
        absent: dict[str, int] = {}
        for section, accessions in self.sequences.by_section().items():
            present[section] = sum(1 for a in accessions
                                   if canonical in domains[a])
            absent[section] = len(accessions) - present[section]
        return {"canonical_present_by_section": dict(sorted(present.items())),
                "canonical_absent_by_section": dict(sorted(absent.items()))}

    def architecture_census(self, domains: Mapping[str, frozenset[str]],
                            accessions: Sequence[str] | None = None) -> dict[str, int]:
        keys = accessions if accessions is not None else sorted(domains)
        counts = collections.Counter(self.architecture(domains[a]) for a in keys)
        return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))

    def family_census(self, domains: Mapping[str, frozenset[str]],
                      accessions: Sequence[str] | None = None) -> dict[str, int]:
        keys = accessions if accessions is not None else sorted(domains)
        counts: collections.Counter = collections.Counter()
        for accession in keys:
            counts.update(domains[accession])
        return {m: counts.get(m, 0) for m in sorted(self.rules.panel)}

    def absent(self, domains: Mapping[str, frozenset[str]]) -> tuple[str, ...]:
        canonical = self.rules.canonical_model
        return tuple(a for a in sorted(domains) if canonical not in domains[a])

    # ------------------------------------------------------------ the fragment
    #                                                              explanation

    def length_profile(self, domains: Mapping[str, frozenset[str]]) -> dict:
        """Length quartiles and a short-sequence count for each side.

        A sequence can miss a 348-column model simply by being a fragment, which
        would make the whole audit an artefact of sequencing completeness rather
        than a statement about annotation. The profile is graded so that the
        explanation is tested rather than assumed either way.
        """
        canonical = self.rules.canonical_model
        groups = {"canonical_present": [], "canonical_absent": []}
        for accession, models in sorted(domains.items()):
            side = "canonical_present" if canonical in models else "canonical_absent"
            groups[side].append(self.sequences.records[accession]["length"])
        out = {}
        for side, lengths in sorted(groups.items()):
            if not lengths:
                raise AuditError(
                    f"{side} is empty, so the campaign has only one side to "
                    "compare and the control is not a control")
            ordered = sorted(lengths)
            out[side] = {
                "count": len(ordered),
                "min": ordered[0],
                "median": _median(ordered),
                "max": ordered[-1],
                "shorter_than": {str(band): sum(1 for v in ordered if v < band)
                                 for band in self.rules.length_bands},
            }
        return out

    # ------------------------------------------------------------ the verdicts

    def misannotated(self, domains: Mapping[str, frozenset[str]],
                     policy: str) -> tuple[str, ...]:
        """Accessions a declared policy counts as misannotated.

        The policy names the models that are ACCEPTED as an architecture for this
        EC number. Everything else is counted. The strict policy accepts only the
        canonical model; a broader one accepts a neighbouring family as well, and
        the point of varying it is that the two answers differ by hundreds of
        sequences on real data.
        """
        accepted = self.rules.misannotation_policies.get(policy)
        if accepted is None:
            raise AuditError(f"{policy!r} is not a declared misannotation policy")
        out = []
        for accession in sorted(domains):
            if not (domains[accession] & frozenset(accepted)):
                out.append(accession)
        return tuple(out)

    def policy_census(self, domains: Mapping[str, frozenset[str]]) -> dict:
        return {name: len(self.misannotated(domains, name))
                for name in sorted(self.rules.misannotation_policies)}


def _median(ordered: Sequence[int]) -> float:
    n = len(ordered)
    if not n:
        raise AuditError("median of an empty length set")
    if n % 2:
        return float(ordered[n // 2])
    return round((ordered[n // 2 - 1] + ordered[n // 2]) / 2, PRECISION)
