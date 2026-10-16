"""Domain architecture of a multidomain biosynthetic protein, and who agrees.

An NRPS or PKS is a single protein carrying a chain of catalytic domains, and
"the architecture" is the ordered list of those domains. Computing it from a
pinned HMMER against a pinned Pfam release sounds like a lookup. It is not: the
answer is a function of three choices nobody writes down, and this campaign
declares all three and measures what each one costs.

Five things are load-bearing.

**Hits overlap, so an architecture needs a resolution rule.** The catalogued
design called for a "non-overlapping" domain table; measured on this corpus that
is simply false -- nested Pfam families (`ketoacyl-synt` inside `Thiolase_N`,
`adh_short` against `KR`, five methyltransferase families against each other) hit
the same residues. Resolution is declared, the baseline resolves, and the
unresolved architecture is one of R4's axes.

**Pfam's own clan rule is NOT available.** Pfam resolves intra-clan overlaps by
score, but the pinned library carries no clan annotation at all -- zero `CL`
lines across 19,632 models -- so the rule here is score and geometry, declared as
such rather than borrowed and misattributed.

**The coordinate convention is a choice.** HMMER reports an alignment span and a
wider envelope span, and ordering by one or the other can reorder the
architecture. The baseline is declared and the alternative is measured.

**A Pfam family is not a domain type.** The curated records name an adenylation
domain; Pfam has `AMP-binding` and `AMP-binding_C`, which are two halves of one
thing. Which family "is" the domain is a declared mapping, and R4 varies it,
because counting both double-counts every module.

**The reconciliation target is independent of the tool.** MIBiG's modules are
curated from the literature and the product's structure, gene by gene; the
architecture is computed from sequence. Comparing their COUNTS is arithmetic. The
curated coordinates are not used, because they are unusable -- the campaign
records how unusable rather than quietly working around it.

Everything outside the HMMER calls is pure stdlib and sorted order throughout.
"""

from __future__ import annotations

import collections
import gzip
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

TEMPLATE_ID = "domain-architecture"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "architecture_rules.json"
PRECISION = 6

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 3600


class ArchitectureError(RuntimeError):
    """The protein set or the environment does not support the audit."""


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


def _digest_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# ----------------------------------------------------------------------- rules


@dataclass(frozen=True)
class Rules:
    module_family: str
    accepted_module_types: tuple[str, ...]
    excluded_statuses: tuple[str, ...]
    panel: tuple[str, ...]
    cutoffs: tuple[str, ...]
    baseline_cutoff: str
    coordinate_spans: tuple[str, ...]
    baseline_span: str
    overlap_fraction: float
    resolution_policies: tuple[str, ...]
    baseline_resolution: str
    domain_mappings: Mapping[str, tuple[str, ...]]
    baseline_mapping: str
    curated_domain_keys: tuple[str, ...]
    primary_domain: str
    secondary_domain: str
    module_pattern: tuple[str, ...]
    strata: Mapping[str, tuple[str, ...]]
    verdicts: tuple[str, ...]
    resources: Mapping[str, Resource]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        a = spec["audit"]
        mappings = {name: tuple(fams)
                    for name, fams in sorted(a["domain_mappings"].items())}
        for field, choices in (("baseline_cutoff", a["cutoffs"]),
                               ("baseline_span", a["coordinate_spans"]),
                               ("baseline_resolution", a["resolution_policies"]),
                               ("baseline_mapping", sorted(mappings)),
                               ("primary_domain", a["curated_domain_keys"]),
                               ("secondary_domain", a["curated_domain_keys"])):
            if a[field] not in choices:
                raise ArchitectureError(
                    f"{field} is {a[field]!r}, which is not among the declared "
                    f"choices {list(choices)}")
        if a["primary_domain"] == a["secondary_domain"]:
            raise ArchitectureError(
                "the primary and secondary reconciled domain types are the same "
                "key, so the two reconciliations would be one column twice")
        # The secondary domain is counted through the mapping whose key is its
        # own name. If that keying drifts, the secondary reconciliation computes
        # zero for every gene and nothing else complains -- which is exactly what
        # happened when this template was generalised, so it is checked here.
        if a["secondary_domain"] not in mappings:
            raise ArchitectureError(
                f"the secondary domain {a['secondary_domain']!r} has no entry in "
                f"domain_mappings (declared: {sorted(mappings)}), so it would be "
                "counted as absent on every gene")
        if a["baseline_mapping"] == a["secondary_domain"]:
            raise ArchitectureError(
                "the baseline mapping is the secondary domain's own mapping, so "
                "the primary reconciliation would count the secondary domain")
        strata = {name: tuple(sorted(types))
                  for name, types in sorted(a["module_strata"].items())}
        covered = {t for types in strata.values() for t in types}
        if covered != set(a["accepted_module_types"]):
            raise ArchitectureError(
                "the declared strata must partition the accepted module types "
                f"exactly; accepted {sorted(a['accepted_module_types'])}, "
                f"covered {sorted(covered)}")
        return cls(
            module_family=a["module_family"],
            accepted_module_types=tuple(sorted(a["accepted_module_types"])),
            excluded_statuses=tuple(a["excluded_statuses"]),
            panel=tuple(a["model_panel"]),
            cutoffs=tuple(a["cutoffs"]),
            baseline_cutoff=a["baseline_cutoff"],
            coordinate_spans=tuple(a["coordinate_spans"]),
            baseline_span=a["baseline_span"],
            overlap_fraction=float(a["overlap_fraction"]),
            resolution_policies=tuple(a["resolution_policies"]),
            baseline_resolution=a["baseline_resolution"],
            domain_mappings=mappings,
            baseline_mapping=a["baseline_mapping"],
            curated_domain_keys=tuple(a["curated_domain_keys"]),
            primary_domain=a["primary_domain"],
            secondary_domain=a["secondary_domain"],
            module_pattern=tuple(a["module_pattern"]),
            strata=strata,
            verdicts=tuple(a["reconciliation_verdicts"]),
            resources=load_declared(spec["external_resources"]),
        )


# ------------------------------------------------------------------- sequences


#: The shipped FASTA names every record `<cluster>|<gene>`, written by this
#: template's own prepare_input, so a header that does not parse is a broken
#: input rather than something to tolerate.
_HEADER = re.compile(r"^(?P<bgc>BGC\d+)\|(?P<gene>\S+)$")


def read_proteins(path: pathlib.Path) -> dict[tuple[str, str], str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        text = handle.read()
    out: dict[tuple[str, str], str] = {}
    key: tuple[str, str] | None = None
    buffer: list[str] = []

    def flush() -> None:
        if key is None:
            return
        if key in out:
            raise ArchitectureError(f"{key} appears twice in {path.name}")
        out[key] = "".join(buffer)

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            match = _HEADER.match(line[1:].strip())
            if match is None:
                raise ArchitectureError(
                    f"unparseable FASTA header {line[:60]!r}; this campaign's "
                    "own input names every record cluster|gene")
            key = (match.group("bgc"), match.group("gene"))
            buffer = []
        else:
            buffer.append(line.strip())
    flush()
    if not out:
        raise ArchitectureError(f"{path.name} holds no sequences")
    return dict(sorted(out.items()))


# --------------------------------------------------------------- curated side


@dataclass(frozen=True)
class CuratedGene:
    bgc: str
    gene: str
    status: str
    modules: int
    module_types: tuple[str, ...]
    domain_counts: Mapping[str, int]
    carriers: int
    modification_domains: int

    @property
    def key(self) -> tuple[str, str]:
        return (self.bgc, self.gene)


def read_curated(path: pathlib.Path, rules: Rules) -> tuple[
        dict[tuple[str, str], CuratedGene], dict]:
    """Every gene a curated module names, plus the census of what was dropped.

    A module naming more than one gene is excluded, because the per-gene counts
    this campaign reconciles cannot be apportioned between them without a rule
    nothing in the record supplies.
    """
    if not path.is_file():
        raise ArchitectureError(f"the annotation archive {path} is not present")
    census: collections.Counter = collections.Counter()
    coordinates: collections.Counter = collections.Counter()
    active: collections.Counter = collections.Counter()
    # Per module TYPE, how many of its modules record each typed domain block.
    # Counted at module level and not per gene, because a gene can carry modules
    # of two types and attributing its totals to both would make a type that
    # never records a domain look as though it sometimes does.
    by_type: dict[str, collections.Counter] = {}
    genes: dict[tuple[str, str], dict] = {}
    statuses: dict[str, str] = {}
    with tarfile.open(path, "r:*") as archive:
        for member in sorted(archive.getmembers(), key=lambda m: m.name):
            if not member.isfile() or not member.name.endswith(".json"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            document = json.loads(handle.read().decode())
            accession = document.get("accession") or pathlib.Path(member.name).stem
            modules = (document.get("biosynthesis") or {}).get("modules") or []
            if not modules:
                continue
            census["entries_with_modules"] += 1
            statuses[accession] = document.get("status")
            for module in modules:
                census["modules"] += 1
                active[repr(module.get("active"))] += 1
                typed = by_type.setdefault(module.get("type"),
                                           collections.Counter())
                for key in rules.curated_domain_keys:
                    block = module.get(key)
                    if not block:
                        continue
                    typed[key] += 1
                    location = block.get("location") or {}
                    usable = (isinstance(location.get("from"), int)
                              and isinstance(location.get("to"), int)
                              and location["from"] >= 0 and location["to"] > 0)
                    coordinates["usable" if usable else "placeholder"] += 1
                named = module.get("genes") or []
                if len(named) != 1:
                    census["modules_spanning_several_genes"] += 1
                    continue
                block = genes.setdefault((accession, named[0]), {
                    "modules": 0, "types": [], "carriers": 0,
                    "modification_domains": 0,
                    "domains": collections.Counter()})
                block["modules"] += 1
                block["types"].append(module.get("type"))
                block["carriers"] += len(module.get("carriers") or [])
                block["modification_domains"] += len(
                    module.get("modification_domains") or [])
                for key in rules.curated_domain_keys:
                    if module.get(key):
                        block["domains"][key] += 1
    if not genes:
        raise ArchitectureError(f"{path.name} holds no curated modules")
    curated = {
        key: CuratedGene(bgc=key[0], gene=key[1],
                         status=statuses.get(key[0]),
                         modules=block["modules"],
                         module_types=tuple(sorted(set(block["types"]))),
                         domain_counts={k: block["domains"].get(k, 0)
                                        for k in rules.curated_domain_keys},
                         carriers=block["carriers"],
                         modification_domains=block["modification_domains"])
        for key, block in sorted(genes.items())}
    census["genes_named_by_a_single_gene_module"] = len(curated)
    return curated, {
        "annotation_census": dict(sorted(census.items())),
        "curated_coordinate_usability": dict(sorted(coordinates.items())),
        "curated_active_flag": dict(sorted(active.items())),
        "domains_by_module_type": {
            module_type: dict(sorted(counts.items()))
            for module_type, counts in sorted(by_type.items())
            if module_type is not None},
    }


# ------------------------------------------------------------------ the corpus


@dataclass(frozen=True)
class Hit:
    model: str
    ali: tuple[int, int]
    env: tuple[int, int]
    score: float

    def span(self, which: str) -> tuple[int, int]:
        return self.ali if which == "ali" else self.env


@dataclass(frozen=True)
class Corpus:
    rules: Rules
    release: str
    proteins: Mapping[tuple[str, str], str]
    curated: Mapping[tuple[str, str], CuratedGene]
    annotation_census: Mapping[str, object]
    keys: tuple[tuple[str, str], ...]
    shipped_proteins: int
    curated_genes: int
    #: memo for `_absent_strata`, which walks the whole curated set.
    _absent_cache: frozenset[str] | None = None

    @classmethod
    def load(cls, campaign: pathlib.Path, annotations: str,
             proteins: str) -> "Corpus":
        campaign = pathlib.Path(campaign)
        rules = Rules.load(campaign)
        provenance = json.loads(
            (campaign / "inputs" / "provenance.json").read_text())
        sequences = read_proteins(campaign / proteins)
        curated, census = read_curated(campaign / annotations, rules)

        # The declared selection: every gene whose curated module types are ALL
        # among the accepted ones, on an entry the status rule admits, with a
        # shipped translation. "All accepted" rather than "any" because a gene
        # carrying both NRPS and PKS modules has an architecture drawn from two
        # vocabularies and belongs to neither instantiation. A family can accept
        # several types -- a modular PKS gene normally carries a loading module
        # of its own type alongside its elongation modules -- and the strata
        # below keep those distinguishable.
        keys = []
        # A closed vocabulary: the selection rule has three clauses and the
        # census names all three whether or not each removed anything. A Counter
        # left to itself omits the clauses that happened to remove nothing, which
        # makes the census's SHAPE depend on the data and a reader cannot tell an
        # absent key from a clause that does not exist.
        drops: collections.Counter = collections.Counter(
            {"excluded_by_status": 0, "other_module_type": 0,
             "no_shipped_translation": 0})
        for key, gene in curated.items():
            if gene.status in rules.excluded_statuses:
                drops["excluded_by_status"] += 1
                continue
            if not set(gene.module_types) <= set(rules.accepted_module_types):
                drops["other_module_type"] += 1
                continue
            if key not in sequences:
                drops["no_shipped_translation"] += 1
                continue
            keys.append(key)
        if not keys:
            raise ArchitectureError(
                f"no gene survived the selection for module family "
                f"{rules.module_family!r}")
        census = {**census, "selection_drops": dict(sorted(drops.items()))}
        return cls(rules=rules, release=provenance["release"],
                   proteins=sequences, curated=curated,
                   annotation_census=census, keys=tuple(sorted(keys)),
                   shipped_proteins=len(sequences), curated_genes=len(curated))

    # ---------------------------------------------------------- the inventory

    def length_summary(self) -> dict:
        lengths = sorted(len(self.proteins[k]) for k in self.keys)
        return {"count": len(lengths), "min": lengths[0],
                "median": _median(lengths), "max": lengths[-1],
                "residues": sum(lengths)}

    def curated_module_census(self) -> dict[str, int]:
        counts = collections.Counter(
            self.curated[k].modules for k in self.keys)
        return {str(n): counts[n] for n in sorted(counts)}

    def curated_domain_totals(self) -> dict[str, int]:
        totals: collections.Counter = collections.Counter()
        for key in self.keys:
            for name, n in self.curated[key].domain_counts.items():
                totals[name] += n
        return {k: totals.get(k, 0) for k in self.rules.curated_domain_keys}

    # ------------------------------------------------------------- the tool

    def hmmer_version(self) -> str:
        done = _run([str(tool_prefix() / "hmmsearch"), "-h"])
        for line in done.stdout.splitlines():
            if line.startswith("# HMMER"):
                return line[2:].strip()
        raise ArchitectureError("hmmsearch did not report a version")

    def fasta(self) -> str:
        return "".join(f">{bgc}|{gene}\n{self.proteins[(bgc, gene)]}\n"
                       for bgc, gene in self.keys)

    def _panel(self, work: pathlib.Path) -> pathlib.Path:
        library = self.rules.resources["pfam_a_hmm"].resolve()
        names = work / "panel.txt"
        names.write_text("\n".join(self.rules.panel) + "\n")
        panel = work / "panel.hmm"
        with open(panel, "w") as handle:
            done = _run([str(tool_prefix() / "hmmfetch"), "-f", str(library),
                         str(names)], stdout=handle)
        if done.returncode != 0:
            raise ArchitectureError(f"hmmfetch failed: {done.stderr[-600:]}")
        found = sum(1 for line in panel.read_text().splitlines()
                    if line.startswith("NAME"))
        if found != len(self.rules.panel):
            raise ArchitectureError(
                f"the declared panel names {len(self.rules.panel)} models and "
                f"{found} came out of {library}; a campaign cannot score a model "
                "the pinned release does not carry")
        return panel

    def hits(self, cutoff: str | None = None) -> dict[tuple[str, str], list[Hit]]:
        """Panel hits on the corpus at a declared cutoff, keyed by gene.

        Memoised when and only when a tool memo is enabled, which the sweep does
        for stub systems alone. The key carries the sequences, the library, the
        tool version, the panel and the cutoff -- a memo missing any of those
        would answer from the wrong world.
        """
        chosen = cutoff or self.rules.baseline_cutoff
        if chosen not in self.rules.cutoffs:
            raise ArchitectureError(f"{chosen!r} is not a declared cutoff")
        library = self.rules.resources["pfam_a_hmm"]
        text = self.fasta()

        def compute() -> list[dict]:
            return [{"bgc": bgc, "gene": gene, "model": h.model,
                     "ali_from": h.ali[0], "ali_to": h.ali[1],
                     "env_from": h.env[0], "env_to": h.env[1], "score": h.score}
                    for (bgc, gene), hs in sorted(self._hits(text, chosen).items())
                    for h in hs]

        if not cache.enabled():
            return self._hits(text, chosen)
        memo = cache.memoize_json(
            ["domain_architecture.hits", _digest_text(text), library.sha256,
             self.hmmer_version(), ",".join(self.rules.panel), chosen],
            compute)
        out: dict[tuple[str, str], list[Hit]] = {}
        for row in memo:
            out.setdefault((row["bgc"], row["gene"]), []).append(
                Hit(model=row["model"],
                    ali=(row["ali_from"], row["ali_to"]),
                    env=(row["env_from"], row["env_to"]),
                    score=row["score"]))
        return {k: v for k, v in sorted(out.items())}

    def _hits(self, text: str, cutoff: str) -> dict[tuple[str, str], list[Hit]]:
        """The tool call itself. Never memoised."""
        flags = cutoff.split()
        if not flags:
            raise ArchitectureError("an empty cutoff is not a declared cutoff")
        with tempfile.TemporaryDirectory(prefix="npbench-arch-") as tmp:
            work = pathlib.Path(tmp)
            panel = self._panel(work)
            fasta = work / "proteins.faa"
            fasta.write_text(text)
            table = work / "hits.domtbl"
            done = _run([str(tool_prefix() / "hmmsearch"), *flags,
                         "--cpu", "2", "--noali", "--domtblout", str(table),
                         str(panel), str(fasta)], stdout=subprocess.DEVNULL)
            if done.returncode != 0:
                raise ArchitectureError(f"hmmsearch failed: {done.stderr[-600:]}")
            out: dict[tuple[str, str], list[Hit]] = {}
            for line in table.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                # hmmsearch --domtblout: column 1 the target (sequence), column 4
                # the query (model); columns 18-19 the alignment span and 20-21
                # the envelope.
                match = _HEADER.match(f[0])
                if match is None:
                    raise ArchitectureError(f"unexpected target name {f[0]!r}")
                out.setdefault((match.group("bgc"), match.group("gene")), []).append(
                    Hit(model=f[3], ali=(int(f[17]), int(f[18])),
                        env=(int(f[19]), int(f[20])), score=float(f[13])))
        unknown = {h.model for v in out.values() for h in v} - set(self.rules.panel)
        if unknown:
            raise ArchitectureError(
                f"hmmsearch reported models outside the declared panel: "
                f"{sorted(unknown)}")
        return {k: sorted(v, key=lambda h: (h.ali, h.model))
                for k, v in sorted(out.items())}

    # ------------------------------------------------------ the architecture

    def resolve(self, hits: Sequence[Hit], span: str | None = None,
                policy: str | None = None) -> list[Hit]:
        """The hits a declared policy keeps, in coordinate order.

        `keep_all` keeps everything. `resolve_by_score` walks the hits best score
        first and drops any that overlaps an already-kept hit by more than the
        declared fraction of the SHORTER span -- the fraction is of the shorter
        one because a short domain nested inside a long one is the case the rule
        exists for, and measuring against the longer span would never trigger.

        The walk order is (-score, model, span) so it is total: two hits with the
        same score must still resolve the same way on every run.
        """
        which = span or self.rules.baseline_span
        chosen = policy or self.rules.baseline_resolution
        if which not in self.rules.coordinate_spans:
            raise ArchitectureError(f"{which!r} is not a declared coordinate span")
        if chosen not in self.rules.resolution_policies:
            raise ArchitectureError(f"{chosen!r} is not a declared resolution policy")
        order = sorted(hits, key=lambda h: (-h.score, h.model, h.span(which)))
        kept: list[Hit] = []
        for hit in order:
            if chosen == "resolve_by_score":
                lo_a, hi_a = hit.span(which)
                clash = False
                for other in kept:
                    lo_b, hi_b = other.span(which)
                    overlap = max(0, min(hi_a, hi_b) - max(lo_a, lo_b) + 1)
                    shorter = min(hi_a - lo_a, hi_b - lo_b) + 1
                    if overlap / shorter > self.rules.overlap_fraction:
                        clash = True
                        break
                if clash:
                    continue
            kept.append(hit)
        return sorted(kept, key=lambda h: (h.span(which), h.model))

    def architecture_table(self, hits: Mapping[tuple[str, str], Sequence[Hit]],
                           span: str | None = None,
                           policy: str | None = None) -> list[dict]:
        which = span or self.rules.baseline_span
        rows = []
        for bgc, gene in self.keys:
            raw = list(hits.get((bgc, gene), ()))
            kept = self.resolve(raw, which, policy)
            rows.append({
                "bgc": bgc, "gene": gene,
                "length": len(self.proteins[(bgc, gene)]),
                "hits": len(raw), "domains": len(kept),
                "resolved_away": len(raw) - len(kept),
                "architecture": "-".join(h.model for h in kept),
                "first_domain_start": kept[0].span(which)[0] if kept else None,
                "last_domain_end": kept[-1].span(which)[1] if kept else None,
            })
        return rows

    def architecture_census(self, table: Sequence[Mapping]) -> dict:
        per_gene = sorted(row["domains"] for row in table)
        return {
            "genes": len(table),
            "genes_with_no_domain": sum(1 for row in table if not row["domains"]),
            "distinct_architectures": len({row["architecture"] for row in table}),
            "domains_total": sum(row["domains"] for row in table),
            "domains_per_gene": {"min": per_gene[0], "median": _median(per_gene),
                                 "max": per_gene[-1]},
        }

    def overlap_census(self, hits: Mapping[tuple[str, str], Sequence[Hit]],
                       span: str | None = None) -> dict:
        """How much of the raw hit set overlaps, before any resolution.

        The number the catalogued "non-overlapping domain table" assumed was
        zero. It is not, and the pairs are counted here so that the resolution
        policy is justified by a measurement rather than by an assertion.
        """
        which = span or self.rules.baseline_span
        pairs = 0
        genes = 0
        total = 0
        for key in self.keys:
            order = sorted(hits.get(key, ()), key=lambda h: h.span(which))
            total += len(order)
            found = 0
            for i, left in enumerate(order):
                lo_a, hi_a = left.span(which)
                for right in order[i + 1:]:
                    lo_b, hi_b = right.span(which)
                    if lo_b > hi_a:
                        break
                    if min(hi_a, hi_b) - max(lo_a, lo_b) + 1 > 0:
                        found += 1
            pairs += found
            if found:
                genes += 1
        return {"hits_total": total, "overlapping_pairs": pairs,
                "genes_with_an_overlap": genes}

    # ------------------------------------------------------ the reconciliation

    def reconcile(self, hits: Mapping[tuple[str, str], Sequence[Hit]],
                  mapping: str | None = None, span: str | None = None,
                  policy: str | None = None) -> dict:
        """Curated module counts against computed domain counts, gene by gene.

        The curated side is a human reading of the literature and the product's
        structure; the computed side is HMMER on sequence. Only COUNTS are
        compared, because the curated coordinates are placeholders on all but a
        handful of records -- which R1 reports rather than works around.

        `span` and `policy` are passed through to the resolution rather than
        taken from the baseline, because the counterfactual asks for the
        agreement a DIFFERENT reading produces. An earlier version took the
        architecture table as an argument and then re-resolved at the baseline,
        so the sweep's agreement column silently ignored two of its three axes;
        the numbers happened to coincide on this corpus, which is exactly why it
        was worth fixing rather than leaving.
        """
        chosen = mapping or self.rules.baseline_mapping
        families = self.rules.domain_mappings.get(chosen)
        if families is None:
            raise ArchitectureError(f"{chosen!r} is not a declared domain mapping")
        secondary = self.rules.domain_mappings.get(self.rules.secondary_domain, ())
        rows = []
        for key in self.keys:
            kept = self.resolve(hits.get(key, ()), span, policy)
            counts = collections.Counter(h.model for h in kept)
            curated = self.curated[key]
            computed = sum(counts.get(f, 0) for f in families)
            declared = curated.domain_counts.get(self.rules.primary_domain, 0)
            s_computed = sum(counts.get(f, 0) for f in secondary)
            s_declared = curated.domain_counts.get(self.rules.secondary_domain, 0)
            rows.append({
                "bgc": key[0], "gene": key[1],
                "stratum": self.stratum(key),
                "curated_modules": curated.modules,
                "curated_primary": declared,
                "computed_primary": computed,
                "primary_delta": computed - declared,
                "primary_verdict": _verdict(computed, declared),
                "curated_secondary": s_declared,
                "computed_secondary": s_computed,
                # A curated count of zero is `not_curated` only where the record
                # simply does not say. Where the module TYPE means the domain is
                # absent -- a trans-AT module has no acyltransferase by
                # definition -- zero is a claim, and the campaign grades it as
                # one; that is the data's own negative control.
                "secondary_verdict": (
                    _verdict(s_computed, s_declared)
                    if s_declared or self.stratum_asserts_absence(key)
                    else "not_curated"),
            })
        a_counts = collections.Counter(r["primary_verdict"] for r in rows)
        c_counts = collections.Counter(r["secondary_verdict"] for r in rows)
        unknown = (set(a_counts) | set(c_counts)) - set(self.rules.verdicts)
        if unknown:
            raise ArchitectureError(
                f"reconciliation verdicts outside the declared vocabulary: "
                f"{sorted(unknown)}")
        deltas = collections.Counter(r["primary_delta"] for r in rows)
        # "mixed" is a stratum `stratum()` can return and the declared list
        # cannot contain -- a gene spanning two strata belongs to neither -- so
        # the census has to make room for it. The NRPS instantiation has one
        # stratum and never produces it, which is why this only surfaced when
        # the PKS family was built.
        by_stratum: dict[str, collections.Counter] = {
            name: collections.Counter() for name in [*self.rules.strata, "mixed"]}
        for row in rows:
            by_stratum[row["stratum"]][row["secondary_verdict"]] += 1
        return {
            "reconciliation_table": rows,
            "primary_verdicts": {v: a_counts.get(v, 0) for v in self.rules.verdicts},
            "secondary_verdicts": {v: c_counts.get(v, 0) for v in self.rules.verdicts},
            "primary_delta_census": {str(d): deltas[d] for d in sorted(deltas)},
            "secondary_verdicts_by_stratum": {
                name: {v: counts.get(v, 0) for v in self.rules.verdicts}
                for name, counts in sorted(by_stratum.items())},
        }

    def stratum(self, key: tuple[str, str]) -> str:
        """The declared stratum a gene's module types put it in.

        A gene whose types fall in one stratum takes that stratum's name; one
        spanning several takes "mixed". The strata partition the accepted types,
        which `Rules.load` checks, so every gene lands somewhere.
        """
        types = set(self.curated[key].module_types)
        matched = sorted(name for name, members in self.rules.strata.items()
                         if types & set(members))
        if not matched:
            raise ArchitectureError(
                f"{key} carries module types {sorted(types)}, which no declared "
                "stratum covers")
        return matched[0] if len(matched) == 1 else "mixed"

    def stratum_asserts_absence(self, key: tuple[str, str]) -> bool:
        """Whether the gene's stratum means the secondary domain is absent.

        A stratum asserts absence when every accepted type in it records the
        secondary domain on none of its modules across the whole annotation set.
        That is read from the data rather than declared, so the campaign cannot
        assert a biological fact the record does not actually carry.
        """
        return self.stratum(key) in self._absent_strata()

    def _absent_strata(self) -> frozenset[str]:
        if self._absent_cache is not None:
            return self._absent_cache
        per_type = self.annotation_census["domains_by_module_type"]
        out = []
        for name, types in sorted(self.rules.strata.items()):
            total = sum(per_type.get(module_type, {}).get(
                self.rules.secondary_domain, 0) for module_type in types)
            if total == 0:
                out.append(name)
        frozen = frozenset(out)
        object.__setattr__(self, "_absent_cache", frozen)
        return frozen

    def strata_census(self) -> dict[str, int]:
        counts = collections.Counter(self.stratum(k) for k in self.keys)
        names = [*sorted(self.rules.strata), "mixed"]
        return {name: counts.get(name, 0) for name in names if counts.get(name)}

    def module_decomposition(self,
                             hits: Mapping[tuple[str, str], Sequence[Hit]],
                             span: str | None = None,
                             policy: str | None = None) -> dict:
        """Complete declared module patterns inside each architecture.

        Counting the pattern's non-overlapping occurrences is arithmetic over
        R2's answer, and comparing that count with the curated module count is a
        different claim from comparing adenylation counts -- a gene can carry the
        right number of adenylation domains and not assemble into modules at all.

        The family list comes from the resolved hits and NOT from splitting the
        architecture string, because a Pfam family name can contain the
        separator: `AMP-binding` holds a hyphen, so splitting the hyphen-joined
        architecture shreds every family name and the pattern matches nothing.
        The string is for reading; this is the same list it was built from.
        """
        pattern = list(self.rules.module_pattern)
        rows = []
        for bgc, gene in self.keys:
            kept = self.resolve(hits.get((bgc, gene), ()), span, policy)
            families = [h.model for h in kept]
            found, i = 0, 0
            while i + len(pattern) <= len(families):
                if families[i:i + len(pattern)] == pattern:
                    found += 1
                    i += len(pattern)
                else:
                    i += 1
            curated = self.curated[(bgc, gene)].modules
            rows.append({"bgc": bgc, "gene": gene,
                         "complete_modules": found, "curated_modules": curated,
                         "verdict": _verdict(found, curated)})
        counts = collections.Counter(r["verdict"] for r in rows)
        return {"module_table": rows,
                "module_verdicts": {v: counts.get(v, 0) for v in self.rules.verdicts},
                "complete_modules_total": sum(r["complete_modules"] for r in rows)}

    # ------------------------------------------------------- the counterfactual

    def policy_sweep(self, hits_by_cutoff: Mapping[str, Mapping]) -> dict:
        """The architecture under each declared reading, against the baseline.

        Three axes, kept apart because they are independent choices: how overlaps
        are resolved, which coordinate span orders the architecture, and which
        cutoff the tool ran at. Each cell reports the same four quantities so the
        cost of a choice is directly readable.
        """
        baseline = {row["bgc"] + "|" + row["gene"]: row["architecture"]
                    for row in self.architecture_table(
                        hits_by_cutoff[self.rules.baseline_cutoff])}
        out: dict[str, dict] = {}
        for cutoff in self.rules.cutoffs:
            hits = hits_by_cutoff[cutoff]
            for policy in self.rules.resolution_policies:
                for span in self.rules.coordinate_spans:
                    table = self.architecture_table(hits, span, policy)
                    census = self.architecture_census(table)
                    agree = sum(
                        1 for row in self.reconcile(
                            hits, span=span, policy=policy)["reconciliation_table"]
                        if row["primary_verdict"] == "agrees")
                    moved = sum(
                        1 for row in table
                        if baseline[row["bgc"] + "|" + row["gene"]]
                        != row["architecture"])
                    out[f"{cutoff}|{policy}|{span}"] = {
                        "domains_total": census["domains_total"],
                        "distinct_architectures": census["distinct_architectures"],
                        "a_domain_agreement": agree,
                        "genes_differing_from_baseline": moved,
                    }
        return out

    def mapping_sweep(self,
                      hits: Mapping[tuple[str, str], Sequence[Hit]]) -> dict:
        """Agreement under each declared mapping for the PRIMARY domain type.

        The axis the catalogued design did not have and the data demands: Pfam
        splits a catalytic domain across several families -- the adenylation
        domain into an N- and a C-terminal half, the ketosynthase into three --
        so whether you count one, another or all of them decides the answer.
        Counting all of them multiplies every module, which is the kind of
        mistake that looks like diligence.

        The secondary domain's mapping is left out: it is the OTHER reconciled
        type, not an alternative reading of the primary one, so sweeping it would
        mix two axes.
        """
        out = {}
        for name in sorted(self.rules.domain_mappings):
            if name == self.rules.secondary_domain:
                continue
            block = self.reconcile(hits, name)
            out[name] = {
                "families": list(self.rules.domain_mappings[name]),
                "verdicts": block["primary_verdicts"],
            }
        return out


def _verdict(computed: int, declared: int) -> str:
    if computed == declared:
        return "agrees"
    return "computed_more" if computed > declared else "computed_fewer"


def _median(ordered: Sequence[int]) -> float:
    n = len(ordered)
    if not n:
        raise ArchitectureError("median of an empty set")
    if n % 2:
        return float(ordered[n // 2])
    return round((ordered[n // 2 - 1] + ordered[n // 2]) / 2, PRECISION)
