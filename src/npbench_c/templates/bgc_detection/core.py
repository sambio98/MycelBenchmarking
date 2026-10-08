"""BGC detection and boundary calling, reconciled against curated boundaries.

The benchmark's first S2 campaign. antiSMASH is run once over a complete,
unsliced chromosome; the regions it calls are read through a declared projection
and compared with the boundaries MIBiG records for the same accession. Then
declared core biosynthetic genes are deleted from the annotation and the search
is repeated.

Five things are load-bearing, and each was settled by measurement rather than
assumption.

**The input is a complete chromosome, not a slice.** Slicing was the obvious way
to keep the run cheap and it is wrong twice over: the slice truncates CDS features
at its edges, which antiSMASH refuses outright ("feature translation extends out
of record"), and a region that reaches the slice edge has its boundary set by the
window rather than by the tool. Measured on a 81 kb window around the actinorhodin
cluster, antiSMASH called one region running to the record end with
`contig_edge=True`. On the whole chromosome all 29 regions come back
`contig_edge=False`, so every boundary is the tool's. A minimal run over 8.67 Mb
takes about 75 seconds, so the slice bought nothing.

**The coordinate conventions differ and the difference was verified, not assumed.**
antiSMASH writes 0-based half-open locations; MIBiG's `from`/`to` are 1-based
inclusive. That is checked against the source annotation rather than inferred:
BGC0000194's `to` is 5,535,091, which is exactly where SCO5092 ends in
AL645882.2. Both sides are converted to one declared convention before anything
is compared.

**Detection is not the question; extent is.** Every one of the 15 MIBiG loci on
this accession is detected. What varies is by how much antiSMASH overcalls, so the
verdict vocabulary is about extent and the Jaccard index is graded.

**Half of what antiSMASH finds has no curated counterpart.** 15 of 29 regions
overlap no MIBiG locus, which is a statement about coverage of the curated
record and not an error in either resource.

**The perturbation targets are named in the rules, not discovered at solve time.**
Each variant lists the genes to delete. The sequence is untouched: antiSMASH reads
annotated input and scans the CDS translations it finds, so removing the feature
removes the protein from detection and perturbs nothing else.

Pure stdlib outside the tool call, sorted order throughout.
"""

from __future__ import annotations

import collections
import gzip
import json
import os
import pathlib
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.tools import cache

TEMPLATE_ID = "bgc-detection"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "detection_rules.json"
PRECISION = 6

DEFAULT_PREFIX = "/opt/npbench-tools/bin"

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 7200

#: GenBank feature table: a feature starts at column 5, qualifiers at column 21.
_FEATURE = re.compile(r"^ {5}(\S+)\s+(.*)$")


class DetectionError(RuntimeError):
    """The input or the environment does not support the detection experiment."""


def tool_prefix() -> pathlib.Path:
    return pathlib.Path(os.environ.get("NPBENCH_TOOL_PREFIX", DEFAULT_PREFIX))


def _run(argv: Sequence[str]) -> subprocess.CompletedProcess:
    """Run a pinned tool with the prefix appended to PATH.

    antiSMASH resolves nine helpers on PATH -- hmmscan, hmmsearch, hmmpress,
    hmmpfam2, blastp, makeblastdb, diamond, prodigal and FastTree -- so the prefix
    has to be there; it goes last so the tool's own interpreter cannot shadow the
    one this code runs under.
    """
    prefix = str(tool_prefix())
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if prefix in parts:
        parts.remove(prefix)
    env = {**os.environ, **FIXED_ENV, "PATH": os.pathsep.join([*parts, prefix])}
    return subprocess.run(list(argv), capture_output=True, text=True, env=env,
                          timeout=TIMEOUT_S)


# -------------------------------------------------------------- coordinates


@dataclass(frozen=True)
class Interval:
    """A half-open 0-based interval, which is the campaign's one convention."""

    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlap(self, other: "Interval") -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))

    def jaccard(self, other: "Interval") -> float:
        union = max(self.end, other.end) - min(self.start, other.start)
        return round(self.overlap(other) / union, PRECISION) if union else 0.0

    @classmethod
    def from_antismash(cls, location: str) -> "Interval":
        """antiSMASH writes '[10992:81283](+)': 0-based, half-open, as-is."""
        digits = re.findall(r"\d+", location)
        if len(digits) < 2:
            raise DetectionError(f"unparseable antiSMASH location {location!r}")
        return cls(int(digits[0]), int(digits[1]))

    @classmethod
    def from_mibig(cls, start: int, end: int) -> "Interval":
        """MIBiG writes 1-based inclusive, verified against the source record.

        BGC0000194's `to` is 5,535,091 and SCO5092 ends at 5,535,091 in
        AL645882.2, so `to` is an inclusive 1-based coordinate and `from` is too.
        Converting gives a half-open interval whose length is to - from + 1, which
        is one more than the difference MIBiG's own numbers suggest at a glance.
        """
        return cls(start - 1, end)

    def as_one_based(self) -> tuple[int, int]:
        """The same interval as 1-based inclusive, for a reader comparing to a
        GenBank record."""
        return (self.start + 1, self.end)


# -------------------------------------------------------------------- rules


@dataclass(frozen=True)
class Rules:
    accession: str
    search_flags: tuple[str, ...]
    cpus: int
    verdict_order: tuple[str, ...]
    congruent_jaccard: float
    generous_jaccard: float
    contained_coverage: float
    partial_coverage: float
    deletion_variants: Mapping[str, tuple[str, ...]]
    rule_fingerprint: str
    rule_paths: tuple[str, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        detection = spec["detection"]
        thresholds = detection["verdict_thresholds"]
        resource = spec["external_resources"]["antismash_detection_rules"]
        return cls(
            accession=detection["accession"],
            search_flags=tuple(detection["search_flags"]),
            cpus=int(detection["cpus"]),
            verdict_order=tuple(detection["verdict_order"]),
            congruent_jaccard=float(thresholds["congruent_jaccard"]),
            generous_jaccard=float(thresholds["generous_jaccard"]),
            contained_coverage=float(thresholds["contained_coverage"]),
            partial_coverage=float(thresholds["partial_coverage"]),
            deletion_variants={
                name: tuple(genes)
                for name, genes in sorted(detection["deletion_variants"].items())},
            rule_fingerprint=resource["fingerprint"],
            rule_paths=tuple(resource["relative_paths"]),
        )


# ------------------------------------------------------------- the GenBank


def read_genbank(path: pathlib.Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.read()


def delete_genes(text: str, names: Sequence[str]) -> tuple[str, list[tuple[str, str]]]:
    """Remove the CDS and gene features for the named genes.

    The sequence is left alone. antiSMASH scans the CDS translations present in an
    annotated record, so removing the feature is what "delete the gene" means for
    a detection run -- and it keeps every other coordinate in the record where it
    was, which cutting the nucleotides out would not.
    """
    wanted = set(names)
    if not wanted:
        return text, []
    lines = text.splitlines(keepends=True)
    out: list[str] = []
    removed: list[tuple[str, str]] = []
    index, in_features = 0, False
    while index < len(lines):
        line = lines[index]
        if line.startswith("FEATURES"):
            in_features = True
            out.append(line)
            index += 1
            continue
        if in_features and (line.startswith("ORIGIN") or line.startswith("//")):
            in_features = False
            out.append(line)
            index += 1
            continue
        match = _FEATURE.match(line) if in_features else None
        if not match:
            out.append(line)
            index += 1
            continue
        block = [line]
        cursor = index + 1
        while cursor < len(lines):
            nxt = lines[cursor]
            if (_FEATURE.match(nxt) or nxt.startswith("ORIGIN")
                    or nxt.startswith("//")):
                break
            block.append(nxt)
            cursor += 1
        body = "".join(block)
        hit = next((n for n in sorted(wanted) if f'/gene="{n}"' in body), None)
        if hit and match.group(1) in ("CDS", "gene"):
            removed.append((match.group(1), hit))
        else:
            out.extend(block)
        index = cursor
    found = {name for _, name in removed}
    missing = wanted - found
    if missing:
        raise DetectionError(
            f"the declared deletion names {sorted(missing)} and the record carries "
            "no such gene; a variant that deletes nothing is not a counterfactual")
    return "".join(out), removed


# ------------------------------------------------------------------ regions


@dataclass(frozen=True)
class Region:
    number: int
    interval: Interval
    products: tuple[str, ...]
    contig_edge: bool

    def key(self) -> tuple:
        return (self.interval.start, self.interval.end, self.products)


@dataclass(frozen=True)
class Locus:
    bgc: str
    interval: Interval
    completeness: str
    classes: tuple[str, ...]


@dataclass(frozen=True)
class Detection:
    accession: str
    genbank_path: pathlib.Path
    loci: tuple[Locus, ...]
    rules: Rules
    record_length: int

    @classmethod
    def load(cls, campaign: pathlib.Path, genbank: str, loci_file: str) -> "Detection":
        campaign = pathlib.Path(campaign)
        rules = Rules.load(campaign)
        path = campaign / genbank
        with (gzip.open(path, "rt") if path.suffix == ".gz" else open(path)) as handle:
            header = handle.readline()
        match = re.match(r"^LOCUS\s+\S+\s+(\d+) bp", header)
        if not match:
            raise DetectionError(f"unreadable GenBank LOCUS line: {header[:70]!r}")
        payload = json.loads((campaign / loci_file).read_text())
        if payload["accession"] != rules.accession:
            raise DetectionError(
                f"the loci table is for {payload['accession']} and the rules "
                f"declare {rules.accession}")
        loci = tuple(sorted(
            (Locus(bgc=item["bgc"],
                   interval=Interval.from_mibig(item["from"], item["to"]),
                   completeness=item["completeness"],
                   classes=tuple(item["classes"])) for item in payload["loci"]),
            key=lambda l: (l.interval.start, l.interval.end, l.bgc)))
        if not loci:
            raise DetectionError("the loci table is empty, so R3 has nothing to "
                                 "reconcile against")
        return cls(accession=rules.accession, genbank_path=path, loci=loci,
                   rules=rules, record_length=int(match.group(1)))

    # ------------------------------------------------------------ tool calls

    def antismash_version(self) -> str:
        done = _run([str(tool_prefix() / "antismash"), "--version"])
        if done.returncode != 0:
            raise DetectionError(f"antismash --version failed: {done.stderr[-400:]}")
        return done.stdout.strip()

    def rule_fingerprint(self) -> str:
        """A hash over the three detection rule files, as the registry declares.

        The rule set IS the gold for a region-detection campaign, so a campaign
        pins what it was built against and a rule change is detectable even when
        the version string is not what moved.
        """
        import hashlib

        root = tool_prefix().parent / "lib"
        digest = hashlib.sha256()
        for relative in self.rules.rule_paths:
            matches = sorted(root.glob(f"python3.*/site-packages/antismash/{relative}"))
            if not matches:
                raise DetectionError(
                    f"antiSMASH detection rule file {relative} is not present "
                    "under the pinned image; this campaign reads it from the "
                    "image rather than shipping it, because it is AGPL package "
                    "data whose licence the benchmark cannot pass on")
            digest.update(matches[-1].read_bytes())
        return digest.hexdigest()[:16]

    def _detect(self, deleted: Sequence[str]) -> list[dict]:
        """One antiSMASH run, projected to the regions it called."""
        text = read_genbank(self.genbank_path)
        text, _ = delete_genes(text, deleted)
        with tempfile.TemporaryDirectory(prefix="npbench-bgcdetect-") as tmp:
            work = pathlib.Path(tmp)
            record = work / "record.gbk"
            record.write_text(text)
            out = work / "out"
            done = _run([str(tool_prefix() / "antismash"),
                         *self.rules.search_flags,
                         "--cpus", str(self.rules.cpus),
                         "--output-dir", str(out),
                         "--output-basename", "record",
                         "--logfile", str(work / "run.log"),
                         str(record)])
            artifact = out / "record.json"
            if done.returncode != 0 or not artifact.is_file():
                tail = (done.stdout or "")[-500:] + (done.stderr or "")[-500:]
                raise DetectionError(f"antismash failed: {tail}")
            document = json.loads(artifact.read_text())
        records = document.get("records")
        if not isinstance(records, list) or len(records) != 1:
            raise DetectionError(
                f"expected one record in the antiSMASH output, got "
                f"{len(records) if isinstance(records, list) else 'none'}")
        projected = []
        for feature in records[0].get("features", []):
            if feature.get("type") != "region":
                continue
            qualifiers = feature.get("qualifiers", {})
            interval = Interval.from_antismash(feature["location"])
            projected.append({
                "number": int(qualifiers["region_number"][0]),
                "start": interval.start,
                "end": interval.end,
                "products": sorted(qualifiers.get("product", [])),
                "contig_edge": qualifiers.get("contig_edge", ["False"])[0] == "True",
            })
        if not projected:
            raise DetectionError(
                "antiSMASH called no region at all on this record, which is a "
                "broken run rather than a finding")
        return sorted(projected, key=lambda r: (r["start"], r["end"]))

    def regions(self, deleted: Sequence[str] = ()) -> tuple[Region, ...]:
        """The declared projection of one antiSMASH run.

        Memoised only where a tool memo is enabled, which is stub sweeps alone.
        The key carries the input digest, the tool version, the detection rule
        fingerprint, the flags and the deleted gene list -- a memo that dropped
        any of those would answer from the wrong world.
        """
        deleted = tuple(sorted(deleted))

        def compute() -> list[dict]:
            return self._detect(deleted)

        if cache.enabled():
            payload = cache.memoize_json(
                ["bgc_detection.regions",
                 cache.file_digest(self.genbank_path),
                 self.antismash_version(),
                 self.rule_fingerprint(),
                 ",".join(self.rules.search_flags) or "none",
                 str(self.rules.cpus),
                 ",".join(deleted) or "none"],
                compute)
        else:
            payload = compute()
        return tuple(Region(number=item["number"],
                            interval=Interval(item["start"], item["end"]),
                            products=tuple(item["products"]),
                            contig_edge=bool(item["contig_edge"]))
                     for item in payload)

    # ------------------------------------------------------------- inventory

    def region_census(self, regions: Sequence[Region]) -> dict:
        products: collections.Counter = collections.Counter()
        for region in regions:
            products.update(region.products)
        covered = sum(r.interval.length for r in regions)
        return {
            "regions": len(regions),
            "distinct_products": len(products),
            "product_census": dict(sorted(products.items())),
            "contig_edge_regions": sum(1 for r in regions if r.contig_edge),
            "region_bases": covered,
            "record_length": self.record_length,
            "fraction_of_record": round(covered / self.record_length, PRECISION),
        }

    def region_table(self, regions: Sequence[Region]) -> list[dict]:
        """Every region, in one declared convention and the reader's convention."""
        rows = []
        for region in sorted(regions, key=lambda r: (r.interval.start, r.interval.end)):
            one_based = region.interval.as_one_based()
            rows.append({
                "number": region.number,
                "start": region.interval.start,
                "end": region.interval.end,
                "start_one_based": one_based[0],
                "end_one_based": one_based[1],
                "length": region.interval.length,
                "products": list(region.products),
                "contig_edge": region.contig_edge,
            })
        return rows

    def span_summary(self, regions: Sequence[Region]) -> dict:
        lengths = sorted(r.interval.length for r in regions)
        if not lengths:
            raise DetectionError("no region to summarise")
        middle = len(lengths) // 2
        median = (float(lengths[middle]) if len(lengths) % 2
                  else round((lengths[middle - 1] + lengths[middle]) / 2, PRECISION))
        return {"min": lengths[0], "median": median, "max": lengths[-1]}

    # --------------------------------------------------------- reconciliation

    def best_region(self, locus: Locus,
                    regions: Sequence[Region]) -> tuple[Region | None, int]:
        """The region overlapping this locus most, under a declared tie-break.

        Ties go to the lower start then the lower end, so a locus sitting in two
        equally overlapping regions does not depend on feature order.
        """
        best, best_overlap = None, 0
        for region in sorted(regions,
                             key=lambda r: (r.interval.start, r.interval.end)):
            overlap = locus.interval.overlap(region.interval)
            if overlap > best_overlap:
                best, best_overlap = region, overlap
        return best, best_overlap

    def verdict(self, locus: Locus, regions: Sequence[Region]) -> dict:
        """The declared verdict for one curated locus against the called regions.

        First match wins over the declared order. The vocabulary is about EXTENT,
        not detection: on this accession every locus is detected, and what varies
        is how much larger the called region is.
        """
        region, overlap = self.best_region(locus, regions)
        if region is None:
            return {"bgc": locus.bgc, "verdict": "not_detected", "region": None,
                    "coverage": 0.0, "jaccard": 0.0}
        coverage = round(overlap / locus.interval.length, PRECISION)
        jaccard = locus.interval.jaccard(region.interval)
        rules = self.rules
        if coverage >= rules.contained_coverage:
            if jaccard >= rules.congruent_jaccard:
                verdict = "congruent"
            elif jaccard >= rules.generous_jaccard:
                verdict = "contained_generous"
            else:
                verdict = "contained_much_larger"
        elif coverage >= rules.partial_coverage:
            verdict = "partial"
        else:
            verdict = "marginal"
        return {"bgc": locus.bgc, "verdict": verdict, "region": region.number,
                "coverage": coverage, "jaccard": jaccard}

    def reconciliation(self, regions: Sequence[Region]) -> dict:
        # Ordered by MIBiG accession, which is a declared key the measurer can
        # reproduce without knowing any coordinate. Coordinate order would read
        # more naturally and would make the row list depend on something the
        # measurer must not need.
        rows = sorted((self.verdict(locus, regions) for locus in self.loci),
                      key=lambda row: row["bgc"])
        counts = collections.Counter(row["verdict"] for row in rows)
        matched = {row["region"] for row in rows if row["region"] is not None}
        shared: collections.Counter = collections.Counter(
            row["region"] for row in rows if row["region"] is not None)
        jaccards = sorted(row["jaccard"] for row in rows)
        middle = len(jaccards) // 2
        median = (jaccards[middle] if len(jaccards) % 2
                  else round((jaccards[middle - 1] + jaccards[middle]) / 2, PRECISION))
        return {
            "loci_total": len(self.loci),
            "verdict_census": {v: counts.get(v, 0) for v in self.rules.verdict_order},
            "per_locus": rows,
            "regions_matching_a_locus": len(matched),
            "regions_without_a_locus": len(regions) - len(matched),
            "regions_holding_several_loci": sorted(
                number for number, n in shared.items() if n > 1),
            "jaccard_median": median,
            "jaccard_min": jaccards[0],
            "jaccard_max": jaccards[-1],
        }

    # ------------------------------------------------------- counterfactual

    def deletion_outcome(self, baseline: Sequence[Region],
                         deleted: Sequence[str]) -> dict:
        """What one declared deletion does to the called region set."""
        after = self.regions(deleted)
        before_keys = {r.key() for r in baseline}
        after_keys = {r.key() for r in after}
        lost = sorted(before_keys - after_keys)
        gained = sorted(after_keys - before_keys)
        return {
            "deleted_genes": sorted(deleted),
            "regions": len(after),
            "regions_lost": len(lost),
            "regions_gained": len(gained),
            "regions_unchanged": len(before_keys & after_keys),
            "lost_products": sorted({p for _, _, products in lost for p in products}),
            "gained_products": sorted({p for _, _, products in gained
                                       for p in products}),
        }
