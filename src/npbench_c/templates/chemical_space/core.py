"""Chemical-space bookkeeping over a pinned chemistry table.

The campaign's subject is the identity key. A corpus of natural products has no
single right answer to "how many compounds is this?": counting full InChIKeys,
connectivity blocks and scaffolds gives three different numbers over the same
records, and every downstream statistic -- duplicates, per-class spread,
scaffold sharing -- inherits the choice. So the keys are declared in
space_rules.json and the ladder grades the consequences.

Scaffold perception needs a chemistry toolkit, so it happens once in
build_reference.py and is pinned into reference/chemistry_table.json. Everything
here is stdlib set arithmetic over that table joined to the release archive, in
sorted order throughout, which is what keeps the campaign S0 and the numbers
reproducible.

One rule in here exists because of what the data turned out to look like. A
Bemis-Murcko scaffold keeps ring systems and discards substituents, so over this
corpus the most widely shared scaffold is plain benzene, in 76 entries. The
informativeness filter is a declared response to that, and both the filtered and
unfiltered counts are graded so the size of the effect stays visible.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

TEMPLATE_ID = "chemical-space"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
TABLE = "chemistry_table.json"
RULES = "space_rules.json"

#: The three declared identity keys, in increasing coarseness.
KEYS = ("compound", "connectivity", "scaffold")


class SpaceError(RuntimeError):
    """The corpus does not support the analysis. A campaign bug on shipped inputs."""


@dataclass(frozen=True)
class Rules:
    status: str
    class_vocabulary: tuple[str, ...]
    accepted_evidence: frozenset[str]
    min_ring_count: int
    min_coverage: float

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads(
            (pathlib.Path(campaign) / "reference" / RULES).read_text())
        informative = spec["informative_scaffold_rule"]
        return cls(
            status=spec["corpus_filter"]["status"],
            class_vocabulary=tuple(spec["cross_class_rule"]["class_vocabulary"]),
            accepted_evidence=frozenset(
                spec["evidence_gate"]["accepted_locus_evidence"]),
            min_ring_count=int(informative["min_ring_count"]),
            min_coverage=float(informative["min_heavy_atom_coverage"]),
        )


def connectivity_key(inchikey: str) -> str:
    """The InChIKey's first block: skeleton and connectivity, nothing else."""
    return inchikey.split("-", 1)[0]


@dataclass(frozen=True)
class ChemicalSpace:
    release: str
    records: tuple[dict, ...]
    corpus_entries: int
    active_entries: int
    classes: Mapping[str, frozenset[str]]
    locus_evidence: Mapping[str, frozenset[str]]
    rules: Rules

    @classmethod
    def load(cls, campaign: pathlib.Path, archive: str) -> "ChemicalSpace":
        from npbench_c.templates.mibig_diff.core import load_release

        campaign = pathlib.Path(campaign)
        rules = Rules.load(campaign)
        table = json.loads((campaign / "reference" / TABLE).read_text())
        corpus = load_release(campaign / archive)

        active = {a: e for a, e in corpus.items() if e.get("status") == rules.status}
        classes, evidence = {}, {}
        for accession, entry in active.items():
            declared = {c.get("class")
                        for c in ((entry.get("biosynthesis") or {}).get("classes") or [])
                        if c.get("class")}
            unknown = declared - set(rules.class_vocabulary)
            if unknown:
                raise SpaceError(
                    f"{accession}: classes outside the declared vocabulary: "
                    f"{sorted(unknown)}")
            classes[accession] = frozenset(declared)
            evidence[accession] = frozenset(
                e.get("method") for locus in (entry.get("loci") or [])
                for e in (locus.get("evidence") or [])
                if isinstance(e, dict) and e.get("method"))

        records = tuple(sorted(table["records"], key=lambda r: r["record_id"]))
        stray = sorted({r["accession"] for r in records} - set(active))
        if stray:
            raise SpaceError(
                f"the chemistry table carries records for entries that are not "
                f"{rules.status} in this archive: {stray[:5]}")
        return cls(release=table["corpus_release"], records=records,
                   corpus_entries=len(corpus), active_entries=len(active),
                   classes=classes, locus_evidence=evidence, rules=rules)

    # ------------------------------------------------------------------ keys

    def key(self, record: Mapping, which: str) -> str | None:
        if which == "compound":
            return record["inchikey"]
        if which == "connectivity":
            return connectivity_key(record["inchikey"])
        if which == "scaffold":
            return record["scaffold_inchikey"]
        raise SpaceError(f"{which!r} is not a declared identity key")

    def distinct(self, records: Sequence[Mapping], which: str) -> int:
        return len({k for r in records if (k := self.key(r, which)) is not None})

    def entries_by_key(self, records: Sequence[Mapping],
                       which: str) -> dict[str, list[str]]:
        """Key -> the distinct accessions carrying it, sorted."""
        out: dict[str, set[str]] = {}
        for record in records:
            k = self.key(record, which)
            if k is not None:
                out.setdefault(k, set()).add(record["accession"])
        return {k: sorted(v) for k, v in sorted(out.items())}

    def cross_entry_duplicates(self, records: Sequence[Mapping], which: str) -> int:
        return sum(1 for accessions in self.entries_by_key(records, which).values()
                   if len(accessions) > 1)

    # ---------------------------------------------------------- informative

    def informative(self, records: Sequence[Mapping], min_rings: int | None = None,
                    min_coverage: float | None = None) -> tuple[dict, ...]:
        """Records whose scaffold clears the declared ring and coverage floors."""
        rings = self.rules.min_ring_count if min_rings is None else min_rings
        coverage = self.rules.min_coverage if min_coverage is None else min_coverage
        out = []
        for record in records:
            if record["scaffold_inchikey"] is None:
                continue
            heavy = record["compound_heavy_atoms"]
            if not heavy:
                raise SpaceError(f"{record['record_id']} has no heavy atoms")
            if record["scaffold_ring_count"] < rings:
                continue
            if record["scaffold_heavy_atoms"] / heavy < coverage:
                continue
            out.append(record)
        return tuple(out)

    # --------------------------------------------------------- class joins

    def shared_scaffolds(self, records: Sequence[Mapping]) -> dict[str, list[str]]:
        return {k: v for k, v in self.entries_by_key(records, "scaffold").items()
                if len(v) > 1}

    def cross_class(self, records: Sequence[Mapping]) -> dict:
        """Scaffolds spanning classes, under the declared rule and the naive one.

        Declared: two of the scaffold's entries have disjoint class sets, so a
        hybrid cluster cannot make a scaffold cross classes on its own. Naive:
        the union of the entries' classes has more than one name, which is what
        most analyses compute and which 546 multi-class entries inflate.
        """
        declared: dict[str, list[str]] = {}
        naive: list[str] = []
        for scaffold, accessions in self.shared_scaffolds(records).items():
            sets = [self.classes.get(a, frozenset()) for a in accessions]
            if len({c for s in sets for c in s}) > 1:
                naive.append(scaffold)
            if any(not (a & b) for i, a in enumerate(sets) for b in sets[i + 1:]):
                declared[scaffold] = accessions
        return {"declared": declared, "naive_count": len(naive)}

    def passes_evidence(self, accession: str) -> bool:
        return bool(self.locus_evidence.get(accession, frozenset())
                    & self.rules.accepted_evidence)

    def evidence_gated(self, records: Sequence[Mapping]) -> tuple[dict, ...]:
        return tuple(r for r in records if self.passes_evidence(r["accession"]))

    # --------------------------------------------------------------- census

    def per_class(self, records: Sequence[Mapping], which: str | None) -> dict:
        """Per class, the record count or the distinct-key count.

        A record contributes to every class its entry carries, so these sum to
        more than the corpus total; the alternative would need a primary-class
        rule MIBiG does not supply.
        """
        buckets: dict[str, set | int] = {}
        for name in self.rules.class_vocabulary:
            buckets[name] = set() if which else 0
        for record in records:
            for name in sorted(self.classes.get(record["accession"], ())):
                if which is None:
                    buckets[name] += 1                      # type: ignore[operator]
                else:
                    k = self.key(record, which)
                    if k is not None:
                        buckets[name].add(k)                # type: ignore[union-attr]
        return {n: (len(v) if which else v) for n, v in sorted(buckets.items())}

    def ring_histogram(self, records: Sequence[Mapping]) -> dict[str, int]:
        """Scaffold ring counts, bucketed. Acyclic compounds are bucket '0'."""
        counts: dict[int, int] = {}
        for record in records:
            counts[record["scaffold_ring_count"]] = (
                counts.get(record["scaffold_ring_count"], 0) + 1)
        return {str(k): counts[k] for k in sorted(counts)}

    def multi_key_connectivity_groups(self, records: Sequence[Mapping]) -> int:
        """Connectivity blocks covering more than one full InChIKey.

        These are the stereoisomers, charge states and labelled variants that a
        connectivity-level count merges. If this were zero the two granularities
        would be the same statistic and R2 would have nothing to say.
        """
        groups: dict[str, set[str]] = {}
        for record in records:
            groups.setdefault(connectivity_key(record["inchikey"]),
                              set()).add(record["inchikey"])
        return sum(1 for v in groups.values() if len(v) > 1)
