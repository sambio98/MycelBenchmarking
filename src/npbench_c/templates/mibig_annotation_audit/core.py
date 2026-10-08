"""Audit how well MIBiG's gene-function annotations are backed by evidence.

This is the benchmark's own G2 doctrine turned on a database: a functional claim
counts only when a structured evidence field says so. MIBiG records three
independent evidence axes, and conflating them is the error the campaign
targets:

  loci[].evidence       how the cluster was linked to its product
  compounds[].evidence  how the compound's structure was determined
  functions[].evidence  how a gene's function was established

Admissibility requires the function axis AND the locus axis AND active status,
which is why the admissible count is far below the merely-backed count.

Pure stdlib, deterministic.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass
from typing import Mapping

TEMPLATE_ID = "mibig-annotation-audit"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"


class AuditError(RuntimeError):
    """The corpus does not support the audit. A campaign bug on shipped inputs."""


@dataclass(frozen=True)
class Policy:
    name: str
    function_evidence: frozenset[str]
    locus_evidence: frozenset[str]
    require_active: bool

    @classmethod
    def from_spec(cls, name: str, spec: Mapping) -> "Policy":
        return cls(
            name=name,
            function_evidence=frozenset(spec["function_evidence_accepted"]),
            locus_evidence=frozenset(spec["locus_evidence_accepted"]),
            require_active=bool(spec.get("require_active_status", True)),
        )


@dataclass(frozen=True)
class AnnotationAudit:
    function_vocabulary: tuple[str, ...]
    function_methods: tuple[str, ...]
    baseline: Policy

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "AnnotationAudit":
        spec = json.loads(
            (pathlib.Path(campaign) / "reference" / "evidence_policy.json").read_text())
        return cls(
            function_vocabulary=tuple(spec["function_vocabulary"]),
            function_methods=tuple(spec["function_evidence_methods"]),
            baseline=Policy.from_spec(spec["baseline_policy"]["name"],
                                      spec["baseline_policy"]),
        )

    # ---------------------------------------------------------- extraction

    @staticmethod
    def pairs(entry: Mapping):
        """Yield (function_name, accepted_methods) per (annotation, function).

        The unit is the pair, not the gene: a gene carrying two functions
        contributes two, which the counting rules state so a correct answer
        cannot differ on the convention.
        """
        for annotation in ((entry.get("genes") or {}).get("annotations") or []):
            for function in (annotation.get("functions") or []):
                if not isinstance(function, dict):
                    continue
                name = ((function.get("function") or {}) or {}).get("name")
                methods = {e.get("method") for e in (function.get("evidence") or [])
                           if isinstance(e, dict)}
                yield name, {m for m in methods if m}

    @staticmethod
    def locus_methods(entry: Mapping) -> set[str]:
        return {e.get("method")
                for locus in (entry.get("loci") or [])
                for e in (locus.get("evidence") or [])
                if isinstance(e, dict) and e.get("method")}

    # ------------------------------------------------------------- census

    def census(self, corpus: Mapping[str, dict]) -> dict:
        """Backed / unbacked counts per function category, plus method totals."""
        backed = {f: 0 for f in self.function_vocabulary}
        unbacked = {f: 0 for f in self.function_vocabulary}
        methods = {m: 0 for m in self.function_methods}
        entries_with_annotations = 0
        total_pairs = 0
        unknown: set[str] = set()

        for accession in sorted(corpus):
            entry = corpus[accession]
            seen = False
            for name, evidence in self.pairs(entry):
                seen = True
                total_pairs += 1
                if name not in backed:
                    unknown.add(str(name))
                    continue
                if evidence:
                    backed[name] += 1
                else:
                    unbacked[name] += 1
                for m in sorted(evidence):
                    if m in methods:
                        methods[m] += 1
            if seen:
                entries_with_annotations += 1

        if unknown:
            raise AuditError(
                f"function names outside the declared vocabulary: {sorted(unknown)}")
        return {
            "entries_with_annotations": entries_with_annotations,
            "total_function_annotations": total_pairs,
            "backed_by_function": dict(sorted(backed.items())),
            "unbacked_by_function": dict(sorted(unbacked.items())),
            "function_evidence_method_totals": dict(sorted(methods.items())),
        }

    def admissible(self, corpus: Mapping[str, dict], policy: Policy) -> dict:
        """Counts admissible under a policy, plus the accessions per category.

        All three conditions apply -- function evidence, locus evidence, active
        status -- so this is strictly smaller than the backed count, and the gap
        is what distinguishes "a curator recorded a function" from "the function
        is experimentally supported for a cluster whose product link is itself
        experimentally supported".
        """
        counts = {f: 0 for f in self.function_vocabulary}
        accessions: dict[str, set[str]] = {f: set() for f in self.function_vocabulary}
        for accession in sorted(corpus):
            entry = corpus[accession]
            if policy.require_active and entry.get("status") != "active":
                continue
            if not (self.locus_methods(entry) & policy.locus_evidence):
                continue
            for name, evidence in self.pairs(entry):
                if name in counts and (evidence & policy.function_evidence):
                    counts[name] += 1
                    accessions[name].add(accession)
        return {
            "policy": policy.name,
            "admissible_by_function": dict(sorted(counts.items())),
            "admissible_entries_by_function": {
                f: sorted(a) for f, a in sorted(accessions.items())
            },
        }
