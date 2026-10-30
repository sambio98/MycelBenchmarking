"""Selectivity ratios from measured activities, and nothing else.

L5's rule is that gold here is measured-or-computed only: a value from a database
with a structured confidence field, or exact arithmetic over such values. Never a
model prediction, never a curated expert call. This engine is the specific form of
that rule -- every admission test reads a structured ChEMBL field, and the only
arithmetic is a median and a quotient.

Two things are load-bearing and easy to get wrong. ChEMBL's confidence_score lives
on the ASSAY, not on the activity, so admission needs a join. And a '>' or '<'
relation is a bound rather than a measurement: treating one as a value is the
commonest way a potency table acquires numbers nobody measured.

Pure stdlib, deterministic, sorted order throughout.
"""

from __future__ import annotations

import gzip
import json
import pathlib
from dataclasses import dataclass
from typing import Mapping, Sequence

TEMPLATE_ID = "chembl-selectivity"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "selectivity_rules.json"
PRECISION = 6


class SelectivityError(RuntimeError):
    """The activity set does not support the comparison. A campaign bug."""


@dataclass(frozen=True)
class Rules:
    standard_type: str
    rejection_order: tuple[str, ...]
    min_confidence: int
    thresholds: tuple[int, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        admission = spec["admission"]
        return cls(
            standard_type=admission["standard_type"],
            rejection_order=tuple(r["reason"]
                                  for r in admission["ordered_rejection_reasons"]),
            min_confidence=int(admission["min_confidence_score"]),
            thresholds=tuple(int(t) for t in spec["selectivity"]["thresholds"]),
        )


def median(values: Sequence[float]) -> float:
    """Declared: the middle value, or the mean of the two middle ones."""
    ordered = sorted(values)
    n = len(ordered)
    if not n:
        raise SelectivityError("median of an empty measurement set")
    if n % 2:
        return round(ordered[n // 2], PRECISION)
    return round((ordered[n // 2 - 1] + ordered[n // 2]) / 2, PRECISION)


def minimum(values: Sequence[float]) -> float:
    return round(min(values), PRECISION)


AGGREGATORS = {"median": median, "minimum": minimum}


@dataclass(frozen=True)
class Admission:
    """Which admission rules are in force. An R4 variant relaxes one of them."""

    require_exact_relation: bool = True
    require_nanomolar: bool = True
    reject_flagged: bool = True
    reject_duplicates: bool = True
    aggregator: str = "median"

    @classmethod
    def from_spec(cls, spec: Mapping) -> "Admission":
        return cls(
            require_exact_relation=bool(spec.get("require_exact_relation", True)),
            require_nanomolar=bool(spec.get("require_nanomolar", True)),
            reject_flagged=bool(spec.get("reject_flagged", True)),
            reject_duplicates=bool(spec.get("reject_duplicates", True)),
            aggregator=str(spec.get("aggregator", "median")),
        )


@dataclass(frozen=True)
class Selectivity:
    release: str
    focal_target: str
    reference_target: str
    activities: tuple[dict, ...]
    assays: Mapping[str, dict]
    rules: Rules

    @classmethod
    def load(cls, campaign: pathlib.Path, relative: str, focal: str,
             reference: str) -> "Selectivity":
        path = pathlib.Path(campaign) / relative
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            payload = json.load(handle)
        activities = tuple(sorted(payload["activities"],
                                  key=lambda a: a["activity_id"]))
        targets = {a["target_chembl_id"] for a in activities}
        missing = {focal, reference} - targets
        if missing:
            raise SelectivityError(
                f"the shipped activities carry nothing for {sorted(missing)}")
        return cls(release=payload["release"], focal_target=focal,
                   reference_target=reference, activities=activities,
                   assays=payload["assays"], rules=Rules.load(campaign))

    # ------------------------------------------------------------- admission

    def reject_reason(self, activity: Mapping,
                      admission: Admission) -> str | None:
        """The FIRST declared reason this activity fails, or None.

        The order is declared so the rejection counts partition the rejected set
        exactly; without it an activity failing two rules would be countable
        twice and the breakdown would not sum.
        """
        if admission.require_exact_relation and activity.get("standard_relation") != "=":
            return "censored_or_missing_relation"
        if admission.require_nanomolar and activity.get("standard_units") != "nM":
            return "unconvertible_units"
        value = activity.get("standard_value")
        if value in (None, ""):
            return "no_value"
        if float(value) <= 0:
            return "non_positive"
        if admission.reject_flagged and activity.get("data_validity_comment"):
            return "flagged_by_chembl"
        if admission.reject_duplicates and activity.get("potential_duplicate"):
            return "potential_duplicate"
        assay = self.assays.get(activity.get("assay_chembl_id") or "", {})
        if (assay.get("confidence_score") or 0) < self.rules.min_confidence:
            return "low_confidence_assay"
        return None

    def of_type(self) -> tuple[dict, ...]:
        return tuple(a for a in self.activities
                     if a.get("standard_type") == self.rules.standard_type)

    def admission_census(self, admission: Admission) -> dict:
        """Per target: admitted, and rejected broken down by first reason."""
        targets = (self.focal_target, self.reference_target)
        admitted = {t: 0 for t in targets}
        rejected = {t: {r: 0 for r in self.rules.rejection_order} for t in targets}
        for activity in self.of_type():
            target = activity["target_chembl_id"]
            if target not in admitted:
                continue
            reason = self.reject_reason(activity, admission)
            if reason is None:
                admitted[target] += 1
            else:
                rejected[target][reason] += 1
        return {
            "admitted_by_target": dict(sorted(admitted.items())),
            "rejected_by_target": {t: dict(sorted(r.items()))
                                   for t, r in sorted(rejected.items())},
        }

    # ----------------------------------------------------------- aggregation

    def aggregated(self, admission: Admission) -> dict[tuple[str, str], float]:
        """(molecule, target) -> one value, under the declared aggregator."""
        buckets: dict[tuple[str, str], list[float]] = {}
        for activity in self.of_type():
            if activity["target_chembl_id"] not in (self.focal_target,
                                                    self.reference_target):
                continue
            if self.reject_reason(activity, admission) is not None:
                continue
            key = (activity["molecule_chembl_id"], activity["target_chembl_id"])
            buckets.setdefault(key, []).append(float(activity["standard_value"]))
        if admission.aggregator not in AGGREGATORS:
            raise SelectivityError(
                f"{admission.aggregator!r} is not a declared aggregator")
        reduce = AGGREGATORS[admission.aggregator]
        return {k: reduce(v) for k, v in sorted(buckets.items())}

    # ----------------------------------------------------------- selectivity

    def ratios(self, admission: Admission) -> dict[str, float]:
        """Molecule -> reference IC50 over focal IC50, for molecules on both."""
        values = self.aggregated(admission)
        molecules = sorted({m for m, t in values if t == self.focal_target}
                           & {m for m, t in values if t == self.reference_target})
        out = {}
        for molecule in molecules:
            focal = values[(molecule, self.focal_target)]
            if focal <= 0:
                continue
            out[molecule] = round(
                values[(molecule, self.reference_target)] / focal, PRECISION)
        if not out:
            raise SelectivityError(
                "no molecule carries an admissible value on both targets")
        return out

    def verdicts(self, ratios: Mapping[str, float]) -> dict:
        """Counts and committed sets at each declared threshold."""
        focal_selective, reference_selective = {}, {}
        for threshold in self.rules.thresholds:
            focal_selective[str(threshold)] = sorted(
                m for m, r in ratios.items() if r > threshold)
            reference_selective[str(threshold)] = sorted(
                m for m, r in ratios.items() if r < 1.0 / threshold)
        ranked = sorted(((-r, m) for m, r in ratios.items()))
        return {
            "eligible_molecules": len(ratios),
            "focal_selective_counts": {k: len(v) for k, v in focal_selective.items()},
            "reference_selective_counts": {k: len(v)
                                           for k, v in reference_selective.items()},
            "focal_selective_at_highest_threshold":
                focal_selective[str(max(self.rules.thresholds))],
            "most_selective_molecule": ranked[0][1],
            "most_selective_ratio": round(-ranked[0][0], PRECISION),
        }
