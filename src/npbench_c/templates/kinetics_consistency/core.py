"""Reconcile BRENDA's three kinetic constants against each other.

BRENDA stores Km, kcat and kcat/Km as independent records, and the third is
arithmetically determined by the first two. So the corpus can be checked against
itself -- no external standard, no curator, no judgement -- and the shape of the
disagreement is the finding.

The units are the whole difficulty. BRENDA's JSON records none: a value is a
string like `0.05 {benzyl alcohol}`, a number and the substrate in braces, and the
unit is a property of the FIELD. Km is in mM, kcat in s^-1, kcat/Km in
mM^-1 s^-1. A system that does not know that cannot relate the three at all, and
one that assumes micromolar is wrong by exactly a factor of 1000 -- which is one
of the bands this campaign measures, because the disagreements in the real data
cluster there and at 60.

Pure stdlib, deterministic, sorted order throughout.
"""

from __future__ import annotations

import json
import pathlib
import re
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

TEMPLATE_ID = "kinetics-consistency"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "kinetic_units.json"

#: The three linked fields, in the order the unit table declares them.
FIELDS = ("km_value", "turnover_number", "kcat_km_value")

#: A number, optionally in exponent form, optionally followed by {substrate}.
VALUE = re.compile(r"^\s*(-?\d+(?:\.\d+)?(?:[eE]-?\d+)?)\s*(?:\{(.*)\})?\s*$")

PRECISION = 6


class KineticsError(RuntimeError):
    """The corpus does not support the reconciliation. A campaign bug."""


@dataclass(frozen=True)
class Band:
    name: str
    low: float
    high: float
    reciprocal: bool

    def matches(self, ratio: float) -> bool:
        if self.low <= ratio <= self.high:
            return True
        return bool(self.reciprocal and ratio > 0
                    and self.low <= 1.0 / ratio <= self.high)


@dataclass(frozen=True)
class Rules:
    units: Mapping[str, str]
    sentinel: str
    bands: tuple[Band, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        bands = tuple(
            Band(b["name"], float(b["low"]), float(b["high"]),
                 bool(b.get("reciprocal")))
            for b in spec["consistency_rule"]["agreement_bands"])
        return cls(
            units={f: spec["fields"][f]["unit"] for f in FIELDS},
            sentinel=str(spec["identity_rule"]["missing_value_sentinel"]),
            bands=bands,
        )

    def band_of(self, ratio: float) -> str:
        """First matching band wins; the order is declared, not incidental."""
        for band in self.bands:
            if band.matches(ratio):
                return band.name
        return "other"


def parse_value(raw: object, sentinel: str) -> tuple[float | None, str, str]:
    """-> (number or None, substrate, outcome).

    outcome is one of parsed / sentinel / unparsed, and nothing is coerced: a
    range like '0.05 - 0.1' is unparsed and counted, never averaged into a number
    nobody wrote down.
    """
    text = raw if isinstance(raw, str) else ""
    match = VALUE.match(text)
    if match is None:
        return None, "", "unparsed"
    substrate = (match.group(2) or "").strip()
    number = float(match.group(1))
    # Numeric comparison, not a string prefix: '-999', '-999.0' and '-9.99e2' are
    # the same sentinel, and '-9990' is not.
    if number == float(sentinel):
        return None, substrate, "sentinel"
    return number, substrate, "parsed"


@dataclass(frozen=True)
class Kinetics:
    release: str
    corpus: Mapping[str, Mapping]
    rules: Rules

    @classmethod
    def load(cls, campaign: pathlib.Path, relative: str) -> "Kinetics":
        import gzip

        path = pathlib.Path(campaign) / relative
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            payload = json.load(handle)
        return cls(release=payload["release"], corpus=payload["data"],
                   rules=Rules.load(pathlib.Path(campaign)))

    # ---------------------------------------------------------------- census

    def census(self) -> dict:
        """Per field: records, parsed, sentinels, unparsed. Never summed across."""
        records = {f: 0 for f in FIELDS}
        parsed = {f: 0 for f in FIELDS}
        sentinel = {f: 0 for f in FIELDS}
        unparsed = {f: 0 for f in FIELDS}
        for ec in sorted(self.corpus):
            for field in FIELDS:
                for entry in (self.corpus[ec].get(field) or []):
                    records[field] += 1
                    _, _, outcome = parse_value(entry.get("value"),
                                                self.rules.sentinel)
                    {"parsed": parsed, "sentinel": sentinel,
                     "unparsed": unparsed}[outcome][field] += 1
        return {
            "records_by_field": records,
            "parsed_by_field": parsed,
            "sentinel_by_field": sentinel,
            "unparsed_by_field": unparsed,
        }

    # ---------------------------------------------------------------- triples

    def triples(self) -> dict[tuple, dict[str, float]]:
        """(ec, proteins, substrate, references) -> the parsed values it carries.

        The key has to be this tight. Dropping the reference would pair values
        from different papers; dropping the protein would pair different enzymes.
        Either way the disagreement rate would stop meaning anything.
        """
        table: dict[tuple, dict[str, float]] = {}
        for ec in sorted(self.corpus):
            for field in FIELDS:
                for entry in (self.corpus[ec].get(field) or []):
                    number, substrate, outcome = parse_value(
                        entry.get("value"), self.rules.sentinel)
                    if outcome != "parsed":
                        continue
                    key = (ec, tuple(entry.get("proteins") or []), substrate,
                           tuple(entry.get("references") or []))
                    table.setdefault(key, {})[field] = number
        return table

    def completeness(self, table: Mapping[tuple, Mapping[str, float]]) -> dict:
        counts = {"1": 0, "2": 0, "3": 0}
        for values in table.values():
            counts[str(len(values))] += 1
        return counts

    # ------------------------------------------------------------ consistency

    def consistency(self, table: Mapping[tuple, Mapping[str, float]]) -> dict:
        """The ratio of computed to stored catalytic efficiency, banded."""
        bands = {band.name: 0 for band in self.rules.bands}
        bands["other"] = 0
        excluded = 0
        complete = 0
        ratios: list[float] = []
        for key in sorted(table):
            values = table[key]
            if len(values) != 3:
                continue
            complete += 1
            km, stored = values["km_value"], values["kcat_km_value"]
            if km <= 0 or stored <= 0:
                excluded += 1
                continue
            ratio = (values["turnover_number"] / km) / stored
            ratios.append(ratio)
            bands[self.rules.band_of(ratio)] += 1
        if not ratios:
            raise KineticsError("no usable complete triple; nothing to reconcile")
        return {
            "complete_triples": complete,
            "excluded_non_positive": excluded,
            "checked": len(ratios),
            "band_counts": dict(sorted(bands.items())),
            "agreement_rate": round(bands["agree"] / len(ratios), PRECISION),
        }

    def band_examples(self, table: Mapping[tuple, Mapping[str, float]],
                      band_name: str, limit: int = 5) -> list[dict]:
        """The first few triples in a band, by sorted key. For audit, not grading."""
        out = []
        for key in sorted(table):
            values = table[key]
            if len(values) != 3:
                continue
            km, stored = values["km_value"], values["kcat_km_value"]
            if km <= 0 or stored <= 0:
                continue
            ratio = (values["turnover_number"] / km) / stored
            if self.rules.band_of(ratio) != band_name:
                continue
            out.append({
                "ec": key[0], "substrate": key[2],
                "km_mM": values["km_value"],
                "kcat_per_s": values["turnover_number"],
                "stored_kcat_km": stored,
                "computed_kcat_km": round(values["turnover_number"] / km, PRECISION),
                "ratio": round(ratio, PRECISION),
            })
            if len(out) >= limit:
                break
        return out

    def under_unit_error(self, table: Mapping[tuple, Mapping[str, float]],
                         km_factor: float, kcat_factor: float) -> dict:
        """The same check with a declared unit mistake applied.

        The point of the counterfactual: if the implicit convention did not
        matter, scaling Km by a thousand would leave the agreement rate alone.
        """
        agree = checked = 0
        for key in sorted(table):
            values = table[key]
            if len(values) != 3:
                continue
            km = values["km_value"] * km_factor
            stored = values["kcat_km_value"]
            if km <= 0 or stored <= 0:
                continue
            ratio = ((values["turnover_number"] * kcat_factor) / km) / stored
            checked += 1
            if self.rules.band_of(ratio) == "agree":
                agree += 1
        return {"checked": checked, "agree": agree,
                "agreement_rate": round(agree / checked, PRECISION) if checked else 0.0}
