"""Reconcile two MIBiG releases across a schema migration.

MIBiG 4.0 restructured the schema: 3.1 nests everything under a ``cluster``
object, 4.0 is flat with twelve top-level keys. A raw field-by-field diff
therefore reports almost every field of almost every entry as changed, which is
true and useless. What is well-posed is reconciliation on a declared field
mapping, separating a change in *representation* from a change in *substance*.

The sharpest instance is in the corpus rather than planted: 3.1 stores the NCBI
taxon id as a JSON string and 4.0 as an integer. Compared raw, all 2442 shared
entries look re-annotated; compared as strings, 8 genuine corrections remain.

Pure stdlib, deterministic.
"""

from __future__ import annotations

import json
import pathlib
import tarfile
from dataclasses import dataclass
from typing import Any, Iterable, Mapping

TEMPLATE_ID = "mibig-release-diff"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"


class DiffError(RuntimeError):
    """The corpora do not support the reconciliation. A campaign bug if it
    fires on shipped inputs."""


def load_release(archive: pathlib.Path) -> dict[str, dict]:
    """Read every BGC record from a release tarball, keyed by accession.

    Members are sorted before reading so the mapping is built in a fixed order
    and nothing depends on the archive's internal ordering.
    """
    out: dict[str, dict] = {}
    with tarfile.open(archive, "r:gz") as tar:
        members = sorted((m for m in tar.getmembers()
                          if m.isfile() and m.name.endswith(".json")
                          and pathlib.PurePath(m.name).name.startswith("BGC")),
                         key=lambda m: m.name)
        for member in members:
            handle = tar.extractfile(member)
            if handle is None:
                continue
            accession = pathlib.PurePath(member.name).stem
            out[accession] = json.loads(handle.read().decode("utf-8"))
    if not out:
        raise DiffError(f"no BGC records found in {archive.name}")
    return out


@dataclass(frozen=True)
class ReleaseDiff:
    fields: Mapping[str, Mapping]
    verdict_vocabulary: tuple[str, ...]
    verdict_precedence: tuple[str, ...]
    class_mapping: Mapping[str, str | None]
    retired_terms: tuple[str, ...]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "ReleaseDiff":
        reference = pathlib.Path(campaign) / "reference"
        mapping = json.loads((reference / "field_mapping.json").read_text())
        classes = json.loads((reference / "class_mapping.json").read_text())
        return cls(
            fields=mapping["fields"],
            verdict_vocabulary=tuple(mapping["verdict_vocabulary"]),
            verdict_precedence=tuple(mapping["verdict_precedence"]),
            class_mapping=classes["mapping"],
            retired_terms=tuple(classes["retired_terms"]),
        )

    # ------------------------------------------------------------ extraction

    @staticmethod
    def _dig(record: Mapping, path: str) -> Any:
        """Resolve a declared path. ``[]`` collects a field across a list."""
        cur: Any = record
        for part in path.split("."):
            if part.endswith("[]"):
                key = part[:-2]
                if not isinstance(cur, Mapping) or key not in cur:
                    return None
                cur = cur[key]
                if not isinstance(cur, list):
                    return None
                return cur
            if not isinstance(cur, Mapping) or part not in cur:
                return None
            cur = cur[part]
        return cur

    def extract(self, record: Mapping, path: str) -> Any:
        """Extract a declared path's value, handling the ``a[].b`` shape."""
        if "[]." in path:
            head, tail = path.split("[].", 1)
            items = self._dig(record, head + "[]")
            if not isinstance(items, list):
                return None
            values = [i.get(tail) for i in items if isinstance(i, Mapping)]
            return [v for v in values if v is not None]
        if path.endswith("[]"):
            return self._dig(record, path)
        return self._dig(record, path)

    # --------------------------------------------------------- normalisation

    def normalise(self, value: Any, rule: str) -> Any:
        if value is None:
            return None
        if rule == "string":
            return str(value)
        if rule == "class_vocabulary":
            terms = value if isinstance(value, list) else [value]
            return sorted({self.class_mapping[t] for t in terms
                           if self.class_mapping.get(t) is not None})
        if isinstance(value, list):
            return sorted({v for v in value if v is not None})
        return value

    def _present(self, value: Any) -> bool:
        if value is None:
            return False
        if isinstance(value, (list, set, dict, str)):
            return len(value) > 0
        return True

    # ----------------------------------------------------------- comparison

    def classify_field(self, field: str, source: Mapping, target: Mapping) -> str:
        """Classify one field of one entry, by the declared precedence."""
        spec = self.fields[field]
        raw_s = self.extract(source, spec["source_path"])
        raw_t = self.extract(target, spec["target_path"])
        rule = spec.get("normalisation", "none")

        if rule == "class_vocabulary":
            terms = raw_s if isinstance(raw_s, list) else ([raw_s] if raw_s else [])
            if any(t in self.retired_terms for t in terms):
                return "no_counterpart_in_target"
        if not self._present(raw_s):
            return "absent_in_source"
        if not self._present(raw_t):
            return "absent_in_target"

        if spec["kind"] == "scalar_in_set":
            targets = raw_t if isinstance(raw_t, list) else [raw_t]
            return "unchanged" if raw_s in targets else "substantive_change"

        if raw_s == raw_t:
            return "unchanged"
        norm_s = self.normalise(raw_s, rule)
        norm_t = self.normalise(raw_t, "none" if rule == "class_vocabulary" else rule)
        if norm_s == norm_t:
            return "representation_only"
        return "substantive_change"

    def raw_differs(self, field: str, source: Mapping, target: Mapping) -> bool:
        """Whether the field differs before any normalisation.

        This is what a direct JSON comparison reports, and the gap between it
        and the classified result is the point of the campaign: for the taxon
        id it is 2442 against 8.
        """
        spec = self.fields[field]
        raw_s = self.extract(source, spec["source_path"])
        raw_t = self.extract(target, spec["target_path"])
        if spec["kind"] == "scalar_in_set":
            targets = raw_t if isinstance(raw_t, list) else [raw_t]
            return raw_s not in targets
        if isinstance(raw_s, list) or isinstance(raw_t, list):
            a = sorted(raw_s) if isinstance(raw_s, list) else [raw_s]
            b = sorted(raw_t) if isinstance(raw_t, list) else [raw_t]
            return a != b
        return raw_s != raw_t

    # ------------------------------------------------------------- corpus

    def corpus_diff(self, source: Mapping[str, dict],
                    target: Mapping[str, dict]) -> dict:
        """Entry-set diff, status transitions, raw and classified field counts."""
        s_keys, t_keys = set(source), set(target)
        shared = sorted(s_keys & t_keys)

        status: dict[str, int] = {}
        for accession in shared:
            value = target[accession].get("status") or "unknown"
            status[value] = status.get(value, 0) + 1

        raw_counts: dict[str, int] = {}
        verdicts: dict[str, dict[str, int]] = {}
        for field in sorted(self.fields):
            raw = 0
            tally: dict[str, int] = {v: 0 for v in self.verdict_vocabulary}
            for accession in shared:
                s, t = source[accession], target[accession]
                if self.raw_differs(field, s, t):
                    raw += 1
                tally[self.classify_field(field, s, t)] += 1
            raw_counts[field] = raw
            verdicts[field] = {k: v for k, v in sorted(tally.items()) if v}

        return {
            "source_entries": len(s_keys),
            "target_entries": len(t_keys),
            "only_in_source": len(s_keys - t_keys),
            "only_in_target": len(t_keys - s_keys),
            "shared_entries": len(shared),
            "shared_status_in_target": dict(sorted(status.items())),
            "raw_difference_counts": dict(sorted(raw_counts.items())),
            "classified_verdict_counts": verdicts,
        }
