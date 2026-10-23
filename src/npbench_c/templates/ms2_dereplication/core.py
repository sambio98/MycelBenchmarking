"""MS2 dereplication against a reference library: the declared pipeline.

T-L4-1. A query spectrum is matched against a library of spectra of known
compounds, and the question is which compound the query is. What makes it
gradeable is the split: the query's own compound is in the library only through
spectra from DIFFERENT contributors, so a correct answer is another lab's
measurement rather than a lookup of the same one. That is a rule the inputs are
built under and a refusal `Corpus.load` enforces, because it first shipped as a
docstring and the audit found 10 of 90 correct answers coming from the query's
own contributor at the query's own instrument type. The INSTRUMENT type is not
guaranteed per query and the campaign does not claim it; `inputs/provenance.json`
reports the instrument support as a census.

Pure stdlib apart from the tool call, which is a subprocess: matchms is a
library pin and importing it would put numpy, scipy and numba inside a grading
core that is stdlib-only by design.

Three things are load-bearing and each was measured before it was pinned.

**The precursor filter is part of the pipeline, not an optimisation.** Candidates
are library spectra whose precursor m/z is within a declared window of the
query's -- which is what dereplication practice does, and which is what an
adduct misassignment breaks. Median candidate set: 20 spectra out of 30,749.

**The tolerance is part of gold.** Across the declared grid the aggregate top-1
rate moves one point (0.877 to 0.890 over 300 queries) and 26 of 300 individual
answers move somewhere in it. A campaign that picked one tolerance and called
the rest a detail would be grading a parameter it never declared.

**The adduct misassignment makes the search silent, not wrong.** Re-reading a
protonated precursor as sodiated empties the candidate window for 246 of 300
queries. So it is graded as a census -- how many annotations go missing, how
many go wrong, how many survive -- and not as a per-query survival bit, which
would have had one answer.
"""

from __future__ import annotations

import gzip
import json
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

from npbench_c.tools import cache
from npbench_c.tools.resources import tool_prefix

TEMPLATE_ID = "ms2-dereplication"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "dereplication_rules.json"
DRIVER = "derep_driver.py"
PRECISION = 6

#: The interpreter that has matchms. Not `python3`: the image's PATH keeps the
#: tool prefix last so the benchmark's own interpreter stays the system one.
INTERPRETER = "python3.11"


def declared_ranking(row: Sequence) -> tuple:
    """Score down, matched peaks down, accession up.

    The third key is not cosmetic: equal scores occur, and two keys would leave
    the answer to whatever order the tool returned. Measured on the shipped
    hundred, neither the second nor the third key changes a single top match --
    8 queries have tied candidate scores and none of them at the top -- so the
    rule is a determinism guarantee rather than a discriminator, and the stub
    ladder's probe gets the order wrong in a way that does bite.
    """
    return (-row[1], -row[2], row[0])


class DereplicationError(RuntimeError):
    """The corpus or the tool cannot support what the campaign grades. A
    campaign bug if it fires on shipped inputs."""


@dataclass(frozen=True)
class Rules:
    similarity: str
    peak_tolerance: float
    precursor_window: float
    tolerance_grid: tuple[float, ...]
    confidence_rule: tuple[tuple[float, int, str], ...]
    confidence_fallback: str
    confidence_vocabulary: tuple[str, ...]
    ambiguity_margin: float
    proton_mass: float
    sodium_mass: float
    assumed_adduct: str
    misread_adduct: str
    tool_version: str

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        document = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        search = document["search"]
        confidence = document["confidence"]
        adduct = document["adduct_misassignment"]
        return cls(
            similarity=str(search["similarity"]),
            peak_tolerance=float(search["peak_tolerance"]),
            precursor_window=float(search["precursor_window"]),
            tolerance_grid=tuple(float(t) for t in search["tolerance_grid"]),
            confidence_rule=tuple(
                (float(floor), int(peaks), str(label))
                for floor, peaks, label in confidence["rule"]),
            confidence_fallback=str(confidence["fallback"]),
            confidence_vocabulary=tuple(confidence["vocabulary"]),
            ambiguity_margin=float(document["ambiguity"]["margin"]),
            proton_mass=float(adduct["proton_mass"]),
            sodium_mass=float(adduct["sodium_mass"]),
            assumed_adduct=str(adduct["assumed"]),
            misread_adduct=str(adduct["misread_as"]),
            tool_version=str(search["tool_version"]),
        )

    def classify(self, score: float, matches: int) -> str:
        """The declared confidence class of one annotation."""
        for floor, peaks, label in self.confidence_rule:
            if score >= floor and matches >= peaks:
                return label
        return self.confidence_fallback


def _read_gz(path: pathlib.Path) -> list[dict]:
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


@dataclass(frozen=True)
class Corpus:
    queries: tuple[dict, ...]
    library: tuple[dict, ...]
    rules: Rules
    campaign: pathlib.Path

    @classmethod
    def load(cls, campaign: pathlib.Path, queries_input: str,
             library_input: str) -> "Corpus":
        campaign = pathlib.Path(campaign)
        queries = _read_gz(campaign / queries_input)
        library = _read_gz(campaign / library_input)
        if not queries or not library:
            raise DereplicationError(
                f"{queries_input} holds {len(queries)} queries and "
                f"{library_input} holds {len(library)} spectra; both are needed")
        accessions = {row["accession"] for row in library}
        overlap = sorted(row["accession"] for row in queries
                         if row["accession"] in accessions)
        if overlap:
            raise DereplicationError(
                f"{len(overlap)} query spectra are also in the library "
                f"({overlap[:3]}). A query that can match itself makes the "
                "retrieval a lookup and every rung above R1 meaningless")
        by_block: dict[str, list[dict]] = {}
        for row in library:
            by_block.setdefault(row["compound_block"], []).append(row)
        for row in queries:
            own = by_block.get(row["compound_block"], [])
            if not own:
                raise DereplicationError(
                    f"{row['accession']}'s compound is not in the library, so "
                    "its correct answer does not exist")
            # The split's claim, checked on every load rather than asserted in a
            # docstring: a library spectrum of the query's compound from the
            # query's own contributor makes a correct answer a replicate lookup.
            same = sorted(other["accession"] for other in own
                          if other["contributor"] == row["contributor"])
            if same:
                raise DereplicationError(
                    f"{row['accession']} is answerable from its own "
                    f"contributor's spectra ({same[:3]}), so a correct retrieval "
                    "would be the same lab's replicate rather than another lab's "
                    "measurement")
        return cls(queries=tuple(sorted(queries, key=lambda r: r["accession"])),
                   library=tuple(sorted(library, key=lambda r: r["accession"])),
                   rules=Rules.load(campaign), campaign=campaign)

    # ------------------------------------------------------------- tool calls

    def matchms_version(self) -> str:
        """The library's own version string, via the interpreter that has it."""
        command = [str(tool_prefix() / INTERPRETER), "-c",
                   "import matchms; print(matchms.__version__)"]

        def query() -> str:
            done = subprocess.run(command, capture_output=True, text=True,
                                  timeout=600)
            if done.returncode != 0:
                raise DereplicationError(
                    f"matchms version query failed: {done.stderr[-400:]}")
            return done.stdout.strip()

        return cache.memoize_version([*command, self.rules.tool_version], query)

    def candidates(self, query: Mapping, centre: float | None = None,
                   window: float | None = None) -> list[dict]:
        """Library spectra whose precursor sits within the declared window."""
        centre = query["precursor_mz"] if centre is None else centre
        window = self.rules.precursor_window if window is None else window
        return [row for row in self.library
                if abs(row["precursor_mz"] - centre) <= window]

    def _score(self, tasks: Sequence[Mapping], tolerance: float) -> dict[str, list]:
        """One driver call for every (query, candidate set) pair.

        Memoised on the tool version, the tolerance and a digest of the job, so
        a stub sweep does not re-run the same comparison once per rung-repeat.
        """
        job = {"peak_tolerance": tolerance, "tasks": list(tasks)}
        body = json.dumps(job, sort_keys=True)
        key = ["ms2_dereplication.score", self.matchms_version(),
               f"tolerance={tolerance}", cache.text_digest(body)]

        def compute() -> dict:
            with tempfile.TemporaryDirectory() as tmp:
                work = pathlib.Path(tmp)
                (work / "job.json").write_text(body)
                done = subprocess.run(
                    [str(tool_prefix() / INTERPRETER), str(CONSTANTS / DRIVER),
                     str(work / "job.json"), str(work / "out.json")],
                    capture_output=True, text=True, timeout=3600)
                if done.returncode != 0:
                    raise DereplicationError(
                        f"the matchms driver failed: {done.stderr[-600:]}")
                return json.loads((work / "out.json").read_text())

        return cache.memoize_json(key, compute)

    # --------------------------------------------------------------- the search

    def rank(self, tolerance: float | None = None, window: float | None = None,
             centres: Mapping[str, float] | None = None,
             ranking: "Callable[[list], tuple] | None" = None) -> dict[str, list]:
        """Every query's candidates, ranked by the declared rule.

        `centres` re-points the precursor filter, which is how the adduct
        misassignment is applied: the spectra do not change, the window does.

        `ranking` exists for one caller: the stub ladder's competence probe, which
        has to get the declared order wrong in a way that is internally
        consistent. Gold never passes it, and `declared_ranking` is the default.
        """
        tolerance = self.rules.peak_tolerance if tolerance is None else tolerance
        tasks = []
        for query in self.queries:
            centre = (centres or {}).get(query["accession"])
            tasks.append({"query": query,
                          "candidates": self.candidates(query, centre, window)})
        scored = self._score(tasks, tolerance)
        ranked: dict[str, list] = {}
        for query in self.queries:
            rows = scored.get(query["accession"]) or []
            ranked[query["accession"]] = sorted(
                ([str(name), float(score), int(matches)]
                 for name, score, matches in rows),
                key=ranking or declared_ranking)
        return ranked

    # ---------------------------------------------------------- the annotations

    def _block_of(self, accession: str) -> str:
        for row in self.library:
            if row["accession"] == accession:
                return row["compound_block"]
        raise DereplicationError(f"{accession} is not a library spectrum")

    def annotations(self, tolerance: float | None = None,
                    window: float | None = None,
                    centres: Mapping[str, float] | None = None,
                    ranking: "Callable[[list], tuple] | None" = None) -> list[dict]:
        """One annotation per query: what it matched and how well."""
        ranked = self.rank(tolerance, window, centres, ranking)
        out = []
        for query in self.queries:
            rows = ranked[query["accession"]]
            record = {"query": query["accession"],
                      "candidates": len(self.candidates(
                          query, (centres or {}).get(query["accession"]), window))}
            if not rows:
                record.update({
                    "verdict": "no_candidate" if not record["candidates"]
                               else "no_match",
                    "top_match": None, "score": None, "matched_peaks": None,
                    "confidence": "no_candidate" if not record["candidates"]
                                  else "no_match",
                    "correct": False, "margin_to_other_compound": None})
                out.append(record)
                continue
            name, score, matches = rows[0]
            block = self._block_of(name)
            competitor = next(
                (row for row in rows[1:] if self._block_of(row[0]) != block), None)
            record.update({
                "verdict": "annotated",
                "top_match": name,
                "top_block": block,
                "score": round(score, PRECISION),
                "matched_peaks": matches,
                "confidence": self.rules.classify(score, matches),
                "correct": block == query["compound_block"],
                "margin_to_other_compound": (None if competitor is None
                                             else round(score - competitor[1],
                                                        PRECISION)),
            })
            out.append(record)
        return out
