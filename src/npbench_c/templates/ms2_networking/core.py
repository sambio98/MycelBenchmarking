"""MS2 molecular networking over one laboratory's spectra: the declared pipeline.

T-L4-2. A set of spectra acquired together is turned into a graph -- nodes are
spectra, an edge joins two whose cosine clears a declared threshold with a
declared number of matched peaks -- and the graph is then annotated from a
reference library of OTHER laboratories' spectra. Most nodes get no annotation of
their own, and what the network does with the few that do is the campaign.

Pure stdlib apart from the tool call, which is a subprocess: matchms is a library
pin and importing it would put numpy, scipy and numba inside a grading core that
is stdlib-only by design.

Two tool calls serve every rung. The pairwise matrix is computed once at the
lowest declared threshold, so the whole R4 grid is a filter over one call; the
annotation layer is a second call over the per-node candidate windows.

Three things are load-bearing.

**The edge PAIRS are graded and the edge SCORES are not.** A cosine is a float
out of a BLAS-backed library. The pair set is a comparison against a declared
threshold, and measured on this node set no pair sits within 1e-4 of any declared
threshold except one, 5.6e-05 below 0.9 -- three orders of magnitude above the
float noise this image's pinned numpy and scipy produce.

**A node carries no identity.** The compound a node is, is the answer. The node
records carry what an experimental spectrum carries and `load` refuses one that
carries a compound block, an InChIKey or a formula; the identities live in
`gold/nodes_truth.json`, outside every sandbox. So `correct` is a gold-side
column and nothing graded asks whether an annotation is right.

**Propagation is a claim about the evidence, not about the truth.** An
unannotated node joined to an annotated one inherits either `same_compound` or
`analogue_with_mass_shift`, decided by the precursor difference. The grounding
measured both against the truth, and the asymmetry it found is the campaign's
subject rather than a graded component.
"""

from __future__ import annotations

import json
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.templates.ms2_dereplication.core import (
    IDENTITY_FIELDS,
    INTERPRETER,
    declared_ranking,
)
from npbench_c.tools import cache
from npbench_c.tools.resources import tool_prefix

TEMPLATE_ID = "ms2-networking"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "networking_rules.json"
DRIVER = "network_driver.py"
PRECISION = 6

#: The nodes' own identities, outside every sandbox.
TRUTH = "gold/nodes_truth.json"


class NetworkError(RuntimeError):
    """The corpus or the tool cannot support what the campaign grades. A
    campaign bug if it fires on shipped inputs."""


@dataclass(frozen=True)
class Rules:
    similarity: str
    peak_tolerance: float
    edge_score_threshold: float
    edge_min_matched_peaks: int
    threshold_grid: tuple[float, ...]
    precursor_window: float
    ranking_rule: str
    confidence_rule: tuple[tuple[float, int, str], ...]
    confidence_fallback: str
    confidence_vocabulary: tuple[str, ...]
    accepted: tuple[str, ...]
    propagation_vocabulary: tuple[str, ...]
    mass_shift_precision: int
    tool_version: str

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        document = json.loads(
            (pathlib.Path(campaign) / "reference" / RULES).read_text())
        network = document["network"]
        annotation = document["annotation"]
        confidence = annotation["confidence"]
        propagation = document["propagation"]
        rules = cls(
            similarity=str(network["similarity"]),
            peak_tolerance=float(network["peak_tolerance"]),
            edge_score_threshold=float(network["edge_score_threshold"]),
            edge_min_matched_peaks=int(network["edge_min_matched_peaks"]),
            threshold_grid=tuple(float(t) for t in network["threshold_grid"]),
            precursor_window=float(annotation["precursor_window"]),
            ranking_rule=str(annotation["ranking_rule"]),
            confidence_rule=tuple((float(floor), int(peaks), str(label))
                                  for floor, peaks, label in confidence["rule"]),
            confidence_fallback=str(confidence["fallback"]),
            confidence_vocabulary=tuple(confidence["vocabulary"]),
            accepted=tuple(annotation["accepted"]),
            propagation_vocabulary=tuple(propagation["vocabulary"]),
            mass_shift_precision=int(propagation["mass_shift_precision"]),
            tool_version=str(network["tool_version"]),
        )
        if rules.edge_score_threshold not in rules.threshold_grid:
            raise NetworkError(
                f"the declared edge threshold {rules.edge_score_threshold} is "
                f"not in the declared grid {rules.threshold_grid}, so R4 has no "
                "identity rung and every movement in it is uncalibrated")
        unknown = [label for label in rules.accepted
                   if label not in rules.confidence_vocabulary]
        if unknown:
            raise NetworkError(
                f"the accepted tiers {unknown} are not in the declared "
                "confidence vocabulary")
        return rules

    def classify(self, score: float, matches: int) -> str:
        """The declared confidence class of one annotation."""
        for floor, peaks, label in self.confidence_rule:
            if score >= floor and matches >= peaks:
                return label
        return self.confidence_fallback


def _read_gz(path: pathlib.Path) -> list[dict]:
    import gzip

    with gzip.open(path, "rt") as handle:
        return json.load(handle)


@dataclass(frozen=True)
class Network:
    nodes: tuple[dict, ...]
    library: tuple[dict, ...]
    rules: Rules
    campaign: pathlib.Path
    #: accession -> compound block, or None where the truth is not readable.
    truth: Mapping[str, str] | None = None

    @classmethod
    def load(cls, campaign: pathlib.Path, nodes_input: str,
             library_input: str) -> "Network":
        campaign = pathlib.Path(campaign)
        nodes = _read_gz(campaign / nodes_input)
        library = _read_gz(campaign / library_input)
        if not nodes or not library:
            raise NetworkError(
                f"{nodes_input} holds {len(nodes)} nodes and {library_input} "
                f"holds {len(library)} spectra; both are needed")

        leaked = sorted({field for row in nodes for field in IDENTITY_FIELDS
                         if field in row})
        if leaked:
            raise NetworkError(
                f"the node records carry {leaked}, which name the compound each "
                "node is -- the answer. A node carries what an experimental "
                "spectrum carries; the identities belong in " + TRUTH)

        contributors = {row["contributor"] for row in nodes}
        if len(contributors) != 1:
            raise NetworkError(
                f"the node set spans {sorted(contributors)}; a network is one "
                "laboratory's spectra, and mixing laboratories mixes their "
                "acquisition into every edge weight")
        instruments = {row["instrument_type"] for row in nodes}
        if len(instruments) != 1:
            raise NetworkError(
                f"the node set spans {sorted(instruments)}; fragmentation "
                "behaviour differs between instrument types, so an edge between "
                "two of them is not the similarity the campaign declares")
        shared = sorted(contributors & {row["contributor"] for row in library})
        if shared:
            raise NetworkError(
                f"the library carries the node contributor's own spectra "
                f"({shared}), so a node could be annotated from its own "
                "laboratory's replicate rather than another's measurement")
        overlap = sorted({row["accession"] for row in nodes}
                         & {row["accession"] for row in library})
        if overlap:
            raise NetworkError(
                f"{len(overlap)} node spectra are also in the library "
                f"({overlap[:3]}); a node that can match itself makes its "
                "annotation a lookup")

        truth = None
        truth_path = campaign / TRUTH
        if truth_path.is_file():
            truth = json.loads(truth_path.read_text())
            missing = sorted(row["accession"] for row in nodes
                             if row["accession"] not in truth)
            if missing:
                raise NetworkError(
                    f"{len(missing)} nodes have no entry in {TRUTH} "
                    f"({missing[:3]}), so their answer is undefined")
        return cls(nodes=tuple(sorted(nodes, key=lambda r: r["accession"])),
                   library=tuple(sorted(library, key=lambda r: r["accession"])),
                   rules=Rules.load(campaign), campaign=campaign, truth=truth)

    # ------------------------------------------------------------- tool calls

    def matchms_version(self) -> str:
        """The library's own version string, via the interpreter that has it."""
        command = [str(tool_prefix() / INTERPRETER), "-s", "-c",
                   "import matchms; print(matchms.__version__)"]

        def query() -> str:
            done = subprocess.run(command, capture_output=True, text=True,
                                  timeout=600)
            if done.returncode != 0:
                raise NetworkError(
                    f"matchms version query failed: {done.stderr[-400:]}")
            return done.stdout.strip()

        return cache.memoize_version([*command, self.rules.tool_version], query)

    def _drive(self, job: dict, label: str) -> object:
        """One driver call, memoised on the tool version and the job's digest."""
        body = json.dumps(job, sort_keys=True)
        key = [f"ms2_networking.{label}", self.matchms_version(),
               cache.text_digest(body)]

        def compute() -> object:
            with tempfile.TemporaryDirectory() as tmp:
                work = pathlib.Path(tmp)
                (work / "job.json").write_text(body)
                done = subprocess.run(
                    # `-s` keeps the surrounding container's user site out of
                    # the import path, so the scores come from the pinned prefix
                    # or not at all.
                    [str(tool_prefix() / INTERPRETER), "-s",
                     str(CONSTANTS / DRIVER),
                     str(work / "job.json"), str(work / "out.json")],
                    capture_output=True, text=True, timeout=7200)
                if done.returncode != 0:
                    raise NetworkError(
                        f"the matchms driver failed: {done.stderr[-600:]}")
                return json.loads((work / "out.json").read_text())

        return cache.memoize_json(key, compute)

    # --------------------------------------------------------------- the graph

    def pairs(self) -> list[list]:
        """Every node pair at or above the LOWEST declared threshold.

        One call serves the whole R4 grid: a higher threshold is a filter over
        this list, and the scores do not depend on which threshold is applied.
        """
        job = {"kind": "pairs", "peak_tolerance": self.rules.peak_tolerance,
               "score_floor": min(self.rules.threshold_grid),
               "spectra": [dict(row) for row in self.nodes]}
        rows = self._drive(job, "pairs")
        return sorted(([str(a), str(b), float(score), int(matches)]
                       for a, b, score, matches in rows),
                      key=lambda row: (row[0], row[1]))

    def edges(self, threshold: float | None = None) -> list[list]:
        """The declared edge rule applied at one threshold."""
        threshold = (self.rules.edge_score_threshold if threshold is None
                     else threshold)
        return [[a, b, round(score, PRECISION), matches]
                for a, b, score, matches in self.pairs()
                if score >= threshold
                and matches >= self.rules.edge_min_matched_peaks]

    def components(self, threshold: float | None = None) -> list[list[str]]:
        """Connected components, largest first and then by first accession."""
        parent = {row["accession"]: row["accession"] for row in self.nodes}

        def find(node: str) -> str:
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node

        for a, b, *_ in self.edges(threshold):
            left, right = find(a), find(b)
            if left != right:
                parent[left] = right
        groups: dict[str, list[str]] = {}
        for row in self.nodes:
            groups.setdefault(find(row["accession"]), []).append(row["accession"])
        return sorted((sorted(members) for members in groups.values()),
                      key=lambda members: (-len(members), members[0]))

    # ---------------------------------------------------------- the annotation

    def candidates(self, node: Mapping) -> list[dict]:
        """Library spectra whose precursor sits inside the declared window."""
        window = self.rules.precursor_window
        return [row for row in self.library
                if abs(row["precursor_mz"] - node["precursor_mz"]) <= window]

    def _scored_windows(self) -> dict[str, list]:
        job = {"kind": "windows", "peak_tolerance": self.rules.peak_tolerance,
               "tasks": [{"query": dict(node),
                          "candidates": [dict(row) for row in self.candidates(node)]}
                         for node in self.nodes]}
        return self._drive(job, "windows")

    def _block_of(self, accession: str) -> str:
        for row in self.library:
            if row["accession"] == accession:
                return row["compound_block"]
        raise NetworkError(f"{accession} is not a library spectrum")

    def annotations(self) -> list[dict]:
        """One row per node: what the library says it is, and how well."""
        scored = self._scored_windows()
        out = []
        for node in self.nodes:
            accession = node["accession"]
            rows = sorted(([str(name), float(score), int(matches)]
                           for name, score, matches in scored.get(accession) or []),
                          key=declared_ranking)
            record = {"node": accession, "candidates": len(self.candidates(node))}
            if not rows:
                record.update({
                    "verdict": "no_candidate" if not record["candidates"]
                               else "no_match",
                    "top_match": None, "block": None, "score": None,
                    "matched_peaks": None,
                    "confidence": "no_candidate" if not record["candidates"]
                                  else "no_match",
                    "correct": None if self.truth is None else False})
                out.append(record)
                continue
            name, score, matches = rows[0]
            block = self._block_of(name)
            record.update({
                "verdict": "annotated", "top_match": name, "block": block,
                "score": round(score, PRECISION), "matched_peaks": matches,
                "confidence": self.rules.classify(score, matches),
                "correct": (None if self.truth is None
                            else block == self.truth[accession])})
            out.append(record)
        return out

    def accepted(self, rows: Sequence[Mapping]) -> dict[str, str]:
        """Node -> compound block, for the tiers the rules accept as a seed."""
        return {row["node"]: row["block"] for row in rows
                if row["confidence"] in self.rules.accepted}

    # --------------------------------------------------------- the propagation

    def claims(self, rows: Sequence[Mapping], threshold: float | None = None,
               seeds: Mapping[str, str] | None = None) -> dict[str, dict]:
        """The declared propagation claim for every node.

        `seeds` is the accepted annotation layer, passed in so that R4 can run
        the same rule over a smaller one.
        """
        seeds = self.accepted(rows) if seeds is None else seeds
        precursor = {row["accession"]: row["precursor_mz"] for row in self.nodes}
        adjacency: dict[str, list[tuple[str, float]]] = {}
        for a, b, score, _matches in self.edges(threshold):
            adjacency.setdefault(a, []).append((b, score))
            adjacency.setdefault(b, []).append((a, score))

        out: dict[str, dict] = {}
        for node in self.nodes:
            accession = node["accession"]
            if accession in seeds:
                out[accession] = {"claim": "annotated_directly",
                                  "block": seeds[accession], "mass_shift": None}
                continue
            near = [(other, score) for other, score
                    in sorted(adjacency.get(accession, ()))
                    if other in seeds]
            if not near:
                out[accession] = {"claim": "no_annotated_neighbour",
                                  "block": None, "mass_shift": None}
                continue
            blocks = {seeds[other] for other, _ in near}
            if len(blocks) > 1:
                out[accession] = {"claim": "conflicting_neighbours",
                                  "block": None, "mass_shift": None}
                continue
            # The best-scoring annotated neighbour decides the shift. Ties go to
            # the lower accession, which `near` is already sorted by.
            best = max(near, key=lambda pair: (pair[1], ))
            shift = round(precursor[accession] - precursor[best[0]],
                          self.rules.mass_shift_precision)
            same = abs(shift) <= self.rules.precursor_window
            out[accession] = {
                "claim": "same_compound" if same else "analogue_with_mass_shift",
                "block": blocks.pop(),
                "mass_shift": None if same else shift}
        unknown = sorted({row["claim"] for row in out.values()}
                         - set(self.rules.propagation_vocabulary))
        if unknown:
            raise NetworkError(
                f"the propagation rule produced {unknown}, which the declared "
                "vocabulary does not contain")
        return out
