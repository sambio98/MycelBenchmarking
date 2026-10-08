"""Gene cluster family membership as a function of the clustering cutoff.

BiG-SCAPE partitions a set of BGCs into families, and the partition it returns is
a function of a distance cutoff the caller supplies. This engine runs the pinned
BiG-SCAPE once over a declared cutoff grid, reconstructs the partition over the
whole corpus at each cutoff, and reconciles it against families derived
independently from chemistry.

Five things are load-bearing, and each was settled by running the tool.

**The clustering table lists only clusters that joined a family.** At the lowest
declared cutoff it holds a third of the corpus; the rest are absent, not
unclustered-and-listed. So the partition over the corpus has to be RECONSTRUCTED:
listed records form their groups, and every record the table omits is a singleton.
A system that reads the table as the partition reports a corpus a third of its
real size and a family structure that happens to look cleaner.

**The join key is the GBK column, not the Record column.** `Record` carries a
region suffix (`BGC0000854.gbk_region_1`); `GBK` is the input stem. Verified at
every declared cutoff: exactly one row per GBK, every row a `region` record. A
corpus with multi-region inputs would break that, so the engine checks rather than
assumes.

**A family label is a name, not a result.** Two runs can agree completely on which
clusters belong together and disagree on what the groups are called. The graded
artifact is the PARTITION -- groups as sorted member sets, sorted among themselves
-- which is the projection the tool registry already declares.

**The cutoff is the answer's parameter, so it is declared, never defaulted.** A
default would make the partition depend on a number nobody published.

**Agreement with the chemistry-derived families is not monotone in the cutoff.**
Low cutoffs fragment the families; high cutoffs absorb a negative control. The
campaign measures where between the two the agreement peaks, which is the whole
point of a cutoff-sensitivity ladder.

Pure stdlib outside the tool call, sorted order throughout.
"""

from __future__ import annotations

import collections
import csv
import itertools
import json
import os
import pathlib
import re
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Mapping, Sequence

from npbench_c.tools import cache
from npbench_c.tools.resources import (
    DEFAULT_PREFIX,
    Resource,
    load_declared,
    tool_prefix,
)

TEMPLATE_ID = "gcf-cutoff"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"
RULES = "cutoff_rules.json"
PRECISION = 6

FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}

TIMEOUT_S = 7200

#: BiG-SCAPE writes one output directory per cutoff, named "<label>_<stamp>_c<cut>".
_CUTOFF_DIR = re.compile(r"_c([0-9.]+)$")


class CutoffError(RuntimeError):
    """The corpus or the environment does not support the clustering experiment."""


def _run(argv: Sequence[str]) -> subprocess.CompletedProcess:
    """Run a pinned tool with the prefix appended to PATH, last."""
    prefix = str(tool_prefix())
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if prefix in parts:
        parts.remove(prefix)
    env = {**os.environ, **FIXED_ENV, "PATH": os.pathsep.join([*parts, prefix])}
    return subprocess.run(list(argv), capture_output=True, text=True, env=env,
                          timeout=TIMEOUT_S)


# ---------------------------------------------------------------------- rules


@dataclass(frozen=True)
class Rules:
    cutoffs: tuple[float, ...]
    baseline_cutoff: float
    cores: int
    include_gbk: str
    verdict_order: tuple[str, ...]
    resources: Mapping[str, Resource]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
        clustering = spec["clustering"]
        resources = load_declared(spec["external_resources"])
        cutoffs = tuple(float(c) for c in clustering["cutoffs"])
        baseline = float(clustering["baseline_cutoff"])
        if baseline not in cutoffs:
            raise CutoffError(
                f"the baseline cutoff {baseline} is not in the declared grid "
                f"{cutoffs}, so R2's partition is not one the sweep computed")
        return cls(cutoffs=cutoffs, baseline_cutoff=baseline,
                   cores=int(clustering["cores"]),
                   include_gbk=clustering["include_gbk"],
                   verdict_order=tuple(clustering["verdict_order"]),
                   resources=resources)

    def label(self, cutoff: float) -> str:
        """The declared string form of a cutoff, used as a report key."""
        return f"{cutoff:g}"


# --------------------------------------------------------------------- corpus


@dataclass(frozen=True)
class Clustering:
    corpus: tuple[str, ...]
    designed: Mapping[str, tuple[str, ...]]
    controls: tuple[str, ...]
    input_dir: pathlib.Path
    rules: Rules

    @classmethod
    def load(cls, campaign: pathlib.Path, gbk_dir: str,
             families_file: str) -> "Clustering":
        campaign = pathlib.Path(campaign)
        rules = Rules.load(campaign)
        directory = campaign / gbk_dir
        if not directory.is_dir():
            raise CutoffError(f"no input directory at {directory}")
        corpus = tuple(sorted(p.stem for p in directory.glob("*.gbk")))
        if not corpus:
            raise CutoffError(f"{directory} holds no .gbk files")
        payload = json.loads((campaign / families_file).read_text())
        designed = {name: tuple(sorted(members))
                    for name, members in sorted(payload["families"].items())}
        controls = tuple(sorted(payload["controls"]))
        members = sorted({m for v in designed.values() for m in v} | set(controls))
        if members != list(corpus):
            missing = sorted(set(corpus) - set(members))
            extra = sorted(set(members) - set(corpus))
            raise CutoffError(
                "the declared families plus controls must partition the shipped "
                f"corpus exactly; missing {missing[:4]}, extra {extra[:4]}")
        overlap = [name for name, v in designed.items() if set(v) & set(controls)]
        if overlap:
            raise CutoffError(
                f"{overlap} share a member with the negative controls, so the "
                "controls are not controls")
        return cls(corpus=corpus, designed=designed, controls=controls,
                   input_dir=directory, rules=rules)

    # ------------------------------------------------------------ tool calls

    def bigscape_version(self) -> str:
        done = _run([str(tool_prefix() / "bigscape"), "--version"])
        if done.returncode != 0:
            raise CutoffError(f"bigscape --version failed: {done.stderr[-400:]}")
        return done.stdout.strip()

    def _cluster(self) -> dict[str, dict[str, list[str]]]:
        """One BiG-SCAPE run over the whole declared cutoff grid.

        Every cutoff comes from the same run, so the distances behind them are
        identical by construction and a difference between two cutoffs is the
        cutoff rather than two searches.
        """
        pfam = self.rules.resources["pfam_a_hmm"].resolve()
        grid = ",".join(self.rules.label(c) for c in self.rules.cutoffs)
        with tempfile.TemporaryDirectory(prefix="npbench-gcf-") as tmp:
            work = pathlib.Path(tmp)
            out = work / "out"
            done = _run([str(tool_prefix() / "bigscape"), "cluster",
                         "-i", str(self.input_dir),
                         "-o", str(out),
                         "-p", str(pfam),
                         "--cores", str(self.rules.cores),
                         "--gcf-cutoffs", grid,
                         "--include-gbk", self.rules.include_gbk,
                         "--label", "campaign"])
            if done.returncode != 0:
                tail = (done.stdout or "")[-600:] + (done.stderr or "")[-600:]
                raise CutoffError(f"bigscape cluster failed: {tail}")
            found: dict[str, dict[str, list[str]]] = {}
            for directory in sorted((out / "output_files").glob("*_c*")):
                match = _CUTOFF_DIR.search(directory.name)
                if not match:
                    continue
                label = self.rules.label(float(match.group(1)))
                rows: dict[str, list[str]] = {}
                tables = sorted(directory.glob("*/*_clustering_c*.tsv"))
                if not tables:
                    raise CutoffError(
                        f"no clustering table under {directory}; the run produced "
                        "a cutoff directory with nothing in it")
                for table in tables:
                    with open(table) as handle:
                        for row in csv.DictReader(handle, delimiter="\t"):
                            # The GBK column is the input stem; Record carries a
                            # region suffix. One row per GBK is checked rather
                            # than assumed: a multi-region input would need a
                            # different reconstruction rule than this one.
                            stem = row["GBK"]
                            if stem in rows:
                                raise CutoffError(
                                    f"{stem} appears twice at cutoff {label}; the "
                                    "corpus has a multi-region input and the "
                                    "declared reconstruction rule does not cover "
                                    "it")
                            if row["Record_Type"] != "region":
                                raise CutoffError(
                                    f"{stem} at cutoff {label} is a "
                                    f"{row['Record_Type']!r} record, not a region")
                            rows[stem] = [table.parent.name, row["Family"]]
                found[label] = rows
            missing = {self.rules.label(c) for c in self.rules.cutoffs} - set(found)
            if missing:
                raise CutoffError(
                    f"the run produced no output for cutoffs {sorted(missing)}")
            return found

    def assignments(self) -> dict[str, dict[str, tuple[str, str]]]:
        """Per cutoff, the class and family of every LISTED record.

        Memoised only where a tool memo is enabled, which is stub sweeps alone.
        """
        def compute() -> dict[str, dict[str, list[str]]]:
            return self._cluster()

        key_parts = [
            "gcf_cutoff.assignments",
            self._corpus_digest(),
            self.bigscape_version(),
            self.rules.resources["pfam_a_hmm"].sha256,
            ",".join(self.rules.label(c) for c in self.rules.cutoffs),
            self.rules.include_gbk,
            str(self.rules.cores),
        ]
        payload = (cache.memoize_json(key_parts, compute) if cache.enabled()
                   else compute())
        return {cut: {stem: (value[0], value[1]) for stem, value in rows.items()}
                for cut, rows in payload.items()}

    def _corpus_digest(self) -> str:
        """One digest over every input file, in sorted name order."""
        import hashlib

        digest = hashlib.sha256()
        for path in sorted(self.input_dir.glob("*.gbk")):
            digest.update(path.name.encode())
            digest.update(cache.file_digest(path).encode())
        return digest.hexdigest()

    # ------------------------------------------------------- the partition

    def partition(self, assignment: Mapping[str, tuple[str, str]]) -> list[list[str]]:
        """The partition over the WHOLE corpus at one cutoff.

        Listed records form their groups; every record the table omits is a
        singleton. That reconstruction is the campaign's first trap: the table
        lists only clusters that joined a family, so reading it as the partition
        silently shrinks the corpus.
        """
        groups: dict[tuple[str, str], list[str]] = {}
        for stem, key in sorted(assignment.items()):
            groups.setdefault(key, []).append(stem)
        unlisted = [[stem] for stem in self.corpus if stem not in assignment]
        out = sorted([sorted(members) for members in groups.values()] + unlisted)
        flat = sorted(m for group in out for m in group)
        if flat != list(self.corpus):
            raise CutoffError(
                "the reconstructed partition does not cover the corpus exactly; "
                "the clustering table names a record the corpus does not hold")
        return out

    def group_census(self, partition: Sequence[Sequence[str]]) -> dict:
        sizes = collections.Counter(len(group) for group in partition)
        return {
            "groups": len(partition),
            "size_census": {str(size): n for size, n in sorted(sizes.items())},
            "largest_group": max((len(g) for g in partition), default=0),
            "singletons": sizes.get(1, 0),
            "clustered_records": sum(len(g) for g in partition if len(g) > 1),
        }

    def class_census(self, assignment: Mapping[str, tuple[str, str]]) -> dict[str, int]:
        counts = collections.Counter(cls for cls, _ in assignment.values())
        return dict(sorted(counts.items()))

    # --------------------------------------------------- the reconciliation

    def family_verdict(self, members: Sequence[str],
                       partition: Sequence[Sequence[str]]) -> dict:
        """How one chemistry-derived family fares against the called partition."""
        wanted = set(members)
        spanning = [sorted(g) for g in partition if wanted & set(g)]
        outsiders = sorted({m for g in spanning for m in g} - wanted)
        exact = len(spanning) == 1 and not outsiders
        if exact:
            verdict = "exact"
        elif len(spanning) > 1 and outsiders:
            verdict = "fragmented_and_absorbed"
        elif len(spanning) > 1:
            verdict = ("dissolved" if all(len(g) == 1 for g in spanning)
                       else "fragmented")
        else:
            verdict = "absorbed"
        return {
            "verdict": verdict,
            "groups_spanned": len(spanning),
            "group_sizes": sorted((len(g) for g in spanning), reverse=True),
            "outsiders": outsiders,
        }

    def reconciliation(self, partition: Sequence[Sequence[str]]) -> dict:
        per_family = {name: self.family_verdict(members, partition)
                      for name, members in sorted(self.designed.items())}
        counts = collections.Counter(v["verdict"] for v in per_family.values())
        absorbed = sorted(
            control for control in self.controls
            if any(len(g) > 1 and control in g for g in partition))
        return {
            "per_family": per_family,
            "verdict_census": {v: counts.get(v, 0)
                               for v in self.rules.verdict_order},
            "families_exact": sum(1 for v in per_family.values()
                                  if v["verdict"] == "exact"),
            "controls_absorbed": absorbed,
        }

    # ------------------------------------------------------ pair statistics

    def designed_pairs(self) -> set[tuple[str, str]]:
        pairs: set[tuple[str, str]] = set()
        for members in self.designed.values():
            pairs |= {tuple(sorted(pair))
                      for pair in itertools.combinations(sorted(members), 2)}
        return pairs

    def pair_statistics(self, partition: Sequence[Sequence[str]]) -> dict:
        """Co-clustering agreement at the level of pairs.

        A pair-counting measure rather than a label comparison, because family
        labels are names: two partitions can agree perfectly and use different
        ones. Recovered pairs are designed pairs the tool co-clusters; false joins
        are pairs it co-clusters that no chemistry group puts together.
        """
        designed = self.designed_pairs()
        called: set[tuple[str, str]] = set()
        for group in partition:
            called |= {tuple(sorted(pair))
                       for pair in itertools.combinations(sorted(group), 2)}
        recovered = len(called & designed)
        false_joins = len(called - designed)
        missed = len(designed - called)
        union = recovered + false_joins + missed
        return {
            "designed_pairs": len(designed),
            "co_clustered_pairs": len(called),
            "recovered_pairs": recovered,
            "missed_pairs": missed,
            "false_join_pairs": false_joins,
            "pair_jaccard": round(recovered / union, PRECISION) if union else 0.0,
        }

    # ------------------------------------------------------- the cutoff sweep

    def sweep(self, assignments: Mapping[str, Mapping[str, tuple[str, str]]]) -> dict:
        """Every declared cutoff: the partition census, the pairs, the verdicts."""
        rows = {}
        for cutoff in self.rules.cutoffs:
            label = self.rules.label(cutoff)
            partition = self.partition(assignments[label])
            reconciliation = self.reconciliation(partition)
            rows[label] = {
                **self.group_census(partition),
                **self.pair_statistics(partition),
                "families_exact": reconciliation["families_exact"],
                "controls_absorbed": reconciliation["controls_absorbed"],
                "family_verdicts": {name: block["verdict"] for name, block
                                    in sorted(reconciliation["per_family"].items())},
            }
        return rows

    def sweep_summary(self, sweep: Mapping[str, Mapping]) -> dict:
        """The three boundary claims the sweep supports.

        The best cutoff under the declared tie-break (highest pair Jaccard, then
        the lowest cutoff), the first cutoff at which a negative control is
        absorbed, and per family the first cutoff at which it is recovered
        exactly -- or null, because a family need not be recoverable at any
        declared cutoff and reporting a number there would be an invention.
        """
        ordered = sorted(sweep, key=float)
        best = max(ordered, key=lambda label: (sweep[label]["pair_jaccard"],
                                               -float(label)))
        first_false = next((label for label in ordered
                            if sweep[label]["false_join_pairs"]), None)
        first_absorbed = next((label for label in ordered
                               if sweep[label]["controls_absorbed"]), None)
        first_exact = {}
        for name in sorted(self.designed):
            first_exact[name] = next(
                (label for label in ordered
                 if sweep[label]["family_verdicts"][name] == "exact"), None)
        return {
            "best_cutoff": best,
            "best_pair_jaccard": sweep[best]["pair_jaccard"],
            "first_cutoff_with_a_false_join": first_false,
            "first_cutoff_absorbing_a_control": first_absorbed,
            "first_cutoff_recovering_each_family": first_exact,
        }
