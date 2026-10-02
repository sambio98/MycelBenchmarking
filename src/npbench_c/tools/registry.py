"""The pinned-tool registry: what is in the image, and under what controls.

Every entry declares three things that a version number alone does not capture:

  the CONTROLS the tool must always be run under, because several of these tools
  are deterministic only when told to be -- a thread count left to the host, a
  reordered output stream, or an implicit result cap will each move an answer;

  the INVOCATION the invariance suite exercises, on a pinned fixture, so the
  determinism claim is about a command this benchmark actually issues rather than
  about the tool in general;

  the NORMALISATION the suite is allowed to apply before comparing, with a reason
  per rule. This matters more than it looks. "Identical output" is cheap to claim
  if you quietly strip whatever differs, so the suite reports the raw verdict and
  the normalised verdict separately and names which rules fired. A tool that is
  byte-identical only after a line is removed is a tool whose output carries
  something unstable, and the benchmark should know which.

`python -m npbench_c.tools.registry --verify` checks the installed versions
against the pins and exits non-zero on drift. The image build runs it, so a pin
cannot drift silently into a score.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

#: The suite runs every tool at these thread counts. A tool with no thread flag
#: is run the same number of times at its only setting, so repeat-determinism is
#: still measured.
THREAD_COUNTS = (1, 8)
REPEATS = 2


@dataclass(frozen=True)
class Normalisation:
    """A declared licence to ignore one kind of line before comparing."""

    pattern: str
    reason: str


@dataclass(frozen=True)
class Canonicalisation:
    """A declared order-fixing transform applied to an artifact before USE.

    This is not the same thing as a normalisation and the difference is the point.
    A normalisation says "this line is not a result, ignore it when comparing". A
    canonicalisation says "the content is right but its order is not guaranteed,
    so fix the order before anybody reads it" -- and it therefore applies to the
    oracle and the agent alike, not only to the comparison inside this suite. A
    tool that needs one is admissible on a gold-producing path; a tool that needs
    one and does not get it is not.

    mode "sort_lines" is only sound for an artifact whose lines are independent
    records with no header, which the declaring tool must state.
    """

    mode: str
    reason: str


@dataclass(frozen=True)
class Invocation:
    """One command the suite runs, plus whatever has to exist first.

    Placeholders: {bin} the tool prefix, {fixtures} the fixture directory,
    {setup} a directory shared by every run of this invocation, {work} a
    directory private to one run, {threads} the thread count under test.

    `setup` runs ONCE, at one thread, and its output is shared. That is
    deliberate: it keeps the comparison about the step under test. Whether an
    index build is itself thread-invariant is a separate question, and a separate
    invocation is the place to ask it.
    """

    name: str
    argv: tuple[str, ...]
    artifact: str
    setup: tuple[tuple[str, ...], ...] = ()
    stdout_artifact: bool = False
    #: Set when the invocation exists to WITNESS a known non-determinism rather
    #: than to pass. The suite reports it as XFAIL and does not fail the run; an
    #: unexpected pass is reported as XPASS, because a tool that quietly became
    #: deterministic is a pin worth revisiting.
    expected_to_fail: str | None = None


@dataclass(frozen=True)
class Tool:
    name: str
    binary: str
    pinned_version: str
    needed_by: tuple[str, ...]
    version_argv: tuple[str, ...]
    version_pattern: str
    controls: tuple[str, ...]
    thread_flag: str | None
    invocations: tuple[Invocation, ...]
    normalisations: tuple[Normalisation, ...] = field(default_factory=tuple)
    canonicalisation: Canonicalisation | None = None

    @property
    def thread_counts(self) -> tuple[int, ...]:
        return THREAD_COUNTS if self.thread_flag else (1,)

    def installed_version(self, prefix: str | None = None) -> str | None:
        binary = f"{prefix}/{self.binary}" if prefix else self.binary
        if shutil.which(binary) is None:
            return None
        try:
            result = subprocess.run([binary, *self.version_argv],
                                    capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.TimeoutExpired):
            return None
        match = re.search(self.version_pattern,
                          (result.stdout or "") + "\n" + (result.stderr or ""))
        return match.group(1) if match else None


# --------------------------------------------------------------- the registry

HMMER = Tool(
    name="hmmer",
    binary="hmmsearch",
    pinned_version="3.4",
    needed_by=("T-L1-2", "T-L2-1", "T-L3-4"),
    version_argv=("-h",),
    version_pattern=r"HMMER (\S+)",
    controls=("pin version", "fixed --cpu"),
    thread_flag="--cpu",
    invocations=(
        Invocation(
            name="hmmsearch_tblout",
            setup=(("{bin}/hmmbuild", "--cpu", "1", "{setup}/halo.hmm",
                    "{fixtures}/halogenases.afa"),),
            argv=("{bin}/hmmsearch", "--cpu", "{threads}", "--tblout",
                  "{work}/hits.tbl", "-o", "/dev/null", "{setup}/halo.hmm",
                  "{fixtures}/halogenases.faa"),
            artifact="hits.tbl",
        ),
        Invocation(
            name="hmmbuild",
            argv=("{bin}/hmmbuild", "--cpu", "{threads}", "{work}/halo.hmm",
                  "{fixtures}/halogenases.afa"),
            artifact="halo.hmm",
        ),
    ),
    normalisations=(
        Normalisation(r"^# (Program|Version|Pipeline mode|Query file|Target file|"
                      r"Option settings|Current dir|Date|hmmsearch|hmmbuild|"
                      r"HMMER|Copyright|Freely distributed|- - -)",
                      "HMMER writes a provenance banner into --tblout and into "
                      "the HMM file: program name, the option string including "
                      "the thread count under test, the working directory and "
                      "the wall-clock date. All four differ between runs by "
                      "design and none is a result."),
        Normalisation(r"^DATE\s", "hmmbuild stamps the HMM with its build date."),
        Normalisation(r"^COM\s", "hmmbuild records its own command line, which "
                                 "contains the thread count under test."),
    ),
)

DIAMOND = Tool(
    name="diamond",
    binary="diamond",
    pinned_version="2.2.8",
    needed_by=("T-L1-3", "T-L2-3"),
    version_argv=("--version",),
    version_pattern=r"diamond version (\S+)",
    controls=("version >= 2.2.7", "fixed --threads",
              "explicit --max-target-seqs, which changes which hits survive",
              "NEVER --no-reorder: measured here, that flag is what BREAKS "
              "determinism rather than what provides it"),
    thread_flag="--threads",
    invocations=(
        Invocation(
            name="blastp_tabular",
            setup=(("{bin}/diamond", "makedb", "--in",
                    "{fixtures}/halogenases.faa", "-d", "{setup}/halo",
                    "--threads", "1", "--quiet"),),
            # --no-reorder is deliberately ABSENT. The suite measured it: with
            # the flag, 2.2.8 emits the same 192 hits in a different query order
            # on each multithreaded run; without it, every run at 1, 4 and 8
            # threads is byte-identical. DIAMOND's default restores query order,
            # and --no-reorder is documented as switching that off for speed --
            # so listing it as a determinism control, as this project's
            # provisioning notes originally did, had it exactly backwards.
            argv=("{bin}/diamond", "blastp", "-q", "{fixtures}/halogenases.faa",
                  "-d", "{setup}/halo", "-o", "{work}/hits.tsv",
                  "--threads", "{threads}",
                  "--max-target-seqs", "25", "--quiet"),
            artifact="hits.tsv",
        ),
        Invocation(
            name="blastp_tabular_no_reorder",
            setup=(("{bin}/diamond", "makedb", "--in",
                    "{fixtures}/halogenases.faa", "-d", "{setup}/halo",
                    "--threads", "1", "--quiet"),),
            # Kept as a STANDING witness, not as a working invocation. It is
            # expected to fail, and a future DIAMOND release that makes it pass
            # should be noticed rather than assumed.
            argv=("{bin}/diamond", "blastp", "-q", "{fixtures}/halogenases.faa",
                  "-d", "{setup}/halo", "-o", "{work}/hits.tsv",
                  "--threads", "{threads}", "--no-reorder",
                  "--max-target-seqs", "25", "--quiet"),
            artifact="hits.tsv",
            expected_to_fail="--no-reorder leaves query output order unspecified "
                             "under multithreading; measured, 2.2.8 gives three "
                             "distinct orderings over four runs with an identical "
                             "hit set",
        ),
    ),
)

MMSEQS2 = Tool(
    name="mmseqs2",
    binary="mmseqs",
    pinned_version="18.8cc5c",
    needed_by=("T-L2-3",),
    controls=("pin post-fix build", "always createindex", "fixed --threads",
              "explicit --max-seqs; upstream issue #277 reports results differing "
              "by core count, so this is the tool the suite exists for",
              "declared canonical line order: measured here, multithreaded runs "
              "emit the SAME hit set in a different order, and no flag changes "
              "that, so the order is fixed by the invocation contract instead"),
    version_argv=("version",),
    version_pattern=r"^(\S+)",
    thread_flag="--threads",
    invocations=(
        Invocation(
            name="search_convertalis",
            setup=(("{bin}/mmseqs", "createdb", "{fixtures}/halogenases.faa",
                    "{setup}/db", "-v", "0"),
                   ("{bin}/mmseqs", "createindex", "{setup}/db", "{setup}/tmp",
                    "--threads", "1", "-v", "0")),
            argv=("{bin}/mmseqs", "easy-search", "{fixtures}/halogenases.faa",
                  "{setup}/db", "{work}/hits.tsv", "{work}/tmp",
                  "--threads", "{threads}", "--max-seqs", "300", "-v", "0"),
            artifact="hits.tsv",
        ),
    ),
    canonicalisation=Canonicalisation(
        mode="sort_lines",
        reason="Tabular MMseqs2 output has no header and each line is an "
               "independent hit record, so sorting the lines changes nothing but "
               "the order. Measured: at one thread every run is byte-identical; "
               "at eight the line MULTISET is identical and only the order moves. "
               "Since no flag fixes the order, the canonical sort is part of the "
               "contract for every gold-producing use -- oracle and agent alike "
               "-- and not a licence this suite takes for itself.",
    ),
)

PRODIGAL = Tool(
    name="prodigal",
    binary="prodigal",
    pinned_version="2.6.3",
    needed_by=("T-L1-1", "T-L1-4"),
    version_argv=("-v",),
    version_pattern=r"Prodigal V(\S+?):",
    controls=("pin version AND mode: single and meta make different calls",
              "single-threaded, so the suite measures repeat determinism only",
              "Prodigal-GV is a different tool and is not in the image"),
    thread_flag=None,
    invocations=(
        Invocation(
            name="single_mode_proteins",
            argv=("{bin}/prodigal", "-i", "{fixtures}/synthetic.fna",
                  "-a", "{work}/proteins.faa", "-o", "{work}/genes.gbk",
                  "-p", "single", "-q"),
            artifact="proteins.faa",
        ),
    ),
)

MAFFT = Tool(
    name="mafft",
    binary="mafft",
    pinned_version="7.526",
    needed_by=("T-L2-3",),
    version_argv=("--version",),
    version_pattern=r"v(\S+)",
    controls=("pin major version", "prefer deterministic modes: fixed --retree "
              "and --maxiterate rather than the --auto strategy chooser, whose "
              "choice depends on the input size it measures"),
    thread_flag="--thread",
    invocations=(
        Invocation(
            name="retree2_noiterate",
            argv=("{bin}/mafft", "--thread", "{threads}", "--anysymbol",
                  "--retree", "2", "--maxiterate", "0", "--quiet",
                  "{fixtures}/halogenases.faa"),
            artifact="aln.afa",
            stdout_artifact=True,
        ),
    ),
)

BLAST = Tool(
    name="blast",
    binary="blastp",
    pinned_version="2.17.0",
    needed_by=("T-L1-3", "T-L2-2"),
    version_argv=("-version",),
    version_pattern=r"blastp: (\d+\.\d+\.\d+)",
    controls=("pin version", "fixed -num_threads", "explicit -max_target_seqs"),
    thread_flag="-num_threads",
    invocations=(
        Invocation(
            name="blastp_tabular",
            setup=(("{bin}/makeblastdb", "-in", "{fixtures}/halogenases.faa",
                    "-dbtype", "prot", "-out", "{setup}/halo"),),
            argv=("{bin}/blastp", "-query", "{fixtures}/halogenases.faa",
                  "-db", "{setup}/halo", "-outfmt", "6", "-num_threads",
                  "{threads}", "-max_target_seqs", "25", "-out",
                  "{work}/hits.tsv"),
            artifact="hits.tsv",
        ),
    ),
)

TOOLS: tuple[Tool, ...] = (BLAST, DIAMOND, HMMER, MAFFT, MMSEQS2, PRODIGAL)
BY_NAME = {tool.name: tool for tool in TOOLS}

#: Declared absences. Named here so "the image does not have it" is a recorded
#: decision rather than something an auditor has to infer from a missing row.
NOT_IN_IMAGE = {
    "antismash": "Needed by T-L3-1 and T-L3-5. Its pin is a PAIR -- binary plus "
                 "ClusterBlast and Pfam database releases -- that has to be "
                 "resolved together, and detection rules went 58 to 88 across "
                 "v5 to v7.1, so results must never be compared across versions. "
                 "Next image task.",
    "bigscape": "Needed by T-L3-2. Its pin depends on the antiSMASH and pyhmmer "
                "versions, so it follows antiSMASH.",
    "fasttreemp": "Banned outright: thread order affects the neighbour-joining "
                  "heuristic, so it can never sit on a gold-producing path.",
    "iqtree": "Not needed by any catalogued template. If it ever is: always "
              "--seed, always fixed -T N, never AUTO, and a single-gene ML tree "
              "is not admissible as gold.",
}


def verify(prefix: str | None = None) -> list[tuple[str, str, str | None, bool]]:
    rows = []
    for tool in TOOLS:
        found = tool.installed_version(prefix)
        rows.append((tool.name, tool.pinned_version, found,
                     found == tool.pinned_version))
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true",
                    help="compare installed versions against the pins")
    ap.add_argument("--prefix", default=None,
                    help="directory holding the binaries; default is PATH")
    args = ap.parse_args(argv)

    rows = verify(args.prefix)
    width = max(len(name) for name, *_ in rows)
    for name, pinned, found, ok in rows:
        mark = "  ok  " if ok else (" MISS " if found is None else " DRIFT")
        print(f"[{mark}] {name:<{width}}  pinned {pinned:<10} "
              f"installed {found or '-'}")
    bad = [name for name, _, _, ok in rows if not ok]
    print(f"\n{len(rows) - len(bad)}/{len(rows)} tools at their pinned version"
          + (f"; wrong or missing: {bad}" if bad else ""))
    for name, reason in sorted(NOT_IN_IMAGE.items()):
        print(f"[ n/a  ] {name:<{width}}  {reason.split('.')[0]}.")
    if not args.verify:
        return 0
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
