"""The thread-invariance suite.

Run every pinned tool on a pinned fixture at 1 and 8 threads, twice each, and
require the gold-producing output to be identical. Anything that fails and cannot
be pinned into determinism is made single-threaded in the image or excluded from
gold-producing paths -- that is the rule, and this is what enforces it.

Two verdicts are reported per invocation, never one:

  RAW         every run byte-identical with no normalisation at all
  NORMALISED  identical after the tool's declared normalisation rules

Reporting only the second would make the suite easy to pass by stripping whatever
differs. Reporting only the first would fail every tool that stamps its own
command line into its output, which several of these do and which says nothing
about determinism. So both are reported, and the rules that actually fired are
named, because a line that has to be removed is a line carrying something
unstable and the benchmark should know which.

A tool's absence is reported as SKIP with its reason, never as a pass. A suite
that quietly measures nothing is worse than one that fails.

Usage:
  python -m npbench_c.tools.invariance --all [--json out.json] [--prefix DIR]
  python -m npbench_c.tools.invariance --tool diamond
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass

from npbench_c.tools.registry import (
    BANNED_BINARIES,
    BY_NAME,
    NOT_IN_IMAGE,
    REPEATS,
    TOOLS,
    Invocation,
    Tool,
    banned_present,
    fingerprint,
)

FIXTURES = pathlib.Path(__file__).resolve().parent / "fixtures"
TIMEOUT_S = 1800

PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"
#: An invocation that exists to witness a known non-determinism. XFAIL does not
#: fail the suite; XPASS does not either, but it is reported loudly, because a
#: tool that quietly became deterministic is a pin worth revisiting.
XFAIL, XPASS = "XFAIL", "XPASS"

#: The environment every run gets. A tool that reads its thread count from the
#: environment instead of its flag would otherwise make the comparison
#: meaningless -- both arms would silently use the host's core count.
FIXED_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0", "LC_ALL": "C", "TZ": "UTC",
}


@dataclass
class RunResult:
    threads: int
    repeat: int
    returncode: int
    raw_digest: str | None
    normalised_digest: str | None
    canonical_digest: str | None
    rules_fired: tuple[str, ...]
    stderr_tail: str


def _normalise(text: str, tool: Tool) -> tuple[str, tuple[str, ...]]:
    """Drop lines matching the tool's declared rules; report which rules fired."""
    if not tool.normalisations:
        return text, ()
    patterns = [(re.compile(rule.pattern), rule.pattern)
                for rule in tool.normalisations]
    kept, fired = [], set()
    for line in text.splitlines():
        dropped = False
        for compiled, label in patterns:
            if compiled.search(line):
                fired.add(label)
                dropped = True
                break
        if not dropped:
            kept.append(line)
    return "\n".join(kept) + "\n", tuple(sorted(fired))


class ProjectionError(RuntimeError):
    """A projection found nothing to project.

    This is a hard failure rather than an empty result, because a projection that
    silently yields nothing turns the determinism test into a tautology: every run
    would agree on the same emptiness. If the tool's output structure changed, the
    suite must say so.
    """


def _project_antismash_regions(raw: bytes) -> str:
    """Per record, each detected region's coordinates and sorted products.

    This is the declared answer surface for a BGC-detection campaign. Everything
    else antiSMASH writes -- the input path, its own version, the record
    timestamp, the full HMM hit table, the HTML -- is by this declaration not a
    result. Region order is sorted, so a reordered features list is not mistaken
    for a different detection.
    """
    document = json.loads(raw)
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise ProjectionError("antiSMASH JSON carries no records list")

    def sort_key(location: str) -> tuple:
        """Numeric where the location parses, lexical where it does not.

        A string sort alone is deterministic, which is all this suite needs, but
        the projection is also what a campaign would grade -- and there
        "[100:2000]" before "[10:90]" is an order no reader would expect. So the
        coordinates are parsed when they can be, with the raw string as the
        tie-break and the fallback.
        """
        digits = re.findall(r"-?\d+", location)
        if len(digits) >= 2:
            return (0, int(digits[0]), int(digits[1]), location)
        return (1, 0, 0, location)

    lines: list[str] = []
    for record in records:
        name = str(record.get("id", "?"))
        regions = []
        for feature in (record.get("features") or []):
            if feature.get("type") != "region":
                continue
            location = str(feature.get("location", ""))
            products = sorted(
                str(x) for x in ((feature.get("qualifiers") or {}).get("product") or []))
            regions.append((location, tuple(products)))
        for location, products in sorted(regions, key=lambda r: sort_key(r[0])):
            lines.append(f"{name}\t{location}\t{','.join(products)}")
    if not lines:
        raise ProjectionError(
            "no region features found; either the fixture has no detectable "
            "cluster or the JSON structure moved, and both must be looked at "
            "rather than passed as 'identical'")
    return "\n".join(lines) + "\n"


def _project_antismash_domains(raw: bytes) -> str:
    """Per record and CDS, the ordered NRPS/PKS domain architecture.

    This is the claim a domain-architecture or substrate-specificity campaign
    would grade: which domains, in which order, over which residues. The e-value
    and bitscore that accompany each hit are deliberately NOT projected -- a
    float's last digits are a reduction-order artefact, and nothing downstream
    needs them to state an architecture.
    """
    document = json.loads(raw)
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise ProjectionError("antiSMASH JSON carries no records list")

    lines: list[str] = []
    for record in records:
        name = str(record.get("id", "?"))
        module = (record.get("modules") or {}).get(
            "antismash.detection.nrps_pks_domains")
        if not module:
            continue
        for cds, result in sorted((module.get("cds_results") or {}).items()):
            hits = sorted(
                (int(h.get("query_start", -1)), int(h.get("query_end", -1)),
                 str(h.get("hit_id", "?")))
                for h in (result.get("domain_hmms") or []))
            if not hits:
                continue
            architecture = ";".join(f"{hid}:{start}-{end}"
                                    for start, end, hid in hits)
            lines.append(f"{name}\t{cds}\t{architecture}")
    if not lines:
        raise ProjectionError(
            "no NRPS/PKS domain results found; either the run did not include the "
            "domain detection module or the JSON structure moved, and both must be "
            "looked at rather than passed as 'identical'")
    return "\n".join(lines) + "\n"


def _project_bigscape_partition(raw: bytes) -> str:
    """The gene-cluster-family PARTITION, with the family labels thrown away.

    `FAM_00001` is a name, not a result. Two runs can agree completely on which
    clusters belong together and disagree on what the groups are called, and
    comparing labels would score that as a difference. So each group becomes its
    sorted member set, the groups are sorted among themselves, and the connected
    component number goes the same way as the label.

    The concatenated input is the per-class clustering tables, each preceded by a
    `# <basename>` line; the basename carries the class, which IS a result.
    """
    text = raw.decode(errors="replace")
    groups: dict[tuple[str, str], set[str]] = {}
    current = "?"
    for line in text.splitlines():
        if line.startswith("# "):
            basename = line[2:].strip()
            current = basename.split("_clustering_")[0]
            continue
        if not line.strip() or line.startswith("Record\t"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        groups.setdefault((current, fields[5]), set()).add(fields[0])
    if not groups:
        raise ProjectionError(
            "no clustering rows found; either the run produced no families or the "
            "table layout moved, and both must be looked at rather than passed as "
            "'identical'")
    rendered = sorted((cls, ",".join(sorted(members)))
                      for (cls, _label), members in groups.items())
    return "\n".join(f"{cls}\t{members}" for cls, members in rendered) + "\n"


PROJECTIONS = {
    "antismash_regions": _project_antismash_regions,
    "antismash_nrps_pks_domains": _project_antismash_domains,
    "bigscape_gcf_partition": _project_bigscape_partition,
}


def _pfam_path(prefix: str) -> str | None:
    """The Pfam HMM library, reused from the antiSMASH databases.

    BiG-SCAPE needs Pfam-A.hmm and antiSMASH already ships one with its pressed
    indexes, so the image carries one copy rather than two -- 1.5 GB saved, and
    one fewer thing whose release has to be pinned separately.
    """
    import glob

    matches = sorted(glob.glob(str(
        pathlib.Path(prefix).parent / "lib" / "python3.*" / "site-packages"
        / "antismash" / "databases" / "pfam" / "*" / "Pfam-A.hmm")))
    return matches[-1] if matches else None


def _canonicalise(text: str, tool: Tool) -> str:
    """Apply the tool's declared order-fixing transform, if it has one.

    Unlike a normalisation, this is part of the invocation contract for every
    gold-producing use of the tool, not a licence this suite takes for itself.
    """
    if tool.canonicalisation is None:
        return text
    if tool.canonicalisation.mode == "sort_lines":
        return "\n".join(sorted(text.splitlines())) + "\n"
    raise ValueError(
        f"{tool.name}: unknown canonicalisation mode "
        f"{tool.canonicalisation.mode!r}")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()[:16]


def _format(argv, **subs) -> list[str]:
    return [token.format(**subs) for token in argv]


def run_invocation(tool: Tool, invocation: Invocation, prefix: str,
                   root: pathlib.Path) -> dict:
    """Run one invocation at every declared thread count and compare."""
    setup_dir = root / "setup"
    setup_dir.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, **FIXED_ENV}
    subs = {"bin": prefix, "fixtures": str(FIXTURES), "setup": str(setup_dir),
            "pfam": _pfam_path(prefix) or ""}
    if "{pfam}" in " ".join(invocation.argv) and not subs["pfam"]:
        return {"invocation": invocation.name, "status": FAIL,
                "detail": "Pfam-A.hmm not found in the antiSMASH databases; this "
                          "invocation reuses it rather than downloading a second "
                          "copy, so the databases must be present"}

    for command in invocation.setup:
        argv = _format(command, **subs, work=str(setup_dir), threads="1")
        done = subprocess.run(argv, capture_output=True, text=True, env=env,
                              timeout=TIMEOUT_S)
        if done.returncode != 0:
            return {"invocation": invocation.name, "status": FAIL,
                    "detail": f"setup failed: {' '.join(argv)}\n"
                              f"{done.stderr[-800:]}"}

    runs: list[RunResult] = []
    for threads in tool.thread_counts:
        for repeat in range(REPEATS):
            work = root / f"t{threads}_r{repeat}"
            work.mkdir(parents=True, exist_ok=True)
            argv = _format(invocation.argv, **subs, work=str(work),
                           threads=str(threads))
            timeout = invocation.timeout_s
            if invocation.stdout_artifact:
                target = work / invocation.artifact
                with target.open("wb") as handle:
                    done = subprocess.run(argv, stdout=handle,
                                          stderr=subprocess.PIPE, env=env,
                                          timeout=timeout)
                stderr = (done.stderr or b"").decode(errors="replace")
            else:
                done = subprocess.run(argv, capture_output=True, text=True,
                                      env=env, timeout=timeout)
                stderr = done.stderr or ""
                target = work / invocation.artifact
            if "*" in invocation.artifact:
                # A glob: the output directory is timestamped, so the artifact is
                # every matching table, concatenated in basename order with a
                # header line naming each. The basename is stable; the directory
                # above it is not.
                import glob as _glob

                matches = sorted(_glob.glob(str(work / invocation.artifact)),
                                 key=lambda m: pathlib.Path(m).name)
                if done.returncode != 0 or not matches:
                    runs.append(RunResult(threads, repeat, done.returncode, None,
                                          None, None, (), stderr[-800:]))
                    continue
                raw = b"".join(
                    f"# {pathlib.Path(m).name}\n".encode()
                    + pathlib.Path(m).read_bytes() for m in matches)
            else:
                if done.returncode != 0 or not target.is_file():
                    runs.append(RunResult(threads, repeat, done.returncode, None,
                                          None, None, (), stderr[-800:]))
                    continue
                raw = target.read_bytes()
            if invocation.projection is not None:
                try:
                    raw = PROJECTIONS[invocation.projection.mode](raw).encode()
                except (ProjectionError, json.JSONDecodeError, KeyError) as exc:
                    runs.append(RunResult(threads, repeat, 0, None, None, None, (),
                                          f"projection failed: {exc}"))
                    continue
            normalised, fired = _normalise(raw.decode(errors="replace"), tool)
            canonical = _canonicalise(normalised, tool)
            runs.append(RunResult(threads, repeat, 0, _digest(raw),
                                  _digest(normalised.encode()),
                                  _digest(canonical.encode()), fired,
                                  stderr[-300:]))

    failed = [r for r in runs if r.raw_digest is None]
    if failed:
        return {"invocation": invocation.name, "status": FAIL,
                "detail": f"{len(failed)}/{len(runs)} runs produced no artifact; "
                          f"first stderr: {failed[0].stderr_tail}"}

    raw_digests = {r.raw_digest for r in runs}
    norm_digests = {r.normalised_digest for r in runs}
    canon_digests = {r.canonical_digest for r in runs}
    fired = sorted({rule for r in runs for rule in r.rules_fired})
    # Repeat determinism at a single thread count is a different property from
    # invariance across thread counts, and a tool can have one without the other;
    # both DIAMOND under --no-reorder and MMseqs2 do.
    per_threads = {t: {r.canonical_digest for r in runs if r.threads == t}
                   for t in tool.thread_counts}
    repeatable = all(len(v) == 1 for v in per_threads.values())

    invariant = len(canon_digests) == 1
    if invocation.expected_to_fail:
        status = XPASS if invariant else XFAIL
    else:
        status = PASS if invariant else FAIL

    # A raw match while a normalisation rule fired is a COINCIDENCE, not a
    # property: hmmbuild stamps a build date, and four runs inside one second
    # agree on it. Saying "raw-identical" there would turn a timing accident into
    # a claim, so it is reported as incidental and the normalised verdict is the
    # one that means anything.
    incidental = len(raw_digests) == 1 and bool(fired)

    # What was compared, said plainly. For a projected invocation the raw
    # artifact is certainly NOT identical -- antiSMASH's JSON carries the input
    # path and a timestamp -- so reporting "raw-identical" there would overstate
    # the result in exactly the way the incidental rule below guards against.
    compared = (f"projection:{invocation.projection.mode}"
                if invocation.projection else "artifact")

    if len(raw_digests) == 1 and not fired:
        detail = ("identical across all runs with no normalisation"
                  if not invocation.projection else
                  f"the {invocation.projection.mode} projection is identical "
                  "across all runs; the raw artifact is not, and by this "
                  "invocation's declaration is not a result")
    elif incidental:
        detail = (f"identical after {len(fired)} declared normalisation rule(s); "
                  "the raw bytes also agreed, but a declared-unstable field is "
                  "present and only happened to match")
    elif len(norm_digests) == 1:
        detail = f"identical after {len(fired)} declared normalisation rule(s)"
    elif invariant:
        detail = ("identical only after the declared canonicalisation "
                  f"({tool.canonicalisation.mode}); the hit set is the same in "
                  "every run and only its order moved")
    else:
        detail = (f"NOT invariant: {len(canon_digests)} distinct outputs across "
                  f"{len(runs)} runs even after canonicalisation; repeatable at "
                  f"fixed threads={repeatable}")
    if invocation.expected_to_fail and not invariant:
        detail = f"expected: {invocation.expected_to_fail}"
    elif invocation.expected_to_fail and invariant:
        detail = ("XPASS -- this invocation is a witness for a KNOWN "
                  "non-determinism and it just passed. Re-read the pin before "
                  "trusting it: " + invocation.expected_to_fail)

    return {
        "invocation": invocation.name,
        "status": status,
        "thread_counts": list(tool.thread_counts),
        "runs": len(runs),
        "compared": compared,
        "raw_identical": (len(raw_digests) == 1 and not fired
                          and not invocation.projection),
        "projection_identical": (len(raw_digests) == 1 and not fired
                                 and bool(invocation.projection)),
        "raw_identical_incidental": incidental,
        "normalised_identical": len(norm_digests) == 1,
        "canonical_identical": invariant,
        "repeatable_at_fixed_threads": repeatable,
        "normalisation_rules_fired": fired,
        "canonicalisation": (tool.canonicalisation.mode
                             if tool.canonicalisation else None),
        "projection": (invocation.projection.mode
                       if invocation.projection else None),
        "digests": {f"t{r.threads}_r{r.repeat}":
                    {"raw": r.raw_digest, "normalised": r.normalised_digest,
                     "canonical": r.canonical_digest}
                    for r in runs},
        "detail": detail,
    }


def run_tool(tool: Tool, prefix: str | None = None,
             root: pathlib.Path | None = None) -> dict:
    binary = f"{prefix}/{tool.binary}" if prefix else tool.binary
    if shutil.which(binary) is None:
        return {"tool": tool.name, "status": SKIP,
                "detail": f"{tool.binary} is not on PATH; the pinned-tool image "
                          "is not present, so nothing was measured"}
    found = tool.installed_version(prefix)
    if found != tool.pinned_version:
        return {"tool": tool.name, "status": SKIP,
                "detail": f"installed {found}, pinned {tool.pinned_version}; "
                          "an invariance result for the wrong version would not "
                          "be about the image"}

    resolved_prefix = prefix or str(pathlib.Path(shutil.which(binary)).parent)
    # A tool whose behaviour lives in data files gets that recorded too, so a
    # campaign can pin what it was built against. antiSMASH's package root is a
    # sibling of its bin directory inside the environment.
    package_root = None
    if tool.fingerprint is not None:
        import glob

        candidates = sorted(glob.glob(str(
            pathlib.Path(resolved_prefix).parent / "lib" / "python3.*"
            / "site-packages" / tool.name)))
        package_root = candidates[-1] if candidates else None
    own = root is None
    base = pathlib.Path(tempfile.mkdtemp(prefix=f"invariance_{tool.name}_")) \
        if own else root / tool.name
    base.mkdir(parents=True, exist_ok=True)
    try:
        results = [run_invocation(tool, inv, resolved_prefix, base / inv.name)
                   for inv in tool.invocations]
    finally:
        if own:
            shutil.rmtree(base, ignore_errors=True)

    failing = [r for r in results if r["status"] == FAIL]
    return {
        "tool": tool.name,
        "status": PASS if not failing else FAIL,
        "version": found,
        "needed_by": list(tool.needed_by),
        "controls": list(tool.controls),
        "slow": tool.slow,
        "fingerprint": fingerprint(tool, package_root),
        "fingerprint_paths": (list(tool.fingerprint.paths)
                              if tool.fingerprint else None),
        "invocations": results,
        "canonicalisation": (tool.canonicalisation.mode
                             if tool.canonicalisation else None),
        "detail": (f"{len([r for r in results if r['status'] in (PASS, XFAIL)])}"
                   f"/{len(results)} invocation(s) invariant at "
                   f"{tool.thread_counts}"
                   + ("" if tool.thread_flag else
                      " (no thread flag; repeat determinism only)")
                   if not failing else
                   f"{len(failing)}/{len(results)} invocation(s) not invariant: "
                   f"{[r['invocation'] for r in failing]}"),
    }


def run_all(prefix: str | None = None) -> dict:
    report: dict = {
        "thread_counts": list(__import__(
            "npbench_c.tools.registry", fromlist=["THREAD_COUNTS"]).THREAD_COUNTS),
        "repeats": REPEATS,
        "fixtures": sorted(p.name for p in FIXTURES.iterdir() if p.is_file()),
        "not_in_image": NOT_IN_IMAGE,
        "banned_binaries": {
            "declared": BANNED_BINARIES,
            "present": banned_present(prefix),
        },
        "tools": {},
    }
    for tool in TOOLS:
        report["tools"][tool.name] = run_tool(tool, prefix)
    statuses = [r["status"] for r in report["tools"].values()]
    invocations = [i for r in report["tools"].values()
                   for i in r.get("invocations", [])]
    report["summary"] = {
        "pass": statuses.count(PASS), "fail": statuses.count(FAIL),
        "skip": statuses.count(SKIP),
        "measured": statuses.count(PASS) + statuses.count(FAIL),
        "invocations": len(invocations),
        "xfail": sum(1 for i in invocations if i["status"] == XFAIL),
        "xpass": sum(1 for i in invocations if i["status"] == XPASS),
        "needed_canonicalisation": sorted(
            r["tool"] for r in report["tools"].values()
            if r.get("canonicalisation")),
    }
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--tool", action="append", default=[])
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)

    if args.tool:
        report = {"tools": {name: run_tool(BY_NAME[name], args.prefix)
                            for name in args.tool}}
        statuses = [r["status"] for r in report["tools"].values()]
        report["summary"] = {"pass": statuses.count(PASS),
                             "fail": statuses.count(FAIL),
                             "skip": statuses.count(SKIP),
                             "measured": statuses.count(PASS) + statuses.count(FAIL)}
    elif args.all:
        report = run_all(args.prefix)
    else:
        ap.error("pass --all or --tool NAME")

    marks = {PASS: "  ok  ", FAIL: " FAIL ", SKIP: " skip ",
             XFAIL: " xfail", XPASS: " XPASS"}
    for name, result in report["tools"].items():
        print(f"[{marks[result['status']]}] {name:<9} {result['detail']}")
        for inv in result.get("invocations", []):
            if inv.get("raw_identical"):
                flag = "raw-identical"
            elif inv.get("projection_identical"):
                flag = f"projection-identical ({inv['projection']})"
            elif inv.get("raw_identical_incidental"):
                flag = (f"normalised ({len(inv['normalisation_rules_fired'])} "
                        "rules; raw match incidental)")
            elif inv.get("normalised_identical"):
                flag = f"normalised ({len(inv['normalisation_rules_fired'])} rules)"
            elif inv.get("canonical_identical"):
                flag = f"canonical ({inv['canonicalisation']})"
            else:
                flag = ""
            print(f"           - {inv['invocation']:<26} "
                  f"{marks.get(inv['status'], inv['status']).strip():<5} "
                  f"{flag or inv['detail'][:70]}")
            if inv["status"] in (FAIL, XFAIL, XPASS):
                print(f"             {inv['detail']}")

    summary = report["summary"]
    print(f"\n{summary['pass']} pass / {summary['fail']} fail / "
          f"{summary['skip']} skip ({summary['measured']} tools measured)")
    if args.json:
        pathlib.Path(args.json).write_text(json.dumps(report, indent=2,
                                                      sort_keys=True) + "\n")
        print(f"wrote {args.json}")
    return 0 if summary["fail"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
