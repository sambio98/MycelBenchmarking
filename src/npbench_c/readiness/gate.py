"""The readiness gate: may this campaign leave the building?

A campaign is not ready because someone believes it is; it is ready when this
returns all-pass. Checks that require internal agent runs cannot be satisfied
mechanically and report PENDING rather than PASS -- the gate never launders an
unmeasured property into a green tick. Collaborator attention is the scarcest
input in the project, so nothing reaches an auditor until the mechanical checks
are green and the agent-run checks have real numbers.
"""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
from dataclasses import dataclass

from npbench_c.grading.grade import grade, grade_json
from npbench_c.sweep.prompt import CONVENTIONS_KEY
from npbench_c.sweep.runner import _absolute_pythonpath
from npbench_c.grading.ladder import MAX_CHANCE_FLOOR
from npbench_c.grading.version import grader_version

PURITY_REPEATS = 100
FORBIDDEN_SOURCES = ("kegg", "metacyc", "biocyc", "drugbank")
VENDOR_TOKENS = ("mycel",)

PASS, FAIL, PENDING = "PASS", "FAIL", "PENDING"

#: The gold-side files each campaign's build may read, with the reason each one
#: is not in the sandbox. Outside campaigns/, where no sandbox can reach it.
GOLD_SIDE_INPUTS = (pathlib.Path(__file__).resolve().parents[3]
                    / "audit" / "gold_side_inputs.json")

#: What a campaign may declare about the licence of the data it ships.
#:
#:   open         no downstream condition beyond attribution (CC BY, CC0, public
#:                domain). Most of the benchmark.
#:   nc           non-commercial source, so unavailable to commercial evaluators.
#:   share_alike  the shipped data is an ADAPTATION of a ShareAlike source, so it
#:                carries that licence onward. The obligation attaches to the
#:                campaign's own data files, not to the benchmark: a collection of
#:                separable campaigns is not an adaptation of any one of them, so
#:                the code, the grader and the other campaigns are unaffected and
#:                such a campaign can be dropped without touching the rest.
LICENSE_CLASSES = ("open", "nc", "share_alike")

AGENT_RUN_CHECKS = (
    "monotonicity_measured",
    "r1_clear_rate",
    "leakage_no_tool_ablation",
    "difficulty_gate",
    "discrimination",
)


@dataclass
class Check:
    name: str
    status: str
    detail: str


def _check(name: str, ok: bool, detail: str) -> Check:
    return Check(name, PASS if ok else FAIL, detail)


def _oracle_and_purity(campaign: pathlib.Path) -> list[Check]:
    import yaml

    checks = []
    sub = campaign / "oracle_submission"
    measurement = sub / "measurement.json"

    # The measure command comes from the campaign's own task.yaml, exactly as it
    # does for the sweep runner. The gate previously globbed for a *.faa input,
    # which is a construct-design assumption that broke the moment a campaign
    # had different inputs.
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    harness = task.get("harness") or {}
    measure_command = harness.get("measure_command")
    if not measure_command:
        return [Check("oracle_runs", FAIL,
                      "task.yaml declares no harness.measure_command")]

    # A campaign built from a shared template carries no oracle directory at
    # all, so the solve step is whatever task.yaml declares. The legacy form --
    # a campaign-local oracle/solve.py -- is still honoured.
    solve_command = harness.get("solve_command") or [sys.executable, "oracle/solve.py",
                                                     "{submission}"]
    env = _absolute_pythonpath(os.environ)
    try:
        subprocess.run(
            [tok.format(campaign=str(campaign), submission=str(sub)) for tok in solve_command],
            cwd=campaign, env=env, check=True, capture_output=True, timeout=1800,
        )
        subprocess.run(
            [tok.format(campaign=str(campaign), submission=str(sub),
                        measurement=str(measurement)) for tok in measure_command],
            cwd=campaign, env=env, check=True, capture_output=True, timeout=1800,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
        return [Check("oracle_runs", FAIL, f"oracle pipeline failed: {exc}")]

    measurement_data = json.loads(measurement.read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    spec = json.loads((campaign / "grading.json").read_text())
    result = grade(measurement_data, gold, spec)

    checks.append(_check("oracle_scores_1.0", result["score"] == 1.0,
                         f"score={result['score']} depth={result['depth']}"))
    checks.append(_check("oracle_ladder_monotonic", result["monotonic"],
                         f"rungs_cleared={result['rungs_cleared']} depth={result['depth']}"))
    # The cap applies to the discriminating rungs (R2 upward). R1 is an
    # execution precondition whose design target is an 85-95% clear rate, so
    # including it would make the cap unsatisfiable for every well-formed ladder.
    checks.append(_check("chance_floor_within_limit",
                         result["discriminating_chance_floor"] <= MAX_CHANCE_FLOOR,
                         f"discriminating={result['discriminating_chance_floor']} "
                         f"(full={result['chance_floor']}) limit={MAX_CHANCE_FLOOR}"))

    reference = grade_json(measurement_data, gold, spec)
    identical = all(grade_json(measurement_data, gold, spec) == reference
                    for _ in range(PURITY_REPEATS))
    checks.append(_check("grader_purity", identical,
                         f"{PURITY_REPEATS} repeat gradings byte-identical={identical}"))
    checks.append(_check("grader_version_recorded",
                         spec.get("campaign_id") is not None,
                         f"grader_version={grader_version()}"))

    # The pinned version must match the grader that just ran. Without this the
    # campaign can silently carry a stale hash: a grader change is a benchmark
    # version bump, and a score whose grader_version does not identify the code
    # that produced it is not a comparable score.
    import yaml

    task = yaml.safe_load((campaign / "task.yaml").read_text())
    pinned, current = task.get("grader_version"), grader_version()
    checks.append(_check("grader_version_matches_pin", pinned == current,
                         f"pinned={pinned} current={current}"))
    return checks


def gold_side_inputs(campaign_name: str,
                     path: pathlib.Path | None = None) -> dict[str, str]:
    """The gold-side files one campaign's build may read, as {file: reason}."""
    source = path or GOLD_SIDE_INPUTS
    if not source.is_file():
        return {}
    declared = json.loads(source.read_text()).get("inputs") or {}
    return {entry["file"]: entry["reason"]
            for entry in declared.get(campaign_name, [])}


def _rebuild_gold(campaign: pathlib.Path, into: pathlib.Path,
                  restore: tuple[str, ...]) -> tuple[bool, str]:
    """Rebuild gold in a tree holding only the sandbox plus `restore`."""
    import shutil

    import yaml

    from npbench_c.sweep.runner import build_sandbox

    task = yaml.safe_load((campaign / "task.yaml").read_text())
    gold_command = ((task.get("harness") or {}).get("gold_command")
                    or [sys.executable, "oracle/generate_gold.py"])
    sandbox = build_sandbox(campaign, into, None)
    (sandbox / "gold").mkdir(exist_ok=True)
    for name in restore:
        shutil.copyfile(campaign / "gold" / name, sandbox / "gold" / name)
    try:
        subprocess.run([tok.format(campaign=str(sandbox)) for tok in gold_command],
                       cwd=sandbox, env=_absolute_pythonpath(os.environ),
                       check=True, capture_output=True, timeout=7200)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
        return False, str(exc)[-200:]
    produced = sandbox / "gold" / "gold.json"
    if not produced.is_file():
        return False, "the build wrote no gold.json"
    if produced.read_bytes() != (campaign / "gold" / "gold.json").read_bytes():
        return False, "the rebuilt gold differs from the committed gold"
    return True, "byte-identical"


def _oracle_determinism(campaign: pathlib.Path) -> Check:
    """Regenerate gold from the SANDBOX and require byte-identity.

    This was a rebuild inside the campaign directory, which has everything --
    so it could say the build is deterministic and could not say whether a
    value a system is asked for depends on a file the system never sees. The
    rebuild happens in a tree holding only the published sandbox, which answers
    both with one build: 21 of 25 campaigns reproduce their gold from the
    sandbox alone, and four need exactly one declared file --
    audit/gold_side_inputs.json names them, and none of the four is graded.

    Each declared file also has to be needed: the build must refuse without it,
    which costs one extra rebuild on those four campaigns and nothing on the
    other twenty-one. A declaration nothing requires is a hole in the check.
    """
    import tempfile

    declared = gold_side_inputs(campaign.name)
    restore = tuple(sorted(declared))
    with tempfile.TemporaryDirectory() as tmp:
        ok, detail = _rebuild_gold(campaign, pathlib.Path(tmp) / "cold", restore)
        if not ok:
            plus = f" plus {', '.join(restore)}" if restore else ""
            return Check("oracle_determinism", FAIL,
                         f"rebuilt from the sandbox{plus}: {detail}")
        if restore:
            bare, _ = _rebuild_gold(campaign, pathlib.Path(tmp) / "bare", ())
            if bare:
                return Check("oracle_determinism", FAIL,
                             f"declared gold-side input(s) {', '.join(restore)} "
                             "are not needed: the sandbox alone rebuilds this gold")
    where = ("the sandbox alone" if not restore
             else f"the sandbox plus {', '.join(restore)}, each of them required")
    return Check("oracle_determinism", PASS,
                 f"gold rebuilt from {where}, byte-identical")


def _contract_and_provenance(campaign: pathlib.Path) -> list[Check]:
    import yaml  # declared dependency; not in the scoring path

    checks = []
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    spec = json.loads((campaign / "grading.json").read_text())

    tiers = [r.get("gold_tier") for r in spec["rungs"]]
    checks.append(_check("gold_tier_declared_per_rung", all(tiers),
                         f"tiers={tiers}"))
    checks.append(_check("at_most_one_g2_rung", tiers.count("G2") <= 1,
                         f"G2 rungs={tiers.count('G2')}"))
    checks.append(_check("no_g3_outside_repro_audit",
                         "G3" not in tiers or task.get("family") == "reproducibility_audit",
                         f"tiers={tiers} family={task.get('family')}"))
    license_class = task.get("license_class")
    checks.append(_check("license_class_declared",
                         license_class in LICENSE_CLASSES,
                         f"license_class={license_class} "
                         f"(allowed: {sorted(LICENSE_CLASSES)})"))

    # A share-alike campaign carries an obligation that travels with its data, so
    # the obligation has to travel with the campaign. A notice naming the source,
    # the licence and what a redistributor must do is the whole of it -- and it is
    # checked, because a licence condition recorded only in a design document is a
    # condition the person redistributing the files will never see.
    notice = task.get("license_notice") or {}
    required = ("source", "license", "attribution", "redistribution")
    missing = [k for k in required if not notice.get(k)]
    checks.append(_check(
        "share_alike_notice_present",
        license_class != "share_alike" or not missing,
        "not a share-alike campaign" if license_class != "share_alike"
        else (f"license_notice missing: {missing}" if missing
              else f"notice names {', '.join(required)}")))

    blob = json.dumps(task).lower()
    bad = [s for s in FORBIDDEN_SOURCES if s in blob]
    checks.append(_check("no_license_landmine_sources", not bad,
                         f"forbidden sources referenced: {bad}" if bad else "none"))

    contract = task.get("output_contract") or {}
    checks.append(_check("output_contract_present", bool(contract),
                         f"{len(contract)} artifacts contracted"))

    # Every pinned table the task requires must be inside a path the sandbox
    # ships, or the campaign is unsolvable and no run can tell you so: the
    # agent simply fails and looks incapable.
    from npbench_c.sweep.runner import PUBLIC_CAMPAIGN_PATHS

    pinned = task.get("constraints", {}).get("pinned_tables") or []
    pinned = [pinned] if isinstance(pinned, str) else pinned
    reachable = [p for p in pinned
                 if pathlib.PurePath(p).parts[0] in PUBLIC_CAMPAIGN_PATHS
                 and (campaign / p).exists()]
    checks.append(_check("pinned_tables_reachable_by_agent",
                         bool(pinned) and len(reachable) == len(pinned),
                         f"{len(reachable)}/{len(pinned)} pinned tables are under a "
                         f"public sandbox path and present on disk"))

    # A campaign may need reference data it cannot ship -- too large, or under a
    # licence the benchmark is not in a position to pass on. Such a resource is
    # declared rather than copied, and the declaration has to be checkable: the
    # file must be where the campaign says, and must hash to what it recorded.
    # An unverifiable external resource is the solvability defect one step out:
    # the sandbox looks complete and the answer depends on bytes nobody pinned.
    external = task.get("constraints", {}).get("external_resources") or []
    external = [external] if isinstance(external, dict) else list(external)
    resource_detail, resources_ok = "no external resources declared", True
    if external:
        named = [e for e in external if e.get("name") and e.get("declared_in")]
        found, mismatched, absent_files = [], [], []
        for entry in named:
            declared_in = campaign / entry["declared_in"]
            if not declared_in.is_file():
                absent_files.append(entry["name"])
                continue
            spec = json.loads(declared_in.read_text())
            block = (spec.get("external_resources") or {}).get(entry["name"])
            if not block:
                absent_files.append(entry["name"])
                continue
            prefix = pathlib.Path(
                os.environ.get("NPBENCH_TOOL_PREFIX", "/opt/npbench-tools/bin"))

            # Two declared shapes, because a resource is not always one file. A
            # single file carries `relative_path` plus a full `sha256`; a set of
            # files carries `relative_paths` plus a `fingerprint` hashed over them
            # in the declared order. Either way the declaration has to be
            # checkable, or it is a comment rather than a pin.
            paths, expected = None, None
            if block.get("relative_path") and block.get("sha256"):
                paths, expected = [block["relative_path"]], block["sha256"]
            elif block.get("relative_paths") and block.get("fingerprint"):
                paths, expected = list(block["relative_paths"]), block["fingerprint"]
            if not paths:
                absent_files.append(entry["name"])
                continue

            digest, resolved = hashlib.sha256(), True
            for relative in paths:
                # A path may name the interpreter version with a glob, since the
                # image's site-packages directory carries it.
                matches = (sorted((prefix.parent).glob(relative))
                           if any(c in relative for c in "*?[")
                           else [prefix.parent / relative])
                matches = [m for m in matches if m.is_file()]
                if not matches:
                    matches = sorted(
                        (prefix.parent / "lib").glob(
                            f"python3.*/site-packages/antismash/{relative}"))
                    matches = [m for m in matches if m.is_file()]
                if not matches:
                    resolved = False
                    break
                with open(matches[-1], "rb") as handle:
                    for chunk in iter(lambda: handle.read(1 << 20), b""):
                        digest.update(chunk)
            if not resolved:
                absent_files.append(entry["name"])
                continue
            got = digest.hexdigest()
            if got == expected or got[:len(expected)] == expected:
                found.append(entry["name"])
            else:
                mismatched.append(entry["name"])
        resources_ok = (len(named) == len(external)
                        and len(found) == len(external))
        parts = [f"{len(found)}/{len(external)} declared external resources "
                 "present and matching their recorded sha256"]
        if mismatched:
            parts.append(f"fingerprint mismatch: {mismatched}")
        if absent_files:
            parts.append(f"not found or not fully declared: {absent_files}")
        if len(named) != len(external):
            parts.append("an entry is missing name or declared_in")
        resource_detail = "; ".join(parts)
    checks.append(_check("external_resources_pinned_and_present", resources_ok,
                         resource_detail))

    # Every rung must carry the key so the declaration is explicit; R1 has no
    # earlier rung to start from, so null is the correct value there. Bundles
    # above R1 must exist on disk: a missing bundle silently turns a warm run
    # into a cold one and the per-rung difficulty profile becomes a fiction.
    rungs = sorted(task.get("rungs", []), key=lambda r: r["ordinal"])
    declared = [("warm_start_bundle" in r) for r in rungs]
    bundles_ok = all(
        (r.get("warm_start_bundle") is None) if i == 0
        else bool(r.get("warm_start_bundle"))
        and (campaign / r["warm_start_bundle"]).is_dir()
        for i, r in enumerate(rungs)
    )
    checks.append(_check("warm_start_declared_per_rung",
                         bool(rungs) and all(declared) and bundles_ok,
                         f"{sum(declared)}/{len(rungs)} rungs declare the key; "
                         f"bundles present and non-empty above R1: {bundles_ok}"))

    checks.append(_check("excluded_fields_justified",
                         all("reason" in e for e in task.get("excluded_from_grading", [])),
                         f"{len(task.get('excluded_from_grading', []))} exclusions, each with a measurement"))
    return checks


import yaml  # declared dependency; never in the scoring path


def _resolve_gold(gold: object, path: str) -> object:
    """Walk a grading spec's dotted gold path. Raises on anything unresolvable."""
    for part in path.split("."):
        if part.endswith("]"):
            name, index = part[:-1].split("[", 1)
            gold = gold[name][int(index)]       # type: ignore[index]
        else:
            gold = gold[part]                   # type: ignore[index]
    return gold


def _warm_bundles_withhold_answers(campaign: pathlib.Path) -> Check:
    """A rung may not be handed the answer to one of its own components.

    A warm start exists to remove the EARLIER rungs' work, so if a rung's bundle
    already contains what that rung is graded on, the component belongs to an
    earlier rung. The failure is invisible in an oracle run -- the oracle scores
    1.0 either way -- and only shows up when a real system clears a rung it
    should not have, by which point the per-rung difficulty profile is fiction.

    Only COMPOSITE gold values are checked. A scalar drawn from a declared
    vocabulary ("direct_repeat", "unchanged") legitimately appears throughout an
    earlier rung's census, so flagging it would be a false positive in every
    enum-scored component.
    """
    import tempfile

    from npbench_c.sweep.runner import build_sandbox

    spec = json.loads((campaign / "grading.json").read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    bundles = {r["rung_id"]: r.get("warm_start_bundle")
               for r in task.get("rungs", [])}

    hits, checked = [], 0
    with tempfile.TemporaryDirectory() as tmp:
        for rung in spec["rungs"]:
            bundle = bundles.get(rung["rung_id"])
            if not bundle:
                continue
            sandbox = build_sandbox(campaign, pathlib.Path(tmp) / rung["rung_id"],
                                    bundle)
            handed = [p for p in sorted(sandbox.iterdir()) if p.suffix == ".json"]
            blob = json.dumps([json.loads(p.read_text()) for p in handed],
                              sort_keys=True)
            for component in rung["components"]:
                reference = component["args"].get("gold")
                if not isinstance(reference, str) or not reference.startswith("gold:"):
                    continue
                try:
                    value = _resolve_gold(gold, reference.split(":", 1)[1])
                except (KeyError, IndexError, TypeError, ValueError):
                    hits.append(f"{rung['rung_id']}:{component['name']}:unresolvable")
                    continue
                if not isinstance(value, (dict, list)):
                    continue
                checked += 1
                if json.dumps(value, sort_keys=True) in blob:
                    hits.append(f"{rung['rung_id']}:{component['name']}")
    return _check("warm_bundles_withhold_rung_answers", not hits,
                  f"{checked} composite components checked against their rung's "
                  f"bundle" + (f"; handed over: {hits}" if hits else "; none handed over"))


def _cold_sandbox_withholds_answers(campaign: pathlib.Path) -> Check:
    """The cold sandbox may not state a graded number in prose.

    The warm-bundle check guards one rung against its own bundle. This guards
    every rung against the files the agent gets before it has done anything:
    `task.yaml` and `reference/`. Those exist to publish the objective, the
    contract and the conventions -- and a note explaining a design decision drifts
    into quoting the measurement that motivated it, at which point R1 can be
    cleared by reading the task description.

    It is easy to do and invisible afterwards: the oracle still scores 1.0, and the
    no-tool ablation still fails because a stub does not read prose. Only a real
    system benefits, and then the per-rung profile is fiction.

    Only DISTINCTIVE scalars are checked. A graded count of 2 or a rate of 0.5
    appears in any prose by coincidence, so the threshold is four or more digits
    for an integer and four or more decimals for a float -- enough that a match is
    the number rather than a collision. Composites are checked whole, as in the
    warm-bundle case.
    """
    import re as _re
    import tempfile

    from npbench_c.sweep.runner import build_sandbox

    spec = json.loads((campaign / "grading.json").read_text())
    gold = json.loads((campaign / "gold" / "gold.json").read_text())

    def distinctive(value) -> list[str]:
        """The string forms of a value worth searching for, or none."""
        if isinstance(value, bool) or value is None:
            return []
        if isinstance(value, int):
            return [str(value), f"{value:,}"] if abs(value) >= 1000 else []
        if isinstance(value, float):
            text = repr(value)
            _, _, decimals = text.partition(".")
            return [text] if len(decimals) >= 4 else []
        if isinstance(value, str):
            return []
        if isinstance(value, (dict, list)):
            # An empty composite serialises to "[]" or "{}", which occurs in any
            # YAML or JSON file by coincidence. It also carries no answer to
            # disclose, so there is nothing to search for.
            return [json.dumps(value, sort_keys=True)] if value else []
        return []

    hits, checked = [], 0
    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(campaign, pathlib.Path(tmp) / "cold", None)
        # task.yaml and reference/ only. inputs/ is the data the campaign asks the
        # agent to read, so a number appearing there is not a leak -- a shipped
        # record's length or a curated coordinate table is the input, and a
        # component whose value can be read straight out of inputs/ without doing
        # the work is a free component rather than a disclosure.
        parts = []
        for path in sorted(sandbox.rglob("*")):
            if not path.is_file() or path.suffix not in (
                    ".yaml", ".yml", ".json", ".md", ".txt"):
                continue
            if path.relative_to(sandbox).parts[0] == "inputs":
                continue
            parts.append(path.read_text(errors="replace"))
        blob = "\n".join(parts)
        # Digit grouping varies, so the haystack is also searched with separators
        # stripped: "4,552" and "4552" are the same leak.
        plain = _re.sub(r"(?<=\d)[,_](?=\d)", "", blob)
        for rung in spec["rungs"]:
            for component in rung["components"]:
                reference = component["args"].get("gold")
                if not isinstance(reference, str) or not reference.startswith("gold:"):
                    continue
                try:
                    value = _resolve_gold(gold, reference.split(":", 1)[1])
                except (KeyError, IndexError, TypeError, ValueError):
                    continue
                for needle in distinctive(value):
                    checked += 1
                    # A bare integer needs word boundaries: 2359 occurs inside
                    # "BGC0002359" in a shipped accession table, which is a
                    # coincidence rather than a disclosure. A JSON composite is
                    # distinctive enough to match as a substring.
                    if needle.startswith(("[", "{")):
                        found = needle in blob or needle in plain
                    else:
                        pattern = r"(?<![0-9.])" + _re.escape(needle) + r"(?![0-9])"
                        found = bool(_re.search(pattern, blob)
                                     or _re.search(pattern, plain))
                    if found:
                        hits.append(f"{rung['rung_id']}:{component['name']}")
                        break
    return _check(
        "cold_sandbox_withholds_answers", not hits,
        f"{checked} distinctive graded values searched in task.yaml and "
        f"reference/" + (f"; stated there: {sorted(set(hits))}" if hits
                         else "; none stated"))


def _cold_sandbox_withholds_small_answers(campaign: pathlib.Path) -> Check:
    """Check 22: the two shapes check 21's distinctiveness threshold cannot reach.

    Check 21 above searches for DISTINCTIVE graded values -- four digits, four
    decimals -- because a graded count of 2 appears in any prose by coincidence.
    That threshold is right for what it does and it leaves a gap, which an audit
    found twelve answers sitting in: a small scalar stated in a sentence
    (`15` as a locus total, `456` as a corpus count, `0.62` as a GC floor), and a
    declared term stated against its own label, where neither half of the pair is
    distinctive but the pair is the answer.

    So this check anchors on the component rather than on the digit count, which
    reaches those at the cost of real collisions. Collisions are declared in
    `audit/cold_sandbox_exemptions.json`, which sits outside `campaigns/` because
    an exemption names the declared value an answer equals and so states the
    answer. Every exemption must still collide: a stale one is a check nobody is
    getting, and this project produced one within a minute of writing its first
    waiver list.
    """
    import tempfile

    from npbench_c.readiness import cold_sandbox
    from npbench_c.sweep.runner import build_sandbox

    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(campaign, pathlib.Path(tmp) / "cold", None)
        hits = cold_sandbox.findings(campaign, sandbox, _resolve_gold)
        stale = cold_sandbox.stale_exemptions(campaign, sandbox, _resolve_gold)
        waived = cold_sandbox.exemptions(campaign.name)

    problems = ([f"stated: {h}" for h in sorted(hits)]
                + [f"stale exemption: {s}" for s in stale])
    detail = (f"{len(waived)} declared collision(s)"
              if waived else "no declared collisions")
    if problems:
        detail += "; " + "; ".join(problems)
    else:
        detail += "; no small graded value or label/term pair stated"
    return _check("cold_sandbox_withholds_small_answers", not problems, detail)


def _cold_sandbox_withholds_other_campaigns_answers(campaign: pathlib.Path) -> Check:
    """Check 23: the cold sandbox may not state ANOTHER campaign's answer either.

    Checks 21 and 22 compare a sandbox against its own gold, so a note naming a
    neighbouring campaign's graded answer passes both every time. Two instances
    shipped: the malleobactin campaign's audit notes stated the sevadicin
    instantiation's verdict while explaining that the field varies across
    instantiations, and a codon-table provenance block added to four
    construct-design campaigns stated the gene-calling campaign's graded GC
    fraction while recording which record the weights came from. Both were true,
    useful sentences about this campaign that happened to answer another one.

    The scope is the set of campaigns shipped together, because that is when the
    disclosure can reach a system: every built campaign beside this one. The
    neighbours are read from `campaigns/` directly rather than through a sandbox
    apiece -- a sandbox copies the inputs, and 22 of them per gate run would cost
    minutes to learn vocabularies that `reference/` already states.

    Which shapes it runs is measured rather than assumed, and the small-scalar
    shape is excluded: see the table in `cold_sandbox`'s docstring.
    """
    import tempfile

    from npbench_c.readiness import cold_sandbox
    from npbench_c.sweep.runner import build_sandbox

    others = [c for c in sorted(campaign.parent.iterdir())
              if c.is_dir() and c.name != campaign.name
              and (c / "grading.json").is_file()
              and (c / "gold" / "gold.json").is_file()]
    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(campaign, pathlib.Path(tmp) / "cold", None)
        hits = cold_sandbox.cross_findings(campaign, sandbox, others, _resolve_gold)
        stale = cold_sandbox.stale_cross_exemptions(
            campaign, sandbox, others, _resolve_gold)
        waived = cold_sandbox.cross_exemptions(campaign.name)

    problems = ([f"states: {h}" for h in sorted(hits)]
                + [f"stale exemption: {s}" for s in stale])
    detail = f"{len(others)} neighbouring campaign(s) checked"
    if waived:
        detail += f", {len(waived)} declared collision(s)"
    detail += ("; " + "; ".join(problems) if problems
               else "; no other campaign's answer stated")
    return _check("cold_sandbox_withholds_other_campaigns_answers",
                  not problems, detail)


#: Design-label fields a campaign may not declare. `structural_elements` was one
#: for twenty-three campaigns and nothing read it; the census that retired it is
#: in docs. The synonyms are listed because the field's defect was its
#: unverifiability, not its name, and a replacement under another name would
#: inherit it.
UNVALIDATED_DESIGN_LABELS = ("structural_elements", "design_elements",
                             "structural_features", "design_features",
                             "design_properties")

#: A campaign directory name: lower-case words joined by - or _, ending in the
#: two-digit instantiation number.
CAMPAIGN_NAME = re.compile(r"\b[a-z][a-z0-9]*(?:[-_][a-z0-9][a-z0-9._]*)+-\d{2}\b")


def _design_claims_are_checkable(campaign: pathlib.Path) -> Check:
    """Check 24: a design claim in task.yaml must be one a reader can check.

    `structural_elements` was a free-text list of design features --
    `verification`, `counterfactual`, `control`, `planted`, `abstention`,
    `abstention_free`, `recovery` -- that no code read and no check validated.
    The census that retired it: `verification` and `counterfactual` were declared
    by 23 of 23 campaigns, which is what the four-rung ladder means, so they
    carried no information. `planted` was declared by 11, and no campaign ships
    an altered input: the two that explained the claim describe something the
    campaign DOES at run time (a declared single-base deletion, a declared query
    shuffle), and every campaign that wrote out its reasoning concluded "nothing
    is planted" while nine others claimed it with no explanation at all.
    `control` was declared by 14 and meant a different thing in each -- a
    negative subset, two chains that must agree, a clean comparison section --
    with no common signature to check. `recovery` was what one campaign measures.

    And one of them was answer-adjacent, which is what settled it: on the
    constraint control pair, the conflict half declared `abstention` and the
    feasible half did not, while the graded R2 answer of each is exactly
    `infeasible_gc_unreachable` against `feasible`. The label told a reader which
    half of the pair they were holding.

    So the field is gone, this check keeps it gone under any of its names, and
    what it gestured at lives where something reads it: the ladder in
    `grading.json`, the controls in gold where a component grades them, the
    abstention verdict in the graded vocabulary, the planting rule in the
    campaign's own rules file.

    The other half of the check is the design claims that DO resolve: a campaign
    name mentioned in `task.yaml` must be a campaign that exists, and
    `control_pair_with` must name one that names this campaign back. A dangling
    reference is how a renamed or retired campaign leaves a sentence pointing at
    nothing, and a one-way pairing is how a control pair becomes a pair of
    unrelated campaigns.
    """
    import yaml

    text = (campaign / "task.yaml").read_text()
    task = yaml.safe_load(text)
    problems = []

    for field in UNVALIDATED_DESIGN_LABELS:
        if field in task:
            problems.append(
                f"{field} is declared ({task[field]!r}); it is a design label no "
                "check can verify, and on a control pair it states which half "
                "the reader is holding")

    known = {c.name for c in campaign.parent.iterdir() if c.is_dir()}
    # A campaign names itself in its own id; what matters is the names of others.
    mentioned = set(CAMPAIGN_NAME.findall(text)) - {campaign.name}
    for name in sorted(mentioned - known):
        problems.append(f"names {name}, which is not a campaign")

    pair = task.get("control_pair_with")
    if pair is not None:
        if pair not in known:
            problems.append(f"control_pair_with={pair} is not a campaign")
        else:
            other = yaml.safe_load((campaign.parent / pair / "task.yaml").read_text())
            if other.get("control_pair_with") != campaign.name:
                problems.append(
                    f"control_pair_with={pair}, which names "
                    f"{other.get('control_pair_with')!r} back")

    detail = (f"{len(mentioned & known)} reference(s) to other campaigns resolve"
              + (f", paired with {pair}" if pair else "")
              + ("; " + "; ".join(problems) if problems
                 else "; no unvalidated design label"))
    return _check("design_claims_are_checkable", not problems, detail)


#: The provenance field that says which antiSMASH processed a shipped record set.
#: A field rather than the archive's name, because a check can read a field.
ANTISMASH_PROCESSING_KEY = "antismash_version_that_processed_them"

#: Matches a snapshot tag or a provenance source/url that names antiSMASH.
NAMES_ANTISMASH = re.compile(r"antismash", re.IGNORECASE)


def _antismash_versions_are_not_mixed(campaign: pathlib.Path) -> Check:
    """Check 25: one antiSMASH per campaign, and it must be written down.

    antiSMASH's own rules forbid comparing regions from two of its versions, and
    this benchmark has the ingredients for exactly that mistake sitting in the
    repository. MIBiG publishes its reference set as
    `mibig_antismash_4.0_gbk_as8b1.tar.bz2`, processed with antiSMASH **8.0 beta
    1**; the image pins **8.0.4**. Four campaigns ship records out of that set and
    one campaign runs antiSMASH, and the open item that led here asked whether to
    reprocess the set or to declare that both sides come from the published one.

    Measured over all 23 campaigns, the two groups do not overlap: the four that
    ship the processed set drive HMMER or BiG-SCAPE and never invoke antiSMASH,
    and the one that invokes antiSMASH reads a raw INSDC chromosome and MIBiG's
    curated loci. So nothing is mixed today, reprocessing would buy no campaign
    anything, and it would cost the provenance of a published, citable set. What
    was missing was anything that notices if that changes.

    So this is the rule, enforced per campaign:

    - a campaign that ships antiSMASH-processed records must RECORD which
      antiSMASH processed them, in `inputs/provenance.json` under
      `antismash_version_that_processed_them`, where a check can read it;
    - a campaign that also runs antiSMASH must run the version that processed
      them. If it does not, the campaign is making the forbidden comparison and
      has two ways out: reprocess the records with the pinned antiSMASH, or take
      both sides from the published set and stop running the tool;
    - a declared antiSMASH pin must be the version the registry pins, so a
      campaign cannot claim a tool the image does not carry.
    """
    import yaml

    from npbench_c.tools.registry import ANTISMASH

    text = (campaign / "task.yaml").read_text()
    task = yaml.safe_load(text)
    snapshots = task.get("snapshots_required") or {}
    pins = [tag for tag in (snapshots.get("image") or [])
            if str(tag).lower().startswith("antismash@")]
    shipped = [tag for tag in (snapshots.get("shipped") or [])
               if NAMES_ANTISMASH.search(str(tag).split("@")[0])]

    # Only `source` and `url` are read, never a free-text note: the campaign that
    # RUNS antiSMASH explains in its note that a sliced record is one antiSMASH
    # refuses, and a scan over prose would read that as shipping processed data.
    processed_fields: list[str] = []
    recorded: list[str] = []
    provenance = campaign / "inputs" / "provenance.json"
    if provenance.is_file():
        def walk(node: object, path: str = "") -> None:
            if isinstance(node, dict):
                for field in ("source", "url"):
                    value = node.get(field)
                    if isinstance(value, str) and NAMES_ANTISMASH.search(value):
                        processed_fields.append(f"{path or '.'}:{field}")
                value = node.get(ANTISMASH_PROCESSING_KEY)
                if isinstance(value, str) and value.strip():
                    recorded.append(value)
                for key, child in node.items():
                    walk(child, f"{path}.{key}")
            elif isinstance(node, list):
                for index, child in enumerate(node):
                    walk(child, f"{path}[{index}]")
        walk(json.loads(provenance.read_text()))

    problems = []
    ships_processed = bool(shipped or processed_fields)
    if ships_processed and not recorded:
        problems.append(
            f"ships antiSMASH-processed records ({sorted(set(shipped + processed_fields))}) "
            f"and records no {ANTISMASH_PROCESSING_KEY} in inputs/provenance.json, "
            "so nothing can check which antiSMASH produced them")

    for pin in pins:
        version = pin.split("@", 1)[1]
        if version != ANTISMASH.pinned_version:
            problems.append(
                f"declares {pin} and the registry pins "
                f"antismash {ANTISMASH.pinned_version}")
        for record in recorded:
            if version not in record:
                problems.append(
                    f"runs antiSMASH {version} and ships records processed by "
                    f"{record.split(',')[0]!r}. That is the cross-version "
                    "comparison antiSMASH's rules forbid: either reprocess the "
                    "records with the pinned antiSMASH, or take both sides from "
                    "the published set and stop running the tool here")

    detail = (f"{len(pins)} antiSMASH pin(s), {len(set(shipped + processed_fields))} "
              f"processed-record reference(s), {len(recorded)} version record(s)")
    detail += ("; " + "; ".join(problems) if problems
               else "; no campaign mixes two antiSMASH versions")
    return _check("antismash_versions_are_not_mixed", not problems, detail)


def _single_term_answers_are_not_stated(campaign: pathlib.Path) -> Check:
    """Check 26: a component whose whole answer is one declared term, stated.

    Measured gap. Check 22's shapes pair a graded term with the label in its gold
    path or inside the graded value; neither reaches a term stated against a
    DIFFERENT label. The constraint control pair shipped `expect: feasible` as a
    template parameter -- R3's entire answer, in the sandbox -- and two more
    campaigns stated their single-term verdict in prose. All three passed checks
    21, 22 and 23 every time they ran.
    """
    import tempfile

    from npbench_c.readiness import cold_sandbox
    from npbench_c.sweep.runner import build_sandbox

    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(campaign, pathlib.Path(tmp) / "cold", None)
        hits = cold_sandbox.single_term_findings(campaign, sandbox, _resolve_gold)
        stale = cold_sandbox.stale_single_term_exemptions(
            campaign, sandbox, _resolve_gold)
        waived = cold_sandbox.single_term_exemptions(campaign.name)

    problems = ([f"stated: {h}" for h in sorted(hits)]
                + [f"stale exemption: {s}" for s in stale])
    detail = (f"{len(waived)} declared collision(s)" if waived
              else "no declared collisions")
    detail += ("; " + "; ".join(problems) if problems
               else "; no single-term answer stated outside its vocabulary")
    return _check("single_term_answers_are_not_stated", not problems, detail)

def _contract_fields(declared: object) -> dict:
    """{field: nested-spec-or-None} for one artefact's contract, recursively.

    A contract key may name several fields at once (`a / b / c`), which is how
    the suite declares a group of fields that share one description.
    """
    out: dict = {}
    if not isinstance(declared, dict):
        return out
    for key, description in declared.items():
        if key == CONVENTIONS_KEY:
            continue
        child = _contract_fields(description) if isinstance(description, dict) else None
        for name in (part.strip() for part in str(key).split("/")):
            if name:
                out[name] = child
    return out


def _project(value: object, spec: dict | None) -> object:
    """`value` with only the fields `spec` declares. A declared leaf keeps all
    of its subtree: the contract stops naming fields where the answer becomes
    data."""
    if spec is None or not isinstance(value, dict):
        return value
    return {k: _project(v, spec[k]) for k, v in value.items() if k in spec}


def _field_disagreements(written: object, spec: dict, where: str) -> list[str]:
    """Every name one side has and the other does not, at every declared depth."""
    if not isinstance(written, dict):
        return [f"{where} is not an object, so its fields cannot be declared"]
    out = [f"{where}.{k} written, not declared"
           for k in sorted(written) if k not in spec]
    out += [f"{where}.{k} declared, not written"
            for k in sorted(spec) if k not in written]
    for name, child in sorted(spec.items()):
        if child and name in written:
            out += _field_disagreements(written[name], child, f"{where}.{name}")
    return out


def _contract_names_what_is_graded(campaign: pathlib.Path) -> Check:
    """Check 27: an answer written to the PUBLISHED contract scores what the
    oracle scores.

    The oracle does not go through the contract. `solve.py` renders a submission
    from gold's own key tuples, so a contract that named a field the measurer
    does not read -- or failed to name one it does -- still scored 1.0 while
    being unusable by any system that followed it. Measured gap, five of
    twenty-five campaigns: two whose field names were never updated when their
    template was generalised (a system following them scored zero at R1), two
    whose R3 and R4 answers live under wrapper keys the contract never mentioned,
    and one whose field name was abbreviated into something unwritable.

    So: project the oracle's own submission onto the field names the contract
    declares, and require the projection to grade exactly as the oracle does.

    The check reaches as deep as the contract names keys. Where a contract stops
    naming keys and describes the shape below in prose -- a table's columns, a
    map keyed by an identifier -- the projection keeps that subtree whole, so a
    column misnamed in prose is not something this can see.

    The baseline is the committed measurement rather than a second solve run:
    `_oracle_and_purity` regenerates it earlier in the same gate run and asserts
    it scores 1.0, so running the oracle again here would measure nothing new.
    """
    import tempfile

    import yaml

    task = yaml.safe_load((campaign / "task.yaml").read_text())
    contract = task.get("output_contract") or {}
    if not contract:
        return Check("contract_names_what_is_graded", FAIL,
                     "task.yaml declares no output_contract")
    measure_command = (task.get("harness") or {}).get("measure_command")
    if not measure_command:
        return Check("contract_names_what_is_graded", FAIL,
                     "task.yaml declares no harness.measure_command")

    oracle = campaign / "oracle_submission"
    written = sorted(p.name for p in oracle.iterdir()
                     if p.is_file() and p.name != "measurement.json")
    problems = [f"{name} written, no contract entry" for name in written
                if name not in contract]
    problems += [f"{name} declared, not written" for name in sorted(contract)
                 if name not in written]

    gold = json.loads((campaign / "gold" / "gold.json").read_text())
    spec = json.loads((campaign / "grading.json").read_text())
    baseline = grade(json.loads((oracle / "measurement.json").read_text()), gold, spec)

    with tempfile.TemporaryDirectory() as tmp:
        projected = pathlib.Path(tmp) / "submission"
        projected.mkdir()
        for name in written:
            fields = _contract_fields(contract.get(name))
            text = (oracle / name).read_text()
            try:
                raw = json.loads(text)
            except json.JSONDecodeError:
                # Not every artefact is JSON -- an ORF ships as FASTA, whose
                # contract is conventions only. Nothing to project.
                (projected / name).write_text(text)
                continue
            if not fields:
                problems.append(f"{name} declares no field, so it cannot be written")
                (projected / name).write_text(text)
                continue
            problems += _field_disagreements(raw, fields, name)
            (projected / name).write_text(
                json.dumps(_project(raw, fields), indent=2) + "\n")

        measurement = pathlib.Path(tmp) / "measurement.json"
        try:
            subprocess.run(
                [tok.format(campaign=str(campaign), submission=str(projected),
                            measurement=str(measurement)) for tok in measure_command],
                cwd=campaign, env=_absolute_pythonpath(os.environ),
                check=True, capture_output=True, timeout=1800,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
            return Check("contract_names_what_is_graded", FAIL,
                         f"measuring the contract-shaped submission failed: {exc}")
        result = grade(json.loads(measurement.read_text()), gold, spec)

    same = (result["depth"], result["score"]) == (baseline["depth"], baseline["score"])
    detail = (f"contract-shaped submission: depth={result['depth']} "
              f"score={result['score']} (oracle depth={baseline['depth']} "
              f"score={baseline['score']})")
    if problems:
        detail += "; " + "; ".join(problems[:6])
        if len(problems) > 6:
            detail += f"; and {len(problems) - 6} more"
    return _check("contract_names_what_is_graded", same and not problems, detail)


def _every_key_name_is_findable(campaign: pathlib.Path) -> Check:
    """Check 28: a key name the answer requires, that no sandbox file states.

    Check 27 reaches as deep as the contract names keys. Below that a contract
    described its shape in prose -- "entries carrying curated modules, modules,
    modules naming more than one gene" -- and never gave the key names, which are
    load-bearing: measured on arch-pks, renaming one key inside
    `selection_census` takes an otherwise perfect answer from depth 4 to depth 0.
    227 such names across 18 campaigns were findable nowhere a system can read.

    A name the campaign's own DATA supplies -- a product name the tool emitted, a
    status the records carry -- is declared in audit/contract_key_names.json
    instead, because a contract that listed those would be listing an answer.
    """
    import tempfile

    from npbench_c.readiness import cold_sandbox
    from npbench_c.sweep.runner import build_sandbox

    with tempfile.TemporaryDirectory() as tmp:
        sandbox = build_sandbox(campaign, pathlib.Path(tmp) / "cold", None)
        required = cold_sandbox.key_names_required(campaign)
        words = cold_sandbox._sandbox_words(sandbox)
    waived = cold_sandbox.key_name_exemptions(campaign.name)
    unfindable = {name for name in required if name not in words}
    undeclared = sorted(unfindable - set(waived))
    stale = sorted(name for name in waived if name in words)

    problems = ([f"unfindable: {name}" for name in undeclared]
                + [f"stale exemption: {name}" for name in stale])
    detail = (f"{len(required)} key name(s) required, "
              f"{len(waived)} declared as supplied by the data")
    detail += ("; " + "; ".join(problems[:6]) if problems
               else "; every one of them stated in the sandbox")
    if len(problems) > 6:
        detail += f"; and {len(problems) - 6} more"
    return _check("every_key_name_is_findable", not problems, detail)


def _integrity(campaign: pathlib.Path) -> list[Check]:
    hits = []
    for path in sorted(campaign.rglob("*")):
        if not path.is_file() or path.suffix not in (".py", ".yaml", ".json", ".md"):
            continue
        text = path.read_text(errors="replace").lower()
        for token in VENDOR_TOKENS:
            if token in text:
                hits.append(f"{path.relative_to(campaign)}:{token}")
    return [_check("no_vendor_tool_names", not hits, f"hits={hits}" if hits else "none")]


def _pending_agent_runs(campaign: pathlib.Path) -> list[Check]:
    """Report measured values from an internal sweep if present, else PENDING."""
    sweep = campaign / "internal_sweep.json"
    if not sweep.is_file():
        return [Check(name, PENDING, "no internal_sweep.json; requires internal agent runs")
                for name in AGENT_RUN_CHECKS]
    data = json.loads(sweep.read_text())

    # A stub-based sweep is a harness self-test. Its systems are built to reach
    # preset depths, so promoting its numbers to PASS would launder a fixture
    # into evidence -- the exact thing this gate exists to prevent.
    stub_based = bool(data.get("_meta", {}).get("stub_based"))
    stubs = data.get("_meta", {}).get("stub_systems", [])

    out = []
    for name in AGENT_RUN_CHECKS:
        if name not in data:
            out.append(Check(name, PENDING, "absent from internal_sweep.json"))
            continue
        entry = data[name]
        detail = str(entry.get("detail", ""))
        if stub_based:
            out.append(Check(name, PENDING,
                             f"stub-based sweep (fixtures: {', '.join(stubs)}); "
                             f"real systems required. Harness result: {detail}"))
        else:
            out.append(Check(name, PASS if entry.get("pass") else FAIL, detail))
    return out


def run(campaign_dir: pathlib.Path) -> list[Check]:
    campaign = pathlib.Path(campaign_dir).resolve()
    checks: list[Check] = []
    checks += _oracle_and_purity(campaign)
    checks.append(_oracle_determinism(campaign))
    checks += _contract_and_provenance(campaign)
    checks.append(_warm_bundles_withhold_answers(campaign))
    checks.append(_cold_sandbox_withholds_answers(campaign))
    checks.append(_cold_sandbox_withholds_small_answers(campaign))
    checks.append(_cold_sandbox_withholds_other_campaigns_answers(campaign))
    checks.append(_design_claims_are_checkable(campaign))
    checks.append(_antismash_versions_are_not_mixed(campaign))
    checks.append(_single_term_answers_are_not_stated(campaign))
    checks.append(_contract_names_what_is_graded(campaign))
    checks.append(_every_key_name_is_findable(campaign))
    checks += _integrity(campaign)
    checks += _pending_agent_runs(campaign)
    return checks


def report(campaign_dir: pathlib.Path) -> dict:
    checks = run(campaign_dir)
    counts = {s: sum(1 for c in checks if c.status == s) for s in (PASS, FAIL, PENDING)}
    return {
        "campaign_dir": str(campaign_dir),
        "grader_version": grader_version(),
        "ready_to_audit": counts[FAIL] == 0 and counts[PENDING] == 0,
        "mechanical_checks_green": counts[FAIL] == 0,
        "summary": counts,
        "checks": [{"name": c.name, "status": c.status, "detail": c.detail} for c in checks],
    }


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print("usage: python -m npbench_c.readiness.gate <campaign_dir>", file=sys.stderr)
        return 2
    rep = report(pathlib.Path(argv[0]))
    width = max(len(c["name"]) for c in rep["checks"])
    for c in rep["checks"]:
        mark = {"PASS": "  ok  ", "FAIL": " FAIL ", "PENDING": " pend "}[c["status"]]
        print(f"[{mark}] {c['name']:<{width}}  {c['detail']}")
    s = rep["summary"]
    print(f"\n{s['PASS']} pass / {s['FAIL']} fail / {s['PENDING']} pending")
    print(f"mechanical checks green: {rep['mechanical_checks_green']}")
    print(f"ready to send for audit: {rep['ready_to_audit']}")
    return 0 if rep["mechanical_checks_green"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
