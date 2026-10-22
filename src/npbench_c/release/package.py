"""Build the release tree, and refuse to build one that carries the answers.

This is the other side of the sandbox. `sweep.runner.build_sandbox` assembles
what ONE system sees for ONE campaign and audits it; this assembles what leaves
the building, and the failure mode is the same defect at a different scale: a
release that ships `gold/` hands over every answer at once.

Three things it excludes that are less obvious than `gold/`.

**`docs/`.** The build log and the catalogue are full of graded values, and that
is deliberate: this project's rule for a measurement that motivated a design
decision is that it belongs in `docs/`, "which the sandbox never sees". True of a
sandbox, and false of a release that shipped `docs/` -- which would publish the
verdict censuses, the recovery rates, the call counts and every other number the
rungs grade. The rule only holds while `docs/` stays out of the release.

**`src/npbench_c/templates/`.** Each template's `analyse.py` and `build.py` ARE
the oracle: they compute gold from the inputs. Shipping them lets a participant
regenerate every answer, tool runs and all. A participant needs the contract,
not the machinery, so no source ships at all.

**`grading.json` and `internal_sweep.json`.** The first says which fields are
graded and with what tolerance, which is a map of where to spend effort rather
than an answer, and the second carries measured per-rung difficulty. Both are
organiser material.

What it adds, and checks: a manifest with a digest per file, and an attribution
file assembled from each campaign's own licence statements. That assembly turned
up something worth stating: attribution lives in `task.yaml`'s `license_notice`
for eleven campaigns and in `inputs/provenance.json` for twelve, and nothing
before this checked that every campaign has it in at least one of the two. A
campaign with neither is a refusal here, because several inputs are CC BY and
attribution is a licence condition rather than a courtesy.

Usage:
  python -m npbench_c.release.package <output_dir> [--campaign NAME ...]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import shutil
import sys

import yaml

from npbench_c.grading.version import grader_version
from npbench_c.readiness import cold_sandbox
from npbench_c.readiness.gate import _resolve_gold
from npbench_c.sweep.runner import (
    NEVER_IN_SANDBOX,
    PUBLIC_CAMPAIGN_PATHS,
    SECRET_CAMPAIGN_PATHS,
    GoldLeak,
)

#: Repository trees that must never appear in a participant release, each for a
#: reason the module docstring states. Named rather than filtered by guesswork,
#: so an auditor can check the list against the reasons.
EXCLUDED_TREES = ("docs", "src", "tests", "image", "audit", "oracle")

MANIFEST = "MANIFEST.json"
ATTRIBUTION = "ATTRIBUTION.md"
README = "README.md"


class ReleaseError(RuntimeError):
    """The release cannot be built as asked, or would carry what it must not."""


def _digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def campaigns(repo: pathlib.Path) -> list[pathlib.Path]:
    """Every campaign a release can ship: one that can also be SCORED.

    A task.yaml is not enough. A campaign under construction has its inputs and
    its rules and no grading spec yet, and shipping that hands a system a task
    with nothing behind it -- no rungs to clear and no score to earn. The
    criterion is mechanical rather than a status flag anyone has to remember to
    set: no grading.json, not a campaign yet. `unfinished` names the ones left
    out so the manifest can say so rather than letting a reader infer it from a
    count.
    """
    return [c for c in sorted((repo / "campaigns").iterdir())
            if (c / "task.yaml").is_file() and (c / "grading.json").is_file()]


def unfinished(repo: pathlib.Path) -> list[pathlib.Path]:
    """Campaign directories with a task.yaml and no grading spec."""
    return [c for c in sorted((repo / "campaigns").iterdir())
            if (c / "task.yaml").is_file() and not (c / "grading.json").is_file()]


#: A key whose NAME states a licence or an attribution. "source" is not on the
#: list and the reason is a false pass this check produced on its first run:
#: construct-ecoli-prna-01 cleared it on the strength of `"source": "PubMed"`
#: inside an evidence record, which says where a claim came from, not what may
#: be redistributed. A rule that accepts that accepts almost any JSON.
LICENCE_KEYS = ("licen", "attribut", "copyright", "citation")

#: Carried along once a licence is found nearby, never sufficient alone.
SUPPORTING_KEYS = ("source", "url", "doi", "accession", "database", "retrieved")


def _licence_subtrees(node, path: str = "") -> dict:
    """Every part of a provenance document that states a licence.

    Walks rather than reading top-level keys, which is the second thing this
    check got wrong: `pocket-rebh-01` and `residues-prna-01` state their
    licences inside `annotation_source` and `structure_sources`, so a shallow
    scan declared two correctly-attributed campaigns unattributed. A licence
    nested two levels down is still a licence.
    """
    found: dict = {}
    if isinstance(node, dict):
        names = {k.lower() for k in node}
        if any(word in name for name in names for word in LICENCE_KEYS):
            keep = {k: v for k, v in sorted(node.items())
                    if not isinstance(v, (dict, list))
                    and any(word in k.lower()
                            for word in LICENCE_KEYS + SUPPORTING_KEYS)}
            if keep:
                found[path or "."] = keep
        for key, value in sorted(node.items()):
            found.update(_licence_subtrees(value, f"{path}/{key}" if path else key))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            found.update(_licence_subtrees(value, f"{path}[{index}]"))
    return found


def _attribution_for(campaign: pathlib.Path) -> dict:
    """Every licence statement this campaign makes, from both places it may live.

    `task.yaml`'s `license_notice` is the documented home and eleven campaigns
    use it. Twelve others state the same thing in `inputs/provenance.json`,
    which is shipped and so carries the attribution with the data. Both are
    read; neither is assumed.
    """
    task = yaml.safe_load((campaign / "task.yaml").read_text())
    out: dict = {"campaign": campaign.name,
                 "license_class": task.get("license_class"),
                 "from_task_yaml": task.get("license_notice"),
                 "from_provenance": None}
    provenance = campaign / "inputs" / "provenance.json"
    if provenance.is_file():
        try:
            data = json.loads(provenance.read_text())
        except json.JSONDecodeError as exc:
            raise ReleaseError(
                f"{campaign.name}: inputs/provenance.json is not readable JSON "
                f"({exc}), so its licence statement cannot be collected") from exc
        out["from_provenance"] = _licence_subtrees(data) or None
    if not out["from_task_yaml"] and not out["from_provenance"]:
        raise ReleaseError(
            f"{campaign.name} states no licence or attribution in either "
            "task.yaml (license_notice) or inputs/provenance.json, and it ships "
            "data. Several inputs here are CC BY, where attribution is a "
            "condition rather than a courtesy, so the release refuses rather "
            "than publishing unattributed material")
    return out


def _render_attribution(entries: list[dict]) -> str:
    lines = [
        "# Attribution",
        "",
        "Every campaign in this release, the licence class it declares, and the",
        "licence statements it carries. Assembled from each campaign's own files",
        "rather than written here, so this file cannot drift from them: a",
        "campaign's statement lives in `task.yaml` (`license_notice`) or in",
        "`inputs/provenance.json` beside the data, and both are collected.",
        "",
    ]
    for entry in entries:
        lines.append(f"## {entry['campaign']}")
        lines.append("")
        lines.append(f"- declared licence class: `{entry['license_class']}`")
        if entry["from_task_yaml"]:
            lines.append("- from `task.yaml`:")
            lines.append("")
            lines.append("```yaml")
            lines.append(yaml.safe_dump({"license_notice": entry["from_task_yaml"]},
                                        sort_keys=True, width=78).rstrip())
            lines.append("```")
        if entry["from_provenance"]:
            lines.append("- from `inputs/provenance.json`:")
            lines.append("")
            lines.append("```json")
            lines.append(json.dumps(entry["from_provenance"], indent=2,
                                    sort_keys=True))
            lines.append("```")
        lines.append("")
    return "\n".join(lines) + "\n"


README_TEXT = """# NPBench-C release

One directory per campaign, each holding exactly what a system is given:

- `inputs/` the task's data
- `reference/` every pinned table the task requires
- `task.yaml` the objective, the output contract, the counting rules and the
  closed vocabularies
- `submission/` empty, for the answer

A campaign is solvable from its own directory: every table, vocabulary, counting
rule and convention its grading depends on is inside `reference/` or `task.yaml`.
If something needed is missing, that is a defect in the campaign rather than a
puzzle to solve.

What is deliberately absent, so nobody looks for it: gold, the grading spec, the
oracle that computes gold, the per-rung difficulty measurements, and this
project's build log. The scoring runs on the organiser's side.

`MANIFEST.json` carries a digest for every file here, the grader version the
campaigns were built against, and the campaign list. `ATTRIBUTION.md` carries
each campaign's licence statements, assembled from the campaigns themselves.
"""


def build_release(repo: pathlib.Path, out: pathlib.Path,
                  only: tuple[str, ...] = ()) -> dict:
    """Assemble the release tree at `out` and return its manifest."""
    repo = pathlib.Path(repo).resolve()
    out = pathlib.Path(out).resolve()
    if out == repo or repo in out.parents:
        raise ReleaseError(
            f"the release would be written inside the build repo ({out}), where "
            "the next packaging run would find it and package it again")
    chosen = campaigns(repo)
    if only:
        names = {c.name for c in chosen}
        unknown = sorted(set(only) - names)
        if unknown:
            raise ReleaseError(f"no such campaign: {unknown}")
        chosen = [c for c in chosen if c.name in only]
    if not chosen:
        raise ReleaseError(
            "no campaign carries both a task.yaml and a grading.json; nothing "
            "to release")

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    for campaign in chosen:
        target = out / "campaigns" / campaign.name
        target.mkdir(parents=True)
        for name in PUBLIC_CAMPAIGN_PATHS:
            source = campaign / name
            if source.is_dir():
                shutil.copytree(source, target / name)
            elif source.is_file():
                shutil.copyfile(source, target / name)
        (target / "submission").mkdir()

    entries = [_attribution_for(campaign) for campaign in chosen]
    (out / ATTRIBUTION).write_text(_render_attribution(entries))
    (out / README).write_text(README_TEXT)

    skipped = [c.name for c in unfinished(repo) if not only or c.name in only]
    manifest = {
        "grader_version": grader_version(),
        "campaigns": [campaign.name for campaign in chosen],
        "omitted_campaigns": {
            name: "under construction: a task.yaml and inputs but no grading "
                  "spec, so nothing in it can be scored yet"
            for name in skipped},
        "excluded_trees": list(EXCLUDED_TREES),
        "withheld_campaign_paths": list(SECRET_CAMPAIGN_PATHS),
        "files": {},
    }
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != MANIFEST:
            manifest["files"][str(path.relative_to(out))] = _digest(path)
    manifest["file_count"] = len(manifest["files"])
    (out / MANIFEST).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    audit_release(repo, out, chosen)
    return manifest


def audit_release(repo: pathlib.Path, out: pathlib.Path,
                  chosen: list[pathlib.Path]) -> None:
    """Refuse a release that carries gold, the oracle, or a graded answer.

    Runs after assembly rather than instead of it, so the check sees the tree
    that would actually ship -- including the manifest and attribution this
    module adds, which are the files a path-based rule would never think to
    examine.
    """
    for path in sorted(out.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(out)
        if rel.name in NEVER_IN_SANDBOX:
            raise GoldLeak(f"{rel} is in the release")
        parts = rel.parts
        if parts[0] in EXCLUDED_TREES:
            raise GoldLeak(f"{rel} is in an excluded tree")
        if parts[0] == "campaigns" and len(parts) > 2:
            if parts[2] in SECRET_CAMPAIGN_PATHS:
                raise GoldLeak(f"{rel} is withheld campaign material")

    for campaign in chosen:
        packaged = out / "campaigns" / campaign.name
        for name in PUBLIC_CAMPAIGN_PATHS:
            if (campaign / name).exists() and not (packaged / name).exists():
                raise ReleaseError(
                    f"{campaign.name}: {name} did not reach the release, so the "
                    "campaign is not solvable from it")
        if not (packaged / "submission").is_dir():
            raise ReleaseError(f"{campaign.name}: no empty submission/ directory")
        if any((packaged / "submission").iterdir()):
            raise ReleaseError(f"{campaign.name}: submission/ is not empty")

        # The same scans the readiness gate runs as checks 22 and 23, pointed
        # at the release instead of a sandbox: gold and the grading spec come
        # from the build repo, the prose from what would ship.
        if not (campaign / "gold" / "gold.json").is_file():
            continue
        hits = cold_sandbox.findings(campaign, packaged, _resolve_gold)
        if hits:
            raise GoldLeak(
                f"{campaign.name}: the release states {len(hits)} graded "
                f"answer(s): {sorted(hits)[:4]}")

        # And check 23, over THIS release's campaigns rather than every built
        # one. A disclosure of a neighbour's answer only matters where both
        # reach the same system, so the set being shipped is the right scope --
        # a release that drops a campaign should not refuse over a sentence
        # about it, and one that adds a campaign must be rechecked against it.
        neighbours = [c for c in chosen
                      if c.name != campaign.name
                      and (c / "gold" / "gold.json").is_file()
                      and (c / "grading.json").is_file()]
        cross = cold_sandbox.cross_findings(
            campaign, packaged, neighbours, _resolve_gold)
        if cross:
            raise GoldLeak(
                f"{campaign.name}: the release states {len(cross)} graded "
                f"answer(s) belonging to a campaign shipped with it: "
                f"{sorted(cross)[:4]}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", help="directory to write the release into")
    ap.add_argument("--repo", default=".", help="build repository root")
    ap.add_argument("--campaign", action="append", default=[],
                    help="package only this campaign; repeatable")
    args = ap.parse_args(argv)

    manifest = build_release(pathlib.Path(args.repo), pathlib.Path(args.out),
                             tuple(args.campaign))
    print(f"release written to {args.out}")
    print(f"  grader {manifest['grader_version']}")
    print(f"  {len(manifest['campaigns'])} campaigns, "
          f"{manifest['file_count']} files")
    if manifest["omitted_campaigns"]:
        print(f"  omitted (no grading spec yet): "
              f"{', '.join(sorted(manifest['omitted_campaigns']))}")
    print(f"  withheld: {', '.join(manifest['withheld_campaign_paths'])}")
    print(f"  excluded trees: {', '.join(manifest['excluded_trees'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
