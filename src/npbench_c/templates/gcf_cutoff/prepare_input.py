"""Ship the clustering corpus and the chemistry-derived families for it.

Usage:
  python -m npbench_c.templates.gcf_cutoff.prepare_input <campaign> \
      --source <dir of antiSMASH-processed .gbk> --provenance <provenance.json>

Run once at instantiation. The corpus and its designed families both come from a
declared selection rule recorded in the source provenance, so the campaign ships
the files rather than re-deriving them, and records the rule and the rule's own
caveats alongside.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import shutil


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--source", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--key", default="bigscape_input")
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    inputs = campaign / "inputs"
    target = inputs / "clusters"
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    source = pathlib.Path(args.source)
    names = sorted(p.name for p in source.glob("*.gbk"))
    if not names:
        raise RuntimeError(f"{source} holds no .gbk files")
    for name in names:
        shutil.copyfile(source / name, target / name)

    block = json.loads(pathlib.Path(args.provenance).read_text())[args.key]
    families = {name: sorted(members)
                for name, members in sorted(block["families"].items())}
    controls = sorted(block["controls"])
    members = sorted({m for v in families.values() for m in v} | set(controls))
    stems = sorted(n[: -len(".gbk")] for n in names)
    if members != stems:
        raise RuntimeError(
            "the declared families plus controls must partition the corpus "
            f"exactly; corpus {len(stems)}, declared {len(members)}")

    (inputs / "designed_families.json").write_text(json.dumps({
        "families": families,
        "controls": controls,
        "selection_rule": block["selection_rule"],
        "excluded_placeholder_names": block["excluded_placeholder_names"],
        "excluded_placeholder_reason": block["excluded_placeholder_reason"],
        "note": ("The families are derived from CHEMISTRY -- entries sharing an "
                 "InChIKey connectivity block in the chemical-space campaign's "
                 "table -- and not from any clustering of these sequences. That "
                 "is what makes them something to reconcile against rather than "
                 "a restatement of the answer. The controls are one entry per "
                 "biosynthetic class that belongs to no such group."),
    }, indent=2, sort_keys=True) + "\n")

    digests = {name: hashlib.sha256((target / name).read_bytes()).hexdigest()
               for name in names}
    provenance = {
        "source": block["source"],
        "url": block["url"],
        "license": block["license"],
        "attribution": ("MIBiG antiSMASH-processed reference set "
                        f"({block['url']}), licensed {block['license']}"),
        "archive_sha256": block["archive_sha256"],
        "antismash_version_that_processed_them":
            block["antismash_version_that_processed_them"],
        "clusters": len(names),
        "files": {"clusters": "clusters/", "designed_families":
                  "designed_families.json"},
        "sha256": digests,
        "note": ("Every input in this corpus was processed by the SAME antiSMASH "
                 "version, which is what the cross-version rule requires. This "
                 "campaign never compares these files against output from the "
                 "image's antiSMASH, so no run mixes two versions."),
    }
    (inputs / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"clusters/: {len(names)} antiSMASH-processed GBKs")
    print(f"designed_families.json: {len(families)} families, "
          f"{len(controls)} controls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
