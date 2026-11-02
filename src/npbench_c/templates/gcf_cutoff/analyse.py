"""Everything the GCF-cutoff ladder grades, from the sandbox plus the image.

One BiG-SCAPE run covers the whole declared cutoff grid, so the whole ladder --
the baseline partition, the reconciliation and the sweep -- rests on a single
tool invocation and the cutoffs are comparable by construction.
"""

from __future__ import annotations

import json
import pathlib

import yaml

from npbench_c.templates.gcf_cutoff.core import RULES, Clustering, CutoffError


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def declared(campaign: pathlib.Path, key: str) -> object:
    spec = json.loads((pathlib.Path(campaign) / "reference" / RULES).read_text())
    return spec["clustering"][key]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    clustering = Clustering.load(campaign, params["cluster_dir"],
                                 params["designed_families_input"])
    rules = clustering.rules

    pfam = rules.resources["pfam_a_hmm"]
    found = pfam.resolve()
    from npbench_c.tools.cache import file_digest

    digest = file_digest(found)
    if digest != pfam.sha256:
        raise CutoffError(
            f"{pfam.name} at {found} hashes to {digest[:16]} and the campaign "
            f"declares {pfam.sha256[:16]}. A partition computed against a "
            "different Pfam release is not this campaign's gold.")

    assignments = clustering.assignments()
    baseline_label = rules.label(rules.baseline_cutoff)
    baseline = clustering.partition(assignments[baseline_label])
    sweep = clustering.sweep(assignments)

    return {
        "corpus_size": len(clustering.corpus),
        "designed_families": {name: len(members) for name, members
                              in sorted(clustering.designed.items())},
        "controls": list(clustering.controls),
        "bigscape_version": clustering.bigscape_version(),
        "pfam_version": pfam.version,
        "pfam_sha256": digest,
        "cutoffs": [rules.label(c) for c in rules.cutoffs],
        "baseline_cutoff": baseline_label,
        "include_gbk": rules.include_gbk,
        "class_census": clustering.class_census(assignments[baseline_label]),
        **clustering.group_census(baseline),
        "baseline_partition": baseline,
        **clustering.reconciliation(baseline),
        **clustering.pair_statistics(baseline),
        "cutoff_sweep": sweep,
        **clustering.sweep_summary(sweep),
    }
