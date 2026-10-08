"""Everything the domain-architecture ladder grades, from the sandbox plus the image.

The calls that are not pure computation are the HMMER scans: one per declared
cutoff, which is four. The declared panel keeps each to about half a minute
rather than the seven minutes a full Pfam-A pass costs on a corpus of this size.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.domain_architecture.core import Corpus, ArchitectureError


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["annotations_input"],
                         params["proteins_input"])
    rules = corpus.rules

    resource = rules.resources["pfam_a_hmm"].verify()
    if not resource["matches"]:
        raise ArchitectureError(
            f"{resource['name']} at {resource['path']} hashes to "
            f"{resource['found_sha256'][:16]} and the campaign declares "
            f"{resource['declared_sha256'][:16]}. A domain architecture against a "
            "different Pfam release is not this campaign's gold.")

    hits_by_cutoff = {cutoff: corpus.hits(cutoff) for cutoff in rules.cutoffs}
    hits = hits_by_cutoff[rules.baseline_cutoff]
    table = corpus.architecture_table(hits)

    return {
        "mibig_release": corpus.release,
        "module_type": rules.module_type,
        "baseline_cutoff": rules.baseline_cutoff,
        "baseline_resolution": rules.baseline_resolution,
        "baseline_span": rules.baseline_span,
        "baseline_mapping": rules.baseline_mapping,
        "panel_size": len(rules.panel),
        "hmmer_version": corpus.hmmer_version(),
        "pfam_version": resource["version"],
        "pfam_sha256": resource["found_sha256"],

        # ---- R1: how the corpus was arrived at, and what the curated side holds
        "selection_census": {
            **corpus.annotation_census["annotation_census"],
            **{f"drop_{k}": v for k, v in
               corpus.annotation_census["selection_drops"].items()},
            "shipped_proteins": corpus.shipped_proteins,
            "corpus_genes": len(corpus.keys),
            "corpus_entries": len({b for b, _ in corpus.keys}),
        },
        "curated_coordinate_usability":
            corpus.annotation_census["curated_coordinate_usability"],
        "curated_active_flag": corpus.annotation_census["curated_active_flag"],
        "curated_module_census": corpus.curated_module_census(),
        "curated_domain_totals": corpus.curated_domain_totals(),
        "protein_length_summary": corpus.length_summary(),
        "panel_census": {
            model: sum(1 for hs in hits.values() for h in hs if h.model == model)
            for model in rules.panel},

        # ---- R2: the architecture itself
        "architecture_table": table,
        "architecture_census": corpus.architecture_census(table),
        "overlap_census": corpus.overlap_census(hits),

        # ---- R3: against the curated modules
        **corpus.reconcile(hits),
        **corpus.module_decomposition(hits),

        # ---- R4: the declared readings
        "policy_sweep": corpus.policy_sweep(hits_by_cutoff),
        "mapping_sweep": corpus.mapping_sweep(hits),
    }
