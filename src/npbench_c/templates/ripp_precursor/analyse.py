"""Everything the RiPP-precursor ladder grades, from the sandbox plus the image.

Two calls are not pure computation, and both are HMMER: the declared panel over
the corpus as shipped, and the same panel over the corpus with each declared
leader removed. The declared panel keeps that to seconds rather than the
half-minute a full Pfam-A pass costs even on sequences this short.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.ripp_precursor.core import Corpus, PrecursorError


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["annotations_input"],
                         params["reference_records_input"])
    rules = corpus.rules

    resource = rules.resources["pfam_a_hmm"].verify()
    if not resource["matches"]:
        raise PrecursorError(
            f"{resource['name']} at {resource['path']} hashes to "
            f"{resource['found_sha256'][:16]} and the campaign declares "
            f"{resource['declared_sha256'][:16]}. A domain verdict against a "
            "different Pfam release is not this campaign's gold.")

    table = corpus.localisation_table()
    before = corpus.domains(corpus.fasta())
    after = corpus.domains(corpus.fasta(table))

    return {
        "mibig_release": corpus.release,
        "biosynthetic_class": rules.biosynthetic_class,
        "baseline_policy": rules.baseline_policy,
        "threshold": rules.threshold,
        "panel_size": len(rules.panel),
        "hmmer_version": corpus.hmmer_version(),
        "pfam_version": resource["version"],
        "pfam_sha256": resource["found_sha256"],

        # ---- R1: the selection, and the shapes the field arrives in
        "selection_census": {
            "entries_in_class": corpus.entries_in_class,
            "entries_with_a_named_precursor": corpus.entries_with_a_named_precursor,
            "entries_excluded_by_status": corpus.entries_excluded_by_status,
            "entries_absent_from_reference_set":
                corpus.entries_absent_from_reference_set,
            "entries_excluded_by_both": corpus.entries_excluded_by_both,
            "reference_records_shipped": corpus.reference_records_shipped,
            "records_named": corpus.records_named,
            "records_with_an_unresolved_gene": corpus.records_unresolved_gene,
            "corpus_entries": len({r.bgc for r in corpus.records}),
            "corpus_records": len(corpus.records),
        },
        "core_sequence_shape_census": corpus.shape_census(),
        "core_sequence_case_census": corpus.case_census(),
        "records_declaring_a_cleavage_location":
            len(corpus.declaring_a_cleavage_location()),
        "ripp_type_census": corpus.ripp_type_census(),
        "status_census": corpus.status_census(),
        "precursor_length_summary": corpus.precursor_length_summary(),

        # ---- R2: where the core sits
        "localisation_table": table,
        "localisation_census": corpus.localisation_census(table),
        "leader_length_summary": corpus.leader_length_summary(table),

        # ---- R3: who says so -- the coordinate, and the domain model
        **corpus.reconcile(table),
        **corpus.domain_audit(table, before),

        # ---- R4: the two counterfactual axes
        "policy_sweep": corpus.policy_sweep(),
        **corpus.deletion_outcome(table, before, after),
    }
