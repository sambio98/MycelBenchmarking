"""Everything the remote-homology ladder grades, from the sandbox plus the image.

The calls that are not pure computation are the searches: one over the declared
query batch, and one per declared shuffle window over the same batch mutated.
Four runs at this ladder's size, each against the same indexed corpus.
"""

from __future__ import annotations

import pathlib

import yaml

from npbench_c.templates.remote_homology.core import Corpus


def load_params(campaign: pathlib.Path) -> dict:
    return yaml.safe_load(
        (pathlib.Path(campaign) / "task.yaml").read_text())["template_params"]


def analyse(campaign: pathlib.Path) -> dict:
    campaign = pathlib.Path(campaign)
    params = load_params(campaign)
    corpus = Corpus.load(campaign, params["corpus_input"], params["labels_input"])
    rules = corpus.rules

    queries = corpus.queries()
    baseline = corpus.search({accession: corpus.sequences[accession]
                              for accession in sorted(queries.values())})

    census = {accession: corpus.band_census(baseline, accession)
              for accession in sorted(queries.values())}

    return {
        "label_resource": "InterPro 110.0",
        "corpus_superfamily": rules.corpus_superfamily,
        "sensitivity": rules.sensitivity,
        "evalue": rules.evalue,
        "max_seqs": rules.max_seqs,
        "threads": rules.threads,
        "format_fields": list(rules.format_fields),
        "mmseqs_version": corpus.mmseqs_version(),

        # ---- R1: the corpus as parsed, the declared queries, and the search
        "corpus_sequences": corpus.size,
        "corpus_sha256": corpus.sequences_sha256(),
        "query_table": {
            entry: {"accession": accession,
                    "length": len(corpus.sequences[accession]),
                    "family_members": sum(
                        1 for entries in corpus.families.values() if entry in entries),
                    "family_name": corpus.entry_names.get(entry, "")}
            for entry, accession in sorted(queries.items())},
        "baseline_hits": {accession: sum(
            1 for h in baseline if h.query == accession and h.target != accession)
            for accession in sorted(queries.values())},
        "baseline_band_sizes": corpus.baseline_band_sizes(baseline),

        # ---- R2: what the search found of what the resource asserts
        "recovery": {accession: corpus.recovery(baseline, accession)
                     for accession in sorted(queries.values())},
        "band_census": census,

        # ---- R3: the declared verdict per query and band
        "band_verdicts": {
            accession: {band: corpus.verdict(cells[band])
                        for band in rules.band_labels}
            for accession, cells in sorted(census.items())},
        "verdict_census": _verdict_census(corpus, census),

        # ---- R4: the shuffle ladder
        "shuffle_ladder": corpus.ladder(baseline),
    }


def _verdict_census(corpus: Corpus, census: dict) -> dict[str, int]:
    counts = {verdict: 0 for verdict in corpus.rules.verdicts}
    for cells in census.values():
        for band in corpus.rules.band_labels:
            counts[corpus.verdict(cells[band])] += 1
    return counts
