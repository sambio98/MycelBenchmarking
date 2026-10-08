"""Deterministic stub systems for the domain-architecture ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each scans the declared Pfam panel with the pinned HMMER.

  reader      clears R1: the selection census, the curated record's own shape,
              the per-family hit census and the tool identity
  architect   clears R2: the architecture of each gene after resolution
  reconcile   clears R3: against MIBiG's curated modules
  complete    clears R4
  unresolved  COMPETENCE probe, tooled: scans correctly and reports every hit as
              a domain, in coordinate order, with no overlap resolution. Clears
              R1 and fails R2.
  noncompute  LEAKAGE ablation: never runs the tool and never opens the inputs.

`unresolved` is the probe this campaign needs, and it is the mistake the
catalogued design itself made. Pfam's nested families hit the same residues, so
listing every hit produces an architecture with a ketosynthase core spliced into
a thiolase and five methyltransferases in a row -- longer than the real one,
internally consistent, and wrong on hundreds of genes. Nothing in the output
says a domain was counted twice.
"""

from __future__ import annotations

import argparse
import json
import pathlib

from npbench_c.templates.domain_architecture import analyse as dar
from npbench_c.templates.domain_architecture.core import Corpus

LEVELS = ("reader", "architect", "reconcile", "complete", "unresolved",
          "noncompute")
REPORT = "architecture_report.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", required=True, choices=LEVELS)
    ap.add_argument("--submission", required=True)
    args = ap.parse_args()

    sandbox = pathlib.Path.cwd()
    out = pathlib.Path(args.submission)
    out.mkdir(parents=True, exist_ok=True)
    target = out / REPORT

    def write(payload: dict) -> int:
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return 0

    if args.level == "noncompute":
        return write({
            # Family-agnostic on purpose: this level serves both
            # instantiations and must not name either one's vocabulary, because
            # it is the ablation that reads nothing at all.
            "mibig_release": "", "module_family": "",
            "primary_domain": "", "secondary_domain": "",
            "baseline_cutoff": "", "baseline_resolution": "",
            "baseline_span": "", "baseline_mapping": "", "panel_size": 0,
            "hmmer_version": "HMMER 3.4", "pfam_version": "Pfam 35.0",
            "pfam_sha256": "0" * 64,
            "selection_census": {}, "curated_coordinate_usability": {},
            "curated_active_flag": {}, "curated_module_census": {},
            "strata_census": {}, "domains_by_module_type": {},
            "curated_domain_totals": {}, "protein_length_summary": {},
            "panel_census": {},
            "architecture_table": [], "architecture_census": {},
            "overlap_census": {},
            "reconciliation_table": [], "primary_verdicts": {},
            "secondary_verdicts": {}, "secondary_verdicts_by_stratum": {},
            "primary_delta_census": {},
            "module_table": [], "module_verdicts": {},
            "complete_modules_total": 0,
            "policy_sweep": {}, "mapping_sweep": {},
        })

    params = dar.load_params(sandbox)
    corpus = Corpus.load(sandbox, params["annotations_input"],
                         params["proteins_input"])
    rules = corpus.rules
    pfam = rules.resources["pfam_a_hmm"]
    from npbench_c.tools.cache import file_digest

    hits = corpus.hits()
    report = {
        "mibig_release": corpus.release,
        "module_family": rules.module_family,
        "primary_domain": rules.primary_domain,
        "secondary_domain": rules.secondary_domain,
        "baseline_cutoff": rules.baseline_cutoff,
        "baseline_resolution": rules.baseline_resolution,
        "baseline_span": rules.baseline_span,
        "baseline_mapping": rules.baseline_mapping,
        "panel_size": len(rules.panel),
        "hmmer_version": corpus.hmmer_version(),
        "pfam_version": pfam.version,
        "pfam_sha256": file_digest(pfam.resolve()),
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
        "strata_census": corpus.strata_census(),
        "domains_by_module_type":
            corpus.annotation_census["domains_by_module_type"],
        "curated_domain_totals": corpus.curated_domain_totals(),
        "protein_length_summary": corpus.length_summary(),
        "panel_census": {
            model: sum(1 for hs in hits.values() for h in hs if h.model == model)
            for model in rules.panel},
    }
    if args.level == "reader":
        return write(report)

    if args.level == "unresolved":
        # The tool ran and every hit is real; none is resolved away.
        table = corpus.architecture_table(hits, policy="keep_all")
        report["architecture_table"] = table
        report["architecture_census"] = corpus.architecture_census(table)
        report["overlap_census"] = corpus.overlap_census(hits)
        return write(report)

    table = corpus.architecture_table(hits)
    report["architecture_table"] = table
    report["architecture_census"] = corpus.architecture_census(table)
    report["overlap_census"] = corpus.overlap_census(hits)
    if args.level == "architect":
        return write(report)

    report.update(corpus.reconcile(hits))
    report.update(corpus.module_decomposition(hits))
    if args.level == "reconcile":
        return write(report)

    hits_by_cutoff = {c: corpus.hits(c) for c in rules.cutoffs}
    report["policy_sweep"] = corpus.policy_sweep(hits_by_cutoff)
    report["mapping_sweep"] = corpus.mapping_sweep(hits)
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
