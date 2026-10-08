"""Deterministic stub systems for the RiPP-precursor ladder.

Harness fixtures. Sandbox-only plus the pinned image, which every level above
`noncompute` needs: each scans the declared Pfam panel with the pinned HMMER.

  reader      clears R1: the selection census, the field's shapes and the tool
              identity
  localise    clears R2: where each declared core sits in its translation
  evidence    clears R3: the elected coordinate convention and the domain audit
  complete    clears R4
  firstshape  COMPETENCE probe, tooled: scans correctly and reads `core_sequence`
              as `core_sequence[0]`, which is right for the list shapes and takes
              the first RESIDUE where the field is a bare string. Clears R1 and
              fails R2.
  noncompute  LEAKAGE ablation: never runs the tool and never opens the inputs.

`firstshape` is the probe this campaign needs, because its symptom is a report
that looks complete rather than broken. Taking element zero of `core_sequence` is
the natural thing to write against the shape the field has most often, it raises
no error on any record, and it turns a third of the corpus into a one-residue
core that still localises -- just ambiguously, and to the wrong place. Nothing in
the output says a parse went wrong.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import pathlib

from npbench_c.templates.ripp_precursor import analyse as rpa
from npbench_c.templates.ripp_precursor.core import Corpus

LEVELS = ("reader", "localise", "evidence", "complete", "firstshape",
          "noncompute")
REPORT = "precursor_report.json"


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
            "mibig_release": "4.0", "biosynthetic_class": "ribosomal",
            "baseline_policy": "p1_case_folded", "threshold": "--cut_ga",
            "panel_size": 8,
            "hmmer_version": "HMMER 3.4", "pfam_version": "Pfam 35.0",
            "pfam_sha256": "0" * 64,
            "selection_census": {}, "core_sequence_shape_census": {},
            "core_sequence_case_census": {},
            "records_declaring_a_cleavage_location": 30,
            "ripp_type_census": {}, "status_census": {},
            "precursor_length_summary": {},
            "localisation_table": [], "localisation_census": {},
            "leader_length_summary": {},
            "convention_tally": {},
            "elected_convention": "to_as_0based_core_start",
            "cleavage_table": [], "domain_table": [], "panel_census": {},
            "records_with_a_panel_hit": 12,
            "boundary_relation_census": {}, "boundary_offset_summary": {},
            "policy_sweep": {}, "deletion_table": [], "deletion_census": {},
        })

    params = rpa.load_params(sandbox)
    corpus = Corpus.load(sandbox, params["annotations_input"],
                         params["reference_records_input"])
    rules = corpus.rules
    pfam = rules.resources["pfam_a_hmm"]
    from npbench_c.tools.cache import file_digest

    report = {
        "mibig_release": corpus.release,
        "biosynthetic_class": rules.biosynthetic_class,
        "baseline_policy": rules.baseline_policy,
        "threshold": rules.threshold,
        "panel_size": len(rules.panel),
        "hmmer_version": corpus.hmmer_version(),
        "pfam_version": pfam.version,
        "pfam_sha256": file_digest(pfam.resolve()),
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
    }
    if args.level == "reader":
        return write(report)

    if args.level == "firstshape":
        # Element zero of the field, whatever the field is. Correct for the list
        # shapes; one residue where the field is a bare string.
        broken = dataclasses.replace(
            corpus, records=tuple(
                dataclasses.replace(r, raw_core_sequence=r.raw_core_sequence[0])
                for r in corpus.records))
        table = broken.localisation_table()
        report["localisation_table"] = table
        report["localisation_census"] = broken.localisation_census(table)
        report["leader_length_summary"] = broken.leader_length_summary(table)
        return write(report)

    table = corpus.localisation_table()
    report["localisation_table"] = table
    report["localisation_census"] = corpus.localisation_census(table)
    report["leader_length_summary"] = corpus.leader_length_summary(table)
    if args.level == "localise":
        return write(report)

    before = corpus.domains(corpus.fasta())
    report.update(corpus.reconcile(table))
    report.update(corpus.domain_audit(table, before))
    if args.level == "evidence":
        return write(report)

    after = corpus.domains(corpus.fasta(table))
    report["policy_sweep"] = corpus.policy_sweep()
    report.update(corpus.deletion_outcome(table, before, after))
    return write(report)


if __name__ == "__main__":
    raise SystemExit(main())
