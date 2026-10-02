# Re-grounding pass — 2026-10-02

Every remaining `build: now` row in the catalog, probed against its actual data
source before any oracle work. Four earlier catalogued designs had already been
corrected this way; this pass covers the rest so the size of the S0 set is known
rather than assumed.

## Verdicts

| Template | Catalog said | Re-grounded | Decided by |
|---|---|---|---|
| T-L1-3 annotation transfer | `now` / likely | **`image`** | finding the closest hit needs alignment |
| T-L2-2 catalytic residues | `now` / likely | **`now` / verified** | UniProt features + PDB coordinates both parse |
| T-L2-5 kinetics harmonisation | `now` / likely | **blocked — credentials** | BRENDA data file 404s |
| T-L4-4 chemical space | `now` / likely | **`now` / verified** | 3,449 SMILES, 100% parseable |
| T-L5-1 selectivity ratios | `now` / likely | **verified, licence-blocked** | ChEMBL fine; share-alike unresolved |
| T-L5-2 mutation → target | `now` / likely | **verified, licence-blocked** | structured variant fields exist |
| T-L5-3 pocket geometry | `now` / likely | **`now` / verified** | mmCIF coordinates + bound ligand |

**The S0 buildable set is 3 templates, not 7.** Two more are data-ready and wait
only on a licensing decision; one needs credentials; one belongs in the
container tier.

## What each probe found

### T-L2-2 catalytic residues — verified, and cheaper than expected

UniProt serves per-residue features with **ECO evidence codes**, so the gold can
be evidence-gated by the same doctrine as the MIBiG G2 rule rather than by
curation:

```
Q8KHZ8  Active site 1, Binding site 17   ECO:0000269 (PubMed) x32, ECO:0007744 (PDB) x31
P95480  Active site 1, Binding site 24   PDB xrefs: 2APG 2AQJ 2AR8 2ARD ...
```

The Active site feature carries an exact position and a PubMed-backed ECO code
(`{"type": "Active site", "location": {"start": {"value": 79}}, "evidences":
[{"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "17260957"}]}`).

Two conveniences worth exploiting: **RebH and PrnA are already in the
benchmark**, and both carry these features plus 7–10 PDB cross-references. So
this template needs no new input provenance. Instantiation pool across UniProt:
**39,232** reviewed entries with `existence:1` and a PDB cross-reference.

### T-L5-3 pocket geometry — verified

`files.rcsb.org` serves mmCIF directly (2E4G, 970 KB, 13,816 lines, 9,338
ATOM+HETATM records, parseable with stdlib string handling). Hetero residues:
833 HOH and **30 TRP** — a bound tryptophan, which is RebH's substrate, so
contact-set and pocket-volume computation has a real ligand to work against.

Pure coordinate arithmetic, so **S0**, no chemistry toolkit needed.

### T-L2-2 and T-L5-3 should share one template

Both read UniProt features plus PDB coordinates for the same proteins. That is
the `construct_design` situation again — one engine, two ladders — and building
them separately would duplicate the mmCIF parser and the feature extractor.
Recommend a `structure_features` template with `residues/` and `pocket/`
ladders.

### T-L4-4 chemical space — verified and rich

Over MIBiG 4.0 active entries:

```
compounds carrying SMILES           3,449
RDKit-parseable                     3,449  (100.0%)
of those, strong loci evidence      1,238
distinct Murcko scaffolds           1,655
distinct InChIKeys                  3,115
```

3,449 compounds against 3,115 InChIKeys means **334 cross-entry duplicates** —
the same compound annotated under more than one BGC. That is campaign content in
its own right, not noise.

RDKit is needed here (scaffold perception), so this is **S1** unless scaffolds
are precomputed into `reference/` the way monomer formulas were for mass
balance — which would make it S0 and is the recommended route.

### T-L5-1 and T-L5-2 — the data is fine, the licence is not

ChEMBL's REST API works and has the structure both templates need:

- 5,265 IC50 activities for a single target, with `standard_value`,
  `standard_units`, `pchembl_value`, `data_validity_comment`,
  `potential_duplicate`
- `confidence_score` lives on the **assay**, not the activity
  (`CHEMBL652769` → 8, "Homologous single protein target assigned"), so an
  activity → assay join is required. Worth stating in the output contract.
- activities carry **`assay_variant_accession`** and
  **`assay_variant_mutation`** — structured mutation fields, which is what makes
  T-L5-2's "which mutation confers resistance" computable from measured IC50
  ratios rather than from prose.

**Blocker: ChEMBL is CC BY-SA 3.0 and share-alike propagates to redistributed
subsets.** Both templates are otherwise ready. This is a decision for whoever
owns licensing, not a data problem, and it gates two of the seven rows.

### T-L2-5 kinetics harmonisation — blocked on credentials

`brenda-enzymes.org/download.php` returns 200, but that is the HTML page; the
actual archive (`/download/brenda_download.tar.gz`) **404s**. The data is behind
registration and an acceptance gate, as §4 anticipated. So this is Tier F *with
credentials*, not `build: now`, and it cannot be authored until someone
registers and the licence acceptance is recorded.

### T-L1-3 annotation transfer — re-tiered to `image`

The task needs the closest homolog of a query to be *found*, which needs
alignment. Shipping a precomputed hit table would ship the answer, and
reimplementing alignment in the oracle would be reinventing a pinned tool. So it
belongs in the container tier alongside the other DIAMOND/MMseqs2 rows.

There **is** an S0-shaped neighbour: UniProt annotates function with ECO codes,
so "which of these annotations were transferred by similarity rather than
measured" is decidable from structured fields — the same move that turned
catalogued T-L3-5 into the gene-function evidence audit. It is a different task
from annotation-transfer *error* detection and should be catalogued separately if
wanted, not substituted silently.

## Consequences for the plan

1. **The S0 set is 3 templates.** With 5 templates built, the near-term ceiling
   without the container or a licensing decision is about 8 templates, not the
   catalog's implied 12.
2. **Phase 1 matters sooner than the catalog implies.** Nine of the remaining
   rows need pinned tools, and T-L1-3 now joins them.
3. **One decision unblocks two templates**: ChEMBL share-alike.
4. **One registration unblocks one**: BRENDA.
5. Build order for the three that are ready: `structure_features` (two ladders,
   shared engine, reuses RebH/PrnA provenance) then chemical space.
