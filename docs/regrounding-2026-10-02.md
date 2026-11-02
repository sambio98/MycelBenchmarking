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
| T-L2-5 kinetics harmonisation | `now` / likely | **wrong on both counts** | see the correction below |
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
distinct InChIKeys                  3,115
distinct InChIKey block 1           2,996
distinct Murcko scaffolds           1,843   (by scaffold InChIKey)
compounds with no ring system         182
```

**Two numbers in this section were wrong when first written and are corrected
above.** The scaffold count was recorded as 1,655, and it does not reproduce
under any definition: 1,843 keyed on the scaffold's InChIKey, 1,850 keyed on its
canonical SMILES over active entries, 2,254 over all entries including retired
ones. The build uses the InChIKey, because canonical SMILES is a toolkit output
and cannot be a key in a benchmark whose scores must stay comparable. And
3,449 records against 3,115 InChIKeys is **334 duplicate records**, not 334
cross-entry duplicates: most of those records repeat a structure inside one
entry. The cross-entry figure — keys appearing under more than one accession — is
**217** at the full key and **259** at the connectivity block. Both are campaign
content; they are not the same statistic and conflating them was the error.

What the full build then found, which no count in this table shows: **Bemis-Murcko
reduction is much coarser than the catalogued design assumed.** The most widely
shared scaffold in this corpus is plain benzene, in 76 entries. See
`docs/catalog.md` under T-L4-4 for the two scope corrections that follow.

RDKit is needed here (scaffold perception), so this is **S1** unless scaffolds
are precomputed into `reference/` the way monomer formulas were for mass
balance — which would make it S0 and is the recommended route. **Taken**: the
built campaign pins the chemistry once at instantiation and is S0.

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

### T-L2-5 kinetics harmonisation — this verdict was wrong twice over

What this section said: the archive 404s, the data is behind registration, the row
cannot be authored until someone registers. **Both halves were wrong**, and the
error was guessing a URL instead of reading the page.

`/download/brenda_download.tar.gz` does 404 — because it is not the download
path. `download.php` serves a form: a licence-acceptance checkbox enables three
buttons, and a POST with `dlfile=dl-json&accept-license=1` returns
`brenda_2026_1.json.tar.gz`, 83 MB, no account and no API key. The gate is a
**licence acceptance**, not a registration wall, and the licence is plain CC BY
4.0 — the page says so in as many words.

The second error was in the task, not the access. "Harmonisation across unit
conventions" is not what the data supports: BRENDA has already normalised units
per field and records **no unit on any value**, so there are no competing
conventions in the file to harmonise. What it does support is the consequence of
the convention being implicit — Km, kcat and kcat/Km are stored independently and
arithmetically linked, so the corpus can be checked against itself, and the
disagreements cluster at exactly the factors a unit mistake would produce. Built
as `campaigns/kinetics-brenda-ec1_1-01`; see `docs/catalog.md` under T-L2-5.

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
4. **One registration unblocks one**: BRENDA. *(Not a registration at all — see
the correction above.)*
5. Build order for the three that are ready: `structure_features` (two ladders,
   shared engine, reuses RebH/PrnA provenance) then chemical space.
   **All three are now built** — `residues-prna-01`, `pocket-rebh-01`,
   `chemspace-mibig-4_0-01` — so the S0 set is exhausted and the three blockers
   above (Phase 1, ChEMBL share-alike, BRENDA credentials) are now the only
   things between the benchmark and its remaining 19 templates.
