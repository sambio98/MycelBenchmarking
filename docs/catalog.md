# NPBench-C campaign catalog — 32 campaigns from 24 templates

Authoring unit is the **template** (one oracle, one contract, one ladder); a
campaign is a template instantiated on specific data. 24 templates, 8 of them
instantiated twice, gives 32 campaigns and 128 graded rungs.

The 8 doubled templates are the ones with a clean instantiation axis, and they
carry the **private split**, because that is where variant regeneration is
needed when the rotation trigger fires.

## Grounding status

Every row carries two flags, because "a campaign we can describe" and "a
campaign we can build" are different things and conflating them is how a
catalog fills up with ladders nobody can implement.

**Data** — is the source confirmed to contain the fields the ladder needs?
- `verified` — checked this session, field paths and counts below
- `likely` — source is open and standard, fields not yet inspected
- `probe` — needs a feasibility probe before authoring

**Build** — can the oracle run before the Phase 1 container exists?
- `now` — oracle is pure computation over downloadable data
- `image` — needs a pinned tool (antiSMASH, HMMER, Prodigal, BiG-SCAPE, matchms)

**12 of 24 templates are `build: now`.** Those should be authored first: they
need no container, so they cannot be blocked by Phase 1, and they exercise the
grading core on real data immediately.

## What was verified this session

MIBiG 4.0 (`dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz`, 954 KB,
3,013 entries, CC BY 4.0):

- **`loci[].evidence` is the BGC–product link evidence** and carries the
  vocabulary the G2 rule needs: Heterologous expression (638), Knock-out studies
  (542), Enzymatic assays (284), Gene expression correlated with compound
  production (219), Correlation of genomic and metabolomic data (82),
  Homology-based prediction (31), In vitro expression (25).
  **`compounds[].evidence` is a different axis** — how the *structure* was
  determined (NMR 1234, Mass spectrometry 550, MS/MS 237, X-ray 39, Total
  synthesis 34). The G2 admission rule reads `loci[].evidence`, not
  `compounds[].evidence`.
- **`quality` is not a usable filter.** 2,710 of 3,013 entries are
  `questionable`, which tracks `completeness: unknown` rather than scientific
  doubt; a `high` entry can carry no loci evidence at all. All 43 entries
  passing the strict mass-balance filter below are `questionable`, so filtering
  on quality as well would yield zero. The defensible filter is evidence.
  *Open: confirm with MIBiG curators what `quality` encodes before citing it.*
- `status` ∈ {active, retired, pending}; 213+ retired entries are diff material
  and must never be gold for a non-diff campaign.
- Classes: PKS 1238, NRPS 1020, ribosomal 447, other 421, terpene 241,
  saccharide 216 — enough spread for the coverage matrix.
- **NRPS A-domain substrates are annotated with name, SMILES and a
  `proteinogenic` flag** (2,879 of 4,034 modules carry substrate annotation).
- **PKS `at_domain.substrates` is empty throughout.** Extender-unit identity
  (malonyl vs methylmalonyl) is not in the data, so **PKS mass balance is not
  buildable from MIBiG** and the mass-balance template is NRPS-only.
- `modification_domains` is a closed 14-type vocabulary (ketoreductase 797,
  dehydratase 533, methyltransferase 203, epimerase 160, enoylreductase 116,
  thioesterase 41, oxidase 29, ligase 17, branching 6, product_template 4,
  hydroxylase 4, thioreductase 1, aminotransferase 1, other 23), so attributing
  a formula gap to tailoring is gradable as an enum.
- Strict mass-balance admission filter (active + compound formula + ≥2 modules +
  all modules NRPS + every module exactly one annotated substrate + strong loci
  evidence) yields **43 candidates**, 2 to 18 modules.

Also verified: UniProt REST serves version-pinned sequences (used by L6-1);
RDKit 2026.03.6 installs and computes molecular formula correctly — used **only**
for formula, never for canonical SMILES, with InChIKey remaining the comparison
key.

---

## L1 — Sequence & gene (4 templates, 5 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L1-1 | Gene calling under assembly fragmentation | probe | image |
| T-L1-2 ×2 | Domain architecture parsing | likely | image |
| T-L1-3 | Annotation transfer error detection | likely | now |
| T-L1-4 | Frameshift / pseudogene detection | probe | image |

**T-L1-1** Prodigal gene calls on a contig, then the same contig fragmented.
R1 emit a GFF-like call table honouring the coordinate contract · R2 calls match
gold within overlap-F1 tolerance · R3 claim: which genes are truncated at a
contig boundary, from a closed enum · R4 contig re-fragmented at new
breakpoints — predict which calls change. G1 throughout. R4 perturbation:
`truncate_contig`.

**T-L1-2** HMMER/Pfam domain architecture of a multidomain biosynthetic protein.
R1 schema-valid domain table (ordered, non-overlapping, coordinate convention
pinned) · R2 architecture matches gold by overlap-F1 · R3 claim: domain order
string and which domain is catalytically essential, closed enum · R4 a domain is
deleted from the sequence — predict the new architecture and the lost function.
G1. Instantiation axis: domain family (NRPS A-PCP-C vs PKS KS-AT-ACP).

**T-L1-3** Annotation transfer error: a protein whose closest BLAST/DIAMOND hit
is a paralog with a different function. R1 hit table · R2 hits match gold ·
R3 claim: is the transferred annotation supported, verdict enum + the
discriminating feature · R4 the paralog is removed from the reference set —
predict the new top hit and whether the verdict flips. `build: now` if the
reference set ships as a fixed FASTA subset rather than a live database.

**T-L1-4** Frameshift / pseudogene detection. R1 ORF table · R2 matches gold ·
R3 claim: which locus is a pseudogene and why, closed enum · R4 an indel is
introduced elsewhere — predict the new pseudogene call. Perturbation:
`shift_coordinates`.

## L2 — Enzyme & function (5 templates, 7 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L2-1 ×2 | EC misannotation triage (EC 1.1.3.15) | likely | image |
| T-L2-2 | Catalytic residue identification | likely | now |
| T-L2-3 | Remote homology twilight zone | probe | image |
| T-L2-4 ×2 | A-domain substrate specificity | verified | image |
| T-L2-5 | Kinetics harmonisation | likely | now |

**T-L2-1** Built on Rembeza & Engqvist 2021 (*PLoS Comput Biol* 17(9):e1009446):
78% of proteins annotated EC 1.1.3.15 in BRENDA 2017.1 lack the canonical FMN-dh
domain. R1 domain table over the sequence set · R2 matches gold · R3 claim: which
sequences are likely misannotated, plus the discriminating domain, closed enum ·
R4 a sequence is mutated to disrupt the FMN-dh domain — predict the
reclassification. Discrimination ECs: 1.13.12.4, 1.1.2.3, 1.1.99.31.
Instantiation axis: EC family. **Gold is the domain presence computed in-container
(G1), never the paper's percentage (which would be G3).**

**T-L2-2** Catalytic residues from PDB coordinates + UniProt annotation.
R1 residue table with numbering convention pinned · R2 matches gold · R3 claim:
which residue is the catalytic nucleophile and its role, closed enum · R4 the
residue is mutated — predict the geometric consequence (distance change to the
substrate analogue, interval-scored). PDB is CC0 and geometry is arithmetic, so
`build: now`.

**T-L2-3** Remote homology in the twilight zone: MMseqs2/DIAMOND over a set
seeded with true-but-distant homologues and plausible non-homologues.
R1 hit table · R2 matches gold · R3 claim: homologous or not per candidate,
with the evidence class · R4 the query region carrying the shared motif is
shuffled — predict which hits survive. Needs the thread-pinned MMseqs2 build
(issue #277: results differ by core count), so `build: image`.

**T-L2-4** NRPS A-domain substrate specificity. **Verified**: 2,879 modules carry
a substrate with SMILES and a `proteinogenic` flag. R1 per-module specificity
table · R2 matches gold by exact monomer identity · R3 claim: the specificity
code residues and the predicted monomer class · R4 specificity residues are
swapped to another annotated A-domain's — predict the new monomer. Instantiation
axis: monomer class (proteinogenic vs non-proteinogenic — the data has
4-hydroxyphenylglycine 80, 2,4-diaminobutyric acid 53, and similar).

**T-L2-5** Kinetics harmonisation across unit conventions. BRENDA is CC BY 4.0
but behind an acceptance gate, so it is Tier F (fetch-on-first-use, fail-closed
hash). R1 normalised kinetics table, units pinned · R2 matches gold · R3 claim:
which reported value is unit-inconsistent, plus the corrected value
(interval-scored) · R4 a value is rescaled by a factor — predict the corrected
table. Perturbation: `rescale_units`. Tag `license_class: open` but Tier F.

## L3 — Cluster & pathway (5 templates, 7 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L3-1 ×2 | BGC detection & boundary calling | likely | image |
| T-L3-2 | GCF clustering cutoff sensitivity | likely | image |
| T-L3-3 ×2 | MIBiG version diff reconciliation | verified | now |
| T-L3-4 | RiPP precursor annotation | likely | image |
| T-L3-5 | Self-resistance target identification | verified | now |

**T-L3-1** antiSMASH detection and boundary calling. R1 region table, coordinate
convention pinned · R2 boundaries match gold by overlap-F1 within measured
tolerance · R3 claim: cluster class and the boundary-defining gene, closed enum ·
R4 a core biosynthetic gene is deleted — predict the new boundary and whether
detection survives. **Pin antiSMASH version and ClusterBlast DB; never compare
across versions** (rule count went 58 → 88 across v5–v7.1). Instantiation axis:
BGC class (PKS, NRPS).

**T-L3-2** BiG-SCAPE GCF membership as a function of cutoff. R1 distance matrix ·
R2 matches gold · R3 claim: family membership at the declared cutoff · R4 the
cutoff is changed — report the range over which membership flips, scored with
`interval_score`. Absorbs the AN-5 robustness family. Pin BiG-SCAPE version,
`--mibig-version`, antiSMASH version, pyhmmer version, classification mode.

**T-L3-3** MIBiG 3.0 → 4.0 diff. **Verified**: `changelog.releases[]` and
`status` ∈ {active, retired, pending} are present, and 4.0 added 557 entries and
modified 590. R1 diff table (created / modified / retired) · R2 matches gold ·
R3 claim: for a given entry, what changed and in which field, closed enum ·
R4 a specific entry's prior version is supplied — predict the 4.0 delta.
Mechanically generable from the two tarballs, so this template can produce many
instantiations and is the natural private-split generator. Instantiation axis:
change class (retirement vs re-annotation).

**T-L3-4** RiPP precursor annotation: core vs leader peptide boundary.
R1 precursor table · R2 matches gold · R3 claim: the core peptide sequence and
the cleavage motif · R4 the leader is altered — predict the new core. 447
ribosomal-class entries available.

**T-L3-5** Self-resistance target identification: a BGC carrying a resistant
copy of the target the product inhibits. R1 gene table with the housekeeping
paralog identified · R2 matches gold · R3 claim: the resistance gene and the
inferred target, closed enum · R4 the resistance gene is removed — predict
whether the target call survives and what evidence is lost. Gold restricted to
entries with strong `loci[].evidence`.

## L4 — Molecule & structure (4 templates, 6 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L4-1 ×2 | MS² dereplication vs reference standards | likely | image |
| T-L4-2 | Molecular networking, mis-seeded annotation | likely | image |
| T-L4-3 ×2 | Mass balance & route validation (NRPS) | verified | now |
| T-L4-4 | Chemical space / scaffold analysis | likely | now |

**T-L4-1** Dereplication against authentic-standard libraries only — MassBank
2026.03 (CC BY 4.0, shippable) and the GNPS reference-standard sets (MSMLS 863
spectra, NIST14 high-confidence matches 5,763, MIADB 172, Bruker MetaboBASE
Plant ~1,300). CASMI 2022 (500 compounds, +/−ESI, solutions public) as the
externally calibrated difficulty anchor. R1 annotation table, InChIKey as the
key · R2 matches gold at the declared InChIKey block · R3 claim: confidence
level per annotation and the discriminating fragment · R4 an adduct shift is
applied — predict which annotations survive. Perturbation:
`adduct_misassignment`. **The nearest-neighbour suspect library is a distractor,
never gold** — it is propagated by construction.

**T-L4-2** Molecular networking with one planted library-backed-looking but
unconfirmed annotation. R1 network edge table · R2 matches gold · R3 claim:
which annotation is unsupported and why, closed enum from the confound
vocabulary · R4 the planted annotation is removed — predict the network change.

**T-L4-3** **Verified and ready to author.** NRPS mass balance: given the module
architecture and the product formula, validate the route.
R1 route table — ordered modules, each with monomer identity and molecular
formula, formulas parse · R2 naive assembly closes: Σ monomer formulas −
(n−1)·H₂O, atom-balanced per element, matched against the product formula ·
R3 claim: the residual between naive assembly and the actual product formula,
attributed to specific `modification_domains` (methyltransferase +CH₂,
ketoreductase/dehydratase H/O changes) and/or a starter unit, with the signed
per-element discrepancy · R4 a module is deleted — predict the new product
formula and the exact mass shift. All G1. 43 admissible entries; instantiation
axis: module count (BGC0001873, 2 modules, C13H20N2O2 is the smallest; orfamide A
BGC0000399, 10 modules, C64H114N10O17 exercises a fatty-acid starter).
Monomer formulas come from a hash-pinned table with ChEBI provenance, computed
from the SMILES MIBiG supplies.

**T-L4-4** Scaffold / chemical-space analysis keyed on InChIKey block 1
(connectivity), which is the right granularity where stereochemistry is not
determinable. R1 scaffold table · R2 matches gold · R3 claim: which compounds
share a scaffold and the distinguishing substitution · R4 a substituent is
changed — predict the new scaffold assignment.

## L5 — Target & mechanism (3 templates, 3 campaigns)

Gold here is **measured-or-computed only**: a value from a database with a
structured assay-confidence field, or exact arithmetic over such values. Never a
model prediction, never a curated expert call. This is the layer where integrity
breaks if the rule is relaxed.

| # | Template | Data | Build |
|---|---|---|---|
| T-L5-1 | Selectivity ratios from measured IC50 | likely | now |
| T-L5-2 | Resistance mutation → target assignment | likely | now |
| T-L5-3 | Binding-site geometry from coordinates | likely | now |

**T-L5-1** ChEMBL measured activities (CC BY-SA 3.0 — **share-alike propagates
to redistributed subsets; confirm licensing before this is load-bearing**).
R1 activity table, assay-confidence field retained · R2 matches gold ·
R3 claim: the selectivity ratio between two targets (exact arithmetic over
measured values) and whether it exceeds the declared threshold · R4 a
low-confidence assay is excluded — predict the new ratio.

**T-L5-2** Resistance mutation to target assignment from structured fields.
R1 mutation table · R2 matches gold · R3 claim: the implicated target and the
evidence class · R4 a mutation is removed from the set — predict whether the
assignment survives.

**T-L5-3** Binding-site geometry from PDB coordinates (CC0, pure arithmetic).
R1 residue-contact table, distance cutoff pinned · R2 matches gold ·
R3 claim: which residues line the pocket and the pocket volume (interval-scored)
· R4 a side chain is substituted in silico — predict the contact-set change.

## L6 — Translation & design (3 templates, 4 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L6-1 ×2 | Construct & sequence design | verified | **BUILT** |
| T-L6-2 | Multi-constraint conflict resolution | verified | now |
| T-L6-3 | Prioritisation under budget | likely | image |

**T-L6-1** **Built**: `campaigns/construct-ecoli-rebh-01`. Oracle 1.0 on all
four rungs, measured discrimination 1/2/3/3/4, binding constraint flips between
hosts. Second instantiation: PrnA (UniProt P95480, 538 aa, reviewed) — the
instantiation axis is the enzyme, with the host swap staying as R4.

**T-L6-2** A variant where the declared constraint set is **unsatisfiable** and
the correct answer is to identify the conflicting pair and abstain on the
design. Abstention gold is **constructed, not curated**: the conflict is created
by tightening two constraints past feasibility, so `insufficient` is true by
construction and regenerable. R4: one constraint is relaxed — predict whether
the design becomes feasible and which constraint now binds.

**T-L6-3** Candidate prioritisation under a hard container query cap.
R1 ranked table · R2 ranking matches gold by Spearman within tolerance ·
R3 claim: the top candidate and the discriminating criterion · R4 the budget is
halved — predict the new top-k. Budget is a hard container limit, not a counted
metric, so a system with coarser tool granularity is not differently penalised.

---

## Coverage check

**NP class** (each ≥2 campaigns): PKS (L3-1, L4-4) · NRPS (L2-4 ×2, L4-3 ×2,
L3-1) · RiPP (L3-4, L4-2) · terpene (L4-4, L5-1) · saccharide (L4-4, L3-2) ·
enzyme-centric non-BGC (L2-1 ×2, L2-2, L2-3, L2-5, L6-1 ×2, L6-2).
*Gap: terpene and saccharide coverage is thin and leans on shared templates —
worth a dedicated template if a terpene cyclase product-prediction task can be
grounded.*

**Task archetype** (each ≥2): detect (L1-1, L3-1, L3-4) · quantify (L2-5, L4-3,
L5-1, L5-3) · reconcile (L1-3, L3-3, L2-1) · adjudicate (L4-2, L2-3, L5-2) ·
counterfactual (every R4 by construction) · design (L6-1, L6-2) · abstain
(L6-2, L2-3).

A CI check should compute both matrices from `task.yaml` files and fail the
build on an empty row or column, so generality is enforced rather than claimed.

## Authoring order

1. **`build: now` + `data: verified`** — T-L4-3 (mass balance), T-L3-3 (MIBiG
   diff), T-L6-1 second instantiation, T-L6-2, T-L3-5. Five campaigns needing no
   container.
2. **`build: now` + `data: likely`** — T-L1-3, T-L2-2, T-L2-5, T-L4-4, T-L5-1,
   T-L5-2, T-L5-3. Probe each source first.
3. **`build: image`** — everything else, blocked on Phase 1 and on the
   thread-invariance suite for MMseqs2, DIAMOND, HMMER, antiSMASH.

Author ~44 candidates to ship 32: the internal sweep will drop some for
non-discrimination and the audit will drop some for gold problems.
