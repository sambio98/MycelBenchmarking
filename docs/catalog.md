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
| T-L1-3 | Annotation transfer error detection | verified | **image** (re-tiered) |
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

## L2 — Enzyme & function (5 templates, 7 catalogued campaigns, **4 reachable**)

T-L2-4 is retired on grounding and T-L2-1's second instantiation is unresolved,
so this layer now offers 4 campaigns rather than the catalogued 7. Both
reductions are measured below.

| # | Template | Data | Build |
|---|---|---|---|
| T-L2-1 | EC misannotation triage (EC 1.1.3.15) | verified | **BUILT** (first S1) |
| T-L2-2 | Catalytic residue identification | verified | **BUILT** |
| T-L2-3 | Remote homology twilight zone | probe | image |
| T-L2-4 ×2 | A-domain substrate specificity | **thin / vocabulary-bound** | re-tiered, see below |
| T-L2-5 | Kinetics consistency (re-scoped) | verified | **BUILT** |

**T-L2-1** Built on Rembeza & Engqvist 2021 (*PLoS Comput Biol* 17(9):e1009446):
78% of proteins annotated EC 1.1.3.15 in BRENDA 2017.1 lack the canonical FMN-dh
domain. R1 domain table over the sequence set · R2 matches gold · R3 claim: which
sequences are likely misannotated, plus the discriminating domain, closed enum ·
R4 a sequence is mutated to disrupt the FMN-dh domain — predict the
reclassification. Discrimination ECs: 1.13.12.4, 1.1.2.3, 1.1.99.31.
Instantiation axis: EC family. **Gold is the domain presence computed in-container
(G1), never the paper's percentage (which would be G3).**

**Built**: `campaigns/ecaudit-fmn-dh-1_1_3_15-01`, the benchmark's **first S1
campaign** — the first that executes a tool, and the first that reads a pinned
resource it is not allowed to ship. UniProt 2026_03, every entry annotated EC
1.1.3.15: **7,675 sequences, 29 reviewed and 7,646 unreviewed**, both sections in
full because the split is the control and sampling it would put a selection
decision inside the input. HMMER 3.4 against Pfam 35.0 at `--cut_ga`, over a
declared 13-model panel.

Measured:

| | FMN_dh present | absent | rate |
|---|---|---|---|
| reviewed | 29 / 29 | **0** | 0.0000 |
| unreviewed | 6,300 / 7,646 | **1,346** | **0.1760** |

**The paper's 78% does not reproduce on UniProt 2026_03, and that is the point of
the G1 rule.** Rembeza & Engqvist measured BRENDA 2017.1; this campaign measures
what the pinned tool says about the current release and gets 17.6%. Quoting the
published figure as gold would have been a literature-recall test *and* wrong.

Two findings make the ladder more than a tool invocation.

**The canonical-absent sequences are not fragments.** That was the explanation
that would have made the whole audit an artefact of sequencing completeness: a
sequence can miss a 348-column model by being short. Measured, the absent side is
**longer** — median 466 residues against 367 — and only 2.4% of it is under 300
residues against 12.6% of the present side. The length profile is a graded
component so an auditor sees the refutation rather than taking it on trust.

**They are a coherent enzyme family, not noise.** 967 of the 1,346 (71.8%) carry
exactly `FAD-oxidase_C + FAD_binding_4` — the FAD-linked glycolate oxidase subunit
architecture, which oxidises the same substrate through a different cofactor. So
absence of the canonical domain is **not** the same as misannotation, and the
campaign refuses to assert a misannotation rate. It grades the count under three
published policies instead: **1,346** accepting only the canonical FMN domain,
**319** admitting the FAD-oxidase architecture, **22** admitting any flavin
architecture in the panel. The p0/p1 gap is the honest statement of where the
uncertainty lives, and it is R4's load-bearing axis.

**Deviation: the catalogued R4 cannot be built, and the discrimination ECs are not
used.** Mutating a sequence to disrupt an HMM match means choosing which residues
to break, and any choice big enough to drop a curated gathering threshold is a
choice about the answer — with gold then describing a sequence nobody deposited.
R4 varies the declared cutoff and the declared policy instead. The cutoff turns
out to be nearly inert (1,342 to 1,371 absent over six cutoffs, with `--cut_ga`,
`--cut_nc` and `--cut_tc` identical), which is a robustness result rather than a
defect, and is admissible only because the policy axis in the same rung moves by a
factor of sixty. As for the catalogued discrimination ECs: 1.13.12.4 has 204
entries and 1.1.99.31 has 234, and **EC 1.1.2.3 is itself an FMN_dh family
member** — all three reviewed L-lactate dehydrogenase (cytochrome) entries hit the
model at the gathering threshold — so the canonical domain does not discriminate
it from EC 1.1.3.15 at all. Measured, not assumed.

**The ×2 instantiation is now an open question, not a plan.** The catalogued axis
was "EC family", and this build used one. A second instantiation needs an EC with
the same shape — a canonical domain the annotation implies, a reviewed set big
enough to be a control, and an unreviewed set large enough for a rate — and the
catalogued discrimination ECs do not supply it: 1.1.2.3 shares the canonical
domain, and 1.13.12.4 and 1.1.99.31 have 204 and 234 entries. Finding one is a
grounding pass of its own, so this row counts as **one** campaign until somebody
does it.

**The panel is declared, not discovered.** It was fixed once at build time by
scanning the canonical-absent set against the whole of Pfam 35.0 and listing every
family that fired; the campaign then publishes it, and an agent scans with exactly
those 13 models in about 11 seconds instead of the 2.5 minutes a full Pfam-A pass
costs per run. It also **partitions this set perfectly**: all 6,329
canonical-present sequences carry FMN_dh and nothing else from the panel, so the
architecture census is a genuine partition rather than a ranking of overlapping
families.

**And the library is read, not shipped.** Pfam-A.hmm is 1.5 GB, so it stays in the
image and the campaign declares the path, the provider, the version and its
sha256; the build refuses to emit gold if the file it finds does not match. Here
the reason is size and Pfam is CC0 — recorded explicitly, because for the next S2
campaign the reason will be licence. This is what the gate's new twentieth check,
`external_resources_pinned_and_present`, verifies.

**T-L2-2** Catalytic residues from PDB coordinates + UniProt annotation.
**Built as the `residues` ladder of the `structure_features` template**:
`campaigns/residues-prna-01`. R1 the annotation census under a declared dedup
rule · R2 the ECO evidence classes and the PDB citation index · R3 each
annotation checked against the geometry of the structure it cites · R4 the same
check on a second entry with a different substrate bound, plus four
evidence-policy variants. PDB is CC0 and geometry is arithmetic, so `build: now`.

**Deviation from the catalogued design, and why.** The catalogued R3 was "which
residue is the catalytic nucleophile, closed enum" and the catalogued R4 was a
predicted distance change on mutation. Both were re-scoped after the data was
read. The nucleophile is a single structured UniProt field (`Active site`), so
asking for it is a lookup, not a claim — it survives as one component of R3
rather than as the rung. And a predicted distance change on mutation is not
computable from deposited coordinates without rebuilding a side chain, which is
modelling, not arithmetic; what IS computable is the contact consequence of
deleting a side chain, which the pocket ladder grades. The claim R3 makes
instead — *does the structure UniProt cites actually support the annotation* —
is the same evidence doctrine the MIBiG audit applies, turned on a second
database, and it is decided entirely by arithmetic.

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

**Re-tiered. Every rung this row can actually carry is either near-constant or
decided by a naming artefact, and the measurements say so.** Grounded against
MIBiG 4.0 and the Stachelhaus signature table antiSMASH 8.0.4 ships
(`databases/nrps_pks/stachelhaus/1.1/signatures.tsv`, 2,319 rows).

What is there: 2,901 A-domains in 583 entries, 1,228 distinct gene accessions,
2,879 with a substrate. What is not:

- **The A-domain's position in its protein is absent.** `location` is `-1/-1` for
  **2,852 of 2,901** (98.3%). The catalogued R3 wants "the specificity code
  residues", which needs the domain located in a sequence; MIBiG does not carry it
  and MIBiG carries no sequences either (no `translation` anywhere in the JSON).
- **The substrate labels are mostly not experimental.** 530 A-domains have
  experimental support (18.3%), 2,076 are inference-only (71.6%) and 295 carry no
  evidence at all (10.2%); 1,049 name "Sequence-based prediction". So "R2 matches
  gold by exact monomer identity" would grade an agent against labels that are
  three-quarters inference.

The signature table joins to MIBiG on protein accession — 714 shared accessions —
but only **304 genes carry exactly one A-domain on both sides**. 180 genes have
equal counts above one (473 A-domains) and are mappable only under the assumption
that MIBiG's module order equals the table's `.A1/.A2` index, which nothing in
either resource confirms; 230 have mismatched counts; 514 are MIBiG-only and 923
signature-table-only. On the 304 unambiguous pairs, four comparison keys give four
different answers:

| comparison key | decidable | agreement | why |
|---|---|---|---|
| raw substrate names | 304 | **0.000** | `Ile` vs `isoleucine`: different vocabularies |
| the `aaSMILES` comment field as a name table | 193 | 0.912 | 17 of its 18 disagreements are artefacts |
| stereo-aware InChIKey | 219 | **0.151** | the resources differ on annotating stereochemistry |
| connectivity InChIKey (block 1) | 219 | **0.995** | 198 agree, 20 within-multi, **1 real disagreement** |

The middle row is the sharpest warning. antiSMASH's `aaSMILES.txt` has a trailing
`#` comment that usually holds a substrate name, and used as a name table it
produces 18 disagreements of which **one** is substantive: the rest are a hyphen
(`4-hydroxy-phenylglycine` against MIBiG's `4-hydroxyphenylglycine`) or a note to
the maintainer (`salicylic acid (not in norine yet)`, `actually aile,
allo-isoleucine`). The comment field is a comment, not a vocabulary.

Compared properly — connectivity-level InChIKey, the key this benchmark already
pinned for the chemical-space campaign — the two resources agree on **218 of 219**
comparable A-domains. The single disagreement in the entire set is ABA59548.1,
where the table says `Val` and MIBiG says pyruvic acid. **84 of the 304 cannot be
compared at all**, because `aaSMILES` names 94 substrates and the signature table
names 280.

The catalogued ×2 axis also fails: on the proteinogenic side agreement is
**1.000** with zero disagreements, and on the non-proteinogenic side **83 of 112**
are unmappable. And the within-MIBiG route — group A-domains by identical
signature and ask whether MIBiG's own labels agree, which needs no bridge at all —
collapses too: the 34-residue extended signature is **unique for all 304**, so it
groups nothing, and the 10-residue signature yields 33 groups of two or more, 10
of them inconsistent, where the inconsistencies are again granularity
(`threonine` against `allo-threonine`, `glutamic acid` against `d-glutamic acid
branched`).

So there is no rung here whose answer varies for a reason an agent can reason
about: either everything agrees, or everything disagrees for a lexical reason.
**A campaign built on it would be a vocabulary quiz with a hidden key**, which is
the thing the no-human-judgement rule exists to keep out. The numbers are recorded
so the next person does not re-probe.

**A licensing rule falls out of this, and it binds every future S2 campaign.**
`aaSMILES.txt` is antiSMASH package data under **AGPL-3.0-or-later**, and
`stachelhaus/1.1/signatures.tsv` comes from the antiSMASH database download with
**no stated licence at all**. Neither can be shipped inside a campaign: you cannot
redistribute data whose licence you cannot name, and pulling AGPL material into
`inputs/` is not a decision to make by accident. An S2 campaign reads such a
resource **from the pinned image at a declared path with a recorded sha256**, and
redistributes nothing.

**T-L2-5** **Built**: `campaigns/kinetics-brenda-ec1_1-01`. BRENDA 2026.1,
EC subclass 1.1.- , 437 EC numbers and 35,188 kinetic records.
R1 the per-field parse census, with BRENDA's -999 missing-value sentinel counted
and excluded rather than averaged · R2 the triple inventory under the declared
(EC, proteins, substrate, references) key · R3 the consistency verdict: compute
kcat/Km from kcat and Km and band the ratio against the stored value · R4 the same
check under three declared unit mistakes, with an identity control.

*Two scope corrections.* **The access was not blocked.** This catalog recorded
T-L2-5 as blocked on BRENDA credentials; the download is a licence-acceptance
checkbox, not a registration wall, and the licence is plain CC BY 4.0. The earlier
verdict came from guessing a URL instead of reading the page.

**And "harmonisation across unit conventions" is not what the data supports.**
BRENDA has already normalised units per field and records no unit on any value, so
there are no competing conventions in the file to harmonise. What it does support
is the consequence of the convention being implicit: the three constants are
arithmetically linked, so the corpus checks itself, and a system that does not
know Km is in mM and kcat in s⁻¹ cannot do the check at all. Measured over 2,528
complete triples: **63.2% agree within 5%**, and the disagreements include **31 at
a factor of ~1000** (millimolar read as micromolar) and **25 at ~60** (per-second
read as per-minute). Those are the catalogued unit conventions, found in the data
rather than postulated. The R4 counterfactual makes the dependence explicit —
applying the micromolar mistake deliberately drops the agreement rate from 63.2%
to 0.36%.

## L3 — Cluster & pathway (5 templates, 7 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L3-1 ×2 | BGC detection & boundary calling | likely | image |
| T-L3-2 | GCF clustering cutoff sensitivity | likely | image |
| T-L3-3 ×2 | MIBiG version diff reconciliation | verified | **BUILT** |
| T-L3-4 | RiPP precursor annotation | likely | image |
| T-L3-5 | Self-resistance target identification | verified | **image** (re-tiered) |
| T-L3-5a | Gene-function evidence audit | verified | **BUILT** |

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

**T-L3-3** **Built**: `campaigns/mibig-diff-3_1-to-4_0-01` (3.1 → 4.0, the
consecutive-release pair; 3.1 is the last 3.x).

*Scope correction made during the build.* "What changed and in which field" is
**ill-posed**: 4.0 restructured the schema entirely — 3.1 nests everything under
a `cluster` object, 4.0 is flat with twelve top-level keys — so a raw
field-by-field diff reports almost every field of almost every entry as changed.
True and useless.

What is well-posed, and better: **reconcile the two releases on a declared field
mapping, separating a change in representation from a change in substance.** The
campaign publishes the mapping, the per-field normalisation rules, the verdict
vocabulary and the precedence over it.

The sharpest case is supplied by the corpus rather than planted: 3.1 stores the
NCBI taxon id as a JSON **string**, 4.0 as an **integer**. Compared raw, all
2442 shared entries look re-annotated; compared as strings, **8** genuine
corrections remain. R2 asks for the raw counts and R3 for the classified ones,
so a system that never normalises clears R2 and fails R3 — which the sweep
demonstrates directly, since `stub-nonormalise` lands at depth 2.

Measured gold, 2502 → 3013 entries, 2442 shared:

| field | raw diff | classified |
|---|---|---|
| ncbi_tax_id | 2442 | 2434 representation-only, **8** substantive |
| biosynthetic_class | 2442 | 2353 representation-only, 73 retired-term, **16** substantive |
| compound_formulas | 181 | 1380 unchanged, 1015 absent-in-source, 39 substantive |
| compound_names | 152 | 2289 unchanged, 152 substantive |
| organism | 8 | 8 substantive |
| locus_accession | 0 | 2442 unchanged |

3.1's `Alkaloid` class has no 4.0 counterpart and 73 of its entries persist
reclassified, which is why `no_counterpart_in_target` precedes
`substantive_change` in the published precedence. The 8 taxon-id and 8
organism-name corrections are nearly disjoint (overlap: BGC0002347), so they are
independent curation events.

R4 applies five **declared operations** to a base record — taxon id retyped,
taxon id revalued, a formula altered, the class set replaced with the retired
term, the organism removed — and asks the resulting verdict. The first is the
campaign in miniature: retyping moves the verdict from `representation_only` to
`unchanged`, which a system that has not modelled the distinction gets wrong.

The build **refuses to emit gold** if no field's raw count differs from its
classified substantive count, since R3 would then merely restate R2.

Still the natural private-split generator: the field mapping and the
perturbation operations generalise to any release pair.

**T-L3-4** RiPP precursor annotation: core vs leader peptide boundary.
R1 precursor table · R2 matches gold · R3 claim: the core peptide sequence and
the cleavage motif · R4 the leader is altered — predict the new core. 447
ribosomal-class entries available.

**T-L3-5** **Re-tiered to `build: image`.** MIBiG 4.0 *does* carry a structured
`Resistance/immunity` gene-function category (43 annotations, 20 entries with
strong `loci[].evidence`), which makes "identify the resistance gene" a field
lookup rather than an inference. The **target** the product inhibits is not in
MIBiG at all: establishing it needs either a free-text product string
("resistant DNA gyrase subunit B" — transcription, not reasoning) or homology
search to show the gene is a paralog of an essential enzyme, which needs
HMMER/BLAST. So the inferential version belongs in the container tier, not at
S0.

**T-L3-5a** **Built**: `campaigns/mibig-gene-function-evidence-01` — the
S0-feasible audit underneath it, which is the benchmark's own G2 gold-admission
doctrine turned on the database that supplies much of its gold.

MIBiG records **three independent evidence axes**, and conflating them is the
error the campaign targets: `loci[].evidence` (how the cluster was linked to its
product), `compounds[].evidence` (how the structure was determined),
`functions[].evidence` (how a gene's function was established). An annotation is
admissible only under all three conditions — accepted function evidence,
accepted locus evidence, active status.

Measured over 3013 entries, 272 with gene annotations, 2566 (annotation,
function) pairs:

| function | backed | unbacked | admissible | gap |
|---|---|---|---|---|
| Scaffold biosynthesis | 355 | 636 | 225 | 130 |
| Tailoring | 303 | 334 | 204 | 99 |
| Precursor biosynthesis | 134 | 269 | 84 | 50 |
| Regulation | 82 | 151 | 63 | 19 |
| Activation / processing | 32 | 51 | 20 | 12 |
| **Transport** | **23** | **153** | 19 | 4 |
| **Resistance/immunity** (focal) | **21** | 22 | **5** | 16 |

950 backed collapses to **620 admissible**. R2 asks for the backed census and R3
for the admissible one, so `stub-conflated` — which reports merely-backed as
admissible — clears R2 and fails R3. R4 changes the evidence policy four ways
and every variant moves the total (595 / 749 / 636 / 640 against 620); the build
**refuses to emit gold** if they do not.

## L4 — Molecule & structure (4 templates, 6 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L4-1 ×2 | MS² dereplication vs reference standards | likely | image |
| T-L4-2 | Molecular networking, mis-seeded annotation | likely | image |
| T-L4-3 ×2 | Mass balance & route validation (NRPS) | verified | **BUILT** |
| T-L4-4 | Chemical space / scaffold analysis | verified | **BUILT** |

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

**T-L4-3** **Built**: `campaigns/massbalance-nrps-malleobactin-01` (BGC0000386,
4 modules, 4 congeners including a dimer). Oracle 1.0 on all four rungs, stub
discrimination 1/2/3/4, all five sweep gates green.

*Scope correction made during the build.* "Validate the route" is **not**
buildable: only 1 of 43 admissible entries is exactly balanced and 4 differ by
one water, because MIBiG annotates the enzymatic module count rather than the
number of monomer incorporations — iterative reuse, fatty-acid starters, prenyl
groups and metal centres all break the correspondence. The campaign instead
computes what the annotated modules *would* produce and quantifies the gap,
which is well-posed; attributing the gap to named chemistry would be judgement
and is excluded with its reason recorded in `task.yaml`.
R1 ordered module route, each with monomer identity and a parsing formula ·
R2 monomer formulas correct, naive assembly Σ monomers − (n−1)·H₂O correct per
element, peptide-bond count and assembly monoisotopic mass correct · R3 the
signed per-element residual **for every annotated congener** plus a mechanical
verdict from a closed vocabulary · R4 for **each** single-module deletion, the
resulting assembly formula and exact mass shift. All G1.

Monomer formulas are resolved once at construction with RDKit and pinned into
`reference/monomer_formulas.json` with their SMILES source, so the oracle needs
no chemistry toolkit at run time and the campaign is solvable at **S0**.
Second instantiation **built**: `campaigns/massbalance-nrps-sevadicin-01`
(BGC0000426, 3 modules, 1 compound, residual exactly zero). It exercises the
`balanced` verdict branch, which the malleobactin instantiation never reaches —
all four of its congeners are `additional_chemistry_required`.

Both campaigns now contain **no code at all**: the oracle is shared at
`npbench_c.templates.mass_balance_nrps`, so an instantiation is an entry file, a
`task.yaml` and two build commands. That is the template architecture the
catalog's authoring model assumes, and it was factored out at two
instantiations rather than at eight.

**T-L4-4** **Built**: `campaigns/chemspace-mibig-4_0-01`. Scaffold and
chemical-space analysis over MIBiG 4.0's active entries, with the identity key as
the subject rather than an implementation detail: 3,449 compound records are
3,115 full InChIKeys, 2,996 connectivity blocks and 1,843 ring scaffolds, and
every statistic downstream inherits the choice.
R1 the corpus and the three distinct-key counts · R2 what each key merges (111
connectivity blocks cover more than one full key, hiding 119 distinctions) and
the per-class spread · R3 the scaffold-sharing claim that survives a declared
informativeness filter, plus which scaffolds genuinely cross biosynthetic classes
· R4 four filter settings and the experimental-evidence gate.

Chemistry is perceived **once** with RDKit at instantiation and pinned into
`reference/chemistry_table.json` — formula, InChIKey, Bemis-Murcko scaffold, atom
and ring counts per record — exactly as monomer formulas were for mass balance.
The campaign itself is then stdlib set arithmetic, so **S0**.

*Two scope corrections made during the build, both forced by what the scaffolds
turned out to be.* **Bemis-Murcko reduction is far coarser than the catalogued
design assumed.** It keeps ring systems and discards substituents, so the most
widely shared scaffold in this corpus is plain benzene across 76 entries and the
second is tetrahydropyran across 26. "These 76 compounds share a scaffold" is
true and says nothing, and the catalogued R3 — the *distinguishing substitution*
between two compounds sharing a scaffold — is then almost the whole molecule, and
is in any case a set difference between molecular graphs that can only be named
in prose. R3 instead grades the reconciliation between the naive sharing
statistic and the one that survives a declared filter (≥2 rings and ≥50% heavy-atom
coverage): 279 shared scaffolds become 190, and 36 cross-class scaffolds become
14. Both sides are graded, so the size of the effect is visible rather than
hidden inside a choice. **And the catalogued R4** — predict the scaffold after a
substituent change — needs scaffold perception at solve time, which would make
the campaign S1 for one rung; R4 varies the declared thresholds and the evidence
gate instead.

One more counting rule the data forced: a scaffold "crosses classes" only when
two of its entries have **disjoint** class sets. 456 active entries carry more
than one biosynthetic class, so the obvious rule — the union of classes has more
than one name — is inflated more than threefold (115 against 36) by single hybrid
clusters crossing classes with nothing to compare against. Both numbers are
reported.

## L5 — Target & mechanism (3 templates, 3 campaigns)

Gold here is **measured-or-computed only**: a value from a database with a
structured assay-confidence field, or exact arithmetic over such values. Never a
model prediction, never a curated expert call. This is the layer where integrity
breaks if the rule is relaxed.

| # | Template | Data | Build |
|---|---|---|---|
| T-L5-1 | Selectivity ratios from measured IC50 | verified | **BUILT** (`share_alike`) |
| T-L5-2 | Resistance mutation → target assignment | **thin on theme** | re-tiered, see below |
| T-L5-3 | Binding-site geometry from coordinates | verified | **BUILT** |

**T-L5-1** ChEMBL measured activities. **Licensing decided**: CC BY-SA 3.0 is
accepted under `license_class: share_alike`, because ShareAlike attaches to an
adaptation — the derived subset and the gold computed from it — and not to a
collection of separable campaigns, so the obligation is local to these two rows
and they can be dropped without touching the rest. The gate requires such a
campaign to carry a `license_notice`. See `docs/environment.md` for the full
reasoning.

**Built**: `campaigns/selectivity-saureus-topoisomerase-01`, the benchmark's first
`share_alike` campaign. *S. aureus* DNA gyrase (CHEMBL3038482) against
topoisomerase IV (CHEMBL3038508), ChEMBL_37, 1,749 shipped activities of which
1,117 are IC50. R1 the activity and assay inventory, including the
confidence-score distribution · R2 matches gold · R3 the admission census, the
aggregated values, the per-threshold selective counts, the committed selective set
at the highest threshold and the extreme molecule under a declared tie-break ·
R4 four admission variants plus an aggregator swap.

Measured: 739 of the 1,117 IC50 records are admitted, and the rejections
partition exactly — 235 censored relations, 81 unconvertible units, 33 flagged by
ChEMBL, 29 potential duplicates. 668 (molecule, target) pairs carry an aggregated
value and **173 molecules are measured against both targets**. Both directions are
populated at twofold (122 favour gyrase, 21 favour topoisomerase IV), so the
verdict is not foregone either way and the build refuses to emit gold if either
side is empty. The most selective compound is CHEMBL4555272 at 16,666.7×. Every
one of the five R4 variants moves something: admitting censored relations takes
eligibility from 173 to 193, and swapping median for minimum leaves the admitted
count identical while moving the focal counts from 122/74/18 to 139/90/25 — a
reduction choice, not an admission choice, and it changes the answer.

**Deviation: the catalogued R4 is a no-op on this pair.** It was "a low-confidence
assay is excluded — predict the new ratio". All 105 assays behind these IC50
values score 7 ("Direct protein complex subunits assigned"), so excluding on
confidence rejects nothing. The rule stays declared and graded — a rule that
happens to be a no-op on one instantiation is not a rule that can be dropped from
the contract — but the counterfactual has to vary an axis that moves, so R4
relaxes the censored-relation rule, the ChEMBL flag, the duplicate rule and the
aggregator instead. **This is the second half of a rule worth stating in general:
a field being present does not make the quantity derived from it informative, and
only computing the distribution shows which.**

The campaign also carries the join that makes the task non-trivial:
`confidence_score` lives on the **assay**, not the activity, so admission needs an
activity-to-assay join and a system that looks for the field on the activity finds
nothing and admits everything.

**T-L5-2** Resistance mutation to target assignment from structured fields.
R1 mutation table · R2 matches gold · R3 claim: the implicated target and the
evidence class · R4 a mutation is removed from the set — predict whether the
assignment survives.

**Re-tiered: the structured field is there and the data behind it is not, for any
on-theme target.** ChEMBL_37 has 34,921 IC50 activities carrying
`assay_variant_mutation`, which is what the earlier probe confirmed and why this
row read `data: verified`. Counting them per target is what the probe did not do:

- *S. aureus* gyrase complex: **45** variant records (D83N 41). GyrA alone: 39
  (S84L 32, D83N 6). ParE and the topoisomerase IV complex: **0**.
- 44 wild-type/mutant key pairs exist, but most mutant `standard_value`s are
  `None`, so the pair cannot be turned into a ratio. The one clean pair in the
  whole set is ciprofloxacin against gyrase: WT 5.0 nM versus D83N 580,000 nM.
- *P. falciparum* DHFR: 244 variant records, **every one** of them
  "UNDEFINED MUTATION" — the field is populated and says nothing.
- The only well-populated variant target is **HIV-1 RT**: 5,423 records with real
  mutations (Y181C, P236L, K103N).

A campaign built on one ciprofloxacin pair is an anecdote, and a campaign built on
HIV-1 RT is off-theme for a natural-product benchmark — it would grade resistance
reasoning on an antiviral target whose chemistry is nothing like this benchmark's.
So this row is **not built and not blocked on licensing**: ChEMBL is accepted, the
data is thin. Named alternatives, left to the owner rather than chosen here:
instantiate on HIV-1 RT and accept the theme drift, re-scope the ladder to the
census itself (how much resistance data exists per target class, which is a real
and gradeable question), or wait for a release that fills the bacterial variant
values. **The numbers are recorded so the next person does not re-probe.**

**T-L5-3** Binding-site geometry from PDB coordinates (CC0, pure arithmetic).
**Built as the `pocket` ladder of the `structure_features` template**:
`campaigns/pocket-rebh-01`. R1 the structure inventory and the two numbering
offsets read from the alignment record · R2 the contact shell of every ligand
copy plus the two-chain symmetry control · R3 the pocket descriptors, a declared
lattice volume, and the shell reconciled against UniProt's annotations for this
entry · R4 a cutoff sweep and the whole shell truncated to alanine one residue
at a time.

**Deviation: the volume is a declared lattice count, not a published
descriptor.** Every published pocket-volume definition carries parameter choices
that would become part of the answer, so the rule is pinned in
`geometry_rules.json` — 1 Å lattice anchored on integers, 1.4 Å probe, Bondi
radii, three stated admission conditions — and the graded quantity is the
admitted point count. Reproducible to the point rather than approximately
comparable. Solvent-accessible surface area and any druggability score are
excluded for the same reason, recorded in the campaign's
`excluded_from_grading`.

## L6 — Translation & design (3 templates, 4 campaigns)

| # | Template | Data | Build |
|---|---|---|---|
| T-L6-1 ×2 | Construct & sequence design | verified | **BUILT** |
| T-L6-2 ×2 | Multi-constraint conflict + feasible control | verified | **BUILT** |
| T-L6-3 | Prioritisation under budget | likely | image |

**T-L6-1** **Built**: `campaigns/construct-ecoli-rebh-01`. Oracle 1.0 on all
four rungs, measured discrimination 1/2/3/3/4, binding constraint flips between
hosts. Now factored onto the shared `construct_design` engine and carrying no
code, so the second instantiation — PrnA (UniProt P95480, 538 aa, reviewed) — is
a data change plus one build command. Gold was verified byte-identical across
the factoring.

**T-L6-2** **Built**: `campaigns/constraint-conflict-ecoli-rebh-01`. The
declared constraint set is unsatisfiable and the correct answer is to abstain on
the design, name the conflicting pair and prove it.

Abstention gold is **constructed, not curated** — and crucially, *provable*. A
greedy optimiser failing to find a design proves nothing, so the conflict is
built where exact arithmetic decides it: for each residue independently the
highest-GC admissible codon is determined, and because codon choices are
independent the resulting whole-ORF bound is **attained** (a test constructs the
attaining sequence). A GC requirement above that bound cannot be met by any
sequence, whatever else applies, since adding constraints only shrinks the
feasible set.

For RebH: forbidding the eight GC3 codons drops the achievable maximum from
0.677338 to 0.584432, and requiring 0.62 is then impossible — while each
constraint alone is satisfiable. Margins of 0.036 and 0.057, so an auditor
recomputing it does not land on a rounding boundary. The build **refuses to emit
gold** if the requirement is still reachable, or if it exceeds even the
unrestricted maximum (which would make the conflict non-pairwise).

R4 asks, for each single relaxation, the new achievable maximum and whether the
requirement becomes reachable. Only the two GC-relevant constraints move it; the
five composition constraints restrict which sequences are admissible without
changing what GC any sequence can reach, and R4 tests exactly that distinction.

**The control pair, now built.** `campaigns/constraint-feasible-ecoli-rebh-01`
is the sufficient-evidence control, and it differs from the conflict campaign by
**one number**: the GC floor is 0.55 rather than 0.62. Same protein, same
forbidden-codon set, same host, same composition constraints. Reflexive
"infeasible" fails there; reflexive "feasible" fails here. Both properties are
tested at the grader rather than inferred, because both ablation stubs happen to
die at R2 before R3 matters.

Abstention precision and recall are reported separately across the pair, never
F1 alone.

**The asymmetry to keep in mind when reusing this template.** Infeasibility is
decidable by arithmetic; feasibility of the *full* set is only provable
**constructively**, by exhibiting a design, since the achievable-GC bound settles
the GC axis but not whether the composition constraints leave any candidate. The
feasible half therefore rests on a verified witness (GC 0.569994, zero
violations, protein preserved) recorded in gold for the audit packet and
deliberately **not graded** — many designs qualify, so grading one would grade
our optimiser rather than the agent's reasoning.

The build checks the **declared expectation** against the arithmetic and refuses
to emit gold when they disagree, in either direction, so a campaign cannot
silently become the opposite of what it was authored to be.

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

## Re-grounding, 2026-10-02

Every remaining `build: now` row has now been probed against its actual data
source. Full findings in `docs/regrounding-2026-10-02.md`; the flags above are
updated. Headline: **the S0 buildable set is 3 templates, not 7.**

- **Ready now**: T-L2-2 catalytic residues, T-L5-3 pocket geometry,
  T-L4-4 chemical space. The first two share their data sources (UniProt
  ECO-coded features plus PDB mmCIF) and the proteins already in the
  benchmark, so they should share one `structure_features` template with two
  ladders rather than duplicate an mmCIF parser.
  **Done**: all three are built. The shared `structure_features` template
  carries both ladders (`residues-prna-01`, `pocket-rebh-01`) and chemical space
  is `chemspace-mibig-4_0-01`. **The S0 set is exhausted**: every remaining row
  needs the Phase 1 container, a licensing decision on ChEMBL share-alike, or
  BRENDA credentials.
- **Blocked on one licensing decision** (ChEMBL CC BY-SA 3.0 share-alike):
  T-L5-1 and T-L5-2. Both are otherwise data-verified, including the structured
  `assay_variant_mutation` field T-L5-2 needs.
  **Decided**: accepted under `license_class: share_alike`. **T-L5-1 is built**
  (`selectivity-saureus-topoisomerase-01`). **T-L5-2 is re-tiered**: the field is
  populated, the values behind it are not for any on-theme target — 45 variant
  records on the *S. aureus* gyrase complex, 0 on topoisomerase IV, and all 244
  *P. falciparum* DHFR records reading "UNDEFINED MUTATION". Counts and the named
  alternatives are in the T-L5-2 entry above. **Correction to this bullet as
  written**: "otherwise data-verified" was wrong for T-L5-2 — the probe verified
  the field, not the data.
- **Blocked on registration**: T-L2-5, since BRENDA's archive 404s without it.
  **Wrong**: a licence-acceptance checkbox, not a registration, and the row is
  built.
- **Re-tiered to `image`**: T-L1-3, because finding a closest homolog needs
  alignment and any S0 route is a lookup.

So nine of the remaining rows need the Phase 1 container, which makes Phase 1
load-bearing sooner than this catalog implied. **Built since**: the image exists
for the six sequence-analysis tools plus antiSMASH and BiG-SCAPE, and clears 10
of the 12 `build: image` rows; matchms is the remaining provisioning work.

## Authoring order

1. **`build: now` + `data: verified`** — T-L4-3 (mass balance), T-L3-3 (MIBiG
   diff), T-L6-1 second instantiation, T-L6-2, T-L3-5. Five campaigns needing no
   container.
2. **`build: now` + `data: likely`** — T-L1-3, T-L2-2, T-L2-5, T-L4-4, T-L5-1,
   T-L5-2, T-L5-3. Probe each source first. **Done**: T-L2-2, T-L2-5, T-L4-4,
   T-L5-1 and T-L5-3 are built; T-L1-3 was re-tiered to `image`; T-L5-2 was
   re-tiered on thin data. This tier is now exhausted.
3. **`build: image`** — everything else. **Phase 1 is now partly done**: the
   pinned-tool image in `image/` carries HMMER 3.4, DIAMOND 2.2.8, MMseqs2
   18.8cc5c, Prodigal 2.6.3, MAFFT 7.526, BLAST+ 2.17.0, antiSMASH 8.0.4 and
   BiG-SCAPE 2.0.3, all eight verified thread-invariant at 1 and 8 threads, which
   **unblocks 10 of the 12 rows**: T-L1-1, T-L1-2, T-L1-3, T-L1-4, T-L2-1,
   T-L2-3, T-L2-4, T-L3-1, T-L3-2 and T-L3-4, plus the catalogued T-L3-5 target
   inference. The other two wait on matchms (T-L4-1, T-L4-2).

   **BiG-SCAPE does not run as its recipe publishes it**: the bioconda bound
   `sqlalchemy >= 2.0.2` resolves to 2.1.x, where 2.0.3 crashes before reading an
   input file. The image pins `sqlalchemy=2.0.54`. And the MIBiG reference set
   ships processed with antiSMASH **8.0 beta 1** while the image pins 8.0.4, so a
   T-L3-2 campaign must either reprocess it or take both sides from the published
   set — the cross-version rule applies to the reference data too.

   **antiSMASH's pin is a pair** and the lock records both halves: the binary at
   8.0.4 and the 9.4 GB database release it was fetched with. On top of that the
   lock carries a fingerprint over the three detection-rule files, because the rule
   count is now **103** (90 strict, 7 relaxed, 6 loose) against the 88 at v7.1 this
   catalog noted — so "never compare across versions" is enforceable rather than
   merely written down.

   Two determinism controls this catalog and the provisioning notes relied on
   turned out to be wrong, and both were caught by measuring rather than by
   reading: DIAMOND's `--no-reorder` *breaks* multithreaded determinism instead of
   providing it, and MMseqs2 cannot be made order-stable by any flag and needs a
   declared canonical line sort as part of its invocation contract. See
   `docs/environment.md` §2.

Author ~44 candidates to ship 32: the internal sweep will drop some for
non-discrimination and the audit will drop some for gold problems.
