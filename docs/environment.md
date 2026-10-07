# Tools, reference data and databases for the 32 campaigns

What the container must contain, what each campaign needs, and — the part that
decides whether the numbers mean anything — **what the agent is given versus
what is withheld**.

Status flags: `verified` = fetched and field-inspected this session ·
`likely` = open and standard, not yet inspected · `probe` = needs a feasibility
check before any campaign depends on it.

---

## 1. Databases and snapshots

| Snapshot | Version pin | License | Tier | Size | Status |
|---|---|---|---|---|---|
| **MIBiG** | 4.0 | CC BY 4.0 | shipped | 954 KB (3,013 JSON) | **verified** |
| **MIBiG** (for diffs) | 3.0 | CC BY 4.0 | shipped | ~1 MB | likely |
| **UniProtKB / Swiss-Prot** | pin release + per-entry `sequence_version` | CC BY 4.0 | shipped subset | small | **verified** |
| Pfam (HMM library) | 38.2 via InterProScan 5.78-109.0 | CC0 | shipped | ~1.5 GB | likely |
| InterPro | 109.0 | CC BY 4.0 | shipped subset | subset | likely |
| MassBank | 2026.03 | CC BY 4.0 | shipped | ~240 MB | likely |
| GNPS reference-standard libraries | frozen pull | CC0 (native) | shipped subset | subset | likely |
| CASMI 2022 | fixed release | public w/ solutions | shipped | probe | probe |
| **PDB (mmCIF)** | per-entry snapshot, entry id is the pin | CC0 | shipped subset | 0.5-1 MB per entry | **verified** |
| antiSMASH DB / ClusterBlast | pin to antiSMASH version | open | shipped subset | large | probe |
| **ChEMBL** | pin release | **CC BY-SA 3.0** → `license_class: share_alike` | shipped subset | subset | **verified, decided** |
| **BRENDA** | 2026.1 (March 2026) | CC BY 4.0 | **shipped subset** | 335 KB (EC 1.1.-) | **verified** |
| **NPAtlas** | 2024_09 | **CC BY-NC 4.0** | **fetched** | — | probe |

### Licence classes, and the ChEMBL decision

A campaign declares one of three classes, and the readiness gate checks it:

| Class | Meaning | Effect |
|---|---|---|
| `open` | attribution at most — CC BY, CC0, public domain | none |
| `nc` | non-commercial source | unavailable to commercial evaluators |
| `share_alike` | the shipped data is an **adaptation** of a ShareAlike source | those data files carry that licence onward |

**ChEMBL is CC BY-SA 3.0, and the decision is to use it** under
`license_class: share_alike`. The reasoning, so it can be overturned on its
merits rather than re-derived:

- ShareAlike attaches to an **Adaptation**, not to a **Collection**. A filtered
  ChEMBL subset is an adaptation, and so is gold computed from ChEMBL values, so
  both carry CC BY-SA 3.0. A benchmark made of separable campaigns is a
  collection, not an adaptation of any one of them — so the code, the grader and
  every other campaign are unaffected.
**In use, and the machinery is exercised end to end**:
`campaigns/selectivity-saureus-topoisomerase-01` (T-L5-1) is the first
`share_alike` campaign. Its `license_notice` names the source (ChEMBL_37 activity
and assay records for two targets), the licence, the attribution string and what a
redistributor must do, and the gate's `share_alike_notice_present` check passes on
it — so the nineteenth check is confirmed against a real campaign rather than a
fixture. **T-L5-2, the other row this decision was taken for, is not built**: the
licence is fine and the data is thin. See `docs/catalog.md`.

- The obligation is therefore satisfiable and **local**: attribution, the same
  licence on those files, and a notice. A share-alike campaign can be dropped
  without touching anything else, which is exactly why the class is per-campaign
  rather than per-benchmark.
- The alternatives do not help. IUPHAR/Guide to Pharmacology is CC BY-SA 4.0 —
  the same condition. PubChem BioAssay is largely a mirror of ChEMBL, so routing
  through it would be laundering rather than compliance. ChEMBL is also the source
  a reviewer expects for measured IC50s.
- **What is left for the owner to decide**: whether any evaluator's legal position
  rules out shipping ShareAlike content at all. If so, T-L5-1 and T-L5-2 come out
  and nothing else changes.

A `share_alike` campaign must carry a `license_notice` naming the source, the
licence, the attribution and what a redistributor must do; the gate fails the
campaign without it, because a licence condition recorded only in a design
document is a condition the person redistributing the files will never see.

**NPAtlas** is non-commercial from 2024_09, so any campaign touching it is tagged
`license_class: nc` and is unavailable to commercial evaluators.

**BRENDA is not behind credentials.** The download is a licence-acceptance
checkbox on `download.php` — no account, no API key — and the licence is plain
**CC BY 4.0**: *"All copyrightable parts of BRENDA are licensed under Creative
Commons Attribution License 4.0."* The project previously recorded this row as
blocked on credentials; that was wrong. One caveat worth carrying: BRENDA's terms
note benefit-sharing obligations for Digital Sequence Information under the UN
Convention on Biological Diversity, including contributions to the Cali Fund from
**commercial** users. That is not a copyright restriction and does not touch the
CC BY 4.0 grant, but a commercial evaluator should know it exists, so it is
recorded in the campaign's provenance for the same reason NPAtlas carries its
non-commercial flag.

**Excluded outright — do not build on these:** KEGG (bulk download is paid),
MetaCyc / BioCyc (subscription), DrugBank (restrictive terms). Pathway and
function tasks pull toward KEGG by reflex; use **Rhea** (CC BY), **ChEBI**
(CC BY) and **Reactome** (CC BY 4.0) instead. A readiness check already greps
`task.yaml` for these names and fails the campaign.

## 2. Tools, with their determinism controls

Every pin below is required, not advisory. Coupling is real: antiSMASH ↔
ClusterBlast DB ↔ MIBiG version ↔ pyhmmer move together.

| Tool | Needed by | Required control |
|---|---|---|
| **Python 3.11 + stdlib** | all | the grading core is stdlib-only by design |
| **RDKit** | L4-3, L4-4, L5-1 | pin exact (2026.03.6 verified installable). **Formula and InChIKey only — never canonical SMILES**, which changed at 2022.09, 2023.03, 2023.09, 2026.03.1 and is not guaranteed idempotent |
| **HMMER** | L1-2, L2-1, L3-4 | pin version, fixed `--cpu` |
| **BLAST+** | L1-3, L2-2 | pin version, fixed `-num_threads` |
| **DIAMOND** | L1-3, L2-3 | **≥ 2.2.7**, fixed `--threads`, explicit `--max-target-seqs` (it changes which hits survive). **Never `--no-reorder`** — measured, that flag is what breaks determinism, not what provides it (see below) |
| **MMseqs2** | L2-3 | pin post-fix build, always `createindex`, fixed `--threads`, explicit `--max-seqs` (issue #277: results differ by core count), **plus a declared canonical line order** — no flag fixes its multithreaded output order |
| **Prodigal** | L1-1, L1-4 | pin version **and mode** (single vs meta changes calls); Prodigal-GV is a different tool |
| **antiSMASH** | L3-1, L3-5 | pin the exact version **and** the database release — they are one pin. Fixed `--cpus`. **Annotated input only** (GenBank/EMBL, never FASTA), or the gene caller's determinism is folded into antiSMASH's. Grade a declared **projection** of the JSON, never the HTML and never the whole JSON. **Never compare across versions** — 58 rules at v5, 88 at v7.1, **103 at 8.0.4** |
| **BiG-SCAPE** | L3-2 | pin the version **and `sqlalchemy < 2.1`** — the recipe's own bound is wrong and 2.0.3 crashes on 2.1 before reading an input. Pin the antiSMASH version that produced the inputs. Fixed `--cores`, explicit `--gcf-cutoffs`, explicit `--include-gbk` when inputs are not named `*.region*.gbk`. **Never `--mibig-version`** on a gold path: it downloads a reference set at run time |
| **MAFFT** | L2-3 | pin major version, prefer deterministic modes |
| **matchms / pyOpenMS** | L4-1, L4-2 | pin versions, record every filter parameter |
| **FastTree** | — | **FastTreeMP is banned from the image** (thread order affects the NJ heuristic) |
| **IQ-TREE** | — | if ever used: always `--seed`, always fixed `-T N`, never `AUTO`. **Single-gene ML trees are not admissible as gold** |

`tests/test_thread_invariance.py` runs every tool at 1 and 8 threads, twice
each, on pinned fixtures and asserts identical output. Anything that fails and
cannot be pinned into determinism is made single-threaded in the image or
excluded from gold-producing paths.

### What measuring it found

Run against the image (`image/thread_invariance.json` records the verdicts), on
a 24-sequence halogenase fixture and a synthetic contig:

| Tool | Verdict | What it took |
|---|---|---|
| BLAST+ 2.17.0 | invariant | nothing; byte-identical at 1 and 8 threads |
| MAFFT 7.526 | invariant | nothing, at fixed `--retree`/`--maxiterate` |
| Prodigal 2.6.3 | invariant | nothing; single-threaded, so repeat determinism only |
| HMMER 3.4 | invariant | one declared normalisation: HMMER stamps its own command line, working directory and the wall-clock date into `--tblout` and into the HMM |
| **DIAMOND 2.2.8** | invariant **only without `--no-reorder`** | dropping the flag this project had declared as a control |
| **MMseqs2 18.8cc5c** | invariant **only after a declared canonical line sort** | no flag exists; the contract fixes the order instead |
| **antiSMASH 8.0.4** | invariant on both invocations | a declared **projection**: the raw JSON cannot be identical, since it holds the input path and a timestamp, so the comparison is on the detected regions and on the NRPS/PKS domain architecture |
| **BiG-SCAPE 2.0.3** | invariant | a declared **projection** to the GCF *partition*, plus an `sqlalchemy < 2.1` pin without which it does not run at all |

**The DIAMOND control was backwards.** This table previously listed
`--no-reorder` as a required determinism control. Measured: *with* the flag,
2.2.8 emits the same 192 hits in three different query orders across four
multithreaded runs; *without* it, every run at 1, 4 and 8 threads is
byte-identical. DIAMOND's default is to restore query order and `--no-reorder`
documents itself as switching that off for speed, so the flag was the cause and
not the cure. It is now forbidden on any gold-producing path, and a standing
witness invocation keeps the fact under test so a future release that fixes it is
noticed rather than assumed.

**MMseqs2 needs a canonical order, and loses nothing to it.** At one thread every
run is byte-identical; at eight the hit *multiset* is identical and only the order
moves. No flag changes that, so the order is fixed by the invocation contract: the
tabular output is sorted by line before anything reads it. That is sound only
because the artifact has no header and each line is an independent record, which
the declaration states — it would not be sound for an aligned or sectioned format.

Worth separating, because the suite reports them separately and they are different
properties:

- a **normalisation** says a line is not a result, so ignore it when comparing. It
  applies inside the suite only, and every rule carries a reason.
- a **canonicalisation** says the content is right but its order is not
  guaranteed, so fix the order before anybody reads it. It applies to the oracle
  and the agent alike, and a tool that needs one and does not get it is not
  admissible as a gold source.

And one reporting rule that only looks pedantic until it bites: a raw byte match
while a normalisation rule fired is reported as **incidental**, not as
raw-identical. `hmmbuild` writes a build date, and four runs inside the same
second agree on it; calling that raw-identical would turn a timing accident into a
claim.

## 3. What the agent is given

This is the section that just caught a real bug, so it is worth stating as a
rule rather than a convention.

**The sandbox contains exactly three things**, and nothing else:

```
inputs/      the task's data (sequences, spectra, entry JSON)
reference/   every pinned table the task requires the agent to use
task.yaml    objective, output contract, counting rules, closed vocabularies
submission/  empty, for the agent to fill
```

Withheld, and audited on every sandbox build: `gold/`, `oracle/`,
`grading.json`, `oracle_submission/`, `internal_sweep.json`. A violation raises
`GoldLeak` rather than warning.

**Warm-start runs** additionally receive a **pre-rendered bundle**
(`gold/warm/<rung>/`) holding only the earlier rungs' deliverables — never a
whole gold file, because copying `gold.json` to supply "the previous answer"
also supplies the answer the rung is being asked to compute.

### The rule the bug taught us

> **A campaign must be solvable from its own sandbox.** Every table, vocabulary,
> counting rule and convention that its gold depends on must be inside
> `reference/` or `task.yaml`.

The construct-design campaign shipped briefly with its codon weight table in
`oracle/`, which the sandbox excludes. R3's gold is a pure function of those
weights, so the campaign was **literally unanswerable** — and no run would have
told us: an agent would simply have failed and looked incapable. The stub
systems hid it because they were handed `--oracle` on the command line.

Two guards now exist: `pinned_tables_reachable_by_agent` in the readiness gate,
and a test that recomputes R3's gold using only files visible in the sandbox.
The oracle loads the same `reference/` JSON the agent gets, so the two cannot
drift.

This is also why **stub sweeps cannot promote a readiness check to PASS**. A
fixture given privileged access proves nothing about solvability.

### "Skills" — the procedural material that ships with a campaign

Agents are given no tool documentation beyond what their own surface provides,
but every campaign ships the **task-specific knowledge its gold assumes**:

| Material | Where | Example |
|---|---|---|
| Output contract | `task.yaml` | coordinate convention, identifier namespace, sort order, precision, units, null-vs-abstain, whether a stop codon is included |
| Counting rules | `task.yaml` | "direct_repeat = distinct 11-mers occurring more than once" |
| Closed vocabularies | `task.yaml` | the 6 constraint names; the 18-entry confound vocabulary; `modification_domains` enum |
| Pinned tables | `reference/` | codon weights, host GC ranges, forbidden-site lists, monomer formula tables |
| Derivation rules | `reference/` | "unconstrained optimum = highest-weight codon, ties to lexicographically smallest" |

If a rung's gold depends on a convention, that convention is published. An
unstated convention is a grader defect, not a hard campaign — it is item 4 on
the auditor's form ("is there a legitimate alternative answer this contract
would mark wrong?") and the single most valuable thing a collaborator can find.

## 4. Tool surfaces (for the model × tool-surface design)

| Surface | Contents |
|---|---|
| **S0 bash-only** | Python 3.11 + stdlib, coreutils. No domain tools, no RDKit |
| **S1 generic bioinformatics** | S0 + BLAST+, HMMER, DIAMOND, MMseqs2, MAFFT, Prodigal, Biopython, RDKit, matchms |
| **S2 NP-specialised** | S1 + antiSMASH, BiG-SCAPE, pyOpenMS |
| **S3 agent frameworks** | Biomni, Mycel — each with its own tool set |

**A surface is an installation, not a set of commands.** T-L2-4's grounding
produced the first campaign design that needed **S2 while invoking no tool at
all**: the thing it wanted from antiSMASH was two pinned tables, not a run.
The scale still reads correctly — a campaign that needs antiSMASH's database
layer needs S2 whether or not it ever calls `antismash` — but the declaration
should say which it is, so the model × surface matrix does not read a table lookup
as a tool-execution result.

**And data that cannot be named cannot be shipped.** Of the two tables that design
wanted, `modules/nrps_pks/data/aaSMILES.txt` is antiSMASH package data under
**AGPL-3.0-or-later** and `databases/nrps_pks/stachelhaus/1.1/signatures.tsv`
arrives from the database download with **no stated licence**. So neither goes into
a campaign's `inputs/` or `reference/`: an S2 campaign reads such a resource from
the pinned image at a declared path with a recorded sha256 and redistributes
nothing. Pulling AGPL material into a campaign directory is a licensing decision
about the whole repository, and it is not one to make by accident while writing an
oracle.

**The integrity rule:** a campaign must declare which surfaces can possibly
solve it. If a campaign needs HMMER and S0 lacks HMMER, S0's failure measures
tooling, not capability — which is a legitimate and interesting measurement, but
only if it is *declared*, so the matrix reads as a finding rather than a
confound. Add `min_tool_surface` to `task.yaml` per campaign.

Of the current catalog, **12 of 24 templates are `build: now`** — their oracles
are pure computation over downloadable data, so they are solvable at **S0** and
measure reasoning rather than tool access. That makes them the most
discriminating campaigns in the set for the question you actually want answered,
and the right ones to author first.

## 5. Per-template requirements

| Template | Databases | Tools | Min surface |
|---|---|---|---|
| T-L1-1 gene calling | — | Prodigal | S1 |
| T-L1-2 domain architecture | Pfam | HMMER | S1 |
| T-L1-3 annotation transfer | UniProt 2026_03 EC 1.14 (shipped, closed) | DIAMOND 2.2.8 | **S1, BUILT** |
| T-L1-4 pseudogene | — | Prodigal | S1 |
| T-L2-1 EC misannotation | UniProt 2026_03 (shipped), Pfam 35.0 (image) | HMMER 3.4 | **S1, BUILT** |
| T-L2-2 catalytic residues | PDB, UniProt | — (stdlib only) | **S0** |
| T-L2-3 remote homology | UniProt subset | MMseqs2/DIAMOND | S1 |
| T-L2-4 A-domain specificity | *retired on grounding — see `docs/catalog.md`* | — | — |
| T-L2-5 kinetics consistency | BRENDA 2026.1 (shipped subset) | — (stdlib only) | **S0** |
| T-L3-1 BGC detection | antiSMASH DB | antiSMASH | S2 |
| T-L3-2 GCF cutoff | MIBiG, antiSMASH DB | BiG-SCAPE | S2 |
| T-L3-3 MIBiG diff | MIBiG 3.0 + 4.0 | — | **S0** |
| T-L3-4 RiPP precursor | MIBiG, Pfam | HMMER | S1 |
| T-L3-5 self-resistance | MIBiG 4.0 | — | **S0** |
| T-L4-1 MS² dereplication | MassBank, GNPS, CASMI | matchms, RDKit | S1 |
| T-L4-2 molecular networking | GNPS | matchms | S1 |
| T-L4-3 mass balance | MIBiG 4.0, ChEBI | RDKit (formula only) | **S0/S1** |
| T-L4-4 chemical space | MIBiG 4.0 | RDKit (instantiation only) | **S0** |
| T-L5-1 selectivity ratios | ChEMBL | — (arithmetic) | **S0** |
| T-L5-2 resistance → target | ChEMBL/UniProt | — | **S0** |
| T-L5-3 pocket geometry | PDB, UniProt | — (stdlib only) | **S0** |
| T-L6-1 construct design | UniProt | — | **S0** |
| T-L6-2 constraint conflict | UniProt | — | **S0** |
| T-L6-3 prioritisation | varies | varies | S1 |

### Chemistry perceived once, then pinned

Two templates need a chemistry toolkit for *perception* and nothing else, and
both handle it the same way: RDKit runs once at instantiation, its output is
written into the campaign's `reference/`, and the campaign thereafter needs no
toolkit at all. That keeps both at **S0**, makes the oracle and the agent read the
same table, and means a toolkit upgrade cannot move a published score.

| Template | Script | Pinned table | What it carries |
|---|---|---|---|
| T-L4-3 mass balance | `mass_balance_nrps.build_reference` | `monomer_formulas.json` | a molecular formula per A-domain substrate |
| T-L4-4 chemical space | `chemical_space.build_reference` | `chemistry_table.json` | per compound record: formula, InChIKey, Bemis-Murcko scaffold SMILES and InChIKey, heavy-atom and ring counts (3,449 records, ~1.2 MB) |

The rule both follow: **RDKit output is admissible as a pinned input only where
it is stable across versions.** A molecular formula is atom counting plus
valence-implied hydrogens, and an InChIKey is a versioned standard. A canonical
SMILES is neither, so `scaffold_smiles` ships as a human-readable label and
`scaffold_inchikey` is the identity key — a distinction that is not academic here,
since the two give 1,843 and 1,850 distinct scaffolds over the same records.

Re-running a `build_reference` script is an instantiation-time act, not part of
the gold build: `chemical_space.build` refuses outright if the table is absent
rather than quietly deriving it.

### What the two structure-features ladders actually need

Worth stating precisely, because "structural biology task" reads like a request
for a toolkit and neither ladder needs one.

**Databases.** Two files per campaign and nothing else: a UniProtKB entry as
JSON (`rest.uniprot.org/uniprotkb/<accession>.json`, CC BY 4.0) and one or two
mmCIF coordinate files (`files.rcsb.org/download/<id>.cif`, CC0). Both are
committed into the campaign's `inputs/`, so the sandbox needs no network and the
pin is the file's own hash, recorded in `inputs/provenance.json`. A PDB entry id
is itself a version pin — deposited coordinates do not change under an id — so
no release snapshot is needed, unlike MIBiG.

**Tools.** None. No Biopython, no gemmi, no PyMOL, no RDKit, no alignment. The
mmCIF parser, the numbering mapper and the geometry are stdlib, which is why
both campaigns declare `min_tool_surface: S0`. That is a deliberate
constraint rather than an accident of what was available: the moment pocket
geometry runs through a toolkit, the toolkit's defaults — hydrogen treatment,
altloc resolution, symmetry expansion, radii — become part of the answer, and
the benchmark grades the toolkit's conventions instead of the system's
reasoning.

**What must be given to the agent.** Four pinned tables, all under
`reference/` and all audited as reachable by the readiness gate:

| Table | What it fixes |
|---|---|
| `eco_policy.json` | the functional feature vocabulary, the ECO code → evidence class map, and the counting rules (record identity, position unit, citation unit, census class vs policy admission) |
| `geometry_rules.json` | atom selection (model, altloc, hydrogens, waters), the contact cutoff and reported precision, the centroid and radius-of-gyration reductions, the lattice volume rule with its probe and radii, the alanine-truncation rule, and the numbering rule |
| `residue_alphabet.json` | three-letter → one-letter, including the modified residues deposited coordinates actually contain |
| `campaign_keys.json` | the campaign's declared join keys: ChEBI ligand id → PDB chemical component id, the focal structure and chain, the policy and cutoff variants |

The ChEBI → component join is the one that is easy to miss. UniProt names a
ligand by ChEBI; coordinates name it by a three-character component id; nothing
in either file states the correspondence. It is a convention, not a finding, so
it is declared — and declaring it gives nothing away, because the positions are
still to be computed.

**What must NOT be given.** The pocket ladder's R4 truncates every residue in
the focal contact shell, and the shell is R3's answer. Listing those positions
in `campaign_keys.json` would hand it over, so the scope is declared by rule
(`truncation_scope: focal_shell`) and a test asserts that no number in the
sandbox's declared keys is a shell position.

## 5b. The Phase 1 image

`image/` holds the pinned-tool image and its lock:

| File | What it is |
|---|---|
| `image/Dockerfile` | Ubuntu 24.04 plus micromamba, creating `/opt/npbench-tools` from exact `version=build` specs. Fails the build if `registry --verify` reports pin drift |
| `image/environment.lock.json` | the resolved 79-package closure with every artifact URL and sha256, so the image is rebuildable byte-for-byte rather than merely version-matched |
| `image/thread_invariance.json` | the measured verdicts, recorded the way a campaign records `internal_sweep.json` |
| `src/npbench_c/tools/registry.py` | per tool: the pin, what needs it, the controls it must run under, the invocation the suite exercises, and any declared normalisation or canonicalisation |
| `src/npbench_c/tools/invariance.py` | the suite |
| `src/npbench_c/tools/fixtures/` | the pinned inputs, with hashes and provenance |

The full suite takes about four minutes, most of it antiSMASH and BiG-SCAPE.
BiG-SCAPE is marked **slow**: `invariance --all` always measures it and its
verdict is committed, but the default test run skips it with a reason naming the
runtime, and `NPBENCH_SLOW_TOOLS=1` re-measures. To skip every invariance test,
point `NPBENCH_TOOL_PREFIX` at a directory that does not exist. In both cases the
test reports SKIP with its reason, which is the honest form of "not measured
here" — and distinct from a green tick.

**The Dockerfile has not been built end to end.** The host this was written on
had the Docker CLI but no daemon. Every step it performs was run directly instead
— the micromamba download, the create with these exact `version=build` specs,
`registry --verify` against the resulting prefix, and the full invariance suite,
whose verdicts are recorded — so the pins and the determinism results are real
while the layer sequence is not yet proven. Build it once before relying on it.

**Ubuntu's archive cannot satisfy the pins**, which is why the image uses
conda-forge and bioconda: `apt` ships DIAMOND 2.1.9 where the controls require
≥ 2.2.7, and BLAST+ 2.12.0 against a 2.17.0 pin. HMMER, Prodigal and MAFFT would
have been fine from `apt`; mixing sources for some tools and not others would make
the lock harder to reason about than it is worth.

**What this image unblocks: 10 of the 12 `build: image` templates** — T-L1-1,
T-L1-2, T-L1-3, T-L1-4, T-L2-1, T-L2-3, T-L2-4, T-L3-1, T-L3-2 and T-L3-4, plus
the catalogued T-L3-5 target inference that was re-tiered to `image` during the
S0 build.

**What it does not**, a recorded decision rather than an omission:

| Missing | Why it is a separate task |
|---|---|
| matchms / pyOpenMS (T-L4-1, T-L4-2) | the tools are pip-installable, but the rows need MS² fixtures from MassBank first, and a spectral-matching fixture is its own grounding job: for a spectral match the invariance question becomes a tolerance question |

### BiG-SCAPE, and a dependency bound the upstream recipe gets wrong

**It does not run as published.** The bioconda recipe for 2.0.3 asks for
`sqlalchemy >= 2.0.2` with no upper bound, so a fresh solve installs 2.1.x, and
BiG-SCAPE then raises `ObjectNotExecutableError` **before it reads a single input
file**: it passes a compiled statement object to `Connection.execute`, which 2.0
tolerated and 2.1 rejects. The image carries `sqlalchemy=2.0.54`, a bound the
upstream recipe does not. The general lesson is cheap to state and was not cheap
to find: **installing is not the same as working**, and only running the tool
tells them apart — which is the whole argument for this suite existing before any
campaign is authored on these tools.

Two smaller things the same run surfaced. BiG-SCAPE filters its input directory by
filename, defaulting to `cluster,region` to match antiSMASH's `*.region001.gbk`
convention; MIBiG reference files are named `BGC0000852.gbk`, so without an
explicit `--include-gbk` the run fails with *no valid input GBKs* rather than with
anything about names. And BiG-SCAPE needs `Pfam-A.hmm`, which the antiSMASH
databases already ship with its pressed indexes — so the image reuses that copy,
saving 1.5 GB and one separately-pinned release.

**The GCF partition is the result; the family labels are not.** `FAM_00001` is a
name. Two runs can agree completely on which clusters belong together and disagree
on what the groups are called, and comparing labels would score that as a
difference. So the projection keeps, per class, the groups as sorted member sets,
sorted among themselves, and drops the family label and the connected-component
number. The cutoff stays in the *invocation*, not the projection, because a
partition at a different cutoff is a different answer rather than a different
rendering of one.

**The fixture has designed family structure.** 32 antiSMASH-processed MIBiG
reference clusters, selected by a rule computed from tables already in this
repository: the four largest groups of 4–10 entries sharing an InChIKey
connectivity block — ectoine (10), ochratoxin A (7), kanamycin (5), aflatoxin B1
(5) — plus one singleton per biosynthetic class as a negative control. BiG-SCAPE
recovers those families at a 0.3 cutoff, which is what makes the fixture a test of
clustering rather than a test of whether 32 singletons stay 32 singletons.

One judgement in that selection is named rather than left implicit: groups whose
MIBiG compound name is a **class placeholder** — *capsular polysaccharide*,
*lipopolysaccharide*, *melanin*, *carotenoid*, *exopolysaccharide* — are excluded,
because MIBiG gives those generic names a representative structure, so entries
sharing one are not necessarily homologous clusters. Two different molecular
formulas appear under *capsular polysaccharide* alone. Including them would have
put apparent family structure in the fixture that the biology does not support.

**A version coupling worth stating before a campaign is built on it.** The MIBiG
reference set is published as `mibig_antismash_4.0_gbk_as8b1.tar.bz2` — processed
with antiSMASH **8.0 beta 1**, while this image pins **8.0.4**. For a determinism
fixture that is harmless, since both arms of every comparison read the same files.
For a T-L3-2 campaign it is not: comparing its own 8.0.4 regions against that
reference set is exactly the cross-version comparison antiSMASH's own rules
forbid. Either the reference set is reprocessed with the pinned antiSMASH, or the
campaign states that both sides come from the published set.

### antiSMASH, and the three things adding it forced

**Its pin is a pair, so the lock records both halves.** The binary is 8.0.4; the
databases are what `download-antismash-databases` fetched — knownclusterblast 4.0,
Pfam 35.0, MITE 1.3, as-js 0.16, plus clusterblast, clustercompare, comparippson,
nrps_pks, resfam and tigrfam — **9.4 GB on disk**, recorded in the lock under
`antismash_databases` with the directory versions. Recording only the binary
version would leave half the pin unrecorded, and the half left out is the half
that decides which regions come out.

On top of the version, the lock records a **fingerprint over the detection rule
files**: `04add3eb0e86a816`, a sha256 over `strict.txt`, `relaxed.txt` and
`loose.txt`. The rule count is **103** here (90 strict, 7 relaxed, 6 loose),
against 88 at v7.1 and 58 at v5. "Never compare across versions" is a standing
instruction, and the fingerprint is its mechanical form: a campaign pins what it
was built against, so a changed rule set is detectable even when the version
string is not what moved.

**A banned binary arrives with it.** antiSMASH depends on the `fasttree` package,
which ships `FastTreeMP` — banned outright, because thread order affects its
neighbour-joining heuristic. The image deletes the binary after install and
`registry --verify` fails if it reappears, which is the difference between a ban
and a note. Single-threaded `FastTree` stays, because antiSMASH needs it.

**It brings its own Python, and that nearly became the benchmark's.** The prefix
now contains a conda interpreter. With the tool prefix first on `PATH`, `python3`
silently resolves to it — and the stdlib-only grading core would then be running
on an interpreter nobody chose. Measured the hard way: `python3 -m pytest` stopped
finding pytest. antiSMASH's console scripts carry an absolute shebang to their own
interpreter, so the prefix does not need to come first at all; the image now puts
it **last** on `PATH`.

### The first campaign that runs a tool, and what it needed from the image

`ecaudit-fmn-dh-1_1_3_15-01` (T-L2-1) is the first campaign whose gold is a tool's
output rather than arithmetic over shipped data, so it is the first real test of
whether the image is usable from a campaign rather than only from the invariance
suite. Two things it needed, both now mechanisms rather than one-offs:

**A resource read from the image, not shipped.** The audit scans against Pfam
35.0, which is 1.5 GB. The campaign declares the path, the provider (the antiSMASH
8.0.4 database layer), the version and the file's sha256 in
`reference/audit_rules.json`, and the build refuses to emit gold if what it finds
does not match — a domain verdict against a different Pfam release is not that
campaign's gold. The readiness gate's twentieth check verifies the declaration.
Here the reason for not shipping is **size**, and Pfam is CC0; the note says so
explicitly, because the T-L2-4 grounding established that the next campaign to use
this mechanism will be withholding for **licence**, and the two cases need to be
told apart by anyone reading the campaign later.

**The declared panel, so a tool campaign stays affordable.** A full Pfam-A pass
over the campaign's 1,346-sequence subset takes **2m36s**; the 13-model panel the
campaign publishes takes **11 seconds** over all 7,675. The panel was fixed once at
build time by the full scan and then published in `reference/`. Any future S1 or S2
campaign wants the same shape: do the open-ended scan once at instantiation,
publish what it found as a declared vocabulary, and have the agent run the bounded
version.

That alone was not enough for the sweep, which re-ran the same call for every stub
level, rung and repeat: **29m32s** for one campaign. `npbench_c.tools.cache`
memoises tool output for **stub systems in tooled mode only** — never a real
system, never a `no_tool` ablation, with the runner raising rather than allowing
it — keyed on the input digest, the library fingerprint, the tool version and the
flags. That brought the sweep to **3m44s** with identical gate results. A cached
sweep's wall clock is consequently not a cost measurement of the campaign, which
the report states in `_meta.tool_cache_note`.

### Two shapes of S1 campaign, and the difference matters

The two built S1 campaigns need the image for different reasons, and the
declaration should say which:

- `ecaudit-fmn-dh-1_1_3_15-01` needs **the binary and a pinned data resource**:
  HMMER plus the 1.5 GB Pfam library it scans against, declared as an external
  resource with a fingerprint and read from the image.
- `transfer-ec1_14-uniprot-01` needs **only the binary**. DIAMOND builds its
  reference database from the shipped FASTA, so nothing outside the sandbox and
  the tool itself is involved and the campaign declares no external resource at
  all.

The second shape is the one to prefer where a campaign can be written either way:
a closed, shipped reference set makes gold a function of the sandbox, whereas a
campaign searching a library in the image inherits that library's version as part
of its gold. Both are legitimate; only one can be re-run on a machine that never
downloaded 9.4 GB.

### The other half of the PATH lesson: a measured run needs the prefix on it

Putting the prefix last was right and not sufficient. The invariance runner invoked
`{bin}/tool` by absolute path and left the prefix off the subprocess `PATH`
entirely — and antiSMASH resolves nine helpers on `PATH`: hmmscan, hmmsearch,
hmmpress, hmmpfam2, blastp, makeblastdb, diamond, prodigal and FastTree. So the
measurement read whatever the measuring shell happened to export. On the shell that
produced the original 8/8 verdict the prefix was exported; on a clean one
`antismash --check-prereqs` reports **44 prerequisite failures** across those nine
executables, every one of them in the prefix, and both antiSMASH invocations fail
4/4 runs with `RuntimeError: Modules failing prerequisites`.

The verdicts were right and the way they were obtained was not, which is the worse
failure of the two: **a determinism measurement that depends on the ambient
environment is not a measurement of the image.** `_image_path()` now builds the
subprocess `PATH` as the ambient one with the prefix appended, de-duplicated — last,
so the earlier lesson still holds — and a test pins both the position and the
uniqueness. `image/thread_invariance.json` is regenerated under it: 8/8, 10
invocations, 1 declared xfail, unchanged verdicts honestly obtained.

### Projections: a third kind of declared transform

antiSMASH forced one more distinction. Its output is a single JSON object holding
the input path, the tool version, a record timestamp and the full HMM hit table in
the same structure as the detected regions. No line-drop rule can reach inside
that, so the comparison is made on a **declared projection** instead:

| Invocation | Projection | What it states is a result |
|---|---|---|
| `minimal_detection` | `antismash_regions` | per record, each region's coordinates and sorted products |
| `default_modules_domains` | `antismash_nrps_pks_domains` | per CDS, the ordered domain architecture — without the e-values and bitscores, whose last digits are a reduction-order artefact |

Like a canonicalisation and unlike a normalisation, a projection is part of the
invocation contract: it declares which fields this benchmark treats as the answer,
so a campaign built on the tool grades the fields the suite proved stable, and
everything outside the projection is by that declaration **not a result**. The
suite reports a projected invocation as `projection-identical`, never as
`raw-identical`, because the raw artifact demonstrably is not.

The fixture is three complete deposited records concatenated — *Vibrio
anguillarum* pJM1 (anguibactin, NRPS), a *Kamptonema* landornamide cluster
(ribosomal) and a *Streptomyces sampsonii* julichrome cluster (PKS) — so the rules
exercised are not all of a kind: the run yields `NRP-metallophore`+`NRPS`,
`lanthipeptide-class-ii`+`proteusin` and `T2PKS`. All three are whole records
whose MIBiG locus starts at position 1, so the provenance is three accessions with
no coordinates to get wrong, and all three are annotated, so the verdict is about
antiSMASH rather than about a gene caller feeding it.

The fixtures are worth one note. They are committed files with hashes, because a
determinism suite whose input moves measures nothing. `halogenases.faa` is 24
reviewed UniProt entries carrying Pfam PF04820 — real homologs spanning the
similarity range, including RebH and PrnA, which the benchmark already ships, plus
one fragment, since a short record is where a tool's tie-breaks show.
`synthetic.fna` is reverse-translated from that set under a declared codon cycle,
so Prodigal has a contig with 24 known ORFs that is reproducible from the repo
with no download. `halogenases.afa` is a frozen MAFFT alignment; that it came from
a tool the suite also tests is not circular, because once written it is a file —
a later MAFFT regression would not change the HMMER results computed from it.

## 6. Open provisioning questions

1. **Pfam / InterProScan is the largest single dependency** (~1.5 GB+). Confirm
   whether a Pfam subset restricted to biosynthetic domains suffices; if so the
   image shrinks substantially.
2. **antiSMASH + ClusterBlast DB size** must be measured before committing to
   the four S2 campaigns — it may dominate the image.
3. **ChEMBL share-alike** needs a licensing decision before L5 is authored.
4. **CASMI 2022** distribution terms and file layout need a probe.
5. **Container parallelism target.** A full sweep at 32 campaigns × cold+warm ×
   n=3 × the chosen surfaces runs into thousands of container-hours; at 100
   parallel containers that is days, at 10 it is weeks. This determines the
   orchestration design and cannot be retrofitted cheaply.
6. **`MIBiG quality` semantics.** 2,710 of 3,013 entries are `questionable`,
   which appears to track curation status rather than scientific doubt. Confirm
   with MIBiG curators before the field is cited anywhere; meanwhile all
   filtering is on `loci[].evidence`, which is explicit.
