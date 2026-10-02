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
| **ChEMBL** | pin release | **CC BY-SA 3.0** | shipped subset | subset | likely |
| **BRENDA** | 2026.1 | CC BY 4.0 **+ acceptance gate** | **fetched** | — | probe |
| **NPAtlas** | 2024_09 | **CC BY-NC 4.0** | **fetched** | — | probe |

**Two licensing cautions.** ChEMBL's **share-alike propagates to redistributed
subsets** — confirm with whoever owns licensing before it becomes load-bearing
for L5. NPAtlas is non-commercial from 2024_09, so any campaign touching it is
tagged `license_class: nc` and is unavailable to commercial evaluators.

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
| **antiSMASH** | L3-1, L3-5 | pin exact version + ClusterBlast DB. **Never compare across versions** — detection rules went 58 → 88 across v5–v7.1 |
| **BiG-SCAPE** | L3-2 | pin version, `--mibig-version`, antiSMASH version, pyhmmer version; set classification mode explicitly |
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
| T-L1-3 annotation transfer | UniProt subset | BLAST+/DIAMOND | S1 |
| T-L1-4 pseudogene | — | Prodigal | S1 |
| T-L2-1 EC misannotation | UniProt, Pfam | HMMER | S1 |
| T-L2-2 catalytic residues | PDB, UniProt | — (stdlib only) | **S0** |
| T-L2-3 remote homology | UniProt subset | MMseqs2/DIAMOND | S1 |
| T-L2-4 A-domain specificity | MIBiG 4.0 | HMMER | S1 |
| T-L2-5 kinetics harmonisation | BRENDA (fetched) | — | **S0** |
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

**What this image unblocks: 8 of the 12 `build: image` templates** — T-L1-1,
T-L1-2, T-L1-3, T-L1-4, T-L2-1, T-L2-3, T-L2-4 and T-L3-4, which between them
need only Prodigal, HMMER, DIAMOND, BLAST+, MMseqs2 and MAFFT.

**What it does not**, each a recorded decision rather than an omission:

| Missing | Why it is a separate task |
|---|---|
| antiSMASH (T-L3-1, T-L3-5) | its pin is a *pair* — binary plus ClusterBlast and Pfam database releases — that has to be resolved together, and detection rules went 58 → 88 across v5–v7.1, so results must never be compared across versions |
| BiG-SCAPE (T-L3-2) | its pin depends on the antiSMASH and pyhmmer versions, so it follows antiSMASH |
| matchms / pyOpenMS (T-L4-1, T-L4-2) | the tools are pip-installable, but the invariance suite would need MS² fixtures from MassBank, and a spectral-matching fixture is its own grounding job |

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
