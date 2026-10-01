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
| PDB | pin weekly snapshot | CC0 | shipped subset | subset | likely |
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
| **DIAMOND** | L1-3, L2-3 | **≥ 2.2.7**, `--no-reorder`, fixed threads; note `-k/--max-target-seqs` changes which hits survive |
| **MMseqs2** | L2-3 | pin post-fix build, always `createindex`, fixed `--threads`, explicit `--max-seqs` (issue #277: results differ by core count) |
| **Prodigal** | L1-1, L1-4 | pin version **and mode** (single vs meta changes calls); Prodigal-GV is a different tool |
| **antiSMASH** | L3-1, L3-5 | pin exact version + ClusterBlast DB. **Never compare across versions** — detection rules went 58 → 88 across v5–v7.1 |
| **BiG-SCAPE** | L3-2 | pin version, `--mibig-version`, antiSMASH version, pyhmmer version; set classification mode explicitly |
| **MAFFT** | L2-3 | pin major version, prefer deterministic modes |
| **matchms / pyOpenMS** | L4-1, L4-2 | pin versions, record every filter parameter |
| **FastTree** | — | **FastTreeMP is banned from the image** (thread order affects the NJ heuristic) |
| **IQ-TREE** | — | if ever used: always `--seed`, always fixed `-T N`, never `AUTO`. **Single-gene ML trees are not admissible as gold** |

`tests/test_thread_invariance.py` (Phase 1) must run every tool at 1 and 8
threads on fixed input and assert identical output. Anything that fails and
cannot be pinned into determinism is made single-threaded in the image or
excluded from gold-producing paths.

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
| T-L2-2 catalytic residues | PDB, UniProt | — (geometry) | **S0** |
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
| T-L4-4 chemical space | MIBiG 4.0 | RDKit | S1 |
| T-L5-1 selectivity ratios | ChEMBL | — (arithmetic) | **S0** |
| T-L5-2 resistance → target | ChEMBL/UniProt | — | **S0** |
| T-L5-3 pocket geometry | PDB | — (geometry) | **S0** |
| T-L6-1 construct design | UniProt | — | **S0** |
| T-L6-2 constraint conflict | UniProt | — | **S0** |
| T-L6-3 prioritisation | varies | varies | S1 |

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
