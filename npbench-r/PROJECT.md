# NPBench-R: a computational reproducibility benchmark for natural products and biosynthetic gene clusters

Status: design specification, pre-implementation.
Owner: Saman, Zhao lab, UIUC.
Name settled: **NPBench-R**, R for reproducibility. `NPBench` alone denotes part one of the suite (the curated Q&A set) and has also been used by an unrelated group for a 50-question type II polyketide QA set in the ChatT2 paper (arXiv 2609.25620), so the suffix is load-bearing. Always write the suffix.

## 1. What this is

A benchmark that measures whether an AI agent can reproduce the computational analyses of a published natural products or BGC study, starting from the raw deposited data and ending at structured data artifacts that are graded programmatically against the artifacts the original authors reported.

An agent receives three things: the research objective, method guidance naming the tools and parameters the original paper used, and the raw data staged on disk. It must produce a specified set of output files. Each file is compared to the corresponding published file with a deterministic grader. The task score is the fraction of artifacts that pass.

This is part two of the NPBench suite. Part one tests domain knowledge in a question-and-answer format. Part two tests whether an agent can actually run the analysis.

## 2. Why it does not exist yet

A survey of the 2025 to 2026 agentic biology benchmark landscape found nothing that does this for natural products.

The closest neighbours:

- **BixBench3** (arXiv 2608.25286, Edison Scientific) is the direct methodological template. Twenty tasks drawn from published studies, 138 graded artifacts, agents scored 0.00 to 0.48. Its twenty tasks cover transcriptomics (16), epigenomics (6), proteomics (2), genomics (1) and microbiome sequencing (1). No mass spectrometry, no cheminformatics, no BGCs.
- **MetaBench** (arXiv 2510.14944, ACL 2026) is a metabolomics benchmark, but it is roughly 8,000 knowledge and reasoning test cases built from HMDB, KEGG, PathBank and MetaboLights. Knowledge QA, not reproduction.
- **CompBioBench**, **GenoTEX**, **BioDSA-1K**, **VariantBench**, **GeneBench-Pro**, **scBench-Long**, **SpatialBench-Long** are all computational biology agent benchmarks with no NP or BGC content.
- **CHAMOIS**, **BGC-MLM** and **BGC-Master** are prediction models with their own held-out sets, not agent benchmarks.

The gap is specific and defensible: no benchmark asks an agent to run antiSMASH, cluster the output with BiG-SCAPE, network a set of MS/MS spectra, or link a gene cluster to a molecule.

## 3. What BixBench3 got right, and what we take from it

Three design decisions are worth copying directly.

**Prescribe the method.** The task prompt names the tools, versions, parameters, contrasts and filtering rules from the original paper, without revealing its results. Without this, an agent could pick a scientifically valid alternative and score zero because its output differs in format or scale. The cost is that the benchmark tests execution, not analysis design. BixBench3 names this in its own limitations section.

**Grade artifacts, not answers.** Each paper is decomposed into a directed acyclic graph of analyses, with raw data at the root. Every node is a structured file. This gives partial credit, makes failures diagnosable, and sidesteps the problem that a single final number rarely reproduces exactly.

**Block the answer path.** Network access runs through an allowlist covering package registries and OS mirrors. Anything else goes through a request tool where an adjudicator denies resources that would leak the answer (the original paper, its processed matrices, its supplementary files) while allowing software, documentation and reference genomes. This is cleaner than perturbing inputs, and it is essential here because several of our candidate papers have their ground truth sitting in a public supplementary XLSX.

One thing to handle differently. Claude Opus 5 ranked seventh in BixBench3 largely because it renamed specified output columns on three tasks, and excluding those three tasks it ranked second. That is instruction-following contaminating a measurement of analytical capability. We normalise identifiers before comparison and report format compliance as a separate axis.

## 4. Scope of the first release

Six to eight tasks, forty to sixty graded artifacts total, covering the eight analysis categories below. BixBench3 got 138 artifacts from 20 tasks, a median of five per task. We should land in the same range per task.

Target analysis categories:

1. BGC detection and annotation from genomes and metagenomes (antiSMASH, GECCO, DeepBGC)
2. Gene cluster family clustering and novelty scoring (BiG-SCAPE, BiG-SLiCE, cblaster, clinker)
3. MS/MS molecular networking and feature-based molecular networking (GNPS, MZmine)
4. Dereplication and structure annotation (SIRIUS, CSI:FingerID, NP Atlas)
5. Paired omics and metabologenomics link ranking (NPLinker, Paired Omics Data Platform)
6. NRPS and PKS domain parsing, adenylation domain substrate specificity
7. Cheminformatics on NP structures (RDKit scaffolds, NP-likeness, similarity)
8. Phylogenetics of biosynthetic enzymes and taxonomy versus chemistry concordance

A task may cover more than one category. Strain prioritisation and bioactivity tables ride along with categories 3 and 5.

## 5. Candidate source studies

Full screening results are in the research report. The six strongest, in build order:

| # | Study | Data | Categories | Depth |
|---|---|---|---|---|
| 1 | Crüsemann et al. 2017, *J Nat Prod* 80:588, 146 marine strains | MassIVE MSV000078836, PoDP project | networking, dereplication, prioritisation | 2 to 3 |
| 2 | Letzel/Jensen 2017, *Environ Microbiol* 19:3660, 119 *Salinispora* genomes, 176 BGCs | JGI IMG per-strain IDs | BGC detection, comparative genomics | 2 |
| 3 | Vertical inheritance in BGCs, 2021, *mBio* 12:e02700-21, 118 genomes, 3,041 BGCs, 305 GCFs | JGI IMG; OSF osf.io/entqc | GCF clustering, phylogenetics, concordance | 3+ |
| 4 | Castro-Falcón et al. 2025, *AEM* 92:e02171-25, *Salinispora* across spatial scales | MassIVE MSV000098188, NCBI genomes | paired omics, GCF clustering, networking | 3+ |
| 5 | Hjörleifsson Eldjárn et al. 2021, *PLoS Comput Biol* 17:e1008920, NPLinker | MSV000078836, MSV000085018, MSV000085038 | paired omics, full four-tool chain | 3+ |
| 6 | Tietz et al. 2017, *Nat Chem Biol* 13:470, RODEO lasso peptides, 1,419 BGCs | NCBI accessions; Supplementary Dataset 1 XLSX | RiPP mining, SSN, phylogenetics | 2 to 3 |

Three things to verify by hand before any of these is locked in:

- Crüsemann supplementary file format (XLSX versus PDF) and total MassIVE download size.
- The mBio vertical inheritance Table S1 and S2 file types, and what is on the OSF deposit.
- The Leão cyanobacteria accession discrepancy. The NPLinker paper cites MSV000085038; the NPClassScore paper cites MSV000084950 for what appears to be the same dataset.

Tier 2 candidates, if more coverage is needed: BiG-SLiCE (GigaScience 2020, scoped to a query subset), the fungal BGC atlas (PNAS 2021, 118:e2020230118, scoped to one clade), IIMN (Nat Commun 2021, 12:3832, which deposits MZmine parameters in Supplementary Table 4), the activated sludge long-read MAG study (mSystems 2022, which names antiSMASH and BiG-SCAPE v1.1.2 with cutoff 0.30), and the Micromonosporaceae study (Microb Genom 2024, 10:001167), which deposits both all BGC GenBank files and the exact command lines used.

## 6. Task bundle contract

Every task is a directory, not a question-answer pair.

```
tasks/<task_id>/
  task.yaml            # metadata, provenance, budget, depth labels
  prompt.md            # objective, method guidance, required outputs
  data/                # staged raw inputs, read-only at runtime
  data.manifest.json   # URL, accession, sha256, bytes per file
  reference/           # frozen database snapshots with versions
  truth/               # published artifacts, one file per graded artifact
  grading.yaml         # per-artifact metrics, keys, tolerances
  solution/            # our reference scripts that produced truth/
  environment/         # Dockerfile or conda lock, pinned
```

The prompt has four parts, following BixBench3: general instructions about the workspace and network policy, the research objective with a description of the staged data, method guidance naming tools and parameters, and a required outputs section giving the exact path and schema for every artifact.

Required outputs specify path, format, row key, expected row universe, and column names. Agents also write `outputs/METHODS.md` describing what they did, which feeds the failure-mode judge.

## 7. Grading

Four metrics, selected per artifact by what the artifact contains.

- **Row and column recovery**: F1 over identifiers, penalising both missing and extra entries. Used for gene lists, BGC inventories, feature IDs.
- **Numerical agreement**: Lin's concordance correlation coefficient, with a log transform where the quantity is log-distributed. CCC rather than Pearson, because CCC penalises a changed slope or offset, and a normalisation error that preserves rank order is a real error.
- **Categorical agreement**: macro F1 across label classes for aligned rows. Used for BGC class assignments, GCF membership, compound class.
- **Interval recovery**: overlap F1 under a minimum reciprocal overlap rule. Used for BGC coordinates on contigs.

The metrics applicable to one artifact are multiplied into a single artifact score. A missing, empty or unreadable file scores zero. An artifact passes at 0.80 or above.

The 0.80 threshold is not arbitrary and must be re-derived for this domain rather than inherited. BixBench3 obtained it by having domain experts rate 25 agent-produced artifacts from 1 to 5 on whether the biological meaning survived, finding that artifacts rated 4 or higher averaged 0.80. Run the same calibration with NP domain experts on our own artifacts. Reviewers will ask about this.

Identifier normalisation happens before comparison: case, whitespace, common accession formatting, and a declared alias map per artifact. Format compliance is recorded separately so that a model which produces correct science with wrong column names is visible as such.

## 8. Difficulty calibration

Difficulty is measured, not designed. Build more candidate artifacts than needed, run two or three strong models on each at least twice, and keep artifacts whose pass rate falls between roughly 0.1 and 0.8. Artifacts that every model passes and artifacts that no model passes both carry zero information about model differences.

Target a top-model task score around 0.3 to 0.5.

The levers that produce that range, in rough order of effect:

**Analysis depth.** BixBench3's clearest signal: mean artifact pass was 0.30 at depth 1, 0.44 at depth 2, and 0.24 at depth 3 or more. Most of our artifacts should sit at depth 3 or beyond.

**Cross-modal joins.** Genome to BGC to predicted product to observed MS/MS feature to known structure. Each join needs an identifier mapping that does not exist cleanly, and that is where agents silently inner-join away half the data.

**Real mess, left in.** Contigs that split a BGC across a boundary. Duplicate NP Atlas entries under different names. Stereochemistry collapsed by canonical SMILES. Isomers sharing a molecular formula. Strain names that disagree between the metadata sheet and the FASTA headers. The Crüsemann dataset is itself an example, described as 146 strains in the title, 124 in the MassIVE record, and 120 or 159 in downstream reprocessing. Do not clean these. The naive pipeline gives a confidently wrong answer and the careful one does not.

**Tool weight.** antiSMASH, BiG-SCAPE and a local networking run all have genuine install and runtime cost. Environment setup failure is one of BixBench3's ten tagged failure modes and was enriched among the worst runs.

**Scale that defeats eyeballing.** Files too large to read into context force scripted analysis. Runtime pressure forces parallelisation.

Two things that make tasks look hard without measuring anything: strict output formats and installation pain. Neither is a difficulty lever. Both are confounds.

## 9. Reproducibility constraints

This is where an NP benchmark differs most from a transcriptomics one, and it is the main technical risk.

**Service-dependent tools.** GNPS library matching runs against a continuously recurated spectral library, so a library hit list produced in 2021 is not reproducible today. SIRIUS and CSI:FingerID call a web service whose model version determines the output; SIRIUS 4, 5.8.x and 6.x disagree. Treat library match tables as soft artifacts, or freeze a library export inside the task bundle and grade against that. Network topology from a fixed .mgf with fixed parameters is reproducible and is the safer artifact.

**Version-sensitive tools.** antiSMASH detects different cluster type counts across major versions (58 in v5, 71 in v6, 101 in v8), so a BGC inventory is only reproducible against the exact version. BiG-SCAPE GCF counts depend on the raw distance cutoff (0.30 is the common choice) and on whether hybrids are included. Both must be written into the method guidance, not left to the agent.

**Database snapshots.** MIBiG, the antiSMASH database, the GNPS libraries and NP Atlas all move. Every reference database ships inside the task bundle with a recorded version and checksum. NP Atlas v2024_09 is a 676 MB download covering 36,545 compounds and is licensed CC BY-NC 4.0, which is compatible with research distribution but should be confirmed against whatever licence the benchmark itself carries.

**Manual curation in source papers.** Manual BGC boundary verification and manual structure assignment are common in this literature. Grade only the tool-generated artifacts.

**Record drift.** RiPP genome-neighbourhood studies note that some NCBI accessions no longer resolve. Frozen input bundles with checksums are not optional.

Tools that are safe to pin and run offline: antiSMASH, GECCO, DeepBGC, BiG-SCAPE (v1 and v2), BiG-SLiCE, RDKit, HMMER, Prodigal, all available through bioconda. The nf-core/funcscan pipeline already wires several of these together with pinned containers and is a useful reference for the environment definitions.

## 10. Compute budget

BixBench3 averaged 67 GB of raw input per task, 6.8 hours, 102 million tokens and about $43 per task attempt, with a worst case of 24 hours and $525. Across 13 models and 20 tasks that is a five-figure evaluation.

Our data is smaller. A few hundred bacterial genomes, a GNPS-scale MS/MS dataset, or a set of NP structure tables are orders of magnitude below 67 GB of RNA-seq. Keep every task well under 100 GB, ideally under 20 GB, and supply assembled genomes or MAGs rather than reads for anything that would otherwise require metagenome assembly. Doing a comparable benchmark an order of magnitude cheaper is worth stating in the paper.

Runtime target: 32 vCPUs, 128 GB RAM, no GPU, 24 hour cap. This matches BixBench3's environment and keeps results comparable.

## 11. Differentiator worth building

BixBench3 states in its own limitations that because each task supplies an explicit methodological plan, success means an agent can execute a specified pipeline, not that it can decide which analyses are worth pursuing.

A second tier where method guidance is withheld tests exactly what they cannot. The grading is harder, probably a rubric over the final claim or a tolerance on a headline number rather than artifact matching. Keep it small, maybe three tasks, and report it separately from the main score. For a copilot like Mycel, this is the capability that actually matters.

## 12. Relationship to Mycel

NPBench part one and part two together are the evaluation substrate for Mycel. Part two tasks exercise the same tool wrappers Mycel ships, so a task that Mycel's open code substrate cannot complete is a gap in the capability set. Keep the benchmark harness independent of Mycel so the benchmark remains credible when evaluating other systems, and so that a Mycel result is not self-graded.

## 13. Open questions

- Whether the expert calibration panel is lab-internal or external. External is more defensible and slower.
- Licence for the benchmark itself, and whether the NP Atlas CC BY-NC term constrains it.
- Whether to host staged data on Hugging Face (BixBench3's choice) or Zenodo.
- How many models to evaluate at release. Thirteen is thorough and expensive. Six is probably enough to show the depth and data-size effects.
