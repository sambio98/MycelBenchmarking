# Attribution

Every campaign in this release, the licence class it declares, and the
licence statements it carries. Assembled from each campaign's own files
rather than written here, so this file cannot drift from them: a
campaign's statement lives in `task.yaml` (`license_notice`) or in
`inputs/provenance.json` beside the data, and both are collected.

## arch-nrps-mibig-4_0-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'MIBiG 4.0 (https://mibig.secondarymetabolites.org/) and the MIBiG
    antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/),
    both licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries the annotation archive unchanged and the translations
    unmodified under a published extraction rule, with the source archive''s sha256
    recorded in inputs/provenance.json. CC BY 4.0 requires attribution and nothing
    further. The Pfam library is NOT redistributed: it is read from the pinned image
    at a declared path with a recorded sha256.

    '
  source: 'MIBiG 4.0 annotation records, and the CDS translations of every gene
    named by a curated module, taken from the antiSMASH-processed reference records

    '
```
- from `inputs/provenance.json`:

```json
{
  "annotations": {
    "attribution": "MIBiG 4.0 (https://mibig.secondarymetabolites.org/), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG"
  },
  "proteins": {
    "attribution": "MIBiG antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG antiSMASH-processed reference set",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2"
  }
}
```

## arch-pks-mibig-4_0-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'MIBiG 4.0 (https://mibig.secondarymetabolites.org/) and the MIBiG
    antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/),
    both licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries the annotation archive unchanged and the translations
    unmodified under a published extraction rule, with the source archive''s sha256
    recorded in inputs/provenance.json. CC BY 4.0 requires attribution and nothing
    further. The Pfam library is NOT redistributed: it is read from the pinned image
    at a declared path with a recorded sha256.

    '
  source: 'MIBiG 4.0 annotation records, and the CDS translations of every gene
    named by a curated module, taken from the antiSMASH-processed reference records

    '
```
- from `inputs/provenance.json`:

```json
{
  "annotations": {
    "attribution": "MIBiG 4.0 (https://mibig.secondarymetabolites.org/), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG"
  },
  "proteins": {
    "attribution": "MIBiG antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG antiSMASH-processed reference set",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2"
  }
}
```

## bgcdetect-scoelicolor-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'Streptomyces coelicolor A3(2) chromosome AL645882.2 from NCBI Nucleotide;
    curated cluster boundaries from MIBiG 4.0 (https://mibig.secondarymetabolites.org/),
    licensed CC BY 4.0.

    '
  license: 'GenBank records carry no copyright claim; MIBiG 4.0 is CC BY 4.0

    '
  redistribution: 'inputs/ carries the complete GenBank record (6.8 MB gzipped)
    and an extract of MIBiG''s locus table for this accession, one row per locus.
    How many rows that is, is R3''s `loci_total` and so is not stated here; an earlier
    version of this notice stated it. The antiSMASH detection rule files are NOT
    redistributed: they are AGPL-3.0-or-later package data, so the campaign reads
    them from the pinned image and records their paths, the provider, the version
    and a fingerprint instead.

    '
  source: 'NCBI Nucleotide AL645882.2 (Streptomyces coelicolor A3(2) complete chromosome)
    and the curated cluster boundaries MIBiG 4.0 maps to it

    '
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "AL645882.2",
    "attribution": "Streptomyces coelicolor A3(2) chromosome AL645882.2 from NCBI Nucleotide; curated cluster boundaries from MIBiG 4.0 (CC BY 4.0, https://mibig.secondarymetabolites.org/)",
    "license": "NCBI GenBank records carry no copyright claim; MIBiG 4.0 is CC BY 4.0",
    "sequence_source": "NCBI Nucleotide efetch, rettype=gbwithparts"
  }
}
```

## chemspace-mibig-4_0-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "license": "CC BY 4.0",
    "retrieved": "2026-10-02",
    "source": "MIBiG",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"
  }
}
```

## constraint-gc55-ecoli-rebh-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "Q8KHZ8",
    "license": "CC BY 4.0",
    "retrieved": "2026-10-01",
    "source": "UniProtKB"
  }
}
```

## constraint-gc62-ecoli-rebh-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "Q8KHZ8",
    "license": "CC BY 4.0",
    "retrieved": "2026-10-01",
    "source": "UniProtKB"
  }
}
```

## construct-ecoli-prna-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "P95480",
    "license": "CC BY 4.0",
    "retrieved": "2026-10-02",
    "sequence_version_source": "FASTA header SV=1",
    "source": "UniProtKB"
  }
}
```

## construct-ecoli-rebh-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "Q8KHZ8",
    "license": "CC BY 4.0",
    "retrieved": "2026-10-02",
    "sequence_version_source": "FASTA header SV=1",
    "source": "UniProtKB",
    "url": "https://rest.uniprot.org/uniprotkb/Q8KHZ8.json"
  }
}
```

## derep-massbank-2026_03-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'MassBank consortium, MassBank-data release 2026.03 (https://github.com/MassBank/MassBank-data),
    licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries spectra from the release unchanged apart from
    the declared reshaping: peaks sorted by m/z, m/z rounded to six decimals and
    intensity to four, and the fields the campaign reads kept. CC BY 4.0 requires
    attribution and nothing further. The selection rule, the census it produced
    and the per-file digests are in inputs/provenance.json.

    '
  source: MassBank-data, release tag 2026.03
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "attribution": "MassBank consortium, MassBank-data release 2026.03 (https://github.com/MassBank/MassBank-data), CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MassBank-data, release tag 2026.03, cloned shallow from https://github.com/MassBank/MassBank-data",
    "url": "https://github.com/MassBank/MassBank-data"
  }
}
```

## ecaudit-fmn-dh-1_1_3_15-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'UniProt Consortium, UniProtKB release 2026_03 (https://www.uniprot.org/),
    licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries the UniProt FASTA unchanged apart from concatenating
    the two sections; CC BY 4.0 requires attribution and nothing further. The Pfam
    library the audit scans against is NOT redistributed: it is read from the pinned
    image at the declared path, and reference/audit_rules.json records its provider,
    version and sha256 instead.

    '
  source: UniProtKB release 2026_03, all entries annotated EC 1.1.3.15
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "attribution": "UniProt Consortium, UniProtKB release 2026_03 (https://www.uniprot.org/), CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "UniProtKB REST"
  }
}
```

## frameshift-scoelicolor-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'Sequence from INSDC accession AL645882.2 (https://www.ebi.ac.uk/ena/browser/view/AL645882).

    '
  license: public domain / no restrictions on use (INSDC)
  redistribution: 'inputs/contig.fna.gz carries the record''s bases unchanged under
    a published extraction rule. INSDC sequence data carries no redistribution restriction.
    No annotation is shipped, and the mutated contigs are not shipped either: they
    follow from the declared rule.

    '
  source: 'The INSDC record for the Streptomyces coelicolor A3(2) chromosome, AL645882.2,
    reduced to its bases

    '
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "AL645882.2",
    "license": "public domain / no restrictions on use (INSDC)",
    "source": "EMBL-EBI/GenBank record for the Streptomyces coelicolor A3(2) chromosome",
    "url": "https://www.ebi.ac.uk/ena/browser/view/AL645882"
  }
}
```

## gcfcutoff-mibig-32-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'MIBiG antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/),
    licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/clusters/ carries the selected GenBank files unchanged,
    and inputs/designed_families.json records the groups and the rule that selected
    them. CC BY 4.0 requires attribution and nothing further. The Pfam library is
    NOT redistributed: it is read from the pinned image at a declared path with
    a recorded sha256.

    '
  source: 'The MIBiG antiSMASH-processed reference set, restricted to a declared
    selection of clusters, plus the chemistry-derived groups for them

    '
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "attribution": "MIBiG antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG, antiSMASH-processed reference set",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2"
  }
}
```

## genecall-scoelicolor-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'Sequence from INSDC accession AL645882.2 (https://www.ebi.ac.uk/ena/browser/view/AL645882).

    '
  license: public domain / no restrictions on use (INSDC)
  redistribution: 'inputs/contig.fna.gz carries the record''s bases unchanged under
    a published extraction rule, with the source archive''s digest in inputs/provenance.json.
    INSDC sequence data carries no redistribution restriction. No annotation is
    shipped.

    '
  source: 'The INSDC record for the Streptomyces coelicolor A3(2) chromosome, AL645882.2,
    reduced to its bases

    '
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "AL645882.2",
    "license": "public domain / no restrictions on use (INSDC)",
    "source": "EMBL-EBI/GenBank record for the Streptomyces coelicolor A3(2) chromosome",
    "url": "https://www.ebi.ac.uk/ena/browser/view/AL645882"
  }
}
```

## homology-twilight-ipr036188-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  redistribution: 'inputs/ carries the corpus sequences and the InterPro assignment
    table, both under licences that permit redistribution with attribution. The
    searched database is not shipped: it is built from the shipped corpus by the
    declared invocation. The mutated queries are not shipped either -- they follow
    from the declared shuffle rule.

    '
  sources:
  - attribution: 'UniProt Consortium, UniProtKB release 2026_03 (https://www.uniprot.org/).
      Sequences are redistributed unchanged under a published selection rule.

      '
    license: CC BY 4.0
    name: UniProtKB reviewed (Swiss-Prot), release 2026_03
  - attribution: 'InterPro release 110.0 (https://www.ebi.ac.uk/interpro/). Entry
      types and per-protein assignments are redistributed as a table.

      '
    license: CC0 1.0
    name: InterPro 110.0
```
- from `inputs/provenance.json`:

```json
{
  "corpus": {
    "attribution": "UniProt Consortium, UniProtKB release 2026_03 (https://www.uniprot.org/).",
    "license": "CC BY 4.0 (UniProt)",
    "source": "UniProtKB reviewed (Swiss-Prot), release 2026_03"
  },
  "labels": {
    "attribution": "InterPro release 110.0 (https://www.ebi.ac.uk/interpro/).",
    "license": "CC0 1.0 (InterPro)",
    "source": "InterPro 110.0, via the InterPro REST API for the entry types and UniProtKB 2026_03 cross-references for the per-protein assignments"
  }
}
```

## kinetics-brenda-ec1_1-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "license": "CC BY 4.0",
    "license_statement": "The BRENDA download page states: \"All copyrightable parts of BRENDA are licensed under Creative Commons Attribution License 4.0 (CC BY 4.0).\"",
    "retrieved": "2026-10-06",
    "source": "BRENDA"
  }
}
```

## massbalance-nrps-malleobactin-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "BGC0000386",
    "license": "CC BY 4.0",
    "retrieved": "2026-10-01",
    "source": "MIBiG",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"
  }
}
```

## massbalance-nrps-sevadicin-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "accession": "BGC0000426",
    "license": "CC BY 4.0",
    "retrieved": "2026-10-01",
    "source": "MIBiG",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"
  }
}
```

## mibig-diff-3_1-to-4_0-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "license": "CC BY 4.0",
    "retrieved": "2026-10-02",
    "source": "MIBiG",
    "url": "https://dl.secondarymetabolites.org/mibig/"
  }
}
```

## mibig-gene-function-evidence-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  ".": {
    "license": "CC BY 4.0",
    "retrieved": "2026-10-02",
    "source": "MIBiG",
    "url": "https://dl.secondarymetabolites.org/mibig/"
  }
}
```

## network-massbank-2026_03-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'MassBank consortium, MassBank-data release 2026.03 (https://github.com/MassBank/MassBank-data),
    licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries spectra from the release unchanged apart from
    the declared reshaping: peaks sorted by m/z, m/z rounded to six decimals and
    intensity to four, and the fields each side of the split may carry. CC BY 4.0
    requires attribution and nothing further. The two selection rules, the census
    they produced and the per-file digests are in inputs/provenance.json.

    '
  source: MassBank-data, release tag 2026.03
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "attribution": "MassBank consortium, MassBank-data release 2026.03 (https://github.com/MassBank/MassBank-data), CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MassBank-data, release tag 2026.03, cloned shallow from https://github.com/MassBank/MassBank-data",
    "url": "https://github.com/MassBank/MassBank-data"
  }
}
```

## pocket-rebh-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  "annotation_source": {
    "accession": "Q8KHZ8",
    "database": "UniProtKB",
    "license": "CC BY 4.0",
    "url": "https://rest.uniprot.org/uniprotkb/Q8KHZ8.json"
  },
  "structure_sources[0]": {
    "database": "RCSB PDB",
    "license": "CC0 1.0",
    "url": "https://files.rcsb.org/download/2E4G.cif"
  }
}
```

## residues-prna-01

- declared licence class: `open`
- from `inputs/provenance.json`:

```json
{
  "annotation_source": {
    "accession": "P95480",
    "database": "UniProtKB",
    "license": "CC BY 4.0",
    "url": "https://rest.uniprot.org/uniprotkb/P95480.json"
  },
  "structure_sources[0]": {
    "database": "RCSB PDB",
    "license": "CC0 1.0",
    "url": "https://files.rcsb.org/download/2AQJ.cif"
  },
  "structure_sources[1]": {
    "database": "RCSB PDB",
    "license": "CC0 1.0",
    "url": "https://files.rcsb.org/download/2AR8.cif"
  }
}
```

## ripp-precursor-mibig-4_0-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'MIBiG 4.0 (https://mibig.secondarymetabolites.org/) and the MIBiG
    antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/),
    both licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries both archives with their contents unchanged,
    and inputs/provenance.json records the selection rule and the fingerprints.
    CC BY 4.0 requires attribution and nothing further. The Pfam library is NOT
    redistributed: it is read from the pinned image at a declared path with a recorded
    sha256.

    '
  source: 'MIBiG 4.0 annotation records, and the antiSMASH-processed reference records
    for every entry of the declared biosynthetic class

    '
```
- from `inputs/provenance.json`:

```json
{
  "annotations": {
    "attribution": "MIBiG 4.0 (https://mibig.secondarymetabolites.org/), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG"
  },
  "reference_records": {
    "attribution": "MIBiG antiSMASH-processed reference set (https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2), licensed CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "MIBiG antiSMASH-processed reference set",
    "url": "https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2"
  }
}
```

## selectivity-saureus-topoisomerase-01

- declared licence class: `share_alike`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'Data from ChEMBL (https://www.ebi.ac.uk/chembl/), EMBL-EBI, release
    ChEMBL_37 of 2026-05-01, licensed CC BY-SA 3.0.

    '
  license: CC BY-SA 3.0
  redistribution: 'inputs/ and gold/ of this campaign are an adaptation of a ShareAlike
    work, so anyone redistributing them must do so under CC BY-SA 3.0 with the attribution
    above. The obligation is local to this campaign: the benchmark is a collection
    of separable campaigns rather than an adaptation of any one of them, so the
    code, the grader and every other campaign are unaffected and this campaign can
    be removed without touching them.

    '
  source: ChEMBL 37 (EMBL-EBI), activities and assay records for two targets
```

## transfer-ec1_14-uniprot-01

- declared licence class: `open`
- from `task.yaml`:

```yaml
license_notice:
  attribution: 'UniProt Consortium, UniProtKB release 2026_03 (https://www.uniprot.org/),
    licensed CC BY 4.0.

    '
  license: CC BY 4.0
  redistribution: 'inputs/ carries the UniProt sequences and their EC annotations.
    The FASTA headers are reduced to bare accessions so the hit table needs no re-parsing;
    nothing else is altered. CC BY 4.0 requires attribution and nothing further.

    '
  source: UniProtKB release 2026_03, reviewed entries annotated under EC 1.14
```
- from `inputs/provenance.json`:

```json
{
  ".": {
    "attribution": "UniProt Consortium, UniProtKB release 2026_03 (https://www.uniprot.org/), CC BY 4.0",
    "license": "CC BY 4.0",
    "source": "UniProtKB REST"
  }
}
```

