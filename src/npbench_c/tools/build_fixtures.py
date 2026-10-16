"""Build the thread-invariance fixtures once and pin them.

A determinism suite is only as good as its inputs: if the fixture changes between
runs, "identical output" means nothing. So the fixtures are committed files with
recorded provenance, and this script regenerates them byte-identically.

Two fixtures:

  halogenases.faa   24 reviewed UniProt proteins carrying Pfam PF04820, fetched
                    once. Real homologs spanning the similarity range, including
                    RebH and PrnA (already in the benchmark) and one fragment --
                    which is the kind of record that makes an aligner's tie-breaks
                    visible.
  halogenases.afa   the MAFFT alignment of that set, frozen.
  synthetic.fna     a nucleotide contig built by reverse-translating the set with
                    a declared codon cycle and joining with declared spacers.
  bigscape_input/   32 antiSMASH-processed MIBiG reference clusters, selected by
                    a declared rule so the set has real family structure rather
                    than 32 unrelated singletons.
  bgc_triplet.embl  three complete deposited records, each carrying one
                    characterised cluster of a different class, concatenated.
                    antiSMASH needs annotated input, and whole records rather
                    than slices mean the provenance is three accessions with no
                    coordinates to get wrong.

The alignment fixture is produced by a tool the suite also tests, which is worth
being explicit about: it is not circular, because once written the alignment is a
frozen file with a hash. If MAFFT later fails invariance, the HMMER results
computed from this file are unaffected -- the file does not change.

The nucleotide fixture is synthetic on purpose. Prodigal needs a contig with open
reading frames in it, and no nucleotide source in this project ships one;
reverse-translation from the protein fixture gives a contig whose gene content is
known exactly and whose bytes are reproducible from the repo alone, with no
further download.

Usage: python -m npbench_c.tools.build_fixtures [--mafft <path>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys

FIXTURES = pathlib.Path(__file__).resolve().parent / "fixtures"
PROTEINS = "halogenases.faa"
ALIGNMENT = "halogenases.afa"
CONTIG = "synthetic.fna"
CLUSTER = "bgc_triplet.embl"
BGC_SET = "bigscape_input"
PROVENANCE = "provenance.json"

#: The standard genetic code, codons sorted so the cycle below is reproducible.
CODONS: dict[str, tuple[str, ...]] = {
    "A": ("GCA", "GCC", "GCG", "GCT"), "C": ("TGC", "TGT"),
    "D": ("GAC", "GAT"), "E": ("GAA", "GAG"), "F": ("TTC", "TTT"),
    "G": ("GGA", "GGC", "GGG", "GGT"), "H": ("CAC", "CAT"),
    "I": ("ATA", "ATC", "ATT"), "K": ("AAA", "AAG"),
    "L": ("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"), "M": ("ATG",),
    "N": ("AAC", "AAT"), "P": ("CCA", "CCC", "CCG", "CCT"),
    "Q": ("CAA", "CAG"), "R": ("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"),
    "S": ("AGC", "AGT", "TCA", "TCC", "TCG", "TCT"),
    "T": ("ACA", "ACC", "ACG", "ACT"), "V": ("GTA", "GTC", "GTG", "GTT"),
    "W": ("TGG",), "Y": ("TAC", "TAT"),
}
STOPS = ("TAA", "TAG", "TGA")

#: Declared spacer: a fixed 60-nt intergenic sequence, repeated between ORFs so
#: Prodigal sees gene boundaries rather than one fused reading frame.
SPACER = "ATGCCTAGGATCCGATTACGAATTCGGTACCAAGCTTGCATGCCTGCAGGTCGACTCTAGA"
LINE_WIDTH = 60


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def reverse_translate(protein: str) -> str:
    """Declared rule: residue i takes codon i % len(synonyms) for its amino acid.

    Cycling rather than always taking the first codon keeps the base composition
    realistic enough for a gene caller to work on, and is still a pure function
    of the sequence. Residues outside the standard twenty (X, B, Z, U, O) are
    dropped, which affects only the fragment record.
    """
    out = []
    for index, residue in enumerate(protein.upper()):
        synonyms = CODONS.get(residue)
        if synonyms is None:
            continue
        out.append(synonyms[index % len(synonyms)])
    return "".join(out)


def build_contig(records: list[tuple[str, str]]) -> str:
    """One contig: spacer, then each ORF with a stop codon, each followed by a spacer."""
    pieces = [SPACER]
    for index, (_, protein) in enumerate(records):
        pieces.append(reverse_translate(protein) + STOPS[index % len(STOPS)])
        pieces.append(SPACER)
    return "".join(pieces)


def wrap(sequence: str) -> str:
    return "\n".join(sequence[i:i + LINE_WIDTH]
                     for i in range(0, len(sequence), LINE_WIDTH))


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mafft", default="mafft",
                    help="path to mafft; the alignment fixture needs it once")
    ap.add_argument("--skip-alignment", action="store_true",
                    help="leave halogenases.afa as committed")
    args = ap.parse_args(argv)

    if not (FIXTURES / BGC_SET).is_dir():
        print(f"{FIXTURES / BGC_SET} is missing; it is a committed fixture, "
              "not derived", file=sys.stderr)
        return 2
    for required in (PROTEINS, CLUSTER):
        if not (FIXTURES / required).is_file():
            print(f"{FIXTURES / required} is missing; it is a committed fixture, "
                  "not derived", file=sys.stderr)
            return 2
    proteins = FIXTURES / PROTEINS
    if not proteins.is_file():
        print(f"{proteins} is missing; it is a committed fixture, not derived",
              file=sys.stderr)
        return 2
    records = read_fasta(proteins)

    contig = build_contig(records)
    (FIXTURES / CONTIG).write_text(
        f">npbench_c_synthetic_contig reverse-translated from {PROTEINS}, "
        f"{len(records)} ORFs, declared codon cycle and spacer\n"
        + wrap(contig) + "\n")

    if not args.skip_alignment:
        # --thread 1 and a fixed strategy: the fixture must not depend on the
        # machine that happened to build it.
        result = subprocess.run(
            [args.mafft, "--thread", "1", "--anysymbol", "--retree", "2",
             "--maxiterate", "0", "--quiet", str(proteins)],
            capture_output=True, text=True, check=True)
        (FIXTURES / ALIGNMENT).write_text(result.stdout)

    aligned = read_fasta(FIXTURES / ALIGNMENT)
    (FIXTURES / PROVENANCE).write_text(json.dumps({
        "description": "Pinned inputs for the thread-invariance suite. A "
                       "determinism claim is only as good as its fixture, so "
                       "these are committed files with hashes, not downloads.",
        "built_on": "2026-10-02",
        PROTEINS: {
            "source": "UniProtKB",
            "url": "https://rest.uniprot.org/uniprotkb/stream"
                   "?query=xref:pfam-PF04820+AND+reviewed:true&format=fasta",
            "query": "xref:pfam-PF04820 AND reviewed:true",
            "license": "CC BY 4.0",
            "records": len(records),
            "residues": sum(len(s) for _, s in records),
            "note": "Flavin-dependent halogenases. Includes RebH (Q8KHZ8) and "
                    "PrnA (P95480), already in the benchmark, and one fragment "
                    "(P81540) whose partial length is deliberate: a short record "
                    "is where an aligner's and a search tool's tie-breaks show.",
            "sha256": sha256(proteins),
        },
        ALIGNMENT: {
            "derived_from": PROTEINS,
            "command": "mafft --thread 1 --anysymbol --retree 2 --maxiterate 0",
            "columns": len(aligned[0][1]) if aligned else 0,
            "note": "Produced by a tool the suite also tests. Not circular: once "
                    "written this is a frozen file, so a later MAFFT invariance "
                    "failure does not change the HMMER results computed from it.",
            "sha256": sha256(FIXTURES / ALIGNMENT),
        },
        CLUSTER: {
            "source": "ENA (INSDC)",
            "url_template": "https://www.ebi.ac.uk/ena/browser/api/embl/"
                            "{accession}?lineLimit=0",
            "assembled_by": "concatenating the three records in the order below; "
                            "EMBL flat files terminate each record with //, so "
                            "concatenation is the whole of it",
            "records": [
                {"accession": "AY312585.1", "length_bp": 65009, "cds": 57,
                 "organism": "Vibrio anguillarum 775",
                 "description": "plasmid pJM1, complete sequence",
                 "mibig": "BGC0002468, NRPS, knock-out evidence, anguibactin"},
                {"accession": "MN724926.1", "length_bp": 16516, "cds": 9,
                 "organism": "Kamptonema sp. PCC 6506",
                 "description": "landornamide biosynthetic gene cluster",
                 "mibig": "BGC0002330, ribosomal, knock-out evidence, "
                          "landornamide A"},
                {"accession": "MT799798.1", "length_bp": 16028, "cds": 15,
                 "organism": "Streptomyces sampsonii",
                 "description": "julichrome biosynthetic gene cluster",
                 "mibig": "BGC0002514, PKS, knock-out evidence, julichrome Q6.6"},
            ],
            "why_these_records": "Three classes rather than one, so the detection "
                                 "rules exercised are not all of a kind: the run "
                                 "yields NRP-metallophore plus NRPS, "
                                 "lanthipeptide-class-ii plus proteusin, and "
                                 "T2PKS. All three are COMPLETE deposited records "
                                 "whose MIBiG locus starts at position 1, so the "
                                 "provenance is three accessions with no "
                                 "coordinates to get wrong. All three are "
                                 "ANNOTATED, so antiSMASH's verdict is about "
                                 "antiSMASH rather than about a gene caller "
                                 "feeding it. All three are bacterial, so one "
                                 "--taxon setting covers the file.",
            "license": "INSDC records are freely available without restriction",
            "sha256": sha256(FIXTURES / CLUSTER),
        },
        BGC_SET: {
            "source": "MIBiG, antiSMASH-processed reference set",
            "url": "https://dl.secondarymetabolites.org/mibig/"
                   "mibig_antismash_4.0_gbk_as8b1.tar.bz2",
            "archive_sha256": "cf7e27514197983c814fbe2b6ade68e7"
                              "23311d10e2f12b2583447192ad6d3310",
            "license": "CC BY 4.0",
            "antismash_version_that_processed_them": "8.0 beta 1, as the archive "
                "name as8b1 states. Worth noticing rather than passing over: the "
                "image pins antiSMASH 8.0.4, so a campaign that compares its own "
                "antiSMASH 8.0.4 regions against this reference set is mixing two "
                "antiSMASH versions -- the comparison antiSMASH's own rules "
                "forbid. For a determinism fixture that is harmless, because both "
                "arms of every comparison use the same files; for a T-L3-2 "
                "campaign it is a constraint that has to be stated.",
            "selection_rule": "From the chemical-space campaign's connectivity "
                "groups: the four largest groups of 4 to 10 entries sharing an "
                "InChIKey connectivity block, after excluding groups whose MIBiG "
                "compound name is a class placeholder rather than a specific "
                "compound, plus one singleton entry per biosynthetic class as a "
                "negative control. Deterministic, and reproducible from tables "
                "already in this repository.",
            "excluded_placeholder_names": ["capsular polysaccharide",
                                           "lipopolysaccharide", "melanin",
                                           "carotenoid", "exopolysaccharide"],
            "excluded_placeholder_reason": "MIBiG gives these generic names a "
                "representative structure, so entries sharing one are not "
                "necessarily homologous clusters. Two different formulas appear "
                "under 'capsular polysaccharide' alone. Including them would have "
                "put apparent family structure in the fixture that the biology "
                "does not support.",
            "families": {
                "ectoine": ["BGC0000852", "BGC0000853", "BGC0000854",
                            "BGC0000855", "BGC0000856", "BGC0000857",
                            "BGC0000858", "BGC0000859", "BGC0000860",
                            "BGC0002052"],
                "ochratoxin A": ["BGC0001030", "BGC0002598", "BGC0002605",
                                 "BGC0002606", "BGC0002607", "BGC0002608",
                                 "BGC0002609"],
                "kanamycin": ["BGC0000702", "BGC0000703", "BGC0000704",
                              "BGC0000705", "BGC0000706"],
                "Aflatoxin B1": ["BGC0000006", "BGC0000007", "BGC0000008",
                                 "BGC0000009", "BGC0000011"],
            },
            "controls": ["BGC0000026", "BGC0000605", "BGC0001082", "BGC0001483",
                         "BGC0002197"],
            "note": "32 clusters over four classes. BiG-SCAPE recovers the "
                    "designed families at a 0.3 cutoff, which is what makes the "
                    "fixture a test of clustering rather than a test of whether "
                    "32 singletons stay 32 singletons.",
        },
        CONTIG: {
            "derived_from": PROTEINS,
            "rule": "Each protein is reverse-translated by taking, for residue i, "
                    "codon i modulo the number of synonyms for that amino acid, "
                    "from the standard genetic code with codons sorted "
                    "alphabetically. A stop codon is appended, cycling TAA, TAG, "
                    "TGA. ORFs are separated and flanked by a fixed 60-nt spacer. "
                    "Residues outside the standard twenty are dropped.",
            "length_nt": len(contig),
            "orfs": len(records),
            "note": "Synthetic on purpose: Prodigal needs a contig and no "
                    "nucleotide source in this project ships one. Gene content is "
                    "known exactly and the bytes are reproducible from the repo "
                    "with no download.",
            "sha256": sha256(FIXTURES / CONTIG),
        },
    }, indent=2, sort_keys=True) + "\n")

    print(f"{PROTEINS}: {len(records)} records, "
          f"{sum(len(s) for _, s in records)} residues")
    print(f"{ALIGNMENT}: {len(aligned)} records, "
          f"{len(aligned[0][1]) if aligned else 0} columns")
    print(f"{CONTIG}: {len(contig)} nt, {len(records)} ORFs")
    cluster = FIXTURES / CLUSTER
    if cluster.is_file():
        text = cluster.read_text()
        cds = sum(1 for line in text.splitlines() if line.startswith("FT   CDS"))
        ids = sum(1 for line in text.splitlines() if line.startswith("ID   "))
        print(f"{CLUSTER}: {cluster.stat().st_size // 1024} KB, "
              f"{ids} records, {cds} CDS")
    bgcs = FIXTURES / BGC_SET
    if bgcs.is_dir():
        files = sorted(bgcs.glob("*.gbk"))
        size = sum(f.stat().st_size for f in files)
        print(f"{BGC_SET}/: {len(files)} clusters, {size // 1024} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
