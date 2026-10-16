"""Recompute the pinned codon tables from the genome records that declare them.

Usage:
  python -m npbench_c.templates.construct_design.prepare_tables \
      --record ecoli_bl21=<CP001509.3 GenBank> \
      --record streptomyces_coelicolor=<AL645882.2 GenBank or .gb.gz> \
      [--adopt-digests]

Run when the manifest changes, not at instantiation. The output is
constants/codon_tables.json -- which params.copy_constants then copies into every
construct-design campaign's reference/ -- plus audit/codon_table_provenance.json,
which carries the per-codon counts and the gene list and is NOT shipped in a
sandbox.

Why this script exists. Gold does not depend on the weights being the single
correct measurement of codon usage: it depends on their being frozen, hashed and
published, so any party recomputes the same gold. But an auditor will ask where
a weight came from, and "representative values for highly-expressed genes" is not
an answer anyone can check. So the weights are derived here, by a declared rule
over a pinned record, and the rule's own output is pinned too.

What is pinned, and what is merely recorded. The enforced pins are the DERIVED
ones: the number of reference genes, and the digest of the reference set itself
(each gene's location, locus tag and length). The record's file digest is
recorded but not enforced, because the same accession re-fetched from a different
INSDC mirror is a different byte stream and the same biology; a refusal on the
transport encoding would fire on a correct re-run and teach the next person to
pass a flag. A refusal on the reference set fires only when the record really did
change.

The genetic code is imported from core, not restated: one table, one place.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pathlib
import re
import sys
from typing import Mapping, Sequence

from .core import AA_CODONS, CODON_TO_AA, CONSTANTS

AUDIT = pathlib.Path(__file__).resolve().parents[4] / "audit"
SOURCES = CONSTANTS / "codon_table_sources.json"
TABLES = CONSTANTS / "codon_tables.json"
PROVENANCE = "codon_table_provenance.json"

# A reference set smaller than this is not a codon-usage measurement, it is a
# handful of genes. The classic CAI sets are 25-40 genes; a bacterial genome
# annotates about 55 ribosomal proteins, so a count this low means the product
# pattern stopped matching the record's annotation style.
MIN_GENES = 40

# A family whose most-used codon is seen fewer times than this has its whole
# row decided by a handful of observations.
MIN_FAMILY_TOP = 20

# Unobserved SENSE codons, which the pseudocount rescues. A record whose
# reference set misses more than this many is too small or too skewed for the
# ratio to mean anything at the bottom of the table.
MAX_UNOBSERVED_SENSE = 12

# Product-matching features dropped by the validity checks. A few are expected
# (a partial location, a selenocysteine readthrough); a tenth of the set means
# the reader and the record disagree about something.
MAX_DROPPED_FRACTION = 0.10

ROUND = 6


class TableError(RuntimeError):
    """The tables cannot be derived from what the manifest declares. Always a
    provenance bug: a table whose rule did not run does not ship."""


# ------------------------------------------------------------------- the record


def _open(path: pathlib.Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def read_record(path: pathlib.Path) -> tuple[str, list[dict]]:
    """The ORIGIN bases and the feature table of the single record in `path`.

    Qualifier continuation is decided by quote state, not by a leading slash: a
    /note value can itself hold a line that starts with '/', and a reader that
    treats that as a new qualifier loses the location it was reading.
    """
    chunks: list[str] = []
    features: list[dict] = []
    in_features = in_origin = False
    current: dict | None = None
    key: str | None = None
    open_quote = False

    with _open(path) as handle:
        for line in handle:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN"):
                in_features, in_origin = False, True
                if current is not None:
                    features.append(current)
                    current = None
                continue
            if line.startswith("//"):
                break
            if in_origin:
                chunks.append("".join(line.split()[1:]))
                continue
            if not in_features:
                continue
            if line[:21].strip() and not line.startswith(" " * 21):
                if current is not None:
                    features.append(current)
                current = {"kind": line[5:21].strip(),
                           "location": line[21:].strip(),
                           "qualifiers": {}}
                key, open_quote = "__location__", False
                continue
            if current is None:
                continue
            body = line[21:].strip()
            start = None if open_quote else re.match(
                r"/([A-Za-z_0-9]+)(?:=(.*))?$", body)
            if start is not None:
                key, value = start.group(1), start.group(2)
                if value is None:
                    current["qualifiers"].setdefault(key, []).append("")
                    open_quote = False
                elif value.startswith('"'):
                    open_quote = not (len(value) > 1 and value.endswith('"'))
                    current["qualifiers"].setdefault(key, []).append(
                        value[1:] if open_quote else value[1:-1])
                else:
                    current["qualifiers"].setdefault(key, []).append(value)
                    open_quote = False
            elif key == "__location__":
                current["location"] += body
            elif key is not None:
                if body.endswith('"'):
                    body, open_quote = body[:-1], False
                # A translation wraps mid-sequence and is joined without a
                # separator; free text is joined with a space, or two words run
                # together across the line break.
                current["qualifiers"][key][-1] += (
                    "" if key == "translation" else " ") + body

    if current is not None:
        features.append(current)
    if not chunks:
        raise TableError(f"{path} carries no ORIGIN bases to read codons from")
    return "".join(chunks).upper(), features


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def simple_span(location: str) -> tuple[int, int, int] | None:
    """(start, end, strand) for one contiguous, complete range; else None.

    A join, a partial bound (< or >), or a remote accession returns None and the
    feature is dropped: this script counts codons, and a feature whose frame is
    spliced or whose end is unknown has no codons to count.
    """
    text = location.strip()
    strand = 1
    inner = re.fullmatch(r"complement\((.*)\)", text)
    if inner is not None:
        strand, text = -1, inner.group(1)
    span = re.fullmatch(r"(\d+)\.\.(\d+)", text)
    if span is None:
        return None
    return int(span.group(1)) - 1, int(span.group(2)), strand


def extract(sequence: str, span: tuple[int, int, int]) -> str:
    start, end, strand = span
    bases = sequence[start:end]
    return bases.translate(_COMPLEMENT)[::-1] if strand < 0 else bases


def translate(nt: str) -> str:
    return "".join(CODON_TO_AA[nt[i:i + 3]] for i in range(0, len(nt), 3))


# ------------------------------------------------------- the reference gene set


def reference_set(sequence: str, features: Sequence[Mapping],
                  pattern: re.Pattern) -> tuple[list[dict], dict[str, int]]:
    """The genes the declared rule selects, and why each candidate was dropped.

    The census is returned rather than logged: the dropped count is a refusal
    threshold, and a threshold needs the number it is checked against.
    """
    genes: list[dict] = []
    census = {"product_match": 0, "pseudo": 0, "spliced_or_partial": 0,
              "codon_start": 0, "transl_table": 0, "not_whole_codons": 0,
              "stop_structure": 0, "translation_mismatch": 0, "no_translation": 0}

    for feature in features:
        if feature["kind"] != "CDS":
            continue
        qualifiers = feature["qualifiers"]
        product = (qualifiers.get("product") or [""])[0].strip()
        if not pattern.match(product):
            continue
        census["product_match"] += 1
        if "pseudo" in qualifiers or "pseudogene" in qualifiers:
            census["pseudo"] += 1
            continue
        span = simple_span(feature["location"])
        if span is None:
            census["spliced_or_partial"] += 1
            continue
        if (qualifiers.get("codon_start") or ["1"])[0].strip() != "1":
            census["codon_start"] += 1
            continue
        if (qualifiers.get("transl_table") or ["11"])[0].strip() != "11":
            census["transl_table"] += 1
            continue
        nt = extract(sequence, span)
        if len(nt) % 3 or set(nt) - set("ACGT"):
            census["not_whole_codons"] += 1
            continue
        protein = translate(nt)
        if not protein.endswith("*") or "*" in protein[:-1]:
            census["stop_structure"] += 1
            continue
        declared = re.sub(r"\s+", "", (qualifiers.get("translation") or [""])[0])
        if not declared:
            census["no_translation"] += 1
            continue
        # The record's own translation is the check that the frame is right. The
        # initiator is compared as methionine whatever its codon, which is how
        # the record writes it.
        if "M" + protein[1:-1] != declared:
            census["translation_mismatch"] += 1
            continue
        genes.append({
            "locus_tag": (qualifiers.get("locus_tag") or [""])[0].strip(),
            "gene": (qualifiers.get("gene") or [""])[0].strip(),
            "product": product,
            "location": feature["location"],
            "length_nt": len(nt),
            "bases": nt,
        })

    genes.sort(key=lambda gene: (gene["location"], gene["locus_tag"]))
    return genes, census


def gene_set_digest(genes: Sequence[Mapping]) -> str:
    """A digest of the reference SET, not of the file it came from.

    Location, locus tag and length: enough that a different record, a different
    annotation release, or a pattern that started matching something new changes
    the digest, and a re-fetch of the same accession from another mirror does
    not.
    """
    return hashlib.sha256("\n".join(
        f"{gene['location']}\t{gene['locus_tag']}\t{gene['length_nt']}"
        for gene in genes).encode()).hexdigest()


def codon_counts(genes: Sequence[Mapping]) -> dict[str, int]:
    """Every codon of every reference gene, the terminal stop included.

    CAI itself excludes stops, so the stop row never reaches a score; it is
    counted anyway, because a table that reports a weight of 1.0 for all three
    stops is reporting the pseudocount, not the record.
    """
    counts = {codon: 0 for codon in CODON_TO_AA}
    for gene in genes:
        nt = gene["bases"]
        for i in range(0, len(nt), 3):
            counts[nt[i:i + 3]] += 1
    return counts


def relative_adaptiveness(counts: Mapping[str, int],
                          pseudocount: float) -> tuple[dict[str, float], list[str]]:
    """w = count / max(count in the synonymous family), zeros rescued first."""
    weights: dict[str, float] = {}
    unobserved: list[str] = []
    for amino_acid, codons in AA_CODONS.items():
        adjusted = {}
        for codon in codons:
            if counts[codon] == 0:
                unobserved.append(codon)
                adjusted[codon] = float(pseudocount)
            else:
                adjusted[codon] = float(counts[codon])
        top = max(adjusted.values())
        for codon in codons:
            weights[codon] = round(adjusted[codon] / top, ROUND)
    return weights, sorted(unobserved)


# ------------------------------------------------------------------ the refusals


def _assert_usable(host: str, genes: Sequence[Mapping], census: Mapping[str, int],
                   counts: Mapping[str, int], weights: Mapping[str, float],
                   unobserved: Sequence[str]) -> None:
    if len(genes) < MIN_GENES:
        raise TableError(
            f"{host}: the declared rule selected {len(genes)} reference genes, "
            f"fewer than {MIN_GENES}. A bacterial genome annotates about 55 "
            "ribosomal proteins, so this record names them in a style the "
            "product pattern does not match; fix the pattern in the manifest "
            "rather than deriving a table from a handful of genes")

    dropped = census["product_match"] - len(genes)
    if census["product_match"] and dropped > MAX_DROPPED_FRACTION * census["product_match"]:
        raise TableError(
            f"{host}: {dropped} of {census['product_match']} product-matching "
            f"features failed the validity checks ({dict(census)}), more than "
            f"{MAX_DROPPED_FRACTION:.0%}. The reader and the record disagree "
            "about something; the weights would be computed from whichever "
            "genes happened to survive")

    sense_unobserved = [c for c in unobserved if CODON_TO_AA[c] != "*"]
    if len(sense_unobserved) > MAX_UNOBSERVED_SENSE:
        raise TableError(
            f"{host}: {len(sense_unobserved)} sense codons are unobserved in the "
            f"reference set ({sense_unobserved}), more than {MAX_UNOBSERVED_SENSE}. "
            "That many rows would be decided by the pseudocount rather than by "
            "the record")

    thin = sorted(
        (max(counts[c] for c in codons), amino_acid)
        for amino_acid, codons in AA_CODONS.items() if len(codons) > 1)
    if thin and thin[0][0] < MIN_FAMILY_TOP:
        count, amino_acid = thin[0]
        raise TableError(
            f"{host}: the most-used codon of the {amino_acid} family is seen only "
            f"{count} times, fewer than {MIN_FAMILY_TOP}. Every weight in that "
            "family is then a ratio of single-digit counts")

    for amino_acid, codons in AA_CODONS.items():
        tops = [c for c in codons if weights[c] == 1.0]
        if len(tops) != 1:
            raise TableError(
                f"{host}: the {amino_acid} family has {len(tops)} codons at weight "
                f"1.0 ({tops}). Relative adaptiveness is a ratio to the family's "
                "most-used codon, so exactly one codon per family is 1.0; a tie "
                "makes the table's scale ambiguous and is resolved in the "
                "manifest, not silently")
    if any(weight <= 0.0 for weight in weights.values()):
        raise TableError(
            f"{host}: a weight is not positive, so CAI's logarithm is undefined. "
            "The pseudocount exists to prevent exactly this")


def _assert_window_holds(host: str, spec: Mapping, sequence: str) -> float:
    """The declared design window is not a measurement, but it must contain one.

    A synthesis window that excludes the host's own genome is a typo, and it
    would make every campaign on that host unsatisfiable in a way the oracle
    discovers late.
    """
    gc = sum(1 for base in sequence if base in "GC") / len(sequence)
    low, high = spec["design_window"]["gc_global"]
    if not low <= gc <= high:
        raise TableError(
            f"{host}: the pinned record's GC fraction is {gc:.5f}, outside the "
            f"declared gc_global window [{low}, {high}]. Either the window is a "
            "typo or the record is not this host's genome")
    return gc


# --------------------------------------------------------------------- the build


def _file_digest(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def build_tables(manifest: Mapping, records: Mapping[str, pathlib.Path],
                 adopt_digests: bool = False) -> tuple[dict, dict]:
    declared = set(manifest["hosts"])
    given = set(records)
    if declared != given:
        raise TableError(
            f"the manifest declares hosts {sorted(declared)} and --record was "
            f"given {sorted(given)}. Every declared host needs its record, and a "
            "record for a host the manifest does not declare has no window, "
            "licence or accession to publish")

    rule = manifest["reference_set"]
    pattern = re.compile(rule["product_pattern"], re.IGNORECASE)
    pseudocount = float(rule["pseudocount"])
    if not 0.0 < pseudocount < 1.0:
        raise TableError(
            f"the manifest's pseudocount is {pseudocount}; it stands in for a "
            "codon seen zero times and must be between 0 and 1, or an "
            "unobserved codon outranks an observed one")

    accessions = [spec["accession"] for spec in manifest["hosts"].values()]
    if len(set(accessions)) != len(accessions):
        raise TableError(
            f"two hosts are pinned to the same accession {accessions}. Their "
            "tables would be identical and the campaign's host swap -- the "
            "whole of R4 -- would compare a table against itself")

    tables: dict = {
        "description": (
            "Pinned codon relative-adaptiveness weights and host design windows. "
            "Derived by constants/prepare_tables.py from the records named in "
            "provenance; see audit/codon_table_provenance.json for the per-codon "
            "counts and the gene list. Campaign data: gold does not depend on "
            "these being the single correct measurement of codon usage, only on "
            "their being frozen, hashed and published so any party recomputes "
            "the same gold."),
        "hosts": {},
        "provenance": {
            "generated_by": "npbench_c.templates.construct_design.prepare_tables",
            "reference_set_rule": rule["rule"],
            "pseudocount": pseudocount,
            "product_pattern": rule["product_pattern"],
            "why_ribosomal": rule["why_ribosomal"],
            "hosts": {},
        },
    }
    audit: dict = {
        "_what_this_is": (
            "The arithmetic behind constants/codon_tables.json: the reference "
            "gene list and the per-codon counts, so the weights can be checked "
            "without re-downloading a genome. Audit-side on purpose -- it is not "
            "copied into any campaign's reference/ and never reaches a sandbox."),
        "generated_by": "npbench_c.templates.construct_design.prepare_tables",
        "hosts": {},
    }

    for host in sorted(manifest["hosts"]):
        spec = manifest["hosts"][host]
        path = pathlib.Path(records[host]).resolve()
        sequence, features = read_record(path)
        genes, census = reference_set(sequence, features, pattern)
        counts = codon_counts(genes)
        weights, unobserved = relative_adaptiveness(counts, pseudocount)
        _assert_usable(host, genes, census, counts, weights, unobserved)
        gc = _assert_window_holds(host, spec, sequence)
        digest = gene_set_digest(genes)

        if adopt_digests:
            spec["genes"] = len(genes)
            spec["gene_set_sha256"] = digest
            spec["record_sha256_when_derived"] = _file_digest(path)
        else:
            for field, found in (("genes", len(genes)),
                                 ("gene_set_sha256", digest)):
                if field not in spec:
                    raise TableError(
                        f"{host}: the manifest declares no {field}, so there is "
                        "nothing to check this run against. Run once with "
                        "--adopt-digests to record what the rule selects, read "
                        "the diff, and commit it")
                if spec[field] != found:
                    raise TableError(
                        f"{host}: the reference set does not match the manifest: "
                        f"{field} is {found!r}, the manifest pins {spec[field]!r}. "
                        "The record, its annotation release, or the product "
                        "pattern changed. Decide what the table should be, then "
                        "re-adopt the digests deliberately")

        tables["hosts"][host] = {
            "gc_global": list(spec["design_window"]["gc_global"]),
            "gc_window": list(spec["design_window"]["gc_window"]),
            "label": spec["label"],
            "weights": {codon: weights[codon] for codon in sorted(weights)},
        }
        tables["provenance"]["hosts"][host] = {
            "accession": spec["accession"],
            "codons_counted": sum(counts.values()),
            "design_window_basis": spec["design_window_basis"],
            "gene_set_sha256": digest,
            "genes": len(genes),
            "license": spec["license"],
            "organism": spec["organism"],
            "record_gc_fraction": round(gc, ROUND),
            "record_sha256_when_derived": spec.get(
                "record_sha256_when_derived", _file_digest(path)),
            "source": spec["source"],
            "unobserved_codons": unobserved,
            "url": spec["url"],
        }
        audit["hosts"][host] = {
            "accession": spec["accession"],
            "census_of_dropped_candidates": census,
            "codon_counts": {codon: counts[codon] for codon in sorted(counts)},
            "gene_set_sha256": digest,
            "genes": [{key: gene[key] for key in
                       ("locus_tag", "gene", "product", "location", "length_nt")}
                      for gene in genes],
            "weights": {codon: weights[codon] for codon in sorted(weights)},
        }

    return tables, audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", action="append", default=[], metavar="HOST=PATH",
                        help="the GenBank record for one declared host")
    parser.add_argument("--adopt-digests", action="store_true",
                        help="record what the rule selects into the manifest "
                             "instead of checking against it")
    args = parser.parse_args(argv)

    records: dict[str, pathlib.Path] = {}
    for item in args.record:
        if "=" not in item:
            raise TableError(f"--record {item!r} is not HOST=PATH")
        host, _, path = item.partition("=")
        if host in records:
            raise TableError(f"--record names {host} twice")
        records[host] = pathlib.Path(path)
        if not records[host].exists():
            raise TableError(f"{host}: no record at {path}")

    manifest = json.loads(SOURCES.read_text())
    tables, audit = build_tables(manifest, records, adopt_digests=args.adopt_digests)

    if args.adopt_digests:
        SOURCES.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    TABLES.write_text(json.dumps(tables, indent=2, sort_keys=True) + "\n")
    AUDIT.mkdir(parents=True, exist_ok=True)
    (AUDIT / PROVENANCE).write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")

    for host in sorted(tables["hosts"]):
        block = tables["provenance"]["hosts"][host]
        print(f"  {host:<24} {block['accession']:<12} "
              f"{block['genes']} genes  {block['codons_counted']} codons  "
              f"GC {block['record_gc_fraction']}  "
              f"unobserved {len(block['unobserved_codons'])}")
        print(f"  {'':<24} gene_set_sha256 {block['gene_set_sha256']}")
    print(f"  wrote {TABLES}")
    print(f"  wrote {AUDIT / PROVENANCE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
