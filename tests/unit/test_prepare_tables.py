"""The codon tables are derived, and the derivation refuses rather than guesses.

Two things are tested here. First, that what is committed is what the rule
produces: regenerating from the two pinned records reproduces
`constants/codon_tables.json` and `audit/codon_table_provenance.json` byte for
byte, and the manifest's pinned reference-set digests still describe the records.
That is the whole provenance claim -- a table nobody can recompute is an
assertion with extra steps.

Second, every refusal. A derivation that silently produced a table from a handful
of genes, or from a record whose annotation style the product pattern stopped
matching, would be worse than the representative values it replaced: it would
look derived. So each `TableError` in `prepare_tables` is raised by a synthetic
record built to trip exactly it.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from npbench_c.templates.construct_design import prepare_tables as pt
from npbench_c.templates.construct_design.core import AA_CODONS

ROOT = pathlib.Path(__file__).resolve().parents[2]
SENSE = [codon for aa, codons in AA_CODONS.items() if aa != "*" for codon in codons]
# One designated most-used codon per family that has a choice. Stops are
# excluded: an extra stop codon inside a gene body is an internal stop, and the
# gene is then dropped as the rule says it should be.
TOPS = {aa: codons[0] for aa, codons in AA_CODONS.items()
        if len(codons) > 1 and aa != "*"}


# --------------------------------------------------------------- the real thing


def test_the_committed_tables_are_what_the_rule_produces():
    """The provenance claim, checked rather than described."""
    manifest = json.loads(pt.SOURCES.read_text())
    records = {host: ROOT / spec["repository_copy"]
               for host, spec in manifest["hosts"].items()}
    for host, path in records.items():
        if not path.is_file():
            pytest.skip(f"{host}: {path} not in this checkout")

    tables, audit = pt.build_tables(manifest, records)
    assert tables == json.loads(pt.TABLES.read_text())
    assert audit == json.loads((pt.AUDIT / pt.PROVENANCE).read_text())


def test_the_reference_sets_are_the_sizes_the_manifest_pins():
    manifest = json.loads(pt.SOURCES.read_text())
    assert sorted(manifest["hosts"]) == ["ecoli_bl21", "streptomyces_coelicolor"]
    for host, spec in manifest["hosts"].items():
        assert spec["genes"] >= pt.MIN_GENES, host
        assert len(spec["gene_set_sha256"]) == 64, host
        assert len(spec["record_sha256_when_derived"]) == 64, host


def test_every_weight_is_a_positive_ratio_with_one_top_per_family():
    tables = json.loads(pt.TABLES.read_text())
    for host, block in tables["hosts"].items():
        weights = block["weights"]
        assert sorted(weights) == sorted(pt.CODON_TO_AA), host
        assert all(0.0 < w <= 1.0 for w in weights.values()), host
        for amino_acid, codons in AA_CODONS.items():
            tops = [c for c in codons if weights[c] == 1.0]
            assert len(tops) == 1, (host, amino_acid, tops)


def test_the_provenance_block_names_a_record_per_host():
    tables = json.loads(pt.TABLES.read_text())
    for host, block in tables["provenance"]["hosts"].items():
        assert block["accession"]
        assert block["url"].startswith("https://")
        assert block["license"]
        assert block["genes"] >= pt.MIN_GENES
        assert block["codons_counted"] > 1000
        low, high = tables["hosts"][host]["gc_global"]
        assert low <= block["record_gc_fraction"] <= high


# ---------------------------------------------------------------- the fixtures


def _translation(bases: str) -> str:
    return "M" + pt.translate(bases)[1:-1]


def _record(tmp_path: pathlib.Path, bodies: list[list[str]],
            products: list[str] | None = None,
            qualifiers: list[str] | None = None,
            wrap_translation: bool = False,
            name: str = "record.gb") -> pathlib.Path:
    """A GenBank record carrying one CDS per body, each on the plus strand.

    The last gene is written on the minus strand so the reverse-complement path
    is exercised by every fixture rather than by one special case.
    """
    genes = [("ATG" + "".join(body) + "TAA") for body in bodies]
    products = products or [f"50S ribosomal protein L{i + 1}"
                            for i in range(len(genes))]
    lines, origin, cursor = [], [], 0
    for index, (bases, product) in enumerate(zip(genes, products)):
        start, end = cursor + 1, cursor + len(bases)
        cursor = end
        minus = index == len(genes) - 1
        stored = bases.translate(pt._COMPLEMENT)[::-1] if minus else bases
        origin.append(stored)
        location = f"complement({start}..{end})" if minus else f"{start}..{end}"
        lines.append(f"     CDS             {location}")
        lines.append(f'                     /locus_tag="T_{index:04d}"')
        lines.append(f'                     /product="{product}"')
        lines.append("                     /codon_start=1")
        lines.append("                     /transl_table=11")
        for extra in (qualifiers or []):
            lines.append(f"                     {extra}")
        translation = _translation(bases)
        if wrap_translation:
            head, tail = translation[:20], translation[20:]
            lines.append(f'                     /translation="{head}')
            lines.append(f'                     {tail}"')
        else:
            lines.append(f'                     /translation="{translation}"')

    sequence = "".join(origin)
    body = [f"{i + 1:>9} " + sequence[i:i + 60].lower()
            for i in range(0, len(sequence), 60)]
    text = ("\n".join([
        f"LOCUS       TEST01{len(sequence):>18} bp    DNA     circular BCT 16-FEB-2017",
        "DEFINITION  Synthetic record for the codon-table refusals.",
        "ACCESSION   TEST01",
        "FEATURES             Location/Qualifiers",
        *lines,
        "ORIGIN",
        *body,
        "//",
    ]) + "\n")
    path = tmp_path / name
    path.write_text(text)
    return path


def _bodies(genes: int = 45, extras: dict[str, int] | None = None,
            omit: tuple[str, ...] = ()) -> list[list[str]]:
    """Every sense codon once per gene, plus extras in the first gene.

    The extras break the tie: without them every codon of a family has the same
    count and the family has no single most-used codon, which is its own refusal.
    """
    base = [c for c in SENSE if c not in omit]
    bodies = [list(base) for _ in range(genes)]
    for codon, count in (extras if extras is not None
                         else {c: 10 for c in TOPS.values()}).items():
        bodies[0].extend([codon] * count)
    return bodies


def _manifest(record: pathlib.Path, gc_global=(0.0, 1.0), **overrides) -> dict:
    host = {
        "accession": "TEST01",
        "design_window": {"gc_global": list(gc_global), "gc_window": [0.3, 0.7]},
        "design_window_basis": "synthetic",
        "label": "Test host",
        "license": "none",
        "organism": "Testus fixturus",
        "source": "synthetic",
        "url": "https://example.invalid/TEST01",
    }
    host.update(overrides)
    return {
        "hosts": {"test_host": host},
        "reference_set": {
            "pseudocount": 0.5,
            "product_pattern": json.loads(pt.SOURCES.read_text())
                ["reference_set"]["product_pattern"],
            "rule": "as declared in the real manifest",
            "why_ribosomal": "as declared in the real manifest",
        },
    }


def _build(record, manifest=None, **kwargs):
    manifest = manifest or _manifest(record)
    return pt.build_tables(manifest, {"test_host": record},
                           adopt_digests=kwargs.pop("adopt_digests", True),
                           **kwargs)


def test_the_fixture_itself_derives_a_table(tmp_path):
    """The happy path, so every refusal below is the only thing that changed."""
    record = _record(tmp_path, _bodies())
    tables, audit = _build(record)
    block = tables["provenance"]["hosts"]["test_host"]
    assert block["genes"] == 45
    assert audit["hosts"]["test_host"]["census_of_dropped_candidates"][
        "product_match"] == 45
    weights = tables["hosts"]["test_host"]["weights"]
    for amino_acid, codon in TOPS.items():
        assert weights[codon] == 1.0, amino_acid
    # Every gene ends in the same stop, so the stop row reports the record
    # rather than the pseudocount.
    assert weights["TAA"] == 1.0
    assert 0.0 < weights["TGA"] < 0.05


def test_a_wrapped_translation_is_joined_without_a_separator(tmp_path):
    """A translation wraps mid-sequence; joining with a space breaks the check
    that the frame is right, and every gene would then be dropped."""
    record = _record(tmp_path, _bodies(), wrap_translation=True)
    tables, _ = _build(record)
    assert tables["provenance"]["hosts"]["test_host"]["genes"] == 45


def test_a_note_holding_a_slash_line_does_not_end_the_location(tmp_path):
    """Quote state decides continuation, not a leading slash.

    The real S. coelicolor record carries notes that wrap onto a line starting
    with '/', and a reader that treats that as a new qualifier loses the feature.
    """
    record = _record(tmp_path, _bodies(), qualifiers=[
        '/note="similar to SCO1234 fasta scores; opt: 526',
        '/unexpected: this line is inside the note',
        'and this one closes it"'])
    tables, _ = _build(record)
    assert tables["provenance"]["hosts"]["test_host"]["genes"] == 45


# ----------------------------------------------------------------- the refusals


def test_a_host_without_a_record_refuses(tmp_path):
    record = _record(tmp_path, _bodies())
    manifest = _manifest(record)
    manifest["hosts"]["second_host"] = dict(manifest["hosts"]["test_host"])
    manifest["hosts"]["second_host"]["accession"] = "TEST02"
    with pytest.raises(pt.TableError, match="Every declared host needs its record"):
        pt.build_tables(manifest, {"test_host": record}, adopt_digests=True)


def test_a_record_for_an_undeclared_host_refuses(tmp_path):
    record = _record(tmp_path, _bodies())
    with pytest.raises(pt.TableError, match="has no window"):
        pt.build_tables(_manifest(record),
                        {"test_host": record, "ghost": record}, adopt_digests=True)


def test_two_hosts_on_one_accession_refuse(tmp_path):
    record = _record(tmp_path, _bodies())
    manifest = _manifest(record)
    manifest["hosts"]["second_host"] = dict(manifest["hosts"]["test_host"])
    with pytest.raises(pt.TableError, match="pinned to the same accession"):
        pt.build_tables(manifest, {"test_host": record, "second_host": record},
                        adopt_digests=True)


@pytest.mark.parametrize("pseudocount", [0.0, 1.0, 2.0, -0.5])
def test_a_pseudocount_outside_the_unit_interval_refuses(tmp_path, pseudocount):
    record = _record(tmp_path, _bodies())
    manifest = _manifest(record)
    manifest["reference_set"]["pseudocount"] = pseudocount
    with pytest.raises(pt.TableError, match="must be between 0 and 1"):
        pt.build_tables(manifest, {"test_host": record}, adopt_digests=True)


def test_too_few_reference_genes_refuses(tmp_path):
    record = _record(tmp_path, _bodies(genes=pt.MIN_GENES - 1))
    with pytest.raises(pt.TableError, match="fewer than 40"):
        _build(record)


def test_a_record_naming_nothing_ribosomal_refuses(tmp_path):
    """The failure mode the gene floor is really for: a record whose annotation
    style the product pattern does not match selects nothing at all."""
    record = _record(tmp_path, _bodies(),
                     products=[f"hypothetical protein {i}" for i in range(45)])
    with pytest.raises(pt.TableError, match="does not match"):
        _build(record)


def test_too_many_dropped_candidates_refuse(tmp_path):
    """Half the product-matching features carry a frame the reader cannot read."""
    record = _record(tmp_path, _bodies(genes=60))
    # Eight of sixty read under a different code, which the rule will not
    # translate -- enough to leave more than MIN_GENES and still over a tenth.
    text = record.read_text().replace("/transl_table=11", "/transl_table=4", 8)
    broken = tmp_path / "dropped.gb"
    broken.write_text(text)
    with pytest.raises(pt.TableError, match="failed the validity checks"):
        _build(broken)


def test_too_many_unobserved_sense_codons_refuse(tmp_path):
    omit = tuple(SENSE[:pt.MAX_UNOBSERVED_SENSE + 1])
    extras = {c: 10 for c in TOPS.values() if c not in omit}
    record = _record(tmp_path, _bodies(extras=extras, omit=omit))
    with pytest.raises(pt.TableError, match="sense codons are unobserved"):
        _build(record)


def test_a_family_decided_by_a_handful_of_counts_refuses(tmp_path):
    """Cysteine has two codons; here they are seen five times between them."""
    cys = AA_CODONS["C"]
    bodies = _bodies(omit=cys)
    for body in bodies[:5]:
        body.append(cys[0])
    record = _record(tmp_path, bodies)
    with pytest.raises(pt.TableError, match="family is seen only"):
        _build(record)


def test_a_tied_family_refuses(tmp_path):
    """Two codons at weight 1.0 make the table's scale ambiguous."""
    extras = {c: 10 for c in TOPS.values() if c != TOPS["A"]}
    record = _record(tmp_path, _bodies(extras=extras))
    with pytest.raises(pt.TableError, match="codons at weight"):
        _build(record)


def test_a_non_positive_weight_refuses():
    """Unreachable through the manifest -- the pseudocount exists to prevent it --
    so the guard is exercised directly rather than left as a comment."""
    counts = {codon: 50 for codon in pt.CODON_TO_AA}
    weights = {codon: 1.0 if codon == AA_CODONS["A"][0] else 0.5
               for codon in pt.CODON_TO_AA}
    for amino_acid, codons in AA_CODONS.items():
        for index, codon in enumerate(codons):
            weights[codon] = 1.0 if index == 0 else 0.5
    weights[AA_CODONS["A"][1]] = 0.0
    genes = [{"locus_tag": f"T_{i}", "location": f"{i}..{i + 9}",
              "length_nt": 9} for i in range(pt.MIN_GENES)]
    census = {"product_match": pt.MIN_GENES}
    with pytest.raises(pt.TableError, match="not positive"):
        pt._assert_usable("test_host", genes, census, counts, weights, [])


def test_a_design_window_that_excludes_the_record_refuses(tmp_path):
    record = _record(tmp_path, _bodies())
    manifest = _manifest(record, gc_global=(0.9, 0.95))
    with pytest.raises(pt.TableError, match="outside the declared gc_global"):
        pt.build_tables(manifest, {"test_host": record}, adopt_digests=True)


def test_a_record_without_bases_refuses(tmp_path):
    record = _record(tmp_path, _bodies())
    text = record.read_text()
    headless = tmp_path / "headless.gb"
    headless.write_text(text[:text.index("ORIGIN")] + "//\n")
    with pytest.raises(pt.TableError, match="no ORIGIN bases"):
        _build(headless)


def test_an_unpinned_manifest_refuses_to_be_checked(tmp_path):
    record = _record(tmp_path, _bodies())
    with pytest.raises(pt.TableError, match="--adopt-digests"):
        pt.build_tables(_manifest(record), {"test_host": record})


def test_a_changed_reference_set_refuses(tmp_path):
    """The pin that is enforced: the derived SET, not the file's bytes."""
    record = _record(tmp_path, _bodies())
    manifest = _manifest(record)
    pt.build_tables(manifest, {"test_host": record}, adopt_digests=True)
    assert manifest["hosts"]["test_host"]["genes"] == 45

    # One more gene: same rule, same record style, different reference set.
    bigger = _record(tmp_path, _bodies(genes=46), name="bigger.gb")
    with pytest.raises(pt.TableError, match="does not match the manifest"):
        pt.build_tables(manifest, {"test_host": bigger})


def test_a_reformatted_record_is_not_a_refusal(tmp_path):
    """The pin is deliberately NOT the file digest: the same accession from
    another INSDC mirror is a different byte stream and the same biology."""
    record = _record(tmp_path, _bodies())
    manifest = _manifest(record)
    pt.build_tables(manifest, {"test_host": record}, adopt_digests=True)
    reformatted = tmp_path / "reformatted.gb"
    reformatted.write_text(record.read_text().replace(
        "DEFINITION  Synthetic record for the codon-table refusals.",
        "DEFINITION  Synthetic record, re-exported by another mirror."))
    tables, _ = pt.build_tables(manifest, {"test_host": reformatted})
    assert tables["provenance"]["hosts"]["test_host"]["genes"] == 45


# --------------------------------------------------------------- the CLI itself


@pytest.mark.parametrize("argv,message", [
    (["--record", "no-equals-sign"], "is not HOST=PATH"),
    (["--record", "test_host=/nonexistent/record.gb"], "no record at"),
])
def test_the_cli_refuses_bad_record_arguments(argv, message):
    with pytest.raises(pt.TableError, match=message):
        pt.main(argv)


def test_the_cli_refuses_a_repeated_host(tmp_path):
    record = _record(tmp_path, _bodies())
    with pytest.raises(pt.TableError, match="names test_host twice"):
        pt.main(["--record", f"test_host={record}",
                 "--record", f"test_host={record}"])
