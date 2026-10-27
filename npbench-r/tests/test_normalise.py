import pytest

from npbench.grading.normalise import (
    NormalisationError,
    build_alias_map,
    normalise,
    normalise_all,
)

ALL = ("strip", "collapse_whitespace", "casefold", "strip_accession_version")


@pytest.mark.parametrize(
    "value,rules,expected",
    [
        ("  BGC0001  ", ("strip",), "BGC0001"),
        ("BGC   0001", ("collapse_whitespace",), "BGC 0001"),
        ("NRPS-Like", ("casefold",), "nrps-like"),
        ("NZ_CP011122.1", ("strip_accession_version",), "NZ_CP011122"),
        ("GCF_000203835.2", ("strip_accession_version",), "GCF_000203835"),
        ("  NZ_CP011122.1 ", ALL, "nz_cp011122"),
        # A decimal must survive: the rule is narrow on purpose.
        ("1.5", ALL, "1.5"),
        # So must an accession with no version suffix.
        ("NZ_CP011122", ALL, "nz_cp011122"),
        # And an identifier whose tail is not a bare integer.
        ("strain.A1", ALL, "strain.a1"),
        (None, ALL, ""),
    ],
)
def test_rules(value, rules, expected):
    assert normalise(value, rules) == expected


def test_rules_apply_in_declared_order():
    # casefold before strip_accession_version must not change the outcome here,
    # but the order is the author's to declare, so it is applied as written.
    assert normalise("NZ_CP011122.1", ("casefold", "strip_accession_version")) == "nz_cp011122"


def test_unknown_rule_is_an_error():
    with pytest.raises(NormalisationError, match="unknown normalisation rule"):
        normalise("x", ("lowercase",))


def test_alias_map_keys_are_normalised_too():
    # Written in the source paper's casing, used against casefolded values.
    aliases = build_alias_map({"NRPS-Like": "NRPS"}, ("casefold",))
    assert normalise_all(["nrps-like", "pks"], ("casefold",), aliases) == ["nrps", "pks"]


def test_alias_is_applied_after_rules_not_before():
    aliases = build_alias_map({"  T1PKS  ": "pks"}, ALL)
    assert normalise_all(["t1pks"], ALL, aliases) == ["pks"]
