"""The task 02 reference solution.

Truth is whatever this script produces, so a bug here is a bug in the ground
truth and every score computed against it. The family assignment is the part
worth guarding: it uses an inverted index and union-find for speed, and a
missed pair would silently change the partition.
"""

import gzip
import importlib.util
import itertools
import random
from pathlib import Path

import pytest

SOLUTION = (
    Path(__file__).resolve().parents[1]
    / "tasks" / "02_salinispora_bgc" / "solution" / "make_truth.py"
)


@pytest.fixture(scope="module")
def ref():
    spec = importlib.util.spec_from_file_location("task02_solution", SOLUTION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def brute_force_families(ref, domains):
    """Single linkage the obvious way: every pair, then merge components."""
    eligible = sorted(c for c, d in domains.items() if len(d) >= ref.MIN_DOMAINS)
    component = {c: {c} for c in eligible}
    owner = {c: c for c in eligible}
    for a, b in itertools.combinations(eligible, 2):
        x, y = domains[a], domains[b]
        if len(x & y) / len(x | y) >= ref.JACCARD_THRESHOLD:
            ca, cb = owner[a], owner[b]
            if ca != cb:
                component[ca] |= component[cb]
                for member in component[cb]:
                    owner[member] = ca
                del component[cb]
    return {m: min(component[owner[m]]) for m in eligible}


class TestDomainSet:
    def test_deduplicates_the_multiset(self, ref):
        # GECCO lists a domain once per occurrence. Jaccard is over sets.
        assert ref.domain_set("PF00165;PF00165;PF00165;PF00106") == {"PF00165", "PF00106"}

    def test_empty_field(self, ref):
        assert ref.domain_set("") == frozenset()


class TestAssignFamilies:
    def test_matches_brute_force_on_random_inputs(self, ref):
        # The check that matters: the inverted index must not miss a pair.
        rng = random.Random(20261027)
        pool = [f"PF{i:05d}" for i in range(60)]
        for trial in range(25):
            domains = {
                f"c{j:03d}": frozenset(rng.sample(pool, rng.randint(1, 12)))
                for j in range(60)
            }
            assert ref.assign_families(domains) == brute_force_families(ref, domains), (
                f"trial {trial} disagreed with brute force"
            )

    def test_thin_clusters_are_unassigned(self, ref):
        domains = {"a": frozenset("AB"), "b": frozenset("ABC"), "c": frozenset("ABC")}
        families = ref.assign_families(domains)
        assert "a" not in families           # 2 distinct domains, below MIN_DOMAINS
        assert families == {"b": "b", "c": "b"}

    def test_family_is_named_for_its_smallest_member(self, ref):
        domains = {"zz": frozenset("ABC"), "aa": frozenset("ABC"), "mm": frozenset("ABC")}
        assert set(ref.assign_families(domains).values()) == {"aa"}

    def test_single_linkage_chains_through_an_intermediate(self, ref):
        # J(a,b) = 3/5 and J(b,c) = 3/5, both above threshold.
        # J(a,c) = 2/6 = 0.33, below it. Single linkage still puts all three
        # in one family; complete linkage would not, which is the distinction
        # an agent has to get right.
        a, b, c = frozenset("ABCD"), frozenset("ABCE"), frozenset("ABEF")
        assert len(a & b) / len(a | b) >= ref.JACCARD_THRESHOLD
        assert len(b & c) / len(b | c) >= ref.JACCARD_THRESHOLD
        assert len(a & c) / len(a | c) < ref.JACCARD_THRESHOLD
        families = ref.assign_families({"a": a, "b": b, "c": c})
        assert families == {"a": "a", "b": "a", "c": "a"}

    def test_disjoint_clusters_are_singletons(self, ref):
        families = ref.assign_families({"a": frozenset("ABC"), "b": frozenset("XYZ")})
        assert families == {"a": "a", "b": "b"}

    def test_is_order_independent(self, ref):
        rng = random.Random(7)
        pool = [f"PF{i:05d}" for i in range(30)]
        domains = {f"c{j}": frozenset(rng.sample(pool, 6)) for j in range(40)}
        shuffled = dict(sorted(domains.items(), key=lambda kv: kv[0][::-1]))
        assert ref.assign_families(domains) == ref.assign_families(shuffled)


class TestHelpers:
    @pytest.mark.parametrize(
        "organism,expected",
        [
            ("Salinispora arenicola CNT849", "Salinispora arenicola"),
            ("Salinispora arenicola", "Salinispora arenicola"),
            ("Salinispora sp.", "Salinispora sp."),
            ("Salinispora tropica CNB-440", "Salinispora tropica"),
        ],
    )
    def test_species_of_collapses_inconsistent_organism_strings(self, ref, organism, expected):
        assert ref.species_of(organism) == expected

    @pytest.mark.parametrize(
        "values,expected",
        [
            ([5], 5),
            ([1, 3], 2),
            ([1, 2, 3], 2),
            ([4, 1, 3, 2], 2),      # even count: mean of the two middles, floored
            ([10, 20, 30, 41], 25),
        ],
    )
    def test_median_is_an_integer(self, ref, values, expected):
        assert ref._median(values) == expected

    def test_contig_lengths_reads_a_gzipped_multi_fasta(self, ref, tmp_path):
        path = tmp_path / "g.fna.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(">c1 some description\nACGT\nACG\n>c2\nAC\n")
        assert ref.contig_lengths(path) == {"c1": 7, "c2": 2}

    def test_contig_lengths_on_an_empty_file(self, ref, tmp_path):
        path = tmp_path / "e.fna.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("")
        assert ref.contig_lengths(path) == {}
