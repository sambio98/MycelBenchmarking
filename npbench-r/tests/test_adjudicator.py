"""The network adjudicator. It fails closed, so the default-deny cases matter most."""

import pytest

from npbench.harness.adjudicator import Adjudicator, Policy

TASK_POLICY = Policy(
    deny_terms=("10.1111/1462-2920.13867", "Letzel", "MSV000078836"),
)


@pytest.fixture
def adj():
    return Adjudicator(TASK_POLICY)


@pytest.fixture
def bare():
    return Adjudicator()


class TestAllowed:
    @pytest.mark.parametrize(
        "url",
        [
            "https://pypi.org/simple/pandas/",
            "https://files.pythonhosted.org/packages/ab/cd/gecco-0.11.0-py3-none-any.whl",
            "https://conda.anaconda.org/bioconda/linux-64/repodata.json",
            "https://deb.debian.org/debian/dists/bookworm/InRelease",
            "https://github.com/zellerlab/GECCO",
            "https://raw.githubusercontent.com/zellerlab/GECCO/master/README.md",
            "https://gecco.readthedocs.io/en/latest/",
            "https://docs.python.org/3/library/csv.html",
            "https://dl.secondarymetabolites.org/releases/7.1.0/",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/425/",
        ],
    )
    def test_software_data_and_docs_are_allowed(self, adj, url):
        assert adj.decide(url).allowed, adj.decide(url)

    def test_subdomains_of_an_allowed_host(self, adj):
        assert adj.decide("https://anything.readthedocs.io/page").allowed

    def test_a_search_with_no_task_specific_term(self, adj):
        decision = adj.decide("gecco cluster type classification defaults")
        assert decision.allowed and decision.rule == "search"


class TestDenied:
    @pytest.mark.parametrize(
        "url",
        [
            "https://doi.org/10.1111/x",
            "https://pubs.acs.org/doi/10.1021/acs.jnatprod.6b00722",
            "https://onlinelibrary.wiley.com/doi/10.1111/1462-2920",
            "https://journals.asm.org/doi/10.1128/mbio.02700-21",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5368041/",
            "https://europepmc.org/article/MED/28628233",
            "https://www.researchgate.net/publication/317",
            "https://scholar.google.com/scholar?q=salinispora",
            "https://osf.io/entqc/",
            "https://zenodo.org/records/1234",
            "https://figshare.com/articles/dataset/1234",
        ],
    )
    def test_papers_preprints_and_deposits(self, bare, url):
        assert not bare.decide(url).allowed, bare.decide(url)

    @pytest.mark.parametrize(
        "url",
        [
            "https://example.org/doi/suppl/10.1021/x/suppl_file/np6b00722_si_001.xlsx",
            "https://example.org/cms/supplementary/table_s1.xlsx",
            "https://example.org/articles/media/4567",
            "https://example.org/data/dataset_s1.zip",
        ],
    )
    def test_supplementary_paths_on_any_host(self, bare, url):
        decision = bare.decide(url)
        assert not decision.allowed
        assert decision.rule in {"deny_path", "not_allowlisted"}

    def test_an_unknown_host_is_denied_by_default(self, bare):
        decision = bare.decide("https://some-random-mirror.example.net/data.tar.gz")
        assert not decision.allowed
        assert decision.rule == "not_allowlisted"

    def test_deny_terms_beat_an_allowed_host(self, adj):
        # The decisive rule: an allowlisted host cannot be used as a carrier for
        # a paper-specific lookup.
        decision = adj.decide("https://github.com/search?q=MSV000078836")
        assert not decision.allowed
        assert decision.rule == "deny_term"

    def test_deny_terms_apply_to_searches(self, adj):
        assert not adj.decide("Letzel Salinispora biosynthetic gene clusters").allowed

    def test_deny_terms_are_case_insensitive(self, adj):
        assert not adj.decide("LETZEL salinispora").allowed

    def test_the_doi_is_denied_bare(self, adj):
        assert not adj.decide("10.1111/1462-2920.13867").allowed

    def test_empty_request(self, adj):
        assert not adj.decide("   ").allowed


class TestPolicyComposition:
    def test_a_task_can_widen_the_denylist(self):
        adj = Adjudicator(Policy(deny_hosts=("ftp.ncbi.nlm.nih.gov",)))
        # Allowlisted benchmark-wide, denied for this task: deny wins.
        assert not adj.decide("https://ftp.ncbi.nlm.nih.gov/genomes/").allowed

    def test_a_task_can_widen_the_allowlist(self):
        adj = Adjudicator(Policy(allow_hosts=("mibig.secondarymetabolites.org",)))
        assert adj.decide("https://mibig.secondarymetabolites.org/repository/").allowed

    def test_a_task_cannot_allowlist_its_way_past_a_deny_term(self):
        adj = Adjudicator(
            Policy(deny_terms=("MSV000078836",), allow_hosts=("massive.ucsd.edu",))
        )
        assert not adj.decide("https://massive.ucsd.edu/MSV000078836/").allowed

    def test_blank_deny_terms_are_ignored(self):
        adj = Adjudicator(Policy(deny_terms=("", "   ")))
        assert adj.decide("https://pypi.org/simple/").allowed


def test_decision_renders_for_a_log(adj):
    assert str(adj.decide("https://pypi.org/simple/")).startswith("ALLOW [allow_host]")
    assert str(adj.decide("https://doi.org/10")).startswith("DENY [deny_host]")
