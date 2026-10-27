"""The network adjudicator.

PROJECT.md §3: routine traffic runs through an allowlist of package registries
and OS mirrors. Everything else goes through a request tool, where this decides
whether the resource would leak the answer.

It fails closed. Anything neither allowed nor recognised is denied, because the
cost of wrongly denying a documentation page is an agent that reads a man page
instead, while the cost of wrongly allowing a supplementary XLSX is a task that
measures nothing. Several of the candidate studies have their ground truth
sitting in a public supplementary file, so this is not a theoretical concern.

Rules are evaluated deny-first. A host on both lists is denied.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from urllib.parse import urlsplit

# Package registries, OS mirrors, tool source and documentation. Fetching a
# wheel or reading a manual cannot reveal what antiSMASH found in our genomes.
ALLOW_HOSTS = frozenset(
    {
        "pypi.org", "files.pythonhosted.org",
        "conda.anaconda.org", "repo.anaconda.com", "anaconda.org",
        "registry.npmjs.org", "cran.r-project.org",
        "deb.debian.org", "security.debian.org",
        "archive.ubuntu.com", "security.ubuntu.com",
        "github.com", "raw.githubusercontent.com", "objects.githubusercontent.com",
        "codeload.github.com",
        "docs.python.org", "readthedocs.io", "readthedocs.org",
        # Tool distributions and their reference databases.
        "dl.secondarymetabolites.org", "antismash.secondarymetabolites.org",
        "ftp.ncbi.nlm.nih.gov", "api.ncbi.nlm.nih.gov", "eutils.ncbi.nlm.nih.gov",
        "ftp.ebi.ac.uk", "www.ebi.ac.uk",
    }
)

# Publisher, preprint and aggregator hosts: the paper itself and its
# supplementary material.
DENY_HOSTS = frozenset(
    {
        "doi.org", "dx.doi.org",
        "pubs.acs.org", "journals.asm.org", "onlinelibrary.wiley.com",
        "academic.oup.com", "journals.plos.org", "www.nature.com", "nature.com",
        "www.sciencedirect.com", "sciencedirect.com", "link.springer.com",
        "www.pnas.org", "pnas.org", "www.cell.com", "cell.com",
        "www.biorxiv.org", "biorxiv.org", "www.medrxiv.org", "arxiv.org",
        "pmc.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov", "europepmc.org",
        "pubmed.ncbi.nlm.nih.gov",
        "scholar.google.com", "www.semanticscholar.org", "semanticscholar.org",
        "www.researchgate.net", "researchgate.net", "sci-hub.se", "sci-hub.st",
        # Deposits that hold processed matrices and supplementary tables. A task
        # that genuinely needs one declares it in its own allow_hosts.
        "osf.io", "zenodo.org", "figshare.com", "datadryad.org",
        "github.io",
    }
)

# Paths that are supplementary material whatever the host.
DENY_PATH_PATTERNS = (
    r"/suppl(ement(ary|al)?)?/",
    r"/doi/suppl/",
    r"_si_\d*\.",
    r"/media/\d",
    r"\bsupp(orting)?[-_ ]?information\b",
    r"\btable[-_ ]?s\d+\b",
    r"\bdataset[-_ ]?s\d+\b",
)

SEARCH_HOSTS = frozenset(
    {"www.google.com", "google.com", "duckduckgo.com", "www.bing.com", "bing.com"}
)


@dataclass(frozen=True)
class Decision:
    allowed: bool
    rule: str
    reason: str

    def __str__(self) -> str:
        return f"{'ALLOW' if self.allowed else 'DENY'} [{self.rule}] {self.reason}"


@dataclass(frozen=True)
class Policy:
    """Per-task additions to the benchmark-wide lists.

    `deny_terms` is the important field: the DOI, the author names, the dataset
    accessions and any phrase distinctive enough to find the paper. It is
    matched case-insensitively against the whole request, URL or search query
    alike.
    """

    deny_terms: tuple[str, ...] = ()
    deny_hosts: tuple[str, ...] = ()
    allow_hosts: tuple[str, ...] = ()


def _host_of(request: str) -> str | None:
    parts = urlsplit(request if "//" in request else f"//{request}")
    return (parts.hostname or "").lower() or None


def _host_matches(host: str, listed: frozenset[str] | tuple[str, ...]) -> str | None:
    """Match a host or any of its parents, so readthedocs.io covers x.readthedocs.io."""
    labels = host.split(".")
    for i in range(len(labels) - 1):
        candidate = ".".join(labels[i:])
        if candidate in listed:
            return candidate
    return None


def _looks_like_a_url(request: str) -> bool:
    return bool(re.match(r"^(https?://|ftp://|www\.)", request.strip(), re.I)) or bool(
        re.match(r"^[a-z0-9.-]+\.[a-z]{2,}(/|$|:)", request.strip(), re.I)
    )


class Adjudicator:
    def __init__(self, policy: Policy | None = None):
        self.policy = policy or Policy()
        self.allow = ALLOW_HOSTS | frozenset(h.lower() for h in self.policy.allow_hosts)
        self.deny = DENY_HOSTS | frozenset(h.lower() for h in self.policy.deny_hosts)
        self.terms = tuple(t.lower() for t in self.policy.deny_terms if t.strip())

    def decide(self, request: str) -> Decision:
        """Allow or deny one web access request: a URL, or a search query."""
        text = request.strip()
        if not text:
            return Decision(False, "empty", "no request given")
        folded = text.lower()

        # A paper-specific term anywhere in the request is decisive, whether it
        # arrives as a URL or as a search query.
        for term in self.terms:
            if term in folded:
                return Decision(
                    False, "deny_term",
                    f"request contains {term!r}, which identifies this task's source study",
                )

        if not _looks_like_a_url(text):
            return Decision(
                True, "search",
                "search query carries no task-specific term; results themselves are "
                "re-adjudicated when fetched",
            )

        host = _host_of(text)
        if host is None:
            return Decision(False, "unparseable", f"cannot read a hostname from {text!r}")

        denied = _host_matches(host, self.deny)
        if denied:
            return Decision(
                False, "deny_host",
                f"{denied} publishes papers, preprints or processed deposits",
            )

        path = urlsplit(text if "//" in text else f"//{text}").path.lower()
        for pattern in DENY_PATH_PATTERNS:
            if re.search(pattern, path):
                return Decision(
                    False, "deny_path",
                    f"path looks like supplementary material (matched {pattern!r})",
                )

        if _host_matches(host, SEARCH_HOSTS):
            return Decision(
                True, "search",
                "search engine with no task-specific term in the query",
            )

        allowed = _host_matches(host, self.allow)
        if allowed:
            return Decision(True, "allow_host", f"{allowed} serves software, data or documentation")

        return Decision(
            False, "not_allowlisted",
            f"{host} is not on the allowlist; request a specific package registry, "
            f"tool source or documentation host instead",
        )
