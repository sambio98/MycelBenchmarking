"""The web access request tool.

PROJECT.md §3: routine traffic goes through an allowlist; anything else is
requested here and adjudicated. The agent sees the decision and, when allowed,
the content. A denial says why, so the agent can ask for something legitimate
instead of guessing.

The fetcher is injected so the decision logic can be tested without network.
"""

from __future__ import annotations

import urllib.request
from typing import Awaitable, Callable

from inspect_ai.tool import Tool, tool

from .adjudicator import Adjudicator, Policy

MAX_BYTES = 256_000

Fetcher = Callable[[str], Awaitable[str]]


async def http_fetch(url: str) -> str:
    """Default fetcher: a capped GET. Runs only after the adjudicator allows."""
    request = urllib.request.Request(url, headers={"User-Agent": "npbench-r/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        body = response.read(MAX_BYTES + 1)
    text = body[:MAX_BYTES].decode("utf-8", errors="replace")
    if len(body) > MAX_BYTES:
        text += f"\n\n[truncated at {MAX_BYTES} bytes]"
    return text


def decide_and_describe(adjudicator: Adjudicator, request: str) -> tuple[bool, str]:
    """Adjudicate one request and render the message the agent will read."""
    decision = adjudicator.decide(request)
    if decision.allowed:
        return True, f"ALLOWED: {decision.reason}"
    return False, (
        f"DENIED: {decision.reason}.\n\n"
        f"The allowlist covers package registries, OS mirrors, tool source and "
        f"tool documentation. Published results for this dataset are out of "
        f"scope: the analysis has to be run, not looked up."
    )


@tool
def request_web_access(
    adjudicator: Adjudicator | None = None, fetcher: Fetcher | None = None
) -> Tool:
    async def execute(request: str) -> str:
        """Request a web resource that is not on the network allowlist.

        Use this for a package index, a tool's source repository or its
        documentation. Requests for the published results of this dataset are
        denied. Pass either a URL or a search query.

        Args:
            request: A URL to fetch, or a search query to run.

        Returns:
            The adjudicator's decision, and the page content if allowed.
        """
        allowed, message = decide_and_describe(adjudicator or Adjudicator(Policy()), request)
        if not allowed:
            return message
        if not request.lower().startswith(("http://", "https://")):
            return f"{message}\n\nThis is a query, not a URL. Fetch a specific URL to read it."
        try:
            body = await (fetcher or http_fetch)(request)
        except Exception as exc:
            return f"{message}\n\nFetch failed: {type(exc).__name__}: {exc}"
        return f"{message}\n\n{body}"

    return execute
