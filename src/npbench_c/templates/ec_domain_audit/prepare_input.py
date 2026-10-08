"""Fetch the EC-annotated sequence set from UniProt and pin it.

Usage: python -m npbench_c.templates.ec_domain_audit.prepare_input <campaign> \
           [--ec 1.1.3.15]

Run once at instantiation. The campaign ships the FASTA, so the gold does not
depend on UniProt being reachable or unchanged at grading time, and the
provenance file records the release, the exact queries and the counts.

Both UniProt sections are taken in full. Sampling the unreviewed set would put a
selection decision inside the input, and the reviewed/unreviewed split is the
campaign's control -- it has to be the whole of each side or it controls nothing.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pathlib
import sys
import urllib.parse
import urllib.request

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
TIMEOUT_S = 600


def _get(url: str, params: dict) -> tuple[bytes, dict]:
    """Fetch, returning the body and a lower-cased header map.

    HTTP header names are case-insensitive and UniProt sends X-Total-Results in
    title case, so the map is normalised here rather than at each use.
    """
    query = urllib.parse.urlencode(params)
    request = urllib.request.Request(f"{url}?{query}",
                                     headers={"Accept-Encoding": "identity"})
    with urllib.request.urlopen(request, timeout=TIMEOUT_S) as response:
        headers = {k.lower(): v for k, v in response.headers.items()}
        return response.read(), headers


def fetch(ec: str) -> tuple[str, dict, dict]:
    """Return the combined FASTA, the per-section counts and the release."""
    parts, counts, release = [], {}, None
    for reviewed, section in (("true", "reviewed"), ("false", "unreviewed")):
        query = f"ec:{ec} AND reviewed:{reviewed}"
        head, headers = _get(SEARCH, {"query": query, "size": "1",
                                      "fields": "accession", "format": "tsv"})
        counts[section] = int(headers["x-total-results"])
        release = headers.get("x-uniprot-release", release)
        body, _ = _get(STREAM, {"query": query, "format": "fasta"})
        text = body.decode()
        # Count header LINES, not '>' characters: seven UniProt descriptions in
        # this set contain 'Delta 5-->4-isomerase', and counting characters
        # overcounts by exactly that and fails the completeness check for the
        # wrong reason.
        got = sum(1 for line in text.splitlines() if line.startswith(">"))
        if got != counts[section]:
            raise RuntimeError(
                f"{section}: the search reports {counts[section]} entries and the "
                f"stream returned {got}; a partial download would silently change "
                "the campaign's denominator")
        parts.append(text if text.endswith("\n") else text + "\n")
    return "".join(parts), counts, release


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--ec", default="1.1.3.15")
    args = parser.parse_args(argv)

    campaign = pathlib.Path(args.campaign).resolve()
    inputs = campaign / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)

    fasta, counts, release = fetch(args.ec)
    name = f"ec_{args.ec.replace('.', '_')}_uniprot_{release}.fasta.gz"
    target = inputs / name
    # mtime=0 so the shipped bytes are a function of the content alone.
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(target, "wb"),
                       mtime=0) as handle:
        handle.write(fasta.encode())

    provenance = {
        "source": "UniProtKB REST",
        "release": release,
        "license": "CC BY 4.0",
        "attribution": ("UniProt Consortium, UniProtKB release "
                        f"{release} (https://www.uniprot.org/), CC BY 4.0"),
        "queries": {section: f"ec:{args.ec} AND reviewed:{reviewed}"
                    for reviewed, section in (("true", "reviewed"),
                                              ("false", "unreviewed"))},
        "counts_by_section": counts,
        "sequences_total": sum(counts.values()),
        "file": name,
        "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "note": ("Both sections in full, no sampling: the reviewed/unreviewed "
                 "split is the campaign's control."),
    }
    (inputs / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"{name}: {provenance['sequences_total']} sequences "
          f"({counts}), UniProt {release}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
