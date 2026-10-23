"""Read a compressed tar archive the way gzip wants to be read: once, in order.

Why this module exists, measured. Two campaigns load their corpus from a
`.tar.gz` of a few thousand small documents, and both did it in a way that reads
naturally and costs a hundredfold:

- the RiPP corpus opened the archive again for every record it wanted, through a
  per-record reader closure -- 54 opens of an 11.5 MB archive;
- the architecture corpus extracted members in NAME order rather than archive
  order, which for a stream is the same defect wearing a sort.

A gzip stream cannot seek backwards. Every out-of-order `extractfile` re-inflates
from the beginning of the file, so a scan that should cost one decompression
costs one per member. Profiled on the RiPP corpus: 5.3 seconds of a 7.6 second
load, 133,873 calls to `zlib.decompress`, against 0.2 seconds for a single pass.
Multiplied by the 72 runs of a stub sweep, that was most of a 15-minute sweep.

So the rule is: one open, one pass, members in archive order, and sort the
RESULT rather than the reading. The whole uncompressed archive is held in memory
-- 10 to 12 MB for these corpora, which is the cost that was being traded for a
hundred re-inflations.
"""

from __future__ import annotations

import pathlib
import tarfile


class ArchiveError(RuntimeError):
    """The archive cannot be read as the campaign declares it."""


def read_members(path: pathlib.Path, suffix: str) -> dict[str, bytes]:
    """Every regular member whose name ends in `suffix`, read in one pass.

    Keyed by member name, in the order the archive stores them -- a dict keeps
    insertion order, so a caller that wants deterministic processing sorts the
    keys and pays nothing for it. A duplicate member name is a malformed archive
    and raises rather than silently keeping the last copy.
    """
    path = pathlib.Path(path)
    if not path.is_file():
        raise ArchiveError(f"the archive {path} is not present")
    out: dict[str, bytes] = {}
    with tarfile.open(path, "r:*") as archive:
        # Iterating the archive object streams it; getmembers() would build the
        # index first, which for a .tar.gz means an extra full inflation.
        for member in archive:
            if not member.isfile() or not member.name.endswith(suffix):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            if member.name in out:
                raise ArchiveError(
                    f"{member.name} appears twice in {path.name}")
            out[member.name] = handle.read()
    if not out:
        raise ArchiveError(f"{path.name} holds no {suffix} members")
    return out
