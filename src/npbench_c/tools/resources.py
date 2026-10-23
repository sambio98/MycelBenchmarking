"""A pinned file a campaign reads from the image and does not redistribute.

Some reference data cannot be shipped inside a campaign: it is too large, or it
carries a licence the benchmark is not in a position to pass on. Such a resource
is declared with the path, the provider, the version and a fingerprint, so a
campaign records exactly what it was built against even though the bytes live
outside it.

This lives in `tools` rather than in a template because three campaigns now read
the same Pfam library from the image, and the readiness gate's
`external_resources_pinned_and_present` check reads the declaration out of the
campaign's rules file. One implementation of the declaration means one place
where the path, the provider and the fingerprint are interpreted.
"""

from __future__ import annotations

import hashlib
import os
import pathlib
from dataclasses import dataclass
from typing import Mapping

#: Where the pinned image puts its binaries. The suite and the templates share
#: the override so a different prefix needs no code change.
DEFAULT_PREFIX = "/opt/npbench-tools/bin"


class ResourceError(RuntimeError):
    """A declared external resource is missing, or is not the declared one."""


def tool_prefix() -> pathlib.Path:
    return pathlib.Path(os.environ.get("NPBENCH_TOOL_PREFIX", DEFAULT_PREFIX))


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Resource:
    name: str
    relative_path: str
    provider: str
    version: str
    sha256: str
    note: str

    def resolve(self) -> pathlib.Path:
        """The file as found under the image, or a refusal naming the pin.

        `relative_path` may be a glob, because a resource shipped inside a tool's
        own database layer is versioned by directory name and a campaign that
        declares the version separately should not have to restate it in the
        path. A glob that matches several files takes the LAST in sorted order,
        which is the highest version, and the fingerprint is what decides whether
        that was the declared one.
        """
        root = tool_prefix().parent
        if any(c in self.relative_path for c in "*?["):
            matches = [p for p in sorted(root.glob(self.relative_path))
                       if p.is_file()]
        else:
            candidate = root / self.relative_path
            matches = [candidate] if candidate.is_file() else []
        if not matches:
            raise ResourceError(
                f"{self.name} is not present under {root}. This campaign reads "
                f"it from the pinned image ({self.provider}) rather than shipping "
                f"it: {self.note}")
        return matches[-1]

    def verify(self) -> dict:
        """Fingerprint the resource as found, for the record."""
        path = self.resolve()
        found = sha256(path)
        return {"name": self.name, "path": str(path), "provider": self.provider,
                "version": self.version, "declared_sha256": self.sha256,
                "found_sha256": found, "matches": found == self.sha256}


def load_declared(spec: Mapping[str, Mapping]) -> dict[str, Resource]:
    """Build the resource table from a rules file's `external_resources` block."""
    return {name: Resource(name=name, relative_path=block["relative_path"],
                           provider=block["provider"], version=block["version"],
                           sha256=block["sha256"], note=block["note"])
            for name, block in sorted(spec.items())}
