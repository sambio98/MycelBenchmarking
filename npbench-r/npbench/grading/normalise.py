"""Identifier normalisation, applied before any comparison.

PROJECT.md §7: case, whitespace, common accession formatting, and a declared
alias map per artifact. Normalisation is declared per artifact in grading.yaml;
there is no implicit default, so a task author always sees which rules a score
was computed under.
"""

from __future__ import annotations

import re

RULES = ("strip", "collapse_whitespace", "casefold", "strip_accession_version")

# A bare accession with a trailing version suffix: NZ_CP011122.1, GCF_000203835.2.
# Deliberately narrow. It must start with a letter and hold only word characters
# before the suffix, so a decimal number is never mangled into an integer.
_VERSIONED_ACCESSION = re.compile(r"[A-Za-z][A-Za-z0-9_]*\.\d+")


class NormalisationError(ValueError):
    """An unknown normalisation rule was declared."""


def normalise(value: object, rules: tuple[str, ...]) -> str:
    """Apply `rules`, in the order given, to a single value."""
    text = "" if value is None else str(value)
    for rule in rules:
        if rule == "strip":
            text = text.strip()
        elif rule == "collapse_whitespace":
            text = " ".join(text.split())
        elif rule == "casefold":
            text = text.casefold()
        elif rule == "strip_accession_version":
            if _VERSIONED_ACCESSION.fullmatch(text):
                text = text.rsplit(".", 1)[0]
        else:
            raise NormalisationError(
                f"unknown normalisation rule {rule!r}; known rules: {', '.join(RULES)}"
            )
    return text


def build_alias_map(aliases: dict[str, str], rules: tuple[str, ...]) -> dict[str, str]:
    """Normalise both sides of a declared alias map.

    The map is written by a human in whatever casing the source used, so its
    keys are normalised too. Otherwise an alias silently fails to fire on an
    artifact that declares `casefold`.
    """
    return {normalise(k, rules): normalise(v, rules) for k, v in aliases.items()}


def normalise_all(
    values: object, rules: tuple[str, ...], alias_map: dict[str, str] | None = None
) -> list[str]:
    """Normalise an iterable of values, then resolve aliases."""
    aliases = alias_map or {}
    out = []
    for value in values:
        text = normalise(value, rules)
        out.append(aliases.get(text, text))
    return out
