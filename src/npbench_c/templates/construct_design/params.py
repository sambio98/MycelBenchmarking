"""Campaign parameters for the construct-design templates.

Read from task.yaml so the ladder modules hold no campaign specifics: which
protein, which hosts, and (for the conflict ladder) which constraint pair was
tightened into infeasibility.
"""

from __future__ import annotations

import json
import pathlib
import shutil
from dataclasses import dataclass

import yaml

from npbench_c.templates.construct_design.core import CONSTANTS


def load_task(campaign: pathlib.Path) -> dict:
    return yaml.safe_load((pathlib.Path(campaign) / "task.yaml").read_text())


def copy_constants(campaign: pathlib.Path) -> None:
    """Copy the template's pinned tables into the campaign's reference/.

    The campaign's sandbox must be self-contained -- the agent needs every table
    its gold depends on -- while the template stays the single source of truth.
    """
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    for name in ("codon_tables.json", "cloning_strategies.json",
                 "composition_limits.json"):
        shutil.copyfile(CONSTANTS / name, reference / name)


def read_protein(campaign: pathlib.Path, relative: str, add_stop: bool = True) -> str:
    path = pathlib.Path(campaign) / relative
    lines = path.read_text().strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"{path.name}: not a FASTA record")
    seq = "".join(l.strip().upper() for l in lines[1:] if not l.startswith(">"))
    return seq + "*" if add_stop else seq


@dataclass(frozen=True)
class ConflictSpec:
    """The constructed conflict.

    Infeasibility is CONSTRUCTED, never curated: the campaign tightens two
    constraints past the point where any sequence can satisfy both, so
    "infeasible" is true by arithmetic and regenerable. Curating "an expert
    judged this infeasible" would put human judgement into gold.
    """

    forbidden_codons: tuple[str, ...]
    gc_global_min: float
    host: str
    strategy: str
    protein_input: str
    expect: str          # "infeasible" | "feasible", verified against the arithmetic

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "ConflictSpec":
        params = load_task(campaign)["template_params"]
        conflict = params["conflict"]
        return cls(
            forbidden_codons=tuple(sorted(conflict["forbidden_codons"])),
            gc_global_min=float(conflict["gc_global_min"]),
            host=params["host"],
            strategy=params["strategy"],
            protein_input=params["protein_input"],
            expect=conflict["expect"],
        )

    def write_reference(self, campaign: pathlib.Path) -> None:
        """Publish the conflict's own constraints to the agent.

        Without this the campaign is unsolvable: the forbidden-codon set and the
        GC requirement ARE the constraints under test, and nothing else in the
        sandbox states them.
        """
        reference = pathlib.Path(campaign) / "reference"
        reference.mkdir(parents=True, exist_ok=True)
        (reference / "conflict_constraints.json").write_text(json.dumps({
            "description": "The constraint set this campaign asks you to assess. "
                           "forbidden_codons and gc_global_min are additional to "
                           "the composition constraints in "
                           "composition_limits.json; gc_global_min REPLACES the "
                           "host's gc_global window for this campaign.",
            "host": self.host,
            "cloning_strategy": self.strategy,
            "forbidden_codons": list(self.forbidden_codons),
            "forbidden_codons_rationale": "Disallowed by the synthesis vendor for "
                                          "this order; treat the list as given.",
            "gc_global_min": self.gc_global_min,
            "constraint_vocabulary": list(CONFLICT_CONSTRAINTS),
            "achievable_gc_rule": "For each residue independently, the lowest and "
                                  "highest GC content reachable with an allowed "
                                  "codon; summed over residues and divided by the "
                                  "ORF length in nucleotides. Codon choices are "
                                  "independent, so these bounds are attained and "
                                  "therefore exact.",
            "verdict_vocabulary": list(VERDICTS),
        }, indent=2, sort_keys=True) + "\n")


# The constraints this ladder reasons about. The host's gc_global window is
# replaced by the explicit gc_global_min requirement, so it is not listed.
CONFLICT_CONSTRAINTS: tuple[str, ...] = (
    "forbidden_codons",
    "gc_global_min",
    "forbidden_sites",
    "gc_window",
    "homopolymer",
    "direct_repeat",
    "internal_rbs",
)

VERDICTS: tuple[str, ...] = (
    "feasible",
    "infeasible_gc_unreachable",
    "infeasible_no_admissible_codon",
    "infeasible_other",
    "insufficient_information",
)
