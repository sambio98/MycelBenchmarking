"""Campaign parameters for the structure-features ladders.

The ladder modules hold no campaign specifics: which protein, which deposited
structure is focal, which structure the counterfactual swaps to, and which
ChEBI ligand joins to which PDB chemical component all come from task.yaml.
"""

from __future__ import annotations

import json
import pathlib
import shutil
from dataclasses import dataclass

import yaml

from npbench_c.templates.structure_features.core import CONSTANTS

CONSTANT_FILES = ("eco_policy.json", "geometry_rules.json", "residue_alphabet.json")


def load_task(campaign: pathlib.Path) -> dict:
    return yaml.safe_load((pathlib.Path(campaign) / "task.yaml").read_text())


def copy_constants(campaign: pathlib.Path) -> None:
    """Publish the template's pinned tables into the campaign's reference/.

    The oracle then loads exactly the files the sandbox ships, so a campaign
    whose gold depends on a table the agent cannot reach is impossible to build
    rather than merely discouraged.
    """
    reference = pathlib.Path(campaign) / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    for name in CONSTANT_FILES:
        shutil.copyfile(CONSTANTS / name, reference / name)


@dataclass(frozen=True)
class StructureSpec:
    """One campaign's protein, its structures and the declared join keys."""

    accession: str
    annotation_input: str
    focal_pdb_id: str
    focal_input: str
    focal_chain: str
    alternate_pdb_id: str | None
    alternate_input: str | None
    ligand_components: dict[str, str]        # ChEBI id -> PDB chemical component id
    focal_ligand_chebi: str | None           # the pocket ladder's graded ligand
    truncation_scope: str | None             # pocket ladder R4; "focal_shell"
    cutoff_variants: tuple[float, ...]
    policy_variants: dict[str, dict]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "StructureSpec":
        params = load_task(campaign)["template_params"]
        structures = params["structures"]
        focal = structures["focal"]
        alternate = structures.get("alternate") or {}
        return cls(
            accession=params["accession"],
            annotation_input=params["annotation_input"],
            focal_pdb_id=focal["pdb_id"].upper(),
            focal_input=focal["input"],
            focal_chain=focal["chain"],
            alternate_pdb_id=(alternate.get("pdb_id") or "").upper() or None,
            alternate_input=alternate.get("input"),
            ligand_components={k: v for k, v in sorted(
                (params.get("ligand_components") or {}).items())},
            focal_ligand_chebi=params.get("focal_ligand_chebi"),
            truncation_scope=params.get("truncation_scope"),
            cutoff_variants=tuple(float(c) for c in (params.get("cutoff_variants") or ())),
            policy_variants={k: v for k, v in sorted(
                (params.get("policy_variants") or {}).items())},
        )

    def write_reference(self, campaign: pathlib.Path) -> None:
        """Publish the campaign's own declared join keys and focal choices.

        Without this the campaign is unsolvable: nothing else in the sandbox
        says that ChEBI:CHEBI:57912 is the ligand deposited as TRP, or which
        chain and which structure the graded numbers are computed over. These
        are conventions, not findings, so handing them over gives nothing away.
        """
        reference = pathlib.Path(campaign) / "reference"
        reference.mkdir(parents=True, exist_ok=True)
        payload = {
            "description": "Declared join keys and focal choices for this "
                           "campaign. ChEBI identifies a ligand in UniProt; a "
                           "three-character chemical component id identifies it "
                           "in deposited coordinates. Nothing in either file "
                           "states the correspondence, so it is declared here "
                           "rather than guessed from the ligand's name.",
            "accession": self.accession,
            "focal_structure": {"pdb_id": self.focal_pdb_id, "chain": self.focal_chain},
            "alternate_structure": ({"pdb_id": self.alternate_pdb_id}
                                    if self.alternate_pdb_id else None),
            "ligand_components": self.ligand_components,
            "focal_ligand_chebi": self.focal_ligand_chebi,
            "truncation_scope": self.truncation_scope,
            "cutoff_variants": list(self.cutoff_variants),
            "policy_variants": self.policy_variants,
            "position_convention": "Every position this campaign grades is a "
                                   "UniProt sequence position, 1-based. Author "
                                   "and entity numbering in the coordinate files "
                                   "are mapped onto it by the rule in "
                                   "geometry_rules.json#numbering_rule.",
        }
        (reference / "campaign_keys.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n")
