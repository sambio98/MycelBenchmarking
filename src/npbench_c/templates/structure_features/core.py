"""Shared engine for the structure-features ladders.

Two ladders sit on this module because they read the same two sources for the
same proteins: UniProt's per-residue functional annotations, which carry ECO
evidence codes, and the deposited mmCIF coordinates those annotations cite.

  residues/  census the annotations, classify their evidence, then check which
             of them the cited structure's geometry actually supports
  pocket/    compute the contact shell and pocket descriptors from coordinates,
             then reconcile them against the annotation set

Building them separately would duplicate the mmCIF parser, the numbering mapper
and the ECO classifier, and the two would drift. Everything here is pure stdlib
and deterministic: no superposition, no symmetry expansion, no minimisation, no
chemistry toolkit. Reductions run in file order and every convention that could
change a number is pinned in constants/geometry_rules.json.
"""

from __future__ import annotations

import json
import math
import pathlib
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

TEMPLATE_ID = "structure-features"
CONSTANTS = pathlib.Path(__file__).resolve().parent / "constants"

#: Distances and descriptors are reported at the precision the deposited
#: coordinates support; see geometry_rules.json.
DISTANCE_DECIMALS = 2
DESCRIPTOR_DECIMALS = 3

NULL_TOKENS = (".", "?", "")


class StructureError(RuntimeError):
    """The structure or the annotation set does not support the task.

    Raised on shipped inputs only, so it is always a campaign bug: an element
    outside the pinned radius table, a residue outside the pinned alphabet, a
    per-chain numbering disagreement. Refusing beats substituting a default,
    because a default would make a wrong number look like a computed one.
    """


# --------------------------------------------------------------------- mmCIF


def _tokens(line: str) -> list[str]:
    """Tokenise one mmCIF line, honouring single and double quoting.

    A quote opens a string only at the start of a token and closes it only when
    followed by whitespace or end of line, which is what lets a residue name
    contain an apostrophe. A '#' starts a comment only at the start of a token.
    """
    out: list[str] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c == "#":
            break
        if c in "'\"":
            j = i + 1
            while j < n and not (line[j] == c and (j + 1 >= n or line[j + 1] in " \t")):
                j += 1
            out.append(line[i + 1:j])
            i = j + 1
            continue
        j = i
        while j < n and line[j] not in " \t":
            j += 1
        out.append(line[i:j])
        i = j
    return out


def read_cif(text: str, categories: Iterable[str]) -> dict[str, list[dict]]:
    """Read the named mmCIF categories into lists of row dicts.

    Both serialisations are handled, because the files in use carry both: a
    multi-row category appears as a `loop_`, and a single-row category appears
    as bare `_category.item value` lines. A parser that handles only the loop
    form silently returns nothing for the second -- which is how a numbering
    offset ends up assumed instead of read.
    """
    want = set(categories)
    data: dict[str, list[dict]] = {c: [] for c in want}
    lines = text.splitlines()
    n = len(lines)

    def semicolon_value(i: int) -> tuple[str, int]:
        buf = [lines[i][1:]]
        i += 1
        while i < n and not lines[i].startswith(";"):
            buf.append(lines[i])
            i += 1
        return "\n".join(buf).strip(), i + 1

    i = 0
    while i < n:
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue

        if stripped.lower() == "loop_":
            i += 1
            names: list[str] = []
            while i < n and lines[i].strip().startswith("_"):
                names.append(lines[i].strip().split()[0])
                i += 1
            if not names:
                continue
            category = names[0][1:].split(".", 1)[0]
            keys = [nm[1:].split(".", 1)[1] for nm in names]
            rows: list[dict] = []
            buf: list[str] = []
            while i < n:
                line = lines[i]
                st = line.strip()
                if not st:
                    i += 1
                    continue
                if (st.startswith("#") or st.startswith("_") or st.startswith("data_")
                        or st.startswith("save_") or st.lower() in ("loop_", "stop_")):
                    break
                if line.startswith(";"):
                    value, i = semicolon_value(i)
                    buf.append(value)
                else:
                    buf.extend(_tokens(line))
                    i += 1
                while len(buf) >= len(keys):
                    rows.append(dict(zip(keys, buf[:len(keys)])))
                    buf = buf[len(keys):]
            if category in want:
                data[category].extend(rows)
            continue

        if stripped.startswith("_"):
            name = stripped.split()[0]
            category, key = name[1:].split(".", 1)
            rest = stripped[len(name):].strip()
            if rest:
                value = (_tokens(rest) or [""])[0]
                i += 1
            else:
                i += 1
                if i < n and lines[i].startswith(";"):
                    value, i = semicolon_value(i)
                else:
                    value = (_tokens(lines[i]) if i < n else [""])
                    value = (value or [""])[0]
                    i += 1
            if category in want:
                if not data[category]:
                    data[category].append({})
                data[category][0][key] = value
            continue

        i += 1
    return data


def _clean(value: str | None) -> str | None:
    return None if value is None or value in NULL_TOKENS else value


# ------------------------------------------------------------------ rules


@dataclass(frozen=True)
class Rules:
    """The pinned conventions, loaded from the campaign's reference/ directory.

    The oracle reads exactly the files the sandbox ships, so a campaign whose
    gold depends on a table the agent cannot see is impossible by construction
    rather than by review.
    """

    functional_types: tuple[str, ...]
    secondary_types: tuple[str, ...]
    eco_classes: Mapping[str, str]
    class_vocabulary: tuple[str, ...]
    accepted_classes: frozenset[str]
    cutoff: float
    grid_spacing: float
    probe_radius: float
    vdw: Mapping[str, float]
    retained_atoms: frozenset[str]
    three_to_one: Mapping[str, str]

    @classmethod
    def load(cls, campaign: pathlib.Path) -> "Rules":
        ref = pathlib.Path(campaign) / "reference"
        eco = json.loads((ref / "eco_policy.json").read_text())
        geo = json.loads((ref / "geometry_rules.json").read_text())
        alphabet = json.loads((ref / "residue_alphabet.json").read_text())
        volume = geo["pocket_volume_rule"]
        return cls(
            functional_types=tuple(eco["functional_feature_types"]),
            secondary_types=tuple(eco["secondary_structure_feature_types"]),
            eco_classes=dict(eco["eco_classes"]),
            class_vocabulary=tuple(eco["eco_class_vocabulary"]),
            accepted_classes=frozenset(eco["accepted_classes_baseline"]),
            cutoff=float(geo["contact_rule"]["cutoff_angstrom"]),
            grid_spacing=float(volume["grid_spacing"]),
            probe_radius=float(volume["probe_radius"]),
            vdw={k.upper(): float(v) for k, v in volume["vdw_radii_angstrom"].items()},
            retained_atoms=frozenset(geo["truncation_rule"]["retained_atom_names"]),
            three_to_one=dict(alphabet["three_to_one"]),
        )

    #: Strongest first. A record's class is the strongest among its citations.
    CLASS_ORDER = ("experimental", "structural", "inferred")

    def eco_class(self, codes: Iterable[str]) -> str:
        seen = set()
        for code in codes:
            if code not in self.eco_classes:
                raise StructureError(
                    f"ECO code {code!r} is outside the pinned policy table")
            seen.add(self.eco_classes[code])
        for name in self.CLASS_ORDER:
            if name in seen:
                return name
        return "none"

    def one_letter(self, comp: str) -> str:
        if comp not in self.three_to_one:
            raise StructureError(
                f"residue {comp!r} is outside the pinned residue alphabet")
        return self.three_to_one[comp]

    def radius(self, element: str) -> float:
        key = element.upper()
        if key not in self.vdw:
            raise StructureError(
                f"element {element!r} is outside the pinned van der Waals table")
        return self.vdw[key]


# ------------------------------------------------------------- annotations


@dataclass(frozen=True)
class FeatureRecord:
    """One UniProt per-residue feature, after the declared dedup key is applied."""

    type: str
    start: int
    end: int
    ligand_name: str | None
    ligand_id: str | None
    original: str | None
    substitution: str | None
    citations: tuple[tuple[str, str, str], ...]

    @property
    def key(self) -> tuple:
        return (self.type, self.start, self.end, self.ligand_id, self.substitution)

    @property
    def codes(self) -> tuple[str, ...]:
        return tuple(c[0] for c in self.citations)

    def cites_structure(self, pdb_id: str) -> bool:
        """True when the record names this PDB entry as structural evidence.

        This is the field that makes structure-confirmation computable: UniProt
        does not merely say "position 348 binds chloride", it says which
        deposited structure that claim rests on.
        """
        upper = pdb_id.upper()
        return any(code == "ECO:0007744" and source == "PDB" and ident.upper() == upper
                   for code, source, ident in self.citations)


@dataclass(frozen=True)
class Annotations:
    accession: str
    sequence: str
    records: tuple[FeatureRecord, ...]
    secondary_counts: Mapping[str, int]
    pdb_entries: tuple[str, ...]

    @classmethod
    def load(cls, path: pathlib.Path, rules: Rules) -> "Annotations":
        entry = json.loads(pathlib.Path(path).read_text())
        sequence = entry["sequence"]["value"]
        functional = set(rules.functional_types)
        secondary = {t: 0 for t in rules.secondary_types}

        seen: dict[tuple, FeatureRecord] = {}
        for feature in entry.get("features") or []:
            ftype = feature["type"]
            if ftype in secondary:
                secondary[ftype] += 1
                continue
            if ftype not in functional:
                continue
            location = feature["location"]
            ligand = feature.get("ligand") or {}
            alternative = feature.get("alternativeSequence") or {}
            alts = alternative.get("alternativeSequences") or []
            citations = tuple(sorted(
                (c["evidenceCode"], str(c.get("source") or ""), str(c.get("id") or ""))
                for c in (feature.get("evidences") or [])))
            record = FeatureRecord(
                type=ftype,
                start=int(location["start"]["value"]),
                end=int(location["end"]["value"]),
                ligand_name=ligand.get("name"),
                ligand_id=ligand.get("id"),
                original=alternative.get("originalSequence"),
                substitution=",".join(alts) if alts else None,
                citations=citations,
            )
            # Declared dedup: byte-identical repeats collapse, but two ligands
            # or two substitutions at one position stay two records.
            existing = seen.get(record.key)
            if existing is None:
                seen[record.key] = record
            elif existing.citations != record.citations:
                raise StructureError(
                    f"two records share the dedup key {record.key} but cite "
                    "different evidence, so the key does not identify a record")

        pdb_entries = tuple(sorted(
            x["id"] for x in (entry.get("uniProtKBCrossReferences") or [])
            if x.get("database") == "PDB"))

        return cls(
            accession=entry["primaryAccession"],
            sequence=sequence,
            records=tuple(sorted(seen.values(), key=lambda r: r.key)),
            secondary_counts=dict(sorted(secondary.items())),
            pdb_entries=pdb_entries,
        )

    def of_type(self, ftype: str) -> tuple[FeatureRecord, ...]:
        return tuple(r for r in self.records if r.type == ftype)

    def ligands(self) -> dict[str, str]:
        """ChEBI id -> ligand name, over binding-site records. ChEBI is the key."""
        out: dict[str, str] = {}
        for record in self.records:
            if record.ligand_id:
                name = out.setdefault(record.ligand_id, record.ligand_name or "")
                if name != (record.ligand_name or ""):
                    raise StructureError(
                        f"ChEBI {record.ligand_id} carries two names: "
                        f"{name!r} and {record.ligand_name!r}")
        return dict(sorted(out.items()))

    def positions_for_ligand(self, chebi: str, pdb_id: str | None = None) -> tuple[int, ...]:
        """Binding-site positions for a ligand, optionally only those citing a PDB entry."""
        out = {r.start for r in self.records
               if r.type == "Binding site" and r.ligand_id == chebi
               and (pdb_id is None or r.cites_structure(pdb_id))}
        return tuple(sorted(out))


# -------------------------------------------------------------- structures


@dataclass(frozen=True)
class Atom:
    group: str
    element: str
    comp: str
    chain: str
    residue: str          # auth_seq_id plus insertion code when present
    label_seq: int | None
    name: str
    x: float
    y: float
    z: float

    @property
    def is_polymer(self) -> bool:
        return self.group == "ATOM"


@dataclass(frozen=True)
class Structure:
    pdb_id: str
    atoms: tuple[Atom, ...]
    water_atoms: int
    excluded_hydrogens: int
    excluded_altloc: int
    ref_seq: tuple[dict, ...]

    @classmethod
    def load(cls, path: pathlib.Path, pdb_id: str) -> "Structure":
        text = pathlib.Path(path).read_text()
        blocks = read_cif(text, ("atom_site", "struct_ref_seq"))
        rows = blocks["atom_site"]
        if not rows:
            raise StructureError(f"{path}: no _atom_site records")

        atoms: list[Atom] = []
        water = hydrogens = altloc = 0
        for row in rows:
            if _clean(row.get("pdbx_PDB_model_num")) not in (None, "1"):
                continue
            if _clean(row.get("label_alt_id")) is not None:
                altloc += 1
                continue
            element = (row.get("type_symbol") or "").upper()
            if element in ("H", "D"):
                hydrogens += 1
                continue
            comp = row["label_comp_id"]
            if comp == "HOH":
                water += 1
                continue
            ins = _clean(row.get("pdbx_PDB_ins_code"))
            label_seq = _clean(row.get("label_seq_id"))
            atoms.append(Atom(
                group=row["group_PDB"],
                element=element,
                comp=comp,
                chain=row["auth_asym_id"],
                residue=row["auth_seq_id"] + (ins or ""),
                label_seq=int(label_seq) if label_seq and label_seq.lstrip("-").isdigit() else None,
                name=row["label_atom_id"],
                x=float(row["Cartn_x"]), y=float(row["Cartn_y"]), z=float(row["Cartn_z"]),
            ))
        return cls(pdb_id=pdb_id.upper(), atoms=tuple(atoms), water_atoms=water,
                   excluded_hydrogens=hydrogens, excluded_altloc=altloc,
                   ref_seq=tuple(blocks["struct_ref_seq"]))

    # ------------------------------------------------------------ inventory

    def polymer(self) -> tuple[Atom, ...]:
        return tuple(a for a in self.atoms if a.is_polymer)

    def chains(self) -> tuple[str, ...]:
        return tuple(sorted({a.chain for a in self.polymer()}))

    def ligand_copies(self) -> dict[tuple[str, str, str], tuple[Atom, ...]]:
        """Non-water hetero groups, keyed by (comp id, chain, residue)."""
        out: dict[tuple[str, str, str], list[Atom]] = {}
        for atom in self.atoms:
            if atom.is_polymer:
                continue
            out.setdefault((atom.comp, atom.chain, atom.residue), []).append(atom)
        return {k: tuple(v) for k, v in sorted(out.items())}

    def resolved_residues(self, chain: str) -> dict[int, str]:
        """Author residue number -> comp id, for one polymer chain.

        Residues carrying an insertion code are skipped: their author number is
        not an integer, so they have no unambiguous position on this axis.
        """
        out: dict[int, str] = {}
        for atom in self.polymer():
            if atom.chain != chain or not atom.residue.lstrip("-").isdigit():
                continue
            out[int(atom.residue)] = atom.comp
        return dict(sorted(out.items()))

    # ------------------------------------------------------------ numbering

    def offsets(self, accession: str) -> dict[str, int]:
        """label_seq_id and auth_seq_id offsets to UniProt numbering.

        Every alignment record naming the accession must agree, because a
        per-chain offset would make "the" UniProt position of a contact
        ambiguous and the reconciliation meaningless.
        """
        found = []
        for row in self.ref_seq:
            if (row.get("pdbx_db_accession") or "").upper() != accession.upper():
                continue
            db = int(row["db_align_beg"])
            found.append((db - int(row["seq_align_beg"]),
                          db - int(row["pdbx_auth_seq_align_beg"])))
        if not found:
            raise StructureError(
                f"{self.pdb_id}: no _struct_ref_seq record cites {accession}")
        if len(set(found)) != 1:
            raise StructureError(
                f"{self.pdb_id}: alignment records for {accession} disagree on the "
                f"numbering offsets: {sorted(set(found))}")
        label, auth = found[0]
        return {"label_to_uniprot": label, "auth_to_uniprot": auth}

    # ------------------------------------------------------------- geometry

    def contacts(self, ligand: Sequence[Atom], cutoff: float,
                 polymer: Sequence[Atom] | None = None) -> dict[tuple[str, str, str], float]:
        """(chain, residue, comp) -> minimum heavy-atom distance within the cutoff."""
        pool = self.polymer() if polymer is None else polymer
        limit = cutoff * cutoff
        out: dict[tuple[str, str, str], float] = {}
        for atom in pool:
            best = None
            for other in ligand:
                d2 = ((atom.x - other.x) ** 2 + (atom.y - other.y) ** 2
                      + (atom.z - other.z) ** 2)
                if d2 <= limit and (best is None or d2 < best):
                    best = d2
            if best is not None:
                key = (atom.chain, atom.residue, atom.comp)
                distance = math.sqrt(best)
                if key not in out or distance < out[key]:
                    out[key] = distance
        return dict(sorted(out.items()))

    def truncated_polymer(self, chain: str, residue: str,
                          rules: Rules) -> tuple[tuple[Atom, ...], int]:
        """Polymer atoms with one residue's side chain deleted beyond CB.

        Exact arithmetic over the deposited coordinates: nothing is rebuilt and
        nothing moves, so the result answers "which contacts did this side chain
        make" and claims nothing about a real mutant.
        """
        target = [a for a in self.polymer() if a.chain == chain and a.residue == residue]
        if not target:
            raise StructureError(
                f"{self.pdb_id}: no polymer residue {chain}:{residue} to truncate")
        comp = target[0].comp
        if comp in ("GLY", "ALA"):
            raise StructureError(
                f"{self.pdb_id}: {chain}:{residue} is {comp}, which has no side "
                "chain to truncate; the campaign names the wrong position")
        kept, removed = [], 0
        for atom in self.polymer():
            if atom.chain == chain and atom.residue == residue:
                if atom.name in rules.retained_atoms:
                    kept.append(atom)
                else:
                    removed += 1
            else:
                kept.append(atom)
        return tuple(kept), removed


# --------------------------------------------------------- free functions


def centroid(atoms: Sequence[Atom]) -> tuple[float, float, float]:
    n = len(atoms)
    if not n:
        raise StructureError("centroid of an empty atom set")
    sx = sy = sz = 0.0
    for atom in atoms:          # file order, so the reduction is reproducible
        sx += atom.x
        sy += atom.y
        sz += atom.z
    return (round(sx / n, DESCRIPTOR_DECIMALS),
            round(sy / n, DESCRIPTOR_DECIMALS),
            round(sz / n, DESCRIPTOR_DECIMALS))


def radius_of_gyration(atoms: Sequence[Atom]) -> float:
    n = len(atoms)
    if not n:
        raise StructureError("radius of gyration of an empty atom set")
    cx = sum(a.x for a in atoms) / n
    cy = sum(a.y for a in atoms) / n
    cz = sum(a.z for a in atoms) / n
    total = 0.0
    for atom in atoms:
        total += (atom.x - cx) ** 2 + (atom.y - cy) ** 2 + (atom.z - cz) ** 2
    return round(math.sqrt(total / n), DESCRIPTOR_DECIMALS)


def pocket_volume(ligand: Sequence[Atom], polymer: Sequence[Atom],
                  rules: Rules) -> dict:
    """Grid-counted cavity volume under the pinned rule.

    The lattice is anchored on integers so the grid does not move when the
    bounding box does, which is what makes the count reproducible rather than
    merely repeatable.
    """
    if not ligand:
        raise StructureError("pocket volume of an empty ligand")
    cutoff, probe, step = rules.cutoff, rules.probe_radius, rules.grid_spacing
    lo = [math.floor(min(getattr(a, axis) for a in ligand) - cutoff)
          for axis in ("x", "y", "z")]
    hi = [max(getattr(a, axis) for a in ligand) + cutoff for axis in ("x", "y", "z")]

    # Only polymer atoms that could touch the box matter; the rest cannot change
    # any point's verdict, so excluding them is exact, not an approximation.
    reach = cutoff + probe + max(rules.vdw.values())
    near = [a for a in polymer
            if lo[0] - reach <= a.x <= hi[0] + reach
            and lo[1] - reach <= a.y <= hi[1] + reach
            and lo[2] - reach <= a.z <= hi[2] + reach]
    radii = [rules.radius(a.element) for a in near]
    clash2 = [(r + probe) ** 2 for r in radii]
    lining2 = (cutoff + probe) ** 2
    ligand2 = cutoff * cutoff

    counts = [0, 0, 0]
    admitted = 0
    axes = []
    for k in range(3):
        steps = int(math.floor((hi[k] - lo[k]) / step)) + 1
        axes.append([lo[k] + i * step for i in range(steps)])
    for px in axes[0]:
        for py in axes[1]:
            for pz in axes[2]:
                near_ligand = False
                for atom in ligand:
                    if ((px - atom.x) ** 2 + (py - atom.y) ** 2
                            + (pz - atom.z) ** 2) <= ligand2:
                        near_ligand = True
                        break
                if not near_ligand:
                    continue
                counts[0] += 1
                clash = False
                lining = False
                for idx, atom in enumerate(near):
                    d2 = ((px - atom.x) ** 2 + (py - atom.y) ** 2 + (pz - atom.z) ** 2)
                    if d2 < clash2[idx]:
                        clash = True
                        break
                    if d2 <= lining2:
                        lining = True
                if clash:
                    counts[1] += 1
                    continue
                if not lining:
                    counts[2] += 1
                    continue
                admitted += 1
    return {
        "grid_points_near_ligand": counts[0],
        "grid_points_rejected_clash": counts[1],
        "grid_points_rejected_bulk": counts[2],
        "pocket_grid_points": admitted,
        "pocket_volume_a3": round(admitted * step ** 3, DESCRIPTOR_DECIMALS),
    }


def contact_label(chain: str, residue: str, comp: str, rules: Rules) -> str:
    """The graded identifier for a contact residue: 'A:357:E'."""
    return f"{chain}:{residue}:{rules.one_letter(comp)}"


def uniprot_positions(contacts: Iterable[tuple[str, str, str]], chain: str,
                      auth_offset: int) -> tuple[int, ...]:
    """Contacts in one chain, expressed as UniProt positions."""
    out = set()
    for ch, residue, _comp in contacts:
        if ch != chain or not residue.lstrip("-").isdigit():
            continue
        out.add(int(residue) + auth_offset)
    return tuple(sorted(out))
