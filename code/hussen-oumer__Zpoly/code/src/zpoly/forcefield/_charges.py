"""User charges (QM / RESP) for selected atoms"""

# rq-cb6ca230

from __future__ import annotations

import logging
import math
import os
from collections.abc import Mapping, Sequence
from pathlib import Path

import yaml
from rdkit import Chem

from ._errors import ChargeFileError
from ._model import ChargeReport, Topology

log = logging.getLogger(__name__)
LARGE_CORRECTION = 0.1  # e per atom
SCALE = 10_000  # charges are written with 4 decimals


# rq-31b6ec98
def load_charge_file(path: str | os.PathLike) -> dict[str, dict[str, float]]:
    """Read a charge file"""
    p = str(path)
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except OSError as exc:
        raise ChargeFileError(f"cannot read the charge file ({exc.strerror or exc})", p) from exc
    except yaml.YAMLError as exc:
        raise ChargeFileError(f"not valid YAML ({exc})", p) from exc
    return check_charge_mapping(data, p)


def check_charge_mapping(data, path: str | None = None) -> dict[str, dict[str, float]]:
    if data is None:
        return {}  # an empty (or fully commented) file: no user charges
    if not isinstance(data, Mapping):
        raise ChargeFileError("a charge file must map residue names to {atom: charge} tables", path)
    out: dict[str, dict[str, float]] = {}
    for residue, atoms in data.items():
        if atoms is None:
            continue  # a template residue with every atom still commented out
        if not isinstance(atoms, Mapping):
            raise ChargeFileError(f"residue {residue}: expected a table of atom: charge", path, str(residue))
        table: dict = {}
        for atom, value in atoms.items():
            where = f"residue {residue}, atom {atom}"
            if isinstance(value, Mapping):
                keys = {str(k) for k in value}
                if not keys or keys - {"charge", "H"}:
                    raise ChargeFileError(f"{where}: use {{charge: q, H: qH}} (one or both keys), got {dict(value)!r}",
                                          path, str(residue), str(atom))
                table[str(atom)] = {k: _number(value[k], f"{where}, {k}", path, residue, atom) for k in keys}
            else:
                table[str(atom)] = _number(value, where, path, residue, atom)
        out[str(residue)] = table
    return out


def round_to_target(values: Sequence[float], target: float, adjustable: Sequence[int]) -> tuple[float, ...]:
    """Round to 4 decimals with an exact sum."""
    exact = [q * SCALE for q in values]
    ints = [round(x) for x in exact]
    adjustable = list(adjustable)
    for i in adjustable:
        ints[i] = math.floor(exact[i])
    missing = round(target * SCALE) - sum(ints)
    if adjustable and 0 <= missing <= len(adjustable):
        by_fraction = sorted(adjustable, key=lambda i: (-(exact[i] - ints[i]), i))
        for i in by_fraction[:missing]:
            ints[i] += 1
    elif adjustable:  # far from target (should not happen after balancing): spread evenly
        for k in range(abs(missing)):
            ints[adjustable[k % len(adjustable)]] += 1 if missing > 0 else -1
    return tuple(i / SCALE for i in ints)


def _number(value, where: str, path, residue, atom) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ChargeFileError(f"{where}: charge must be a number, got {value!r}", path, str(residue), str(atom))
    return float(value)


def expand_entries(residue: str, entries: Mapping, names: Sequence[str],
                   hydrogens_of: Mapping[str, Sequence[str]], path: str | None = None) -> dict[str, float]:
    """Turn {atom"""
    known = set(names)
    out: dict[str, float] = {}

    def put(atom: str, q: float, source: str) -> None:
        if atom in out and abs(out[atom] - q) > 1e-12:
            raise ChargeFileError(f"residue {residue}: {atom} is given two different charges "
                                  f"({out[atom]} and {q}, via {source})", path, residue, atom)
        out[atom] = q

    for atom, value in entries.items():
        if atom not in known:
            raise ChargeFileError(f"residue {residue} has no atom named {atom} (atoms: {', '.join(names)})",
                                  path, residue, atom)
        if isinstance(value, Mapping):
            if "charge" in value:
                put(atom, value["charge"], atom)
            if "H" in value:
                hs = hydrogens_of.get(atom, ())
                if not hs:
                    raise ChargeFileError(f"residue {residue}: {atom} has no hydrogens, so 'H' cannot be used",
                                          path, residue, atom)
                for h in hs:
                    put(h, value["H"], f"{atom}'s H")
        else:
            put(atom, value, atom)
    return out


# rq-3541c593
def apply_user_charges(residue: str, names: Sequence[str], backend: Sequence[float], target: int,
                       user: Mapping[str, float], path: str | None = None):
    """Keep user charges, balance the rest to the integer target"""
    unknown = sorted(set(user) - set(names))
    if unknown:
        raise ChargeFileError(f"residue {residue} has no atom named {unknown[0]} "
                              f"(atoms: {', '.join(names)})", path, residue, unknown[0])
    fixed = [i for i, nm in enumerate(names) if nm in user]
    free = [i for i, nm in enumerate(names) if nm not in user]
    values = [user[nm] if nm in user else q for nm, q in zip(names, backend, strict=True)]
    if not free:
        if abs(sum(values) - target) > 1e-4:
            raise ChargeFileError(f"residue {residue}: all atoms are given but they sum to {sum(values):.4f}, "
                                  f"not {target}", path, residue)
        correction = 0.0
    else:
        correction = (target - sum(values)) / len(free)
        for i in free:
            values[i] += correction
    # free atoms absorb the rounding
    rounded = round_to_target(values, target, free or fixed)
    large = abs(correction) > LARGE_CORRECTION
    if large:
        log.warning("User charges in %s needed a correction of %+.4f e on each other atom; were they fitted "
                    "for the same fragment and charge state?", residue, correction)
    return rounded, ChargeReport(residue, len(fixed), correction, large)


def _neighbour_label(mol: Chem.Mol, atom_idx: int, resnum: int) -> str:
    atom = mol.GetAtomWithIdx(atom_idx)
    r = atom.GetIntProp("resnum")
    prefix = "" if r == resnum else ("-" if r < resnum else "+")
    return prefix + atom.GetProp("name")


# rq-194370d8
def charge_template(structure, topology: Topology | None = None) -> str:
    """A ready-to-edit charge file"""
    mol = structure.mol
    by_res: dict[int, list[int]] = {}
    for a in mol.GetAtoms():
        by_res.setdefault(a.GetIntProp("resnum"), []).append(a.GetIdx())
    unit_of = {}
    head = 1 if structure.chain.head_group is not None else 0
    for k, z in enumerate(structure.units):
        unit_of[k + 1 + head] = z.name
    seen: set[str] = set()
    lines = ["# Charge file for zpoly. Uncomment the atoms you have QM/RESP (e.g. DFT) charges for;",
             "# the other atoms of each residue are balanced so it keeps its integer charge.",
             "# {charge: q, H: qH} sets an atom and all its hydrogens; {H: qH} only the hydrogens.",
             "# Hydrogens can also be given one by one by name, e.g. H6: 0.2042.", ""]
    roles = {"F": "front unit", "R": "repeat unit", "E": "end unit"}
    for resnum in sorted(by_res):
        idx = by_res[resnum]
        res = mol.GetAtomWithIdx(idx[0]).GetProp("res")
        if res in seen:
            continue
        seen.add(res)
        if resnum in unit_of:
            desc = f"{roles.get(res[-1], 'unit')} of {unit_of[resnum]}"
        else:
            desc = "initiator end group"
        lines.append(f"# {res}: {desc}.")
        lines.append(f"{res}:")
        for i in idx:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() == 1:
                continue  # hydrogens are written with their parent atom
            q = topology.charges[i] if topology is not None else 0.0
            hs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
            heavy_nbrs = " ".join(_neighbour_label(mol, n.GetIdx(), resnum)
                                  for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)
            name = atom.GetProp("name")
            if hs:
                qh = sum(topology.charges[h] for h in hs) / len(hs) if topology is not None else 0.0
                h_names = " ".join(mol.GetAtomWithIdx(h).GetProp("name") for h in hs)
                lines.append(f"  # {name}: {{charge: {q:.4f}, H: {qh:.4f}}}   # {atom.GetSymbol()}, {len(hs)} H "
                             f"({h_names}), bonded to {heavy_nbrs}")
            else:
                lines.append(f"  # {name}: {q:.4f}   # {atom.GetSymbol()}, bonded to {heavy_nbrs}")
        lines.append("")
    return "\n".join(lines)
