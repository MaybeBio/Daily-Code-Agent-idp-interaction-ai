"""PDB writer for PolymerStructure (GROMACS-friendly atom order, CONECT, CRYST1)"""

# rq-18964b06

from __future__ import annotations

import os
from pathlib import Path

from ._build import PolymerStructure
from ._errors import PdbLimitError

MAX_ATOMS = 99_999
MAX_RESIDUES = 9_999


def _atom_name_field(name: str, element: str) -> str:
    # PDB atom-name columns
    if len(element) == 1 and len(name) < 4:
        return f" {name:<3s}"
    return f"{name:<4s}"


def _charge_field(q: int) -> str:
    return "  " if q == 0 else f"{abs(q)}{'+' if q > 0 else '-'}"


# rq-671d8f4c
def write_pdb(structure: PolymerStructure, path: str | os.PathLike[str], padding: float = 10.0) -> Path:
    """Write the structure as PDB (ATOM, CONECT, CRYST1, TER, END) and return the path"""
    mol = structure.mol
    n_atoms = mol.GetNumAtoms()
    if n_atoms > MAX_ATOMS:
        raise PdbLimitError(f"{n_atoms} atoms exceed the PDB limit of {MAX_ATOMS}.")
    residues = {a.GetIntProp("resnum") for a in mol.GetAtoms()}
    if len(residues) > MAX_RESIDUES:
        raise PdbLimitError(f"{len(residues)} residues exceed the PDB limit of {MAX_RESIDUES}.")
    if padding < 0:
        raise ValueError("padding must be >= 0")

    xyz = mol.GetConformer().GetPositions()
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    box = (hi - lo) + 2 * padding
    xyz = xyz - lo + padding

    lines = [f"CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1"]
    for i, atom in enumerate(mol.GetAtoms()):
        x, y, z = xyz[i]
        element = atom.GetSymbol()
        lines.append(
            f"ATOM  {i + 1:5d} {_atom_name_field(atom.GetProp('name'), element)} "
            f"{atom.GetProp('res'):>3s} A{atom.GetIntProp('resnum'):4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          "
            f"{element.upper():>2s}{_charge_field(atom.GetFormalCharge())}"
        )
    lines.append("TER")
    for atom in mol.GetAtoms():
        partners = sorted(nb.GetIdx() + 1 for nb in atom.GetNeighbors())
        for start in range(0, len(partners), 4):
            chunk = "".join(f"{p:5d}" for p in partners[start:start + 4])
            lines.append(f"CONECT{atom.GetIdx() + 1:5d}{chunk}")
    lines.append("END")

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="ascii")
    return out


__all__ = ["write_pdb"]
