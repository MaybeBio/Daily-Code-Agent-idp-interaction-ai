"""Shared multi-chain assembly"""

# rq-b2c69832 rq-88cb8d61

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from ..forcefield import Topology
from ..forcefield._gromacs import forcefield_files, moleculetype_text
from ..structure import PolymerStructure
from ._errors import SystemBuildError
from ._pack import Piece
from ._solvent import (
    COUNTERIONS,
    Molecule,
    atomtype_lines,
    check_defaults,
    locate_solvent_ff,
    madrid2019,
    read_molecules,
    scaled,
    type_mass,
)

AMU_G = 1.66053906660e-24  # g per atomic mass unit
NAME = "polymer"


@dataclass
class SolventParts:
    folder: Path
    molecules: dict[str, Molecule]  # name -> block, for every molecule the system uses
    type_lines: dict[str, str]  # atom type -> [ atomtypes ] line to append
    masses: dict[str, float]  # molecule name -> total mass (u)
    pairs: str = ""  # [ nonbond_params ] lines (Madrid-2019)


def counterions_for(structure: PolymerStructure, chains: int) -> tuple[str | None, int]:
    """(GROMACS ion molecule, count) that neutralizes `chains` copies of the chain"""
    q = structure.chain.net_charge
    if q == 0:
        return None, 0
    ion = structure.chain.counterion
    if ion not in COUNTERIONS:
        raise SystemBuildError(f"Counterion {ion!r} is not supported (supported: {', '.join(COUNTERIONS)}); "
                               "the chain's net charge cannot be neutralized.")
    if (q > 0) == ion.endswith("+"):
        raise SystemBuildError(f"Counterion {ion} has the same sign as the chain charge {q:+d}.")
    return COUNTERIONS[ion], abs(q) * chains


def solvent_parts(topology: Topology, molecule_names: Sequence[str], solvent_ff=None,
                  water_file: str | None = None, solvent: str = "tip3p") -> SolventParts:
    if solvent == "tip4p2005":  # rq-96c58164: built-in TIP4P/2005 + Madrid-2019 (charges already +-0.85)
        with_bt = any(at.bonded_key != at.name for at in topology.parameters.atomtypes.values())
        molecules, lines, pairs = madrid2019(molecule_names, with_bt)
        masses = {n: sum(m.masses) for n, m in molecules.items()}
        return SolventParts(Path("madrid2019"), molecules, lines, masses, pairs)
    if solvent != "tip3p":
        raise ValueError(f"solvent must be 'tip3p' or 'tip4p2005', got {solvent!r}")
    folder = locate_solvent_ff(topology.parameters.defaults[1], solvent_ff)
    check_defaults(folder, topology.parameters.defaults)
    molecules = read_molecules(folder / "ions.itp")
    if water_file is not None:
        if not (folder / water_file).is_file():
            raise SystemBuildError(f"{folder} has no {water_file}")
        molecules.update(read_molecules(folder / water_file))
    chosen = {}
    for name in molecule_names:
        if name not in molecules:
            raise SystemBuildError(f"{folder} defines no molecule {name}")
        chosen[name] = molecules[name]
        if topology.charge_scale != 1 and len(chosen[name].atom_types) == 1:  # ions, not water
            chosen[name] = scaled(chosen[name], topology.charge_scale)
    types = {t for m in chosen.values() for t in m.atom_types}
    lines = atomtype_lines(folder, types) if types else {}
    clash = sorted(t for t in lines if t in topology.parameters.atomtypes)
    if clash:
        raise SystemBuildError(f"Water/ion atom type(s) {', '.join(clash)} clash with polymer atom types.")
    masses = {}
    for name, m in chosen.items():
        masses[name] = sum(mass if mass > 0 else type_mass(lines[t]) for t, mass in zip(m.atom_types, m.masses,
                                                                                         strict=True))
    return SolventParts(folder, chosen, lines, masses)


def box_angstrom(box_nm) -> tuple[float, float, float]:
    """User box in nm (one edge for a cube, or three edges) -> Å"""
    edges = (box_nm,) * 3 if isinstance(box_nm, (int, float)) else tuple(box_nm)
    if len(edges) != 3 or not all(isinstance(e, (int, float)) and e > 0 for e in edges):
        raise ValueError(f"box_nm must be a positive edge or three positive edges (nm), got {box_nm!r}")
    return tuple(float(e) * 10.0 for e in edges)


def density_g_cm3(mass_u: float, volume_A3: float) -> float:
    return mass_u * AMU_G / (volume_A3 * 1e-24)


def chain_mass(topology: Topology) -> float:
    return float(sum(topology.parameters.atomtypes[t].mass for t in topology.atom_types))


def chain_piece(structure: PolymerStructure, emc: bool = False) -> Piece:
    """The chain as a packing piece"""
    mol = structure.mol
    piece = Piece(mol.GetConformer().GetPositions().copy(), [a.GetSymbol() for a in mol.GetAtoms()])
    if emc:
        piece.pivots, piece.near = _pivots_and_near(mol)
    return piece


def _pivots_and_near(mol):
    nbrs = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {a: k for k, a in enumerate(heavy)}
    near = np.zeros((len(heavy), len(heavy)), dtype=bool)
    for a in heavy:  # heavy atoms within 3 bonds
        frontier, seen = {a}, {a}
        for _ in range(3):
            frontier = {n for f in frontier for n in nbrs[f] if n not in seen}
            seen |= frontier
        for b in seen:
            if b in pos:
                near[pos[a], pos[b]] = True
    pivots = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        heavy_deg = [sum(mol.GetAtomWithIdx(n).GetAtomicNum() > 1 for n in nbrs[i]) for i in (a, b)]
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing() or min(heavy_deg) < 2:
            continue
        side, todo = {b}, [b]  # atoms on b's side of the bond
        while todo:
            for n in nbrs[todo.pop()]:
                if n not in side and n != a:
                    side.add(n)
                    todo.append(n)
        if 2 * len(side) > mol.GetNumAtoms():  # rotate the smaller side
            side = set(range(mol.GetNumAtoms())) - side
            a, b = b, a
        pivots.append((a, b, np.array(sorted(side - {b}))))
    return tuple(pivots), near


def ion_piece(molecule: Molecule) -> Piece:
    """A one-atom ion; it counts as a heavy atom for the distance checks"""
    return Piece(np.zeros((1, 3)), [molecule.atom_names[0]])


def write_forcefield(out: Path, topology: Topology, parts: SolventParts | None) -> None:
    files = forcefield_files(topology, NAME)
    if parts is not None and parts.type_lines:
        files[f"{NAME}.ff/ffnonbonded.itp"] += ("\n; water and ion types from " + parts.folder.name + "\n"
                                                + "\n".join(parts.type_lines[t] for t in parts.type_lines) + "\n")
    if parts is not None and parts.pairs:
        files[f"{NAME}.ff/ffnonbonded.itp"] += ("\n[ nonbond_params ]\n; i  j  func  sigma (nm)  epsilon (kJ/mol)\n"
                                                + parts.pairs + "\n")
    if parts is not None:
        ions = [m for m in parts.molecules.values() if len(m.atom_types) == 1]
        if ions:
            files[f"{NAME}.ff/ions.itp"] = f"; ions from {parts.folder.name}\n\n" + "\n".join(m.text for m in ions)
        water = [m for m in parts.molecules.values() if len(m.atom_types) > 1]
        if water:
            files[f"{NAME}.ff/water.itp"] = f"; water from {parts.folder.name}\n\n" + "\n".join(m.text for m in water)
    (out / f"{NAME}.ff").mkdir(parents=True, exist_ok=True)
    for rel, text in files.items():
        (out / rel).write_text(text, encoding="ascii")
    (out / f"{NAME}.itp").write_text("; polymer chain (generated by Zpoly)\n\n"
                                     + moleculetype_text(topology, NAME), encoding="ascii")


def write_system_top(out: Path, title: str, chains: int, molecules: Sequence[tuple[str, int]],
                     has_ions: bool, has_water: bool, charge_scale: float = 1.0) -> None:
    lines = [f"; {title} (generated by Zpoly)"]
    if charge_scale != 1:
        lines.append(f"; ionic groups and ions: charges scaled by {charge_scale:g}")
    lines += [f'#include "{NAME}.ff/forcefield.itp"',
             f'#include "{NAME}.itp"']
    if has_water:
        lines.append(f'#include "{NAME}.ff/water.itp"')
    if has_ions:
        lines.append(f'#include "{NAME}.ff/ions.itp"')
    lines += ["", "[ system ]", f"  {title}", "", "[ molecules ]", "; compound   nmols", f"  {NAME:<12s} {chains}"]
    lines += [f"  {name:<12s} {count}" for name, count in molecules if count]
    (out / "system.top").write_text("\n".join(lines) + "\n", encoding="ascii")


def write_gro(path: Path, title: str, structure: PolymerStructure, chain_xyz: Sequence[np.ndarray],
              extra: Sequence[tuple[Molecule, np.ndarray]], box_A: Sequence[float]) -> int:
    """Write a .gro (nm); returns the number of atoms."""
    mol = structure.mol
    atoms = list(mol.GetAtoms())
    nres = max(a.GetIntProp("resnum") for a in atoms)
    rows = []
    for k, xyz in enumerate(chain_xyz):
        for a, (x, y, z) in zip(atoms, xyz, strict=True):
            rows.append((k * nres + a.GetIntProp("resnum"), a.GetProp("res"), a.GetProp("name"), x, y, z))
    resnum = len(chain_xyz) * nres
    for molecule, xyz in extra:
        resnum += 1
        for name, (x, y, z) in zip(molecule.atom_names, xyz, strict=True):
            rows.append((resnum, molecule.residue, name, x, y, z))
    out = [title, f"{len(rows):5d}"]
    for i, (rn, res, name, x, y, z) in enumerate(rows, 1):
        out.append(f"{rn % 100000:5d}{res:<5s}{name:>5s}{i % 100000:5d}{x / 10:8.3f}{y / 10:8.3f}{z / 10:8.3f}")
    out.append("".join(f"{b / 10:10.5f}" for b in box_A))
    path.write_text("\n".join(out) + "\n", encoding="ascii")
    return len(rows)


def total_charge(topology: Topology, chains: int, parts: SolventParts | None, counts: dict[str, int]) -> float:
    q = sum(topology.charges) * chains
    if parts is not None:
        q += sum(sum(parts.molecules[n].charges) * c for n, c in counts.items())
    return q
