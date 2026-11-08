"""Ionic groups and their capped fragments."""

# rq-d4bd3142

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem

CAP_BOND = 1.09  # Å, C-H / X-H cap length
_CENTRAL = {"S", "P", "C"}
_TERMINAL = {"O", "S"}


@dataclass(frozen=True)
class IonGroup:
    atoms: tuple[int, ...]  # heavy atoms and their hydrogens, sorted
    charge: int


def _core(mol: Chem.Mol, atom: Chem.Atom) -> set[int]:
    if atom.GetIsAromatic():
        ri = mol.GetRingInfo()
        rings = [set(r) for r in ri.AtomRings() if atom.GetIdx() in r and all(
            mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
        system = set().union(*rings) if rings else {atom.GetIdx()}
        grown = True
        while grown:  # fused aromatic rings
            grown = False
            for r in ri.AtomRings():
                if set(r) & system and not set(r) <= system and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
                    system |= set(r)
                    grown = True
        return system
    if atom.GetFormalCharge() < 0 and atom.GetSymbol() in _TERMINAL and atom.GetDegree() == 1:
        centre = atom.GetNeighbors()[0]
        if centre.GetSymbol() in _CENTRAL:
            terminal = {n.GetIdx() for n in centre.GetNeighbors()
                        if n.GetSymbol() in _TERMINAL and n.GetDegree() == 1}
            return {centre.GetIdx()} | terminal
    return {atom.GetIdx()}


def ionic_groups(mol: Chem.Mol) -> list[IonGroup]:
    """Charged groups of a molecule"""
    groups: list[tuple[set[int], int]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1 or atom.GetFormalCharge() == 0:
            continue
        core = _core(mol, atom)
        shell: set[int] = set()
        for i in core:
            for n in mol.GetAtomWithIdx(i).GetNeighbors():
                if n.GetAtomicNum() > 1 and n.GetIdx() not in core:
                    shell.add(n.GetIdx())
        for i in list(shell):
            a = mol.GetAtomWithIdx(i)
            if a.GetSymbol() != "C":
                shell |= {n.GetIdx() for n in a.GetNeighbors()
                          if n.GetSymbol() == "C" and n.GetIdx() not in core}
        heavy = core | shell
        heavy |= {n.GetIdx() for i in shell for n in mol.GetAtomWithIdx(i).GetNeighbors()
                  if n.GetAtomicNum() > 1 and n.GetDegree() == 1}  # =O, F, ... stay on their atom
        groups.append(heavy)
    # merge groups that overlap, or same-sign groups that touch
    merged: list[set[int]] = []
    for heavy in groups:
        heavy = set(heavy)
        hit = [m for m in merged if m & heavy or (_sign(mol, m) == _sign(mol, heavy) and _bonded(mol, m, heavy))]
        for m in hit:
            merged.remove(m)
            heavy |= m
        merged.append(heavy)
    out = []
    for heavy in merged:
        charge = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in heavy)
        hs = {n.GetIdx() for i in heavy for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1}
        out.append(IonGroup(tuple(sorted(heavy | hs)), charge))
    return sorted(out, key=lambda g: (-g.charge, g.atoms))


def _sign(mol: Chem.Mol, atoms: set[int]) -> int:
    q = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in atoms)
    return (q > 0) - (q < 0)


def _bonded(mol: Chem.Mol, a: set[int], b: set[int]) -> bool:
    return any(n.GetIdx() in b for i in a for n in mol.GetAtomWithIdx(i).GetNeighbors())


@dataclass(frozen=True)
class Fragment:
    mol: Chem.Mol  # with explicit H and one conformer
    parent_atoms: tuple[int, ...]  # parent index per fragment atom, -1 for a cap
    caps: tuple[tuple[int, int], ...]  # (cap atom, fragment atom it is bonded to)
    charge: int


def fragment(mol: Chem.Mol, positions: np.ndarray, group: IonGroup) -> Fragment:
    """Cut a group out of a molecule, capping cut bonds with H."""
    keep = list(group.atoms)
    index = {p: k for k, p in enumerate(keep)}
    rw = Chem.RWMol()
    for p in keep:
        a = mol.GetAtomWithIdx(p)
        new = Chem.Atom(a.GetAtomicNum())
        new.SetFormalCharge(a.GetFormalCharge())
        new.SetIsAromatic(a.GetIsAromatic())
        new.SetNoImplicit(True)
        rw.AddAtom(new)
    xyz = [positions[p] for p in keep]
    caps = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in index and b in index:
            rw.AddBond(index[a], index[b], bond.GetBondType())
        elif (a in index) != (b in index):
            inside, outside = (a, b) if a in index else (b, a)
            h = rw.AddAtom(Chem.Atom(1))
            rw.AddBond(index[inside], h, Chem.BondType.SINGLE)
            direction = positions[outside] - positions[inside]
            xyz.append(positions[inside] + direction / np.linalg.norm(direction) * CAP_BOND)
            caps.append((h, index[inside]))
    frag = rw.GetMol()
    Chem.SanitizeMol(frag)
    conf = Chem.Conformer(frag.GetNumAtoms())
    for k, (x, y, z) in enumerate(xyz):
        conf.SetAtomPosition(k, (float(x), float(y), float(z)))
    frag.AddConformer(conf, assignId=True)
    parents = tuple(keep) + (-1,) * len(caps)
    return Fragment(frag, parents, tuple(caps), group.charge)


_LINKER = Chem.MolFromSmarts("[CX4:1][CX3:2](=[OX1:3])[OX2,NX3:4][#6:5]")


def linker_group(mol: Chem.Mol, backbone_carbon: int | None) -> tuple[IonGroup, tuple[int, ...], str] | None:
    """Ester/amide linker group for a DFT fragment."""
    if backbone_carbon is None:
        return None
    for alpha, carbonyl, oxo, x, side in mol.GetSubstructMatches(_LINKER):
        if alpha != backbone_carbon:
            continue
        core = {carbonyl, oxo, x, side}
        backbone = {alpha} | {n.GetIdx() for n in mol.GetAtomWithIdx(alpha).GetNeighbors()
                              if n.GetSymbol() == "C" and n.GetIdx() != carbonyl}
        hs = lambda atoms: {n.GetIdx() for i in atoms for n in mol.GetAtomWithIdx(i).GetNeighbors()  # noqa: E731
                            if n.GetAtomicNum() == 1}
        apply = core | hs({x, side})
        kind = "amide" if mol.GetAtomWithIdx(x).GetSymbol() == "N" else "ester"
        return IonGroup(tuple(sorted(core | backbone | hs(core | backbone))), 0), tuple(sorted(apply)), kind
    return None
