"""RDKit helpers"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

from ._errors import (
    AmbiguousPolymerizableGroupError,
    ChargeError,
    IncompatibleBlocksError,
    InvalidStructureError,
)
from ._model import Block, BlockCategory

# terminal vinyl C=C
_TERMINAL_VINYL = Chem.MolFromSmarts("[CH2]=[#6]")


def parse_smiles(smiles: str) -> Chem.Mol:
    """Parse and sanitize a SMILES string, raising InvalidStructureError on failure"""
    if not isinstance(smiles, str) or not smiles.strip():
        raise InvalidStructureError(str(smiles), "empty SMILES")
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidStructureError(smiles, "RDKit could not parse or sanitize the SMILES")
    return mol


def parse_smarts(smarts: str) -> Chem.Mol:
    with rdBase.BlockLogs():
        patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise InvalidStructureError(smarts, "RDKit could not parse the SMARTS")
    return patt


def net_charge(mol: Chem.Mol) -> int:
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def check_zwitterion(mol: Chem.Mol, smiles: str) -> None:
    """Require net charge 0 with at least one positive and one negative centre"""
    charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
    total = sum(charges)
    if total != 0:
        raise ChargeError(smiles, total, f"net formal charge is {total:+d}, expected 0")
    if not any(c > 0 for c in charges) or not any(c < 0 for c in charges):
        raise ChargeError(smiles, 0, "a zwitterion needs at least one positive and one negative centre")


# rq-ddcf8658
def check_comonomer(mol: Chem.Mol, smiles: str, counterion: str | None) -> None:
    """A non-zwitterionic comonomer"""
    total = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    if total == 0 and counterion is not None:
        raise ChargeError(smiles, 0, f"a neutral comonomer takes no counterion (got {counterion})")
    if total != 0:
        if counterion not in COMONOMER_COUNTERIONS:
            raise ChargeError(smiles, total, f"net charge {total:+d} needs a counterion "
                                             f"({', '.join(COMONOMER_COUNTERIONS)})")
        if (total > 0) == counterion.endswith("+"):
            raise ChargeError(smiles, total, f"counterion {counterion} has the same sign as the charge {total:+d}")


COMONOMER_COUNTERIONS = ("Cl-", "Br-", "Na+", "K+")


def block_attachment_labels(mol: Chem.Mol) -> list[int]:
    return sorted(a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomicNum() == 0)


def assemble(blocks: Sequence[Block]) -> Chem.Mol:
    """Join blocks in order; block i's [*:2] bonds to block i+1's [*:1]"""
    keys = [b.key for b in blocks]
    if not blocks:
        raise IncompatibleBlocksError(keys, 0, "no blocks given")
    if blocks[0].category is not BlockCategory.POLYMERIZABLE:
        raise IncompatibleBlocksError(keys, 0, f"'{keys[0]}' is not a polymerizable block")
    last = len(blocks) - 1
    for i, block in enumerate(blocks):
        if i > 0 and not block.has_inner:
            raise IncompatibleBlocksError(keys, i, f"'{block.key}' has no backbone-side attachment")
        if i < last and not block.has_outer:
            raise IncompatibleBlocksError(keys, i, f"'{block.key}' is terminal but is not last")
        if i == last and block.has_outer:
            raise IncompatibleBlocksError(keys, i, f"'{block.key}' is not terminal but is last")
    pieces = []
    for i, block in enumerate(blocks):
        mol = Chem.MolFromSmiles(block.smiles)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                # [*:1] -> i, [*:2] -> i+1
                atom.SetAtomMapNum(i if atom.GetAtomMapNum() == 1 else i + 1)
        pieces.append(mol)
    combined = pieces[0]
    for piece in pieces[1:]:
        combined = Chem.CombineMols(combined, piece)
    zipped = Chem.molzip(combined)
    if any(a.GetAtomicNum() == 0 for a in zipped.GetAtoms()):
        raise IncompatibleBlocksError(keys, last, "unmatched attachment points remain")
    Chem.SanitizeMol(zipped)
    return zipped


def canonical(mol: Chem.Mol) -> tuple[Chem.Mol, str]:
    """Return (molecule renumbered to canonical SMILES order, canonical SMILES)"""
    smiles = Chem.MolToSmiles(mol)
    order = list(map(int, mol.GetProp("_smilesAtomOutputOrder").strip("[],").split(",")))
    return Chem.RenumberAtoms(mol, order), smiles


def matching_keys(mol: Chem.Mol, blocks: Iterable[Block]) -> tuple[str, ...]:
    return tuple(b.key for b in blocks if b.smarts and mol.HasSubstructMatch(parse_smarts(b.smarts)))


def find_polymerizable_group(
    mol: Chem.Mol,
    smiles: str,
    polymerizable_blocks: Sequence[Block],
    vinyl_atoms: tuple[int, int] | None = None,
) -> tuple[str | None, tuple[int, int] | None]:
    """Find the polymerizable C=C group."""
    candidates = sorted({(m[0], m[1]) for m in mol.GetSubstructMatches(_TERMINAL_VINYL)})
    if vinyl_atoms is not None:
        pair = tuple(int(i) for i in vinyl_atoms)
        if pair not in candidates and pair[::-1] not in candidates:
            raise InvalidStructureError(smiles, f"atoms {pair} are not a terminal C=C group")
        chosen = pair if pair in candidates else pair[::-1]
    elif not candidates:
        return None, None
    elif len(candidates) == 1:
        chosen = candidates[0]
    else:
        chosen = _single_group_covering(mol, candidates, polymerizable_blocks) or _equivalent(mol, candidates)
        if chosen is None:
            raise AmbiguousPolymerizableGroupError(smiles, candidates)
    for block in polymerizable_blocks:
        for match in mol.GetSubstructMatches(parse_smarts(block.smarts)):
            if chosen[0] in match and chosen[1] in match:
                return block.key, chosen
    return None, None


def _equivalent(mol, candidates):
    """The first C=C if all candidates are symmetry-equivalent (e.g."""
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    if len({(ranks[a], ranks[b]) for a, b in candidates}) == 1:
        return candidates[0]
    return None


def _single_group_covering(mol, candidates, polymerizable_blocks):
    """If one block match contains every candidate C=C (e.g."""
    for block in polymerizable_blocks:
        for match in mol.GetSubstructMatches(parse_smarts(block.smarts)):
            atoms = set(match)
            if all(a in atoms and b in atoms for a, b in candidates):
                return candidates[0]
    return None


def formula_and_weight(mol: Chem.Mol) -> tuple[str, float]:
    return rdMolDescriptors.CalcMolFormula(mol), round(Descriptors.MolWt(mol), 3)
