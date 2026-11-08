"""Data types of the zwitterion library"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType

from rdkit import Chem

from ._errors import NoPolymerizableGroupError


def _norm_token(value: str) -> str:
    return value.replace("-", "").replace("_", "").replace(" ", "").casefold()


# rq-a09ee60d
class ForceField(str, Enum):
    """Force fields known to the project"""

    OPLSAA = "oplsaa"
    GAFF2 = "gaff2"

    @classmethod
    def _missing_(cls, value: object) -> ForceField | None:
        if isinstance(value, str):
            wanted = _norm_token(value)
            for member in cls:
                if member.value == wanted:
                    return member
        return None


# rq-158c925e
class BlockCategory(str, Enum):
    """Role of a building block in a monomer"""

    POLYMERIZABLE = "polymerizable"
    LINKER = "linker"
    CATION = "cation"
    SPACER = "spacer"
    ANION = "anion"


# rq-e896b7cf
@dataclass(frozen=True)
class Block:
    """A building block with attachment points."""

    key: str
    category: BlockCategory
    name: str
    smiles: str
    smarts: str | None
    formal_charge: int
    terminal: bool
    notes: str = ""

    @property
    def has_inner(self) -> bool:
        """True if the block has a [*:1] attachment (toward the backbone)"""
        return "[*:1]" in self.smiles

    @property
    def has_outer(self) -> bool:
        """True if the block has a [*:2] attachment (away from the backbone)"""
        return "[*:2]" in self.smiles


_EMPTY: Mapping = MappingProxyType({})


# rq-5e9d599c
@dataclass(frozen=True)
class Zwitterion:
    """A validated zwitterionic monomer or small-molecule zwitterion"""

    name: str
    smiles: str
    full_name: str = ""
    family: str = ""
    status: str = "user"
    aliases: tuple[str, ...] = ()
    blocks: tuple[str, ...] | None = None
    polymerizable_group: str | None = None
    vinyl_atoms: tuple[int, int] | None = None
    cation_classes: tuple[str, ...] = ()
    anion_classes: tuple[str, ...] = ()
    formula: str = ""
    molecular_weight: float = 0.0
    source: str = ""
    notes: str = ""
    residue: str | None = None
    zwitterionic: bool = True
    counterion: str | None = None
    properties: Mapping[str, float | str] = field(default=_EMPTY)
    parameters: Mapping[ForceField, str] = field(default=_EMPTY)

    def __hash__(self) -> int:
        return hash((self.name, self.smiles))

    @property
    def polymerizable(self) -> bool:
        return self.polymerizable_group is not None

    def to_mol(self, add_hs: bool = True) -> Chem.Mol:
        """Return a new RDKit molecule for this structure"""
        mol = Chem.MolFromSmiles(self.smiles)
        return Chem.AddHs(mol) if add_hs else mol

    def repeat_unit_smiles(self) -> str:
        """Repeat-unit SMILES with [*:1] and [*:2] attachment points."""
        if self.vinyl_atoms is None:
            raise NoPolymerizableGroupError(f"'{self.name}' has no polymerizable group.")
        head, tail = self.vinyl_atoms
        mol = Chem.RWMol(Chem.MolFromSmiles(self.smiles))
        mol.GetBondBetweenAtoms(head, tail).SetBondType(Chem.BondType.SINGLE)
        for atom_idx, label in ((head, 1), (tail, 2)):
            dummy = Chem.Atom(0)
            dummy.SetAtomMapNum(label)
            new_idx = mol.AddAtom(dummy)
            mol.AddBond(atom_idx, new_idx, Chem.BondType.SINGLE)
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)


# rq-4711c4d5 rq-fbb9e221
class Ends(str, Enum):
    """Which chain ends carry initiator end groups"""

    BOTH = "both"
    HEAD = "head"
    TAIL = "tail"
    NONE = "none"


# rq-6114615c
@dataclass(frozen=True)
class Initiator:
    """A radical initiator"""

    name: str
    head_smiles: str | None
    tail_smiles: str | None
    full_name: str = ""
    aliases: tuple[str, ...] = ()
    head_charge: int = 0
    tail_charge: int = 0
    counterion: str | None = None
    source: str = ""
    notes: str = ""
    parameters: Mapping[ForceField, str] = field(default=_EMPTY)

    def __hash__(self) -> int:
        return hash((self.name, self.head_smiles, self.tail_smiles))


# rq-f3793f01
@dataclass(frozen=True)
class Chain:
    """A linear chain of repeat units with optional initiator end groups (2D structure)"""

    smiles: str
    units: tuple[str, ...]
    initiator: str | None
    ends: Ends
    head_group: str | None
    tail_group: str | None
    net_charge: int
    counterion: str | None
    formula: str
    molecular_weight: float
    cap: str = "CH3"  # group on ends without an end group: "CH3" or "H"
