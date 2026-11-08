"""Data types shared by force-field backends and writers"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, runtime_checkable

from rdkit import Chem

from ..structure import PolymerStructure


# rq-a68ca1c4
@dataclass(frozen=True)
class AtomType:
    name: str
    atomic_number: int
    mass: float  # u
    sigma: float  # nm
    epsilon: float  # kJ/mol
    bond_type: str = ""  # bonded-parameter class (OPLS 'CT'); empty means the type name

    @property
    def bonded_key(self) -> str:
        return self.bond_type or self.name


# rq-840650f8
@dataclass(frozen=True)
class ParameterSet:
    """GROMACS-unit parameters"""

    defaults: tuple[int, int, str, float, float]
    rtp_bondedtypes: tuple[int, ...]
    atomtypes: Mapping[str, AtomType] = field(default_factory=dict)
    bondtypes: Mapping[tuple[str, str], tuple[float, float]] = field(default_factory=dict)
    angletypes: Mapping[tuple[str, str, str], tuple[float, float]] = field(default_factory=dict)
    dihedraltypes: Mapping[tuple[str, str, str, str], tuple[tuple[float, float, int], ...]] = field(
        default_factory=dict)
    impropertypes: Mapping[tuple[str, str, str, str], tuple[float, float, int]] = field(default_factory=dict)
    dihedral_funct: int = 9  # 9: Fourier terms (phase, k, mult); 3: six Ryckaert-Bellemans coefficients


# rq-81766163
@dataclass(frozen=True)
class FragmentParameters:
    atom_types: tuple[str, ...]
    charges: tuple[float, ...]
    parameters: ParameterSet


# rq-9616f2d1
@dataclass(frozen=True)
class ResidueBlock:
    name: str
    atom_names: tuple[str, ...]
    atom_types: tuple[str, ...]
    charges: tuple[float, ...]
    bonds: tuple[tuple[str, str], ...]
    impropers: tuple[tuple[str, str, str, str], ...]
    charge: float  # formal charge x charge_scale


# rq-853c44c5
@dataclass(frozen=True)
class ChargeReport:
    """How user charges were balanced in one residue block"""

    residue: str
    user_atoms: int
    correction_per_atom: float
    large: bool


@dataclass(frozen=True, eq=False)
class Topology:
    name: str
    force_field: str
    structure: PolymerStructure
    blocks: tuple[ResidueBlock, ...]
    atom_types: tuple[str, ...]
    charges: tuple[float, ...]
    parameters: ParameterSet
    bonds: tuple[tuple[int, int], ...]
    pairs: tuple[tuple[int, int], ...]
    angles: tuple[tuple[int, int, int], ...]
    propers: tuple[tuple[int, int, int, int], ...]
    impropers: tuple[tuple[int, int, int, int], ...]
    net_charge: int
    charge_report: tuple[ChargeReport, ...] = ()
    user_parameters: tuple[tuple[str, tuple[str, ...]], ...] = ()
    missing_parameters: Mapping[str, tuple] = field(default_factory=dict)
    dft_blocks: tuple[str, ...] = ()  # residue blocks whose charges came from DFT
    dft_missing: tuple[str, ...] = ()  # molecules (monomer / end group) without DFT charges
    ion_blocks: tuple[str, ...] = ()  # residues whose ionic atoms carry DFT ion-fragment charges
    charge_scale: float = 1.0  # ionic groups scaled by this factor (rq-12060285)
    builtin_blocks: tuple[str, ...] = ()  # residues with built-in RESP charges
    explicit: Mapping = field(default_factory=dict)  # atom indices -> values written on that bonded line


# rq-54f83eba
@runtime_checkable
class ForceFieldBackend(Protocol):
    name: str

    def check_available(self) -> None: ...

    def parametrize(self, mol: Chem.Mol, net_charge: int, work_dir: Path) -> FragmentParameters: ...
    # optional keyword: charges=False (types/parameters only)
