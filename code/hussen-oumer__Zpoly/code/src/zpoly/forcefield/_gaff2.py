"""GAFF2 backend"""

# rq-85e6069f rq-c7ea04ff

from __future__ import annotations

import functools
import os
import re
import shutil
import subprocess
from dataclasses import replace
from pathlib import Path

from rdkit import Chem

from ._amber import read_amber_parameters
from ._errors import ForceFieldUnavailableError, ParametrizationError
from ._model import FragmentParameters, ParameterSet

INSTALL_HINT = "Install AmberTools with: conda install -c conda-forge ambertools"
_ELEMENT = re.compile(r"[A-Za-z]+")


@functools.lru_cache(maxsize=4)
def _cached_dat(path: str, mtime: float) -> ParameterSet:
    return read_amber_parameters(path)


# rq-a54a6990
def read_mol2_atoms(path: Path) -> list[tuple[str, str, float]]:
    """(atom name, atom type, charge) for each atom of a Tripos mol2 file"""
    atoms, in_atoms = [], False
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            continue
        if in_atoms and line.strip():
            tok = line.split()
            atoms.append((tok[1], tok[5], float(tok[8])))
    return atoms


def mol2_types_and_charges(path: Path, mol: Chem.Mol) -> tuple[tuple[str, ...], tuple[float, ...]]:
    """Types and charges from antechamber's mol2."""
    atoms = read_mol2_atoms(path)
    if len(atoms) != mol.GetNumAtoms():
        raise ParametrizationError("antechamber", 0, "", f"antechamber returned {len(atoms)} atoms, "
                                   f"expected {mol.GetNumAtoms()}")
    for i, (aname, _, _) in enumerate(atoms):
        symbol = mol.GetAtomWithIdx(i).GetSymbol()
        if not _element_matches(aname, symbol):
            raise ParametrizationError("antechamber", 0, "", f"atom {i + 1} is '{aname}' in the mol2 "
                                       f"file but {symbol} in the input")
    return tuple(t for _, t, _ in atoms), tuple(q for _, _, q in atoms)


def _element_matches(atom_name: str, symbol: str) -> bool:
    letters = _ELEMENT.match(atom_name)
    if not letters:
        return False
    text = letters.group(0)
    return text[:2].casefold() == symbol.casefold() or text[:1].casefold() == symbol.casefold()


# rq-2702f408 rq-1b967d27 rq-8c905ad1
class Gaff2Backend:
    """GAFF2 atom types and AM1-BCC charges via AmberTools"""

    name = "gaff2"

    def __init__(self, timeout: float = 3600.0, antechamber: str = "antechamber", parmchk2: str = "parmchk2",
                 gaff2_dat: str | os.PathLike | None = None) -> None:
        self.timeout = timeout
        self.antechamber = antechamber
        self.parmchk2 = parmchk2
        self.gaff2_dat = gaff2_dat

    def _locate(self) -> tuple[str, str, Path]:
        ac = shutil.which(self.antechamber)
        if ac is None:
            raise ForceFieldUnavailableError(self.name, f"'{self.antechamber}' was not found on PATH", INSTALL_HINT)
        pc = shutil.which(self.parmchk2)
        if pc is None:
            raise ForceFieldUnavailableError(self.name, f"'{self.parmchk2}' was not found on PATH", INSTALL_HINT)
        if self.gaff2_dat is not None:
            dat = Path(self.gaff2_dat)
        elif os.environ.get("AMBERHOME"):
            dat = Path(os.environ["AMBERHOME"]) / "dat" / "leap" / "parm" / "gaff2.dat"
        else:
            dat = Path(ac).resolve().parent.parent / "dat" / "leap" / "parm" / "gaff2.dat"
        if not dat.is_file():
            raise ForceFieldUnavailableError(
                self.name, f"gaff2.dat was not found at {dat}",
                f"Set AMBERHOME to your AmberTools installation. {INSTALL_HINT}")
        return ac, pc, dat

    def check_available(self) -> None:
        self._locate()

    def full_parameters(self) -> ParameterSet:
        """Every parameter of gaff2.dat (for same_as lookups in parameter files)."""
        _, _, dat = self._locate()
        return _cached_dat(str(dat), dat.stat().st_mtime)

    def _run(self, tool: str, argv: list[str], cwd: Path, expected: Path) -> None:
        try:
            proc = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, timeout=self.timeout)
        except subprocess.TimeoutExpired as exc:
            out = (exc.stdout or "") + (exc.stderr or "")
            out = out.decode(errors="replace") if isinstance(out, bytes) else out
            raise ParametrizationError(tool, -1, out, f"{tool} timed out after {self.timeout:g} s") from None
        log = (proc.stdout or "") + (proc.stderr or "")
        if proc.returncode != 0:
            raise ParametrizationError(tool, proc.returncode, log)
        if not expected.is_file():
            raise ParametrizationError(tool, proc.returncode, log, f"{tool} did not write {expected.name}")

    def parametrize(self, mol: Chem.Mol, net_charge: int, work_dir: Path, *, charges: bool = True
                    ) -> FragmentParameters:
        ac, pc, dat = self._locate()
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        sdf, mol2, frcmod = work_dir / "frag.sdf", work_dir / "frag.mol2", work_dir / "frag.frcmod"
        Chem.MolToMolFile(mol, str(sdf), kekulize=True)
        # rq-8a4726a9: no AM1-BCC with DFT charges
        charge_args = ["-c", "bcc", "-nc", str(net_charge)] if charges else []
        self._run("antechamber", [ac, "-i", sdf.name, "-fi", "mdl", "-o", mol2.name, "-fo", "mol2",
                                  *charge_args, "-at", "gaff2", "-pf", "y", "-s", "0"],
                  work_dir, mol2)
        self._run("parmchk2", [pc, "-i", mol2.name, "-f", "mol2", "-o", frcmod.name, "-s", "gaff2"],
                  work_dir, frcmod)
        types, mol2_charges = mol2_types_and_charges(mol2, mol)
        params = self._parameters(dat, frcmod, mol, types)
        return FragmentParameters(types, mol2_charges if charges else tuple(0.0 for _ in types), params)

    def _parameters(self, dat: Path, frcmod: Path, mol: Chem.Mol, types: tuple[str, ...]) -> ParameterSet:
        base = _cached_dat(str(dat), dat.stat().st_mtime)
        mod = read_amber_parameters(frcmod, frcmod=True)
        present = set(types)

        def keep(key) -> bool:
            return all(t in present or t == "X" for t in key)

        def combine(attr):
            merged = {**getattr(base, attr), **getattr(mod, attr)}  # frcmod wins
            return {k: v for k, v in merged.items() if keep(k)}

        atomic = {}
        for atom, t in zip(mol.GetAtoms(), types, strict=True):
            atomic.setdefault(t, atom.GetAtomicNum())
        atomtypes = {}
        for t in sorted(present):
            at = mod.atomtypes.get(t) or base.atomtypes.get(t)
            if at is None:
                raise ParametrizationError("parmchk2", 0, "", f"no mass/Lennard-Jones parameters for GAFF2 type {t}")
            atomtypes[t] = replace(at, atomic_number=atomic[t])
        return ParameterSet(base.defaults, base.rtp_bondedtypes, atomtypes, combine("bondtypes"),
                            combine("angletypes"), combine("dihedraltypes"), combine("impropertypes"))
