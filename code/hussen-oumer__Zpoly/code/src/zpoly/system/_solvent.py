"""Water and ion parameters."""

# rq-39bee73a

from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from ._errors import SystemBuildError

COUNTERIONS = {"K+": "K", "Na+": "NA", "Cl-": "CL"}  # chain counterion -> GROMACS ion molecule
SOLVENT_FF = {2: "amber99sb-ildn.ff", 3: "oplsaa.ff"}  # combination rule -> matching folder


@dataclass(frozen=True)
class Molecule:
    name: str
    text: str  # the [ moleculetype ] block
    atom_types: tuple[str, ...]
    residue: str
    atom_names: tuple[str, ...]
    charges: tuple[float, ...]
    masses: tuple[float, ...]


def locate_solvent_ff(comb_rule: int, solvent_ff: str | os.PathLike | None = None) -> Path:
    if solvent_ff is not None:
        folder = Path(solvent_ff)
    else:
        name = SOLVENT_FF.get(comb_rule)
        if name is None:
            raise SystemBuildError(f"No matching water/ion force field for combination rule {comb_rule}; "
                                   "pass solvent_ff.")
        candidates = []
        if os.environ.get("GMXDATA"):
            candidates.append(Path(os.environ["GMXDATA"]) / "top" / name)
        for exe in ("gmx", "gmx_mpi"):
            found = shutil.which(exe)
            if found:
                candidates.append(Path(found).resolve().parent.parent / "share" / "gromacs" / "top" / name)
        folder = next((c for c in candidates if (c / "ffnonbonded.itp").is_file()), None)
        if folder is None:
            raise SystemBuildError(f"Cannot find {name} for water and ion parameters; install GROMACS or pass "
                                   "solvent_ff=/path/to/" + name)
    for f in ("forcefield.itp", "ffnonbonded.itp", "ions.itp"):
        if not (folder / f).is_file():
            raise SystemBuildError(f"{folder} is not a usable water/ion force field: {f} is missing")
    return folder


def _active_lines(path: Path):
    """Data lines with #ifdef blocks evaluated."""
    stack: list[bool] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        w = raw.split()[:1]
        if w == ["#ifdef"]:
            stack.append(False)
            continue
        if w == ["#ifndef"]:
            stack.append(True)
            continue
        if w == ["#else"] and stack:
            stack[-1] = not stack[-1]
            continue
        if w == ["#endif"] and stack:
            stack.pop()
            continue
        if all(stack):
            yield raw


def solvent_defaults(folder: Path) -> tuple:
    for raw in _active_lines(folder / "forcefield.itp"):
        data = raw.split(";")[0].split()
        if len(data) >= 5 and data[0].isdigit():
            return int(data[0]), int(data[1]), data[2], float(data[3]), float(data[4])
    raise SystemBuildError(f"{folder / 'forcefield.itp'} has no [ defaults ] line")


def check_defaults(folder: Path, polymer_defaults: tuple) -> None:
    sd = solvent_defaults(folder)
    ok = (sd[0] == polymer_defaults[0] and sd[1] == polymer_defaults[1]
          and str(sd[2]).lower() == str(polymer_defaults[2]).lower()
          and abs(sd[3] - polymer_defaults[3]) < 1e-3 and abs(sd[4] - polymer_defaults[4]) < 1e-3)
    if not ok:
        raise SystemBuildError(f"The water/ion force field {folder.name} has [ defaults ] {sd}, which does not match "
                               f"the polymer force field's {tuple(polymer_defaults)}; pass a matching solvent_ff.")


def read_molecules(path: Path) -> dict[str, Molecule]:
    """[ moleculetype ] blocks of an .itp (ions.itp, tip3p.itp) by molecule name."""
    return read_molecules_text(_active_lines(path))


def read_molecules_text(lines) -> dict[str, Molecule]:
    """[ moleculetype ] blocks of .itp lines (or text) by molecule name."""
    if isinstance(lines, str):
        lines = lines.splitlines()
    blocks: list[list[str]] = []
    for raw in lines:
        if re.match(r"\s*\[\s*moleculetype\s*\]", raw):
            blocks.append([])
        if blocks:
            blocks[-1].append(raw)
    out = {}
    for lines in blocks:
        section, name = None, None
        types, names, charges, masses, residue = [], [], [], [], ""
        for raw in lines:
            m = re.match(r"\s*\[\s*(\S+)\s*\]", raw)
            if m:
                section = m.group(1)
                continue
            data = raw.split(";")[0].split()
            if not data:
                continue
            if section == "moleculetype" and name is None:
                name = data[0]
            elif section == "atoms":
                types.append(data[1])
                residue = data[3]
                names.append(data[4])
                charges.append(float(data[6]))
                masses.append(float(data[7]) if len(data) > 7 else 0.0)
        if name:
            out[name] = Molecule(name, "\n".join(lines).strip() + "\n", tuple(types), residue, tuple(names),
                                 tuple(charges), tuple(masses))
    return out


def atomtype_lines(folder: Path, names: set[str]) -> dict[str, str]:
    """[ atomtypes ] lines for the given types."""
    found = {}
    for raw in _active_lines(folder / "ffnonbonded.itp"):
        data = raw.split(";")[0].split()
        if data and data[0] in names:
            found[data[0]] = raw.rstrip()
    missing = names - set(found)
    if missing:
        raise SystemBuildError(f"{folder / 'ffnonbonded.itp'} has no atom type(s) {', '.join(sorted(missing))}")
    return found


def type_mass(line: str) -> float:
    """Mass column of an [ atomtypes ] line (7- or 8-column layout)"""
    tok = line.split(";")[0].split()
    p = next(i for i, t in enumerate(tok) if t in {"A", "D", "S", "V"} and i >= 4)
    return float(tok[p - 2])


# rq-12060285
def scaled(molecule: Molecule, s: float) -> Molecule:
    """The molecule with its charges scaled."""
    lines, section = [], None
    for raw in molecule.text.splitlines():
        m = re.match(r"\s*\[\s*(\S+)\s*\]", raw)
        if m:
            section = m.group(1)
        body, sep, comment = raw.partition(";")
        data = body.split()
        if section == "atoms" and not m and len(data) >= 7:
            data[6] = f"{float(data[6]) * s:.5f}"
            raw = "  ".join(data) + (f"  ;{comment}" if sep else "")
        lines.append(raw)
    return Molecule(molecule.name, "\n".join(lines) + "\n", molecule.atom_types, molecule.residue,
                    molecule.atom_names, tuple(q * s for q in molecule.charges), molecule.masses)


# rq-96c58164
def madrid2019(names, with_bond_types: bool) -> tuple[dict[str, Molecule], dict[str, str], str]:
    """TIP4P/2005 water and Madrid-2019 ions"""
    import math
    from importlib import resources

    import yaml

    data = yaml.safe_load(resources.files("zpoly.system").joinpath("data/madrid2019.yaml").read_text())
    types, w = data["types"], data["water"]
    a = w["d_om"] / (2 * w["d_oh"] * math.cos(math.radians(w["angle_hoh"] / 2)))  # M = O + a (H1-O) + a (H2-O)
    d_hh = 2 * w["d_oh"] * math.sin(math.radians(w["angle_hoh"] / 2))
    texts = {ion: f"[ moleculetype ]\n{ion}  1\n\n[ atoms ]\n  1  {t}  1  {ion}  {ion}  1  {types[t][2]:.4f}  "
                  f"{types[t][1]:.4f}\n" for ion, t in data["ions"].items()}
    texts["SOL"] = (
        "[ moleculetype ]\nSOL  2\n\n[ atoms ]\n"
        f"  1  OW_tip4p  1  SOL  OW   1   0.0000  {types['OW_tip4p'][1]:.4f}\n"
        f"  2  HW_tip4p  1  SOL  HW1  1  {w['q_h']:.4f}  {types['HW_tip4p'][1]:.4f}\n"
        f"  3  HW_tip4p  1  SOL  HW2  1  {w['q_h']:.4f}  {types['HW_tip4p'][1]:.4f}\n"
        f"  4  MW_tip4p  1  SOL  MW   1  {w['q_m']:.4f}  0.0000\n\n"
        f"[ settles ]\n  1  1  {w['d_oh']:.5f}  {d_hh:.5f}\n\n"
        f"[ virtual_sites3 ]\n  4  1  2  3  1  {a:.9f}  {a:.9f}\n\n"
        "[ exclusions ]\n  1  2  3  4\n  2  1  3  4\n  3  1  2  4\n  4  1  2  3\n")
    molecules = {}
    for name in names:
        if name not in texts:
            raise SystemBuildError(f"Madrid-2019 has no molecule {name} (it has {', '.join(texts)})")
        molecules[name] = read_molecules_text(texts[name])[name]
    used = sorted({t for m in molecules.values() for t in m.atom_types})
    lines = {}
    for t in used:
        z, mass, _q, sigma, eps = types[t]
        bt = f"{t:<10s} " if with_bond_types else ""
        lines[t] = (f"  {t:<10s} {bt}{z:6d}{mass:12.4f}{0.0:10.4f}     {'D' if z == 0 else 'A'}"
                    f"{sigma:15.6e}{eps:15.6e}   ; Madrid-2019 / TIP4P/2005")
    lj = [t for t in used if types[t][4] > 0]
    pairs = []
    for i, ti in enumerate(lj):
        for tj in lj[i:]:
            given = data["pairs"].get(f"{ti}-{tj}") or data["pairs"].get(f"{tj}-{ti}")
            sigma, eps = given or ((types[ti][3] + types[tj][3]) / 2, math.sqrt(types[ti][4] * types[tj][4]))
            note = "Table III/IV" if given or ti == tj else "Lorentz-Berthelot"
            pairs.append(f"  {ti:<10s} {tj:<10s} 1 {sigma:12.6f} {eps:12.6f}   ; {note}")
    return molecules, lines, "\n".join(pairs)
