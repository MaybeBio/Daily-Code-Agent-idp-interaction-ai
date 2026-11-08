"""Read Amber parameter files (gaff2.dat, frcmod) into GROMACS units."""

# rq-e2b2a86f

from __future__ import annotations

import os
from pathlib import Path

from ._errors import ForceFieldError
from ._model import AtomType, ParameterSet

KCAL = 4.184
AMBER_DEFAULTS = (1, 2, "yes", 0.5, 0.8333)
AMBER_RTP_BONDEDTYPES = (1, 1, 9, 4, 1, 3, 1, 0)
SIGMA_FROM_RSTAR = 2 ** (5 / 6) / 10  # sigma (nm) = R* (Å) * 2^(5/6) / 10

_SECTIONS = ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB", "HBON", "CMAP")


def _key(text: str, n: int) -> tuple[str, ...]:
    parts = tuple(p.strip() for p in text.split("-"))
    if len(parts) != n or not all(parts):
        raise ValueError(f"bad type key {text!r}")
    return parts


class _Reader:
    def __init__(self) -> None:
        self.masses: dict[str, float] = {}
        self.nonbonded: dict[str, tuple[float, float]] = {}
        self.bonds: dict = {}
        self.angles: dict = {}
        self.dihedrals: dict = {}
        self.impropers: dict = {}
        self._continuing: tuple[str, ...] | None = None

    def mass(self, line: str) -> None:
        tok = line.split()
        self.masses[tok[0]] = float(tok[1])

    def bond(self, line: str) -> None:
        k, r0 = (float(x) for x in line[5:].split()[:2])
        self.bonds[_key(line[:5], 2)] = (r0 / 10, 2 * k * KCAL * 100)

    def angle(self, line: str) -> None:
        k, theta = (float(x) for x in line[8:].split()[:2])
        self.angles[_key(line[:8], 3)] = (theta, 2 * k * KCAL)

    def dihedral(self, line: str) -> None:
        key = _key(line[:11], 4)
        idivf, pk, phase, pn = (float(x) for x in line[11:].split()[:4])
        term = (phase, pk / idivf * KCAL, int(abs(round(pn))))
        if self._continuing == key:
            self.dihedrals[key] = self.dihedrals[key] + (term,)
        else:
            self.dihedrals[key] = (term,)  # a new occurrence replaces an earlier definition
        self._continuing = key if pn < 0 else None

    def improper(self, line: str) -> None:
        pk, phase, pn = (float(x) for x in line[11:].split()[:3])
        self.impropers[_key(line[:11], 4)] = (phase, pk * KCAL, int(abs(round(pn))))

    def nonbond(self, line: str) -> None:
        tok = line.split()
        self.nonbonded[tok[0]] = (float(tok[1]), float(tok[2]))

    def result(self) -> ParameterSet:
        atomtypes = {
            t: AtomType(t, 0, self.masses[t], r * SIGMA_FROM_RSTAR, eps * KCAL)
            for t, (r, eps) in self.nonbonded.items() if t in self.masses
        }
        return ParameterSet(AMBER_DEFAULTS, AMBER_RTP_BONDEDTYPES, atomtypes, self.bonds, self.angles,
                            self.dihedrals, self.impropers)


def _blocks(lines: list[str], start: int) -> tuple[list[str], int]:
    """Lines from `start` up to (not including) the next blank line"""
    out = []
    i = start
    while i < len(lines) and lines[i].strip():
        out.append(lines[i])
        i += 1
    return out, i + 1


# rq-455dff93
def read_amber_parameters(path: str | os.PathLike, *, frcmod: bool = False) -> ParameterSet:
    """Parse gaff2.dat or an frcmod file (GROMACS units)."""
    try:
        lines = Path(path).read_text(encoding="ascii", errors="replace").splitlines()
    except OSError as exc:
        raise ForceFieldError(f"Cannot read Amber parameter file {path}: {exc}") from exc
    r = _Reader()
    try:
        if frcmod:
            _read_frcmod(lines, r)
        else:
            _read_dat(lines, r)
    except (ValueError, IndexError) as exc:
        raise ForceFieldError(f"Cannot parse Amber parameter file {path}: {exc}") from exc
    return r.result()


def _read_dat(lines: list[str], r: _Reader) -> None:
    i = 1  # title
    block, i = _blocks(lines, i)
    for line in block:
        r.mass(line)
    i += 1  # hydrophilic atom types (one line)
    for handler in (r.bond, r.angle, r.dihedral, r.improper):
        block, i = _blocks(lines, i)
        for line in block:
            handler(line)
        r._continuing = None
    _, i = _blocks(lines, i)  # 10-12 H-bond terms
    equivalences, i = _blocks(lines, i)
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i < len(lines) and lines[i].split()[:1] == ["MOD4"]:
        block, i = _blocks(lines, i + 1)
        for line in block:
            if line.strip() == "END":
                break
            r.nonbond(line)
    for line in equivalences:
        tok = line.split()
        if tok and tok[0] in r.nonbonded:
            for other in tok[1:]:
                r.nonbonded.setdefault(other, r.nonbonded[tok[0]])


def _read_frcmod(lines: list[str], r: _Reader) -> None:
    handlers = {"MASS": r.mass, "BOND": r.bond, "ANGL": r.angle, "DIHE": r.dihedral,
                "IMPR": r.improper, "NONB": r.nonbond}
    i = 1  # title
    while i < len(lines):
        head = lines[i].strip()[:4].upper()
        if head in _SECTIONS:
            block, i = _blocks(lines, i + 1)
            handler = handlers.get(head)
            if handler is not None:
                for line in block:
                    handler(line)
            r._continuing = None
        else:
            i += 1
