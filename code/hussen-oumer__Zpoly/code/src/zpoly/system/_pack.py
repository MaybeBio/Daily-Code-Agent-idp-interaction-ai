"""Packing molecules into a box"""

# rq-f3321ea0

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._errors import SystemBuildError

MIN_DISTANCE = 2.0  # Å
MAX_TRIES = 1000
PACKMOL_HINT = "Install Packmol (conda install -c conda-forge packmol) or use packer=\"builtin\"."


@dataclass
class Piece:
    """One molecule to place"""

    xyz: np.ndarray
    elements: list[str]
    region: tuple[float, float] | None = None
    pivots: tuple = ()  # (a, b, moving atom indices) per rotatable bond, for the EMC packer
    near: np.ndarray | None = None  # heavy x heavy: True if at most 3 bonds apart (EMC self-overlap check)

    @property
    def heavy(self) -> np.ndarray:
        # names starting with H are hydrogens
        return np.array([not e.upper().startswith("H") for e in self.elements])


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a - b*b + c*c - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a - b*b - c*c + d*d]])


def pack_builtin(pieces: Sequence[Piece], box: Sequence[float], seed: int = 0) -> list[np.ndarray]:
    """Random placement with periodic minimum-distance checks between heavy atoms"""
    rng = np.random.default_rng(seed)
    box = np.asarray(box, dtype=float)
    placed: list[np.ndarray] = []
    placed_heavy = np.zeros((0, 3))
    for k, piece in enumerate(pieces):
        centred = piece.xyz - piece.xyz.mean(axis=0)
        heavy = piece.heavy
        for _ in range(MAX_TRIES):
            rotated = centred @ _random_rotation(rng).T
            pos = rng.uniform(0, box)
            if piece.region is not None:
                # keep every atom inside the z region
                lo = piece.region[0] - rotated[:, 2].min()
                hi = piece.region[1] - rotated[:, 2].max()
                if lo > hi:
                    continue  # this orientation is taller than the region
                pos[2] = rng.uniform(lo, hi)
            xyz = rotated + pos
            h = xyz[heavy] if heavy.any() else xyz
            if len(placed_heavy):
                d = h[:, None, :] - placed_heavy[None, :, :]
                d -= box * np.round(d / box)
                if (np.einsum("ijk,ijk->ij", d, d) < MIN_DISTANCE ** 2).any():
                    continue
            placed.append(xyz)
            placed_heavy = np.vstack([placed_heavy, h])
            break
        else:
            raise SystemBuildError(f"The built-in packer could not place molecule {k + 1} of {len(pieces)} after "
                                   f"{MAX_TRIES} tries; lower the start density or use Packmol.")
    return placed


def _write_pdb(path: Path, piece: Piece) -> None:
    lines = [f"HETATM{i + 1:5d} {e:<4s} MOL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {e:>2s}"
             for i, (e, (x, y, z)) in enumerate(zip(piece.elements, piece.xyz, strict=True))]
    path.write_text("\n".join(lines) + "\nEND\n", encoding="ascii")


def pack_packmol(pieces: Sequence[Piece], box: Sequence[float], seed: int = 0, timeout: float = 3600.0
                 ) -> list[np.ndarray]:
    """Pack with Packmol"""
    exe = shutil.which("packmol")
    if exe is None:
        raise SystemBuildError(f"packmol was not found on PATH. {PACKMOL_HINT}")
    lx, ly, lz = (float(b) for b in box)
    margin = 1.0  # Å from the box faces, so periodic images do not overlap
    with tempfile.TemporaryDirectory(prefix="packmol_") as tmp:
        tmpd = Path(tmp)
        groups: list[tuple[Piece, int]] = []
        for piece in pieces:
            last = groups[-1][0] if groups else None
            if last is not None and last.xyz.shape == piece.xyz.shape and np.array_equal(last.xyz, piece.xyz) \
                    and last.region == piece.region:
                groups[-1] = (last, groups[-1][1] + 1)
            else:
                groups.append((piece, 1))
        lines = [f"tolerance {MIN_DISTANCE}", "filetype pdb", "output packed.pdb", f"seed {seed + 1}", ""]
        for k, (piece, count) in enumerate(groups):
            _write_pdb(tmpd / f"m{k}.pdb", piece)
            z0, z1 = piece.region if piece.region else (0.0, lz)
            lines += [f"structure m{k}.pdb", f"  number {count}",
                      f"  inside box {margin} {margin} {z0 + margin} {lx - margin} {ly - margin} {z1 - margin}",
                      "end structure", ""]
        (tmpd / "packmol.inp").write_text("\n".join(lines), encoding="ascii")
        try:
            with open(tmpd / "packmol.inp") as fh:
                proc = subprocess.run([exe], stdin=fh, cwd=tmpd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise SystemBuildError(f"Packmol did not finish within {timeout:g} s.") from None
        out = tmpd / "packed.pdb"
        if proc.returncode != 0 or not out.is_file() or "ERROR" in proc.stdout:
            tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-15:])
            raise SystemBuildError(f"Packmol failed (exit {proc.returncode}):\n{tail}")
        coords = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                           for line in out.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))])
    result, start = [], 0
    for piece, count in groups:
        n = len(piece.elements)
        for _ in range(count):
            result.append(coords[start:start + n])
            start += n
    if start != len(coords):
        raise SystemBuildError("Packmol returned a different number of atoms than requested.")
    return result


def _emc(pieces, box, seed=0):
    from ._emc import emc_pack

    return emc_pack(pieces, box, seed)


# name -> packer(pieces, box, seed)
PACKERS = {"packmol": pack_packmol, "builtin": pack_builtin, "emc": _emc}


def pack(pieces: Sequence[Piece], box: Sequence[float], packer: str = "packmol", seed: int = 0) -> list[np.ndarray]:
    if packer not in PACKERS:
        raise ValueError(f"packer must be one of {', '.join(PACKERS)}, got {packer!r}")
    return PACKERS[packer](pieces, box, seed)
