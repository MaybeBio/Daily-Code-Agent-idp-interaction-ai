"""EMC-style packer"""

# rq-fdfb3cb2

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np

from ._errors import SystemBuildError
from ._pack import Piece, _random_rotation


@dataclass(frozen=True)
class EmcSettings:
    min_distance: float = 2.0  # Å, heavy atoms of different molecules (or > 3 bonds apart) must be farther
    cell: float = 3.0  # Å, cell-list edge (>= min_distance)
    insertion_trials: int = 200  # random placements per molecule; the one with fewest overlaps is kept
    pivots_per_copy: int = 20  # pivot moves that give each chain copy its own conformation
    max_sweeps: int = 500  # relaxation sweeps per inflation step, over the overlapping molecules
    inflation_steps: int = 6  # overlap radius grown from start_radius to min_distance in this many steps
    start_radius: float = 0.8  # Å
    max_shift: float = 2.0  # Å, translate move
    max_angle: float = 30.0  # degrees, rotate and pivot moves
    move_weights: dict = field(default_factory=lambda: {"translate": 1.0, "rotate": 1.0, "pivot": 1.0,
                                                        "local_pivot": 4.0, "reinsert": 0.5})


# Cell list (periodic)
class CellList:
    """Periodic cell list of placed heavy atoms."""

    def __init__(self, box: np.ndarray, cell: float) -> None:
        self.box = box
        self.n = np.maximum(1, np.floor(box / cell).astype(int))
        self.xyz = np.zeros((0, 3))
        self.owner = np.zeros(0, dtype=int)
        self._sorted = None

    def set(self, xyz: np.ndarray, owner: np.ndarray) -> None:
        self.xyz, self.owner, self._sorted = xyz, owner, None

    def _cells(self, xyz: np.ndarray) -> np.ndarray:
        return np.floor((xyz % self.box) / self.box * self.n).astype(int) % self.n

    def _index(self):
        if self._sorted is None:
            c = self._cells(self.xyz)
            lin = (c[:, 0] * self.n[1] + c[:, 1]) * self.n[2] + c[:, 2]
            order = np.argsort(lin, kind="stable")
            starts = np.searchsorted(lin[order], np.arange(self.n.prod() + 1))
            self._sorted = (order, starts)
        return self._sorted

    def overlaps(self, xyz: np.ndarray, exclude: int, cutoff: float) -> tuple[int, np.ndarray]:
        """Overlapping pairs with other molecules, and the atoms involved."""
        none = (0, np.zeros(0, dtype=int))
        if not len(self.xyz) or not len(xyz):
            return none
        order, starts = self._index()
        c = self._cells(xyz)
        shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
        nb = (c[:, None, :] + shifts[None, :, :]) % self.n
        lin = ((nb[..., 0] * self.n[1] + nb[..., 1]) * self.n[2] + nb[..., 2])
        query = np.repeat(np.arange(len(xyz)), 27)
        lin = lin.ravel()
        if (self.n < 3).any():  # few cells per side: the 27 neighbours repeat, keep each cell once per atom
            pairs = np.unique(np.stack([query, lin], axis=1), axis=0)
            query, lin = pairs[:, 0], pairs[:, 1]
        counts = starts[lin + 1] - starts[lin]
        total = int(counts.sum())
        if total == 0:
            return none
        offsets = np.repeat(starts[lin] - np.cumsum(counts) + counts, counts) + np.arange(total)
        other = order[offsets]
        q = np.repeat(query, counts)
        keep = self.owner[other] != exclude
        d = xyz[q[keep]] - self.xyz[other[keep]]
        d -= self.box * np.round(d / self.box)
        close = np.einsum("ij,ij->i", d, d) < cutoff * cutoff
        return int(close.sum()), np.unique(q[keep][close])


# Molecules being packed
@dataclass
class _Molecule:
    piece: Piece
    xyz: np.ndarray  # all atoms, Å
    heavy: np.ndarray  # bool per atom
    clash: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))  # atoms (all-atom index) overlapping
    region: tuple[float, float] | None = None  # z region the molecule must stay in (slab)


def _inside(mol: _Molecule) -> bool:
    if mol.region is None:
        return True
    z = mol.xyz[:, 2]
    return bool(z.min() >= mol.region[0] and z.max() <= mol.region[1])


def _intra_ok(piece: Piece, xyz: np.ndarray, moved: np.ndarray, settings: EmcSettings) -> bool:
    """No self-overlap after a pivot."""
    if piece.near is None:
        return True
    heavy_idx = np.flatnonzero(piece.heavy)
    pos = {a: k for k, a in enumerate(heavy_idx)}
    m = np.array([pos[a] for a in moved if a in pos], dtype=int)
    if not len(m):
        return True
    rest = np.setdiff1d(np.arange(len(heavy_idx)), m)
    d = xyz[heavy_idx[m]][:, None, :] - xyz[heavy_idx[rest]][None, :, :]
    close = np.einsum("ijk,ijk->ij", d, d) < settings.min_distance ** 2
    return not (close & ~piece.near[np.ix_(m, rest)]).any()


def _rotate_about(xyz: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    v = xyz - origin
    c, s = math.cos(angle), math.sin(angle)
    return origin + v * c + np.cross(axis, v) * s + np.outer(v @ axis, axis) * (1 - c)


def _step(rng) -> float:
    """Step-size factor between 0.1 and 1 (log-uniform)"""
    return 10 ** rng.uniform(-1, 0)


def pivot(piece: Piece, xyz: np.ndarray, rng, settings: EmcSettings, max_angle: float | None = None,
          choices=None):
    """Rotate the atoms on one side of a rotatable bond about it"""
    pool = piece.pivots if choices is None else choices
    if not pool:
        return None
    a, b, moving = pool[rng.integers(len(pool))]
    limit = math.radians(max_angle if max_angle is not None else settings.max_angle * _step(rng))
    new = xyz.copy()
    new[moving] = _rotate_about(xyz[moving], xyz[a], xyz[b] - xyz[a], rng.uniform(-limit, limit))
    return new if _intra_ok(piece, new, moving, settings) else None


def _translate(mol: _Molecule, rng, settings: EmcSettings, world):
    step = settings.max_shift * _step(rng)
    return mol.xyz + rng.uniform(-step, step, 3)


def _rotate(mol: _Molecule, rng, settings: EmcSettings, world):
    centre = mol.xyz.mean(axis=0)
    axis = rng.normal(size=3)
    limit = math.radians(settings.max_angle) * _step(rng)
    return _rotate_about(mol.xyz, centre, axis, rng.uniform(-limit, limit))


def _pivot(mol: _Molecule, rng, settings: EmcSettings, world):
    return pivot(mol.piece, mol.xyz, rng, settings)


def _local_pivot(mol: _Molecule, rng, settings: EmcSettings, world):
    """Pivot near an overlapping atom."""
    if not len(mol.clash) or not mol.piece.pivots:
        return None
    atom = mol.clash[rng.integers(len(mol.clash))]
    near = sorted((len(m), k) for k, (_a, _b, m) in enumerate(mol.piece.pivots) if atom in m)[:3]
    return pivot(mol.piece, mol.xyz, rng, settings, choices=[mol.piece.pivots[k] for _, k in near])


def _reinsert(mol: _Molecule, rng, settings: EmcSettings, world):
    """Move the molecule to the best of new random spots."""
    best, best_n = None, None
    for _ in range(max(1, settings.insertion_trials // 4)):
        trial = _random_placement(mol.xyz, mol.region, world.box, rng)
        if trial is not None:
            n = world.overlap(trial)
            if best_n is None or n < best_n:
                best, best_n = trial, n
    return best


# name -> move(mol, rng, settings, world) -> coordinates | None
MOVES = {"translate": _translate, "rotate": _rotate, "pivot": _pivot, "local_pivot": _local_pivot,
         "reinsert": _reinsert}


def _copy_conformation(piece: Piece, rng, settings: EmcSettings) -> np.ndarray:
    xyz = piece.xyz.copy()
    for _ in range(settings.pivots_per_copy if piece.pivots else 0):
        new = pivot(piece, xyz, rng, settings, max_angle=180.0)
        if new is not None:
            xyz = new
    return xyz


def _random_placement(xyz: np.ndarray, region, box: np.ndarray, rng) -> np.ndarray | None:
    rotated = (xyz - xyz.mean(axis=0)) @ _random_rotation(rng).T
    pos = rng.uniform(0, box)
    if region is not None:
        lo, hi = region[0] - rotated[:, 2].min(), region[1] - rotated[:, 2].max()
        if lo > hi:
            return None
        pos[2] = rng.uniform(lo, hi)
    return rotated + pos


# rq-53277b50 rq-41f479b2 rq-df4bcf3b
def emc_pack(pieces: Sequence[Piece], box: Sequence[float], seed: int = 0,
             settings: EmcSettings | None = None) -> list[np.ndarray]:
    """Pack pieces (chains with pivots, ions) into a periodic box"""
    s = settings or EmcSettings()
    unknown = set(s.move_weights) - set(MOVES)
    if unknown:
        raise ValueError(f"No Monte Carlo move named {', '.join(sorted(unknown))} (moves: {', '.join(MOVES)})")
    rng = np.random.default_rng(seed)
    box = np.asarray(box, dtype=float)
    grid = CellList(box, max(s.cell, s.min_distance))
    mols: list[_Molecule] = []

    def rebuild():
        xyz = np.vstack([m.xyz[m.heavy] for m in mols]) if mols else np.zeros((0, 3))
        owner = np.concatenate([np.full(int(m.heavy.sum()), k) for k, m in enumerate(mols)]) if mols else []
        grid.set(xyz, np.asarray(owner, dtype=int))

    state = {"radius": s.start_radius}  # the overlap radius, inflated below

    def overlap(k: int, xyz: np.ndarray, record: bool = False) -> int:
        n, atoms = grid.overlaps(xyz[mols[k].heavy], k, state["radius"])
        if record:
            mols[k].clash = np.flatnonzero(mols[k].heavy)[atoms]
        return n

    # configurational-bias insertion
    for k, piece in enumerate(pieces):
        conf = _copy_conformation(piece, rng, s)
        best, best_n = None, None
        mols.append(_Molecule(piece, conf, piece.heavy, region=piece.region))
        for _ in range(s.insertion_trials):
            trial = _random_placement(conf, mols[k].region, box, rng)
            if trial is None:
                continue
            n = overlap(k, trial)
            if best_n is None or n < best_n:
                best, best_n = trial, n
                if n == 0:
                    break
        if best is None:
            raise SystemBuildError(f"Molecule {k + 1} does not fit into the slab region in any orientation; "
                                   "make the slab thicker.")
        mols[k].xyz = best
        rebuild()

    # MC relaxation with an inflating overlap radius
    names = list(s.move_weights)
    weights = np.array([s.move_weights[n] for n in names], dtype=float)
    weights /= weights.sum()

    def relax() -> bool:
        for _ in range(s.max_sweeps):
            counts = [overlap(k, m.xyz, record=True) for k, m in enumerate(mols)]
            bad = [k for k, c in enumerate(counts) if c]
            if not bad:
                return True
            for k in rng.permutation(bad):
                mol = mols[k]
                name = names[rng.choice(len(names), p=weights)]
                world = SimpleNamespace(box=box, overlap=lambda xyz, k=k: overlap(k, xyz))
                new = MOVES[name](mol, rng, s, world)
                if new is None:
                    continue
                old_xyz = mol.xyz
                mol.xyz = new
                if _inside(mol) and overlap(k, new) <= overlap(k, old_xyz):
                    rebuild()
                else:
                    mol.xyz = old_xyz
        return False

    for radius in np.linspace(min(s.start_radius, s.min_distance), s.min_distance, max(1, s.inflation_steps)):
        state["radius"] = radius
        done = relax()
    if done:
        return [m.xyz for m in mols]
    left = sum(overlap(k, m.xyz) for k, m in enumerate(mols)) // 2
    raise SystemBuildError(f"The EMC packer left {left} overlapping atom pairs after {s.max_sweeps} sweeps; "
                           "enlarge the box (or the slab), or use fewer chains.")
