"""3D construction of one polymer chain by unit-by-unit growth"""

# rq-ee0a0744

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from ..zwitterions import Chain, Ends, Initiator, Library, Zwitterion, build_chain, default_library
from ._errors import BuildError

log = logging.getLogger(__name__)

N_CONFORMERS = 20
N_ANGLES = 24
ENOUGH_FREE = 48  # coil: stop scanning conformers once this many clash-free placements exist
SOFT_CONTACT = 3.2  # clash score counts (3.2 - d)^2 below this distance (Å)
CLASH_FREE = 2.6  # a placement is clash-free when every scored distance is >= this (Å)
BACKTRACK_BELOW = 1.8  # backtrack when the chosen placement is closer than this (Å)
BACKTRACK_UNITS = 3
MAX_BACKTRACKS = 30  # per growth attempt
STRICT_ATTEMPTS = 3  # growth attempts with the strict backtracking limit
LOOSE_ATTEMPTS = 2  # then attempts that accept closer contacts and need a UFF cleanup
LOOSE_BACKTRACK_BELOW = 1.0
FALLBACK_RELAX_STEPS = 500
FINAL_MIN_DISTANCE = 1.5  # final structure must have heavy atoms > 3 bonds apart at least this far
CH_BOND = 1.09
CONFORMATIONS = ("coil", "extended")


# Templates
@dataclass
class _Template:
    """An embedded fragment whose attachment points are temporary carbon caps"""

    mol: Chem.Mol
    positions: list[np.ndarray]  # one (n_atoms, 3) array per conformer
    cap_in: int | None
    cap_out: int | None
    anchor_in: int | None
    anchor_out: int | None
    keep: list[int]  # copied atoms: each heavy atom followed by its hydrogens
    heavy_keep: list[int]
    junction_out: list[int]  # heavy atoms within 2 bonds of anchor_out
    junction_depth: np.ndarray  # their bond distance from anchor_out
    scored: list[int]  # heavy_keep plus cap_out (where the next unit's head atom will sit)
    scored_depth: np.ndarray  # bond distance of each scored point from anchor_in
    reach: float  # max distance from anchor_in to any scored point, over conformers


# rq-57e067c6
def _make_template(smiles: str, seed: int) -> _Template:
    """smiles carries [*:1] (incoming) and/or [*:2] (outgoing) attachment points"""
    rw = Chem.RWMol(Chem.MolFromSmiles(smiles))
    caps: dict[int, int] = {}
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            caps[atom.GetAtomMapNum()] = atom.GetIdx()
            atom.SetAtomicNum(6)
            atom.SetAtomMapNum(0)
            atom.SetNoImplicit(False)
    Chem.SanitizeMol(rw)
    mol = Chem.AddHs(rw)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    with rdBase.BlockLogs():
        cids = list(AllChem.EmbedMultipleConfs(mol, N_CONFORMERS, params))
        if not cids:
            params.useRandomCoords = True
            cids = list(AllChem.EmbedMultipleConfs(mol, N_CONFORMERS, params))
        if not cids:
            raise BuildError(f"RDKit could not embed the fragment {smiles}.")
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
        elif AllChem.UFFHasAllMoleculeParams(mol):
            AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=500)
    positions = [mol.GetConformer(c).GetPositions() for c in cids]

    cap_in, cap_out = caps.get(1), caps.get(2)
    cap_atoms = {c for c in (cap_in, cap_out) if c is not None}
    anchor = {k: mol.GetAtomWithIdx(c).GetNeighbors()[0].GetIdx() for k, c in caps.items()}
    positions = _with_second_torsion(mol, positions, cap_in, anchor.get(1))
    keep: list[int] = []
    heavy_keep: list[int] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 or atom.GetIdx() in cap_atoms:
            continue
        keep.append(atom.GetIdx())
        heavy_keep.append(atom.GetIdx())
        keep.extend(nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1)
    topo = Chem.GetDistanceMatrix(mol)
    # heavy atoms within 2 bonds of the outgoing anchor
    junction_out: list[int] = []
    junction_depth: list[int] = []
    if 2 in anchor:
        for j in heavy_keep:
            depth = int(topo[anchor[2], j])
            if depth <= 2:
                junction_out.append(j)
                junction_depth.append(depth)
    scored = heavy_keep + ([cap_out] if cap_out is not None else [])
    # bond distance from the head atom
    scored_depth = np.array([int(topo[anchor[1], s]) for s in scored]) if 1 in anchor else np.full(len(scored), 99)
    origin = anchor.get(1, heavy_keep[0])
    reach = max(float(np.linalg.norm(P[scored] - P[origin], axis=1).max()) for P in positions)
    return _Template(mol, positions, cap_in, cap_out, anchor.get(1), anchor.get(2), keep, heavy_keep,
                     junction_out, np.array(junction_depth, dtype=int), scored, scored_depth, reach)


SECOND_TORSION_ANGLES = (2 * np.pi / 3, 4 * np.pi / 3)


def _with_second_torsion(mol: Chem.Mol, positions: list[np.ndarray], cap_in: int | None,
                         head: int | None) -> list[np.ndarray]:
    """Add conformers rotated about the second torsion."""
    if head is None or cap_in is None:
        return positions
    others = [nb.GetIdx() for nb in mol.GetAtomWithIdx(head).GetNeighbors()
              if nb.GetAtomicNum() > 1 and nb.GetIdx() != cap_in]
    if len(others) != 1 or mol.GetBondBetweenAtoms(head, others[0]).IsInRing():
        return positions
    pivot = others[0]
    side, frontier = {pivot}, [pivot]
    while frontier:
        nxt = []
        for idx in frontier:
            for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
                j = nb.GetIdx()
                if j != head and j not in side:
                    side.add(j)
                    nxt.append(j)
        frontier = nxt
    moving = np.array(sorted(side))
    out = list(positions)
    for P in positions:
        axis = P[pivot] - P[head]
        for angle in SECOND_TORSION_ANGLES:
            Q = P.copy()
            Q[moving] = (P[moving] - P[head]) @ _rotation(axis, angle).T + P[head]
            out.append(Q)
    return out


# Geometry helpers
def _rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * (k @ k)


def _align(points: np.ndarray, src_a, src_b, dst_a, dst_b) -> np.ndarray:
    """Move points rigidly so src_a -> dst_a and direction src_a->src_b -> dst_a->dst_b"""
    u = (src_b - src_a) / np.linalg.norm(src_b - src_a)
    v = (dst_b - dst_a) / np.linalg.norm(dst_b - dst_a)
    cross = np.cross(u, v)
    s, c = float(np.linalg.norm(cross)), float(u @ v)
    if s < 1e-9:
        if c > 0:
            rot = np.eye(3)
        else:
            perp = np.cross(u, [1.0, 0.0, 0.0] if abs(u[0]) < 0.9 else [0.0, 1.0, 0.0])
            rot = _rotation(perp, np.pi)
    else:
        rot = _rotation(cross, np.arctan2(s, c))
    return (points - src_a) @ rot.T + dst_a


# Growth
@dataclass
class _Placed:
    positions: np.ndarray  # template coordinates after placement
    heavy: np.ndarray  # positions of heavy_keep atoms


def _place_unit(i, templates, placed, rng, conformation, first_atom):
    """Return (positions, min_scored_distance) for unit i given units 0..i-1."""
    tpl = templates[i]
    if i == 0:
        cid = int(rng.integers(len(tpl.positions))) if conformation == "coil" else 0
        return tpl.positions[cid].copy(), float("inf")

    prev_tpl, prev = templates[i - 1], placed[i - 1]
    cap_pos = prev.positions[prev_tpl.cap_out]
    tail_pos = prev.positions[prev_tpl.anchor_out]
    # pairs within 3 bonds are not scored
    junction = set(prev_tpl.junction_out)
    prev_far = [prev.positions[j] for j in prev_tpl.heavy_keep if j not in junction]
    pool_parts = [p.heavy for p in placed[: i - 1]] + ([np.array(prev_far)] if prev_far else [])
    pool = np.vstack(pool_parts) if pool_parts else np.zeros((0, 3))
    if len(pool):
        pool = pool[np.linalg.norm(pool - cap_pos, axis=1) <= tpl.reach + SOFT_CONTACT]
    near = prev.positions[prev_tpl.junction_out] if prev_tpl.junction_out else np.zeros((0, 3))
    bonded = (tpl.scored_depth[:, None] + 1 + prev_tpl.junction_depth[None, :]) <= 3

    axis = tail_pos - cap_pos
    angles = 2 * np.pi * np.arange(N_ANGLES) / N_ANGLES + rng.uniform(-0.1, 0.1)
    rots = np.stack([_rotation(axis, a) for a in angles])  # (A, 3, 3)

    candidates = []  # (score, min_d, positions)
    order = rng.permutation(len(tpl.positions)) if conformation == "coil" else range(len(tpl.positions))
    n_free = 0
    for c in order:
        if conformation == "coil" and n_free >= ENOUGH_FREE:
            break  # enough clash-free options to choose from; skip the remaining conformers
        conf_positions = tpl.positions[c]
        aligned = _align(conf_positions, conf_positions[tpl.anchor_in], conf_positions[tpl.cap_in],
                         cap_pos, tail_pos)
        moved = np.einsum("aij,nj->ani", rots, aligned - cap_pos) + cap_pos  # (A, n_atoms, 3)
        pts = moved[:, tpl.scored, :]
        min_d = np.full(len(angles), np.inf)
        score = np.zeros(len(angles))
        for targets, mask in ((pool, None), (near, bonded)):
            if not len(targets):
                continue
            d = np.linalg.norm(pts[:, :, None, :] - targets[None, None, :, :], axis=3)
            if mask is not None:
                d = np.where(mask[None, :, :], np.inf, d)
            min_d = np.minimum(min_d, d.reshape(len(angles), -1).min(axis=1))
            score += (np.clip(SOFT_CONTACT - d, 0, None) ** 2).reshape(len(angles), -1).sum(axis=1)
        for a in range(len(angles)):
            candidates.append((float(score[a]), float(min_d[a]), moved[a]))
        n_free += int((min_d >= CLASH_FREE).sum())

    free = [c for c in candidates if c[1] >= CLASH_FREE]
    if free and conformation == "coil":
        # pick from the roomier half
        cutoff = float(np.median([c[0] for c in free]))
        roomy = [c for c in free if c[0] <= cutoff]
        choice = roomy[int(rng.integers(len(roomy)))]
    elif free:
        def reach_of(cand):
            pos = cand[2]
            tip = pos[tpl.cap_out] if tpl.cap_out is not None else pos[tpl.heavy_keep].mean(axis=0)
            return float(np.linalg.norm(tip - first_atom))
        choice = max(free, key=reach_of)
    else:
        choice = min(candidates, key=lambda c: c[0])
    return choice[2], choice[1]


# rq-7414a764
def _grow(templates: list[_Template], conformation: str, rng: np.random.Generator,
          backtrack_below: float = BACKTRACK_BELOW) -> list[_Placed]:
    placed: list[_Placed] = []
    backtracks = 0
    failures: dict[int, int] = {}  # unit index -> consecutive failures there
    first_atom = None
    i = 0
    while i < len(templates):
        positions, min_d = _place_unit(i, templates, placed, rng, conformation, first_atom)
        if min_d < backtrack_below and i > 1:
            if backtracks >= MAX_BACKTRACKS:
                raise BuildError(
                    f"Could not place unit {i + 1} without clashes after {MAX_BACKTRACKS} retries "
                    f"(closest contact {min_d:.2f} Å).", min_d)
            backtracks += 1
            # back up twice as far on repeated failure
            depth = BACKTRACK_UNITS * 2 ** failures.get(i, 0)
            failures[i] = failures.get(i, 0) + 1
            back_to = max(1, i - depth)
            del placed[back_to:]
            i = back_to
            continue
        tpl = templates[i]
        placed.append(_Placed(positions, positions[tpl.heavy_keep]))
        if first_atom is None:
            first_atom = positions[tpl.heavy_keep[0]]
        i += 1
    return placed


# Residue names
_ALNUM = re.compile(r"[^A-Za-z0-9]")
_SUFFIXES = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _base_code(name: str, width: int) -> str:
    code = _ALNUM.sub("", name).upper()[:width]
    return code.ljust(width, "X")


# rq-65ea2a82
def _assign_prefixes(items: Sequence[tuple[str, str]]) -> dict[str, str]:
    """items"""
    used: dict[str, str] = {}  # prefix -> monomer
    out: dict[str, str] = {}
    for key, code in items:
        if key in out:
            continue
        if code in used and used[code] != key:
            for ch in _SUFFIXES:
                alt = code[0] + ch
                if alt not in used:
                    code = alt
                    break
            else:
                raise BuildError(f"Ran out of residue prefixes for '{key}'.")
        used[code] = key
        out[key] = code
    return out


def _position_letter(i: int, n: int) -> str:
    """F, R or E for a unit's chain position."""
    if i == 0:
        return "F"
    return "E" if i == n - 1 else "R"


# Public result type
# rq-12389ffb
@dataclass(frozen=True, eq=False)
class PolymerStructure:
    """A 3D polymer chain"""

    mol: Chem.Mol
    chain: Chain
    sequence: tuple[str, ...]
    units: tuple[Zwitterion, ...] = ()
    initiator: Initiator | None = None
    residue_names: Mapping[str, str] = field(default_factory=dict)
    templates: tuple[str, ...] = ()
    conformation: str = "coil"
    seed: int = 0
    min_nonbonded_distance: float = float("inf")
    radius_of_gyration: float = 0.0
    end_to_end_distance: float = 0.0
    relaxed: bool = False  # True when a UFF cleanup ran (requested, or needed after loose growth)


def _min_nonbonded_distance(mol: Chem.Mol, cutoff: float = 3.0) -> float:
    """Closest heavy-atom contact beyond 3 bonds."""
    heavy = np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    xyz = mol.GetConformer().GetPositions()[heavy]
    neighbours = {int(i): {nb.GetIdx() for nb in mol.GetAtomWithIdx(int(i)).GetNeighbors()} for i in heavy}

    def within3(a: int, b: int) -> bool:
        frontier, seen = {a}, {a}
        for _ in range(3):
            frontier = {n for f in frontier for n in neighbours.get(f, ()) if n not in seen}
            if b in frontier:
                return True
            seen |= frontier
        return False

    best = float("inf")
    chunk = 1024
    for start in range(0, len(heavy), chunk):
        block = xyz[start:start + chunk]
        d = np.linalg.norm(block[:, None, :] - xyz[None, :, :], axis=2)
        rows, cols = np.nonzero(d < cutoff)
        for r, c in zip(rows, cols, strict=True):
            a, b = int(heavy[start + r]), int(heavy[c])
            if a < b and d[r, c] < best and not within3(a, b):
                best = float(d[r, c])
    return best


# rq-b9773ded
def build_polymer(
    sequence: Sequence[str | Zwitterion],
    *,
    initiator: Initiator | str | None = None,
    ends: Ends | str | None = None,
    conformation: str = "coil",
    seed: int = 0,
    relax_iterations: int = 0,
    library: Library | None = None,
    cap: str = "CH3",
) -> PolymerStructure:
    """Build a 3D chain for the given sequence, with optional initiator end groups"""
    if conformation not in CONFORMATIONS:
        raise ValueError(f"conformation must be one of {CONFORMATIONS}, got {conformation!r}")
    if relax_iterations < 0:
        raise ValueError("relax_iterations must be >= 0")
    lib = library if library is not None else default_library()
    units = [lib.get(u) if isinstance(u, str) else u for u in sequence]
    chain = build_chain(units, initiator=initiator, ends=ends, library=lib, cap=cap)  # validates everything
    init = lib.initiator(initiator) if isinstance(initiator, str) else initiator

    # residues: prefix + F/R/E; end groups: prefix + H/T
    prefixes = _assign_prefixes([(z.name, z.residue or _base_code(z.name, 2)) for z in units])
    names: dict[str, str] = dict(prefixes)
    fragments: list[tuple[str, str]] = []  # (residue name, template smiles) in chain order
    if chain.head_group is not None:
        names[f"{init.name}:head"] = _base_code(init.name, 2) + "H"
        fragments.append((names[f"{init.name}:head"], chain.head_group.replace("[*:1]", "[*:2]")))
    for i, z in enumerate(units):
        fragments.append((prefixes[z.name] + _position_letter(i, len(units)), z.repeat_unit_smiles()))
    if chain.tail_group is not None:
        names[f"{init.name}:tail"] = _base_code(init.name, 2) + "T"
        fragments.append((names[f"{init.name}:tail"], chain.tail_group))

    cache: dict[str, _Template] = {}
    template_ids: dict[str, int] = {}  # template SMILES -> index in PolymerStructure.templates
    templates = []
    for _, smi in fragments:
        if smi not in cache:
            cache[smi] = _make_template(smi, seed)
            template_ids[smi] = len(template_ids)
        templates.append(cache[smi])

    placed, needs_cleanup = _grow_with_restarts(templates, conformation, seed)
    mol = _assemble(templates, placed, [res for res, _ in fragments],
                    [template_ids[smi] for _, smi in fragments], cap)

    if Chem.MolToSmiles(Chem.RemoveHs(mol)) != chain.smiles:
        raise BuildError("Internal error: assembled structure does not match the chain connectivity.")
    steps = max(relax_iterations, FALLBACK_RELAX_STEPS) if needs_cleanup else relax_iterations
    relaxed = False
    if steps:
        if AllChem.UFFHasAllMoleculeParams(mol):
            AllChem.UFFGetMoleculeForceField(mol).Minimize(maxIts=steps)
            relaxed = True
        elif needs_cleanup:
            raise BuildError("Growth needed a force-field cleanup, but UFF has no parameters for some atoms.")
        else:
            log.warning("UFF has no parameters for some atoms; skipping relaxation.")

    min_d = _min_nonbonded_distance(mol)
    if min_d < FINAL_MIN_DISTANCE:
        raise BuildError(f"Built structure has overlapping atoms ({min_d:.2f} Å apart).", min_d)

    xyz = mol.GetConformer().GetPositions()
    heavy = np.array([a.GetAtomicNum() > 1 for a in mol.GetAtoms()])
    hx = xyz[heavy]
    rg = float(np.sqrt(((hx - hx.mean(axis=0)) ** 2).sum(axis=1).mean()))
    first = int(mol.GetIntProp("_first_backbone"))
    last = int(mol.GetIntProp("_last_backbone"))
    return PolymerStructure(
        mol=mol, chain=chain, sequence=tuple(z.name for z in units), units=tuple(units),
        initiator=init if (chain.head_group or chain.tail_group) else None,
        templates=tuple(template_ids),
        residue_names=MappingProxyType(dict(names)),
        conformation=conformation, seed=seed, min_nonbonded_distance=min_d, radius_of_gyration=rg,
        end_to_end_distance=float(np.linalg.norm(xyz[first] - xyz[last])), relaxed=relaxed,
    )


def _grow_with_restarts(templates: list[_Template], conformation: str, seed: int) -> tuple[list[_Placed], bool]:
    """Strict growth with seeded restarts"""
    last: BuildError | None = None
    phases = ((STRICT_ATTEMPTS, BACKTRACK_BELOW, False), (LOOSE_ATTEMPTS, LOOSE_BACKTRACK_BELOW, True))
    attempt = 0
    for n_attempts, limit, loose in phases:
        for _ in range(n_attempts):
            try:
                return _grow(templates, conformation, np.random.default_rng([seed, attempt]), limit), loose
            except BuildError as exc:
                last = exc
                log.info("Growth attempt %d failed (%s).", attempt + 1, exc)
            attempt += 1
    raise BuildError(f"{last} Gave up after {attempt} growth attempts; try another seed or "
                     "conformation='extended'.", last.min_distance if last else None)


class _ResidueWriter:
    """Adds one residue's atoms with names."""

    def __init__(self, out: Chem.RWMol, coords: list[np.ndarray], res_name: str, resnum: int,
                 template: int = -1) -> None:
        self.out, self.coords, self.res_name, self.resnum = out, coords, res_name, resnum
        self.template = template
        self.counters: dict[str, int] = {}

    def add(self, atom: Chem.Atom, xyz: np.ndarray, template_atom: int = -1) -> int:
        symbol = atom.GetSymbol()
        self.counters[symbol] = self.counters.get(symbol, 0) + 1
        name = f"{symbol}{self.counters[symbol]}"
        if len(name) > 4:
            raise BuildError(f"Residue {self.resnum} has too many {symbol} atoms for 4-character names.")
        new = Chem.Atom(atom)
        new.SetNoImplicit(True)
        new.SetNumExplicitHs(0)
        new.SetProp("res", self.res_name)
        new.SetIntProp("resnum", self.resnum)
        new.SetProp("name", name)
        new.SetIntProp("template", self.template)  # rq-8c171443: DFT charges map through these
        new.SetIntProp("template_atom", template_atom)
        self.coords.append(xyz)
        return self.out.AddAtom(new)


def _assemble(templates: list[_Template], placed: list[_Placed], res_names: list[str],
              template_ids: list[int] | None = None, cap: str = "H") -> Chem.Mol:
    """Join placed templates into one molecule."""
    out = Chem.RWMol()
    coords: list[np.ndarray] = []
    prev_anchor_out: int | None = None
    first_backbone = last_backbone = None
    n = len(templates)
    for r, (tpl, pl) in enumerate(zip(templates, placed, strict=True)):
        add = _ResidueWriter(out, coords, res_names[r], r + 1, template_ids[r] if template_ids else -1).add
        index: dict[int, int] = {}

        for t_idx in tpl.keep:
            index[t_idx] = add(tpl.mol.GetAtomWithIdx(t_idx), pl.positions[t_idx], t_idx)
            # free ends: H after the anchor's hydrogens (tail CH3: below)
            for c, anchor, is_end in ((tpl.cap_in, tpl.anchor_in, r == 0),
                                      (tpl.cap_out, tpl.anchor_out, r == n - 1 and cap == "H")):
                if is_end and c is not None and anchor is not None and _last_h_of(tpl, anchor) == t_idx:
                    vec = pl.positions[c] - pl.positions[anchor]
                    h_pos = pl.positions[anchor] + CH_BOND * vec / np.linalg.norm(vec)
                    h = add(Chem.Atom(1), h_pos)
                    out.AddBond(index[anchor], h, Chem.BondType.SINGLE)
        for bond in tpl.mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in index and b in index:
                out.AddBond(index[a], index[b], bond.GetBondType())
        if r == n - 1 and cap == "CH3" and tpl.cap_out is not None and tpl.anchor_out is not None:
            # rq-8994396e: template methyl cap -> chain-end CH3
            c = add(tpl.mol.GetAtomWithIdx(tpl.cap_out), pl.positions[tpl.cap_out], tpl.cap_out)
            out.AddBond(index[tpl.anchor_out], c, Chem.BondType.SINGLE)
            for nb in tpl.mol.GetAtomWithIdx(tpl.cap_out).GetNeighbors():
                if nb.GetAtomicNum() == 1:
                    h = add(nb, pl.positions[nb.GetIdx()], nb.GetIdx())
                    out.AddBond(c, h, Chem.BondType.SINGLE)
        if prev_anchor_out is not None and tpl.anchor_in is not None:
            out.AddBond(prev_anchor_out, index[tpl.anchor_in], Chem.BondType.SINGLE)
        prev_anchor_out = index[tpl.anchor_out] if tpl.anchor_out is not None else None
        # backbone end atoms
        if tpl.anchor_in is not None and tpl.anchor_out is not None:
            if first_backbone is None:
                first_backbone = index[tpl.anchor_in]
            last_backbone = index[tpl.anchor_out]
    mol = out.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(coords):
        conf.SetAtomPosition(i, xyz.tolist())
    mol.AddConformer(conf, assignId=True)
    mol.SetIntProp("_first_backbone", first_backbone if first_backbone is not None else 0)
    mol.SetIntProp("_last_backbone", last_backbone if last_backbone is not None else mol.GetNumAtoms() - 1)
    return mol


def _last_h_of(tpl: _Template, anchor: int) -> int:
    """Template index after which an end hydrogen for this anchor is inserted"""
    hs = [nb.GetIdx() for nb in tpl.mol.GetAtomWithIdx(anchor).GetNeighbors() if nb.GetAtomicNum() == 1]
    return hs[-1] if hs else anchor
