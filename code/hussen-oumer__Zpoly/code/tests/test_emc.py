"""EMC-style packer and density mode (rqm/system/emc_packer.md, bulk.md, slab.md)."""

import dataclasses

import numpy as np
import pytest
import yaml
from test_forcefield_core import FakeBackend
from test_system_bulk import read_gro

from zpoly.forcefield import register_backend
from zpoly.structure import build_polymer, homopolymer
from zpoly.system import SystemBuildError, _emc, bulk_system, slab_system
from zpoly.system._assemble import chain_piece
from zpoly.system._pack import Piece

BOX = (3.4, 3.4, 10.2)
MADRID = {"solvent": "tip4p2005", "charge_scale": 0.85}  # built-in water and ions


@pytest.fixture(scope="module")
def dmaps5():
    s = build_polymer(homopolymer("DMAPS", 5), initiator="KPS")
    return s, chain_piece(s, emc=True)


def min_intermolecular(coords, pieces, box):
    heavy = [xyz[p.heavy] for xyz, p in zip(coords, pieces, strict=True)]
    best = np.inf
    for i in range(len(heavy)):
        for j in range(i + 1, len(heavy)):
            d = heavy[i][:, None, :] - heavy[j][None, :, :]
            d -= np.asarray(box) * np.round(d / np.asarray(box))
            best = min(best, float(np.sqrt((d ** 2).sum(-1)).min()))
    return best


# rq-fb226d8e
def test_pivots_keep_geometry(dmaps5):
    s, piece = dmaps5
    rng = np.random.default_rng(1)
    xyz = piece.xyz.copy()
    for _ in range(50):
        new = _emc.pivot(piece, xyz, rng, _emc.EmcSettings(), max_angle=180)
        xyz = new if new is not None else xyz
    assert not np.allclose(xyz, piece.xyz)
    mol = s.mol
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        assert np.linalg.norm(xyz[i] - xyz[j]) == pytest.approx(np.linalg.norm(piece.xyz[i] - piece.xyz[j]), abs=1e-6)
    def cos_angle(c, p, q, r):
        u, v = c[q] - c[p], c[r] - c[p]
        return np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)

    for a in mol.GetAtoms():
        nb = [n.GetIdx() for n in a.GetNeighbors()]
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                t = (a.GetIdx(), nb[x], nb[y])
                assert cos_angle(xyz, *t) == pytest.approx(cos_angle(piece.xyz, *t), abs=1e-6)


# rq-53277b50
def test_slab_without_overlaps(tmp_path):
    register_backend(FakeBackend("fake"))
    out = slab_system("DMAPS", n=6, chains=4, box_nm=BOX, slab_nm=4.0, initiator="KPS", ff="fake", packer="emc",
                      **MADRID,
                      out=tmp_path / "s")
    atoms, box = read_gro(out / "packed.gro")
    z = np.array([a[3][2] for a in atoms])
    assert (z >= box[2] / 2 - 20.01).all() and (z <= box[2] / 2 + 20.01).all()
    assert "$GMX solvate -cp packed.gro" in (out / "run.sh").read_text()


def test_emc_no_overlaps(dmaps5):
    _, piece = dmaps5
    box = (40.0, 40.0, 40.0)
    pieces = [dataclasses.replace(piece) for _ in range(6)] + [Piece(np.zeros((1, 3)), ["K"]) for _ in range(6)]
    coords = _emc.emc_pack(pieces, box, seed=3)
    assert min_intermolecular(coords, pieces, box) >= 2.0 - 1e-6


# rq-758b974f
def test_chain_built_once(tmp_path, monkeypatch):
    register_backend(FakeBackend("fake"))
    from zpoly.system import _bulk

    calls = []
    real = _bulk.build_polymer
    monkeypatch.setattr(_bulk, "build_polymer", lambda *a, **k: calls.append(1) or real(*a, **k))
    bulk_system("DMAPS", n=3, chains=5, box_nm=4.0, initiator="KPS", ff="fake", packer="emc",
                      **MADRID, out=tmp_path / "b")
    assert len(calls) == 1


# rq-41f479b2
def test_reproducible(dmaps5):
    _, piece = dmaps5
    pieces = [dataclasses.replace(piece) for _ in range(3)]
    a = _emc.emc_pack(pieces, (35.0, 35.0, 35.0), seed=5)
    b = _emc.emc_pack(pieces, (35.0, 35.0, 35.0), seed=5)
    assert all(np.array_equal(x, y) for x, y in zip(a, b, strict=True))


# rq-df4bcf3b
def test_too_dense(dmaps5):
    _, piece = dmaps5
    pieces = [dataclasses.replace(piece) for _ in range(20)]
    with pytest.raises(SystemBuildError, match="overlapping atom pairs"):
        _emc.emc_pack(pieces, (15.0, 15.0, 15.0), seed=0, settings=_emc.EmcSettings(max_sweeps=5))


# rq-fdfb3cb2
def test_new_move(dmaps5, monkeypatch):
    _, piece = dmaps5
    used = []

    def jiggle(mol, rng, settings, world):
        used.append(1)
        return mol.xyz + rng.normal(scale=0.1, size=3)

    monkeypatch.setitem(_emc.MOVES, "jiggle", jiggle)
    pieces = [dataclasses.replace(piece) for _ in range(8)]
    try:
        _emc.emc_pack(pieces, (30.0, 30.0, 30.0), seed=0,
                      settings=_emc.EmcSettings(move_weights={"jiggle": 1.0}, max_sweeps=3, insertion_trials=2))
    except SystemBuildError:
        pass
    assert used
    with pytest.raises(ValueError, match="No Monte Carlo move"):
        _emc.emc_pack(pieces, (30.0, 30.0, 30.0), settings=_emc.EmcSettings(move_weights={"nope": 1.0}))


# rq-7c268316
def test_slab_density_mode(tmp_path):
    register_backend(FakeBackend("fake"))
    out = slab_system("DMAPS", n=3, chains=4, density=1.0, aspect=3, initiator="KPS", ff="fake", packer="emc",
                      **MADRID,
                      out=tmp_path / "d")
    proto = yaml.safe_load((out / "protocol.yaml").read_text())
    lx, ly, lz = proto["box_nm"]
    mass = proto["slab_density_g_cm3"] * (lx * 1e-7) ** 3 / 1.66053906660e-24  # u
    assert lx == ly and lz == pytest.approx(3 * lx, abs=2e-5) and proto["density_g_cm3"] == 1.0
    _, cube = read_gro(out / "packed.gro")
    cube_nm = cube[0] / 10
    assert cube_nm == pytest.approx((mass * 1.66053906660e-24 / 0.35) ** (1 / 3) * 1e7, rel=1e-3)
    mdp = (out / "00_compress.mdp").read_text()
    rate = float(next(line for line in mdp.splitlines() if line.startswith("deform")).split("=")[1].split()[0])
    assert rate * 1000 == pytest.approx(lx - cube_nm, rel=1e-4) and "gen-vel              = no" in mdp
    run = (out / "run.sh").read_text()
    order = ["00_em.mdp", "00_heat.mdp", "-t 00_heat.cpt", "trjconv", f"-box {lx:.5f} {ly:.5f} {lz:.5f} -center",
             "solvate -cp slab_box.gro", "genion", "-f 01_em.mdp -c system.gro"]
    positions = [run.index(x) for x in order]
    assert positions == sorted(positions)
    with pytest.raises(ValueError):
        slab_system("DMAPS", n=3, chains=4, density=1.0, box_nm=BOX, ff="fake", out=tmp_path / "x")
    with pytest.raises(ValueError):
        slab_system("DMAPS", n=3, chains=4, density=1.0, aspect=2.5, ff="fake", out=tmp_path / "y")


# rq-da3d8837
def test_bulk_density_mode(tmp_path):
    register_backend(FakeBackend("fake"))
    out = bulk_system("DMAPS", n=3, chains=4, density=1.0, initiator="KPS", ff="fake", packer="emc",
                      **MADRID,
                      out=tmp_path / "b")
    proto = yaml.safe_load((out / "protocol.yaml").read_text())
    edge = proto["box_nm"][0]
    assert proto["box_nm"] == [edge] * 3 and proto["start_density_g_cm3"] == pytest.approx(1.0, rel=1e-3)
    _, packed = read_gro(out / "system.gro")
    assert packed[0] / 10 == pytest.approx(proto["packing_box_nm"][0], rel=1e-4) and packed[0] / 10 > edge
    run = (out / "run.sh").read_text()
    assert run.index("01_heat.mdp") < run.index("-f 01_deform.mdp -c 01_heat.gro -t 01_heat.cpt") < run.index(
        "02_nvt_hot.mdp")
    with pytest.raises(ValueError):
        bulk_system("DMAPS", n=3, chains=4, density=1.0, box_nm=4.0, ff="fake", out=tmp_path / "x")
