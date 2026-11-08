"""Tests for the LLPS slab system builder (rqm/system/slab.md)."""

import shutil
import sys

import numpy as np
import pytest
import yaml
from test_forcefield_core import FakeBackend
from test_system_bulk import AMBER_LIKE, mdp, molecules, read_gro

from zpoly.forcefield import register_backend
from zpoly.system import SystemAnalysisError, analyze_slab, phase_concentrations, slab_system

TIP3P = """[ moleculetype ]
; molname nrexcl
SOL  2

[ atoms ]
1  OW  1  SOL  OW   1  -0.834  16.00000
2  HW  1  SOL  HW1  1   0.417   1.00800
3  HW  1  SOL  HW2  1   0.417   1.00800

#ifndef FLEXIBLE
[ settles ]
1  1  0.09572  0.15139
#else
[ bonds ]
1  2  1  0.09572  502416.0
1  3  1  0.09572  502416.0
#endif

[ exclusions ]
1  2  3
2  1  3
3  1  2
"""
WATER_TYPES = ("OW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01\n"
               "HW           1      1.008    0.0000  A   0.00000e+00  0.00000e+00\n")


BOX = (3.4, 3.4, 10.2)


@pytest.fixture(scope="module")
def solvent(tmp_path_factory):
    path = tmp_path_factory.mktemp("solv") / "amber99sb-ildn.ff"
    path.mkdir(parents=True)
    for name, text in AMBER_LIKE.items():
        (path / name).write_text(text + (WATER_TYPES if name == "ffnonbonded.itp" else ""))
    (path / "tip3p.itp").write_text(TIP3P)
    return path


@pytest.fixture(scope="module")
def slab(tmp_path_factory, solvent):
    register_backend(FakeBackend("fake"))
    out = tmp_path_factory.mktemp("slab") / "s"
    slab_system("DMAPS", n=3, chains=6, box_nm=BOX, slab_nm=3.4, initiator="KPS", ff="fake", packer="builtin",
                solvent_ff=solvent, out=out)
    return out


def make(solvent, tmp_path, **kw):
    register_backend(FakeBackend("fake"))
    args = dict(n=3, chains=6, box_nm=BOX, slab_nm=3.4, initiator="KPS", ff="fake", packer="builtin",
                solvent_ff=solvent, out=tmp_path / "x")
    args.update(kw)
    return slab_system("DMAPS", **args)


# rq-70bb800e
def test_box_is_users(slab):
    _, box = read_gro(slab / "packed.gro")
    assert np.allclose(box, (34.0, 34.0, 102.0), atol=0.01)


# rq-185b551d
def test_box_checks(solvent, tmp_path):
    for kw in ({"box_nm": (3.4, 3.4, 8.0)}, {"box_nm": (3.4, 4.0, 12.0)}, {"slab_nm": 11}):
        with pytest.raises(ValueError):
            make(solvent, tmp_path, **kw)


# rq-9b09c717
def test_centred_slab(slab):
    atoms, box = read_gro(slab / "packed.gro")
    proto = yaml.safe_load((slab / "protocol.yaml").read_text())
    t = proto["slab_thickness_nm"] * 10
    assert t == pytest.approx(34.0)
    z = np.array([a[3][2] for a in atoms])
    assert (z >= box[2] / 2 - t / 2 - 0.01).all() and (z <= box[2] / 2 + t / 2 + 0.01).all()


# rq-a6cd2763
def test_missing_user_inputs(solvent, tmp_path):
    for missing in ("chains", "box_nm", "slab_nm"):
        args = dict(n=3, chains=6, box_nm=BOX, slab_nm=3.4, ff="fake", packer="builtin", solvent_ff=solvent,
                    out=tmp_path / "m")
        del args[missing]
        with pytest.raises(TypeError):
            slab_system("DMAPS", **args)


# rq-fee4e59e
def test_builtin_prelude(slab):
    run = (slab / "run.sh").read_text()
    assert run.index("$GMX solvate") < run.index("-f 01_em.mdp -c system.gro")
    assert "genion -s genion.tpr -o system.gro -p system.top -pname NA -nname CL -conc 0.15 -neutral" in run
    assert molecules(slab / "system.top") == [("polymer", 6), ("K", 12)]


# rq-9dfe556a
@pytest.mark.skipif(shutil.which("packmol") is None, reason="packmol not installed")
def test_packmol_packs_everything(solvent, tmp_path):
    out = make(solvent, tmp_path, packer="packmol", chains=2)
    mols = molecules(out / "system.top")
    names = [m for m, _ in mols]
    assert names[:2] == ["polymer", "K"] and "SOL" in names and "NA" in names and "CL" in names
    proto = yaml.safe_load((out / "protocol.yaml").read_text())
    atoms, _ = read_gro(out / "system.gro")
    assert len(atoms) == proto["atoms"]
    sol = dict(mols)["SOL"]
    assert dict(mols)["NA"] == round(0.15 * sol / 55.5)
    run = (out / "run.sh").read_text()
    assert "solvate" not in run and "genion" not in run


# rq-1223ec9e
def test_no_salt(solvent, tmp_path):
    run = (make(solvent, tmp_path, salt=None) / "run.sh").read_text()
    assert "solvate" in run and "genion" not in run


# rq-cb81b1ce
def test_kcl(solvent, tmp_path):
    run = (make(solvent, tmp_path, salt="KCl") / "run.sh").read_text()
    assert "-pname K -nname CL" in run


# rq-7bc134f2
def test_production_temperatures(solvent, tmp_path):
    out = make(solvent, tmp_path, temperatures=(290, 310))
    for t in (290, 310):
        m = mdp(out / f"prod_{t}K.mdp")
        assert m["gen-vel"] == "yes" and m["ref-t"] == str(t) and m["pcoupl"] == "no"
    run = (out / "run.sh").read_text()
    assert "-f prod_290K.mdp -c 03_npt_z.gro" in run and "-f prod_310K.mdp -c 03_npt_z.gro" in run
    assert (out / "run_production.sh").is_file()


# rq-133ab185
def test_z_only_pressure(slab):
    m = mdp(slab / "03_npt_z.mdp")
    assert m["pcoupltype"] == "semiisotropic" and m["compressibility"] == "0 4.5e-5"


# rq-03fc94cd
def test_index_group(slab):
    text = (slab / "index.ndx").read_text().splitlines()
    assert text[0] == "[ Polymer ]"
    ids = [int(x) for line in text[1:] for x in line.split()]
    atoms, _ = read_gro(slab / "packed.gro")
    n_poly = len(atoms) - 12
    assert ids == list(range(1, n_poly + 1))


def tanh_profile(z):
    return 20 + 980 * 0.5 * (np.tanh((z - 12) / 0.5) - np.tanh((z - 18) / 0.5))


# rq-d067b409
def test_phase_concentrations():
    z = np.arange(0, 30.0, 0.1)
    dense, dilute, width = phase_concentrations(z, tanh_profile(z))
    assert dense == pytest.approx(1000, rel=0.02) and dilute == pytest.approx(20, rel=0.02)
    assert abs(width - 6.0) < 0.3


FAKE_GMX = r'''
import sys, numpy as np
args = sys.argv[1:]
out = args[args.index("-o") + 1]
z = np.arange(0, 30.0, 0.1)
rho = 20 + 980 * 0.5 * (np.tanh((z - 12) / 0.5) - np.tanh((z - 18) / 0.5))
with open(out, "w") as fh:
    fh.write("# fake gmx density\n@ title \"Polymer\"\n")
    for a, b in zip(z, rho):
        fh.write(f"{a:10.4f} {b:12.4f}\n")
'''


# rq-4b2d638c
def test_analyze_slab(slab, tmp_path):
    gmx = tmp_path / "gmx"
    gmx.write_text(f"#!{sys.executable}\n{FAKE_GMX}")
    gmx.chmod(0o755)
    prods = yaml.safe_load((slab / "protocol.yaml").read_text())["productions"]
    for p in prods:
        (slab / f"{p['name']}.xtc").write_text("fake")
    results = analyze_slab(slab, gmx=str(gmx))
    rows = (slab / "phase_diagram.csv").read_text().splitlines()
    assert len(rows) == 1 + len(prods) == 1 + len(results)
    for p in prods:
        assert (slab / f"density_profile_{p['name'][5:]}.csv").is_file()
        (slab / f"{p['name']}.xtc").unlink()
    assert results[0].dense == pytest.approx(1000, rel=0.02)


# rq-f0a89383
def test_missing_trajectory(slab):
    (slab / "prod_280K.xtc").write_text("fake")  # the first temperature is present, the second (300 K) is not
    with pytest.raises(SystemAnalysisError, match="prod_300K"):
        analyze_slab(slab, gmx="true")
    (slab / "prod_280K.xtc").unlink()


# rq-42bf805c
def test_cli_slab(solvent, run_inp, tmp_path):
    register_backend(FakeBackend("fake"))
    assert run_inp(f"""JOB = slab
                      MONOMER = DMAPS
                      N = 3
                      CHAINS = 4
                      BOX = 3.4 3.4 10.2
                      SLAB = 3.4
                      INITIATOR = KPS
                      FF = fake
                      PACKER = builtin
                      TEMPERATURES = 300 320
                      SOLVENT_FF = {solvent}
                      OUTPUT = s""") == 0
    assert (tmp_path / "s/prod_300K.mdp").is_file() and (tmp_path / "s/prod_320K.mdp").is_file()


# rq-96c58164
def test_tip4p2005_madrid(tmp_path):
    register_backend(FakeBackend("fake"))
    args = dict(n=3, chains=2, box_nm=BOX, slab_nm=3.4, initiator="KPS", ff="fake", packer="builtin",
                salt="NaCl", solvent="tip4p2005")
    out = slab_system("DMAPS", charge_scale=0.85, out=tmp_path / "m", **args)
    nb = (out / "polymer.ff/ffnonbonded.itp").read_text()
    for t in ("OW_tip4p", "HW_tip4p", "MW_tip4p", "NA_mad", "K_mad", "CL_mad"):
        assert f"  {t} " in nb
    pairs = {tuple(line.split()[:2]): line.split()[3:5] for line in nb.split("[ nonbond_params ]")[1].splitlines()
             if line.strip() and not line.startswith(";")}
    assert pairs[("K_mad", "OW_tip4p")] == ["0.289040", "1.400430"]
    assert pairs[("CL_mad", "NA_mad")] == ["0.300512", "1.438894"]
    assert pairs[("K_mad", "NA_mad")] == [f"{(0.230140 + 0.221737) / 2:.6f}", f"{(1.985740 * 1.472356) ** 0.5:.6f}"]
    ions = (out / "polymer.ff/ions.itp").read_text()
    assert "K_mad  1  K  K  1  0.8500" in ions and "CL_mad  1  CL  CL  1  -0.8500" in ions
    water = (out / "polymer.ff/water.itp").read_text()
    assert "[ virtual_sites3 ]" in water and "0.131937768" in water and "-1.1128" in water
    assert "-cs tip4p.gro" in (out / "run.sh").read_text()
    assert yaml.safe_load((out / "protocol.yaml").read_text())["solvent"] == "tip4p2005"
    with pytest.raises(ValueError, match="CHARGE_SCALE = 0.85"):
        slab_system("DMAPS", charge_scale=1.0, out=tmp_path / "x", **args)
