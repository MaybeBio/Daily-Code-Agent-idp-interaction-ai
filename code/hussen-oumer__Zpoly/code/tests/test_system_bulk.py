"""Tests for the bulk density/Tg system builder (rqm/system/bulk.md)."""

import dataclasses
import itertools
import shutil
import sys

import numpy as np
import pytest
import yaml
from test_forcefield_core import FakeBackend

from zpoly.forcefield import register_backend
from zpoly.system import SystemAnalysisError, SystemBuildError, analyze_tg, bulk_system, fit_tg

AMBER_LIKE = {
    "forcefield.itp": "[ defaults ]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n1 2 yes 0.5 0.8333\n",
    "ffnonbonded.itp": ("[ atomtypes ]\n#ifdef HEAVY_H\nK  19  99.99  0.0000  A  4.73602e-01  1.37235e-03\n#else\n"
                        "K           19      39.10    0.0000  A   4.73602e-01  1.37235e-03\n#endif\n"
                        "Na          11      22.99    0.0000  A   3.32840e-01  1.15897e-02\n"
                        "Cl          17      35.45    0.0000  A   4.40104e-01  4.18400e-01\n"),
    "ions.itp": "".join(f"[ moleculetype ]\n; molname nrexcl\n{n}  1\n\n[ atoms ]\n"
                        f"1  {t}  1  {n}  {n}  1  {q}\n\n" for n, t, q in
                        (("K", "K", "1.00000"), ("NA", "Na", "1.00000"), ("CL", "Cl", "-1.00000"))),
}


def write_solvent_ff(path):
    path.mkdir(parents=True, exist_ok=True)
    for name, text in AMBER_LIKE.items():
        (path / name).write_text(text)
    return path


@pytest.fixture(scope="module")
def solvent(tmp_path_factory):
    return write_solvent_ff(tmp_path_factory.mktemp("solv") / "amber99sb-ildn.ff")


@pytest.fixture(scope="module")
def bulk(tmp_path_factory, solvent):
    register_backend(FakeBackend("fake"))
    out = tmp_path_factory.mktemp("bulk") / "b"
    bulk_system("DMAPS", n=3, chains=4, box_nm=3.0, initiator="KPS", ff="fake", packer="builtin",
                solvent_ff=solvent, out=out)
    return out


def read_gro(path):
    lines = path.read_text().splitlines()
    n = int(lines[1])
    atoms = [(int(line[0:5]), line[5:10].strip(), line[10:15].strip(),
              np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])]) * 10) for line in lines[2:2 + n]]
    box = np.array([float(x) for x in lines[2 + n].split()]) * 10
    return atoms, box


def molecules(top):
    rows, on = [], False
    for line in top.read_text().splitlines():
        if line.strip().startswith("[ molecules ]"):
            on = True
            continue
        if on and line.strip() and not line.strip().startswith(";"):
            name, count = line.split()
            rows.append((name, int(count)))
    return rows


def chains_from(bulk_dir):
    atoms, box = read_gro(bulk_dir / "system.gro")
    n_chain = sum(1 for line in (bulk_dir / "polymer.itp").read_text().split("[ bonds ]")[0].splitlines()
                  if line.strip() and line.split()[0].isdigit())
    return atoms, box, n_chain


# rq-7f93393d
def test_bulk_folder_contents(bulk):
    names = {p.name for p in bulk.iterdir()}
    assert {"system.gro", "system.top", "polymer.itp", "polymer.ff", "01_em.mdp", "02_nvt_hot.mdp",
            "03_npt_compress.mdp", "04_npt_anneal.mdp", "run.sh", "run_slurm.sh", "protocol.yaml"} <= names
    assert sorted(n for n in names if n.startswith("cool_")) == sorted(f"cool_{t}K.mdp" for t in range(200, 601, 20))


# rq-a09233c5
def test_counterions_neutralize(bulk):
    assert molecules(bulk / "system.top") == [("polymer", 4), ("K", 8)]
    itp_q = sum(float(line.split()[6]) for line in (bulk / "polymer.itp").read_text().split("[ bonds ]")[0].splitlines()
                if line.strip() and line.split()[0].isdigit())
    assert abs(4 * itp_q + 8 * 1.0) < 1e-3
    atoms, _, n_chain = chains_from(bulk)
    assert len(atoms) == 4 * n_chain + 8 and all(a[1] == "K" for a in atoms[4 * n_chain:])


# rq-f413df4d
def test_chains_differ(bulk):
    atoms, _, n_chain = chains_from(bulk)
    mats = []
    for k in range(4):
        xyz = np.array([a[3] for a in atoms[k * n_chain:(k + 1) * n_chain]])
        mats.append(np.linalg.norm(xyz[:, None] - xyz[None], axis=2))
    for a, b in itertools.combinations(mats, 2):
        assert not np.allclose(a, b, atol=0.01)


# rq-1dba016a
def test_builtin_min_distance(bulk):
    atoms, box, n_chain = chains_from(bulk)
    mol_of = [k // n_chain if k < 4 * n_chain else 4 + (k - 4 * n_chain) for k in range(len(atoms))]
    heavy = [i for i, a in enumerate(atoms) if not a[2].startswith("H")]
    xyz = np.array([atoms[i][3] for i in heavy])
    owner = np.array([mol_of[i] for i in heavy])
    d = xyz[:, None] - xyz[None]
    d -= box * np.round(d / box)
    dist = np.linalg.norm(d, axis=2)
    other = owner[:, None] != owner[None]
    assert dist[other].min() >= 2.0 - 0.01  # .gro has 0.001 nm precision


# rq-b96476fb
def test_box_is_users(bulk, solvent, tmp_path):
    _, box, _ = chains_from(bulk)
    assert np.allclose(box, 30.0, atol=0.01)
    out = tmp_path / "rect"
    bulk_system("DMAPS", n=3, chains=2, box_nm=(3.0, 3.5, 4.0), ff="fake", packer="builtin", solvent_ff=solvent,
                out=out)
    _, rect = read_gro(out / "system.gro")
    assert np.allclose(rect, (30.0, 35.0, 40.0), atol=0.01)
    with pytest.raises(TypeError):
        bulk_system("DMAPS", n=3, ff="fake", out=tmp_path / "y")
    with pytest.raises(ValueError):
        bulk_system("DMAPS", n=3, chains=2, box_nm=-1.0, ff="fake", packer="builtin", solvent_ff=solvent,
                    out=tmp_path / "z")


# rq-db5c6f04
@pytest.mark.skipif(shutil.which("packmol") is None, reason="packmol not installed")
def test_packmol_packer(bulk, solvent, tmp_path):
    out = tmp_path / "pm"
    bulk_system("DMAPS", n=3, chains=4, box_nm=3.0, initiator="KPS", ff="fake", packer="packmol",
                solvent_ff=solvent, out=out)
    a, box = read_gro(out / "system.gro")
    b, _ = read_gro(bulk / "system.gro")
    assert [(x[1], x[2]) for x in a] == [(x[1], x[2]) for x in b]
    xyz = np.array([x[3] for x in a])
    assert (xyz.min(axis=0) > -5).all() and (xyz.max(axis=0) < box + 5).all()


# rq-af2c836d
def test_packmol_missing(solvent, tmp_path, monkeypatch):
    register_backend(FakeBackend("fake"))
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(SystemBuildError) as exc:
        bulk_system("DMAPS", n=3, chains=2, box_nm=3.0, initiator="KPS", ff="fake", packer="packmol",
                    solvent_ff=solvent, out=tmp_path / "x")
    assert "conda install -c conda-forge packmol" in str(exc.value) and 'packer="builtin"' in str(exc.value)


# rq-72fe414e
def test_unsupported_counterion(solvent, tmp_path):
    register_backend(FakeBackend("fake"))
    with pytest.raises(SystemBuildError, match=r"NH4\+"):
        bulk_system("DMAPS", n=3, chains=2, box_nm=3.0, initiator="APS", ff="fake", packer="builtin",
                    solvent_ff=solvent, out=tmp_path / "x")


def mdp(path):
    out = {}
    for line in path.read_text().splitlines():
        if "=" in line and not line.startswith(";"):
            k, v = line.split("=", 1)
            out[k.strip()] = v.split(";")[0].strip()
    return out


# rq-9ae6f9ab
def test_protocol_parameters(solvent, tmp_path):
    register_backend(FakeBackend("fake"))
    out = tmp_path / "p"
    bulk_system("DMAPS", n=3, chains=2, box_nm=3.0, ff="fake", packer="builtin", solvent_ff=solvent, out=out,
                t_high=500, t_low=300, t_step=50, times_ns={"cool": 0.5})
    cool = sorted(p.name for p in out.glob("cool_*.mdp"))
    assert cool == sorted(f"cool_{t}K.mdp" for t in (500, 450, 400, 350, 300))
    m = mdp(out / "cool_450K.mdp")
    assert m["ref-t"] == "450" and m["nsteps"] == "250000" and m["dt"] == "0.002"
    assert m["pcoupl"] == "C-rescale" and m["ref-p"] == "1"
    proto = yaml.safe_load((out / "protocol.yaml").read_text())
    assert [s["temperature"] for s in proto["cooling"]["steps"]] == [500, 450, 400, 350, 300]


# rq-9c761dab
def test_run_script_order(bulk):
    lines = [line for line in (bulk / "run.sh").read_text().splitlines() if "grompp" in line]
    names = [line.split("-f ")[1].split(".mdp")[0] for line in lines]
    expected = ["01_em", "02_nvt_hot", "03_npt_compress", "04_npt_anneal"]
    expected += [f"cool_{t}K" for t in range(600, 199, -20)]
    assert names == expected
    for prev, line in zip(["system"] + expected[:-1], lines, strict=True):
        assert f"-c {prev}.gro" in line


# rq-950969b7
def test_bilinear_fit():
    t = np.arange(600, 199, -20, dtype=float)
    d = np.where(t < 400, 1.30 - 0.0002 * (t - 400), 1.30 - 0.0006 * (t - 400))
    assert abs(fit_tg(t, d)["tg"] - 400) < 1.0


FAKE_GMX = r'''
import re, sys
args = sys.argv[1:]
edr, out = args[args.index("-f") + 1], args[args.index("-o") + 1]
t = float(re.search(r"cool_([0-9.]+)K", edr).group(1))
d = 1300 - (0.2 if t < 400 else 0.6) * (t - 400)
with open(out, "w") as fh:
    fh.write("# fake gmx energy\n@ title \"Density\"\n")
    for k in range(5):
        fh.write(f"{k * 100.0:8.1f} {d + (k - 2) * 0.1:10.4f}\n")
'''


# rq-1e5b84fd
def test_analyze_tg(bulk, tmp_path):
    gmx = tmp_path / "gmx"
    gmx.write_text(f"#!{sys.executable}\n{FAKE_GMX}")
    gmx.chmod(0o755)
    for step in yaml.safe_load((bulk / "protocol.yaml").read_text())["cooling"]["steps"]:
        (bulk / f"{step['name']}.edr").write_text("fake")
    result = analyze_tg(bulk, gmx=str(gmx))
    rows = (bulk / "tg_density.csv").read_text().splitlines()
    assert len(rows) == 1 + 21
    saved = yaml.safe_load((bulk / "tg_result.yaml").read_text())
    assert abs(result.tg - 400) < 1.0 and saved["tg_K"] == round(result.tg, 2)
    assert result.density_300 == pytest.approx(1.32, abs=1e-3)
    assert saved["density_300K_g_cm3"] == pytest.approx(1.32, abs=1e-3)
    for step in yaml.safe_load((bulk / "protocol.yaml").read_text())["cooling"]["steps"]:
        (bulk / f"{step['name']}.edr").unlink()


# rq-a0aaaf7b
def test_missing_energy_file(bulk, tmp_path):
    for step in yaml.safe_load((bulk / "protocol.yaml").read_text())["cooling"]["steps"]:
        if step["name"] != "cool_400K":
            (bulk / f"{step['name']}.edr").write_text("fake")
    with pytest.raises(SystemAnalysisError, match="cool_400K"):
        analyze_tg(bulk, gmx="true")
    for p in bulk.glob("cool_*.edr"):
        p.unlink()


class OplsLikeBackend(FakeBackend):
    def parametrize(self, mol, net_charge, work_dir):
        frag = super().parametrize(mol, net_charge, work_dir)
        return dataclasses.replace(frag, parameters=dataclasses.replace(frag.parameters,
                                                                         defaults=(1, 3, "yes", 0.5, 0.5)))


# rq-cca7d81d
def test_solvent_ff_mismatch(solvent, tmp_path):
    register_backend(OplsLikeBackend("oplslike"))
    with pytest.raises(SystemBuildError, match="defaults"):
        bulk_system("DMAPS", n=3, chains=2, box_nm=3.0, initiator="KPS", ff="oplslike", packer="builtin",
                    solvent_ff=solvent, out=tmp_path / "x")


# rq-f75e9890
def test_cli_bulk_and_analyze(solvent, run_inp, tmp_path, monkeypatch, capsys):
    register_backend(FakeBackend("fake"))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("GMXDATA", str(solvent.parent.parent))
    (solvent.parent.parent / "top").mkdir(exist_ok=True)
    if not (solvent.parent.parent / "top" / "amber99sb-ildn.ff").exists():
        shutil.copytree(solvent, solvent.parent.parent / "top" / "amber99sb-ildn.ff")
    assert run_inp("JOB = bulk\nMONOMER = DMAPS\nN = 3\nCHAINS = 2\nBOX = 3.0\nINITIATOR = KPS\nFF = fake\n"
                   "PACKER = builtin\nOUTPUT = b") == 0
    out = capsys.readouterr().out
    assert "b:" in out and "atoms" in out
    assert run_inp("JOB = analyze\nOUTPUT = b", name="a.inp") == 1
    assert "cool_600K" in capsys.readouterr().err
