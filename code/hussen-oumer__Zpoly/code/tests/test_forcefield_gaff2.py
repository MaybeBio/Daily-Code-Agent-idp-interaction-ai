"""Tests for the GAFF2 backend (rqm/forcefield/gaff2.md) with fake AmberTools programs."""

import shutil
import sys
import textwrap

import pytest
from rdkit import Chem

from zpoly.forcefield import (
    ForceFieldUnavailableError,
    Gaff2Backend,
    ParametrizationError,
    parametrize,
    read_amber_parameters,
    write_gromacs,
)
from zpoly.forcefield._gaff2 import mol2_types_and_charges
from zpoly.structure import build_polymer, homopolymer

FAKE_ANTECHAMBER = r'''
import os, sys, time
from rdkit import Chem
args = sys.argv[1:]
with open(os.environ["FAKE_LOG"], "a") as fh:
    fh.write("antechamber " + " ".join(args) + "\n")
mode = os.environ.get("FAKE_MODE", "")
if mode == "fail":
    print("sqm failed"); sys.exit(1)
if mode == "hang":
    time.sleep(60)
inp, out = args[args.index("-i") + 1], args[args.index("-o") + 1]
mol = Chem.MolFromMolFile(inp, removeHs=False, sanitize=False)
conf = mol.GetConformer()
lines = ["@<TRIPOS>MOLECULE", "frag", f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL", "bcc", "",
         "@<TRIPOS>ATOM"]
for a in mol.GetAtoms():
    p = conf.GetAtomPosition(a.GetIdx())
    sym = a.GetSymbol()
    lines.append(f"{a.GetIdx()+1:7d} {sym}{a.GetIdx()+1:<5d} {p.x:10.4f} {p.y:10.4f} {p.z:10.4f} "
                 f"{sym.lower()}g {1:4d} MOL {0.05 * a.GetFormalCharge():10.6f}")
lines.append("@<TRIPOS>BOND")
for k, b in enumerate(mol.GetBonds(), 1):
    lines.append(f"{k:6d} {b.GetBeginAtomIdx()+1:5d} {b.GetEndAtomIdx()+1:5d} 1")
open(out, "w").write("\n".join(lines) + "\n")
'''

FAKE_PARMCHK2 = r'''
import os, sys
from rdkit import Chem
args = sys.argv[1:]
with open(os.environ["FAKE_LOG"], "a") as fh:
    fh.write("parmchk2 " + " ".join(args) + "\n")
out = args[args.index("-o") + 1]
mol = Chem.MolFromMolFile("frag.sdf", removeHs=False, sanitize=False)
t = [a.GetSymbol().lower() + "g" for a in mol.GetAtoms()]
masses = {"cg": 12.01, "hg": 1.008, "ng": 14.01, "og": 16.0, "sg": 32.06, "pg": 30.97}
nb = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
bonds = {tuple(sorted((t[b.GetBeginAtomIdx()], t[b.GetEndAtomIdx()]))) for b in mol.GetBonds()}
angles, dihes = set(), set()
for j in range(len(t)):
    for i in nb[j]:
        for k in nb[j]:
            if i < k:
                angles.add(min((t[i], t[j], t[k]), (t[k], t[j], t[i])))
for b in mol.GetBonds():
    j, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    dihes.add(tuple(sorted((t[j], t[k]))))
L = ["fake frcmod", "MASS"] + [f"{x:<2s} {masses[x]:.3f}" for x in sorted(set(t))] + [""]
L += ["BOND"] + [f"{a:<2s}-{b:<2s}  300.0   1.500" for a, b in sorted(bonds)] + [""]
L += ["ANGLE"] + [f"{a:<2s}-{b:<2s}-{c:<2s}   50.0   109.50" for a, b, c in sorted(angles)] + [""]
L += ["DIHE"] + [f"X -{a:<2s}-{b:<2s}-X     1    0.150       0.0             3." for a, b in sorted(dihes)] + [""]
L += ["IMPROPER", "", "NONBON"] + [f"  {x:<2s}  1.9000  0.1000" for x in sorted(set(t))] + ["", ""]
open(out, "w").write("\n".join(L))
'''

GAFF2_DAT = """AMBER General Force Field (fake for tests)
cg 12.01 0.878
hg 1.008 0.135

cg hg
cg-hg  347.0   1.093

cg-cg-hg   46.4   110.07

X -cg-cg-X    1    0.156       0.0             3.


MOD4      RE
  cg  1.9080  0.1094
  hg  1.4870  0.0157

END
"""


@pytest.fixture
def fake_amber(tmp_path, monkeypatch):
    prefix = tmp_path / "amber"
    bindir = prefix / "bin"
    bindir.mkdir(parents=True)
    for name, body in (("antechamber", FAKE_ANTECHAMBER), ("parmchk2", FAKE_PARMCHK2)):
        path = bindir / name
        path.write_text(f"#!{sys.executable}\n{body}")
        path.chmod(0o755)
    parm = prefix / "dat" / "leap" / "parm"
    parm.mkdir(parents=True)
    (parm / "gaff2.dat").write_text(GAFF2_DAT)
    log = tmp_path / "calls.log"
    monkeypatch.setenv("PATH", f"{bindir}:{__import__('os').environ['PATH']}")
    monkeypatch.delenv("AMBERHOME", raising=False)
    monkeypatch.setenv("FAKE_LOG", str(log))
    return prefix, log


def methane():
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return mol


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(textwrap.dedent(text))
    return path


# rq-8d5a0d31
def test_ambertools_missing(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ForceFieldUnavailableError) as exc:
        Gaff2Backend().check_available()
    assert "conda install -c conda-forge ambertools" in exc.value.hint


# rq-5392e74f
def test_gaff2_dat_from_prefix(fake_amber):
    prefix, _ = fake_amber
    backend = Gaff2Backend()
    backend.check_available()
    assert backend._locate()[2] == prefix.resolve() / "dat" / "leap" / "parm" / "gaff2.dat"


# rq-6c330b7c
def test_gaff2_dat_missing(fake_amber):
    prefix, _ = fake_amber
    (prefix / "dat" / "leap" / "parm" / "gaff2.dat").unlink()
    with pytest.raises(ForceFieldUnavailableError) as exc:
        Gaff2Backend().check_available()
    assert "AMBERHOME" in exc.value.hint


# rq-bd8eb0c7
def test_tool_command_lines(fake_amber, tmp_path):
    _, log = fake_amber
    mol = Chem.AddHs(Chem.MolFromSmiles("CS(=O)(=O)[O-]"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=1)
    Gaff2Backend().parametrize(mol, -1, tmp_path / "work")
    calls = log.read_text()
    ac = next(line for line in calls.splitlines() if line.startswith("antechamber"))
    pc = next(line for line in calls.splitlines() if line.startswith("parmchk2"))
    for frag in ("-c bcc", "-nc -1", "-at gaff2", "-fi mdl"):
        assert frag in ac
    assert "-s gaff2" in pc


MOL2 = """@<TRIPOS>MOLECULE
m
5 4 1 0 0
SMALL
bcc

@<TRIPOS>ATOM
      1 C1      0.0 0.0 0.0 c3 1 MOL -0.1000
      2 {second}      1.0 0.0 0.0 hc 1 MOL  0.0250
      3 H2      0.0 1.0 0.0 hc 1 MOL  0.0250
      4 H3      0.0 0.0 1.0 hc 1 MOL  0.0250
      5 H4     -1.0 0.0 0.0 hc 1 MOL  0.0250
@<TRIPOS>BOND
"""


# rq-9b94cbfd
def test_mol2_types_and_charges(tmp_path):
    path = write(tmp_path, "m.mol2", MOL2.format(second="H1"))
    types, charges = mol2_types_and_charges(path, methane())
    assert types == ("c3", "hc", "hc", "hc", "hc")
    assert charges == (-0.1, 0.025, 0.025, 0.025, 0.025)


# rq-b018fb93
def test_mol2_atom_mismatch(tmp_path):
    path = write(tmp_path, "m.mol2", MOL2.format(second="O1"))
    with pytest.raises(ParametrizationError):
        mol2_types_and_charges(path, methane())


def dat_with(bonds="", angles="", dihedrals="", impropers="", masses="c3 12.01", nonbond="  c3  1.9080  0.1094"):
    """A gaff2.dat-shaped file: each section ends with exactly one blank line, like the real file."""
    parts = ["title", masses, "", "c3"]
    for section in (bonds, angles, dihedrals, impropers):
        if section:
            parts.append(section)
        parts.append("")
    parts += ["", "", "MOD4      RE", nonbond, "", "END"]  # empty H-bond and equivalence sections
    return "\n".join(parts) + "\n"


# rq-f5a523d0
def test_unit_conversion(tmp_path):
    path = tmp_path / "p.dat"
    path.write_text(dat_with(bonds="c3-c3  300.9   1.5375", angles="c3-c3-c3   62.9   111.51"))
    ps = read_amber_parameters(path)
    b0, kb = ps.bondtypes[("c3", "c3")]
    assert b0 == pytest.approx(0.15375) and kb == pytest.approx(251793.12)
    th, k = ps.angletypes[("c3", "c3", "c3")]
    assert th == pytest.approx(111.51) and k == pytest.approx(526.3472)
    at = ps.atomtypes["c3"]
    assert round(at.sigma, 6) == 0.339967 and at.epsilon == pytest.approx(0.4577296) and at.mass == 12.01


# rq-13487d5c
def test_multi_term_dihedral(tmp_path):
    path = tmp_path / "p.dat"
    path.write_text(dat_with(dihedrals="X -c -c -X    1    0.5       180.0            -2.\n"
                                       "X -c -c -X    1    0.3         0.0             3."))
    terms = read_amber_parameters(path).dihedraltypes[("X", "c", "c", "X")]
    assert [(p, round(k, 6), n) for p, k, n in terms] == [(180.0, 2.092, 2), (0.0, 1.2552, 3)]


# rq-ed920f0f
def test_improper_conversion(tmp_path):
    path = tmp_path / "p.dat"
    path.write_text(dat_with(impropers="c -o -o -o          1.1          180.          2."))
    phase, k, n = read_amber_parameters(path).impropertypes[("c", "o", "o", "o")]
    assert (phase, round(k, 6), n) == (180.0, 4.6024, 2)


def methylammonium():
    mol = Chem.AddHs(Chem.MolFromSmiles("C[NH3+]"))
    types = tuple("c3" if a.GetSymbol() == "C" else "n4" if a.GetSymbol() == "N" else
                  ("hc" if a.GetNeighbors()[0].GetSymbol() == "C" else "hn") for a in mol.GetAtoms())
    return mol, types


# rq-2708fa8a
def test_frcmod_overrides_dat(tmp_path):
    dat = tmp_path / "gaff2.dat"
    dat.write_text(dat_with(bonds="c3-n4  300.0   1.500", masses="c3 12.01\nn4 14.01\nhc 1.008\nhn 1.008",
                            nonbond="  c3  1.9  0.1\n  n4  1.8  0.17\n  hc  1.4  0.01\n  hn  0.6  0.01"))
    frcmod = write(tmp_path, "f.frcmod", "frcmod\nBOND\nc3-n4  350.0   1.480\n\n")
    mol, types = methylammonium()
    ps = Gaff2Backend(gaff2_dat=dat)._parameters(dat, frcmod, mol, types)
    assert ps.bondtypes[("c3", "n4")] == pytest.approx((0.148, 2 * 350.0 * 418.4))


# rq-8aa9cd8f
def test_only_fragment_types_kept(tmp_path):
    dat = tmp_path / "gaff2.dat"
    dat.write_text(dat_with(bonds="c3-hc  337.3   1.092\nc3-os  301.5   1.439", angles="hc-c3-hc   39.4   107.58",
                            masses="c3 12.01\nhc 1.008\nos 16.0",
                            nonbond="  c3  1.9  0.1\n  hc  1.4  0.01\n  os  1.7  0.17"))
    frcmod = write(tmp_path, "f.frcmod", "frcmod\n\n")
    mol = methane()
    ps = Gaff2Backend(gaff2_dat=dat)._parameters(dat, frcmod, mol, ("c3", "hc", "hc", "hc", "hc"))
    assert "os" not in ps.atomtypes
    for table in (ps.bondtypes, ps.angletypes, ps.dihedraltypes, ps.impropertypes):
        assert all(t in ("c3", "hc", "X") for key in table for t in key)


# rq-d50c2e10
def test_tool_failure(fake_amber, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_MODE", "fail")
    with pytest.raises(ParametrizationError) as exc:
        Gaff2Backend().parametrize(methane(), 0, tmp_path / "w")
    assert exc.value.tool == "antechamber" and exc.value.returncode == 1 and "sqm failed" in exc.value.log


# rq-d03cde6e
def test_tool_timeout(fake_amber, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_MODE", "hang")
    with pytest.raises(ParametrizationError) as exc:
        Gaff2Backend(timeout=1).parametrize(methane(), 0, tmp_path / "w")
    assert exc.value.returncode == -1


# rq-17e6c2e5
def test_end_to_end_with_fake_ambertools(fake_amber, tmp_path):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3), initiator="KPS"), "gaff2", name="g3")
    path = write_gromacs(top, tmp_path / "out")
    assert path.is_file() and (tmp_path / "out/g3.ff/g3.rtp").is_file()
    for block in top.blocks:
        if block.name.startswith("DM"):
            assert round(sum(round(q * 10_000) for q in block.charges)) == 0
    assert set(top.atom_types) <= {"cg", "hg", "ng", "og", "sg"}


# rq-bf118372
@pytest.mark.skipif(shutil.which("antechamber") is None, reason="AmberTools not installed")
def test_real_ambertools(tmp_path):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3), initiator="KPS"), "gaff2", name="real")
    write_gromacs(top, tmp_path)
    assert abs(sum(top.charges) + 2) < 1e-4


# rq-98ce435e
def test_no_charges_skips_am1bcc(fake_amber, tmp_path):
    _, log = fake_amber
    frag = Gaff2Backend().parametrize(methane(), 0, tmp_path / "w", charges=False)
    ac = next(line for line in log.read_text().splitlines() if line.startswith("antechamber"))
    assert " -c " not in f" {ac} " and "bcc" not in ac
    assert set(frag.charges) == {0.0}
