"""Tests for the OPLS-AA backend (rqm/forcefield/oplsaa.md) with a fixture oplsaa.ff folder."""

import itertools
import os
import re
import shutil
from pathlib import Path

import pytest
from rdkit import Chem

from zpoly.forcefield import (
    ForceFieldUnavailableError,
    OplsaaBackend,
    TypingError,
    load_typing_rules,
    parametrize,
    read_gromacs_forcefield,
    register_backend,
    write_gromacs,
)
from zpoly.structure import build_polymer, homopolymer

# name, bond type, atomic number, mass, charge (fixture values, GROMACS layout)
FIXTURE_TYPES = [
    ("opls_103", "N3", 7, 14.0067, 0.0), ("opls_135", "CT", 6, 12.011, -0.18),
    ("opls_136", "CT", 6, 12.011, -0.12), ("opls_137", "CT", 6, 12.011, -0.06),
    ("opls_139", "CT", 6, 12.011, 0.0), ("opls_140", "HC", 1, 1.008, 0.06),
    ("opls_465", "C", 6, 12.011, 0.51), ("opls_466", "O", 8, 15.9994, -0.43),
    ("opls_467", "OS", 8, 15.9994, -0.33), ("opls_490", "CT", 6, 12.011, 0.19),
    ("opls_493", "SY", 16, 32.06, 1.2), ("opls_494", "OY", 8, 15.9994, -0.6),
    ("opls_235", "C", 6, 12.011, 0.5), ("opls_236", "O", 8, 15.9994, -0.5),
    ("opls_238", "N", 7, 14.0067, -0.5), ("opls_241", "H", 1, 1.008, 0.3),
    ("opls_145", "CA", 6, 12.011, -0.115), ("opls_146", "HA", 1, 1.008, 0.115),
    ("opls_154", "OH", 8, 15.9994, -0.683), ("opls_155", "HO", 1, 1.008, 0.418),
    ("opls_180", "OS", 8, 15.9994, -0.4), ("opls_271", "C", 6, 12.011, 0.7),
    ("opls_272", "O2", 8, 15.9994, -0.8), ("opls_440", "P", 15, 30.97376, 1.4),
    ("opls_441", "O2", 8, 15.9994, -0.8), ("opls_442", "OS", 8, 15.9994, -0.6),
    ("opls_142", "CM", 6, 12.011, -0.115), ("opls_143", "CM", 6, 12.011, -0.23),
    ("opls_144", "HC", 1, 1.008, 0.115),
]
# atomtypes.atp descriptions (fixture wording, containing the expect keywords where relevant)
FIXTURE_DESCRIPTIONS = {
    "opls_145": "Benzene C - 12 site", "opls_146": "Benzene H - 12 site",
    "opls_154": "all-atom O: mono alcohols", "opls_155": "all-atom H(O): mono alcohols",
    "opls_180": "O: dialkyl ether", "opls_271": "C in COO- carboxylate", "opls_272": "O in COO- carboxylate",
    "opls_440": "P in Me2PO4-", "opls_441": "O= in Me2PO4-", "opls_442": "OMe in Me2PO4-",
    "opls_142": "alkene C (RH-C=)", "opls_143": "alkene C (H2-C=)", "opls_144": "alkene H (H-C=)",
}
BOND_TYPES = sorted({t[1] for t in FIXTURE_TYPES})


def write_fixture_ff(folder: Path, descriptions: dict | None = None) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    if descriptions is not False:
        desc = {**FIXTURE_DESCRIPTIONS, **(descriptions or {})}
        (folder / "atomtypes.atp").write_text("".join(
            f"  {name:<10s}{mass:10.5f}  ; {desc.get(name, 'custom type')}\n" for name, _, _, mass, _ in FIXTURE_TYPES))
    (folder / "forcefield.itp").write_text(
        "; fixture OPLS-AA\n[ defaults ]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n"
        "1  3  yes  0.5  0.5\n\n#include \"ffnonbonded.itp\"\n#include \"ffbonded.itp\"\n")
    lines = ["[ atomtypes ]", "; name  bond_type  at.num  mass  charge  ptype  sigma  epsilon"]
    for name, bt, z, mass, q in FIXTURE_TYPES:
        lines.append(f" {name}   {bt:<3s} {z:3d} {mass:10.4f} {q:8.3f}  A  3.50000e-01  2.76144e-01 ; {name}")
    (folder / "ffnonbonded.itp").write_text("\n".join(lines) + "\n")
    out = ["#ifdef HEAVY_H", "; nothing here", "#endif", "[ bondtypes ]"]
    for a, b in itertools.combinations_with_replacement(BOND_TYPES, 2):
        b0, kb = (0.1529, 224262.4) if (a, b) == ("CT", "CT") else (0.14, 250000.0)
        out.append(f"  {a:<3s} {b:<3s}  1  {b0:.4f}  {kb:.1f}   ; bond comment")
    out.append("[ angletypes ]")
    for a, c in itertools.combinations_with_replacement(BOND_TYPES, 2):
        for b in BOND_TYPES:
            out.append(f"  {a:<3s} {b:<3s} {c:<3s}  1  109.500  400.000")
    out.append("[ dihedraltypes ]")
    for b, c in itertools.combinations_with_replacement(BOND_TYPES, 2):
        out.append(f"  X  {b:<3s} {c:<3s} X   3   0.62760  1.88280  0.00000 -2.51040  0.00000  0.00000")
    out += ["#define improper_O_C_X_Y   180.0   43.93200   2",
            "#define improper_Z_N_X_Y   180.0    4.18400   2",
            "#define improper_Z_CA_X_Y  180.0    4.60240   2",
            "#define some_flag"]
    (folder / "ffbonded.itp").write_text("\n".join(out) + "\n")
    return folder


@pytest.fixture
def ff_dir(tmp_path):
    return write_fixture_ff(tmp_path / "oplsaa.ff")


@pytest.fixture
def opls(ff_dir):
    backend = OplsaaBackend(oplsaa_ff=ff_dir)
    register_backend(backend)
    return backend


def typed(backend, structure, tmp_path):
    frag = backend.parametrize(structure.mol, structure.chain.net_charge, tmp_path)
    return structure.mol, frag


def types_where(mol, frag, smarts):
    patt = Chem.MolFromSmarts(smarts)
    return {frag.atom_types[m[0]] for m in mol.GetSubstructMatches(patt)}


# rq-c593d1ee
def test_not_found(monkeypatch, tmp_path):
    for var in ("ZPOLY_OPLSAA_FF", "GMXDATA"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ForceFieldUnavailableError) as exc:
        OplsaaBackend().check_available()
    assert "OPLSAA_FF" in exc.value.hint


# rq-a641bad1
def test_found_from_environment(ff_dir, monkeypatch):
    monkeypatch.setenv("ZPOLY_OPLSAA_FF", str(ff_dir))
    OplsaaBackend().check_available()


# rq-de0e6dd9
def test_read_gromacs_files(ff_dir):
    ff = read_gromacs_forcefield(ff_dir)
    bt, z, mass, q, _, _ = ff.atomtypes["opls_135"]
    assert (bt, z, mass, q) == ("CT", 6, 12.011, -0.18)
    assert ff.bondtypes[("CT", "CT")] == (0.1529, 224262.4)
    assert len(ff.rb_dihedraltypes[("X", "CT", "CT", "X")]) == 6
    assert ff.defines["improper_O_C_X_Y"] == (180.0, 43.932, 2)
    assert ff.defaults == (1, 3, "yes", 0.5, 0.5)


# rq-daf83e1c
def test_preprocessor_and_comments(ff_dir):
    text = (ff_dir / "ffbonded.itp").read_text()
    assert "#ifdef" in text and "; bond comment" in text
    ff = read_gromacs_forcefield(ff_dir)
    assert ("C", "OS") in ff.bondtypes or ("OS", "C") in ff.bondtypes
    assert "some_flag" not in ff.defines


# rq-3b60921d
def test_typing_dmaps(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("DMAPS", 3)), tmp_path)
    assert types_where(mol, frag, "[NX4+]") == {"opls_103"}
    assert types_where(mol, frag, "[SX4]") == {"opls_493"}
    assert types_where(mol, frag, "[OX1]~[SX4]") == {"opls_494"}
    assert types_where(mol, frag, "[CX3](=O)O") == {"opls_465"}
    assert types_where(mol, frag, "[OX1]=[CX3]") == {"opls_466"}
    assert types_where(mol, frag, "[OX2]([CX3]=O)C") == {"opls_467"}
    assert types_where(mol, frag, "[CH2][OX2][CX3]=O") == {"opls_490"}
    assert types_where(mol, frag, "[CX4;H0]") == {"opls_139"}
    assert types_where(mol, frag, "[#1][#6]") == {"opls_140"}


# rq-7fac93b5
def test_typing_acrylamide(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("A3361", 3)), tmp_path)
    assert types_where(mol, frag, "[CX3](=O)N") == {"opls_235"}
    assert types_where(mol, frag, "[OX1]=[CX3]N") == {"opls_236"}
    assert types_where(mol, frag, "[NX3;H1]C=O") == {"opls_238"}
    assert types_where(mol, frag, "[#1][NX3]") == {"opls_241"}
    assert types_where(mol, frag, "[CX4;H1](C)C(=O)N") == {"opls_137"}


# rq-cffa5145
def test_kps_end_groups(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("DMAPS", 2), initiator="KPS"), tmp_path)
    assert types_where(mol, frag, "[SX4]([OX2])(=O)(=O)[O-]") == {"opls_493"}
    assert types_where(mol, frag, "[OX1][SX4][OX2]") == {"opls_494"}
    assert types_where(mol, frag, "[OX2]([#6])[SX4]") == {"opls_467"}


# rq-08faba23
def test_default_charges_are_type_charges(opls, tmp_path):
    charges = {name: q for name, _, _, _, q in FIXTURE_TYPES}
    _, frag = typed(opls, build_polymer(homopolymer("DMAPS", 2)), tmp_path)
    assert all(q == charges[t] for t, q in zip(frag.atom_types, frag.charges, strict=True))


# rq-e02bf89b
def test_untyped_atoms_reported(opls, tmp_path):
    with pytest.raises(TypingError) as exc:
        typed(opls, build_polymer(homopolymer("SPVI", 2)), tmp_path)
    assert exc.value.atoms
    assert any(" N" in a and "bonded to" in a for a in exc.value.atoms)
    assert "OPLS_RULES" in str(exc.value)


# rq-ea516f0e
def test_user_rules_override(ff_dir, tmp_path):
    rules = tmp_path / "rules.yaml"
    rules.write_text('version: 1\nrules:\n  - {name: N-methyl, smarts: "[CX4;H3][NX4+]", type: opls_136}\n')
    backend = OplsaaBackend(oplsaa_ff=ff_dir, rules=rules)
    mol, frag = typed(backend, build_polymer(homopolymer("DMAPS", 2)), tmp_path)
    assert types_where(mol, frag, "[CH3][NX4+]") == {"opls_136"}
    assert types_where(mol, frag, "[CH3][CX4;H0]") == {"opls_135"}


# rq-c7ee0327
def test_rule_type_missing_from_ff(ff_dir, tmp_path):
    rules = tmp_path / "rules.yaml"
    rules.write_text('version: 1\nrules:\n  - {name: odd methyl, smarts: "[CX4;H3]", type: opls_999}\n')
    with pytest.raises(TypingError) as exc:
        OplsaaBackend(oplsaa_ff=ff_dir, rules=rules).check_available()
    assert exc.value.rule == "odd methyl"


# rq-a2b7965f
def test_invalid_rule_file(tmp_path):
    rules = tmp_path / "bad.yaml"
    rules.write_text('version: 1\nrules:\n  - {name: no type, smarts: "[CX4]"}\n')
    with pytest.raises(TypingError) as exc:
        load_typing_rules(rules)
    assert str(rules) in str(exc.value)


# rq-91413f6c
def test_impropers_from_defines(opls, tmp_path):
    _, frag = typed(opls, build_polymer(homopolymer("DMAPS", 2)), tmp_path)
    assert frag.parameters.impropertypes[("X", "X", "C", "X")] == (180.0, 43.932, 2)


# rq-010b9f97
def test_rb_dihedrals_in_output(opls, tmp_path):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "oplsaa", name="o3")
    write_gromacs(top, tmp_path)
    bonded = (tmp_path / "o3.ff/ffbonded.itp").read_text()
    rb_lines = [line.split() for line in bonded.splitlines() if re.match(r"^\s+\S+\s+\S+\s+\S+\s+\S+\s+3\s", line)]
    assert rb_lines and all(len(r) == 11 for r in rb_lines)
    top_text = (tmp_path / "o3.top").read_text()
    propers = top_text.split("; propers")[1].split("[ ")[0]
    rows = [r.split() for r in propers.splitlines() if r.strip() and not r.strip().startswith(";")]
    assert rows and all(r[4] == "3" for r in rows)
    rtp = (tmp_path / "o3.ff/o3.rtp").read_text().splitlines()
    assert rtp[2].split() == ["1", "1", "3", "4", "1", "3", "1", "0"]
    nb = [line.split() for line in (tmp_path / "o3.ff/ffnonbonded.itp").read_text().splitlines()
          if line.strip().startswith("opls_")]
    assert nb and all(len(r) == 8 and r[1] in BOND_TYPES for r in nb)


# rq-34b90470
def test_bonded_matched_by_bond_type(opls):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "oplsaa")
    assert all(not k[0].startswith("opls_") for k in top.parameters.bondtypes)


# rq-b2b4afa0
def test_charge_file_with_opls(opls):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "oplsaa", charges={"DMR": {"N1": 0.1573}})
    mol = top.structure.mol
    n1 = [top.charges[a.GetIdx()] for a in mol.GetAtoms() if a.GetProp("res") == "DMR" and a.GetProp("name") == "N1"]
    assert n1 and all(q == 0.1573 for q in n1)
    for block in top.blocks:
        assert sum(round(q * 10_000) for q in block.charges) == 0


# rq-0e9e4781
def test_cli_with_opls(ff_dir, run_inp, tmp_path):
    assert run_inp(f"MONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = oplsaa\nOPLSAA_FF = {ff_dir}\nOUTPUT = out/d") == 0
    assert (tmp_path / "out/d.top").is_file() and (tmp_path / "out/d.ff").is_dir()


# rq-db5f1df5
def test_cli_opls_options_need_oplsaa(ff_dir, run_inp, capsys):
    assert run_inp(f"MONOMER = DMAPS\nN = 3\nOPLSAA_FF = {ff_dir}") == 1
    assert "FF = oplsaa" in capsys.readouterr().err


def _real_oplsaa():
    if os.environ.get("GMXDATA"):
        cand = Path(os.environ["GMXDATA"]) / "top" / "oplsaa.ff"
        if (cand / "ffnonbonded.itp").is_file():
            return cand
    gmx = shutil.which("gmx") or shutil.which("gmx_mpi")
    if gmx:
        cand = Path(gmx).resolve().parent.parent / "share" / "gromacs" / "top" / "oplsaa.ff"
        if (cand / "ffnonbonded.itp").is_file():
            return cand
    return None


# rq-f129ffeb
@pytest.mark.skipif(_real_oplsaa() is None, reason="no GROMACS oplsaa.ff found")
def test_real_oplsaa():
    ff = read_gromacs_forcefield(_real_oplsaa())
    assert len(ff.atomtypes) >= 800 and ff.atomtypes["opls_135"][0] == "CT"


# rq-e2e215ea
def test_typing_carboxybetaine(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("CBMA", 2)), tmp_path)
    assert types_where(mol, frag, "[CX3](=O)[O-]") == {"opls_271"}
    assert types_where(mol, frag, "[OX1]~[CX3]~[OX1]") == {"opls_272"}


# rq-87d16db3
def test_typing_alcohol_and_ether(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("HSBMA", 2)), tmp_path / "a")
    assert types_where(mol, frag, "[OX2;H1]") == {"opls_154"}
    assert types_where(mol, frag, "[#1][OX2]") == {"opls_155"}
    mol, frag = typed(opls, build_polymer(homopolymer("ZB4", 2)), tmp_path / "b")
    assert types_where(mol, frag, "[OX2;R]([CX4])[CX4]") == {"opls_180"}


# rq-8f9bad27
def test_typing_phosphorylcholine(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("MPC", 2)), tmp_path)
    assert types_where(mol, frag, "[PX4]") == {"opls_440"}
    assert types_where(mol, frag, "[OX1]~[PX4]") == {"opls_441"}
    assert types_where(mol, frag, "[OX2]([#6])[PX4]") == {"opls_442"}
    patt = Chem.MolFromSmarts("[OX1]~[PX4]")
    assert len({m[0] for m in mol.GetSubstructMatches(patt)}) == 4  # two terminal O per unit


# rq-c823b8ce
def test_typing_benzene_ring(opls, tmp_path):
    mol, frag = typed(opls, build_polymer(homopolymer("SBVB", 1)), tmp_path)
    ring_c = [i for i, a in enumerate(mol.GetAtoms()) if a.GetIsAromatic()]
    ring_h = [n.GetIdx() for i in ring_c for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
    assert len(ring_c) == 6 and {frag.atom_types[i] for i in ring_c} == {"opls_145"}
    assert len(ring_h) == 4 and {frag.atom_types[i] for i in ring_h} == {"opls_146"}


# rq-3d69e50c
def test_expect_mismatch(tmp_path):
    folder = write_fixture_ff(tmp_path / "bad.ff", {"opls_271": "alkane CH3"})
    with pytest.raises(TypingError) as exc:
        OplsaaBackend(oplsaa_ff=folder).check_available()
    assert exc.value.rule == "carboxylate C" and "alkane CH3" in str(exc.value)


# rq-2106839f
def test_no_atomtypes_atp(tmp_path, caplog):
    folder = write_fixture_ff(tmp_path / "noatp.ff", descriptions=False)
    assert not (folder / "atomtypes.atp").exists()
    with caplog.at_level("WARNING"):
        OplsaaBackend(oplsaa_ff=folder).check_available()
    assert "expect" in caplog.text and "skipped" in caplog.text


# rq-9dd1d715
def test_conditional_atomtype_variants(ff_dir):
    nb = ff_dir / "ffnonbonded.itp"
    heavy = " opls_140   HC    1     4.0320    0.060  A  2.5e-01  1.2552e-01"
    text = nb.read_text().replace(" opls_140   HC ", f"#ifdef HEAVY_H\n{heavy}\n#else\n opls_140   HC ", 1)
    lines = text.splitlines()
    i = next(k for k, line in enumerate(lines) if line.startswith(" opls_140   HC ") and "1.0080" in line)
    lines.insert(i + 1, "#endif")
    nb.write_text("\n".join(lines) + "\n")
    assert read_gromacs_forcefield(ff_dir).atomtypes["opls_140"][2] == 1.008
