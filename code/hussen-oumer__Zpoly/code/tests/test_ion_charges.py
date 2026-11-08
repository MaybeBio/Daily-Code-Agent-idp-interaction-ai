"""Ion fragments for DFT and charge scaling (rqm/forcefield/dft_charges.md)."""

import pytest
import yaml
from rdkit import Chem
from test_dft_charges import manifest, pattern_charges, rec, write_plain  # noqa: F401  (rec: fixture)
from test_forcefield_core import FakeBackend
from test_system_slab import solvent  # noqa: F401  (fixture)

from zpoly.forcefield import ionic_groups, load_dft_charges, parametrize, register_backend
from zpoly.structure import build_polymer, homopolymer
from zpoly.structure._build import _make_template
from zpoly.zwitterions import default_library


def template(name):
    return _make_template(default_library().get(name).repeat_unit_smiles(), seed=0).mol


def heavy(mol, group):
    return sorted(mol.GetAtomWithIdx(i).GetSymbol() for i in group.atoms if mol.GetAtomWithIdx(i).GetAtomicNum() > 1)


# rq-013e1066
def test_dmaps_groups():
    mol = template("DMAPS")
    cation, anion = ionic_groups(mol)
    assert cation.charge == 1 and heavy(mol, cation) == ["C", "C", "C", "C", "N"]
    assert anion.charge == -1 and heavy(mol, anion) == ["C", "O", "O", "O", "S"]


# rq-39505fe6
def test_merged_groups():
    (vbd,) = ionic_groups(template("VBD"))
    assert vbd.charge == 2
    one_ch2 = Chem.AddHs(Chem.MolFromSmiles("C=CC(=O)OCC[N+](C)(C)CS(=O)(=O)[O-]"))
    (pair,) = ionic_groups(one_ch2)
    assert pair.charge == 0


# rq-cefbfada
def test_export_fragments(run_inp, tmp_path):
    assert run_inp("JOB = dft\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nDFT = q\nDFT_PARTS = ions") == 0
    m = manifest(tmp_path / "q")
    assert set(m) == {"DMAPS", "DMAPS_cation", "DMAPS_anion", "DMAPS_ester", "KPS_end"}
    cat, an = m["DMAPS_cation"], m["DMAPS_anion"]
    assert cat["net_charge"] == 1 and len(cat["atoms"]) == 17 and sorted(set(cat["atoms"])) == ["C", "H", "N"]
    assert an["net_charge"] == -1 and sorted(an["atoms"]) == ["C", "H", "H", "H", "O", "O", "O", "S"]
    parent = m["DMAPS"]["atoms"]
    for frag in (cat, an):
        caps = {c for c, _ in frag["caps"]}
        for k, (element, p) in enumerate(zip(frag["atoms"], frag["parent_atoms"], strict=True), 1):
            assert (p == 0 and k in caps and element == "H") or (p > 0 and parent[p - 1] == element)
    assert "charge=1" in (tmp_path / "q/DMAPS_cation.xyz").read_text().splitlines()[1]
    assert "*id" not in (tmp_path / "q/dft_manifest.yaml").read_text()
    ester = m["DMAPS_ester"]  # methyl pivalate; only C(=O)O-CH2 (+2 H) take its charges
    assert ester["net_charge"] == 0 and len(ester["atoms"]) == 20 and len(ester["apply"]) == 6
    assert sorted(ester["atoms"][k - 1] for k in ester["apply"]) == ["C", "C", "H", "H", "O", "O"]


def write_all(folder, names):
    m = manifest(folder)
    for name in names:
        write_plain(folder, name, m[name]["atoms"], pattern_charges(m[name]["atoms"], m[name]["net_charge"]))


@pytest.fixture
def ion_folder(run_inp, tmp_path):
    assert run_inp("JOB = dft\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nDFT = q\nDFT_PARTS = ions") == 0
    return tmp_path / "q"


# rq-f8346d02
def test_fragment_charges_fixed(rec, ion_folder):  # noqa: F811
    write_all(ion_folder, ["DMAPS", "DMAPS_cation", "DMAPS_anion", "DMAPS_ester", "KPS_end"])
    structure = build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    top = parametrize(structure, "dftfake", dft_charges=ion_folder, dft_parts="ions")
    assert set(top.ion_blocks) == {"DMF", "DMR", "DME"} and not top.dft_missing
    assert not set(top.builtin_blocks) & {"DMF", "DMR", "DME"}  # DFT takes precedence over built-in charges
    ions = load_dft_charges(ion_folder).ions_by_template[default_library().get("DMAPS").repeat_unit_smiles()]
    assert len(ions) == 17 + 8 - 3 + 6  # cation + anion atoms (caps folded in) + 6 linker atoms
    checked = 0
    for atom in structure.mol.GetAtoms():
        k = atom.GetIntProp("template_atom")
        if atom.GetProp("res").startswith("DM") and k in ions:
            assert top.charges[atom.GetIdx()] == pytest.approx(ions[k], abs=1.5e-4)
            checked += 1
    assert checked == 3 * len(ions)
    for block in top.blocks:
        assert sum(block.charges) == pytest.approx(block.charge, abs=1e-9)


# rq-7670d71d
def test_missing_fragments(rec, ion_folder):  # noqa: F811
    write_all(ion_folder, ["DMAPS", "KPS_end"])
    structure = build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    top = parametrize(structure, "dftfake", dft_charges=ion_folder, dft_parts="ions")
    assert {"DMAPS_cation", "DMAPS_anion"} <= set(top.dft_missing) and not top.ion_blocks
    whole = parametrize(structure, "dftfake", dft_charges=ion_folder)
    assert not whole.dft_missing


# rq-e14e60c1
def test_scaling_chain():
    register_backend(FakeBackend("fake"))
    structure = build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    plain = parametrize(structure, "fake")
    top = parametrize(structure, "fake", charge_scale=0.8)
    assert top.charge_scale == 0.8
    sums = {b.name: sum(b.charges) for b in top.blocks}
    assert sums["KPH"] == pytest.approx(-0.8) and sums["KPT"] == pytest.approx(-0.8)
    assert all(sums[r] == pytest.approx(0.0, abs=1e-9) for r in ("DMF", "DMR", "DME"))
    ionic = {i for g in ionic_groups(structure.mol) for i in g.atoms}
    for atom in structure.mol.GetAtoms():
        i = atom.GetIdx()
        if i in ionic and atom.GetProp("res").startswith("DM"):
            assert top.charges[i] == pytest.approx(0.8 * plain.charges[i], abs=1.5e-4)
    with pytest.raises(Exception, match="CHARGE_SCALE"):
        parametrize(structure, "fake", charge_scale=1.5)


# rq-db0975d0
def test_scaling_system(run_inp, solvent, tmp_path):  # noqa: F811
    register_backend(FakeBackend("fake"))
    assert run_inp(f"JOB = bulk\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = fake\nCHAINS = 2\nBOX = 3\n"
                   f"PACKER = builtin\nSOLVENT_FF = {solvent}\nCHARGE_SCALE = 0.8\nOUTPUT = b") == 0
    ions = (tmp_path / "b/polymer.ff/ions.itp").read_text()
    k_line = next(line for line in ions.splitlines() if line.split()[1:2] == ["K"])
    assert float(k_line.split()[6]) == pytest.approx(0.8)
    assert "charges scaled by 0.8" in (tmp_path / "b/system.top").read_text()
    assert yaml.safe_load((tmp_path / "b/protocol.yaml").read_text())["charge_scale"] == 0.8


# rq-c7326013
def test_zwitterionic_liquid_export(run_inp, tmp_path):
    assert run_inp("JOB = dft\nMOLECULE = TEA-PS\nDFT_PARTS = ions\nDFT = z") == 0
    assert set(manifest(tmp_path / "z")) == {"TEA-PS", "TEA-PS_cation", "TEA-PS_anion"}
    assert run_inp("JOB = dft\nMOLECULE = TEA-PS\nN = 3", name="bad.inp") == 1


# rq-a3491212
def test_bad_scale(run_inp, capsys):
    for value in ("0", "1.5"):
        assert run_inp(f"MONOMER = DMAPS\nN = 3\nFF = none\nCHARGE_SCALE = {value}") == 1
        assert "CHARGE_SCALE must be in (0, 1]" in capsys.readouterr().err
