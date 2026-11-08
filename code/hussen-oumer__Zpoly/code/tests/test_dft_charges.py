"""Tests for DFT charge export/import (rqm/forcefield/dft_charges.md)."""

import pytest
import yaml
from test_forcefield_core import FakeBackend

from zpoly.forcefield import (
    DftChargeError,
    export_dft_inputs,
    load_dft_charges,
    parametrize,
    register_backend,
)
from zpoly.structure import build_polymer, homopolymer
from zpoly.structure._build import _make_template
from zpoly.zwitterions import default_library


class RecordingBackend(FakeBackend):
    """Fake backend that also records the `charges` flag."""

    def __init__(self, name="dftfake"):
        super().__init__(name)
        self.charge_flags = []

    def parametrize(self, mol, net_charge, work_dir, *, charges=True):
        self.charge_flags.append(charges)
        return super().parametrize(mol, net_charge, work_dir)


@pytest.fixture
def rec():
    backend = RecordingBackend()
    register_backend(backend)
    return backend


def manifest(folder):
    return {m["name"]: m for m in yaml.safe_load((folder / "dft_manifest.yaml").read_text())["molecules"]}


def pattern_charges(atoms, net):
    """Distinct per-atom charges (0.001 * index) adjusted to sum exactly to `net`."""
    q = [0.001 * (i + 1) for i in range(len(atoms))]
    q[-1] += net - sum(q)
    return q


def write_plain(folder, name, atoms, charges, with_elements=True):
    lines = [f"{e} {q:.6f}" if with_elements else f"{q:.6f}" for e, q in zip(atoms, charges, strict=True)]
    (folder / f"{name}.charges").write_text("# test charges\n" + "\n".join(lines) + "\n")


def exported(tmp_path, initiator=None):
    folder = tmp_path / "dft"
    export_dft_inputs("DMAPS", folder, initiator=initiator)
    return folder, manifest(folder)


# rq-0a9d9a0e
def test_export_monomer(tmp_path):
    folder, m = exported(tmp_path)
    assert {p.name for p in folder.iterdir()} >= {"DMAPS.xyz", "DMAPS.pdb", "dft_manifest.yaml"}
    tpl = _make_template(default_library().get("DMAPS").repeat_unit_smiles(), seed=0)
    lines = (folder / "DMAPS.xyz").read_text().splitlines()
    assert int(lines[0]) == tpl.mol.GetNumAtoms() and "charge=0" in lines[1]
    caps = m["DMAPS"]["caps"]
    assert len(caps) == 8
    assert sorted(m["DMAPS"]["atoms"][i - 1] for i in caps) == ["C", "C"] + ["H"] * 6


# rq-ca2b6f77
def test_export_initiator(tmp_path):
    folder, m = exported(tmp_path, "KPS")
    lines = (folder / "KPS_end.xyz").read_text().splitlines()
    assert int(lines[0]) == 9 and "charge=-1" in lines[1]
    assert len(m["KPS_end"]["templates"]) == 2 and "KPS_head" not in m and "KPS_tail" not in m


# rq-dea7e8ed
def test_export_keeps_entries(tmp_path):
    folder, _ = exported(tmp_path)
    export_dft_inputs("CBMA", folder)
    assert {"DMAPS", "CBMA"} <= set(manifest(folder))


# rq-5abcd130
def test_plain_charge_list(tmp_path):
    folder, m = exported(tmp_path)
    q = pattern_charges(m["DMAPS"]["atoms"], 0)
    write_plain(folder, "DMAPS", m["DMAPS"]["atoms"], q)
    dft = load_dft_charges(folder)
    assert dft.by_template[m["DMAPS"]["templates"][0]] == pytest.approx(q)
    assert dft.sources["DMAPS"] == "DMAPS.charges"


# rq-1132056a
def test_wrong_count(tmp_path):
    folder, m = exported(tmp_path)
    atoms = m["DMAPS"]["atoms"][:-1]
    write_plain(folder, "DMAPS", atoms, [0.0] * len(atoms))
    with pytest.raises(DftChargeError) as exc:
        load_dft_charges(folder)
    assert exc.value.path.endswith("DMAPS.charges")


# rq-0bf14ca4
def test_element_mismatch(tmp_path):
    folder, m = exported(tmp_path)
    atoms = list(m["DMAPS"]["atoms"])
    assert atoms[0] == "C"
    write_plain(folder, "DMAPS", ["O"] + atoms[1:], [0.0] * len(atoms))
    with pytest.raises(DftChargeError) as exc:
        load_dft_charges(folder)
    assert exc.value.atom == 1


# rq-e24ad818
def test_wrong_total(tmp_path):
    folder, m = exported(tmp_path)
    atoms = m["DMAPS"]["atoms"]
    write_plain(folder, "DMAPS", atoms, pattern_charges(atoms, 1))
    with pytest.raises(DftChargeError, match="net charge 0"):
        load_dft_charges(folder)


# rq-d2d065e8
def test_mol2_charges(tmp_path):
    folder, m = exported(tmp_path)
    atoms = m["DMAPS"]["atoms"]
    q = pattern_charges(atoms, 0)
    rows = [f"{i:7d} {e}{i:<4d} 0.0 0.0 0.0 {e.lower()}x 1 MOL {c:10.6f}"
            for i, (e, c) in enumerate(zip(atoms, q, strict=True), 1)]
    (folder / "DMAPS.mol2").write_text("@<TRIPOS>MOLECULE\nm\n\n@<TRIPOS>ATOM\n" + "\n".join(rows)
                                       + "\n@<TRIPOS>BOND\n")
    assert load_dft_charges(folder).by_template[m["DMAPS"]["templates"][0]] == pytest.approx(q)


def qm_rows(atoms, charges, fmt):
    return "\n".join(fmt.format(i=i, e=e, q=q) for i, (e, q) in enumerate(zip(atoms, charges, strict=True), 1))


# rq-108f0cc1
def test_gaussian_esp(tmp_path):
    folder, m = exported(tmp_path)
    atoms = m["DMAPS"]["atoms"]
    mull, esp = [0.0] * len(atoms), pattern_charges(atoms, 0)
    row = "{i:6d}  {e:<2s}  {q:10.6f}"
    text = (" Mulliken charges:\n               1\n" + qm_rows(atoms, mull, row)
            + "\n Sum of Mulliken charges =   0.00000\n"
            " ESP charges:\n               1\n" + qm_rows(atoms, esp, row) + "\n Sum of ESP charges =   0.00000\n"
            " ESP charges with hydrogens summed into heavy atoms:\n               1\n     1  C   0.100000\n")
    (folder / "DMAPS.log").write_text(text)
    assert load_dft_charges(folder).by_template[m["DMAPS"]["templates"][0]] == pytest.approx(esp)


# rq-b599705d
def test_orca_chelpg(tmp_path):
    folder, m = exported(tmp_path)
    atoms = m["DMAPS"]["atoms"]
    mull, chelpg = [0.0] * len(atoms), pattern_charges(atoms, 0)
    text = ("-----------------------\nMULLIKEN ATOMIC CHARGES\n-----------------------\n"
            + qm_rows(atoms, mull, "{i:4d} {e:<2s}:{q:12.6f}") + "\nSum of atomic charges:   0.0000000\n\n"
            "CHELPG Charges            \n--------------------------------\n"
            + qm_rows(atoms, chelpg, "{i:4d}   {e:<2s}  :{q:14.6f}") + "\n--------------------------------\n"
            "Total charge:     0.000000\n")
    (folder / "DMAPS.out").write_text(text)
    assert load_dft_charges(folder).by_template[m["DMAPS"]["templates"][0]] == pytest.approx(chelpg)


# rq-b90af318
def test_missing_charge_file(tmp_path):
    folder, _ = exported(tmp_path)
    assert "DMAPS" in load_dft_charges(folder).missing


def full_dft(tmp_path):
    folder, m = exported(tmp_path, "KPS")
    for name in ("DMAPS", "KPS_end"):
        write_plain(folder, name, m[name]["atoms"], pattern_charges(m[name]["atoms"], m[name]["net_charge"]))
    return folder, m


# rq-90737175
def test_dft_charges_mapped(rec, tmp_path):
    folder, m = full_dft(tmp_path)
    structure = build_polymer(homopolymer("DMAPS", 5), initiator="KPS")
    top = parametrize(structure, "dftfake", dft_charges=folder)
    dft = load_dft_charges(folder).by_template[m["DMAPS"]["templates"][0]]
    mol = structure.mol
    diffs = [top.charges[a.GetIdx()] - dft[a.GetIntProp("template_atom")]
             for a in mol.GetAtoms() if a.GetProp("res") == "DMR"]
    assert max(diffs) - min(diffs) <= 2e-4  # one common shift (plus 4-decimal rounding)
    for block in top.blocks:
        assert sum(round(q * 10_000) for q in block.charges) == block.charge * 10_000
    assert rec.charge_flags and not any(rec.charge_flags)
    assert set(top.dft_blocks) == {"KPH", "DMF", "DMR", "DME", "KPT"}


# rq-fa62c105
def test_partial_dft_falls_back(rec, tmp_path):
    folder, m = exported(tmp_path, "KPS")
    write_plain(folder, "DMAPS", m["DMAPS"]["atoms"], pattern_charges(m["DMAPS"]["atoms"], 0))
    top = parametrize(build_polymer(homopolymer("DMAPS", 3), initiator="KPS"), "dftfake", dft_charges=folder)
    assert set(top.dft_blocks) == {"DMF", "DMR", "DME"}
    assert top.dft_missing == ("KPS_end",)
    assert all(rec.charge_flags)


# rq-a8df24ac
def test_charge_file_after_dft(rec, tmp_path):
    folder, m = exported(tmp_path)
    write_plain(folder, "DMAPS", m["DMAPS"]["atoms"], pattern_charges(m["DMAPS"]["atoms"], 0))
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "dftfake", dft_charges=folder,
                      charges={"DMR": {"N1": 0.2}})
    dmr = next(b for b in top.blocks if b.name == "DMR")
    assert dict(zip(dmr.atom_names, dmr.charges, strict=True))["N1"] == 0.2


# rq-684b318c
def test_cli_export_check_build(rec, run_inp, tmp_path, capsys):
    assert run_inp("JOB = dft\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nDFT = dft") == 0
    m = manifest(tmp_path / "dft")
    for name in ("DMAPS", "KPS_end"):
        write_plain(tmp_path / "dft", name, m[name]["atoms"], pattern_charges(m[name]["atoms"], m[name]["net_charge"]))
    capsys.readouterr()
    assert run_inp("MONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = dftfake\nDFT = dft\nOUTPUT = out/d") == 0
    assert "DFT charges: 5 residue blocks" in capsys.readouterr().out


# rq-4d98df6a
def test_cli_notice_incomplete(rec, run_inp, tmp_path, capsys):
    folder, m = exported(tmp_path, "KPS")
    write_plain(folder, "DMAPS", m["DMAPS"]["atoms"], pattern_charges(m["DMAPS"]["atoms"], 0))
    assert run_inp(f"MONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = dftfake\nDFT = {folder}\nOUTPUT = out/d") == 0
    assert "note: no DFT charges for KPS_end" in capsys.readouterr().err


# rq-10a45821
def test_gaff2_is_default(run_inp, tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PATH", str(tmp_path))  # no AmberTools
    assert run_inp("MONOMER = DMAPS\nN = 3") == 1
    assert "AmberTools" in capsys.readouterr().err
    assert run_inp("MONOMER = DMAPS\nN = 3\nFF = none", name="b.inp") == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == ["DMAPS_n3.pdb", "b.inp", "job.inp"]
