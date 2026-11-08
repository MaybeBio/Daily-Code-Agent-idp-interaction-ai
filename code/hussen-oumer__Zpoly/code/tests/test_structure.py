"""Tests for the polymer chain builder (rqm/polymer_structure.md)."""

import textwrap
from collections import Counter

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from zpoly.cli import main
from zpoly.structure import (
    PdbLimitError,
    PolymerStructure,
    SequenceError,
    alternating,
    block,
    build_polymer,
    homopolymer,
    pattern,
    random_copolymer,
    write_pdb,
)
from zpoly.zwitterions import UnknownEntryError, build_chain, load_library


def residues(structure):
    seen = {}
    for atom in structure.mol.GetAtoms():
        seen.setdefault(atom.GetIntProp("resnum"), atom.GetProp("res"))
    return [seen[k] for k in sorted(seen)]


def heavy_smiles(structure):
    return Chem.MolToSmiles(Chem.RemoveHs(structure.mol))


@pytest.fixture(scope="module")
def kps_trimer():
    return build_polymer(homopolymer("DMAPS", 3), initiator="KPS")


# rq-35976189
def test_homopolymer_sequence():
    assert homopolymer("DMAPS", 50) == ["DMAPS"] * 50


# rq-95c60185
def test_alternating_odd():
    assert alternating("DMAPS", "CBMA", 5) == ["DMAPS", "CBMA", "DMAPS", "CBMA", "DMAPS"]


# rq-2aa4244e
def test_block_sequence():
    assert block([("DMAPS", 2), ("CBMA", 3)]) == ["DMAPS"] * 2 + ["CBMA"] * 3


# rq-c073d054
def test_random_copolymer_counts_and_seed():
    a = random_copolymer({"ZB2": 0.56, "SB": 0.44}, 30, seed=1)
    assert a == random_copolymer({"ZB2": 0.56, "SB": 0.44}, 30, seed=1)
    assert Counter(a) == {"ZB2": 17, "SB": 13}
    assert Counter(random_copolymer({"ZB2": 0.56, "SB": 0.44}, 30, seed=2)) == {"ZB2": 17, "SB": 13}


# rq-d87d9e74
def test_pattern_sequence():
    assert pattern("AAB", {"A": "DMAPS", "B": "CBMA"}, 7) == [
        "DMAPS", "DMAPS", "CBMA", "DMAPS", "DMAPS", "CBMA", "DMAPS"]


# rq-826f003e
@pytest.mark.parametrize("call", [
    lambda: homopolymer("DMAPS", 0),
    lambda: block([]),
    lambda: block([("DMAPS", 0)]),
    lambda: random_copolymer({"DMAPS": 0.5, "CBMA": 0.4}, 10),
    lambda: random_copolymer({"DMAPS": 1.2, "CBMA": -0.2}, 10),
    lambda: pattern("", {"A": "DMAPS"}, 3),
    lambda: pattern("AB", {"A": "DMAPS"}, 3),
])
def test_invalid_sequences(call):
    with pytest.raises(SequenceError):
        call()


# rq-06424eaf
def test_build_50_unit_homopolymer_with_kps():
    s = build_polymer(homopolymer("DMAPS", 50), initiator="KPS")
    assert residues(s) == ["KPH", "DMF"] + ["DMR"] * 48 + ["DME", "KPT"]
    ref = build_chain(["DMAPS"] * 50, initiator="KPS")
    assert heavy_smiles(s) == ref.smiles
    assert s.chain.net_charge == -2
    assert s.min_nonbonded_distance >= 1.5


# rq-6d35271c
def test_build_30_unit_random_copolymer():
    seq = random_copolymer({"DMAPS": 0.5, "CBMA": 0.5}, 30, seed=3)
    s = build_polymer(seq)
    names = residues(s)
    assert len(names) == 30 and {n[:2] for n in names} == {"DM", "CB"}
    assert names[0][2] == "F" and names[-1][2] == "E" and all(n[2] == "R" for n in names[1:-1])
    assert heavy_smiles(s) == build_chain(seq).smiles
    assert s.chain.net_charge == 0


# rq-0660a451
def test_coordinates_and_hydrogens():
    s = build_polymer(homopolymer("DMAPS", 5))
    mol = s.mol
    assert mol.GetNumConformers() == 1
    assert rdMolDescriptors.CalcMolFormula(mol) == s.chain.formula
    assert heavy_smiles(s) == s.chain.smiles
    xyz = mol.GetConformer().GetPositions()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() > 1 and b.GetAtomicNum() > 1:
            d = np.linalg.norm(xyz[a.GetIdx()] - xyz[b.GetIdx()])
            assert 1.1 <= d <= 2.0, (a.GetSymbol(), b.GetSymbol(), d)


# rq-d497f481
def test_same_seed_same_coordinates():
    a = build_polymer(homopolymer("DMAPS", 10), seed=5)
    b = build_polymer(homopolymer("DMAPS", 10), seed=5)
    assert np.array_equal(a.mol.GetConformer().GetPositions(), b.mol.GetConformer().GetPositions())


# rq-3e9f3397
def test_extended_longer_than_coil():
    ext = build_polymer(homopolymer("DMAPS", 20), conformation="extended")
    coil = build_polymer(homopolymer("DMAPS", 20), conformation="coil")
    assert ext.end_to_end_distance > coil.end_to_end_distance


# rq-7f3fc73b
def test_unknown_conformation():
    with pytest.raises(ValueError):
        build_polymer(homopolymer("DMAPS", 3), conformation="helix")


# rq-61c06a74
def test_unknown_monomer_in_sequence():
    with pytest.raises(UnknownEntryError):
        build_polymer(["DMAPS", "NOPE"])


# rq-050546f1
def test_residue_name_collision():
    assert residues(build_polymer(["CBMA", "CBMA-1"])) == ["CBF", "C1E"]


# rq-c1a4dfaa
def test_library_residue_field(tmp_path):
    path = tmp_path / "lib.yaml"
    path.write_text(textwrap.dedent("""
        version: 1
        monomers:
          - name: MYSB
            full_name: my sulfobetaine
            family: Sulfobetaine
            status: specialty
            residue: MS
            blocks: [MA, C2, NMe2, C4, SO3]
    """), encoding="utf-8")
    lib = load_library(path)
    assert residues(build_polymer(["MYSB"] * 3, library=lib)) == ["MSF", "MSR", "MSE"]


# rq-8994396e
def test_uncapped_ends_hydrogen_terminated():
    s = build_polymer(homopolymer("DMAPS", 3))
    first_res = [a.GetIdx() for a in s.mol.GetAtoms() if a.GetIntProp("resnum") == 1]
    patt = Chem.MolFromSmarts("[CH3;!$(C[N+]);!$(CO)][CX4H0](C)(C(=O)O)")
    mol_h = Chem.RemoveHs(s.mol)
    # RemoveHs keeps atom order of heavy atoms; map by residue property.
    hits = [m for m in mol_h.GetSubstructMatches(patt) if mol_h.GetAtomWithIdx(m[0]).GetIntProp("resnum") == 1]
    assert hits and first_res
    assert s.chain.net_charge == 0


def _parse_pdb(path):
    atoms, conect, cryst = [], {}, None
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            atoms.append({"serial": int(line[6:11]), "name": line[12:16].strip(), "res": line[17:20],
                          "resnum": int(line[22:26]),
                          "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                          "element": line[76:78].strip()})
        elif line.startswith("CONECT"):
            fields = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
            conect.setdefault(fields[0], []).extend(fields[1:])
        elif line.startswith("CRYST1"):
            cryst = (float(line[6:15]), float(line[15:24]), float(line[24:33]))
    return atoms, conect, cryst


# rq-cfb82d2e
def test_write_pdb(tmp_path, kps_trimer):
    path = write_pdb(kps_trimer, tmp_path / "out" / "dmaps.pdb")
    assert path.exists()
    text = path.read_text()
    atoms, conect, cryst = _parse_pdb(path)
    assert len(atoms) == kps_trimer.mol.GetNumAtoms()
    assert cryst is not None and conect
    assert "\nTER\n" in text and text.rstrip().endswith("END")
    resnums = [a["resnum"] for a in atoms]
    assert resnums == sorted(resnums)
    by_serial = {a["serial"]: a for a in atoms}
    for i, a in enumerate(atoms):
        if a["element"] == "H":
            parent = conect[a["serial"]][0]
            prev = atoms[i - 1]
            assert prev["serial"] == parent or (prev["element"] == "H" and conect[prev["serial"]][0] == parent)
            assert by_serial[parent]["element"] != "H"
    back = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False, proximityBonding=False)
    assert back.GetNumAtoms() == kps_trimer.mol.GetNumAtoms()
    assert back.GetNumBonds() == kps_trimer.mol.GetNumBonds()


# rq-92413a81
def test_pdb_box_and_centring(tmp_path, kps_trimer):
    path = write_pdb(kps_trimer, tmp_path / "box.pdb", padding=10.0)
    atoms, _, cryst = _parse_pdb(path)
    xyz = np.array([a["xyz"] for a in atoms])
    extent = kps_trimer.mol.GetConformer().GetPositions()
    extent = extent.max(axis=0) - extent.min(axis=0)
    assert np.allclose(cryst, extent + 20.0, atol=2e-3)
    assert np.all(xyz >= 10.0 - 2e-3) and np.all(xyz <= np.array(cryst) - 10.0 + 2e-3)


# rq-8b61e864
def test_pdb_atom_limit(tmp_path, kps_trimer):
    big = Chem.RWMol()
    for _ in range(100_000):
        big.AddAtom(Chem.Atom(2))
    fake = PolymerStructure(mol=big.GetMol(), chain=kps_trimer.chain, sequence=())
    target = tmp_path / "big.pdb"
    with pytest.raises(PdbLimitError):
        write_pdb(fake, target)
    assert not target.exists()


# rq-08692f3c
def test_cli_homopolymer(run_inp, tmp_path, capsys):
    assert run_inp("MONOMER = DMAPS\nN = 50\nINITIATOR = KPS\nFF = none\nOUTPUT = dmaps") == 0
    assert (tmp_path / "dmaps.pdb").exists()
    summary = capsys.readouterr().out
    assert "50 units" in summary and "net charge -2" in summary


# rq-b5e5f3b3
def test_cli_default_output_name(run_inp, tmp_path):
    assert run_inp("COPOLYMER = alternating DMAPS CBMA\nN = 6\nFF = none") == 0
    assert (tmp_path / "DMAPS-alt-CBMA_n6.pdb").exists()


# rq-e83ffc40
def test_cli_random_and_block(run_inp, tmp_path):
    assert run_inp("COPOLYMER = random DMAPS:0.5 CBMA:0.5\nN = 10\nSEED = 2\nFF = none") == 0
    assert (tmp_path / "DMAPS-ran-CBMA_n10_s2.pdb").exists()
    assert run_inp("COPOLYMER = block DMAPS:3 CBMA:2\nFF = none") == 0
    assert (tmp_path / "DMAPS3-CBMA2.pdb").exists()


# rq-a405b58f
def test_cli_library_error(run_inp, capsys):
    assert run_inp("MONOMER = NOPE\nN = 5\nFF = none") == 1
    err = capsys.readouterr().err
    assert err.startswith("error:") and "NOPE" in err


# rq-94f765ba
def test_cli_list_initiators(capsys):
    assert main(["list", "initiators"]) == 0
    assert any(line.startswith("KPS") for line in capsys.readouterr().out.splitlines())



# rq-35c1237e
def test_cli_spacer_and_alkyl(run_inp, tmp_path):
    assert run_inp("MONOMER = DMAPS\nN = 5\nSPACER = 5\nALKYL = 2\nFF = none") == 0
    path = tmp_path / "DMAPS_n5_sp5_alk2.pdb"
    assert path.exists()
    mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False, proximityBonding=False)
    five = Chem.MolFromSmarts("[#7][#6][#6][#6][#6][#6][#16]")
    # The read-back is unsanitized, so hydrogens are explicit atoms: spell out CH3-CH2-N.
    ethyl = Chem.MolFromSmarts("[#6]([#1])([#1])([#1])[#6]([#1])([#1])[#7]")
    assert len(mol.GetSubstructMatches(five)) == 5
    assert len(mol.GetSubstructMatches(ethyl)) == 10  # two ethyl groups on each of the 5 nitrogens


# rq-32089525
def test_cli_variant_error(run_inp, tmp_path, capsys):
    assert run_inp("COPOLYMER = alternating DMAPS SPVI\nN = 4\nALKYL = 2\nFF = none") == 1
    assert capsys.readouterr().err.startswith("error:")
    assert not list(tmp_path.glob("*.pdb"))
