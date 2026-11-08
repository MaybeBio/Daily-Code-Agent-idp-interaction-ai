"""Tests for the force-field core (rqm/forcefield/core.md), using a deterministic fake backend."""

import filecmp
import re
from pathlib import Path

import pytest
from rdkit import Chem

from zpoly.cli import main
from zpoly.forcefield import (
    AtomType,
    ForceFieldError,
    ForceFieldUnavailableError,
    FragmentParameters,
    MissingParameterError,
    ParameterSet,
    UnknownForceFieldError,
    available_backends,
    get_backend,
    parametrize,
    register_backend,
    write_gromacs,
)
from zpoly.forcefield._core import bonded_terms
from zpoly.structure import build_polymer, homopolymer

MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974}


class FakeBackend:
    """Types atoms by element, charge 0.1 each, parameters for every combination present."""

    def __init__(self, name="fake", *, wildcard_dihedrals=False, drop_angle=None, conflict_on_call=None,
                 available=True, drop_bond=None):
        self.name = name
        self.wildcard_dihedrals = wildcard_dihedrals
        self.drop_angle = drop_angle
        self.conflict_on_call = conflict_on_call
        self.available = available
        self.drop_bond = drop_bond
        self.calls = []

    def check_available(self):
        if not self.available:
            raise ForceFieldUnavailableError(self.name, "fake tools missing", "install the fake tools")

    def parametrize(self, mol, net_charge, work_dir):
        self.calls.append(net_charge)
        types = tuple(a.GetSymbol() for a in mol.GetAtoms())
        bonds, angles, propers, _ = bonded_terms(mol)

        def key(term):
            return tuple(types[i] for i in term)

        atomtypes = {}
        for a in mol.GetAtoms():
            mass = MASSES[a.GetSymbol()]
            if self.conflict_on_call == len(self.calls) and a.GetSymbol() == "C":
                mass += 1.0
            atomtypes[a.GetSymbol()] = AtomType(a.GetSymbol(), a.GetAtomicNum(), mass, 0.3, 0.4)
        dropped = {self.drop_bond, tuple(reversed(self.drop_bond or ()))}
        bondtypes = {key(b): (0.15, 250000.0) for b in bonds if key(b) not in dropped}
        angletypes = {key(a): (109.5, 400.0) for a in angles if key(a) not in (self.drop_angle,
                                                                               tuple(reversed(self.drop_angle or ())))}
        if self.wildcard_dihedrals:
            dihedraltypes = {("X", key(p)[1], key(p)[2], "X"): ((0.0, 1.0, 3),) for p in propers}
        else:
            dihedraltypes = {key(p): ((0.0, 1.0, 3),) for p in propers}
        impropertypes = {("X", "X", "C", "O"): (180.0, 4.6, 2)}
        params = ParameterSet((1, 2, "yes", 0.5, 0.8333), (1, 1, 9, 4, 1, 3, 1, 0), atomtypes, bondtypes,
                              angletypes, dihedraltypes, impropertypes)
        return FragmentParameters(types, tuple(0.1 for _ in types), params)


@pytest.fixture
def fake():
    backend = FakeBackend()
    register_backend(backend)
    return backend


@pytest.fixture(scope="module")
def kps5():
    return build_polymer(homopolymer("DMAPS", 5), initiator="KPS")


@pytest.fixture(scope="module")
def kps3():
    return build_polymer(homopolymer("DMAPS", 3), initiator="KPS")


def residue_atoms(top, resnum):
    return [a.GetIdx() for a in top.structure.mol.GetAtoms() if a.GetIntProp("resnum") == resnum]


# rq-be5f3279
def test_blocks_of_kps_homopolymer(fake, kps5):
    top = parametrize(kps5, "fake")
    assert [b.name for b in top.blocks] == ["KPH", "DMF", "DMR", "DME", "KPT"]


# rq-e2efc78d
def test_repeat_residues_identical(fake, kps5):
    top = parametrize(kps5, "fake")
    dmr = next(b for b in top.blocks if b.name == "DMR")
    mol = top.structure.mol
    resnums = sorted({a.GetIntProp("resnum") for a in mol.GetAtoms() if a.GetProp("res") == "DMR"})
    assert len(resnums) == 3
    for r in resnums:
        idx = residue_atoms(top, r)
        assert tuple(mol.GetAtomWithIdx(i).GetProp("name") for i in idx) == dmr.atom_names
        assert tuple(top.atom_types[i] for i in idx) == dmr.atom_types
        assert tuple(top.charges[i] for i in idx) == dmr.charges


# rq-998720ac
def test_block_charges_normalized(fake, kps5):
    top = parametrize(kps5, "fake")
    for b in top.blocks:
        total = round(sum(round(q * 10_000) for q in b.charges))
        assert total == (-10_000 if b.name in ("KPH", "KPT") else 0), b.name
    assert round(sum(top.charges), 4) == -2


# rq-882e553a
def test_context_oligomers(fake, kps5):
    parametrize(kps5, "fake")
    assert sorted(fake.calls) == [-1, -1, 0]


# rq-71d21eb1
def test_one_unit_chain(fake):
    top = parametrize(build_polymer(["DMAPS"]), "fake")
    assert len(fake.calls) == 1
    assert [b.name for b in top.blocks] == ["DMF"]


# rq-c3c66893
def test_copolymer_blocks(fake):
    top = parametrize(build_polymer(["DMAPS", "CBMA", "DMAPS", "CBMA"]), "fake")
    assert [b.name for b in top.blocks] == ["DMF", "CBR", "DMR", "CBE"]


# rq-5ddf514c
def test_conflicting_atom_types(kps5):
    register_backend(FakeBackend("conflict", conflict_on_call=2))
    with pytest.raises(ForceFieldError, match="'C'"):
        parametrize(kps5, "conflict")


# rq-26506f0a
def test_missing_bonded_parameter():
    register_backend(FakeBackend("noangle", drop_angle=("C", "S", "O")))
    with pytest.raises(MissingParameterError) as exc:
        parametrize(build_polymer(homopolymer("DMAPS", 3)), "noangle")
    assert exc.value.kind == "angle"
    assert ("C", "S", "O") in exc.value.types or ("O", "S", "C") in exc.value.types


# rq-37d656d6
def test_wildcard_dihedrals_match():
    register_backend(FakeBackend("wild", wildcard_dihedrals=True))
    parametrize(build_polymer(homopolymer("DMAPS", 3)), "wild")


@pytest.fixture
def files3(fake, kps3, tmp_path):
    top = parametrize(kps3, "fake", name="dmaps3")
    write_gromacs(top, tmp_path)
    return top, tmp_path


def sections(text, name):
    """Lines of every '[ name ]' section (comments and blanks removed)."""
    out, current = [], None
    for line in text.splitlines():
        m = re.match(r"\s*\[\s*(\S+)\s*\]", line)
        if m:
            current = m.group(1)
            continue
        if current == name and line.strip() and not line.strip().startswith(";"):
            out.append(line.split())
    return out


# rq-ba87a649
def test_write_gromacs_files(files3):
    top, d = files3
    for rel in ["dmaps3.pdb", "dmaps3.top", "dmaps3.ff/forcefield.itp", "dmaps3.ff/forcefield.doc",
                "dmaps3.ff/ffnonbonded.itp", "dmaps3.ff/ffbonded.itp", "dmaps3.ff/dmaps3.rtp"]:
        assert (d / rel).is_file(), rel
    ff_itp = (d / "dmaps3.ff/forcefield.itp").read_text()
    assert ff_itp.index("[ defaults ]") < ff_itp.index('#include "ffnonbonded.itp"') < ff_itp.index(
        '#include "ffbonded.itp"')
    top_text = (d / "dmaps3.top").read_text()
    assert '#include "dmaps3.ff/forcefield.itp"' in top_text
    atoms = sections(top_text, "atoms")
    pdb_atoms = [line for line in (d / "dmaps3.pdb").read_text().splitlines() if line.startswith("ATOM")]
    assert len(atoms) == len(pdb_atoms)
    for row, line in zip(atoms, pdb_atoms, strict=True):
        assert row[3] == line[17:20].strip() and row[4] == line[12:16].strip()


def rtp_blocks(text):
    blocks, current, sub = {}, None, None
    for line in text.splitlines():
        m = re.match(r"\[\s*(\S+)\s*\]", line)
        if m:
            current, sub = m.group(1), None
            blocks.setdefault(current, {})
            continue
        m = re.match(r"\s+\[\s*(\S+)\s*\]", line)
        if m:
            sub = m.group(1)
            blocks[current].setdefault(sub, [])
            continue
        if sub and line.strip() and not line.strip().startswith(";"):
            blocks[current][sub].append(line.split())
    return blocks


# rq-fecefe8b
def test_rtp_blocks_and_links(files3):
    _, d = files3
    blocks = rtp_blocks((d / "dmaps3.ff/dmaps3.rtp").read_text())
    assert [k for k in blocks if k != "bondedtypes"] == ["KPH", "DMF", "DMR", "DME", "KPT"]
    assert not any(a.startswith("-") for row in blocks["KPH"]["bonds"] for a in row)
    links = [row for row in blocks["DMF"]["bonds"] if any(a.startswith("-") for a in row)]
    assert len(links) == 1 and links[0][0].startswith("-O")
    for name, block in blocks.items():
        if name == "bondedtypes":
            continue
        for row in block["atoms"]:
            assert len(row) == 4
            float(row[2])
            int(row[3])


# rq-e3e93468
def test_top_bonded_sections(fake, tmp_path):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake", name="d3")
    write_gromacs(top, tmp_path)
    text = (tmp_path / "d3.top").read_text()
    bonds = {tuple(sorted(map(int, r[:2]))) for r in sections(text, "bonds")}
    assert len(bonds) == top.structure.mol.GetNumBonds()
    close = set(bonds) | {tuple(sorted((a + 1, c + 1))) for a, _, c in top.angles}
    assert not any(tuple(sorted(map(int, r[:2]))) in close for r in sections(text, "pairs"))
    funcs = {r[4] for r in sections(text, "dihedrals")}
    assert funcs <= {"9", "4"} and "9" in funcs
    assert "4" in funcs  # the carbonyl carbon gets an improper from the fake ("X","X","C","O") type


# rq-0e0b7ea5
def test_net_charge_comment(files3):
    _, d = files3
    assert "; net charge -2: add 2 K+" in (d / "dmaps3.top").read_text()


# rq-7348427c
def test_deterministic_output(files3, tmp_path):
    top, d = files3
    other = tmp_path / "again"
    write_gromacs(top, other)
    for path in d.rglob("*"):
        if path.is_file() and other not in path.parents:
            twin = other / path.relative_to(d)
            assert filecmp.cmp(path, twin, shallow=False), path.name


# rq-23d1bfb8
def test_unknown_force_field(kps3):
    with pytest.raises(UnknownForceFieldError) as exc:
        parametrize(kps3, "charmm")
    assert "gaff2" in exc.value.available


# rq-d876e1b0
def test_register_backend_at_runtime():
    mine = FakeBackend("mine")
    register_backend(mine)
    assert "mine" in available_backends()
    assert get_backend("MINE") is mine


# rq-ba2adece
def test_cli_writes_ff_files(fake, run_inp, tmp_path, capsys):
    assert run_inp("MONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = fake\nOUTPUT = out/dm") == 0
    assert (tmp_path / "out/dm.pdb").exists() and (tmp_path / "out/dm.top").exists()
    assert (tmp_path / "out/dm.ff").is_dir()
    assert "force field fake" in capsys.readouterr().out


# rq-07b1387b
def test_cli_unavailable_writes_nothing(run_inp, tmp_path, capsys):
    register_backend(FakeBackend("absent", available=False))
    assert run_inp("MONOMER = DMAPS\nN = 3\nFF = absent\nOUTPUT = out/dm") == 1
    assert capsys.readouterr().err.startswith("error:")
    assert not (tmp_path / "out").exists() or not any(Path(tmp_path / "out").iterdir())


# rq-dca9b6e9
def test_cli_lists_forcefields(capsys):
    assert main(["list", "forcefields"]) == 0
    assert any(line.startswith("gaff2") for line in capsys.readouterr().out.splitlines())


def test_rdkit_reads_generated_pdb(files3):
    _, d = files3
    assert Chem.MolFromPDBFile(str(d / "dmaps3.pdb"), removeHs=False, sanitize=False) is not None


# User charges (QM / RESP)
def block(top, name):
    return next(b for b in top.blocks if b.name == name)


def raw_dmr_charge(atom_name):
    return 0.1  # the fake backend gives every atom 0.1


# rq-56f6f132
def test_user_charges_fixed_and_balanced(fake, kps5):
    top = parametrize(kps5, "fake", charges={"DMR": {"N1": 0.5, "S1": 1.2}}, builtin_charges=False)
    dmr = block(top, "DMR")
    q = dict(zip(dmr.atom_names, dmr.charges, strict=True))
    assert q["N1"] == 0.5 and q["S1"] == 1.2
    assert sum(round(c * 10_000) for c in dmr.charges) == 0
    others = [c for n, c in q.items() if n not in ("N1", "S1")]
    assert max(others) - min(others) <= 0.0001 + 1e-12
    report = next(r for r in top.charge_report if r.residue == "DMR")
    assert report.user_atoms == 2


# rq-d3724017
def test_user_charges_in_every_repeat(fake, kps5):
    top = parametrize(kps5, "fake", charges={"DMR": {"N1": 0.5}})
    mol = top.structure.mol
    hits = [top.charges[a.GetIdx()] for a in mol.GetAtoms() if a.GetProp("res") == "DMR" and a.GetProp("name") == "N1"]
    assert len(hits) == 3 and all(q == 0.5 for q in hits)


# rq-e3a90ed3
def test_unknown_residue_ignored(fake, kps3):
    top = parametrize(kps3, "fake", charges={"XXR": {"C1": 0.1}, "DMR": {"N1": 0.5}})
    dmr = block(top, "DMR")
    assert dict(zip(dmr.atom_names, dmr.charges, strict=True))["N1"] == 0.5


# rq-ab5245a5
def test_unknown_atom_name(fake, kps3):
    from zpoly.forcefield import ChargeFileError
    with pytest.raises(ChargeFileError) as exc:
        parametrize(kps3, "fake", charges={"DMR": {"N99": 0.5}})
    assert exc.value.residue == "DMR" and exc.value.atom == "N99"


# rq-1ea0268f
def test_full_residue_wrong_total(fake, kps3):
    from zpoly.forcefield import ChargeFileError
    kph = {"O1": -0.1, "S1": -0.1, "O2": -0.1, "O3": -0.1, "O4": -0.1}  # sums to -0.5, not -1
    with pytest.raises(ChargeFileError) as exc:
        parametrize(kps3, "fake", charges={"KPH": kph})
    assert exc.value.residue == "KPH"


# rq-e4d525a8
def test_large_correction_flagged(fake, kps3):
    top = parametrize(kps3, "fake", charges={"DMR": {"N1": 1.0, "S1": 2.0}})
    assert next(r for r in top.charge_report if r.residue == "DMR").large


# rq-58d9d01e
def test_invalid_charge_file(tmp_path):
    from zpoly.forcefield import ChargeFileError, load_charge_file
    path = tmp_path / "bad.yaml"
    path.write_text("DMR: [1, 2]\n")
    with pytest.raises(ChargeFileError) as exc:
        load_charge_file(path)
    assert str(path) in str(exc.value)


# rq-4700787b
def test_charge_template(kps3, tmp_path):
    from zpoly.forcefield import charge_template, load_charge_file
    text = charge_template(kps3)
    order = [line[:-1] for line in text.splitlines() if re.match(r"^[A-Z0-9]{3}:$", line)]
    assert order == ["KPH", "DMF", "DMR", "DME", "KPT"]
    mol = kps3.mol
    heavy, hydrogens = set(), set()
    for a in mol.GetAtoms():
        key = (a.GetProp("res"), a.GetProp("name"))
        (hydrogens if a.GetAtomicNum() == 1 else heavy).add(key)
    commented = re.findall(r"^  # (\S+): .*?   # (\w+),", text, re.M)
    assert len(commented) == len(heavy)
    named_h = re.findall(r"\d+ H \(([^)]*)\)", text)
    assert sum(len(group.split()) for group in named_h) == len(hydrogens)
    uncommented = re.sub(r"^  # ", "  ", text, flags=re.M)
    path = tmp_path / "all.yaml"
    path.write_text(uncommented)
    loaded = load_charge_file(path)
    assert set(loaded) == {"KPH", "DMF", "DMR", "DME", "KPT"}


def test_fully_commented_template_means_no_user_charges(kps3, tmp_path):
    from zpoly.forcefield import charge_template, load_charge_file
    path = tmp_path / "t.yaml"
    path.write_text(charge_template(kps3))
    assert load_charge_file(path) == {}


# rq-83e95ea1
def test_cli_charges_option(fake, run_inp, tmp_path, capsys):
    (tmp_path / "q.yaml").write_text("DMR:\n  N1: 0.5\n")
    assert run_inp("MONOMER = DMAPS\nN = 3\nFF = fake\nCHARGES = q.yaml\nBUILTIN_CHARGES = no\nOUTPUT = out/d") == 0
    assert "user charges: DMR 1 atoms" in capsys.readouterr().out
    rows = [line.split() for line in (tmp_path / "out/d.top").read_text().splitlines()
            if re.match(r"^\s+\d+\s+\S+\s+\d+\s+DMR\s+N1\s", line)]
    assert rows and all(r[6] == "0.5000" for r in rows)


# rq-3618722f
def test_cli_charges_template(run_inp, tmp_path):
    assert run_inp("JOB = charges\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = none") == 0
    out = (tmp_path / "DMAPS_n3_KPS_charges.yaml").read_text()
    for name in ("KPH:", "DMF:", "DMR:", "DME:", "KPT:"):
        assert f"\n{name}\n" in out or out.startswith(name)



# A heavy atom and its hydrogens in one entry
def hydrogens_of(top, res, atom):
    mol = top.structure.mol
    r = next(a.GetIntProp("resnum") for a in mol.GetAtoms() if a.GetProp("res") == res)
    parent = next(a for a in mol.GetAtoms() if a.GetIntProp("resnum") == r and a.GetProp("name") == atom)
    return [n.GetProp("name") for n in parent.GetNeighbors() if n.GetAtomicNum() == 1]


# rq-0bd185f8
def test_compact_hydrogen_entry(fake):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake", builtin_charges=False,
                      charges={"DMR": {"C7": {"charge": -0.4018, "H": 0.2042}, "N1": 0.1573}})
    dmr = block(top, "DMR")
    q = dict(zip(dmr.atom_names, dmr.charges, strict=True))
    hs = hydrogens_of(top, "DMR", "C7")
    assert len(hs) == 3
    assert q["C7"] == -0.4018 and q["N1"] == 0.1573 and all(q[h] == 0.2042 for h in hs)
    assert next(r for r in top.charge_report if r.residue == "DMR").user_atoms == 5
    assert sum(round(c * 10_000) for c in dmr.charges) == 0


# rq-27bf349e
def test_only_hydrogens_of_atom(fake):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake", charges={"DMR": {"C7": {"H": 0.1}}},
                      builtin_charges=False)
    dmr = block(top, "DMR")
    q = dict(zip(dmr.atom_names, dmr.charges, strict=True))
    hs = hydrogens_of(top, "DMR", "C7")
    assert all(q[h] == 0.1 for h in hs)
    others = [c for n, c in q.items() if n not in hs]
    assert max(others) - min(others) <= 0.0001 + 1e-12  # C7 balanced like every other unlisted atom


# rq-c5f83c02
def test_hydrogen_given_twice(fake):
    from zpoly.forcefield import ChargeFileError
    structure = build_polymer(homopolymer("DMAPS", 3))
    probe = parametrize(structure, "fake")
    h = hydrogens_of(probe, "DMR", "C7")[0]
    with pytest.raises(ChargeFileError) as exc:
        parametrize(structure, "fake", charges={"DMR": {"C7": {"H": 0.2}, h: 0.3}})
    assert exc.value.atom == h


# rq-d93f2f21
def test_h_on_atom_without_hydrogens(fake):
    from zpoly.forcefield import ChargeFileError
    with pytest.raises(ChargeFileError) as exc:
        parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake",
                    charges={"DMR": {"N1": {"charge": 0.15, "H": 0.2}}})
    assert exc.value.atom == "N1"


# rq-a2452483
def test_invalid_compact_entry(tmp_path):
    from zpoly.forcefield import ChargeFileError, load_charge_file
    path = tmp_path / "bad.yaml"
    path.write_text("DMR:\n  C7: {q: 0.1}\n")
    with pytest.raises(ChargeFileError) as exc:
        load_charge_file(path)
    assert exc.value.atom == "C7"


# User bonded parameters (typically the ionic parts) and the ionic template
from zpoly.forcefield import (  # noqa: E402
    ParameterFileError,
    ionic_parameter_template,
    load_parameter_file,
)


def write_params(tmp_path, text, name="p.yaml"):
    path = tmp_path / name
    path.write_text(text)
    return path


# rq-38038e4b
def test_missing_parameters_template(tmp_path):
    register_backend(FakeBackend("gaps", drop_angle=("C", "S", "O"), drop_bond=("C", "S")))
    with pytest.raises(MissingParameterError) as exc:
        parametrize(build_polymer(homopolymer("DMAPS", 3)), "gaps")
    assert set(exc.value.missing) >= {"bond", "angle"}
    loaded = load_parameter_file(write_params(tmp_path, exc.value.template))
    assert loaded["bond"][0]["types"] in (("C", "S"), ("S", "C"))


# rq-3072b091
def test_user_parameters_fill_gap(tmp_path):
    register_backend(FakeBackend("noangle2", drop_angle=("C", "S", "O")))
    path = write_params(tmp_path, "version: 1\nangles:\n  - {types: [C, S, O], theta0: 108.9, k: 619.232}\n")
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "noangle2", params=path)
    assert top.parameters.angletypes[("C", "S", "O")] == (108.9, 619.232)
    assert ("angle", ("C", "S", "O")) in top.user_parameters


# rq-96dc51b6
def test_same_as_copies(tmp_path):
    register_backend(FakeBackend("nobond", drop_bond=("C", "S")))
    path = write_params(tmp_path, "version: 1\nbonds:\n  - {types: [C, S], same_as: [C, C]}\n")
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "nobond", params=path)
    assert top.parameters.bondtypes[("C", "S")] == (0.15, 250000.0)


# rq-74842904
def test_user_parameters_override(fake, tmp_path):
    path = write_params(tmp_path, "version: 1\nbonds:\n  - {types: [C, C], b0: 0.2, kb: 1000.0}\n")
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake", params=path, name="o")
    assert top.parameters.bondtypes[("C", "C")] == (0.2, 1000.0)
    write_gromacs(top, tmp_path / "out")
    line = next(x for x in (tmp_path / "out/o.ff/ffbonded.itp").read_text().splitlines() if x.split()[:2] == ["C", "C"])
    assert line.endswith("; user")


# rq-ba3008a1
def test_wrong_dihedral_style(fake, tmp_path):
    path = write_params(tmp_path, "version: 1\ndihedrals:\n  - {types: [X, C, C, X], rb: [0, 0, 0, 0, 0, 0]}\n")
    with pytest.raises(ParameterFileError) as exc:
        parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake", params=path)
    assert exc.value.entry == "dihedrals[1]"


# rq-db609516
def test_same_as_missing_key(fake, tmp_path):
    path = write_params(tmp_path, "version: 1\nbonds:\n  - {types: [C, S], same_as: [Q, Q]}\n")
    with pytest.raises(ParameterFileError) as exc:
        parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake", params=path)
    assert exc.value.entry == "bonds[1]"


# rq-ff8c4f34
def test_invalid_parameter_file(tmp_path):
    path = write_params(tmp_path, "version: 1\nbonds:\n  - {types: [C], b0: 0.1}\n")
    with pytest.raises(ParameterFileError) as exc:
        load_parameter_file(path)
    assert str(path) in str(exc.value) and exc.value.entry == "bonds[1]"


# rq-0d2848e7
def test_cli_parameter_file(fake, run_inp, tmp_path, capsys):
    path = write_params(tmp_path, "version: 1\nbonds:\n  - {types: [C, C], b0: 0.2, kb: 1000.0}\n")
    assert run_inp(f"MONOMER = DMAPS\nN = 3\nFF = fake\nPARAMS = {path}\nOUTPUT = out/d") == 0
    assert "user parameters: 1" in capsys.readouterr().out


# rq-38b4e442
def test_ionic_template_lists_ionic_parts(fake, tmp_path):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "fake")
    text = ionic_parameter_template(top)
    for needle in ("{types: [C, N],", "{types: [O, S],", "{types: [O, S, O],"):
        assert needle in text, needle
    assert "{types: [O, C, O]," not in text  # the ester carbonyl is not part of an ionic group
    reapplied = parametrize(top.structure, "fake", params=write_params(tmp_path, text))
    from zpoly.forcefield._core import match_exact
    for attr in ("bondtypes", "angletypes"):  # same value for every key, in either orientation
        for key, value in getattr(top.parameters, attr).items():
            assert match_exact(getattr(reapplied.parameters, attr), key) == pytest.approx(value, rel=1e-6), key


# rq-0df95ac2
def test_template_marks_missing(tmp_path):
    register_backend(FakeBackend("nobondcs", drop_bond=("C", "S")))
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "nobondcs", strict=False)
    text = ionic_parameter_template(top)
    line = next(x for x in text.splitlines() if "{types: [C, S]," in x)
    assert "b0: 0.0, kb: 0.0" in line and "# missing in the force field" in line
    with pytest.raises(MissingParameterError):
        write_gromacs(top, tmp_path / "w")


# rq-d76b560f
def test_cli_params_template(fake, run_inp, tmp_path):
    assert run_inp("JOB = params\nMONOMER = DMAPS\nN = 3\nFF = fake") == 0
    out = (tmp_path / "DMAPS_n3_params.yaml").read_text()
    assert "bonds:" in out and "angles:" in out
