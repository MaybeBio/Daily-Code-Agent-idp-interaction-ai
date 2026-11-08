"""Built-in RESP charges, sulfonate parameters and methyl chain caps (rqm/forcefield/builtin_charges.md)."""

import re

import pytest
from test_dft_charges import manifest, pattern_charges, rec, write_plain  # noqa: F401  (rec: fixture)
from test_forcefield_oplsaa import write_fixture_ff

from zpoly.forcefield import OplsaaBackend, export_dft_inputs, parametrize, register_backend
from zpoly.forcefield._gromacs import moleculetype_text
from zpoly.structure import build_polymer, homopolymer


@pytest.fixture(scope="module")
def opls(tmp_path_factory):
    backend = OplsaaBackend(oplsaa_ff=write_fixture_ff(tmp_path_factory.mktemp("ff") / "oplsaa.ff"))
    register_backend(backend)
    return backend


def charges_of(top, res):
    block = next(b for b in top.blocks if b.name == res)
    return dict(zip(block.atom_names, block.charges, strict=True)), block


# rq-5508e9fc
def test_reproduce_dmaps(opls):
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "oplsaa", charge_scale=0.85)
    q, _ = charges_of(top, "DMR")
    expected = {"N1": 0.1337, "C7": -0.3415, "C6": -0.1680, "C11": -0.3987, "S1": 1.1919, "O3": -0.6398,
                "C4": 0.2170, "O1": -0.2918, "O2": -0.2421, "C5": 0.1934}
    for name, value in expected.items():
        assert q[name] == pytest.approx(value, abs=1.5e-4), name
    assert set(top.builtin_blocks) == {"DMF", "DMR", "DME"}
    for block in top.blocks:
        assert sum(block.charges) == pytest.approx(0.0, abs=1e-9)


# rq-df8c138b
def test_kps_end_groups(opls):
    structure = build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    top = parametrize(structure, "oplsaa", charge_scale=0.85)
    for res in ("KPH", "KPT"):
        q, block = charges_of(top, res)
        assert q["O1"] == pytest.approx(-0.3614, abs=1.5e-4) and q["S1"] == pytest.approx(0.5956, abs=1.5e-4)
        assert all(q[o] == pytest.approx(-0.3614, abs=1.5e-4) for o in ("O2", "O3", "O4"))
        assert sum(block.charges) == pytest.approx(-0.85, abs=1e-9)
    unscaled = parametrize(structure, "oplsaa")
    for atom in structure.mol.GetAtoms():  # the backbone carbon bonded to O1 is not scaled
        if atom.GetProp("res") == "KPH" and atom.GetProp("name") == "O1":
            c = next(n for n in atom.GetNeighbors() if n.GetProp("res") == "DMF").GetIdx()
            assert top.charges[c] == pytest.approx(unscaled.charges[c], abs=1.5e-4)


# rq-f9f878b0
def test_switch_off(opls):
    structure = build_polymer(homopolymer("DMAPS", 3))
    off = parametrize(structure, "oplsaa", builtin_charges=False)
    on = parametrize(structure, "oplsaa")
    assert not off.builtin_blocks and on.builtin_blocks
    assert charges_of(off, "DMR")[0]["N1"] != pytest.approx(charges_of(on, "DMR")[0]["N1"])


# rq-44135a00
def test_dft_takes_precedence(rec, tmp_path):  # noqa: F811
    folder = tmp_path / "dft"
    export_dft_inputs("DMAPS", folder)
    m = manifest(folder)
    write_plain(folder, "DMAPS", m["DMAPS"]["atoms"], pattern_charges(m["DMAPS"]["atoms"], 0))
    top = parametrize(build_polymer(homopolymer("DMAPS", 3)), "dftfake", dft_charges=folder)
    assert set(top.dft_blocks) == {"DMF", "DMR", "DME"} and not top.builtin_blocks


def stock_like_ff(tmp_path):
    """Fixture oplsaa.ff with the sulfonate S typed SY2 and, as in the stock files, no SY2 bonds or angles."""
    folder = write_fixture_ff(tmp_path / "oplsaa.ff")
    nb = (folder / "ffnonbonded.itp").read_text().replace(" opls_493   SY ", " opls_493   SY2")
    (folder / "ffnonbonded.itp").write_text(nb)
    kept, section = [], None
    for line in (folder / "ffbonded.itp").read_text().splitlines():
        m = re.match(r"\s*\[\s*(\w+)\s*\]", line)
        section = m.group(1) if m else section
        words = line.split()
        if section in ("bondtypes", "angletypes") and "SY" in words:
            continue  # the stock files have no sulfonate bonds and angles
        kept.append(" ".join("SY2" if w == "SY" else w for w in words) if section == "dihedraltypes" else line)
    (folder / "ffbonded.itp").write_text("\n".join(kept) + "\n")
    return folder


# rq-07139b20
def test_sulfonate_parameters(tmp_path):
    register_backend(OplsaaBackend(oplsaa_ff=stock_like_ff(tmp_path)))
    s = build_polymer(homopolymer("DMAPS", 3))
    top = parametrize(s, "oplsaa")
    assert not top.missing_parameters
    cs = [b for b in top.bonds if {s.mol.GetAtomWithIdx(i).GetSymbol() for i in b} == {"C", "S"}]
    assert len(cs) == 3 and all(top.explicit[b] == pytest.approx((0.181, 185769.6)) for b in cs)
    assert ("CT", "SY2") not in top.parameters.bondtypes  # written on the bond lines, not in ffbonded.itp
    line = next(x for x in moleculetype_text(top, "m").splitlines() if x.split()[:2] == [str(cs[0][0] + 1),
                                                                                         str(cs[0][1] + 1)])
    assert line.split()[2:] == ["1", "0.181", "185769.6"]


# rq-b1961917
def test_methyl_caps():
    s = build_polymer(homopolymer("DMAPS", 3))
    heavy = {}
    for a in s.mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            heavy.setdefault(a.GetProp("res"), []).append(a)
    last = heavy["DME"][-1]
    assert last.GetSymbol() == "C" and last.GetTotalNumHs(includeNeighbors=True) == 3
    assert any(n.GetProp("name") == "C2" for n in last.GetNeighbors())  # the substituted backbone carbon
    names = {r: [a.GetProp("name") for a in atoms] for r, atoms in heavy.items()}
    assert names["DMF"] == names["DMR"] == names["DME"][:-1]
    head = build_polymer(homopolymer("DMAPS", 3), initiator="KPS", ends="head")
    assert head.chain.head_group is not None and head.chain.tail_group is None
    assert head.mol.GetNumAtoms() == s.mol.GetNumAtoms() - 1 + 5  # H of the head -> SO4 (5 atoms)


# rq-c612302d
def test_kps_sulfate_gaff2(tmp_path):
    register_backend(OplsaaBackend(oplsaa_ff=stock_like_ff(tmp_path)))
    s = build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    top = parametrize(s, "oplsaa")
    assert not top.missing_parameters
    kph = {a.GetProp("name"): a.GetIdx() for a in s.mol.GetAtoms() if a.GetProp("res") == "KPH"}
    assert {top.atom_types[i] for i in kph.values()} == {"opls_467", "opls_493", "opls_494"}  # stock types
    o1, s1, o2 = kph["O1"], kph["S1"], kph["O2"]
    c = next(n.GetIdx() for n in s.mol.GetAtomWithIdx(o1).GetNeighbors() if n.GetSymbol() == "C")

    def value(*atoms):
        return top.explicit.get(atoms) or top.explicit.get(atoms[::-1])

    assert value(s1, o2) == pytest.approx((0.14657, 514205.232))  # GAFF2 o-s6
    assert value(o1, s1) == pytest.approx((0.16433, 301072.272))  # GAFF2 os-s6
    assert value(o1, s1, o2) == pytest.approx((108.58, 1307.3326))  # GAFF2 o-s6-os
    assert value(c, o1, s1) == pytest.approx((115.73, 782.9938))  # GAFF2 c3-os-s6
    k = 3.6 / 3 * 4.184  # GAFF2 X-os-s6-X: 180 deg, n = 2  ->  RB (2k, 0, -2k, 0, 0, 0)
    assert value(c, o1, s1, o2) == pytest.approx((2 * k, 0, -2 * k, 0, 0, 0), abs=1e-5)
    sulfonate_s = next(a.GetIdx() for a in s.mol.GetAtoms() if a.GetSymbol() == "S" and a.GetProp("res") == "DMR")
    cs = next(n.GetIdx() for n in s.mol.GetAtomWithIdx(sulfonate_s).GetNeighbors() if n.GetSymbol() == "C")
    assert value(cs, sulfonate_s) == pytest.approx((0.181, 185769.6))  # the DMAPS sulfonate keeps its terms


# rq-a67a239e
def test_periodic_to_rb():
    import math

    from zpoly.forcefield._core import periodic_to_rb
    for n in range(1, 6):
        for phase in (0.0, 180.0):
            rb = periodic_to_rb([(phase, 2.5, n), (0.0, 0.7, 1)])
            for deg in range(0, 360, 7):
                phi = math.radians(deg)
                periodic = 2.5 * (1 + math.cos(n * phi - math.radians(phase))) + 0.7 * (1 + math.cos(phi))
                c = math.cos(phi - math.pi)
                assert sum(x * c ** i for i, x in enumerate(rb)) == pytest.approx(periodic, abs=1e-5)


# rq-4503b8e8
def test_spacer_variant(opls):
    from zpoly.zwitterions import make_variant

    s = build_polymer([make_variant("DMAPS", spacer=4)] * 3)
    top = parametrize(s, "oplsaa", charge_scale=0.85)
    assert len(top.builtin_blocks) == 3
    q, block = charges_of(top, top.builtin_blocks[1])
    assert q["N1"] == pytest.approx(0.1337, abs=1.5e-4) and q["S1"] == pytest.approx(1.1919, abs=1.5e-4)
    assert q["C12"] == pytest.approx(-0.3987, abs=1.5e-4) and q["C4"] == pytest.approx(0.2170, abs=1.5e-4)
    assert q["C10"] == pytest.approx(-0.12, abs=1e-3) and q["C11"] == pytest.approx(-0.12, abs=1e-3)  # OPLS CH2
    assert all(sum(b.charges) == pytest.approx(0.0, abs=1e-9) for b in top.blocks)
    alkyl = parametrize(build_polymer([make_variant("DMAPS", alkyl=2)] * 3), "oplsaa")
    assert not alkyl.builtin_blocks
