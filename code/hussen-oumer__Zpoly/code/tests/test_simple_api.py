"""Tests for the simple Python API (rqm/simple_api.md)."""

import pytest
import yaml
from test_dft_charges import pattern_charges, write_plain
from test_forcefield_core import FakeBackend

import zpoly as pp
from zpoly.forcefield import DftChargeError, register_backend


@pytest.fixture
def fake():
    register_backend(FakeBackend("fake"))


def atoms_of(folder):
    return {m["name"]: m for m in yaml.safe_load((folder / "dft_manifest.yaml").read_text())["molecules"]}


def residues(pdb):
    seen = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM"):
            seen.setdefault(int(line[22:26]), line[17:20])
    return [seen[k] for k in sorted(seen)]


# rq-a3b3b51e
def test_top_level_imports():
    for name in ("build", "dft_inputs", "alternating", "block", "random_copolymer", "pattern"):
        assert callable(getattr(pp, name))


# rq-621a8723
def test_dft_inputs(tmp_path):
    paths = pp.dft_inputs("DMAPS", initiator="KPS", folder=tmp_path / "dft")
    assert [p.name for p in paths] == ["DMAPS.xyz", "KPS_end.xyz"]
    assert all(p.exists() for p in paths)


# rq-bd003416
def test_build_with_dft_folder(fake, tmp_path):
    folder = tmp_path / "dft"
    pp.dft_inputs("DMAPS", initiator="KPS", folder=folder)
    m = atoms_of(folder)
    for name in ("DMAPS", "KPS_end"):
        write_plain(folder, name, m[name]["atoms"], pattern_charges(m[name]["atoms"], m[name]["net_charge"]))
    result = pp.build("DMAPS", n=4, initiator="KPS", dft=folder, ff="fake", out=tmp_path / "o" / "d4")
    assert {p.name for p in (tmp_path / "o").iterdir()} == {"d4.pdb", "d4.top", "d4.ff"}
    assert set(result.dft_blocks) == {"KPH", "DMF", "DMR", "DME", "KPT"} and result.net_charge == -2


# rq-a7d55590
def test_build_with_pasted_values(fake, tmp_path):
    folder = tmp_path / "ref"
    pp.dft_inputs("DMAPS", initiator="KPS", folder=folder)
    m = atoms_of(folder)
    values = {"DMAPS": pattern_charges(m["DMAPS"]["atoms"], 0), "KPS": pattern_charges(m["KPS_end"]["atoms"], -1)}
    result = pp.build("DMAPS", n=3, initiator="KPS", dft=values, ff="fake", out=tmp_path / "d3")
    assert {"DMR", "KPH"} <= set(result.dft_blocks)


# rq-1fe6b25c
def test_pasted_wrong_count(fake, tmp_path):
    with pytest.raises(DftChargeError):
        pp.build("DMAPS", n=3, dft={"DMAPS": [0.0, 0.0]}, ff="fake", out=tmp_path / "x")
    assert not list(tmp_path.iterdir())


# rq-da585695
def test_pasted_unknown_monomer(fake, tmp_path):
    with pytest.raises(DftChargeError, match="CBMA"):
        pp.build("DMAPS", n=3, dft={"CBMA": [0.0] * 40}, ff="fake", out=tmp_path / "x")


# rq-78cb44ef
def test_copolymer_from_helper(fake, tmp_path):
    result = pp.build(pp.alternating("DMAPS", "CBMA", 4), ff="fake", out=tmp_path / "alt")
    assert residues(result.pdb) == ["DMF", "CBR", "DMR", "CBE"]


# rq-7896b345
def test_pdb_only(tmp_path):
    result = pp.build("DMAPS", n=3, ff="none", out=tmp_path / "p")
    assert [p.name for p in tmp_path.iterdir()] == ["p.pdb"] and result.top is None


# rq-ca771be1
def test_length_mismatch(tmp_path):
    with pytest.raises(ValueError):
        pp.build(["DMAPS", "CBMA"], n=5, ff="none", out=tmp_path / "m")
