"""Tests for the one-input-file interface (rqm/input_file.md)."""

import pytest
from test_forcefield_core import FakeBackend
from test_system_slab import solvent  # noqa: F401  (fixture)

from zpoly.cli import main
from zpoly.forcefield import register_backend
from zpoly.inputfile import KEYS, InputError, parse_input


@pytest.fixture
def fake():
    register_backend(FakeBackend("fake"))


# rq-b4454a7b
def test_help_lists_keywords(capsys):
    assert main(["help"]) == 0
    out = capsys.readouterr().out
    for key in ("JOB", "MONOMER", "N", "INITIATOR", "FF", "DFT", "CHAINS", "BOX", "SLAB", "OUTPUT"):
        assert f"\n  {key} " in out
    assert max(len(line) for line in out.splitlines()) <= 80
    listed = [line.split()[0] for line in out.splitlines() if line.startswith("  ") and line.split()[0] in KEYS]
    assert sorted(listed) == sorted(KEYS)  # every keyword exactly once
    assert main([]) == 0


# rq-0af106e4
def test_help_one_key(capsys):
    assert main(["help", "box"]) == 0
    out = capsys.readouterr().out
    assert "Lx Ly Lz" in out and "bulk, slab" in out and "example: BOX = " in out


# rq-73efd1ac
def test_help_unknown_key(capsys):
    assert main(["help", "boxx"]) == 1
    assert "BOX" in capsys.readouterr().err


# rq-1c00331e
def test_list(capsys):
    assert main(["list"]) == 0
    assert any(line.startswith("DMAPS") for line in capsys.readouterr().out.splitlines())
    assert main(["list", "initiators"]) == 0
    assert any(line.startswith("KPS") for line in capsys.readouterr().out.splitlines())


# rq-4ae24f96
def test_template_runs(fake, run_inp, tmp_path, capsys):
    assert main(["template", "build"]) == 0
    text = capsys.readouterr().out.replace("FF           = oplsaa", "FF           = fake")
    assert "FF           = fake" in text
    assert run_inp(text) == 0
    assert (tmp_path / "DMAPS_n12_KPS.pdb").is_file() and (tmp_path / "DMAPS_n12_KPS.top").is_file()
    assert (tmp_path / "DMAPS_n12_KPS.ff").is_dir()
    capsys.readouterr()
    for job in ("dft", "bulk", "slab", "analyze"):
        assert main(["template", job]) == 0
        parse_input(capsys.readouterr().out)


# rq-ff28afbc
def test_build_job(fake, run_inp, tmp_path):
    assert run_inp("MONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = fake\nOUTPUT = d3") == 0
    assert (tmp_path / "d3.pdb").is_file() and (tmp_path / "d3.top").is_file() and (tmp_path / "d3.ff").is_dir()


# rq-c68591e1
def test_pdb_only(run_inp, tmp_path):
    assert run_inp("monomer = DMAPS   # comment\n\nn = 3 ! other comment\nFF = none") == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == ["DMAPS_n3.pdb", "job.inp"]


# rq-244b2f07
def test_copolymer_forms(run_inp, tmp_path):
    for line, n, name in (("alternating DMAPS CBMA", 4, "DMAPS-alt-CBMA_n4"), ("block DMAPS:2 CBMA:2", None,
                          "DMAPS2-CBMA2"), ("random DMAPS:0.5 CBMA:0.5", 4, "DMAPS-ran-CBMA_n4_s0"),
                          ("pattern AAB A=DMAPS B=CBMA", 3, "AAB_n3")):
        text = f"COPOLYMER = {line}\nFF = none" + (f"\nN = {n}" if n else "")
        assert run_inp(text) == 0, line
        assert (tmp_path / f"{name}.pdb").is_file()


# rq-ff7e6eb6
def test_spacer_alkyl(run_inp, tmp_path):
    assert run_inp("MONOMER = DMAPS\nN = 3\nSPACER = 5\nALKYL = 2\nFF = none") == 0
    assert (tmp_path / "DMAPS_n3_sp5_alk2.pdb").is_file()


# rq-eda6a6c5
def test_dft_job(run_inp, tmp_path):
    assert run_inp("JOB = dft\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nDFT = q") == 0
    for name in ("DMAPS.xyz", "KPS_end.xyz", "dft_manifest.yaml"):
        assert (tmp_path / "q" / name).is_file()


# rq-6eca734c
def test_bulk_job(fake, solvent, run_inp, tmp_path):  # noqa: F811
    assert run_inp(f"JOB = bulk\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = fake\nCHAINS = 2\nBOX = 3\n"
                   f"PACKER = builtin\nSOLVENT_FF = {solvent}\nTG = 500 300 50") == 0
    folder = tmp_path / "DMAPS_n3_KPS_bulk"
    assert (folder / "system.top").is_file() and (folder / "run.sh").is_file()
    assert len(list(folder.glob("cool_*K.mdp"))) == 5


# rq-33ff1048
def test_slab_job(fake, solvent, run_inp, tmp_path):  # noqa: F811
    assert run_inp(f"JOB = slab\nMONOMER = DMAPS\nN = 3\nINITIATOR = KPS\nFF = fake\nCHAINS = 4\n"
                   f"BOX = 3.4 3.4 10.2\nSLAB = 3.4\nSALT = KCl 0.1\nTEMPERATURES = 300, 320\n"
                   f"PACKER = builtin\nSOLVENT_FF = {solvent}") == 0
    folder = tmp_path / "DMAPS_n3_KPS_slab"
    assert (folder / "prod_300K.mdp").is_file() and (folder / "prod_320K.mdp").is_file()
    assert "-pname K -nname CL -conc 0.1" in (folder / "run.sh").read_text()


# rq-1c3e7626
def test_analyze_detects_kind(fake, solvent, run_inp, tmp_path, capsys):  # noqa: F811
    assert run_inp(f"JOB = bulk\nMONOMER = DMAPS\nN = 3\nFF = fake\nCHAINS = 2\nBOX = 3\nPACKER = builtin\n"
                   f"SOLVENT_FF = {solvent}\nOUTPUT = b") == 0
    capsys.readouterr()
    assert run_inp("JOB = analyze\nOUTPUT = b", name="a.inp") == 1
    assert "cool_600K" in capsys.readouterr().err
    assert run_inp("JOB = analyze\nOUTPUT = nowhere", name="a.inp") == 1


# rq-390e287d
def test_charges_and_params_jobs(fake, run_inp, tmp_path):
    assert run_inp("JOB = charges\nMONOMER = DMAPS\nN = 3\nFF = fake") == 0
    assert run_inp("JOB = params\nMONOMER = DMAPS\nN = 3\nFF = fake") == 0
    assert (tmp_path / "DMAPS_n3_charges.yaml").is_file() and (tmp_path / "DMAPS_n3_params.yaml").is_file()


# rq-d8db093e
def test_unknown_key(run_inp, capsys):
    assert run_inp("MONMER = DMAPS\nN = 3") == 1
    err = capsys.readouterr().err
    assert "line 1" in err and "MONOMER" in err


# rq-1438d114
def test_repeated_key():
    with pytest.raises(InputError) as exc:
        parse_input("MONOMER = DMAPS\nN = 3\nN = 3")
    assert exc.value.line == 3 and exc.value.key == "N"


# rq-20746e76
def test_bad_value(run_inp, capsys):
    assert run_inp("MONOMER = DMAPS\nN = three") == 1
    assert "N needs an integer" in capsys.readouterr().err
    for text in ("MONOMER = DMAPS\nN = 3\nENDS = left", "JOB = bulk\nMONOMER = DMAPS\nN = 3\nCHAINS = 2\nBOX = 3 3",
                 "garbage line"):
        with pytest.raises(InputError):
            parse_input(text)


# rq-f5c6245f
def test_key_for_other_job():
    with pytest.raises(InputError, match="SLAB is for JOB = slab"):
        parse_input("JOB = bulk\nMONOMER = DMAPS\nN = 3\nCHAINS = 2\nBOX = 3\nSLAB = 3")


# rq-3d1ea632
def test_missing_required():
    with pytest.raises(InputError) as exc:
        parse_input("JOB = bulk\nMONOMER = DMAPS\nN = 3")
    assert "CHAINS" in str(exc.value) and "BOX" in str(exc.value)
    for text in ("MONOMER = DMAPS\nSEQUENCE = DMAPS CBMA\nN = 2", "N = 3", "MONOMER = DMAPS",
                 "COPOLYMER = block DMAPS:2 CBMA:2\nN = 4", "MONOMER = DMAPS\nN = 3\nFF = none\nDFT = dft"):
        with pytest.raises(InputError):
            parse_input(text)


# rq-7474fcd4
def test_errors_reported(run_inp, tmp_path, capsys):
    assert run_inp("MONOMER = NOPE\nN = 3\nFF = none") == 1
    err = capsys.readouterr().err
    assert err.startswith("error:") and "NOPE" in err
    assert sorted(p.name for p in tmp_path.iterdir()) == ["job.inp"]
    assert main(["missing.inp"]) == 1
    assert main(["a.inp", "b.inp"]) == 2


def test_defaults_filled():
    opts = parse_input("MONOMER = DMAPS\nN = 3")
    assert opts["job"] == "build" and opts["ff"] == "gaff2" and opts["seed"] == 0 and opts["conformation"] == "coil"
    assert "chains" not in opts


def test_long_output_name(fake, run_inp, tmp_path):
    name = "DMAPS_sp4_n50_KPS_long_name"
    assert run_inp(f"MONOMER = DMAPS\nN = 2\nFF = fake\nOUTPUT = {name}") == 0
    lines = (tmp_path / f"{name}.top").read_text().splitlines()
    assert lines[lines.index("[ moleculetype ]") + 2].split() == [name, "3"]
    assert lines[-1].split() == [name, "1"]
