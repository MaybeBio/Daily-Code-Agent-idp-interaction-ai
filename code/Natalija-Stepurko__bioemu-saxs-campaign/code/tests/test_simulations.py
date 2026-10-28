"""analysis/simulations/simlib.py: radius of gyration from coordinates, the Guinier comparison, the
convergence check, the buffer rule, the input table and the compute plan."""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from bsc import saxs

_dir = Path(__file__).resolve().parents[1] / "analysis" / "simulations"
sys.path.insert(0, str(_dir))
spec = importlib.util.spec_from_file_location("simlib", _dir / "simlib.py")
simlib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(simlib)


def test_rg_of_two_beads_and_mass_weighting():
    xyz = np.array([[[0.0, 0, 0], [2.0, 0, 0]]])
    assert abs(simlib.rg_frames(xyz, np.array([1.0, 1.0]))[0] - 1.0) < 1e-12
    # masses 1 and 3: centre at x = 1.5, Rg^2 = (1 * 1.5^2 + 3 * 0.5^2) / 4 = 0.75
    assert abs(simlib.rg_frames(xyz, np.array([1.0, 3.0]))[0] - math.sqrt(0.75)) < 1e-12


def test_rg_matches_mdtraj_for_equal_masses():
    import mdtraj as md

    rng = np.random.default_rng(1)
    xyz = rng.normal(size=(4, 30, 3)).astype(np.float32)
    top = md.Topology()
    ch = top.add_chain()
    for _ in range(30):
        top.add_atom("CA", md.element.carbon, top.add_residue("ALA", ch))
    ref = md.compute_rg(md.Trajectory(xyz, top))
    assert np.allclose(simlib.rg_frames(xyz.astype(float), np.ones(30)), ref, atol=1e-5)


def test_ensemble_rg_is_the_root_mean_square():
    assert abs(simlib.ensemble_rg(np.array([3.0, 4.0])) - math.sqrt(12.5)) < 1e-12


def test_residue_masses_add_the_termini():
    m = simlib.residue_masses("GAG", {"G": 57.05, "A": 71.07})
    assert np.allclose(m, [59.05, 71.07, 73.05])


def test_unwrap_joins_a_chain_split_by_the_box():
    box = np.array([10.0, 10.0, 10.0])
    whole = np.array([[[4.0, 5, 5], [8.0, 5, 5], [12.0, 5, 5], [16.0, 5, 5]]])
    wrapped = whole % box
    out = simlib.unwrap_chain(wrapped, box)
    assert np.allclose(np.diff(out, axis=1), np.diff(whole, axis=1))
    assert abs(simlib.rg_frames(out, np.ones(4))[0] - simlib.rg_frames(whole, np.ones(4))[0]) < 1e-12


def _random_coil(n, rng, b=3.8):
    steps = rng.normal(size=(n - 1, 3))
    steps = b * steps / np.linalg.norm(steps, axis=1, keepdims=True)
    return np.vstack([np.zeros(3), np.cumsum(steps, axis=0)])


def _debye(confs, q):
    out = []
    for x in confs:
        d = np.linalg.norm(x[:, None] - x[None], axis=2)
        qd = q[:, None, None] * d[None]
        out.append(np.where(qd > 0, np.sin(qd) / np.where(qd > 0, qd, 1), 1.0).sum(axis=(1, 2)))
    return saxs.ensemble_curve(np.array(out))


def test_guinier_rg_of_a_computed_curve_matches_the_coordinate_rg():
    """The Guinier comparison: for compact particles (points filling a sphere) a Debye curve fitted as the
    archive's curves are returns the coordinate radius of the ensemble (root mean square over conformers)."""
    rng = np.random.default_rng(7)
    confs = []
    for _ in range(10):
        p = rng.uniform(-1, 1, size=(400, 3))
        confs.append(15.0 * p[np.linalg.norm(p, axis=1) <= 1][:150])
    q = np.linspace(0.005, 0.4, 300)
    g = saxs.guinier(q, _debye(confs, q), 0.01 * _debye(confs, q) + 1e-9, 1.3)
    rg_coord = simlib.ensemble_rg(simlib.rg_frames(np.array(confs), np.ones(150)))
    assert g["valid"] and abs(g["rg"] / rg_coord - 1) < 0.02


def test_guinier_rg_of_a_chain_is_below_its_coordinate_rg():
    """For chain-like ensembles the same fit (q Rg up to 1.3) returns a radius a few per cent below the
    coordinate radius, the Debye-function curvature; simlib.guinier_bias_gaussian_chain gives it
    analytically."""
    rng = np.random.default_rng(7)
    confs = [_random_coil(60, rng) for _ in range(20)]
    q = np.linspace(0.005, 0.3, 300)
    raw = _debye(confs, q)
    g = saxs.guinier(q, raw, 0.01 * raw + 1e-9, 1.3)
    ratio = g["rg"] / simlib.ensemble_rg(simlib.rg_frames(np.array(confs), np.ones(60)))
    assert 0.85 < ratio < 0.97
    b = simlib.guinier_bias_gaussian_chain(saxs.guinier)
    assert 0.85 < b < 0.97 and abs(b - ratio) < 0.04


def test_convergence_flags_a_drift():
    rng = np.random.default_rng(0)
    flat = 30 + rng.normal(0, 3, 1000)
    drift = np.linspace(20, 40, 1000) + rng.normal(0, 1, 1000)
    assert simlib.convergence(flat)["converged"]
    c = simlib.convergence(drift)
    assert not c["converged"] and c["halves_rel_diff"] > 0.05 and len(c["block_rg"]) == 5


def test_closer_to_one_on_the_log_scale():
    assert simlib.closer_to_one(np.array([0.9, 1.3]), np.array([1.2, 0.8])).tolist() == [True, False]


def test_buffer_rule():
    sp = simlib.parse_buffer("50 mM HEPES, 100 mM NaCl, 0.5 mM tris(2-carboxyethyl)phosphine (TCEP),")
    assert set(sp) == {"PBS", "HEPES", "NaCl"} and abs(sp["NaCl"] - 0.1) < 1e-12  # TCEP is not Tris
    sp = simlib.parse_buffer("150 mM NaCl, 50 mM BisTris pH 6.8, and 0.5 mM EDTA")
    assert "Tris" not in sp and abs(sp["Bis-Tris"] - 0.05) < 1e-12
    sp = simlib.parse_buffer("25 mM Hepes, 150 NaCl, 1 mM DTT, 5% Glycerol")  # no unit -> mM
    assert abs(sp["NaCl"] - 0.15) < 1e-12
    sp = simlib.parse_buffer("PBS, 10 mM DTT, 6 M urea, 0.3 M KCl")
    assert sp["PBS"] and sp["urea"] == 6.0 and abs(sp["KCl"] - 0.3) < 1e-12
    assert abs(simlib.ionic_strength(sp, 7.4) - 0.47) < 1e-12  # urea adds nothing
    sp = simlib.parse_buffer("10 mM Tris HCl, 150 mM NaCl, 10 mM CaCl₂")
    assert abs(simlib.ionic_strength(sp, 8.1) - (0.005 + 0.15 + 0.03)) < 1e-12  # Tris half protonated
    # phosphate at its pKa: half HPO4(2-), I = c (1 + 2 * 0.5)
    assert abs(simlib.ionic_strength({"PBS": False, "phosphate": 0.02}, 7.2) - 0.04) < 1e-12


def test_conditions_defaults_and_sources():
    c = simlib.conditions("Sodium phosphate buffer", None, float("nan"), float("nan"))
    assert (c["temp_k"], c["ph"], c["ionic_m"]) == (293.15, 7.0, 0.15)
    assert (c["temp_source"], c["ph_source"], c["ionic_source"]) == ("default", "default", "default")
    c = simlib.conditions("20 mM HEPES 150 mM NaCl", "150 mM NaCl", 10.0, 7.5)
    assert c["temp_k"] == 283.15 and c["ionic_source"] == "parsed"
    assert abs(c["ionic_m"] - round(0.15 + 0.01, 4)) < 1e-9  # NaCl counted once
    c = simlib.conditions("20 mM HEPES, 150 mM NaCl, 1 M urea, 5 mM DTT", None, 20.0, 7.0)
    assert c["denaturant"] and c["urea_m"] == 1.0


def test_run_settings_follow_the_idrome_protocol():
    s = simlib.run_settings(51)
    assert s["n_save"] == 7000 and s["n_frames"] == 1010 and abs(s["box_nm"] - 23.0) < 1e-9
    s = simlib.run_settings(475)
    assert s["n_save"] == 68000 and abs(s["box_nm"] - (474 * 0.38 + 4)) < 1e-9


def test_prepare_inputs_reads_the_right_sequence_and_conditions():
    resolv = pd.DataFrame(
        {
            "model": ["bioemu", "bioemu", "bioemu", "peptron"],
            "label": ["SASD2", "SASD1", "SASD3", "SASD4"],
            "clean": [True, True, False, True],
            "disorder_class": ["disordered", "disordered", "disordered", "disordered"],
        }
    )
    entries = pd.DataFrame(
        {
            "label": ["SASD1", "SASD2", "SASD3", "SASD4"],
            "sequence": ["MKKE", "GSHMAA", "AAAA", "PPPP"],
            "rg_exp": [10.0, 12.0, 13.0, 14.0],
        }
    )
    meta = pd.DataFrame(
        {
            "label": ["SASD1", "SASD2", "SASD3", "SASD4"],
            "buffer": ["20 mM Tris, 50 mM NaCl", "PBS", "x", "x"],
            "buffer_additive": [None] * 4,
            "temp": [4.0, float("nan"), 20.0, 20.0],
            "ph": [8.1, 7.4, 7.0, 7.0],
        }
    )
    t = simlib.prepare_inputs(resolv, entries, meta)
    assert t["label"].tolist() == ["SASD1", "SASD2"]  # clean BioEmu-1 rows only, sorted
    assert t["sequence"].tolist() == ["MKKE", "GSHMAA"] and t["length"].tolist() == [4, 6]
    assert t["rg_exp_A"].tolist() == [10.0, 12.0]
    assert (
        t.loc[0, "temp_k"] == 277.15 and t.loc[1, "temp_k"] == 293.15 and t.loc[1, "temp_source"] == "default"
    )
    assert abs(t.loc[0, "ionic_m"] - 0.06) < 1e-9 and t.loc[1, "ionic_m"] == 0.17
    assert t["seed"].tolist() == [simlib.SEED_BASE, simlib.SEED_BASE + 1]


def test_speed_fit_recovers_a_power_law():
    n = np.array([50.0, 180.0, 475.0])
    a, b = simlib.fit_speed(n, 3000.0 * n**-0.8)
    assert abs(a - 3000.0) < 1e-6 and abs(b - 0.8) < 1e-9


def test_makespan_longest_first():
    assert simlib.makespan([5, 4, 3, 3, 3], 2) == 10.0  # 5+3+... on one worker, 4+3+3 on the other
    assert simlib.makespan([10, 1, 1], 6) == 10.0


def test_compute_plan_reduces_long_runs_then_cuts_by_length():
    t = pd.DataFrame(
        {
            "label": list("abcd"),
            "length": [50, 100, 300, 500],
            "n_save": [1000] * 4,
            "n_frames": [1000] * 4,
            "n_discard": [10] * 4,
            "steps_done": [0] * 4,
        }
    )
    # 1000 / 3.6 steps per second at every length: 1000 frames x 1000 steps take 1 hour
    p = simlib.compute_plan(t, a=1000 / 3.6, b=0.0, workers=2, budget_h=24)
    assert p["included"].all() and (p["frames_planned"] == 1000).all()
    t["n_save"] = [1000, 1000, 20000, 60000]  # 20 h and 60 h at full length
    p = simlib.compute_plan(t, a=1000 / 3.6, b=0.0, workers=2, budget_h=24)
    assert (
        p.set_index("label").loc["d", "frames_planned"] == 505
        and p.set_index("label").loc["c", "frames_planned"] == 505
    )
    # d at 505 frames is 30.3 h > 24 h: left out; c at 10.1 h fits
    assert p["included"].tolist() == [True, True, True, False]
