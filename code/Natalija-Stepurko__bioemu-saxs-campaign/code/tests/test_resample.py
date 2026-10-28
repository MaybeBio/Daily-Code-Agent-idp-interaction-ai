"""analysis/resample/rslib.py: sample inventory, RMSD, Cα radius of gyration and the χ² wiring, on
synthetic data."""
import importlib.util
import sys
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

from bsc import saxs

HERE = Path(__file__).resolve().parents[1] / "analysis" / "resample"
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location("rslib", HERE / "rslib.py")
rslib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rslib)


def _batches(folder: Path, ranges: list[tuple[int, int]], seq: str = "GSAG") -> None:
    for lo, hi in ranges:
        n = hi - lo
        np.savez(folder / f"batch_{lo:07d}_{hi:07d}.npz", pos=np.zeros((n, len(seq), 3), np.float32),
                 node_orientations=np.zeros((n, len(seq), 3, 3), np.float32), sequence=seq)


def test_inventory_counts_samples_across_batches(tmp_path):
    _batches(tmp_path, [(0, 4), (4, 8), (8, 10)])
    (tmp_path / "samples.xtc").write_bytes(b"")
    c = rslib.count_npz_samples(tmp_path)
    assert c == {"n_npz_files": 3, "n_samples_npz": 10, "npz_contiguous": True, "npz_sequences": ["GSAG"]}


def test_inventory_flags_a_gap_between_batches(tmp_path):
    _batches(tmp_path, [(0, 4), (5, 8)])
    c = rslib.count_npz_samples(tmp_path)
    assert c["n_samples_npz"] == 7 and not c["npz_contiguous"]


def test_inventory_orders_batches_as_bioemu_reads_them(tmp_path):
    _batches(tmp_path, [(10, 20), (0, 10), (20, 25)])
    assert [p.name for p in rslib.npz_batches(tmp_path)] == [
        "batch_0000000_0000010.npz", "batch_0000010_0000020.npz", "batch_0000020_0000025.npz"]


def test_rg_matches_mdtraj_and_the_ensemble_value_is_the_rms():
    rng = np.random.default_rng(1)
    xyz = rng.normal(size=(5, 30, 3)).astype(np.float32)
    top = md.Topology()
    ch = top.add_chain()
    for _ in range(30):
        top.add_atom("CA", md.element.carbon, top.add_residue("GLY", ch))
    ref = md.compute_rg(md.Trajectory(xyz, top))
    assert np.allclose(rslib.rg_per_conformer(xyz), ref, rtol=1e-5)
    assert abs(rslib.ensemble_rg(np.array([3.0, 4.0])) - np.sqrt(12.5)) < 1e-12


def test_rg_of_a_known_geometry():
    # two points 2 apart: every atom 1 from the centre
    xyz = np.array([[[0.0, 0, 0], [2.0, 0, 0]]])
    assert np.allclose(rslib.rg_per_conformer(xyz), [1.0])


def _exp(n=60):
    q = np.linspace(0.01, 0.3, n)
    I = np.exp(-(q * 20) ** 2 / 3) + 0.01
    s = 0.02 * I
    return pd.DataFrame({"q": q, "I": I, "sigma": s})


def test_raw_fit_drops_nonfinite_conformers_and_nonpositive_errors():
    exp = _exp()
    exp.loc[[3, 7], "sigma"] = 0.0
    rng = np.random.default_rng(2)
    curves = np.exp(-(exp["q"].to_numpy()[None, :] * rng.uniform(15, 25, (6, 1))) ** 2 / 3)
    curves[2, 5] = np.nan
    r = rslib.raw_fit(curves, exp)
    keep = exp["sigma"].to_numpy() > 0
    avg = curves[[0, 1, 3, 4, 5]][:, keep].mean(axis=0)
    want = saxs.chi2(avg, exp["I"].to_numpy()[keep], exp["sigma"].to_numpy()[keep])
    assert r["n_conformers"] == 5 and r["n_dropped"] == 1
    assert abs(r["chi2"] - want) < 1e-12


def test_raw_fit_fits_scale_and_constant_background():
    exp = _exp()
    model = (exp["I"].to_numpy() - 0.01) / 3.0          # the data are 3 x model + 0.01
    r = rslib.raw_fit(np.vstack([model, model]), exp)
    assert r["chi2"] < 1e-12


def test_raw_fit_rejects_a_curve_on_another_grid():
    exp = _exp()
    try:
        rslib.raw_fit(np.ones((2, len(exp) - 1)), exp)
    except ValueError:
        return
    raise AssertionError("expected a ValueError")


def test_paired_statistics():
    a = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
    b = 1.1 * a
    p = rslib.paired(a, b, seed=0, n_boot=200)
    assert p["n"] == 6 and p["n_a_below_b"] == 6
    assert abs(p["median_ratio_a_over_b"] - 1 / 1.1) < 1e-12
    assert abs(p["median_a_minus_median_b"] - (np.median(a) - np.median(b))) < 1e-12
    assert p["wilcoxon_two_sided"]["p"] < 0.05


def test_pepsi_log_values():
    text = ("Best d_rho found................ : 0.0123 e/A^3\nBest r0 found........ : 1.62 A\n"
            "Chi^2 ......................... : 3.21\n")
    assert rslib.pepsi_log_values(text) == {"d_rho": 0.0123, "r0": 1.62, "chi2_frame": 3.21}


def test_kabsch_rmsd_is_zero_for_a_rotated_copy_and_exact_for_a_shift():
    rng = np.random.default_rng(3)
    a = rng.normal(scale=60.0, size=(2, 2000, 3))          # far-flung atoms, where float32 loses 0.1 Å
    t = np.radians(40.0)
    rot = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
    b = a @ rot.T + 5.0
    assert np.all(rslib.kabsch_rmsd(a, b) < 1e-9)
    # one atom moved by 1: RMSD sqrt(1 / n) before superposition, and no larger after it
    c = a.copy()
    c[:, 0, 0] += 1.0
    r = rslib.kabsch_rmsd(a, c)
    assert np.all(r <= np.sqrt(1 / a.shape[1]) + 1e-12) and np.all(r > 0)


def test_raw_fit_reproduces_the_study_scoring(monkeypatch):
    """The same conformer curves give the same raw chi2 and NRMSD through rslib.raw_fit as through the
    study's own bsc.score.score_entry (path C of the design must reproduce results/scores.csv)."""
    from bsc import score

    exp = _exp(80)
    exp.loc[[10], "sigma"] = -1.0
    rng = np.random.default_rng(4)
    curves = np.exp(-(exp["q"].to_numpy()[None, :] * rng.uniform(15, 25, (20, 1))) ** 2 / 3)
    curves[7, 3] = np.inf
    monkeypatch.setattr(score, "read_prediction", lambda model, label: curves)
    rec, _ = score.score_entry("bioemu", "X", exp)
    r = rslib.raw_fit(curves, exp)
    assert r["n_conformers"] == rec["n_conformers"] == 19 and r["n_dropped"] == rec["n_dropped"] == 1
    assert abs(r["chi2"] - rec["chi2_raw"]) < 1e-12 * rec["chi2_raw"]
    assert abs(r["nrmsd"] - rec["nrmsd_raw"]) < 1e-12
