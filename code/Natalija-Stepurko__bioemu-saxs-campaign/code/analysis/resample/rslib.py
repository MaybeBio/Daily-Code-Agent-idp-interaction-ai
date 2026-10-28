"""Pure helpers of the v1.2 re-sample analysis: sample counting, RMSD, Cα radius of gyration, the raw-fit
scoring of a set of conformer curves and the paired statistics. Tested in tests/test_resample.py.

Also imported by 02_convert.py in the BioEmu environment, which has numpy, pandas and scipy but not this
package; the bsc modules are therefore imported where they are used."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BATCH = re.compile(r"batch_(\d+)_(\d+)\.npz$")


def npz_batches(folder: Path) -> list[Path]:
    """The sample batches in BioEmu's order (sorted file names, as bioemu.sample.main reads them)."""
    return sorted(p for p in folder.glob("batch_*.npz") if BATCH.search(p.name))


def count_npz_samples(folder: Path) -> dict:
    """Samples in the npz batches of one protein: files, samples (first axis of `pos`), whether the index
    ranges in the file names tile 0..n without gap or overlap, and the distinct sequences stored."""
    files = npz_batches(folder)
    n, seqs, expect, contiguous = 0, set(), 0, True
    for f in files:
        lo, hi = (int(x) for x in BATCH.search(f.name).groups())
        with np.load(f) as d:
            k = int(d["pos"].shape[0])
            seqs.add(str(d["sequence"]))
        contiguous &= (lo == expect) and (hi - lo == k)
        expect = hi
        n += k
    return {"n_npz_files": len(files), "n_samples_npz": n, "npz_contiguous": bool(contiguous),
            "npz_sequences": sorted(seqs)}


def rg_per_conformer(xyz: np.ndarray) -> np.ndarray:
    """Radius of gyration of each conformer with equal atom weights, xyz shape (n_conformers, n_atoms, 3);
    the same quantity as mdtraj.compute_rg without masses."""
    c = xyz - xyz.mean(axis=1, keepdims=True)
    return np.sqrt((c**2).sum(axis=2).mean(axis=1))


def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """RMSD of each pair of conformers a[i], b[i] after optimal superposition, in float64; shape
    (n, atoms, 3). mdtraj.rmsd works in float32, which for a few thousand atoms far from the centre leaves
    an error of 0.1 Å between identical conformers, larger than the xtc precision this is compared with."""
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    out = np.empty(len(a))
    for i in range(len(a)):
        u, _, vt = np.linalg.svd(a[i].T @ b[i])
        d = np.sign(np.linalg.det(u @ vt))
        r = u @ np.diag([1.0, 1.0, d]) @ vt
        out[i] = np.sqrt(np.mean(np.sum((a[i] @ r - b[i]) ** 2, axis=1)))
    return out


def ensemble_rg(rg: np.ndarray) -> float:
    """The ensemble value used by the coordrg stage: sqrt(mean Rg^2) over conformers."""
    rg = np.asarray(rg, float)
    return float(np.sqrt(np.mean(rg**2))) if len(rg) else float("nan")


def raw_fit(curves: np.ndarray, exp: pd.DataFrame) -> dict:
    """The raw fit as bsc.score.score_entry computes it: conformer curves with a non-finite value are dropped,
    points with sigma <= 0 are left out, the uniform average is scaled to the data with one factor and one
    constant background. curves: (n_conformers, n_q) on the measured q grid; exp: columns q, I, sigma."""
    from bsc import saxs
    q, I, s = exp["q"].to_numpy(), exp["I"].to_numpy(), exp["sigma"].to_numpy()
    if curves.ndim != 2 or curves.shape[1] != len(q):
        raise ValueError(f"curves of shape {curves.shape} for {len(q)} measured points")
    s = np.where(s > 0, s, np.nan)
    keep = np.isfinite(s)
    ok = np.isfinite(curves).all(axis=1)
    out = {"n_conformers": int(ok.sum()), "n_dropped": int((~ok).sum())}
    if ok.sum() == 0:
        return {**out, "chi2": float("nan"), "nrmsd": float("nan")}
    avg = saxs.ensemble_curve(curves[ok][:, keep])
    return {**out, "chi2": saxs.chi2(avg, I[keep], s[keep]), "nrmsd": saxs.nrmsd_log(avg, I[keep], s[keep])}


def median_ci(x: np.ndarray, seed: int, n_boot: int) -> dict:
    from bsc import robust
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    lo, hi = robust.bootstrap_ci(x, np.median, n_boot=n_boot, seed=seed)
    return {"n": int(len(x)), "median": float(np.median(x)) if len(x) else float("nan"), "ci95": [lo, hi]}


def paired(a: np.ndarray, b: np.ndarray, seed: int, n_boot: int) -> dict:
    """Paired comparison of positive quantities a (v1.2) and b (archive) on the same proteins: medians of
    each with bootstrap intervals; d = log a - log b with its median (and exp of it, the median ratio) and a
    bootstrap interval; the difference of the medians with a paired bootstrap interval (proteins resampled
    together); a two-sided Wilcoxon signed-rank test on d; counts of a < b."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    a, b = a[ok], b[ok]
    d = np.log(a) - np.log(b)
    out = {"n": int(len(a)), "a": median_ci(a, seed, n_boot), "b": median_ci(b, seed, n_boot),
           "log_diff": median_ci(d, seed, n_boot), "n_a_below_b": int((a < b).sum())}
    out["median_ratio_a_over_b"] = float(np.exp(out["log_diff"]["median"])) if len(d) else float("nan")
    out["median_ratio_ci95"] = [float(np.exp(v)) for v in out["log_diff"]["ci95"]]
    if len(a) >= 2:
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, len(a), size=(n_boot, len(a)))
        dm = np.median(a[idx], axis=1) - np.median(b[idx], axis=1)
        out["median_a_minus_median_b"] = float(np.median(a) - np.median(b))
        out["median_a_minus_median_b_ci95"] = [float(np.quantile(dm, 0.025)), float(np.quantile(dm, 0.975))]
    if len(d) >= 1 and np.any(d != 0):
        w = stats.wilcoxon(d, alternative="two-sided")
        out["wilcoxon_two_sided"] = {"statistic": float(w.statistic), "p": float(w.pvalue)}
    return out


def pepsi_log_values(text: str) -> dict:
    """d_rho, r0 and chi2 from a Pepsi-SAXS 3.0 log (as analysis/simulations/06_saxs.py reads them)."""
    rec = {}
    for line in text.splitlines():
        if line.startswith("Best d_rho found"):
            rec["d_rho"] = float(line.split()[-2])
        elif line.startswith("Best r0 found"):
            rec["r0"] = float(line.split()[-2])
        elif line.startswith("Chi^2"):
            rec["chi2_frame"] = float(line.split()[-1])
    return rec
