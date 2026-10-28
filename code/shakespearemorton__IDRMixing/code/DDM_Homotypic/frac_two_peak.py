"""
Fraction of frames whose per-frame voxel concentrations resolve two peaks,
for Figure 1E. Only needed at the voxel size Figure1.py plots.

    python frac_two_peak.py <index>    writes peaks_<S>/<name>.npz

Same definition as Revisions/analyze_log_v1.py, which reproduces the n_peaks
in the original results_20_log npz files frame for frame: each frame is
slab-centred as in DDM.py, every voxel (empty ones included) is converted to
mM, histogrammed on 1024 points from 0 to 1.05 x the maximum, smoothed with
a Gaussian filter, and its local maxima counted.

Two versions are written from the same pass over the trajectory:
  frac_two_peak      all voxels, as in the original analysis
  frac_two_peak_occ  empty voxels dropped, so the spike at 0 mM cannot be
                     counted as a phase
run_ddm.py --collect merges the result into slab_results_<S>.csv.
"""

import argparse
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

import run_ddm
from run_ddm import DDM

HERE = Path(__file__).resolve().parent

S = 20                  # voxel side length plotted in Figure 1E
PEAK_CHECK_GRID = 1024  # grid for the per-frame peak count


# --- per-frame peak count ------------------------------------------------

def smoothed_density(v, n_grid=PEAK_CHECK_GRID, lo=None, hi=None):
    lo = v.min() if lo is None else lo
    hi = v.max() if hi is None else hi
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return None, None
    hist, edges = np.histogram(v, bins=n_grid, range=(lo, hi))
    bin_width = edges[1] - edges[0]
    if hist.sum() == 0:
        return None, None
    density = hist / (hist.sum() * bin_width)
    std = v.std(ddof=1) if v.size > 1 else 0.0
    if std <= 0:
        return None, None
    sigma_bins = max(std * v.size ** (-1 / 5) / bin_width, 1.0)
    return (edges[:-1] + 0.5 * bin_width,
            gaussian_filter1d(density, sigma=sigma_bins, mode='nearest'))


def frame_peak_count(s):
    """Number of phases resolved in one frame's voxel concentrations."""
    smax = s.max()
    if smax <= 0:
        return 0
    x, y = smoothed_density(s, lo=0.0, hi=smax * 1.05)
    if x is None:
        return 1 if s.size else 0
    padded = np.concatenate(([y[0] - 1.0], y, [y[-1] - 1.0]))
    pk, _ = find_peaks(padded, height=1e-8)
    pk = pk - 1
    return int(((pk >= 0) & (pk < y.size)).sum())


def peak_counts(u, start, chunk_size=S):
    """Peak count of every frame from `start`, on the DDM.py voxel grid,
    with all voxels and with the empty ones dropped."""
    box = u.dimensions[:3]
    n_chunks = np.maximum(np.round(box / chunk_size).astype(int), 1)
    chunk_dims = box / n_chunks
    ny, nz = n_chunks[1], n_chunks[2]
    n_voxels = int(np.prod(n_chunks))
    conv = DDM.ATOMS_A3_TO_mM / np.prod(chunk_dims) / len(u.segments[0].atoms)

    ag = u.atoms
    lz = u.dimensions[2]
    edges = np.arange(0, lz + 1, 1)
    z = edges[:-1] + 0.5

    n_frames = len(u.trajectory) - start
    n_peaks = np.empty(n_frames, dtype=np.int16)
    n_peaks_occ = np.empty(n_frames, dtype=np.int16)
    for i in range(n_frames):
        DDM.centre_slab(ag, u.trajectory[i + start], z, edges, lz)
        idx = np.clip((ag.positions // chunk_dims).astype(int), 0, n_chunks - 1)
        flat_idx = idx[:, 0] * ny * nz + idx[:, 1] * nz + idx[:, 2]
        counts = np.bincount(flat_idx, minlength=n_voxels)
        n_peaks[i] = frame_peak_count(counts * conv)
        n_peaks_occ[i] = frame_peak_count(counts[counts > 0] * conv)
    return n_peaks, n_peaks_occ


def two_peak_fraction(n_peaks):
    return float((n_peaks >= 2).mean()) if n_peaks.size else np.nan


def run(index):
    s = run_ddm.sample(index)
    if s is None:
        return
    name, pdb, dcd, start = s
    out = HERE / f"peaks_{S}" / f"{name}.npz"
    if out.exists():
        print(f"{name}, S={S}: already done, skipping")
        return

    n_peaks, n_peaks_occ = peak_counts(mda.Universe(str(pdb), str(dcd)), start)
    frac = two_peak_fraction(n_peaks)
    frac_occ = two_peak_fraction(n_peaks_occ)
    print(f"{name}: two peaks in {frac:.1%} of {n_peaks.size} frames "
          f"({frac_occ:.1%} with empty voxels dropped)")

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, frac_two_peak=frac, frac_two_peak_occ=frac_occ,
                        n_peaks=n_peaks, n_peaks_occ=n_peaks_occ,
                        chunk_size=S, start_frame=start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Per-frame two-peak fraction for one homotypic slab.")
    parser.add_argument("index", type=int,
                        help="Index of the sequence directory to analyse.")
    run(parser.parse_args().index)
