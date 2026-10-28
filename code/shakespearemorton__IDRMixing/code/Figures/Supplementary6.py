#!/usr/bin/env python3
"""Supplementary figure: equilibration of 50:50 mixtures, tracked by the
degree of mixing.

Every simulation starts completely demixed, so the degree of mixing (d.o.m.)
of Chen and Jacobs (Eq. S25 of their SI) starts near zero in all systems.
Where a pair mixes it climbs and plateaus; where a pair demixes it stays
low.  A low value is therefore not by itself evidence of a system that has
not equilibrated, which is why each trace is judged against its own value at
the end of the run rather than against an absolute threshold.

    d.o.m. = 1 - (1/2) * sum_[dense bins] (p2 - p1)^2 / (p1 + p2)

with p1 and p2 the normalised probabilities of finding a bead of each
species at z, restricted to the dense phase and normalised over it.  Because
the two species were placed in separate slabs at t = 0, "initialised in
phase alpha" and "is species 1" are the same label here.

This reads whole trajectories frame by frame and is slow, so:

  * one .npz per system is written as soon as that system finishes, and
    systems already on disk are skipped -- an interrupted run resumes where
    it stopped;
  * only the trace itself is kept (two float32 arrays per system, a few kB),
    not the per-frame profiles;
  * a job index can be passed for a batch array, exactly as in
    DDM_Homotypic/run_ddm.py.

Usage
-----
    python Supplementary6.py                 # all systems, then plot
    python Supplementary6.py 17              # one system, for a batch array
    python Supplementary6.py --n-systems     # print the job count and exit
    python Supplementary6.py --plot-only     # re-plot from the cached traces
    python Supplementary6.py --force         # recompute every trace
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from DDM import centre_slab, DEFAULT_START

from Supplementary2 import system_rows, split_species, TRAJ_DIR

plt.rcParams.update({
    'font.size': 18, 'axes.linewidth': 2, 'axes.spines.top': False,
    'axes.spines.right': False, 'xtick.major.size': 6, 'xtick.major.width': 2,
    'ytick.major.size': 6, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300,
})

# ---------------------------------------------------------------- paths
TRACE_DIR = Path("../DDM_Heterotypic/dom_traces")
SUMMARY_CSV = Path("../DATASETS/dom_equilibration.csv")
OUT_PDF = Path("Supplementary_Figure6.pdf")

# ---------------------------------------------------------------- settings
LAST_FRAME = 4000       # frames beyond this are not read
STRIDE = 10             # analyse every Nth frame
N_Z_BINS = 100
DENSE_FRAC = 0.5        # a bin is dense above this fraction of the frame's peak
PLATEAU_FRAMES = 500    # window at the end of the run defining the final value
TOL = 0.05              # |d.o.m. - final| within this counts as equilibrated
EQUIL_LINE = DEFAULT_START      # the equilibration cut used everywhere else


# ======================================================================
# per-system trace
# ======================================================================

def trace_path(base_name):
    return TRACE_DIR / f"{base_name}_dom.npz"


def degree_of_mixing(h1, h2, dense):
    """Eq. S25 over the dense bins of one frame.

    Zero when the two species occupy disjoint regions, one when both are
    everywhere present in their global ratio.
    """
    if not dense.any():
        return np.nan
    p1, p2 = h1[dense].astype(float), h2[dense].astype(float)
    s1, s2 = p1.sum(), p2.sum()
    if s1 <= 0 or s2 <= 0:
        return np.nan
    p1, p2 = p1 / s1, p2 / s2          # each normalised over the dense phase
    tot = p1 + p2
    ok = tot > 0
    return float(1.0 - 0.5 * np.sum((p2[ok] - p1[ok]) ** 2 / tot[ok]))


def compute_trace(row, last=LAST_FRAME, stride=STRIDE):
    """d.o.m. against frame number for one system."""
    import MDAnalysis as mda

    base_name = row['base_name']
    pdb = TRAJ_DIR / f"{base_name}.pdb"
    dcds = [TRAJ_DIR / f"{base_name}.dcd"]
    dcds += sorted(TRAJ_DIR.glob(f"{base_name}_resume*.dcd"))
    dcds += sorted(TRAJ_DIR.glob(f"{base_name}_[0-9]*.dcd"))
    dcds = [str(d) for d in dcds if Path(d).is_file()]
    if not pdb.is_file() or not dcds:
        raise FileNotFoundError(f"no trajectory for {base_name}")

    u = mda.Universe(str(pdb), dcds)
    n_total = len(u.trajectory)
    ag1, ag2 = split_species(u, int(row['len1']), int(row['len2']))
    ag = u.atoms

    lz = u.dimensions[2]
    cedges = np.arange(0, lz + 1, 1)
    cz = cedges[:-1] + 0.5
    pedges = np.linspace(0, lz, N_Z_BINS + 1)

    frames = np.arange(0, min(last, n_total), stride)
    dom = np.full(frames.size, np.nan)

    for i, f in enumerate(tqdm(frames, desc=base_name[:38], leave=False,
                               unit='frame')):
        centre_slab(ag, u.trajectory[int(f)], cz, cedges, lz)
        h1, _ = np.histogram(ag1.positions[:, 2], bins=pedges)
        h2, _ = np.histogram(ag2.positions[:, 2], bins=pedges)
        tot = h1 + h2
        # the condensate, as the bins carrying most of the material in this
        # frame; a fixed mask cannot be used because the slab is still
        # rearranging early on
        dense = tot > DENSE_FRAC * tot.max() if tot.max() > 0 else np.zeros_like(tot, bool)
        dom[i] = degree_of_mixing(h1, h2, dense)

    return frames.astype(np.int32), dom.astype(np.float32), n_total


def ensure_trace(row, force=False):
    """Compute and store one system's trace unless it is already on disk."""
    base_name = row['base_name']
    path = trace_path(base_name)
    if path.exists() and not force:
        return 'cached'

    frames, dom, n_total = compute_trace(row)
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    # written immediately, so an interrupted run resumes at the next system
    np.savez_compressed(path, frames=frames, dom=dom,
                        n_frames_total=n_total,
                        stride=STRIDE, last_frame=LAST_FRAME)
    return f"{frames.size} points"


# ======================================================================
# summary
# ======================================================================

def equilibration_frame(frames, dom, tol=TOL, plateau=PLATEAU_FRAMES):
    """First frame after which the trace stays within tol of its final value.

    The final value is the mean over the last `plateau` frames, so a system
    that never settles simply reports the end of the run.
    """
    ok = np.isfinite(dom)
    if ok.sum() < 3:
        return np.nan, np.nan
    f, d = frames[ok], dom[ok]
    final = float(d[f >= f.max() - plateau].mean())

    within = np.abs(d - final) <= tol
    # walk back from the end to the first point that breaks the tolerance
    idx = len(within)
    while idx > 0 and within[idx - 1]:
        idx -= 1
    return (float(f[idx]) if idx < len(f) else float(f[-1])), final


def build_summary(df, out_csv=SUMMARY_CSV):
    rows = []
    for _, row in df.iterrows():
        path = trace_path(row['base_name'])
        if not path.exists():
            continue
        with np.load(path) as t:
            frames, dom = t['frames'], t['dom'].astype(float)
        f_eq, final = equilibration_frame(frames, dom)
        rows.append(dict(base_name=row['base_name'],
                         seq_name1=row['seq_name1'], seq_name2=row['seq_name2'],
                         D_norm=row['D_norm'], dGij=row['dGij'],
                         dom_initial=float(dom[0]) if np.isfinite(dom[0]) else np.nan,
                         dom_final=final, equil_frame=f_eq,
                         n_points=int(np.isfinite(dom).sum())))
    if not rows:
        sys.exit("no traces found -- run the trace stage first")
    out = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    print(f"Wrote {len(out)} rows to {out_csv}")
    return out


# ======================================================================
# plot
# ======================================================================

def make_plot(df, summary, out_pdf=OUT_PDF):
    # Adjusted to 2 subplots for Panel A and Panel B only
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    norm = Normalize(0, 1)
    cmap = plt.get_cmap('viridis')

    lookup = df.set_index('base_name')

    # A: every trace, coloured by how demixed the system ends up
    ax = axes[0]
    for _, s in summary.iterrows():
        with np.load(trace_path(s['base_name'])) as t:
            frames, dom = t['frames'], t['dom'].astype(float)
        ax.plot(frames, dom, lw=0.8, alpha=0.5,
                color=cmap(norm(s['D_norm'])))
    ax.axvline(EQUIL_LINE, color='k', ls='--', lw=2)
    ax.text(EQUIL_LINE, 1.02, f'  analysis starts\n  at frame {EQUIL_LINE}',
            fontsize=12, va='bottom', ha='left')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Degree of mixing')
    ax.set_ylim(0, 1.05)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                 label='Demixing index $D$')

    # B: single aggregate line representation of degree of mixing / final value
    ax = axes[1]
    
    # Collect all traces into a DataFrame aligned by frame
    all_traces = {}
    for _, s in summary.iterrows():
        if not np.isfinite(s['dom_final']) or s['dom_final'] <= 0:
            continue
        with np.load(trace_path(s['base_name'])) as t:
            frames, dom = t['frames'], t['dom'].astype(float)
        
        y_raw = dom / s['dom_final']
        # Align data by frame index
        all_traces[s['base_name']] = pd.Series(y_raw, index=frames)
        
    df_all = pd.DataFrame(all_traces)
    
    # Calculate mean and standard deviation across all systems per frame
    y_mean = df_all.mean(axis=1)
    y_std = df_all.std(axis=1)
    
    # Apply rolling mean smoothing (window size 15)
    y_mean_smooth = y_mean.rolling(window=15, min_periods=1, center=True).mean()
    y_std_smooth = y_std.rolling(window=15, min_periods=1, center=True).mean()
    
    common_frames = df_all.index
    
    # Plot the aggregated mean line and standard deviation error band
    ax.plot(common_frames, y_mean_smooth, lw=2.5, color='#3d405b', label='Mean')
    ax.fill_between(common_frames, 
                    y_mean_smooth - y_std_smooth, 
                    y_mean_smooth + y_std_smooth, 
                    color='#3d405b', alpha=0.25, lw=0, label='±1 Std Dev')

    ax.axhline(1.0, color='k', ls=':', lw=2)
    ax.axvline(EQUIL_LINE, color='k', ls='--', lw=2)
    ax.set_xlabel('Frame')
    ax.set_ylabel('Degree of mixing / final value')
    ax.set_ylim(0, 1.6)
    
    # Add a legend for clarity
    ax.legend(loc='upper left', fontsize=14)

    # Label axes (A and B)
    for i, ax in enumerate(axes):
        ax.text(-0.12, 1.10, chr(65 + i), transform=ax.transAxes,
                fontsize=22, fontweight='bold', va='top', ha='right')

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f"Saved -> {out_pdf}")


# ======================================================================

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('job_index', nargs='?', type=int, default=None)
    p.add_argument('--n-systems', action='store_true')
    p.add_argument('--plot-only', action='store_true')
    p.add_argument('--force', action='store_true')
    args = p.parse_args()

    # 50:50 mixtures that formed a dense phase
    df = system_rows()
    df = df[(df['comp1'] == 50) & df['two_phase'].astype(bool)]
    df = df.reset_index(drop=True)

    if args.n_systems:
        print(len(df))
        return

    if args.plot_only:
        make_plot(df, build_summary(df))
        return

    if args.job_index is not None:
        if not 0 <= args.job_index < len(df):
            sys.exit(f"index {args.job_index} outside 0-{len(df) - 1}")
        row = df.iloc[args.job_index]
        print(f"[{args.job_index}] {row['base_name']}")
        try:
            print(f"  {ensure_trace(row, force=args.force)}")
        except Exception as e:
            sys.exit(f"  ! {type(e).__name__}: {e}")
        return

    print(f"{len(df)} equimolar systems with a dense phase; "
          f"frames 0-{LAST_FRAME} every {STRIDE}")
    done = failed = 0
    for _, row in tqdm(df.iterrows(), total=len(df), desc='systems',
                       unit='sys'):
        try:
            ensure_trace(row, force=args.force)
            done += 1
        except Exception as e:
            tqdm.write(f"  ! {row['base_name']}: {type(e).__name__}: {e}")
            failed += 1
    print(f"{done} traces available, {failed} failed")

    make_plot(df, build_summary(df))


if __name__ == '__main__':
    main()