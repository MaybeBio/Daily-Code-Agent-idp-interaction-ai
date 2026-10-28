#!/usr/bin/env python3
"""Supplementary Figure 2: the demixing index against the local mixing order
parameter, across three composition regimes.

The order parameter q of Chen and Jacobs (JCTC 2024) is the cosine
similarity between the local and the global composition, averaged over the
dense phase.  Unlike the demixing index it needs the dense phase to be
delimited, which is what the three panels are there to show: the two agree
closely at equimolar compositions and diverge as the mixture becomes
asymmetric and the minority component's dense phase becomes narrow.

A second file, Supplementary_Figure2_profiles.pdf, shows the density,
composition and q(z) profiles of two representative systems with the
dense-phase boundary marked, which is where that sensitivity comes from.

Three stages, so the slow one only runs once:

  profiles  centre every frame on the slab and histogram each species along
            z, writing one npz per system to DDM_Heterotypic/z_profiles/
  metrics   reduce each profile to q, cached in DATASETS/mixing_metrics.csv
  plot      D against q, split by composition regime, plus the example
            profiles

The centring is imported from DDM.py, so the slab is positioned exactly as
it is everywhere else in the paper.

Usage
-----
    python Supplementary2.py                 # profiles (missing only), then plot
    python Supplementary2.py 17              # one system, for a batch array
    python Supplementary2.py --n-systems     # print the job count and exit
    python Supplementary2.py --plot-only     # re-plot from the cached csv
    python Supplementary2.py --force         # recompute every profile
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# The slab centring lives with the production analysis.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from DDM import calc_zpatch, centre_slab, DEFAULT_START

plt.rcParams.update({
    'font.size': 18, 'axes.linewidth': 2, 'axes.spines.top': False,
    'axes.spines.right': False, 'xtick.major.size': 6, 'xtick.major.width': 2,
    'ytick.major.size': 6, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300,
})

# ---------------------------------------------------------------- paths
HET_CSV = Path("../DDM_Heterotypic/heterotypic_main.csv")
TRAJ_DIR = Path("../TRAJECTORIES")
PROFILE_DIR = Path("../DDM_Heterotypic/z_profiles")
METRICS_CSV = Path("../DATASETS/mixing_metrics.csv")
OUT_PDF = Path("Supplementary_Figure2.pdf")
OUT_PROFILE_PDF = Path("Supplementary_Figure2_profiles.pdf")

# ---------------------------------------------------------------- settings
N_Z_BINS = 100          # Methods sec. G: the z axis is discretised into 100 bins
START_FRAME = DEFAULT_START
STRIDE = 1              # analyse every Nth frame after equilibration

MIN_DENSE_FRAC = 0.05   # two phases only if >= 5% of voxels sit at or above c_den
DG_CUTOFF = -3.0        # stable heterotypic condensate only
MAX_COMP = 98           # neither species may exceed this % of residues

COLOR = 'grey'
BINS = [('minority', 'IDR-1 < 15% or > 85%'),
        ('unequal',  'IDR-1 15-35% or 65-85%'),
        ('equimolar', 'IDR-1 35-65%')]

# Two representative systems for the profile figure: one asymmetric mixture
# where the minority component sits at the interface, and one where it forms
# a gradient through the dense phase.  (base_name, colour 1, colour 2)
EXAMPLES = [('Q8WTT2_1_202_Q3LI63_1_56_19_81', '#e74c3c', '#3498db'),
            ('P18089_200_362_Q3LI63_1_56_75_25', '#DCC670', '#3498db')]

A3_TO_NM3 = 1000.0      # residues/A^3 -> residues/nm^3


# ======================================================================
# system table
# ======================================================================

def composition_bin(comp1):
    """Composition regime, as defined in Methods sec. G."""
    if comp1 < 15 or comp1 > 85:
        return 'minority'
    if comp1 < 35 or comp1 > 65:
        return 'unequal'
    return 'equimolar'


def apply_two_phase_gate(df):
    """A voxel map resolves two phases only when at least MIN_DENSE_FRAC of
    its voxels sit at or above c_den. Otherwise c_dil is the mean over all
    voxels and no dense value (and so no dGij) is defined."""
    df = df.copy()
    df['two_phase'] = df['frac_vox_ge_cden'] >= MIN_DENSE_FRAC
    one = ~df['two_phase']
    for dil in ('cdil', 'cdil1', 'cdil2'):
        df.loc[one, dil] = df.loc[one, f'{dil}_onephase']
    df.loc[one, ['cden', 'cden1', 'cden2', 'dGij', 'dGij_err']] = np.nan
    return df


def system_rows(het_csv=HET_CSV):
    """Every system that qualifies, in a stable order.

    Sorted by base_name so the index handed to a batch job means the same
    system between submissions.
    """
    df = apply_two_phase_gate(pd.read_csv(het_csv))
    df = df[df['two_phase']]          # gated-out (single-phase) maps are not analysed
    df = df[(df['dGij'] < DG_CUTOFF)
            & (df['comp1'] <= MAX_COMP) & (df['comp2'] <= MAX_COMP)]
    return df.sort_values('base_name').reset_index(drop=True)


def profile_path(base_name):
    return PROFILE_DIR / f"{base_name}_z_profiles.npz"


# ======================================================================
# profiles
# ======================================================================

def split_species(u, len1, len2):
    """Atom groups for the two species.

    Chains differing in length are split on that; equal-length chains are
    split on their residue-name signature, keeping the order in which the
    two signatures first appear (species 1 is written first by the setup).
    """
    segs = list(u.segments)
    if len1 != len2:
        idx1 = [a.index for s in segs if len(s.atoms) == len1 for a in s.atoms]
        idx2 = [a.index for s in segs if len(s.atoms) == len2 for a in s.atoms]
    else:
        order, groups = [], {}
        for s in segs:
            sig = tuple(r.name for r in s.atoms)
            if sig not in groups:
                groups[sig] = []
                order.append(sig)
            groups[sig].append(s)
        if len(order) != 2:
            raise ValueError(f"expected 2 species, found {len(order)}")
        idx1 = [a.index for s in groups[order[0]] for a in s.atoms]
        idx2 = [a.index for s in groups[order[1]] for a in s.atoms]

    if not idx1 or not idx2:
        raise ValueError("could not split the two species")
    return u.atoms[idx1], u.atoms[idx2]


def compute_profile(base_name, len1, len2, start=START_FRAME, stride=STRIDE):
    """Time-averaged density and composition profiles along z."""
    import MDAnalysis as mda

    pdb = TRAJ_DIR / f"{base_name}.pdb"
    dcds = [TRAJ_DIR / f"{base_name}.dcd"]
    dcds += sorted(TRAJ_DIR.glob(f"{base_name}_resume*.dcd"))
    dcds += sorted(TRAJ_DIR.glob(f"{base_name}_[0-9]*.dcd"))
    dcds = [str(d) for d in dcds if Path(d).is_file()]
    if not pdb.is_file() or not dcds:
        raise FileNotFoundError(f"no trajectory for {base_name}")

    u = mda.Universe(str(pdb), dcds)
    n_total = len(u.trajectory)
    if n_total - start < 2:
        raise ValueError(f"only {n_total - start} frame(s) after equilibration")

    ag1, ag2 = split_species(u, len1, len2)
    ag = u.atoms

    box = u.dimensions[:3].astype(float)
    lz = box[2]
    # 1 A bins for the centring, N_Z_BINS for the profiles themselves
    cedges = np.arange(0, lz + 1, 1)
    cz = cedges[:-1] + 0.5
    pedges = np.linspace(0, lz, N_Z_BINS + 1)
    z_centers = 0.5 * (pedges[:-1] + pedges[1:])
    bin_volume = float(box[0] * box[1] * (pedges[1] - pedges[0]))

    frames = range(start, n_total, stride)
    acc_s1 = np.zeros(N_Z_BINS)
    acc_s2 = np.zeros(N_Z_BINS)

    n_used = 0
    for i in frames:
        centre_slab(ag, u.trajectory[i], cz, cedges, lz)
        acc_s1 += np.histogram(ag1.positions[:, 2], bins=pedges)[0]
        acc_s2 += np.histogram(ag2.positions[:, 2], bins=pedges)[0]
        n_used += 1
        print(f"    frame {n_used}/{len(frames)}", end="\r")
    print(" " * 40, end="\r")

    dens_s1 = acc_s1 / (n_used * bin_volume)      # residues / A^3
    dens_s2 = acc_s2 / (n_used * bin_volume)
    dens_tot = dens_s1 + dens_s2

    with np.errstate(divide='ignore', invalid='ignore'):
        comp_s1 = np.where(dens_tot > 0, dens_s1 / dens_tot, 0.0)
        comp_s2 = np.where(dens_tot > 0, dens_s2 / dens_tot, 0.0)

    return dict(z_centers=z_centers, bin_volume=bin_volume,
                density_total=dens_tot, density_s1=dens_s1, density_s2=dens_s2,
                composition_s1=comp_s1, composition_s2=comp_s2,
                box=box, len1=len1, len2=len2, n_frames=n_used,
                n_frames_total=n_total, start_frame=start)


def write_profile(base_name, prof):
    PROFILE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(profile_path(base_name), **prof)


def ensure_profile(row, force=False):
    """Compute this system's profile unless it is already on disk."""
    path = profile_path(row['base_name'])
    if path.exists() and not force:
        return 'cached'
    prof = compute_profile(row['base_name'], int(row['len1']), int(row['len2']))
    write_profile(row['base_name'], prof)
    return f"{prof['n_frames']} frames"


# ======================================================================
# metric
# ======================================================================

def dense_mask(density):
    """Bins belonging to the dense phase.

    The simplest available definition, as stated in Methods sec. G: every
    bin whose total density exceeds the mean density along z, restricted to
    the largest contiguous run so a stray dilute-phase bin cannot join it.
    """
    idx = np.arange(density.size, dtype=float)
    patch, _ = calc_zpatch(idx, density, cutoff=float(density.mean()))
    mask = np.zeros(density.size, dtype=bool)
    if patch.size:
        mask[patch.astype(int)] = True
    return mask


def order_parameter_profile(comp_s1, comp_s2, phi_target):
    """q(z): cosine similarity between local and global composition."""
    loc = np.stack([comp_s1, comp_s2], axis=1)
    norms = np.linalg.norm(loc, axis=1) * np.linalg.norm(phi_target)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(norms > 0, loc @ phi_target / norms, 0.0)


def order_parameter(comp_s1, comp_s2, phi_target, mask):
    """Mean cosine similarity over the dense phase (Chen and Jacobs)."""
    if not mask.any():
        return np.nan
    return float(np.mean(order_parameter_profile(
        comp_s1, comp_s2, phi_target)[mask]))


def target_composition(row):
    """Global residue fractions, from phi if present and comp otherwise."""
    phi = np.array([row['phi1'], row['phi2']], dtype=float)
    if not np.isfinite(phi).all() or phi.sum() <= 0:
        phi = np.array([row['comp1'], row['comp2']], dtype=float) / 100.0
    return phi


def metrics_for(row):
    """q for one system, from its cached profile."""
    path = profile_path(row['base_name'])
    if not path.exists():
        return None
    with np.load(path, allow_pickle=True) as p:
        dens_tot = p['density_total']
        comp_s1, comp_s2 = p['composition_s1'], p['composition_s2']
        n_frames = int(p['n_frames'])

    if not np.isfinite(dens_tot).all() or dens_tot.sum() <= 0:
        return None

    mask = dense_mask(dens_tot)
    phi = target_composition(row)

    return dict(base_name=row['base_name'],
                seq_name1=row['seq_name1'], seq_name2=row['seq_name2'],
                comp1=row['comp1'], comp2=row['comp2'],
                composition_bin=composition_bin(row['comp1']),
                D_norm=row['D_norm'], dGij=row['dGij'],
                order_parameter=order_parameter(comp_s1, comp_s2, phi, mask),
                n_dense_bins=int(mask.sum()), n_frames=n_frames)


def build_metrics(df, out_csv=METRICS_CSV):
    rows, missing = [], 0
    for _, row in df.iterrows():
        m = metrics_for(row)
        if m is None:
            missing += 1
            continue
        rows.append(m)
    if not rows:
        sys.exit("no profiles found -- run the profile stage first")
    out = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    print(f"Wrote {len(out)} rows to {out_csv}"
          + (f" ({missing} systems had no profile)" if missing else ""))
    return out


# ======================================================================
# plots
# ======================================================================

def make_plot(df, out_pdf=OUT_PDF):
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    print("\n" + "=" * 60)
    print("ORDER PARAMETER vs DEMIXING INDEX")
    print("=" * 60)

    for col, (key, title) in enumerate(BINS):
        ax = axes[col]
        ax.text(-0.12, 1.12, chr(65 + col), transform=ax.transAxes,
                fontsize=22, fontweight='bold', va='top', ha='right')

        d = df[df['composition_bin'] == key].dropna(
            subset=['D_norm', 'order_parameter'])
        if len(d) > 1:
            x, y = d['D_norm'].values, d['order_parameter'].values
            ax.scatter(x, y, s=90, alpha=0.7, color=COLOR,
                       edgecolors='black', lw=0.8, rasterized=True)
            r, _ = pearsonr(x, y)
            fit = np.poly1d(np.polyfit(x, y, 1))
            xl = np.linspace(0, 1, 100)
            ax.plot(xl, fit(xl), 'k--', lw=2)
            ax.text(0.95, 0.95, f'$R^2$ = {r ** 2:.2f}\nn = {len(d)}',
                    transform=ax.transAxes, va='top', ha='right', fontsize=15,
                    bbox=dict(boxstyle='round', fc='white', alpha=0.85))
            print(f"  {key:<10} R2 = {r ** 2:.3f}  (n = {len(d)})")
        else:
            ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                    ha='center', va='center')

        ax.set_xlim(0, 1.05)
        ax.set_ylim(0.7, 1.01)
        ax.grid(True, alpha=0.3)
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.set_xlabel('Demixing index')
        if col == 0:
            ax.set_ylabel('Order parameter $q$')

    d = df.dropna(subset=['D_norm', 'order_parameter'])
    if len(d) > 1:
        r, _ = pearsonr(d['D_norm'], d['order_parameter'])
        print(f"  {'overall':<10} R2 = {r ** 2:.3f}  (n = {len(d)})")
    print("=" * 60)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f"Saved -> {out_pdf}")


def shade_dense(ax, z, mask):
    """Mark the dense-phase boundary used to average q."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return None, None
    dz = z[1] - z[0]
    lo, hi = z[idx[0]] - dz / 2, z[idx[-1]] + dz / 2
    ax.axvline(lo, color='gray', ls='--', lw=3, alpha=0.7)
    ax.axvline(hi, color='gray', ls='--', lw=3, alpha=0.7)
    ax.axvspan(lo, hi, alpha=0.1, color='gray')
    return lo, hi


def plot_profile_comparison(df, out_pdf=OUT_PROFILE_PDF, examples=EXAMPLES):
    """Density, composition and q(z) for two representative systems.

    The shaded band is the dense phase over which q is averaged; the point
    of the figure is how much of q depends on exactly where that band ends.
    """
    lookup = df.set_index('base_name')
    fig, axes = plt.subplots(3, 2, figsize=(11.2, 15))

    for col, (base_name, c1, c2) in enumerate(examples):
        if base_name not in lookup.index:
            axes[0, col].text(0.5, 0.5, f"{base_name}\nnot in dataset",
                              ha='center', va='center', fontsize=10,
                              transform=axes[0, col].transAxes)
            for r in range(3):
                axes[r, col].set_axis_off()
            print(f"  ! {base_name} not in {HET_CSV}")
            continue

        row = lookup.loc[base_name]
        row = row.iloc[0] if isinstance(row, pd.DataFrame) else row
        try:
            ensure_profile({**row.to_dict(), 'base_name': base_name})
        except Exception as e:
            print(f"  ! {base_name}: {type(e).__name__}: {e}")
            continue

        with np.load(profile_path(base_name), allow_pickle=True) as p:
            z = p['z_centers']
            d1, d2 = p['density_s1'], p['density_s2']
            dt = p['density_total']
            comp1, comp2 = p['composition_s1'], p['composition_s2']

        phi = target_composition(row)
        q = order_parameter_profile(comp1, comp2, phi)
        mask = dense_mask(dt)
        name1, name2 = str(row['seq_name1']), str(row['seq_name2'])

        # Row 0: density
        ax = axes[0, col]
        ax.plot(z, d1 * A3_TO_NM3, color=c1, lw=5, label=name1)
        ax.plot(z, d2 * A3_TO_NM3, color=c2, lw=5, label=name2)
        shade_dense(ax, z, mask)
        ax.set_ylabel('Density (Res./nm$^3$)', fontsize=18)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=2,
                  fontsize=16)
        ax.text(-0.1, 1.15, chr(68 + col * 3), transform=ax.transAxes,
                fontsize=24, fontweight='bold', va='top', ha='right')

        # Row 1: composition
        ax = axes[1, col]
        ax.plot(z, comp1, color=c1, lw=5)
        ax.plot(z, comp2, color=c2, lw=5)
        shade_dense(ax, z, mask)
        ax.set_ylabel('Composition', fontsize=18)
        ax.set_ylim(0, 1.05)
        ax.text(-0.1, 1.15, chr(69 + col * 3), transform=ax.transAxes,
                fontsize=24, fontweight='bold', va='top', ha='right')

        # Row 2: order parameter
        ax = axes[2, col]
        ax.plot(z, q, color='#27ae60', lw=5)
        ax.axhline(1.0, color='gray', ls='--', lw=3, alpha=0.5)
        shade_dense(ax, z, mask)
        ax.set_ylabel('Order Parameter', fontsize=18)
        ax.set_xlabel('Z position (Å)', fontsize=18)
        ax.set_ylim(0, 1.1)
        ax.text(-0.1, 1.15, chr(70 + col * 3), transform=ax.transAxes,
                fontsize=24, fontweight='bold', va='top', ha='right')

        print(f"  {base_name}: {name1} {row['comp1']:.0f}% / "
              f"{name2} {row['comp2']:.0f}%, "
              f"{int(mask.sum())} dense bins, "
              f"<q>den = {np.mean(q[mask]):.3f}, D = {row['D_norm']:.3f}")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f"Saved -> {out_pdf}")


# ======================================================================

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('job_index', nargs='?', type=int, default=None,
                   help='analyse a single system and exit (for batch arrays)')
    p.add_argument('--n-systems', action='store_true',
                   help='print the number of systems and exit')
    p.add_argument('--plot-only', action='store_true',
                   help='re-plot from the cached metrics csv')
    p.add_argument('--force', action='store_true',
                   help='recompute profiles even if they exist')
    args = p.parse_args()

    if args.plot_only:
        make_plot(pd.read_csv(METRICS_CSV))
        plot_profile_comparison(system_rows())
        return

    df = system_rows()

    if args.n_systems:
        print(len(df))
        return

    if args.job_index is not None:
        if not 0 <= args.job_index < len(df):
            sys.exit(f"index {args.job_index} outside 0-{len(df) - 1}")
        row = df.iloc[args.job_index]
        print(f"[{args.job_index}] {row['base_name']}")
        try:
            print(f"  {ensure_profile(row, force=args.force)}")
        except Exception as e:
            sys.exit(f"  ! {type(e).__name__}: {e}")
        return

    # no index: fill in whatever is missing, then reduce and plot
    print(f"{len(df)} systems qualify (dGij < {DG_CUTOFF}, "
          f"comp <= {MAX_COMP}%)")
    done = failed = 0
    for i, row in df.iterrows():
        try:
            status = ensure_profile(row, force=args.force)
        except Exception as e:
            print(f"  ! {row['base_name']}: {type(e).__name__}: {e}")
            failed += 1
            continue
        if status != 'cached':
            print(f"  {row['base_name']}: {status}")
        done += 1
    print(f"{done} profiles available, {failed} failed")

    make_plot(build_metrics(df))
    plot_profile_comparison(df)


if __name__ == '__main__':
    main()