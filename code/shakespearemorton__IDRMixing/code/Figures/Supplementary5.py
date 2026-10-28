#!/usr/bin/env python3
"""Supplementary figure: P30536 paired with four partners at 50:50.

Three rows, four columns, one column per system:

  row 1   deliberately blank, sized and aligned for a trajectory snapshot
          to be dropped in later.  The final frame of each system is
          written out as a centred pdb so the render lines up with the
          profile beneath it.
  row 2   time-averaged 1D concentration profiles along z, in mM for each
          species, sharing one y axis so the four columns are comparable.
  row 3   the concentration space: one point per voxel, the two species'
          residue concentrations plotted against each other, annotated with
          the dense and dilute concentrations and the demixing index.

Profiles come from the cache written by Supplementary2.py and the voxel
maps from DDM_Heterotypic/results/main, so nothing here re-derives a quantity that the
production analysis already produced.

Usage
-----
    python Supplementary5.py              # figure + snapshots
    python Supplementary5.py --no-pdb     # skip the snapshot export
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from DDM import centre_slab

from Supplementary2 import (system_rows, ensure_profile, profile_path,
                            TRAJ_DIR)

plt.rcParams.update({
    'font.size': 18, 'axes.linewidth': 2, 'axes.spines.top': False,
    'axes.spines.right': False, 'xtick.major.size': 6, 'xtick.major.width': 2,
    'ytick.major.size': 6, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300,
})

# ---------------------------------------------------------------- paths
RESULTS_DIR = Path("../DDM_Heterotypic/results/main")
HET_CSV = Path("../DDM_Heterotypic/heterotypic_main.csv")
DATA_DIR = Path("supplementary5_data")
OUT_PDF = Path("Supplementary_Figure5.pdf")

# ---------------------------------------------------------------- systems
SYSTEMS = ['P30536_1_169_Q96RR1_1_152_50_50',
           'P30536_1_169_P0DMR3_1_200_50_50',
           'P30536_1_169_P18089_200_362_50_50',
           'P30536_1_169_Q96AA8_412_477_50_50']

# P30536 is the shared component, so it keeps the darkest colour in every
# column; the partners take the remaining four in order.
ANCHOR = 'P30536'
ANCHOR_COLOR = '#3d405b'
PARTNER_COLORS = ['#e07a5f', '#81b29a', '#f2cc8f', '#f4f1de']

# The palette's lightest colour is invisible on white, so every partner line
# gets the same thin dark outline: applied to all of them, not just the pale
# one, so the four columns look like a set.
STROKE = [pe.withStroke(linewidth=6.5, foreground='#3d405b', alpha=0.55)]

ATOMS_A3_TO_mM = 1.660539e6
SHARE_Y = True          # one concentration axis across all four columns
LOG_FLOOR = 1e-3        # never-occupied voxels are clipped here for log axes
KDE_GRID = 512
PEAK_PROMINENCE = 0.05  # of the KDE maximum; absolute-height cuts are
                        # swamped by the dilute spike


# ======================================================================
# data
# ======================================================================

def load_row(df, base_name):
    hit = df[df['base_name'] == base_name]
    if hit.empty:
        raise KeyError(f"{base_name} not in {HET_CSV}")
    return hit.iloc[0]


def load_profile(row):
    """1D profiles in mM per species."""
    ensure_profile(row)
    with np.load(profile_path(row['base_name']), allow_pickle=True) as p:
        z = p['z_centers']
        d1, d2 = p['density_s1'], p['density_s2']
        len1, len2 = int(p['len1']), int(p['len2'])
    # c[mM] = n_res / (V * N) * 1.6605e6, with the density already n_res/V
    return z, d1 / max(len1, 1) * ATOMS_A3_TO_mM, d2 / max(len2, 1) * ATOMS_A3_TO_mM


def load_voxels(base_name):
    """Per-voxel concentrations of the two species, converted to mM.

    The voxel maps -- and every concentration column in
    heterotypic_main.csv -- are in residues/voxel, because that is the unit
    in which the two species sum to the total.  Converting to mM divides
    each species by its own chain length, so the two axes are rescaled by
    different factors; k1 and k2 are returned so the reported
    concentrations and the composition line can be converted the same way.
    """
    path = RESULTS_DIR / f"{base_name}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=True) as p:
        tot = p['avg_tot'].astype(np.float64).ravel()
        s1 = p['avg_s1'].astype(np.float64).ravel()
        voxel_vol = float(np.prod(p['chunk_dims']))       # A^3
        len1, len2 = int(p['len1']), int(p['len2'])
    s2 = tot - s1

    # c[mM] = n_res / (V_voxel * N_chain) * 1.6605e6
    k1 = ATOMS_A3_TO_mM / (voxel_vol * max(len1, 1))
    k2 = ATOMS_A3_TO_mM / (voxel_vol * max(len2, 1))
    return tot, s1 * k1, s2 * k2, k1, k2


def short_name(full):
    """P30536_1_169 -> P30536, for legends and annotations."""
    return str(full).split('_')[0]


def write_snapshot(row, out_dir=DATA_DIR):
    """Final frame, centred the same way the profiles were, as a pdb."""
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
    ts = u.trajectory[-1]
    lz = u.dimensions[2]
    edges = np.arange(0, lz + 1, 1)
    centre_slab(u.atoms, ts, edges[:-1] + 0.5, edges, lz)

    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{base_name}_final.pdb"
    u.atoms.write(str(out))
    return out, len(u.trajectory)


# ======================================================================
# panels
# ======================================================================

def plot_profile(ax, z, c1, c2, name1, name2, colors):
    ax.plot(z, c1, color=colors[0], lw=3.5, label=name1, solid_capstyle='round',
            zorder=3)
    ax.plot(z, c2, color=colors[1], lw=3.5, label=name2, solid_capstyle='round',
            zorder=2, path_effects=STROKE)
    ax.set_xlim(z[0], z[-1])
    ax.set_xlabel('Z position (Å)', fontsize=16)
    ax.legend(loc='upper left', fontsize=13, handlelength=1.4,
              borderaxespad=0.2, labelspacing=0.25)


def plot_concentration_space(ax, d1, d2, row, name1, name2, colors,
                             lim, tot_max, k1, k2):
    """Voxel concentrations on log axes, both species in mM.

    Empty voxels sit at exactly zero, which a log axis cannot show, so they
    are clipped up to LOG_FLOOR -- the dilute phase then piles onto the
    bottom-left corner rather than vanishing.

    The dashed line is the constant-total contour.  The two species sum to a
    fixed residue count, so in mM it is c1/k1 + c2/k2 = const: the
    conversion tilts it off the diagonal by the ratio of the two chain
    lengths, which is the point of showing it.
    """
    x = np.clip(d1, LOG_FLOOR, None)
    y = np.clip(d2, LOG_FLOOR, None)
    ax.scatter(x, y, s=26, alpha=0.45, color=colors[1],
               edgecolors=ANCHOR_COLOR, lw=0.35, rasterized=True, zorder=2)

    t = np.linspace(0, tot_max, 500)          # residues/voxel taken by s1
    ax.plot(np.clip(t * k1, LOG_FLOOR, None),
            np.clip((tot_max - t) * k2, LOG_FLOOR, None),
            ls='--', color='0.45', lw=2, zorder=1)

    lim = tot_max * max(k1, k2) * 1.15
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(f'{name1} (mM)', fontsize=16)
    ax.set_ylabel(f'{name2} (mM)', fontsize=16)

    # the csv is in residues/voxel, so the reported values are converted
    # here with the same per-species factors used for the scatter
    def fmt(v, k):
        return '--' if not np.isfinite(v) else f'{v * k:.1f}'

    txt = (f"{name1}   $c_{{den}}$ {fmt(row['cden1'], k1)}  "
           f"$c_{{dil}}$ {fmt(row['cdil1'], k1)}\n"
           f"{name2}   $c_{{den}}$ {fmt(row['cden2'], k2)}  "
           f"$c_{{dil}}$ {fmt(row['cdil2'], k2)}\n"
           f"$D$ = {row['D_norm']:.3f}")
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, va='top', ha='right',
            fontsize=11, linespacing=1.5, zorder=5,
            bbox=dict(boxstyle='round,pad=0.45', fc='white', ec='#3d405b',
                      lw=1.0, alpha=0.95))


# ======================================================================

def plot_kde_check(panels, out_dir=DATA_DIR):
    """Diagnostic only: the voxel-concentration KDEs behind each column.

    Combined distribution first, then each IDR's contribution.

    Two things to keep straight when reading this:

    * Peaks are selected by PROMINENCE, not by height relative to the
      maximum.  The dilute spike is an order of magnitude taller than the
      dense peak, so a threshold scaled to y.max() sits below the noise and
      every ripple registers.  The combined panel then marks only the two
      phases the analysis itself reports: the lowest-concentration maximum
      is dilute, the tallest of the rest is dense.

    * cden1/cden2 and cdil1/cdil2 are each species' own KDE peaks, found by
      DDM.py from that species' voxel map, so they need not sum to cden and
      cdil.  They are drawn on the species panels as reference lines.

    Not part of the manuscript figure.
    """
    from scipy.stats import gaussian_kde
    from scipy.signal import find_peaks

    n = len(panels)
    fig, axes = plt.subplots(n, 3, figsize=(18, 3.6 * n), squeeze=False)

    for r, p in enumerate(panels):
        row = p['row']
        # v1/v2 are already mM; the combined map is residues/voxel and has
        # no single conversion, so each panel carries its own axis
        series = [(p['tot'], 'Combined', 'black', 'Res./Voxel'),
                  (p['v1'], p['n1'], p['colors'][0], 'mM'),
                  (p['v2'], p['n2'], p['colors'][1], 'mM')]

        for c, (data, label, color, unit) in enumerate(series):
            ax = axes[r, c]
            x = np.linspace(0, data.max() * 1.05, KDE_GRID)
            if np.ptp(data) <= 0:
                ax.text(0.5, 0.5, 'no variation', ha='center', va='center',
                        transform=ax.transAxes, fontsize=11)
                continue
            y = gaussian_kde(data, bw_method='scott')(x)
            ax.fill_between(x, 0, y, color=color, alpha=0.35)
            ax.plot(x, y, color=color, lw=2,
                    path_effects=STROKE if c == 2 else None)

            peaks, props = find_peaks(y, prominence=y.max() * PEAK_PROMINENCE)

            if c == 0 and peaks.size:
                # the analysis's own rule: lowest maximum is the dilute
                # phase, tallest of the remainder is the dense phase
                order = peaks[np.argsort(x[peaks])]
                marks = [(order[0], 'dilute')]
                rest = order[1:]
                if rest.size:
                    marks.append((rest[np.argmax(y[rest])], 'dense'))
            else:
                marks = [(pk, '') for pk in peaks]

            for pk, tag in marks:
                ax.axvline(x[pk], color='0.4', ls='--', lw=1)
                ax.text(x[pk], y.max() * 1.02,
                        f'{x[pk]:.2f}' + (f' ({tag})' if tag else ''),
                        ha='center', va='bottom', fontsize=9, rotation=90)
                print(f"  {row['base_name']} {label} peak {x[pk]:.2f} "
                      f"(prom {props['prominences'][list(peaks).index(pk)]:.4f})"
                      + (f" -> {tag}" if tag else ""))

            # species panels: where DDM.py's per-species peaks sit
            if c > 0:
                k = p['k1'] if c == 1 else p['k2']
                den = (row['cden1'] if c == 1 else row['cden2']) * k
                dil = (row['cdil1'] if c == 1 else row['cdil2']) * k
                for v, ls in ((den, '-'), (dil, ':')):
                    if np.isfinite(v):
                        ax.axvline(v, color=ANCHOR_COLOR, ls=ls, lw=1.6,
                                   alpha=0.8)

            ax.set_yticks([])
            ax.set_title(label, fontsize=13, fontweight='bold')
            ax.set_xlabel(f'Concentration ({unit})', fontsize=13)

        axes[r, 0].text(
            0.97, 0.95,
            f"total   den {row['cden']:.2f} / dil {row['cdil']:.2f}  "
            f"(Res./Voxel)\n"
            f"{p['n1']}   den {row['cden1'] * p['k1']:.1f} / "
            f"dil {row['cdil1'] * p['k1']:.2f} mM\n"
            f"{p['n2']}   den {row['cden2'] * p['k2']:.1f} / "
            f"dil {row['cdil2'] * p['k2']:.2f} mM\n"
            f"D = {row['D_norm']:.3f}",
            transform=axes[r, 0].transAxes, va='top', ha='right', fontsize=9,
            bbox=dict(boxstyle='round', fc='white', ec='0.6', alpha=0.95))
        axes[r, 0].set_ylabel(row['base_name'], fontsize=10)

    fig.suptitle('KDE check (diagnostic).  Species panels: solid line = '
                 'reported $c_{den,k}$, dotted = $c_{dil,k}$ -- each species\' '
                 'own DDM.py peaks.',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / 'kde_check.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  KDE check -> {out}")


# ======================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-pdb', action='store_true',
                    help='skip writing the final-frame snapshots')
    args = ap.parse_args()

    df = system_rows(HET_CSV) if HET_CSV.exists() else system_rows()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # ---- gather everything first, so the shared axes can be sized ----
    panels = []
    for base_name, partner_color in zip(SYSTEMS, PARTNER_COLORS):
        row = load_row(df, base_name)
        n1, n2 = short_name(row['seq_name1']), short_name(row['seq_name2'])
        # the anchor keeps the dark colour whichever slot it occupies
        colors = ((ANCHOR_COLOR, partner_color) if n1 == ANCHOR
                  else (partner_color, ANCHOR_COLOR))
        z, c1, c2 = load_profile(row)
        tot, v1, v2, k1, k2 = load_voxels(base_name)
        panels.append(dict(row=row, n1=n1, n2=n2, colors=colors,
                           z=z, c1=c1, c2=c2, v1=v1, v2=v2, tot=tot,
                           k1=k1, k2=k2))
        print(f"{base_name}: {n1}/{n2}, {tot.size} voxels, "
              f"D = {row['D_norm']:.3f}")

    ymax = max(max(p['c1'].max(), p['c2'].max()) for p in panels) * 1.12
    # both panel rows are in mM, so one limit serves the profiles and the
    # concentration space alike
    vlim = max(max(p['v1'].max(), p['v2'].max()) for p in panels) * 1.15

    # ---- layout: blank snapshot row, profiles, concentration space ----
    fig = plt.figure(figsize=(23, 14.5))
    gs = gridspec.GridSpec(3, 4, figure=fig,
                           height_ratios=[0.52, 0.80, 1.30],
                           hspace=0.42, wspace=0.30)

    for col, p in enumerate(panels):
        # row 1: reserved for the snapshot. The axes exist (so the column
        # widths match exactly) but nothing is drawn in them.
        ax_snap = fig.add_subplot(gs[0, col])
        ax_snap.set_axis_off()
        ax_snap.text(0.5, 1.06, f"{p['n1']} + {p['n2']}",
                     transform=ax_snap.transAxes, ha='center', va='bottom',
                     fontsize=19, fontweight='bold', color='black')
        ax_snap.text(-0.09, 1.15, chr(65 + col), transform=ax_snap.transAxes,
                     fontsize=24, fontweight='bold', va='top', ha='right')

        ax_prof = fig.add_subplot(gs[1, col])
        plot_profile(ax_prof, p['z'], p['c1'], p['c2'], p['n1'], p['n2'],
                     p['colors'])
        if SHARE_Y:
            ax_prof.set_ylim(0, ymax)
        if col == 0:
            ax_prof.set_ylabel('Concentration (mM)', fontsize=16)
        else:
            ax_prof.set_yticklabels([])

        ax_conc = fig.add_subplot(gs[2, col])
        plot_concentration_space(ax_conc, p['v1'], p['v2'], p['row'],
                                 p['n1'], p['n2'], p['colors'], vlim,
                                 p['tot'].max(), p['k1'], p['k2'])

    fig.savefig(OUT_PDF, bbox_inches='tight')
    print(f"\nSaved -> {OUT_PDF}")

    # ---- source data + snapshots ----
    for p in panels:
        base_name = p['row']['base_name']
        pd.DataFrame({'z_A': p['z'], f"{p['n1']}_mM": p['c1'],
                      f"{p['n2']}_mM": p['c2']}).to_csv(
            DATA_DIR / f"{base_name}_profile.csv", index=False)
        pd.DataFrame({f"{p['n1']}_mM": p['v1'],
                      f"{p['n2']}_mM": p['v2']}).to_csv(
            DATA_DIR / f"{base_name}_voxels.csv", index=False)

    plot_kde_check(panels)

    if not args.no_pdb:
        for p in panels:
            try:
                out, n = write_snapshot(p['row'])
                print(f"  snapshot {out.name}  (frame {n} of {n}, centred)")
            except Exception as e:
                print(f"  ! {p['row']['base_name']}: {type(e).__name__}: {e}")

    print(f"Data and snapshots in {DATA_DIR}/")


if __name__ == '__main__':
    main()