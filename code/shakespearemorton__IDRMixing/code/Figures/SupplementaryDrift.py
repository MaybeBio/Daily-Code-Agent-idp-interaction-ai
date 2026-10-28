"""
python SupplementaryDrift.py ../../Extra_Jobs/A1_7b2459a9
"""
import argparse
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import MDAnalysis as mda

RC_BASE = {
    'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2,
    'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300, 'font.size': 20,
}
plt.rcParams.update(RC_BASE)

COLORS = {'dense': '#e02b35', 'mid': '#EBECED', 'dilute': '#2066a8'}
NA = 6.02214076e23  # Avogadro's number

# Same convention as analyze_a1_phase.py: backup_<name>_<YYYYDDMM>_<HH>h<MM>m<SS>s.dcd
BACKUP_RE = re.compile(r'^backup_.*_(\d{8})_(\d{2})h(\d{2})m(\d{2})s\.dcd$')


def dcd_sort_key(path):
    """Backups first (chronological, oldest to newest, by the timestamp in
    their filename), then the live trajectory file last -- it's the one
    still being appended to, so its frames come after every backup.
    (Identical to the sort key in analyze_a1_phase.py.)
    """
    m = BACKUP_RE.match(path.name)
    if not m:
        return (1, path.name)
    date_str, hh, mm, ss = m.groups()
    stamp = datetime.strptime(f"{date_str}{hh}{mm}{ss}", "%Y%d%m%H%M%S")
    return (0, stamp)


def find_trajectory_files(directory, pdb_name='top.pdb'):
    """Locate the topology and every trajectory chunk in `directory`,
    ordered so MDAnalysis concatenates them correctly."""
    directory = Path(directory)
    top = directory / pdb_name
    if not top.exists():
        raise FileNotFoundError(f"no {pdb_name} in {directory}")
    dcds = sorted(directory.glob("*.dcd"), key=dcd_sort_key)
    if not dcds:
        raise FileNotFoundError(f"no .dcd files in {directory}")
    return top, dcds


def compute_z_time_density(pdb_path, dcd_paths, voxel_size=20.0, n_chains=None):
    """Bin atom positions along Z into `voxel_size` (Angstrom) slabs for
    every frame of the (possibly multi-file) trajectory, averaging over
    X/Y.

    Parameters
    ----------
    dcd_paths : path-like or list of path-like
        One or more trajectory files, in the order they should be
        concatenated (MDAnalysis chains them automatically).
    n_chains : int or None
        Number of protein chains. If None, auto-detected as
        len(u.segments); pass an int to override.

    Returns
    -------
    density : (n_frames, n_bins) array. In mM if a chain count is
        available (given or auto-detected), otherwise beads/nm^3.
    z_centers_nm : bin centers along Z, in nm, shifted so 0 = box center.
    times_us : per-frame time, in microseconds.
    is_concentration : bool, whether `density` is in mM (True) or
        beads/nm^3 (False).
    """
    if isinstance(dcd_paths, (str, Path)):
        dcd_paths = [dcd_paths]
    u = mda.Universe(str(pdb_path), [str(d) for d in dcd_paths])
    n_atoms = len(u.atoms)

    if n_chains is None:
        n_chains = len(u.segments) or None

    box = u.dimensions[:3]
    n_bins = max(int(round(box[2] / voxel_size)), 1)
    bin_width = box[2] / n_bins
    voxel_vol_nm3 = (bin_width * box[0] * box[1]) / 1000.0  # A^3 -> nm^3

    n_frames = len(u.trajectory)
    counts = np.zeros((n_frames, n_bins))
    times_us = np.zeros(n_frames)

    for f_idx, ts in enumerate(u.trajectory):
        z = u.atoms.positions[:, 2]
        z_idx = np.clip((z // bin_width).astype(int), 0, n_bins - 1)
        counts[f_idx] = np.bincount(z_idx, minlength=n_bins)
        times_us[f_idx] = ts.time / 1e6  # ps -> us

    z_edges_nm = np.linspace(0, box[2], n_bins + 1) / 10.0
    z_centers_nm = 0.5 * (z_edges_nm[:-1] + z_edges_nm[1:]) - (box[2] / 20.0)

    if n_chains:
        atoms_per_chain = n_atoms / n_chains
        density = counts / atoms_per_chain / (voxel_vol_nm3 * 1e-24) / NA * 1000.0
        is_concentration = True
    else:
        density = counts / voxel_vol_nm3
        is_concentration = False

    return density, z_centers_nm, times_us, is_concentration


def plot_drift_heatmap(density, z_centers_nm, times_us, is_concentration,
                        output_path="SupplementaryDrift.pdf"):
    cmap = LinearSegmentedColormap.from_list(
        'custom', [(0.0, 'white'), (1.0, COLORS['dilute'])])

    fig, ax = plt.subplots(figsize=(4.5, 5.5))

    positive = density[density > 0]
    vmin = positive.min() if positive.size else 1e-2
    vmax = density.max() if density.max() > 0 else 1.0

    z_width = z_centers_nm[1] - z_centers_nm[0] if len(z_centers_nm) > 1 else 1.0
    extent = [z_centers_nm[0] - z_width / 2, z_centers_nm[-1] + z_width / 2,
              times_us[0], times_us[-1]]

    im = ax.imshow(density, origin='lower', aspect='auto', cmap=cmap,
                    norm=LogNorm(vmin=vmin, vmax=vmax), extent=extent)

    ax.set_xlabel('z (nm)')
    ax.set_ylabel(r'Time ($\mu$s)')

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Concentration (mM)' if is_concentration else r'Density (beads/nm$^3$)')

    fig.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate a Z vs. time density heatmap (SupplementaryDrift.pdf) "
                     "from a job directory containing top.pdb + one or more .dcd chunks")
    parser.add_argument('directory', type=Path,
                         help='Job directory, e.g. Extra_Jobs/A1_7b2459a9')
    parser.add_argument('--pdb-name', default='top.pdb',
                         help='Topology filename to look for (default: top.pdb)')
    parser.add_argument('--voxel-size', type=float, default=20.0,
                         help='Z-slab thickness in Angstrom (default: 20)')
    parser.add_argument('--n-chains', type=int, default=None,
                         help='Number of protein chains, to convert to mM. '
                              'Auto-detected from the topology if omitted.')
    parser.add_argument('--output', type=Path, default=Path('SupplementaryDrift.pdf'))
    args = parser.parse_args()

    top, dcds = find_trajectory_files(args.directory, pdb_name=args.pdb_name)
    print(f"Topology: {top}")
    print(f"Trajectory chunks ({len(dcds)}), in concatenation order:")
    for d in dcds:
        print(f"    {d.name}")

    density, z_centers_nm, times_us, is_conc = compute_z_time_density(
        top, dcds, voxel_size=args.voxel_size, n_chains=args.n_chains)
    out = plot_drift_heatmap(density, z_centers_nm, times_us, is_conc, output_path=args.output)
    print(f"Wrote {out}")


if __name__ == '__main__':
    main()