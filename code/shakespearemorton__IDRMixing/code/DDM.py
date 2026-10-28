import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

# sys.path.append('PATH/TO/BLOCKING')   # <---- Change this
sys.path.append('/storage/brno14-ceitec/shared/softmatter/shakespearem/BLOCKING')
from main import BlockAnalysis


ATOMS_A3_TO_mM = 1.660539e6    # residues/A^3 -> mM, per chain
KDE_GRID = 4096                # grid for the trajectory-averaged KDEs
DEFAULT_START = 1000           # frames discarded as equilibration
CHUNK_SIZE = 20.0              # voxel side length S, in Angstrom
EPSILON = 1e-9


# --- species detection ---------------------------------------------------

def detect_species(u):
    """Split the topology into species by residue-name signature.

    Returns (ag1, ag2, n_species). ag2 is empty for a single-component
    system. Segments are grouped by their sequence, so chains of the same
    protein land together regardless of order in the file.
    """
    groups = defaultdict(list)
    for seg in u.segments:
        groups[tuple(res.name for res in seg.atoms)].append(seg)
    sigs = list(groups)

    def collect(sig):
        return u.atoms[[a.index for s in groups[sig] for a in s.atoms]]

    if len(sigs) == 1:
        return collect(sigs[0]), u.atoms[[]], 1
    if len(sigs) == 2:
        return collect(sigs[0]), collect(sigs[1]), 2

    print(f"  ! {len(sigs)} species found; analysing as a single component")
    return u.atoms, u.atoms[[]], 1


# --- slab centring -------------------------------------------------------

def calc_zpatch(z, h, cutoff=0):
    """Find the largest contiguous region above cutoff."""
    ct, ct_max, zpatch, hpatch, zwin, hwin = 0., 0., [], [], [], []
    for ix, x in enumerate(h):
        if x > cutoff:
            ct += x; zwin.append(z[ix]); hwin.append(x)
        else:
            if ct > ct_max: ct_max, zpatch, hpatch = ct, zwin, hwin
            ct, zwin, hwin = 0., [], []
    if ct > ct_max: zpatch, hpatch = zwin, hwin
    return np.array(zpatch), np.array(hpatch)


def centre_slab(ag, ts, z, edges, lz):
    """Densest bin to the box centre, then the density-weighted centre of the
    largest contiguous patch.

    Centring uses the total density, so when two condensates coexist and
    diffuse independently the second is smeared across the time-averaged
    map and its density underestimated.
    """
    h, _ = np.histogram(ag.positions[:, 2], bins=edges)
    ag.translate([0, 0, -z[np.argmax(h)] + 0.5 * lz])
    ts = transformations.wrap(ag)(ts)

    h, _ = np.histogram(ag.positions[:, 2], bins=edges)
    zpatch, hpatch = calc_zpatch(z, h)
    if len(zpatch):
        ag.translate([0, 0, -np.average(zpatch, weights=hpatch) + 0.5 * lz])
        ts = transformations.wrap(ag)(ts)
    return ts


# --- peak finding --------------------------------------------------------

def kde_peaks(x, y):
    pk, props = find_peaks(y, height=0)
    if len(pk) == 0:
        return np.array([]), np.array([])
    return x[pk], props['peak_heights']


def peak_fwhm(x, y, pos):
    i = int(np.argmin(np.abs(x - pos)))
    half = 0.5 * y[i]
    lo = i
    while lo > 0 and y[lo] > half:
        lo -= 1
    hi = i
    while hi < y.size - 1 and y[hi] > half:
        hi += 1
    return float(x[lo]), float(x[hi])


def dense_concentration(s):
    """Dense phase from the linear KDE: tallest peak above the lowest one."""
    try:
        kde = gaussian_kde(s)
    except Exception:
        return np.nan, np.nan, np.nan
    pad = 0.05 * np.ptp(s)
    x = np.linspace(s.min() - pad, s.max() + pad, KDE_GRID)
    y = kde(x)
    pos, height = kde_peaks(x, y)
    if len(pos) < 2:
        return np.nan, np.nan, np.nan
    rest = np.argsort(pos)[1:]
    c = float(pos[rest[np.argmax(height[rest])]])
    lo, hi = peak_fwhm(x, y, c)
    return c, lo, hi


def dilute_concentration(concs):
    """Dilute phase from the KDE in ln(c)."""
    occupied = concs[concs > 0]
    if occupied.size < 2:
        return np.nan, np.nan, np.nan
    ls = np.log(occupied)
    try:
        kde = gaussian_kde(ls)
    except Exception:
        return np.nan, np.nan, np.nan
    pad = 0.05 * np.ptp(ls)
    x = np.linspace(ls.min() - pad, ls.max() + pad, KDE_GRID)
    y = kde(x)
    pos, _ = kde_peaks(x, y)
    if not len(pos):
        return np.nan, np.nan, np.nan
    c = float(pos.min())
    lo, hi = peak_fwhm(x, y, c)
    return float(np.exp(c)), float(np.exp(lo)), float(np.exp(hi))


def phase_concentrations(voxel_map):
    return dilute_concentration(voxel_map), dense_concentration(voxel_map)


def dense_fraction(tot, c_den):
    """Fraction of voxels at or above the dense concentration."""
    tot = np.asarray(tot, float).ravel()
    if not np.isfinite(c_den) or tot.size == 0:
        return 0.0
    return float((tot >= c_den).sum() / tot.size)


# --- demixing ------------------------------------------------------------

def demixing(s1, s2, phi1, phi2):
    """Demixing index from the variances along and perpendicular to the
    composition line, normalised by Dmax and capped at 1."""
    pts = np.column_stack((np.asarray(s1, float).ravel(),
                           np.asarray(s2, float).ravel()))
    norm = np.hypot(phi1, phi2)
    p = phi1 * phi2
    if len(pts) < 10 or norm < EPSILON or p < EPSILON:
        return np.nan
    u_comp = np.array([phi1, phi2]) / norm
    u_ortho = np.array([-u_comp[1], u_comp[0]])
    l1 = np.var(pts @ u_comp)
    l2 = np.var(pts @ u_ortho)
    if l1 + l2 < EPSILON:
        return np.nan
    dmax = 2.0 * p / (1.0 - 2.0 * p)
    d = 2.0 * min(l1, l2) / (l1 + l2)
    return float(np.clip(d / dmax, 0.0, 1.0))


# --- domain decomposition ------------------------------------------------

def analyze_3d_chunks(u, ag1, ag2, start=DEFAULT_START, end=None,
                      chunk_size=CHUNK_SIZE):
    """Time-averaged voxel maps, per-frame total counts, and per-frame D.

    Counts are residues per voxel. Species 2 is accumulated directly rather
    than as (avg_tot - avg_s1): the subtraction leaves ~1e-15 residues in
    genuinely empty voxels, which survive the `concs > 0` filter and drag
    the log-space dilute KDE down.
    """
    box = u.dimensions[:3]
    n_chunks = np.maximum(np.round(box / chunk_size).astype(int), 1)
    chunk_dims = box / n_chunks
    ny, nz = n_chunks[1], n_chunks[2]
    n_voxels = int(np.prod(n_chunks))

    ag = u.atoms
    binary = len(ag2) > 0
    s1_indices = ag1.indices if binary else np.empty(0, dtype=int)
    phi1 = len(ag1) / len(ag) if binary else 1.0
    phi2 = 1.0 - phi1

    lz = u.dimensions[2]
    edges = np.arange(0, lz + 1, 1)
    z = edges[:-1] + 0.5

    stop = len(u.trajectory) if end is None else min(len(u.trajectory), end)
    n_frames = stop - start
    if n_frames <= 0:
        raise ValueError(f"no frames to analyse: start={start}, stop={stop}, "
                         f"trajectory has {len(u.trajectory)} frames")

    avg_tot = np.zeros(n_voxels, dtype=np.float64)
    avg_s1 = np.zeros(n_voxels, dtype=np.float64)
    avg_s2 = np.zeros(n_voxels, dtype=np.float64)
    tot_ts = np.zeros((n_frames, n_voxels), dtype=np.uint16)
    d_series = np.full(n_frames, np.nan, dtype=np.float64)

    for i in range(n_frames):
        centre_slab(ag, u.trajectory[i + start], z, edges, lz)

        idx = np.clip((ag.positions // chunk_dims).astype(int), 0, n_chunks - 1)
        flat_idx = idx[:, 0] * ny * nz + idx[:, 1] * nz + idx[:, 2]
        counts = np.bincount(flat_idx, minlength=n_voxels)
        tot_ts[i] = counts
        avg_tot += (counts - avg_tot) / (i + 1)

        if binary:
            counts1 = np.bincount(flat_idx[s1_indices], minlength=n_voxels)
            counts2 = counts - counts1          # exact: integer difference
            avg_s1 += (counts1 - avg_s1) / (i + 1)
            avg_s2 += (counts2 - avg_s2) / (i + 1)
            d_series[i] = demixing(counts1, counts2, phi1, phi2)

    return dict(avg_tot=avg_tot, avg_s1=avg_s1, avg_s2=avg_s2,
                tot_ts=tot_ts, d_series=d_series,
                n_frames=n_frames, chunk_dims=chunk_dims, n_chunks=n_chunks,
                phi1=phi1, phi2=phi2)


# --- phase labelling and block averaging ---------------------------------

def calc_dG(c_dil, e_dil, c_den, e_den, ndraws=10000, dGmin=-10, rng=None):
    if rng is None:
        rng = np.random

    if not (np.isfinite(c_dil) and np.isfinite(c_den)) or c_dil <= 0 or c_den <= 0:
        return np.nan, np.nan

    dG = np.log(c_dil / c_den)

    if not (np.isfinite(e_dil) and np.isfinite(e_den)):
        dG_error = np.nan
    else:
        spread_dil = rng.normal(c_dil, e_dil, size=ndraws)
        spread_den = rng.normal(c_den, e_den, size=ndraws)
        keep = (spread_dil > 0) & (spread_den > 0)
        if not keep.any():
            dG_error = np.nan
        else:
            if keep.mean() < 0.9:
                print(f"  ! {100 * (1 - keep.mean()):.0f}% of dG draws rejected "
                      f"as non-positive; error bar is unreliable")
            dG_error = float(np.std(np.log(spread_dil[keep] / spread_den[keep])))

    if dG < dGmin:
        dG = dGmin
        dG_error = 0.
        print(f'dG extremely small, setting dG to {dGmin:.1f}')
    return dG, dG_error


def masked_series(counts_ts, mask, chunk=25):
    """Per-frame mean residue count over a fixed set of voxels."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return None
    n = counts_ts.shape[0]
    out = np.empty(n, dtype=np.float64)
    for i in range(0, n, chunk):
        out[i:i + chunk] = counts_ts[i:i + chunk][:, idx].mean(axis=1)
    return out


def phase_masks(s, dil, den):
    # dil_lo/dil_hi come back exponentiated from the log KDE, den_lo/den_hi
    # are linear. Both are concentrations, so the comparison below is valid.
    c_dil, dil_lo, dil_hi = dil
    c_den, den_lo, den_hi = den

    if (np.isfinite(c_dil) and np.isfinite(c_den)
            and np.isfinite(dil_hi) and np.isfinite(den_lo)
            and dil_hi > den_lo):
        cut = float(np.sqrt(max(c_dil, 0.0) * max(c_den, 0.0)))
        dil_hi, den_lo = min(dil_hi, cut), max(den_lo, cut)

    m_dil = np.zeros(s.shape, dtype=bool)
    m_den = np.zeros(s.shape, dtype=bool)
    if np.isfinite(dil_lo) and np.isfinite(dil_hi):
        m_dil = (s > 0) & (s >= dil_lo) & (s <= dil_hi)
    if np.isfinite(den_lo) and np.isfinite(den_hi):
        m_den = (s >= den_lo) & (s <= den_hi) & ~m_dil
    return m_dil, m_den


def block_sem(series):
    """Standard error of the mean of a correlated time series."""
    if series is None:
        return np.nan, np.nan
    s = np.asarray(series, float).ravel()
    s = s[np.isfinite(s)]
    if s.size < 2:
        return np.nan, np.nan
    block = BlockAnalysis(s)
    block.SEM()
    return float(block.sem), float(block.bs)


# --- output --------------------------------------------------------------

def process_sample(pdb, dcds, name='IDR', start=DEFAULT_START, end=None,
                   chunk_size=CHUNK_SIZE, out=None, names=None):
    """Domain decomposition of one simulation at voxel side length S."""
    dcds = [str(d) for d in (dcds if isinstance(dcds, (list, tuple)) else [dcds])]
    u = mda.Universe(str(pdb), dcds if len(dcds) > 1 else dcds[0])
    total_frames = len(u.trajectory)

    ag1, ag2, n_species = detect_species(u)
    binary = n_species == 2
    len1 = len(ag1.segments[0].atoms) if len(ag1) else 0
    len2 = len(ag2.segments[0].atoms) if len(ag2) else 0
    n1, n2 = len(ag1.segments), len(ag2.segments)
    sp1, sp2 = (names if names else (['species1', 'species2'] if binary
                                     else [name, '']))

    print(f"{name}: {total_frames} frames, {n_species} species, "
          f"starting at {start}")
    if binary:
        print(f"  {sp1}: {n1} x {len1} res | {sp2}: {n2} x {len2} res")

    data = analyze_3d_chunks(u, ag1, ag2, start=start, end=end,
                             chunk_size=chunk_size)
    avg_tot, avg_s1, avg_s2 = data['avg_tot'], data['avg_s1'], data['avg_s2']
    phi1, phi2, n_frames = data['phi1'], data['phi2'], data['n_frames']

    vox_vol = float(np.prod(data['chunk_dims']))
    conv = ATOMS_A3_TO_mM / vox_vol / len1 if len1 else np.nan
    scale = 1.0 if binary else conv
    units = 'res/voxel' if binary else 'mM'

    dil, den = phase_concentrations(avg_tot)
    dil1, den1 = phase_concentrations(avg_s1) if binary else ((np.nan,) * 3,) * 2
    dil2, den2 = phase_concentrations(avg_s2) if binary else ((np.nan,) * 3,) * 2
    cdil_onephase = float(avg_tot.mean()) * scale

    m_dil, m_den = phase_masks(avg_tot, dil, den)
    cdil_sem, _ = block_sem(masked_series(data['tot_ts'], m_dil))
    cden_sem, _ = block_sem(masked_series(data['tot_ts'], m_den))
    cdil_sem, cden_sem = cdil_sem * scale, cden_sem * scale

    cdil, cden = dil[0] * scale, den[0] * scale
    dG, dG_err = calc_dG(cdil, cdil_sem, cden, cden_sem)

    frac_vox_ge_cden = dense_fraction(avg_tot, den[0])

    # --- demixing: value from the averaged map, error by blocking
    if binary:
        D = demixing(avg_s1, avg_s2, phi1, phi2)
        D_err, _ = block_sem(data['d_series'])
    else:
        D = D_err = np.nan

    # --- report
    print(f"Dilute: {cdil:.4g} +/- {cdil_sem:.2g} {units}, "
          f"Dense: {cden:.4g} +/- {cden_sem:.2g} {units}, "
          f"dG: {dG:.3f} +/- {dG_err:.3f} kT")
    print(f"  windows: dilute [{dil[1] * scale:.4g}, {dil[2] * scale:.4g}] "
          f"({int(m_dil.sum())} voxels), "
          f"dense [{den[1] * scale:.4g}, {den[2] * scale:.4g}] "
          f"({int(m_den.sum())} voxels)")
    print(f"  frac voxels >= cden: {frac_vox_ge_cden:.3f}")
    if binary:
        print(f"  {sp1}: dil {dil1[0]:.4g} den {den1[0]:.4g} | "
              f"{sp2}: dil {dil2[0]:.4g} den {den2[0]:.4g} {units}")
        print(f"  D = {D:.3f} +/- {D_err:.3f}")
    if not np.isfinite(den[0]):
        print(f"One-phase concentration: {cdil_onephase:.4g} {units}")

    save_path = Path(out) if out else Path(f"{name}_data.npz")
    save_path.parent.mkdir(parents=True, exist_ok=True)

    fields = dict(
        # voxel maps, residues/voxel
        avg_tot=avg_tot.astype(np.float32),
        chunk_dims=data['chunk_dims'], n_chunks=data['n_chunks'],
        chunk_size=chunk_size, conv=conv, conc_units=units,
        # concentrations, in `conc_units`
        cdil=cdil, cden=cden, cdil_onephase=cdil_onephase,
        cdil_sem=cdil_sem, cden_sem=cden_sem,
        dG=dG, dG_err=dG_err,
        frac_vox_ge_cden=frac_vox_ge_cden,
        # system
        base_name=name, n_species=n_species,
        n_frames=n_frames, total_frames=total_frames, start_frame=start)

    if binary:
        fields.update(
            avg_s1=avg_s1.astype(np.float32), avg_s2=avg_s2.astype(np.float32),
            cdil1=dil1[0], cden1=den1[0], cdil2=dil2[0], cden2=den2[0],
            D_norm=D, D_err=D_err,
            dGij=dG,                       # alias, for the mixture scripts
            seq_name1=sp1, seq_name2=sp2,
            len1=len1, len2=len2, n1=n1, n2=n2,
            comp1=int(round(phi1 * 100)), comp2=int(round(phi2 * 100)),
            phi1=phi1, phi2=phi2)
    else:
        fields.update(
            active_dens=(avg_tot * conv).astype(np.float32),  # mM, for Figure1
            seq_len=len1)

    np.savez_compressed(save_path, **fields)
    print(f"Wrote {save_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Domain-decomposition analysis of a slab (phase-separation) "
            "simulation. Centres the slab each frame, bins residues into 3D "
            "voxels, and estimates dilute- and dense-phase concentrations "
            "(with block-averaged errors) and the transfer free energy dG. "
            "One or two species are detected from the topology; for two, the "
            "demixing index D is also reported."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("pdb", help="Topology file (.pdb) matching the trajectory.")
    parser.add_argument("dcd", nargs="+",
                        help="Trajectory file(s) (.dcd), in concatenation order.")
    parser.add_argument("--name", default=None,
                        help="Sample label for console output and the default "
                             "output filename. Defaults to the .pdb stem.")
    parser.add_argument("--out", default=None,
                        help="Output path. Defaults to '<name>_data.npz'.")
    parser.add_argument("--names", nargs=2, metavar=('SEQ1', 'SEQ2'),
                        default=None,
                        help="Labels for the two species in a mixture, in the "
                             "order they first appear in the topology.")
    parser.add_argument("--start", type=int, default=DEFAULT_START,
                        help="Number of initial frames discarded as equilibration.")
    parser.add_argument("--end", type=int, default=None,
                        help="Frame index (exclusive) at which to stop. "
                             "By default all frames after 'start' are used.")
    parser.add_argument("--voxel-size", type=float, default=CHUNK_SIZE,
                        help="Voxel side length S, in Angstrom.")

    args = parser.parse_args()
    process_sample(args.pdb, args.dcd,
                   name=args.name or Path(args.pdb).stem,
                   start=args.start, end=args.end,
                   chunk_size=args.voxel_size, out=args.out, names=args.names)