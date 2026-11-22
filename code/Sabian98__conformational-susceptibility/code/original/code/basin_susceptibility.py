"""
Do IDPs hop between conformational basins, and does heat move chains between basins or reshape them?

Landscape per chain (300 K): per-frame log distances between the centroids of S chain segments (S = N // 8, 6..16),
standardised, projected on the first two principal components (whitened). Free energy F = -ln P on an 80 x 80 grid,
P = histogram smoothed with a Gaussian of width 0.5 n^(-1/6) (half Scott's rule; n_eff = frames / integrated
autocorrelation time of PC1 is reported, but smoothing with it erased even well-separated synthetic basins).
Basins: minima of F separated by a barrier of at least DF_MIN = 0.5 kT (persistence clustering), each holding at
least W_MIN of the frames (smaller ones merge across their lowest saddle). A chain is multi-basin only if both halves
of its trajectory also give at least two basins. Calibration on synthetic two-state switching data (30 features, AR(1)
noise): no false basins on single-well data (0/10); two wells about 5 within-SD apart in the projection found in 6 to
9 of 10 runs (2 to 8 of 10 with the both-halves rule). The detector is conservative: a single basin can hide wells
shallower than this.

Hops: changes of basin label along the trajectory after removing visits shorter than MIN_RUN frames.

Susceptibility split (exact, law of total covariance), A = Rg, G = d(beta H)/dT per frame:
  chi = -Cov(A, G) = between + within
  between = -sum_k p_k (<A>_k - <A>)(<G>_k - <G>)   heat moving chains between basins
  within  = -sum_k p_k Cov_k(A, G)                   heat reshaping each basin
Basin populations: dp_k/dT = -Cov(1_k, G) = -p_k (<G>_k - <G>).
Ceiling for the between share: the same split with k equal-population slices of Rg itself.
Validation: frames at 278, 288, 310 and 320 K are projected with the 300 K standardisation and PCA and assigned on
the 300 K grid; the measured change p_k(320) - p_k(278) is compared with dp_k/dT integrated over the five
temperatures (trapezoid).

G: CALVADOS-2 from the FDT npz (G_T); Mpipi-GG -PotEng/(R T^2) (its Hamiltonian does not depend on T).
Sample: the proteins simulated in both models (common set), the same chains in each model.
Output: Traj Level Analysis/Analysis_FDT/basin_susceptibility/{per_protein,per_basin,summary}.csv
Usage: python3 basin_susceptibility.py [--workers 8] [--limit N]
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import ndimage  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import TODO_T2_trajectory_loader as L  # noqa: E402
import mpipi_temperature_hpc as mh  # noqa: E402

bc = L.bc
TA = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
OUT_DIR = os.path.join(L.OUT_ROOT, 'basin_susceptibility')
PER = os.path.join(L.OUT_ROOT, 'energy_mechanism_split', 'per_protein.csv')
MP_DIR = os.path.join(TA, 'Analysis_FDT_MPIPI')
MP_RUNS = {T: f'/Volumes/User Homes/Nu Project/BENDER-BIO/MPIPI/mpipi_sim/{T}k_runs' for T in (278, 288, 310, 320)}
MP_RUNS[300] = '/Volumes/User Homes/Nu Project/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k '
MP_COMPACT_RUNS = '/Volumes/User Homes/Nu Project/BENDER-MT/mpipi_compact'
TEMPS = [278, 288, 300, 310, 320]
R_KCAL = 1.98720425864083e-3
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
GRID, DF_MIN, W_MIN, MIN_RUN = 80, 0.5, 0.05, 3
BW_FACTOR, BW_FROM_NEFF = 0.5, False  # smoothing width = BW_FACTOR * n^(-1/6), n = frames (n_eff over-smooths)


# ---------------------------------------------------------------- frames, G and Rg per temperature

def _rg(xyz):
    return np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))


def _check_rg(xyz, rg, tag):
    r = _rg(xyz)
    if len(r) != len(rg) or np.corrcoef(r, rg)[0, 1] < 0.999 or np.abs(r - rg).max() > 0.02 * rg.mean():
        raise ValueError(f'{tag}: trajectory frames do not match the npz (Rg check)')


def frames_calvados(uid, T, hint):
    base = L.TRAJ_ROOT_FMT.format(T=T)
    d = L.find_traj_dir(uid, base, hint, flat=None if T == 300 else True)
    if d is None:
        return None
    xyz, _, _ = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
    p = os.path.join(L.OUT_ROOT, 'npz', f'{T}K', f'{uid}.npz')
    if not os.path.exists(p):
        return None
    z = np.load(p)
    G, rg = z['G_T'].astype(float), z['rg'].astype(float)
    _check_rg(xyz, rg, f'CALVADOS-2 {uid} {T} K')
    return xyz, G, rg


def frames_mpipi(uid, T, compact=False):
    d = os.path.join(MP_COMPACT_RUNS if compact else MP_RUNS[T], uid)
    npz = os.path.join(MP_DIR, 'compact', 'npz_300K', f'{uid}.npz') if compact else os.path.join(MP_DIR, f'npz_{T}K', f'{uid}.npz')
    if not (os.path.isfile(os.path.join(d, 'traj.xtc')) and os.path.exists(npz)):
        return None
    seq = bc.read_mpipi_sequence(os.path.join(d, 'seq.dat'))
    xyz, _ = bc.load_mpipi(os.path.join(d, 'traj.xtc'), n_res=len(seq))
    steps = mh.xtc_steps(os.path.join(d, 'traj.xtc'))
    z = np.load(npz)
    pos = {int(s): i for i, s in enumerate(steps)}
    idx = np.array([pos[int(s)] for s in z['step']])
    xyz = np.asarray(xyz, dtype=float)[idx]
    G = -z['PotEng_kcal'].astype(float) / (R_KCAL * T ** 2)
    rg = z['rg'].astype(float)
    _check_rg(xyz, rg, f'Mpipi-GG {uid} {T} K')
    return xyz, G, rg


# ---------------------------------------------------------------- landscape and basins

def seg_features(xyz):
    n = xyz.shape[1]
    S = int(np.clip(n // 8, 6, 16))
    e = np.linspace(0, n, S + 1).astype(int)
    c = np.stack([xyz[:, a:b].mean(1) for a, b in zip(e[:-1], e[1:])], 1)
    i, j = np.triu_indices(S, 1)
    return np.log(np.linalg.norm(c[:, i] - c[:, j], axis=2))


def tau_int(x):
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    tau = 1.0
    for t in range(1, n // 5):
        tau += 2 * acf[t]
        if t >= 5 * tau:
            break
    return max(1.0, tau)


def _neighbour_pairs(c, F):
    """(c_a, c_b, F_a, F_b) for right, down, down-right and down-left neighbours."""
    for sa, sb in (((slice(None), slice(None, -1)), (slice(None), slice(1, None))),
                   ((slice(None, -1), slice(None)), (slice(1, None), slice(None))),
                   ((slice(None, -1), slice(None, -1)), (slice(1, None), slice(1, None))),
                   ((slice(None, -1), slice(1, None)), (slice(1, None), slice(None, -1)))):
        yield c[sa], c[sb], F[sa], F[sb]


class Landscape:
    """Projection and grid fitted at 300 K; basins by persistence of free-energy minima."""

    def __init__(self, X):
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-9
        Z = (X - self.mu) / self.sd
        _, s, vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
        self.V = vt[:2].T
        self.scale = s[:2] / np.sqrt(len(Z))
        self.var2 = float((s[:2] ** 2).sum() / (s ** 2).sum())
        Y = self.project(X)
        self.lo, self.hi = Y.min(0) - 1.0, Y.max(0) + 1.0
        self.cell = (self.hi - self.lo) / GRID

    def project(self, X):
        return ((X - self.mu) / self.sd) @ self.V / self.scale

    def cells(self, Y):
        ij = np.floor((Y - self.lo) / self.cell).astype(int)
        return np.clip(ij, 0, GRID - 1)

    def basins(self, Y):
        """Frame labels (0..k-1) and the cell -> basin map for frames Y."""
        n = len(Y)
        neff = max(10.0, n / tau_int(Y[:, 0]))
        ij = self.cells(Y)
        H = np.zeros((GRID, GRID))
        np.add.at(H, (ij[:, 0], ij[:, 1]), 1.0)
        h = BW_FACTOR * (neff if BW_FROM_NEFF else n) ** (-1 / 6)
        P = ndimage.gaussian_filter(H, sigma=h / self.cell, mode='constant')
        P /= P.sum()
        active = P >= 1e-3 * P.max()
        F = np.where(active, -np.log(np.where(active, P, 1.0)), np.inf)
        comp = self._persistence(F, active)
        # cells outside the active region take the basin of the nearest active cell
        _, (ii, jj) = ndimage.distance_transform_edt(~active, return_indices=True)
        comp = comp[ii, jj]
        lab = comp[ij[:, 0], ij[:, 1]]
        comp, lab = self._merge_small(F, active, comp, lab)
        return lab, comp, neff

    @staticmethod
    def _persistence(F, active):
        parent = -np.ones(GRID * GRID, dtype=int)
        rootF = {}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        flat = F.ravel()
        for c in np.argsort(flat):
            if not np.isfinite(flat[c]):
                break
            i, j = divmod(c, GRID)
            roots = set()
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    a, b = i + di, j + dj
                    if (di or dj) and 0 <= a < GRID and 0 <= b < GRID and parent[a * GRID + b] >= 0:
                        roots.add(find(a * GRID + b))
            if not roots:
                parent[c] = c
                rootF[c] = flat[c]
                continue
            best = min(roots, key=lambda r: rootF[r])
            parent[c] = best
            for r in roots:
                if r != best and flat[c] - rootF[r] < DF_MIN:      # shallow basin: merge across this saddle
                    parent[r] = best
        comp = np.full(GRID * GRID, -1)
        for c in range(GRID * GRID):
            if parent[c] >= 0:
                comp[c] = find(c)
        ids = {r: k for k, r in enumerate(sorted(set(comp[comp >= 0])))}
        return np.array([ids.get(c, -1) for c in comp]).reshape(GRID, GRID)

    @staticmethod
    def _merge_small(F, active, comp, lab):
        while True:
            k = comp.max() + 1
            w = np.bincount(lab, minlength=k) / len(lab)
            small = [b for b in range(k) if w[b] < W_MIN and k > 1]
            if not small:
                break
            b = min(small, key=lambda x: w[x])
            # lowest saddle to another basin across neighbouring cells (4 directions cover all 8 neighbours)
            best, bestF = None, np.inf
            for A, B, FA, FB in _neighbour_pairs(comp, F):
                m = ((A == b) & (B != b)) | ((B == b) & (A != b))
                if m.any():
                    sF = np.maximum(FA, FB)[m]
                    other = np.where(A[m] == b, B[m], A[m])
                    i = int(np.argmin(sF))
                    if sF[i] < bestF:
                        best, bestF = int(other[i]), float(sF[i])
            if best is None:          # not adjacent: merge into the most populated basin
                best = int(np.argmax(np.where(np.arange(k) == b, -1, w)))
            comp = np.where(comp == b, best, comp)
            lab = np.where(lab == b, best, lab)
            ids = {v: i for i, v in enumerate(sorted(set(comp.ravel())))}
            comp = np.vectorize(ids.get)(comp)
            lab = np.vectorize(ids.get)(lab)
        return comp, lab


def hops(lab):
    """Basin changes after folding visits shorter than MIN_RUN frames into the previous visit."""
    runs = []
    for x in lab:
        if runs and runs[-1][0] == x:
            runs[-1][1] += 1
        else:
            runs.append([x, 1])
    clean = []
    for x, n in runs:
        if clean and (n < MIN_RUN or clean[-1][0] == x):
            clean[-1][1] += n
        else:
            clean.append([x, n])
    return len(clean) - 1, float(np.mean([n for _, n in clean]))


def split(A, G, lab, k):
    """chi = -Cov(A, G) (population covariance) split into between- and within-basin parts; dp_k/dT."""
    a, g = A - A.mean(), G - G.mean()
    total = -float((a * g).mean())
    p = np.bincount(lab, minlength=k) / len(lab)
    mA = np.bincount(lab, weights=a, minlength=k) / np.maximum(p * len(lab), 1)
    mG = np.bincount(lab, weights=g, minlength=k) / np.maximum(p * len(lab), 1)
    between = -float((p * mA * mG).sum())
    return total, between, p, -p * mG


# ---------------------------------------------------------------- one chain

def analyse(task):
    uid, model, nu, hint, compact = task
    try:
        temps = [300] if compact else TEMPS
        data = {}
        for T in temps:
            f = frames_calvados(uid, T, hint) if model == 'CALVADOS-2' else frames_mpipi(uid, T, compact)
            if f is not None:
                data[T] = f
        if 300 not in data:
            return None, None, f'{model} {uid}: no 300 K frames'
        xyz, G, rg = data[300]
        ls = Landscape(seg_features(xyz))
        Y = ls.project(seg_features(xyz))
        lab, comp, neff = ls.basins(Y)
        k = int(lab.max()) + 1
        h = len(Y) // 2
        k1 = int(ls.basins(Y[:h])[0].max()) + 1
        k2 = int(ls.basins(Y[h:])[0].max()) + 1
        nh, dwell = hops(lab)
        total, between, p, dpdT = split(rg, G, lab, k)
        # ceiling: k equal-population slices of Rg itself
        if k > 1:
            q = np.searchsorted(np.quantile(rg, np.linspace(0, 1, k + 1)[1:-1]), rg)
            ceil = split(rg, G, q, k)[1] / total if total else np.nan
        else:
            ceil = np.nan
        row = {'UniProt_ID': uid, 'model': model, 'compact_set': compact, 'nu': nu, 'N': xyz.shape[1],
               'var_pc12': ls.var2, 'n_eff': neff, 'n_basins': k, 'n_basins_half1': k1, 'n_basins_half2': k2,
               'multi_basin': bool(k > 1 and k1 > 1 and k2 > 1), 'hops_per_1000': 1000 * nh / len(lab),
               'mean_dwell_frames': dwell, 'p_max': float(p.max()),
               'chi_rg_rel': 100 * total / rg.mean(), 'between_share': between / total if total else np.nan,
               'between_share_ceiling': ceil}
        basins = []
        per_T = {}
        for T, (x, g, r) in data.items():
            labT = lab if T == 300 else comp[tuple(ls.cells(ls.project(seg_features(x))).T)]
            tot, btw, pT, dpT = split(r, g, labT, k)
            per_T[T] = (pT, dpT)
            row[f'between_share_{T}'] = btw / tot if tot else np.nan
        for b in range(k):
            bb = {'UniProt_ID': uid, 'model': model, 'compact_set': compact, 'basin': b, 'n_basins': k,
                  'multi_basin': row['multi_basin']}
            for T, (pT, dpT) in per_T.items():
                bb[f'p_{T}'] = float(pT[b])
                bb[f'dpdT_{T}'] = float(dpT[b])
            if all(T in per_T for T in TEMPS):
                bb['dp_measured'] = float(per_T[320][0][b] - per_T[278][0][b])
                bb['dp_predicted'] = float(np.trapezoid([per_T[T][1][b] for T in TEMPS], TEMPS))
            basins.append(bb)
        row['n_T'] = len(data)
        return row, basins, None
    except Exception as e:
        return None, None, f'{model} {uid}: {e!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


# ---------------------------------------------------------------- sample, run, summarise

def sample():
    """The common set: proteins with five-temperature susceptibility in both models; both models on the same chains."""
    per = pd.read_csv(PER, index_col=0).dropna(subset=['chi_C', 'chi_M', 'nu_C', 'nu_M'])
    kingdom = pd.read_csv(bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    tasks = []
    for model, m in (('CALVADOS-2', 'C'), ('Mpipi-GG', 'M')):
        tasks += [(u, model, float(per.at[u, f'nu_{m}']), kingdom.get(u), False) for u in per.index]
    return tasks


def main(workers, limit):
    os.makedirs(OUT_DIR, exist_ok=True)
    tasks = sample()
    if limit:
        tasks = tasks[:limit]
    print(f'{len(tasks)} chains, {workers} workers', flush=True)
    rows, basins, fails = [], [], []
    with mp.Pool(workers) as pool:
        for i, (r, b, err) in enumerate(pool.imap_unordered(analyse, tasks), 1):
            if err:
                fails.append(err)
            else:
                rows.append(r)
                basins += b
            if i % 25 == 0:
                print(f'  {i}/{len(tasks)} done, {len(fails)} failed', flush=True)
    P = pd.DataFrame(rows)
    B = pd.DataFrame(basins)
    P['b'] = pd.cut(P.nu, BINS).astype(str)
    P.to_csv(os.path.join(OUT_DIR, 'per_protein.csv'), index=False)
    B.to_csv(os.path.join(OUT_DIR, 'per_basin.csv'), index=False)
    with open(os.path.join(OUT_DIR, 'failed.txt'), 'w') as fh:
        fh.write('\n'.join(fails))
    g = P.groupby(['model', 'b'])
    S = pd.DataFrame({'n': g.size(), 'frac_multi_basin': g.multi_basin.mean(), 'n_basins_median': g.n_basins.median(),
                      'hops_per_1000_median': g.hops_per_1000.median(), 'p_max_median': g.p_max.median(),
                      'chi_rg_rel_median': g.chi_rg_rel.median(),
                      'between_share_median_multi': g.apply(lambda x: x.loc[x.multi_basin, 'between_share'].median(), include_groups=False),
                      'between_ceiling_median_multi': g.apply(lambda x: x.loc[x.multi_basin, 'between_share_ceiling'].median(), include_groups=False)})
    S.to_csv(os.path.join(OUT_DIR, 'summary.csv'))
    bc.write_provenance(os.path.join(OUT_DIR, 'summary.csv'), __file__,
                        {'grid': GRID, 'barrier_kT': DF_MIN, 'min_weight': W_MIN, 'min_run_frames': MIN_RUN,
                         'features': 'log distances between centroids of N//8 (6..16) segments, PCA 2D whitened',
                         'bandwidth': f'{BW_FACTOR} x n^(-1/6)', 'set': 'common to both models', 'temps': TEMPS}, [PER])
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 30)
    print(S.round(3).to_string())
    print(f'failed: {len(fails)}')
    for f in fails[:5]:
        print('  ', f[:300])
    v = B.dropna(subset=['dp_measured'])
    v = v[v.multi_basin]
    for model, x in v.groupby('model'):
        if len(x) > 2:
            r = np.corrcoef(x.dp_predicted, x.dp_measured)[0, 1]
            slope = float((x.dp_predicted * x.dp_measured).sum() / (x.dp_predicted ** 2).sum())
            print(f'validation {model}: {len(x)} basins in multi-basin chains, r = {r:.3f}, measured/predicted slope {slope:.3f}')
    print(f'saved {OUT_DIR}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()
    main(a.workers, a.limit)
