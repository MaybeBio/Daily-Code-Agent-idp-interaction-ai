"""
pre_reweight.py — experiment-informed hub occupancy: reweight conformational ensembles against PRE data
(PeptoneDB-Integrative, $BSS_DATA/Downloads/PeptoneDBs) and compute weighted hub occupancy. Standard methods, nothing new:

  forward model  Calpha-based PRE (spin label and amide proton both placed at Calpha):
                 Gamma2 = K <r^-6> J(tau_c),  K = 1.23e16 A^6 s^-2 (nitroxide), J = 4 tau_c + 3 tau_c / (1 + w_H^2 tau_c^2)
                 one tau_c per ensemble, fitted on the unreweighted prior (grid 1-20 ns)
  data           Gamma2 rates (PRE-<site>.dat) are used directly; intensity ratios (PRE-I-<site>.dat) are converted by
                 inverting I = R2 exp(-Gamma2 t) / (R2 + Gamma2) with R2 = 10 s^-1, t = 10 ms; ratios below 0.05
                 (broadened beyond use) are dropped and ratios above 1 are set to 1; residues within 3 of the label
                 are dropped. Error: max(propagated experimental error, 10% of Gamma2, 1 s^-1).
  reweighting    BME (Bottaro, Bengtsen & Lindorff-Larsen 2020) on the Gamma2 values; theta from a fixed grid, chosen
                 by leave-one-label-out cross-validated chi^2 (ties -> larger theta).
  hubs           per frame: Calpha contact graph (0.8 nm, |i-j| >= 3), exact betweenness, top-10% hub set (hub_frames.py);
                 occupancy = sum over frames of weight x [residue is a hub].
These rules were fixed before any reweighted hub occupancy was looked at (methods_paper/RECOMB_PLAN.md).

  ~/.venvs/bender/bin/python pre_reweight.py validate   # PED unrestrained pools -> do they move toward the experts'
                                                         # PRE-restrained ensembles? (alpha-syn, Sic1, drkN)
  ~/.venvs/bender/bin/python pre_reweight.py apply      # our CALVADOS-2 / MPIPI-GG runs of the PRE constructs + ACTR MD
"""
import math
import os
import pickle
import re
import sys
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import ped_null as pn  # noqa: E402

PEPTONE = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Downloads/PeptoneDBs/PeptoneDB-Integrative'))
CACHE = os.path.join(HERE, 'pre_cache')
K = 1.23e16                      # A^6 s^-2
R2, T_DELAY = 10.0, 0.010
THETAS = [0.5, 1, 2, 5, 10, 20, 50, 100, 200]
TAUS = [x * 1e-9 for x in (1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20)]
MAX_FRAMES = 1000
EXCLUDE = int(os.environ.get('PRE_EXCLUDE', '3'))   # drop residues with |i - s| <= EXCLUDE; pre-declared 3, sensitivity 9
# PeptoneDB 'PRE-<s>.dat' values for Sic1 contain the diamagnetic baseline (about 30 s^-1 far from every label and
# label-independent), so they are not Gamma2 and are left out (found in code review, 2026-10-03)
TOTAL_R2 = {'Sic1'}


# ── experimental data ────────────────────────────────────────────────────────
def info(prot):
    out = {}
    for line in open(os.path.join(PEPTONE, prot, 'info.csv')):
        t = line.strip().split(',')
        if t[0].startswith('PRE'):
            out[t[0]] = float(t[4])                       # spectrometer MHz
    return out


def load_pre(prot, n, offset=0):
    """[(label index, residue index, Gamma2, sigma, MHz)] in the ensemble's 0-based numbering (offset added)."""
    mhz = info(prot)
    rows = []
    for f in sorted(os.listdir(os.path.join(PEPTONE, prot))):
        m = re.match(r'(PRE-(I-)?(\d+))\.dat$', f)
        if not m:
            continue
        name, ratio, s = m.group(1), bool(m.group(2)), int(m.group(3)) - 1 + offset
        for line in open(os.path.join(PEPTONE, prot, f)):
            t = line.split()
            if len(t) < 2 or t[1].lower() == 'nan':
                continue
            i = int(float(t[0])) - 1 + offset
            if not (0 <= i < n) or abs(i - s) <= 3 or not (0 <= s < n):
                continue
            v = float(t[1])
            e = float(t[2]) if len(t) > 2 and t[2].lower() != 'nan' else 0.0
            if ratio:
                if v < 0.05:
                    continue
                v = min(v, 1.0)
                g = solve_gamma(v)
                dI = v * (T_DELAY + 1.0 / (R2 + g))        # |dI/dGamma2|
                e = e / dI if dI > 0 else 0.0
            else:
                g = v
            rows.append((s, i, g, max(e, 0.1 * g, 1.0), mhz.get(name, 800.0)))
    return rows


def solve_gamma(I):
    """Invert I = R2 exp(-G t) / (R2 + G) for G >= 0 (monotone decreasing in G)."""
    lo, hi = 0.0, 5000.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if R2 * math.exp(-mid * T_DELAY) / (R2 + mid) > I:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ── ensembles: per-frame r^-6 to the label sites and per-frame hub sets ────
def frame_features(seq, frames, sites):
    """frames: iterable of (xyz, box). Returns r6[f][site][i] (A^-6) and hub sets per frame."""
    n = len(seq)
    n_top = max(1, round(0.10 * n))
    r6, hubs = [], []
    rng = hf.tie_rng(seq)
    for xyz, box in frames:
        hf.unwrap(xyz, box)
        X, Y, Z = xyz
        r6.append([[0.0 if i == s else ((X[i] - X[s]) ** 2 + (Y[i] - Y[s]) ** 2 + (Z[i] - Z[s]) ** 2) ** -3
                    for i in range(n)] for s in sites])
        bc = hf.betweenness(hf.contact_lists(xyz, n), n)
        hubs.append(sorted(hf.top_hubs(bc, n_top, rng)))
    return np.array(r6, dtype=np.float64), hubs


def frames_ped(path):
    seq, models = pn.load(path)
    if len(models) > MAX_FRAMES:
        models = models[::math.ceil(len(models) / MAX_FRAMES)]
    return seq, [([list(m[0]), list(m[1]), list(m[2])], None) for m in models]


def frames_dcd(d, uid):
    seq = ''.join(hf.AA3.get(line[17:20], 'X') for line in open(os.path.join(d, 'top.pdb'))
                  if line.startswith('ATOM') and line[12:16].strip() == 'CA')
    fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(d, f'{uid}.dcd'))
    start = n_frames // 2
    idx = [start + (k * (n_frames - start)) // MAX_FRAMES for k in range(MAX_FRAMES)]
    return seq, (hf.dcd_frame(fh, off, frame, natom, has_cell, k) for k in idx)


def frames_xtc(d):
    import hub_frames_mpipi as hm
    import xtc_stdlib as xtc
    seq = hm.read_sequence(os.path.join(d, 'seq.dat'))
    data = open(os.path.join(d, 'traj.xtc'), 'rb').read()
    offs = xtc.index(data)
    idx = [(k * len(offs)) // MAX_FRAMES for k in range(MAX_FRAMES)]

    def gen():
        for k in idx:
            natoms, step, t, box, X, Y, Z = xtc.decode(data, offs[k])
            yield [[v * 10.0 for v in X], [v * 10.0 for v in Y], [v * 10.0 for v in Z]], \
                (box[0] * 10.0, box[4] * 10.0, box[8] * 10.0)
    return seq, gen()


def build(job):
    """Compute (and cache) per-frame features for one ensemble."""
    name, kind, src, prot, offset = job
    os.makedirs(CACHE, exist_ok=True)
    cp = os.path.join(CACHE, f'{name}.pkl')
    if os.path.exists(cp):
        return name
    if kind == 'ped':
        seq, frames = frames_ped(src)
    elif kind == 'dcd':
        seq, frames = frames_dcd(*src)
    else:
        seq, frames = frames_xtc(src)
    pre = load_pre(prot, len(seq), offset) if prot else []
    sites = sorted({s for s, *_ in pre})
    r6, hubs = frame_features(seq, frames, sites)
    pickle.dump({'seq': seq, 'pre': pre, 'sites': sites, 'r6': r6, 'hubs': hubs, 'prot': prot, 'offset': offset},
                open(cp, 'wb'))
    return name


# ── forward model, tau_c fit and BME ─────────────────────────────────────────
def J(tau, mhz):
    w = 2 * math.pi * mhz * 1e6
    return 4 * tau + 3 * tau / (1 + (w * tau) ** 2)


def calc_matrix(c, tau):
    """F[f, k] = per-frame Gamma2 for every observable k, and y, sigma, label of k."""
    site_ix = {s: j for j, s in enumerate(c['sites'])}
    F = np.stack([K * J(tau, mhz) * c['r6'][:, site_ix[s], i] for s, i, g, e, mhz in c['pre']], axis=1)
    y = np.array([g for _, _, g, _, _ in c['pre']])
    sd = np.array([e for _, _, _, e, _ in c['pre']])
    lab = np.array([s for s, *_ in c['pre']])
    return F, y, sd, lab


def bme(F, y, sd, theta, w0=None):
    n = F.shape[0]
    w0 = np.full(n, 1.0 / n) if w0 is None else w0
    Fs, ys = F / sd, y / sd

    def dual(lam):
        a = -Fs @ lam
        amax = a.max()
        e = w0 * np.exp(a - amax)
        Z = e.sum()
        w = e / Z
        val = math.log(Z) + amax + lam @ ys + 0.5 * theta * lam @ lam
        grad = ys - w @ Fs + theta * lam
        return val, grad
    r = minimize(dual, np.zeros(F.shape[1]), jac=True, method='L-BFGS-B')
    a = -Fs @ r.x
    w = w0 * np.exp(a - a.max())
    return w / w.sum()


def chi2(F, y, sd, w):
    return float(np.mean(((w @ F - y) / sd) ** 2))


def phi_eff(w):
    n = len(w)
    s = -np.sum(w[w > 0] * np.log(w[w > 0] * n))
    return float(math.exp(s))


def restrict(c):
    """Apply the near-label exclusion window (the cache already drops |i - s| <= 3)."""
    return dict(c, pre=[x for x in c['pre'] if abs(x[1] - x[0]) > EXCLUDE])


def reweight(c):
    """Fit tau_c on the prior, choose theta by leave-one-label-out CV, return weights and diagnostics.
    If no theta improves the held-out chi2 over the prior, the prior (uniform weights) is kept."""
    c = restrict(c)
    n = c['r6'].shape[0]
    u = np.full(n, 1.0 / n)
    best = min(TAUS, key=lambda t: chi2(*calc_matrix(c, t)[:3], u))
    F, y, sd, lab = calc_matrix(c, best)
    labels = sorted(set(lab))
    cv = {}
    for th in THETAS:
        scores = []
        for L in labels:
            tr, te = lab != L, lab == L
            if tr.sum() == 0:
                continue
            w = bme(F[:, tr], y[tr], sd[tr], th)
            scores.append(chi2(F[:, te], y[te], sd[te], w))
        cv[th] = float(np.mean(scores)) if scores else float('nan')
    theta = min(THETAS, key=lambda t: (round(cv[t], 6), -t))
    prior_cv = float(np.mean([chi2(F[:, lab == L], y[lab == L], sd[lab == L], u) for L in labels]))
    if not cv[theta] < prior_cv:                   # reweighting does not help on held-out labels: keep the prior
        theta, w = float('inf'), u
        cv[theta] = prior_cv
    else:
        w = bme(F, y, sd, theta)
    return w, {'tau_ns': best * 1e9, 'theta': theta, 'chi2_prior': chi2(F, y, sd, u), 'chi2_rew': chi2(F, y, sd, w),
               'cv_prior': prior_cv, 'cv_rew': cv[theta], 'phi_eff': phi_eff(w), 'n_obs': len(y),
               'n_labels': len(labels)}


def occupancy(c, w=None):
    n = len(c['seq'])
    w = np.full(len(c['hubs']), 1.0 / len(c['hubs'])) if w is None else w
    occ = np.zeros(n)
    for wf, hs in zip(w, c['hubs']):
        occ[hs] += wf
    return occ


def load(name):
    return pickle.load(open(os.path.join(CACHE, f'{name}.pkl'), 'rb'))


# ── drivers ──────────────────────────────────────────────────────────────────
VALIDATE = [  # protein, PeptoneDB folder, restraint-free pool, PRE-restrained PED ensemble (reference)
    ('α-synuclein', 'aSyn', 'PED00425/pdbfile.pdb', 'PED00024/PED00024e001.pdb'),
    ('Sic1', 'Sic1', 'PED00423/pdbfile.pdb', 'PED00160/PED00160e002.pdb'),
    ('drkN', 'drkN', 'PED00427/pdbfile.pdb', 'PED00157/PED00157e002.pdb')]


def offset_of(seq, prot):
    ref = ''.join(line.strip() for line in open(os.path.join(PEPTONE, prot, 'seq.fasta')) if not line.startswith('>'))
    if seq == ref:
        return 0
    for off in range(-10, 11):
        if all(0 <= i + off < len(seq) and seq[i + off] == ref[i] for i in range(10, min(len(ref), 60))):
            return off
    raise ValueError(f'{prot}: cannot align ensemble sequence to the PRE construct')


def validate():
    jobs = []
    for name, prot, pool, ref in VALIDATE:
        seq, _ = pn.load(pool)
        jobs.append((f'pool_{prot}', 'ped', pool, prot, offset_of(seq, prot)))
        jobs.append((f'ref_{prot}', 'ped', ref, None, 0))
    with Pool(6) as p:
        list(p.imap_unordered(build, jobs))
    lines = ['PRE reweighting of restraint-free PED pools (IDPConformerGenerator) with PeptoneDB PRE: does the reweighted '
             'hub occupancy move toward the experts\' PRE-restrained ensemble?\n']
    lines.append(f'  near-label exclusion |i - s| <= {EXCLUDE}; theta from leave-one-label-out CV, prior kept if CV does not improve')
    for name, prot, pool, ref in VALIDATE:
        if prot in TOTAL_R2:
            lines.append(f'  {name:12s} left out: PeptoneDB values include the diamagnetic baseline (not Gamma2)')
            continue
        c, r = load(f'pool_{prot}'), load(f'ref_{prot}')
        w, d = reweight(c)
        o0, o1, oref = occupancy(c), occupancy(c, w), occupancy(r)
        lines.append(f'  {name:12s} labels {d["n_labels"]:2d}, obs {d["n_obs"]:3d}, tau_c {d["tau_ns"]:.0f} ns, theta {d["theta"]:g}, '
                     f'phi_eff {d["phi_eff"]:.2f} | chi2 prior {d["chi2_prior"]:.2f} -> {d["chi2_rew"]:.2f}; held-out label '
                     f'chi2 {d["cv_prior"]:.2f} -> {d["cv_rew"]:.2f}')
        lines.append(f'  {"":12s} hub occupancy vs restrained ensemble: pool {spearmanr(o0, oref).correlation:+.2f} -> '
                     f'reweighted pool {spearmanr(o1, oref).correlation:+.2f}')
    print('\n'.join(lines))
    tag = '' if EXCLUDE == 3 else f'_ex{EXCLUDE}'
    open(os.path.join(HERE, f'pre_reweight_validation{tag}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    {'validate': validate}[sys.argv[1]]()
