"""
phospho_rauh.py — does the NMR result hold with the published CALVADOS phosphoresidue model (Rauh, Hedemark, Tesei &
Lindorff-Larsen, Biophys. J. 2026, 125, 396) in place of the simulated phosphate? Exploratory robustness check,
requested 2026-10-27 after all NMR results were known. Design fixed here before any reweighted value was computed.

Simulated phosphate (residues_phospho_ptyr.csv): the parent bead (lambda, sigma of S, T, Y) with q = -2 at every pH.
Published model (CALVADOS package examples/single_pIDR/input/residues_pCALVADOS2.csv, read 2026-10-27):
  pSer  lambda 0.0925  sigma 0.601 nm  pKa 6.01
  pThr  lambda 0.0013  sigma 0.635 nm  pKa 6.30
  q = -1 - 1 / (1 + 10^(pKa - pH))  (examples/single_pIDR/prepare.py), evaluated at each sample's pH.
pTyr is not in the package. The study repository (KULL-Centre/_2025_rauh_phosphorylation, residue_files/residues.csv)
lists pTyr as lambda 0.977 (unchanged), sigma 0.703, pKa 5.83; used here for the three pTyr pairs, reported separately.

Reweighting, no new simulations: each frame x of the three phosphorylated (P) CALVADOS-2 runs of a pair gets the weight
exp(-beta dU(x)), with dU the exact energy change of the site beads from the simulated to the target parameters
(Ashbaugh-Hatch with the CALVADOS-2 shift and 2 nm cutoff, Debye-Hueckel with the 4 nm cutoff, all pairs |i - j| >= 2,
site-site pairs included once), at each sample's T, ionic strength and pH. Moving to a larger bead is the safe direction
(the target excludes configurations the simulation samples, not the reverse). Two targets:
  charge  pH-dependent charge only (parent lambda and sigma; the 'charge only' baseline of Rauh et al.)
  full    the published model (charge, lambda and sigma)
The unmodified runs (W) are unchanged. E_i is the noise-corrected rewiring of phospho_nmr_test with the P contact
frequencies replaced by their reweighted values; rho_p is its partial Spearman given log distance (csp_table and
partial of phospho_nmr_test, window 3). With uniform weights E_i must reproduce the reported rho_p (check).
Reported per pair and target: rho_p, minimum effective sample fraction over the three P runs (a pair is flagged below
0.10), and the reweighted change in mean Rg (P against W, %). Over pairs: number positive, mean, one-sided Wilcoxon,
for the seven pSer/pThr pairs, the three pTyr pairs and all ten.

  BENDER_REANALYSIS_DIR=/Users/taseef/Downloads/Scripts/bender_reanalysis ~/.venvs/finches/bin/python phospho_rauh.py
      ->  summary.txt, per_pair.csv
"""
import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
os.environ.setdefault('BENDER_REANALYSIS_DIR', '/Users/taseef/Downloads/Scripts/bender_reanalysis')
sys.path[:0] = ['/Volumes/User Homes/Nu Project/BENDER-MT/code', EV, os.path.join(EV, 'phospho_nmr_amend2'),
                os.path.join(EV, 'phospho_nmr_window')]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

import phospho_nmr_test as T  # noqa: E402
import phospho_nmr_amend2_test as A2  # noqa: E402
import phospho_nmr_window as W  # noqa: E402
import salt_ph_phospho as sp  # noqa: E402

fdt = sp.fdt
RES = pd.read_csv('/Users/taseef/Downloads/Scripts/CALVADOS IDP/input_nmr2/residues_phospho_ptyr.csv').set_index('one')
PHOS = {'S': 'B', 'T': 'J', 'Y': 'Z'}
RAUH = {'S': (0.0925, 0.601, 6.01), 'T': (0.0013, 0.635, 6.30), 'Y': (0.977461144934, 0.703, 5.83)}
TARGETS = ('charge', 'full')
ESS_MIN = 0.10


def ah(r, s, lam):
    """CALVADOS-2 Ashbaugh-Hatch pair energy (kJ/mol), shifted at and zero beyond RC_AH."""
    x6 = (s / r) ** 6
    lj = x6 * x6 - x6
    sh = (s / fdt.RC_AH) ** 12 - (s / fdt.RC_AH) ** 6
    e = fdt.EPS_AH
    v = np.where(r < 2.0 ** (1.0 / 6.0) * s, 4 * e * (lj - lam * sh) + e * (1 - lam), 4 * e * lam * (lj - sh))
    return np.where(r < fdt.RC_AH, v, 0.0)


def rg(xyz):
    c = xyz - xyz.mean(1, keepdims=True)
    return np.sqrt((c ** 2).sum(2).mean(1))


def p_terms(job):
    """Per-frame dU for each target, the contact stack and Rg of one phosphorylated run."""
    name, dcd, top, wt, pseq, sites, T_K, ionic, ph = job
    xyz = T.load(dcd, top)
    Cst = T.contact_stack(xyz)
    q, lam, sig, _ = fdt.calvados_residue_params(pseq, RES, ph=ph)
    P = fdt.get_dielectric_params(T_K, ionic)
    k, eps, rc = P['kappa'], P['eps_yu'], fdt.RC_YU
    dh = lambda r: eps * np.where(r < rc, np.exp(-k * r) / r - np.exp(-k * rc) / rc, 0.0)  # noqa: E731
    n, S = len(pseq), set(sites)
    new = {}
    for t in TARGETS:
        q1, l1, s1 = q.copy(), lam.copy(), sig.copy()
        for s in sites:
            l_, s_, pka = RAUH[wt[s]]
            q1[s] = -1.0 - 1.0 / (1.0 + 10.0 ** (pka - ph))
            if t == 'full':
                l1[s], s1[s] = l_, s_
        new[t] = (q1, l1, s1)
    dU = {t: np.zeros(len(xyz)) for t in TARGETS}
    for s in sites:
        js = np.array([j for j in range(n) if j not in S and abs(j - s) >= 2])
        r = np.linalg.norm(xyz[:, js] - xyz[:, [s]], axis=2)
        u0 = ah(r, 0.5 * (sig[s] + sig[js])[None], 0.5 * (lam[s] + lam[js])[None]) + q[s] * q[js][None] * dh(r)
        for t, (q1, l1, s1) in new.items():
            u1 = ah(r, 0.5 * (s1[s] + sig[js])[None], 0.5 * (l1[s] + lam[js])[None]) + q1[s] * q[js][None] * dh(r)
            dU[t] += (u1 - u0).sum(1)
    for ia, a in enumerate(sites):
        for b in sites[ia + 1:]:
            if abs(a - b) >= 2:
                r = np.linalg.norm(xyz[:, a] - xyz[:, b], axis=1)
                u0 = ah(r, 0.5 * (sig[a] + sig[b]), 0.5 * (lam[a] + lam[b])) + q[a] * q[b] * dh(r)
                for t, (q1, l1, s1) in new.items():
                    u1 = ah(r, 0.5 * (s1[a] + s1[b]), 0.5 * (l1[a] + l1[b])) + q1[a] * q1[b] * dh(r)
                    dU[t] += u1 - u0
    return name, dU, Cst, rg(xyz), P['beta'], {t: [float(new[t][0][s]) for s in sites] for t in TARGETS}


def w_terms(job):
    name, dcd, top = job
    xyz = T.load(dcd, top)
    return name, T.contact_stack(xyz).mean(0), float(rg(xyz).mean())


def main():
    T.EXCL = 3
    sets = {'orig': (list(csv.DictReader(open(T.PAIRS))),
                     pd.read_csv(T.SHIFTS, usecols=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])),
            'new': (list(csv.DictReader(open(A2.PAIRS))), pd.read_csv(A2.SHIFTS))}
    pjobs, wjobs, meta = [], [], {}
    for s in ('orig', 'new'):
        root, prefix = W.RUNS[(s, 'cal')]
        pairs, cs = sets[s]
        for p in pairs:
            sites = [int(x[1:]) - 1 for x in p['sites'].split()]
            wt = p['wt_seq']
            pseq = ''.join(PHOS[c] if i in sites else c for i, c in enumerate(wt))
            amide, _ = T.csp_table(p, cs)
            kind = 'pTyr' if all(wt[i] == 'Y' for i in sites) else 'pSer/pThr'
            meta[p['protein']] = (prefix, amide, kind, p['pH'])
            for k in (1, 2, 3):
                for form in ('P', 'WT'):
                    name = f'{prefix}-{p["protein"]}-{form}-r{k}'
                    d = os.path.join(root, f'{p["T_K"]}K', name)
                    if form == 'P':
                        pjobs.append((name, os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), wt, pseq, sites,
                                      float(p['T_K']), float(p['ionic_M']), float(p['pH'])))
                    else:
                        wjobs.append((name, os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb')))
    with Pool(8) as pool:
        wout = {r[0]: r[1:] for r in pool.map(w_terms, wjobs)}
        pout = {r[0]: r[1:] for r in pool.map(p_terms, pjobs)}
    obs = {}
    for f in ('phospho_nmr/per_protein.csv', 'phospho_nmr_amend2/results_cal/per_protein.csv'):
        for r in csv.DictReader(open(os.path.join(EV, f))):
            if r['analysed'] == 'True':
                obs[r['protein']] = float(r['rho_p'])
    rows = []
    for prot, (prefix, amide, kind, ph) in meta.items():
        idx, dd, csp = amide.pos.to_numpy() - 1, amide.d.to_numpy(), amide.CSP.to_numpy()
        Wr = [wout[f'{prefix}-{prot}-WT-r{k}'] for k in (1, 2, 3)]
        Wm = [w[0] for w in Wr]
        noise = np.mean([np.abs(Wm[a] - Wm[b]).sum(1) for a, b in ((0, 1), (0, 2), (1, 2))], 0)
        rg_w = np.mean([w[1] for w in Wr])
        Ps = [pout[f'{prefix}-{prot}-P-r{k}'] for k in (1, 2, 3)]
        row = {'protein': prot, 'kind': kind, 'pH': float(ph), 'n_res': len(amide), 'rho_p_reported': obs.get(prot)}
        for t in ('simulated',) + TARGETS:
            Pm, ess, rgs = [], [], []
            for dU, Cst, rgf, beta, qs in Ps:
                x = np.zeros(len(rgf)) if t == 'simulated' else -beta * dU[t]
                w = np.exp(x - x.max())
                w /= w.sum()
                ess.append(float(1 / (w ** 2).sum() / len(w)))
                Pm.append(np.tensordot(w, Cst.astype(np.float32), axes=1))
                rgs.append(float(w @ rgf))
            E = np.mean([np.abs(p_ - w_).sum(1) for p_ in Pm for w_ in Wm], 0) - noise
            row[f'rho_p_{t}'] = T.partial(E[idx], csp, dd)[0]
            row[f'ess_min_{t}'] = min(ess)
            row[f'dRg_pct_{t}'] = (np.mean(rgs) / rg_w - 1) * 100
        row['q_site'] = ' '.join(f'{v:.2f}' for v in Ps[0][4]['full'])
        rows.append(row)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(HERE, 'per_pair.csv'), index=False)
    L = ['NMR test with the published CALVADOS phosphoresidue model (Rauh et al. 2026), by reweighting the P runs '
         '(exploratory)', '',
         'check, uniform weights against the reported rho_p: ' + '  '.join(
             f'{r.protein} {r.rho_p_simulated:+.3f}/{r.rho_p_reported:+.3f}' for _, r in R.iterrows()), '',
         'per pair: rho_p (min ESS fraction, reweighted dRg %) for simulated | charge only | published model; site q']
    for _, r in R.iterrows():
        L.append(f'  {r.protein:9s} {r.kind:9s} pH {r.pH:4.2f}  ' + ' | '.join(
            f'{r[f"rho_p_{t}"]:+.2f} ({r[f"ess_min_{t}"]:.2f}, {r[f"dRg_pct_{t}"]:+.1f}%)'
            for t in ('simulated',) + TARGETS) + f'   q {r.q_site}'
            + ('   LOW ESS' if r.ess_min_full < ESS_MIN else ''))
    L.append('')
    for label, g in (('pSer/pThr (7)', R[R.kind == 'pSer/pThr']), ('pTyr (3)', R[R.kind == 'pTyr']), ('all (10)', R)):
        for t in ('simulated',) + TARGETS:
            v = g[f'rho_p_{t}']
            pw = wilcoxon(v, alternative='greater').pvalue if len(v) > 1 else np.nan
            L.append(f'  {label:14s} {t:9s}: {int((v > 0).sum())}/{len(v)} positive, mean {v.mean():+.3f}, '
                     f'median {v.median():+.3f}, one-sided Wilcoxon P {pw:.3g}')
        L.append('')
    txt = '\n'.join(L)
    open(os.path.join(HERE, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
