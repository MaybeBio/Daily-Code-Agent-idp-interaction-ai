"""
phospho_nmr_test.py — analysis of phospho_nmr_plan.md (hash 1851cf7c). Written and hashed before any run of the plan was
labelled and before any chemical shift difference was computed. The PAGE4 prediction from the existing BENDER-MT runs
was fixed beforehand (phospho_nmr_page4_prediction.py, 45d347b3; outputs hashed).

Measurement   amide CSP_i = sqrt(dH^2 + (0.14 dN)^2) between the phosphorylated and the unmodified BMRB entry (ZOIDBERG's
              mined table), residues with H and N in both and the same residue type, phosphosites and residues within 3
              of a site excluded, median dH and dN over the analysed residues subtracted first; |dCA| likewise
              (secondary). A protein is analysed if at least 20 residues remain.
Prediction    E_i = mean of 9 r_i(P_a, W_b) - mean of 3 r_i(W_a, W_b), r_i(X, Y) = sum_j |p_ij(X) - p_ij(Y)|, contact
              frequencies p_ij (Calpha < 0.8 nm, |i - j| >= 3) over 1,000 evenly spaced second-half frames
              (bender_common.load_calvados) of each of the 3 + 3 runs at the NMR conditions.
Primary       per protein the partial Spearman rho_p of E_i with CSP_i given log d_i (d_i = distance to the nearest site;
              Spearman of the residuals of the ranks after least squares on the rank of log d_i, with intercept);
              one-sided Wilcoxon signed-rank over the analysed proteins that rho_p > 0; success if P < 0.05.
Secondary     per-protein P from 10,000 random circular shifts of the CSP residuals (numpy default_rng(20261013)); plain
              Spearman of E_i and of -log d_i with CSP_i; |d deg_i| and |d occ_i| in place of E_i; |dCA_i| in place of
              CSP_i; charge-model reweighting of the unmodified runs (all sites, mutual terms between non-bonded sites,
              Debye-Hueckel at each sample's T and ionic strength, His charge at its pH) in place of E_i; the primary
              without H1.0; PAGE4 with BENDER-MT's existing runs at 278 K (single pair r_i, no noise correction; 300 K
              also reported, descriptive).
Added before any shift difference was computed (declared here, secondary): PAGE4 N-motif test. Among the analysed PAGE4
residues, the mean CSP residual (after the rank regression on log d_i) of residues 4-12 against every contiguous block of
9 analysed residues placed circularly along the chain; P = fraction of block positions with a mean at least as large.
The fixed PAGE4 prediction says the N-motif is the most rewired region away from T51.

  ~/.venvs/finches/bin/python phospho_nmr_test.py --page4-existing   (only the PAGE4 secondaries; needs no new runs)
  ~/.venvs/finches/bin/python phospho_nmr_test.py                    (everything; needs all 36 runs)
      ->  phospho_nmr/{summary*.txt, per_residue*.csv, page4_existing.png}
"""
import argparse
import csv
import os
import sys
from multiprocessing import Pool

os.environ.setdefault('BENDER_REANALYSIS_DIR', '/Users/taseef/Downloads/Scripts/bender_reanalysis')
DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import rankdata, spearmanr, wilcoxon  # noqa: E402

import phospho_network as pn  # noqa: E402
from phospho_direct_test import check_frozen  # noqa: E402

SHIFTS = '/Volumes/User Homes/ZOIDBERG-ICML 2027/experimental_data/processed/bmrb_chemical_shifts.csv'
PAIRS = os.path.join(DYN, 'gaivi', 'phospho_nmr', 'pairs.csv')
RUNS = '/Volumes/User Homes/Nu Project/BENDER-DYN/phospho_nmr'
MT_PAGE4 = '/Volumes/User Homes/Nu Project/BENDER-MT/page4_CALVADOS'
OUT = os.path.join(DYN, 'phospho_nmr')
SEED, NSHIFT, EXCL, MIN_RES, WN = 20261013, 10000, 3, 20, 0.14
NMOTIF = list(range(4, 13))
fdt, L, sp = pn.sp.fdt, pn.sp.L, pn.sp


# ---------------------------------------------------------------- measurement
def csp_table(pair, cs):
    sites = [int(s[1:]) for s in pair['sites'].split()]
    wt = pair['wt_seq']
    three = {'SEP': 'S', 'TPO': 'T'}
    u = cs[cs.bmrb_id == int(pair['bmrb_unmod'])].set_index('seq_id')
    p = cs[cs.bmrb_id == int(pair['bmrb_phos'])].set_index('seq_id')
    rows = []
    for s in sorted(set(u.index) & set(p.index)):
        au = three.get(u.loc[s, 'comp_id'], u.loc[s, 'aa_1letter'])
        ap = three.get(p.loc[s, 'comp_id'], p.loc[s, 'aa_1letter'])
        if au != ap or not (1 <= s <= len(wt)) or wt[s - 1] != au:
            continue
        d = min(abs(s - t) for t in sites)
        if d <= EXCL:
            continue
        rows.append({'pos': s, 'aa': au, 'd': d, 'dH': p.loc[s, 'cs_H'] - u.loc[s, 'cs_H'],
                     'dN': p.loc[s, 'cs_N'] - u.loc[s, 'cs_N'], 'dCA': p.loc[s, 'cs_CA'] - u.loc[s, 'cs_CA']})
    t = pd.DataFrame(rows)
    amide = t.dropna(subset=['dH', 'dN']).copy()
    amide['dH'] -= amide.dH.median()
    amide['dN'] -= amide.dN.median()
    amide['CSP'] = np.sqrt(amide.dH ** 2 + (WN * amide.dN) ** 2)
    ca = t.dropna(subset=['dCA']).copy()
    ca['absCA'] = (ca.dCA - ca.dCA.median()).abs()
    return amide[['pos', 'aa', 'd', 'CSP']], ca[['pos', 'd', 'absCA']]


# ---------------------------------------------------------------- statistics
def resid_rank(v, d):
    r = rankdata(v)
    X = np.column_stack([np.ones(len(r)), rankdata(np.log(d))])
    return r - X @ np.linalg.lstsq(X, r, rcond=None)[0]


def partial(pred, meas, d, rng=None):
    re_, rc_ = resid_rank(pred, d), resid_rank(meas, d)
    a, b = rankdata(re_), rankdata(rc_)                      # Spearman of the residuals = Pearson of their ranks
    a, b = a - a.mean(), b - b.mean()
    rho = float(a @ b / np.sqrt((a @ a) * (b @ b)))
    if rng is None:
        return rho, None
    ks = rng.integers(1, len(b), NSHIFT)
    null = np.array([a @ np.roll(b, k) for k in ks]) / np.sqrt((a @ a) * (b @ b))
    return rho, (1 + int((null >= rho).sum())) / (NSHIFT + 1)


def block_test(meas, d, pos, block):
    """Mean CSP residual of the residues in block against every circular contiguous block of the same size."""
    rc_ = resid_rank(meas, d)
    inb = np.isin(pos, block)
    k = int(inb.sum())
    obs = float(rc_[inb].mean())
    order = np.argsort(pos)
    rr = rc_[order]
    null = np.array([np.roll(rr, -s)[:k].mean() for s in range(len(rr))])
    return obs, k, float((null >= obs).mean())


# ---------------------------------------------------------------- prediction from runs
def load(dcd, top):
    return np.asarray(sp.bc.load_calvados(dcd, top, k=1000)[0], float)


def contact_stack(xyz):
    n = xyz.shape[1]
    far = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    C = np.zeros((len(xyz), n, n), bool)
    for f in range(len(xyz)):
        C[f] = (np.linalg.norm(xyz[f][:, None] - xyz[f][None], axis=2) < 0.8) & far
    return C


def phi_cond(xyz, seq, T, ionic, ph):
    q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params(), ph=ph)
    i, j, *_ = fdt.pair_setup(q, lam, sig)
    M = sp.site_matrix(i, j, q, len(seq))
    P = fdt.get_dielectric_params(T, ionic)
    k, eps, rc = P['kappa'], P['eps_yu'], fdt.RC_YU
    ec = np.exp(-k * rc)
    phi = np.zeros((len(xyz), len(seq)))
    step = max(1, int(sp.MAX_ELEMS // max(1, len(i))))
    for a in range(0, len(xyz), step):
        r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
        phi[a:a + step] = (M @ (eps * np.where(r < rc, np.exp(-k * r) / r - ec / rc, 0.0)).T).T
    return phi, (k, eps, rc, ec, P['beta'])


def run_one(job):
    name, dcd, top, key = job
    xyz = load(dcd, top)
    C = contact_stack(xyz)
    D, H = pn.networks(xyz, xyz.shape[1], key)
    return name, xyz, C.mean(0), D.mean(0), H.mean(0), C


def predictions(pair, results):
    prot, T, ph, ionic = pair['protein'], float(pair['T_K']), float(pair['pH']), float(pair['ionic_M'])
    W = [results[f'NMR-{prot}-WT-r{k}'] for k in (1, 2, 3)]
    P = [results[f'NMR-{prot}-P-r{k}'] for k in (1, 2, 3)]
    rr = lambda x, y: np.abs(x[2] - y[2]).sum(1)  # noqa: E731
    E = np.mean([rr(p, w) for p in P for w in W], 0) - np.mean([rr(W[a], W[b]) for a, b in ((0, 1), (0, 2), (1, 2))], 0)
    ddeg = np.abs(np.mean([p[3] for p in P], 0) - np.mean([w[3] for w in W], 0))
    docc = np.abs(np.mean([p[4] for p in P], 0) - np.mean([w[4] for w in W], 0))
    sites = [int(s[1:]) - 1 for s in pair['sites'].split()]
    rw = []
    for w in W:
        xyz, C = w[1], w[5]
        phi, (k, eps, rc, ec, beta) = phi_cond(xyz, pair['wt_seq'], T, ionic, ph)
        dH = pn.DQ * phi[:, sites].sum(1)
        for a in range(len(sites)):
            for b in range(a + 1, len(sites)):
                if abs(sites[a] - sites[b]) >= 2:                    # CALVADOS excludes bonded pairs only
                    r = np.linalg.norm(xyz[:, sites[a]] - xyz[:, sites[b]], axis=1)
                    dH += pn.DQ * pn.DQ * eps * np.where(r < rc, np.exp(-k * r) / r - ec / rc, 0.0)
        x = -beta * dH
        wt = np.exp(x - x.max())
        wt /= wt.sum()
        rw.append(np.abs(np.tensordot(wt, C.astype(float), axes=1) - w[2]).sum(1))
    return {'E': E, 'ddeg': ddeg, 'docc': docc, 'rw': np.mean(rw, 0)}


# ---------------------------------------------------------------- PAGE4 with the existing BENDER-MT runs
def page4_existing(cs, pairs, rng, lines):
    pair = next(p for p in pairs if p['protein'] == 'PAGE4')
    amide, ca = csp_table(pair, cs)
    lines += ['== PAGE4 with BENDER-MT existing runs (secondary; single pair r_i, no noise correction)',
              f'  analysed residues: {len(amide)} amide, {len(ca)} Calpha (T51 and residues 48-54 excluded)']
    out = {}
    for T in ('278K', '300K'):
        d = os.path.join(MT_PAGE4, T)
        pw = contact_stack(load(*[os.path.join(d, 'PAGE4-WT-I150-pH75', f) for f in ('PAGE4-WT-I150-pH75.dcd', 'top.pdb')])).mean(0)
        pp = contact_stack(load(*[os.path.join(d, 'PAGE4-HIPK1T51-I150-pH75', f)
                                  for f in ('PAGE4-HIPK1T51-I150-pH75.dcd', 'top.pdb')])).mean(0)
        r = np.abs(pp - pw).sum(1)
        out[T] = r
        tag = 'plan secondary' if T == '278K' else 'descriptive'
        idx = amide.pos.to_numpy() - 1
        rho, pshift = partial(r[idx], amide.CSP.to_numpy(), amide.d.to_numpy(), rng)
        plain = spearmanr(r[idx], amide.CSP).correlation
        idc = ca.pos.to_numpy() - 1
        rca, pca = partial(r[idc], ca.absCA.to_numpy(), ca.d.to_numpy(), rng)
        lines.append(f'  {T} ({tag}): partial Spearman (CSP | log d) {rho:+.3f}, circular-shift P {pshift:.4g}; plain '
                     f'Spearman {plain:+.3f}; Calpha partial {rca:+.3f} (P {pca:.4g})')
    dist = spearmanr(-np.log(amide.d), amide.CSP).correlation
    lines.append(f'  distance-only baseline: Spearman(-log d, CSP) {dist:+.3f}')
    obs, k, p = block_test(amide.CSP.to_numpy(), amide.d.to_numpy(), amide.pos.to_numpy(), NMOTIF)
    lines.append(f'  N-motif test (added before any shift difference): residues 4-12 analysed {k}/9; mean CSP residual '
                 f'{obs:+.2f} (ranks); circular block P {p:.4g}')
    nm = amide[amide.pos.isin(NMOTIF)]
    lines.append(f'  N-motif CSP median {nm.CSP.median():.4f} ppm vs other analysed residues {amide[~amide.pos.isin(NMOTIF)].CSP.median():.4f} ppm')
    lines.append('')
    per = amide.copy()
    per['r_278K'], per['r_300K'] = out['278K'][per.pos - 1], out['300K'][per.pos - 1]
    per.to_csv(os.path.join(OUT, 'per_residue_page4_existing.csv'), index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(10, 5.5), sharex=True)
    ax[0].plot(np.arange(1, 103), out['278K'], color='#1f77b4', lw=1.3, label='278 K')
    ax[0].plot(np.arange(1, 103), out['300K'], color='#d62728', lw=1.0, alpha=0.7, label='300 K')
    ax[0].set_ylabel('predicted rewiring  sum |dp|')
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].bar(amide.pos, amide.CSP, color='0.3', width=0.8)
    ax[1].set_ylabel('measured amide CSP (ppm)')
    ax[1].set_xlabel('PAGE4 residue  (green: N-motif 4-12; grey: T51 +-3, excluded)')
    for a in ax:
        a.axvspan(3.5, 12.5, color='#2ca02c', alpha=0.12, lw=0)
        a.axvspan(47.5, 54.5, color='0.85', lw=0)
        a.axvline(51, color='0.3', lw=0.8)
    ax[0].set_title('PAGE4 T51 phosphorylation: model (BENDER-MT runs) against NMR (BMRB 26544/26545)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'page4_existing.png'), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--page4-existing', action='store_true')
    a = ap.parse_args()
    for f in (PAIRS, os.path.join(OUT, 'page4_prediction_existing_runs.npz')):
        check_frozen(f)
    cs = pd.read_csv(SHIFTS, usecols=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])
    pairs = list(csv.DictReader(open(PAIRS)))
    rng = np.random.default_rng(SEED)
    lines = [f'phospho_nmr test (plan 1851cf7c); mode: {"PAGE4 existing runs only" if a.page4_existing else "full"}', '']
    page4_existing(cs, pairs, rng, lines)
    if a.page4_existing:
        txt = '\n'.join(lines) + '\n'
        open(os.path.join(OUT, 'summary_page4_existing.txt'), 'w').write(txt)
        print(txt)
        return

    jobs = []
    for p in pairs:
        for form in ('WT', 'P'):
            for k in (1, 2, 3):
                name = f'NMR-{p["protein"]}-{form}-r{k}'
                d = os.path.join(RUNS, f'{p["T_K"]}K', name)
                jobs.append((name, os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), name))
    missing = [j[0] for j in jobs if not os.path.isfile(j[1])]
    if missing:
        sys.exit(f'missing runs: {missing}')
    with Pool(8) as pool:
        results = {r[0]: r for r in pool.map(run_one, jobs)}
    rows, per = [], []
    for p in pairs:
        amide, ca = csp_table(p, cs)
        pred = predictions(p, results)
        if len(amide) < MIN_RES:
            rows.append({'protein': p['protein'], 'n_res': len(amide), 'analysed': False})
            continue
        idx, d, c = amide.pos.to_numpy() - 1, amide.d.to_numpy(), amide.CSP.to_numpy()
        row = {'protein': p['protein'], 'sites': p['sites'], 'n_res': len(amide), 'analysed': True}
        row['rho_p'], row['P_shift'] = partial(pred['E'][idx], c, d, rng)
        row['rho_plain'] = spearmanr(pred['E'][idx], c).correlation
        row['rho_dist'] = spearmanr(-np.log(d), c).correlation
        for k in ('ddeg', 'docc', 'rw'):
            row[f'rho_p_{k}'] = partial(pred[k][idx], c, d)[0]
        idc = ca.pos.to_numpy() - 1
        row['rho_p_CA'] = partial(pred['E'][idc], ca.absCA.to_numpy(), ca.d.to_numpy())[0] if len(ca) >= MIN_RES else np.nan
        rows.append(row)
        per.append(amide.assign(protein=p['protein'], E=pred['E'][idx], ddeg=pred['ddeg'][idx], docc=pred['docc'][idx],
                                rw=pred['rw'][idx]))
    R = pd.DataFrame(rows)
    pd.concat(per).to_csv(os.path.join(OUT, 'per_residue.csv'), index=False)
    R.to_csv(os.path.join(OUT, 'per_protein.csv'), index=False)
    A = R[R.analysed]
    w = wilcoxon(A.rho_p, alternative='greater')
    lines += ['== Primary (36 runs at the NMR conditions)',
              f'  proteins analysed: {len(A)} ({list(R[~R.analysed].protein)} below {MIN_RES} residues)']
    for _, r in A.sort_values('protein', key=lambda s: s.map({'PAGE4': 0, 'p53TAD': 1})).iterrows():
        lines.append(f'  {r.protein:9s} ({r.sites}; {r.n_res} res): rho_p {r.rho_p:+.3f} (shift P {r.P_shift:.3g}); plain '
                     f'{r.rho_plain:+.3f}; distance-only {r.rho_dist:+.3f}; |ddeg| {r.rho_p_ddeg:+.3f}, |docc| '
                     f'{r.rho_p_docc:+.3f}, reweighted {r.rho_p_rw:+.3f}, Calpha {r.rho_p_CA:+.3f}')
    lines.append(f'  PRIMARY: Wilcoxon one-sided rho_p > 0 over {len(A)} proteins, {int((A.rho_p > 0).sum())} positive, '
                 f'P {w.pvalue:.4g} -> {"SUCCESS" if w.pvalue < 0.05 else "NOT SUCCESS"}')
    B = A[A.protein != 'H1-0CTD']
    lines.append(f'  without H1.0: {int((B.rho_p > 0).sum())}/{len(B)} positive, P {wilcoxon(B.rho_p, alternative="greater").pvalue:.4g}')
    for k in ('rho_p_ddeg', 'rho_p_docc', 'rho_p_rw', 'rho_p_CA'):
        v = A[k].dropna()
        lines.append(f'  {k}: {int((v > 0).sum())}/{len(v)} positive, Wilcoxon P {wilcoxon(v, alternative="greater").pvalue:.4g}')
    txt = '\n'.join(lines) + '\n'
    open(os.path.join(OUT, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
