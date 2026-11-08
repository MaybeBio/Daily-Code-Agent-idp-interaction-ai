"""
phospho_nmr_amend2_test.py — Amendment 2 of phospho_nmr_plan.md (phospho_nmr_amend2_new_pairs.md): four new pairs
(alpha-synuclein pS129 and pY125, CD79b, Tir) in CALVADOS-2 or Mpipi-Recharged. Measurement, prediction and partial
correlation as in phospho_nmr_test.py (functions imported unchanged). Primary: mean rho_p over the analysed new pairs
against a joint circular-shift null (each pair's CSP residual ranks shifted independently), 10,000 draws,
default_rng(20261015), one-sided; CONFIRMED (CALVADOS-2) / REPLICATED (Mpipi-Recharged) if P < 0.05.

  ~/.venvs/finches/bin/python phospho_nmr_amend2_test.py --ff cal --qc      (run checks only)
  ~/.venvs/finches/bin/python phospho_nmr_amend2_test.py --ff cal | rech     ->  results_<ff>/{summary.txt, per_*.csv}
"""
import argparse
import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
sys.path.insert(0, EV)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import rankdata, spearmanr, wilcoxon  # noqa: E402

import phospho_nmr_test as T  # noqa: E402
import phospho_nmr_recharged_test as RT  # noqa: E402

PAIRS = os.path.join(HERE, 'pairs_nmr2.csv')
SHIFTS = os.path.join(HERE, 'shifts_nmr2.csv')
RUNS = {'cal': '/Volumes/User Homes/Nu Project/BENDER-DYN/phospho_nmr2',
        'rech': '/Volumes/User Homes/Nu Project/BENDER-DYN/phospho_nmr2_recharged'}
OLD = {'cal': os.path.join(EV, 'phospho_nmr', 'per_protein.csv'),
       'rech': os.path.join(EV, 'phospho_nmr_recharged', 'per_protein.csv')}
SITE_NAME = {'cal': {'S': 'SEP', 'T': 'TPO', 'Y': 'PTR'}, 'rech': {'S': 'GLU', 'T': 'GLU', 'Y': 'PTR'}}
THREE = {'S': 'SER', 'T': 'THR', 'Y': 'TYR'}
SEED, NDRAW = 20261015, 10000


def run_path(ff, p, form, k):
    name = f'NMR2-{p["protein"]}-{form}-r{k}'
    d = os.path.join(RUNS[ff], f'{p["T_K"]}K', name)
    return name, os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb')


def ca_resnames(top):
    """C-alpha residue names from ATOM and HETATM records (2026-10-15 erratum: OpenMM writes non-standard residues such
    as PTR as HETATM; the ATOM-only reader of phospho_nmr_recharged_test.py skipped them and failed QC falsely)."""
    return [line[17:20] for line in open(top) if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA']


def qc(ff, pairs):
    bad = []
    for p in pairs:
        wt, sites = p['wt_seq'], [int(s[1:]) for s in p['sites'].split()]
        for form in ('WT', 'P'):
            for k in (1, 2, 3):
                name, dcd, top = run_path(ff, p, form, k)
                if not (os.path.isfile(dcd) and os.path.isfile(top)):
                    bad.append(f'{name}: missing')
                    continue
                rn = ca_resnames(top)
                if len(rn) != len(wt):
                    bad.append(f'{name}: {len(rn)} residues, expected {len(wt)}')
                    continue
                for s in sites:
                    want = SITE_NAME[ff][wt[s - 1]] if form == 'P' else THREE[wt[s - 1]]
                    if rn[s - 1] != want:
                        bad.append(f'{name}: residue {s} is {rn[s - 1]}, expected {want}')
    return bad


def ranks(pred, meas, d):
    a, b = rankdata(T.resid_rank(pred, d)), rankdata(T.resid_rank(meas, d))
    return a - a.mean(), b - b.mean()


def joint_test(vectors, rng):
    """Mean rho_p over pairs against independent circular shifts of each pair's CSP residual ranks."""
    obs = np.mean([a @ b / np.sqrt((a @ a) * (b @ b)) for a, b in vectors])
    null = np.zeros(NDRAW)
    for a, b in vectors:
        ks = rng.integers(1, len(b), NDRAW)
        null += np.array([a @ np.roll(b, k) for k in ks]) / np.sqrt((a @ a) * (b @ b))
    null /= len(vectors)
    return float(obs), (1 + int((null >= obs).sum())) / (NDRAW + 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ff', choices=['cal', 'rech'], required=True)
    ap.add_argument('--qc', action='store_true')
    a = ap.parse_args()
    pairs = list(csv.DictReader(open(PAIRS)))
    bad = qc(a.ff, pairs)
    if bad:
        sys.exit('QC failed:\n  ' + '\n  '.join(bad))
    print(f'QC passed ({a.ff}): {6 * len(pairs)} runs present, phosphosites exactly at the listed residues')
    if a.qc:
        return
    out = os.path.join(HERE, f'results_{a.ff}')
    os.makedirs(out, exist_ok=True)
    cs = pd.read_csv(SHIFTS)
    jobs = [run_path(a.ff, p, f, k) for p in pairs for f in ('WT', 'P') for k in (1, 2, 3)]
    with Pool(8) as pool:
        results = {r[0]: r for r in pool.map(RT.run_one, jobs)}
    results = {k.replace('NMR2-', 'NMR-'): v for k, v in results.items()}       # RT.predictions expects NMR-<prot>-...
    rng = np.random.default_rng(SEED)
    rows, per, vec = [], [], {}
    for p in pairs:
        amide, _ = T.csp_table(p, cs)
        pred = RT.predictions(p['protein'], results)
        if len(amide) < T.MIN_RES:
            rows.append({'protein': p['protein'], 'n_res': len(amide), 'analysed': False})
            continue
        idx, d, c = amide.pos.to_numpy() - 1, amide.d.to_numpy(), amide.CSP.to_numpy()
        row = {'protein': p['protein'], 'sites': p['sites'], 'n_res': len(amide), 'analysed': True}
        row['rho_p'], row['P_shift'] = T.partial(pred['E'][idx], c, d, rng)
        row['rho_plain'] = spearmanr(pred['E'][idx], c).correlation
        row['rho_dist'] = spearmanr(-np.log(d), c).correlation
        for k in ('ddeg', 'docc'):
            row[f'rho_p_{k}'] = T.partial(pred[k][idx], c, d)[0]
        vec[p['protein']] = {k: ranks(pred[k][idx], c, d) for k in ('E', 'ddeg', 'docc')}
        rows.append(row)
        per.append(amide.assign(protein=p['protein'], E=pred['E'][idx], ddeg=pred['ddeg'][idx], docc=pred['docc'][idx]))
    R = pd.DataFrame(rows)
    pd.concat(per).to_csv(os.path.join(out, 'per_residue.csv'), index=False)
    R.to_csv(os.path.join(out, 'per_protein.csv'), index=False)
    A = R[R.analysed]
    name = {'cal': 'CALVADOS-2', 'rech': 'Mpipi-Recharged'}[a.ff]
    verdict = {'cal': ('CONFIRMED', 'NOT CONFIRMED'), 'rech': ('REPLICATED', 'NOT REPLICATED')}[a.ff]
    L = [f'phospho_nmr Amendment 2, {name}, {6 * len(pairs)} runs', '',
         f'pairs analysed: {len(A)} ({list(R[~R.analysed].protein)} below {T.MIN_RES} residues)']
    for _, r in A.iterrows():
        L.append(f'  {r.protein:9s} ({r.sites}; {r.n_res} res): rho_p {r.rho_p:+.3f} (shift P {r.P_shift:.3g}); plain '
                 f'{r.rho_plain:+.3f}; distance-only {r.rho_dist:+.3f}; |ddeg| {r.rho_p_ddeg:+.3f}, |docc| {r.rho_p_docc:+.3f}')
    obs, P = joint_test([vec[k]['E'] for k in A.protein], rng)
    L.append(f'PRIMARY: mean rho_p {obs:+.3f} over {len(A)} pairs, joint circular-shift P {P:.4g} -> '
             f'{verdict[0] if P < 0.05 else verdict[1]}')
    L.append(f'  positive pairs: {int((A.rho_p > 0).sum())}/{len(A)}')
    for lab, keep in (('without Tir', lambda s: s != 'Tir'), ('without aSynY125', lambda s: s != 'aSynY125'),
                      ('tyrosine pairs (aSynY125, CD79b, Tir)', lambda s: s in ('aSynY125', 'CD79b', 'Tir'))):
        ks = [k for k in A.protein if keep(k)]
        if ks:
            o, pp = joint_test([vec[k]['E'] for k in ks], rng)
            L.append(f'  {lab}: mean rho_p {o:+.3f} over {len(ks)}, P {pp:.4g}')
    for k in ('ddeg', 'docc'):
        o, pp = joint_test([vec[x][k] for x in A.protein], rng)
        L.append(f'  |d{k[1:]}| in place of E: mean rho_p {o:+.3f}, P {pp:.4g}')
    old = pd.read_csv(OLD[a.ff])
    allr = np.r_[old[old.analysed].rho_p.to_numpy(float), A.rho_p.to_numpy(float)]
    L.append(f'  descriptive, all {len(allr)} pairs (6 already known): {int((allr > 0).sum())} positive, one-sided '
             f'Wilcoxon P {wilcoxon(allr, alternative="greater").pvalue:.4g}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(out, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
