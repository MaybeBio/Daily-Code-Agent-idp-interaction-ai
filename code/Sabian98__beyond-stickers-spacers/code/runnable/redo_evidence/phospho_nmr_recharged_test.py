"""
phospho_nmr_recharged_test.py — Amendment 1 to phospho_nmr_plan.md (phospho_nmr_amend1_recharged.md): the NMR case
series with Mpipi-Recharged runs of the same constructs. Reuses the measurement, prediction and statistics functions of
phospho_nmr_test.py unchanged; the charge-model reweighting and the PAGE4 BENDER-MT check are not repeated.

  ~/.venvs/finches/bin/python phospho_nmr_recharged_test.py --qc      (run checks only; no shift difference computed)
  ~/.venvs/finches/bin/python phospho_nmr_recharged_test.py           ->  phospho_nmr_recharged/{summary.txt, per_*.csv}
"""
import argparse
import csv
import os
import sys
from multiprocessing import Pool

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr, wilcoxon  # noqa: E402

import phospho_nmr_test as T  # noqa: E402
from phospho_direct_test import check_frozen  # noqa: E402

RUNS = (__import__('os').environ['BSS_DATA'] + '/BENDER-DYN/phospho_nmr_recharged')
OUT = os.path.join(DYN, 'phospho_nmr_recharged')
CAL_PER = os.path.join(DYN, 'phospho_nmr', 'per_residue.csv')
SEED = 20261014
THREE = {'S': 'SER', 'T': 'THR'}


def run_path(p, form, k):
    name = f'NMR-{p["protein"]}-{form}-r{k}'
    d = os.path.join(RUNS, f'{p["T_K"]}K', name)
    return name, os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb')


def ca_resnames(top):
    return [line[17:20] for line in open(top) if line.startswith('ATOM') and line[12:16].strip() == 'CA']


def qc(pairs):
    """Every run present; GLU at exactly the sites in phosphorylated runs, SER/THR there in unmodified runs."""
    bad = []
    for p in pairs:
        wt, sites = p['wt_seq'], [int(s[1:]) for s in p['sites'].split()]
        for form in ('WT', 'P'):
            for k in (1, 2, 3):
                name, dcd, top = run_path(p, form, k)
                if not (os.path.isfile(dcd) and os.path.isfile(top)):
                    bad.append(f'{name}: missing')
                    continue
                rn = ca_resnames(top)
                if len(rn) != len(wt):
                    bad.append(f'{name}: {len(rn)} residues, expected {len(wt)}')
                    continue
                for s in sites:
                    want = 'GLU' if form == 'P' else THREE[wt[s - 1]]
                    if rn[s - 1] != want:
                        bad.append(f'{name}: residue {s} is {rn[s - 1]}, expected {want}')
    return bad


def run_one(job):
    name, dcd, top = job
    xyz = T.load(dcd, top)
    C = T.contact_stack(xyz)
    D, H = T.pn.networks(xyz, xyz.shape[1], name)
    return name, None, C.mean(0), D.mean(0), H.mean(0)


def predictions(prot, results):
    W = [results[f'NMR-{prot}-WT-r{k}'] for k in (1, 2, 3)]
    P = [results[f'NMR-{prot}-P-r{k}'] for k in (1, 2, 3)]
    rr = lambda x, y: np.abs(x[2] - y[2]).sum(1)  # noqa: E731
    E = np.mean([rr(p, w) for p in P for w in W], 0) - np.mean([rr(W[a], W[b]) for a, b in ((0, 1), (0, 2), (1, 2))], 0)
    ddeg = np.abs(np.mean([p[3] for p in P], 0) - np.mean([w[3] for w in W], 0))
    docc = np.abs(np.mean([p[4] for p in P], 0) - np.mean([w[4] for w in W], 0))
    return {'E': E, 'ddeg': ddeg, 'docc': docc}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--qc', action='store_true')
    a = ap.parse_args()
    check_frozen(T.PAIRS)
    check_frozen(CAL_PER)
    pairs = list(csv.DictReader(open(T.PAIRS)))
    bad = qc(pairs)
    if bad:
        sys.exit('QC failed:\n  ' + '\n  '.join(bad))
    print('QC passed: 36 runs present, phosphate beads exactly at the listed sites')
    if a.qc:
        return
    os.makedirs(OUT, exist_ok=True)
    cs = pd.read_csv(T.SHIFTS, usecols=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])
    jobs = [run_path(p, f, k) for p in pairs for f in ('WT', 'P') for k in (1, 2, 3)]
    with Pool(8) as pool:
        results = {r[0]: r for r in pool.map(run_one, jobs)}
    cal = pd.read_csv(CAL_PER)
    rng = np.random.default_rng(SEED)
    rows, per = [], []
    for p in pairs:
        amide, _ = T.csp_table(p, cs)
        pred = predictions(p['protein'], results)
        if len(amide) < T.MIN_RES:
            rows.append({'protein': p['protein'], 'n_res': len(amide), 'analysed': False})
            continue
        idx, d, c = amide.pos.to_numpy() - 1, amide.d.to_numpy(), amide.CSP.to_numpy()
        row = {'protein': p['protein'], 'sites': p['sites'], 'n_res': len(amide), 'analysed': True}
        row['rho_p'], row['P_shift'] = T.partial(pred['E'][idx], c, d, rng)
        row['rho_plain'] = spearmanr(pred['E'][idx], c).correlation
        for k in ('ddeg', 'docc'):
            row[f'rho_p_{k}'] = T.partial(pred[k][idx], c, d)[0]
        ce = cal[cal.protein == p['protein']].set_index('pos').E.reindex(amide.pos).to_numpy()
        row['rho_E_cal'] = spearmanr(pred['E'][idx], ce).correlation
        rows.append(row)
        per.append(amide.assign(protein=p['protein'], E=pred['E'][idx], ddeg=pred['ddeg'][idx], docc=pred['docc'][idx],
                                E_cal=ce))
    R = pd.DataFrame(rows)
    pd.concat(per).to_csv(os.path.join(OUT, 'per_residue.csv'), index=False)
    R.to_csv(os.path.join(OUT, 'per_protein.csv'), index=False)
    A = R[R.analysed]
    w = wilcoxon(A.rho_p, alternative='greater')
    npos = int((A.rho_p > 0).sum())
    reading = ('REPLICATED' if w.pvalue < 0.05 else
               'DIRECTION REPLICATED, NOT SIGNIFICANT' if npos >= 4 else 'NOT REPLICATED')
    L = ['phospho_nmr Amendment 1 (Mpipi-Recharged), 36 runs at the NMR conditions', '',
         f'proteins analysed: {len(A)} ({list(R[~R.analysed].protein)} below {T.MIN_RES} residues)']
    for _, r in A.iterrows():
        L.append(f'  {r.protein:9s} ({r.sites}; {r.n_res} res): rho_p {r.rho_p:+.3f} (shift P {r.P_shift:.3g}); plain '
                 f'{r.rho_plain:+.3f}; |ddeg| {r.rho_p_ddeg:+.3f}, |docc| {r.rho_p_docc:+.3f}; E vs CALVADOS-2 E '
                 f'{r.rho_E_cal:+.3f}')
    L.append(f'PRIMARY: Wilcoxon one-sided rho_p > 0 over {len(A)} proteins, {npos} positive, P {w.pvalue:.4g} -> '
             f'{"SUCCESS" if w.pvalue < 0.05 else "NOT SUCCESS"}')
    L.append(f'READING: {reading}')
    for k in ('rho_p_ddeg', 'rho_p_docc'):
        v = A[k].dropna()
        L.append(f'  {k}: {int((v > 0).sum())}/{len(v)} positive, Wilcoxon P {wilcoxon(v, alternative="greater").pvalue:.4g}')
    L.append(f'  agreement with CALVADOS-2 E (per-protein Spearman): median {A.rho_E_cal.median():+.3f}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
