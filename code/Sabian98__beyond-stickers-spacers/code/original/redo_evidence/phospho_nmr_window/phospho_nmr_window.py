"""
phospho_nmr_window.py — POST HOC sensitivity of the NMR case series to the exclusion window around each phosphosite.
Requested 2026-10-15, after all results of the original plan and Amendments 1 and 2 were known; descriptive only, no
decision attached. The plans exclude residues within 3 positions of a site; this reruns the identical analysis with
windows of 2 and 5 (3 repeated as a check that the pipeline reproduces the hashed values).

Everything except the window is unchanged: CSP table, E_i, partial Spearman given log distance and the circular-shift
statistics are the hashed functions of phospho_nmr_test.py (T.csp_table reads the module-level T.EXCL, which is set per
window). E_i and |d deg_i| are computed from the contact maps only (hub occupancy is not recomputed). Same runs, same
shift tables, same pairs (original 6 + Amendment 2's 4) in CALVADOS-2 and Mpipi-Recharged.
Tests per window and force field: one-sided Wilcoxon over the 6 original pairs, the joint circular-shift test of
Amendment 2 over the 4 new pairs, and the Wilcoxon over all 10. Seeds as in the plans (original 20261013 for
CALVADOS-2 and 20261014 for Mpipi-Recharged per-protein shifts; 20261015 for the joint test).

  ~/.venvs/finches/bin/python phospho_nmr_window.py   ->   summary.txt, per_protein.csv
"""
import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
sys.path[:0] = [EV, os.path.join(EV, 'phospho_nmr_amend2')]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

import phospho_nmr_test as T  # noqa: E402
import phospho_nmr_amend2_test as A2  # noqa: E402

BASE = '/Volumes/User Homes/Nu Project/BENDER-DYN'
RUNS = {('orig', 'cal'): (os.path.join(BASE, 'phospho_nmr'), 'NMR'),
        ('orig', 'rech'): (os.path.join(BASE, 'phospho_nmr_recharged'), 'NMR'),
        ('new', 'cal'): (os.path.join(BASE, 'phospho_nmr2'), 'NMR2'),
        ('new', 'rech'): (os.path.join(BASE, 'phospho_nmr2_recharged'), 'NMR2')}
SEED_ORIG = {'cal': 20261013, 'rech': 20261014}
WINDOWS = (2, 3, 5)


def maps(job):
    name, dcd, top = job
    C = T.contact_stack(T.load(dcd, top))
    return name, C.mean(0), C.sum(2).mean(0)


def predict(res, prefix, prot):
    W = [res[f'{prefix}-{prot}-WT-r{k}'] for k in (1, 2, 3)]
    P = [res[f'{prefix}-{prot}-P-r{k}'] for k in (1, 2, 3)]
    rr = lambda x, y: np.abs(x[1] - y[1]).sum(1)  # noqa: E731
    E = np.mean([rr(p, w) for p in P for w in W], 0) - np.mean([rr(W[a], W[b]) for a, b in ((0, 1), (0, 2), (1, 2))], 0)
    ddeg = np.abs(np.mean([p[2] for p in P], 0) - np.mean([w[2] for w in W], 0))
    return E, ddeg


def main():
    sets = {'orig': (list(csv.DictReader(open(T.PAIRS))),
                     pd.read_csv(T.SHIFTS, usecols=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])),
            'new': (list(csv.DictReader(open(A2.PAIRS))), pd.read_csv(A2.SHIFTS))}
    jobs = []
    for (s, ff), (root, prefix) in RUNS.items():
        for p in sets[s][0]:
            for form in ('WT', 'P'):
                for k in (1, 2, 3):
                    name = f'{prefix}-{p["protein"]}-{form}-r{k}'
                    d = os.path.join(root, f'{p["T_K"]}K', name)
                    jobs.append((f'{ff}|{name}', os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb')))
    with Pool(8) as pool:
        out = pool.map(maps, jobs)
    res = {ff: {} for ff in ('cal', 'rech')}
    for key, C, deg in out:
        ff, name = key.split('|')
        res[ff][name] = (name, C, deg)
    rows, L = [], ['phospho_nmr exclusion-window sensitivity (POST HOC, descriptive; plans use a window of 3)', '']
    for ff in ('cal', 'rech'):
        for w in WINDOWS:
            T.EXCL = w
            rng_o, rng_j = np.random.default_rng(SEED_ORIG[ff]), np.random.default_rng(20261015)
            vec = {}
            for s in ('orig', 'new'):
                pairs, cs = sets[s]
                prefix = RUNS[(s, ff)][1]
                for p in pairs:
                    amide, _ = T.csp_table(p, cs)
                    row = {'ff': ff, 'window': w, 'set': s, 'protein': p['protein'], 'n_res': len(amide)}
                    if len(amide) >= T.MIN_RES:
                        E, ddeg = predict(res[ff], prefix, p['protein'])
                        idx, d, c = amide.pos.to_numpy() - 1, amide.d.to_numpy(), amide.CSP.to_numpy()
                        row['rho_p'], row['P_shift'] = T.partial(E[idx], c, d, rng_o if s == 'orig' else rng_j)
                        row['rho_p_ddeg'] = T.partial(ddeg[idx], c, d)[0]
                        if s == 'new':
                            vec[p['protein']] = A2.ranks(E[idx], c, d)
                    rows.append(row)
            R = pd.DataFrame([r for r in rows if r['ff'] == ff and r['window'] == w and 'rho_p' in r])
            o, n_ = R[R.set == 'orig'], R[R.set == 'new']
            jo, jp = A2.joint_test([vec[k] for k in n_.protein], rng_j)
            L.append(f'{"CALVADOS-2" if ff == "cal" else "Mpipi-Recharged"}, window {w}'
                     + ('  (as pre-registered)' if w == 3 else ''))
            L.append('  ' + '  '.join(f'{r.protein} {r.rho_p:+.2f} (n {r.n_res})' for _, r in R.iterrows()))
            L.append(f'  original 6: {int((o.rho_p > 0).sum())}/{len(o)} positive, Wilcoxon P '
                     f'{wilcoxon(o.rho_p, alternative="greater").pvalue:.3g}; '
                     f'new 4: mean {jo:+.3f}, joint P {jp:.3g}, {int((n_.rho_p > 0).sum())}/{len(n_)} positive; '
                     f'all {len(R)}: {int((R.rho_p > 0).sum())} positive, Wilcoxon P '
                     f'{wilcoxon(R.rho_p, alternative="greater").pvalue:.3g}, median {R.rho_p.median():+.3f}')
            L.append(f'  |d deg| in place of E, all {len(R)}: {int((R.rho_p_ddeg > 0).sum())} positive, Wilcoxon P '
                     f'{wilcoxon(R.rho_p_ddeg, alternative="greater").pvalue:.3g}')
            L.append('')
    pd.DataFrame(rows).to_csv(os.path.join(HERE, 'per_protein.csv'), index=False)
    txt = '\n'.join(L)
    open(os.path.join(HERE, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
