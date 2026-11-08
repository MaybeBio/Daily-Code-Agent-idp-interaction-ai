"""
phospho_nmr_ceiling.py — POST HOC estimate of the attainable partial correlation in the NMR case series (requested
2026-10-21, after all results were known; descriptive, no decision attached).

The observed partial Spearman rho_p between predicted rewiring E_i and the amide CSP_i (given log distance to the nearest
site) is attenuated by noise on both sides, so its ceiling is below 1. Classical attenuation: the attainable correlation
is about sqrt(R_model * R_meas), with R the reliability of each profile (here of the distance-adjusted rank profile).

Measurement reliability R_meas: Gaussian noise is added independently to every assigned H and N shift of both BMRB entries
(sigma_H, sigma_N), the CSP table is recomputed with the hashed phospho_nmr_test.csp_table (same matching, offsets and
+-3 exclusion), and R_meas is the mean partial Spearman between two independently noised CSP profiles (200 pairs of draws).
Treating the observed shifts as the truth overstates the signal, so R_meas and the ceiling are generous (upper bounds).
Two noise settings: (0.01, 0.05) ppm and (0.02, 0.10) ppm for (H, N).

Model reliability R_model: r_i(P_a, W_a) for the three disjoint replicate pairs (a = 1, 2, 3); R_1 = mean partial Spearman
between two such single-pair profiles over the analysed residues; the 3-replicate prediction is taken as an average of
three independent pairs, R_model = 3 R_1 / (1 + 2 R_1) (Spearman-Brown). Negative R_1 is reported and the ceiling left blank.

Same pairs, shift tables and runs as phospho_nmr_window.py (original 6 + Amendment 2's 4, CALVADOS-2 and Mpipi-Recharged).
numpy default_rng(20261021).
  ~/.venvs/finches/bin/python phospho_nmr_ceiling.py   ->   summary.txt, per_pair.csv
"""
import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
sys.path[:0] = [EV, os.path.join(EV, 'phospho_nmr_amend2'), os.path.join(EV, 'phospho_nmr_window')]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import phospho_nmr_test as T  # noqa: E402
import phospho_nmr_amend2_test as A2  # noqa: E402
import phospho_nmr_window as W  # noqa: E402

SEED, NDRAW = 20261021, 200
NOISE = ((0.01, 0.05), (0.02, 0.10))


def noisy(cs, ids, sh, sn, rng):
    c = cs.copy()
    m = c.bmrb_id.isin(ids).to_numpy()
    c.loc[m, 'cs_H'] = c.loc[m, 'cs_H'] + rng.normal(0, sh, m.sum())
    c.loc[m, 'cs_N'] = c.loc[m, 'cs_N'] + rng.normal(0, sn, m.sum())
    return c


def meas_reliability(p, cs, sh, sn, rng):
    ids = [int(p['bmrb_unmod']), int(p['bmrb_phos'])]
    rs = []
    for _ in range(NDRAW):
        a, _ = T.csp_table(p, noisy(cs, ids, sh, sn, rng))
        b, _ = T.csp_table(p, noisy(cs, ids, sh, sn, rng))
        assert (a.pos.to_numpy() == b.pos.to_numpy()).all()
        rs.append(T.partial(a.CSP.to_numpy(), b.CSP.to_numpy(), a.d.to_numpy())[0])
    return float(np.mean(rs))


def main():
    T.EXCL = 3
    sets = {'orig': (list(csv.DictReader(open(T.PAIRS))),
                     pd.read_csv(T.SHIFTS, usecols=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])),
            'new': (list(csv.DictReader(open(A2.PAIRS))), pd.read_csv(A2.SHIFTS))}
    jobs = []
    for (s, ff), (root, prefix) in W.RUNS.items():
        for p in sets[s][0]:
            for form in ('WT', 'P'):
                for k in (1, 2, 3):
                    name = f'{prefix}-{p["protein"]}-{form}-r{k}'
                    d = os.path.join(root, f'{p["T_K"]}K', name)
                    jobs.append((f'{ff}|{name}', os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb')))
    with Pool(8) as pool:
        out = pool.map(W.maps, jobs)
    res = {'cal': {}, 'rech': {}}
    for key, C, deg in out:
        ff, name = key.split('|')
        res[ff][name] = C
    obs = {}
    for ff, files in (('cal', ('phospho_nmr/per_protein.csv', 'phospho_nmr_amend2/results_cal/per_protein.csv')),
                      ('rech', ('phospho_nmr_recharged/per_protein.csv', 'phospho_nmr_amend2/results_rech/per_protein.csv'))):
        for f in files:
            for r in csv.DictReader(open(os.path.join(EV, f))):
                if r['analysed'] == 'True':
                    obs[(ff, r['protein'])] = float(r['rho_p'])
    rng = np.random.default_rng(SEED)
    rows = []
    for s in ('orig', 'new'):
        pairs, cs = sets[s]
        for p in pairs:
            amide, _ = T.csp_table(p, cs)
            idx, dd = amide.pos.to_numpy() - 1, amide.d.to_numpy()
            rm = {f'R_meas_{sh}_{sn}': meas_reliability(p, cs, sh, sn, rng) for sh, sn in NOISE}
            for ff in ('cal', 'rech'):
                prefix = W.RUNS[(s, ff)][1]
                single = [np.abs(res[ff][f'{prefix}-{p["protein"]}-P-r{k}'] - res[ff][f'{prefix}-{p["protein"]}-WT-r{k}']).sum(1)[idx]
                          for k in (1, 2, 3)]
                r1 = float(np.mean([T.partial(single[a], single[b], dd)[0] for a, b in ((0, 1), (0, 2), (1, 2))]))
                rmod = 3 * r1 / (1 + 2 * r1) if r1 > 0 else np.nan
                row = {'protein': p['protein'], 'set': s, 'ff': ff, 'n_res': len(amide), 'R1_single_pair': r1,
                       'R_model': rmod, 'rho_p_observed': obs[(ff, p['protein'])], **rm}
                for sh, sn in NOISE:
                    k_ = f'{sh}_{sn}'
                    row[f'ceiling_{k_}'] = float(np.sqrt(rmod * rm[f'R_meas_{k_}'])) if rmod > 0 and rm[f'R_meas_{k_}'] > 0 else np.nan
                rows.append(row)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(HERE, 'per_pair.csv'), index=False)
    L = ['phospho_nmr attainable-correlation estimate (POST HOC, descriptive)', '',
         'per pair: R_meas (H 0.01/N 0.05 ppm | H 0.02/N 0.10 ppm), R_model, ceiling (both noise settings), observed rho_p']
    for ff in ('cal', 'rech'):
        L.append(f'{"CALVADOS-2" if ff == "cal" else "Mpipi-Recharged"}')
        g = R[R.ff == ff]
        for _, r in g.iterrows():
            L.append(f'  {r.protein:9s} n {r.n_res:3d}  R_meas {r["R_meas_0.01_0.05"]:+.2f} | {r["R_meas_0.02_0.1"]:+.2f}  '
                     f'R1 {r.R1_single_pair:+.2f} R_model {r.R_model:+.2f}  ceiling {r["ceiling_0.01_0.05"]:.2f} | '
                     f'{r["ceiling_0.02_0.1"]:.2f}  observed {r.rho_p_observed:+.2f}')
        for k_ in ('0.01_0.05', '0.02_0.1'):
            c = g[f'ceiling_{k_}']
            L.append(f'  noise {k_}: median ceiling {c.median():.2f} (range {c.min():.2f}-{c.max():.2f}, {c.notna().sum()} pairs); '
                     f'median observed {g.rho_p_observed.median():+.2f}; median observed/ceiling '
                     f'{(g.rho_p_observed / c).median():.2f}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(HERE, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
