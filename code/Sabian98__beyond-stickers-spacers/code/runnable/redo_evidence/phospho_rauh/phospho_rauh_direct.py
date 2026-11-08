"""
phospho_rauh_direct.py — scores the 12 direct CALVADOS-2 runs with the published phosphoresidue model (Rauh et al.,
Biophys. J. 2026; inputs from make_inputs_rauh.py) for the four NMR pairs that the reweighting check (phospho_rauh.py)
left with low ESS or a changed result: Sic1, H1-0 CTD, NHE1 tail, Tir. Exploratory robustness check, written before the
runs were read.

For each pair, E_i is the noise-corrected rewiring of phospho_nmr_test (mean of the 9 r_i(P_a, W_b) minus the mean of
the 3 r_i(W_a, W_b), contact frequencies over 1,000 second-half frames) with P = the three RAUH runs and W = the three
existing unmodified runs, and rho_p is the partial Spearman of E_i with the amide CSP given log distance (csp_table and
partial of phospho_nmr_test, window 3), with the per-pair circular-shift P of the original analysis (10,000 shifts,
numpy default_rng(20261028)). Reported next to rho_p of the simulated phosphate (original P runs, same code) and of the
reweighted estimate (phospho_rauh/per_pair.csv); also the mean Rg of each form, the Spearman correlation between the
E_i profiles of the two phosphate models, and the 10-pair summary with the six other pairs taken from the reweighted
check (descriptive, labelled as such).
  BENDER_REANALYSIS_DIR=$BSS_REPO/bender_reanalysis ~/.venvs/finches/bin/python phospho_rauh_direct.py
      ->  summary_direct.txt, per_pair_direct.csv
"""
import csv
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
os.environ.setdefault('BENDER_REANALYSIS_DIR', (__import__('os').environ['BSS_REPO'] + '/bender_reanalysis'))
sys.path[:0] = [(__import__('os').environ['BSS_DATA'] + '/BENDER-MT/code'), EV, os.path.join(EV, 'phospho_nmr_amend2'),
                os.path.join(EV, 'phospho_nmr_window')]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr, wilcoxon  # noqa: E402

import phospho_nmr_test as T  # noqa: E402
import phospho_nmr_amend2_test as A2  # noqa: E402
import phospho_nmr_window as W  # noqa: E402

RAUH = (__import__('os').environ['BSS_DATA'] + '/BENDER-DYN/phospho_rauh')
PAIRS = {'Sic1': 'orig', 'H1-0CTD': 'orig', 'NHE1tail': 'orig', 'Tir': 'new'}
PHOS3 = {'SEP': 'S', 'TPO': 'T', 'PTR': 'Y'}


def run_stats(job):
    name, dcd, top = job
    xyz = T.load(dcd, top)
    seq3 = [l[17:20] for l in open(top) if l.startswith('ATOM') and l[12:16].strip() == 'CA']
    c = xyz - xyz.mean(1, keepdims=True)
    rg = float(np.sqrt((c ** 2).sum(2).mean(1)).mean())
    return name, T.contact_stack(xyz).mean(0), rg, seq3, len(xyz)


def main():
    T.EXCL = 3
    sets = {'orig': (list(csv.DictReader(open(T.PAIRS))),
                     pd.read_csv(T.SHIFTS, usecols=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])),
            'new': (list(csv.DictReader(open(A2.PAIRS))), pd.read_csv(A2.SHIFTS))}
    jobs, meta = [], {}
    for prot, s in PAIRS.items():
        root, prefix = W.RUNS[(s, 'cal')]
        p = next(r for r in sets[s][0] if r['protein'] == prot)
        meta[prot] = (p, s, prefix)
        for k in (1, 2, 3):
            for form, base, pre in (('P', RAUH, 'RAUH'), ('P', root, prefix), ('WT', root, prefix)):
                name = f'{pre}-{prot}-{form}-r{k}'
                d = os.path.join(base, f'{p["T_K"]}K', name)
                jobs.append((name, os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb')))
    with Pool(8) as pool:
        out = {r[0]: r[1:] for r in pool.map(run_stats, jobs)}
    rw = {r['protein']: r for r in csv.DictReader(open(os.path.join(HERE, 'per_pair.csv')))}
    rng = np.random.default_rng(20261028)
    rows, L = [], ['NMR pairs with the published phosphoresidue model, direct CALVADOS-2 runs (exploratory)', '']
    for prot, (p, s, prefix) in meta.items():
        sites = [int(x[1:]) - 1 for x in p['sites'].split()]
        amide, _ = T.csp_table(p, sets[s][1])
        idx, dd, csp = amide.pos.to_numpy() - 1, amide.d.to_numpy(), amide.CSP.to_numpy()
        Wr = [out[f'{prefix}-{prot}-WT-r{k}'] for k in (1, 2, 3)]
        noise = np.mean([np.abs(Wr[a][0] - Wr[b][0]).sum(1) for a, b in ((0, 1), (0, 2), (1, 2))], 0)
        res = {}
        for lab, pre in (('rauh', 'RAUH'), ('sim', prefix)):
            Pr = [out[f'{pre}-{prot}-P-r{k}'] for k in (1, 2, 3)]
            for r_ in Pr:   # phosphoresidues exactly at the listed sites, full frame count
                ph = [i for i, x in enumerate(r_[2]) if x in PHOS3]
                assert ph == sites, (pre, prot, ph, sites)
                assert r_[3] == 1000, (pre, prot, r_[3])
            E = np.mean([np.abs(a[0] - b[0]).sum(1) for a in Pr for b in Wr], 0) - noise
            rho, pv = T.partial(E[idx], csp, dd, rng)
            res[lab] = (E, rho, pv, np.mean([r_[1] for r_ in Pr]))
        rg_w = np.mean([w[1] for w in Wr])
        row = {'protein': prot, 'n_res': len(amide), 'rho_p_sim': res['sim'][1], 'P_sim': res['sim'][2],
               'rho_p_rauh_direct': res['rauh'][1], 'P_rauh_direct': res['rauh'][2],
               'rho_p_rauh_reweighted': float(rw[prot]['rho_p_full']), 'ess_min_reweighted': float(rw[prot]['ess_min_full']),
               'E_profile_spearman': spearmanr(res['sim'][0][idx], res['rauh'][0][idx]).correlation,
               'Rg_W_nm': rg_w, 'dRg_pct_sim': (res['sim'][3] / rg_w - 1) * 100,
               'dRg_pct_rauh': (res['rauh'][3] / rg_w - 1) * 100}
        rows.append(row)
        L.append(f'  {prot:9s} n {len(amide):3d}  simulated phosphate {row["rho_p_sim"]:+.3f} (P {row["P_sim"]:.3f})  '
                 f'published model direct {row["rho_p_rauh_direct"]:+.3f} (P {row["P_rauh_direct"]:.3f})  '
                 f'reweighted {row["rho_p_rauh_reweighted"]:+.3f} (ESS {row["ess_min_reweighted"]:.2f})  '
                 f'E profiles agree {row["E_profile_spearman"]:+.2f}  dRg {row["dRg_pct_sim"]:+.1f}% / {row["dRg_pct_rauh"]:+.1f}%')
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(HERE, 'per_pair_direct.csv'), index=False)
    allp = {k: float(v['rho_p_full']) for k, v in rw.items()}
    sim10 = {k: float(v['rho_p_simulated']) for k, v in rw.items()}
    for _, r in R.iterrows():
        allp[r.protein] = r.rho_p_rauh_direct
    v = pd.Series(allp)
    L += ['', f'  four pairs: published model direct mean {R.rho_p_rauh_direct.mean():+.3f}, simulated phosphate '
              f'{R.rho_p_sim.mean():+.3f}, reweighted estimate {R.rho_p_rauh_reweighted.mean():+.3f}',
          f'  all ten pairs with the published model (these four direct, the other six reweighted, min ESS 0.37 to 0.88): '
          f'{int((v > 0).sum())}/10 positive, mean {v.mean():+.3f}, median {v.median():+.3f}, one-sided Wilcoxon P '
          f'{wilcoxon(v, alternative="greater").pvalue:.3g}',
          f'  simulated phosphate, same ten: {int((pd.Series(sim10) > 0).sum())}/10 positive, mean {pd.Series(sim10).mean():+.3f}']
    txt = '\n'.join(L) + '\n'
    open(os.path.join(HERE, 'summary_direct.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
