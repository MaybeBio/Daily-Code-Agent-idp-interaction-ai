"""
heat_crosses_crossmodel.py — the compaction rule in CALVADOS-2 and Mpipi-GG on the same 211 chains (METHODS_solvation §13),
from heat_crosses.py (--ids heat_crosses_shared/ids.txt) and heat_crosses_gg.py outputs. Replaces the inline calculation
that first wrote analysis_crossmodel.csv (2026-10-10), so the numbers can be regenerated.
Per model: usable sites (effective sample > 0.1), median |change of the total| for phosphorylation and doubled salt,
share beyond 2 SE, Spearman of the change with the size change and (salt) with |NCPR|; then site-by-site and chain-by-chain
agreement between the models (both usable).
    python heat_crosses_crossmodel.py [--dir .../heat_crosses_shared]
Output: <dir>/analysis_crossmodel.csv, crossmodel_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

TR = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis/Analysis_FDT/total_response')
BIO = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/BENDER_BIO.csv')


def main(d):
    ncpr = pd.read_csv(BIO, usecols=['UniProt_ID', 'NCPR'], low_memory=False).set_index('UniProt_ID').NCPR.abs()
    F = {'CALVADOS-2': ('chains.csv', 'sites.csv'), 'Mpipi-GG': ('chains_gg.csv', 'sites_gg.csv')}
    C, S, rows = {}, {}, []
    ids = set(open(os.path.join(d, 'ids.txt')).read().split())
    for m, (fc, fs) in F.items():
        c = pd.read_csv(os.path.join(d, fc)).drop_duplicates('UniProt_ID').set_index('UniProt_ID')
        s = pd.read_csv(os.path.join(d, fs)).drop_duplicates(['UniProt_ID', 'pos'])
        c, s = c[c.index.isin(ids)], s[s.UniProt_ID.isin(ids)]
        C[m], S[m] = c, s
        su, cu = s[s.ess > 0.1], c[c.salt_ess > 0.1]
        rows.append({'model': m, 'chains': len(c), 'sites': len(s), 'site_ess_ok': (s.ess > 0.1).mean(),
                     'phospho |d| median': su.phospho_dtotal.abs().median(),
                     'phospho >2SE': (su.phospho_dtotal.abs() > 2 * su.phospho_dtotal_se).mean(),
                     'rho(phospho d, dRg)': stats.spearmanr(su.phospho_dtotal, su.phospho_dRg_pct)[0],
                     'salt |d| median': cu.salt_dtotal.abs().median(),
                     'salt >2SE': (cu.salt_dtotal.abs() > 2 * cu.salt_dtotal_se).mean(),
                     'rho(salt d, dRg)': stats.spearmanr(cu.salt_dtotal, cu.salt_dRg_pct)[0],
                     'rho(salt d, |NCPR|)': stats.spearmanr(cu.salt_dtotal, ncpr.reindex(cu.index))[0]})
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(d, 'analysis_crossmodel.csv'), index=False)
    a, b = 'CALVADOS-2', 'Mpipi-GG'
    js = S[a].merge(S[b], on=['UniProt_ID', 'pos'], suffixes=('_C', '_M'))
    js = js[(js.ess_C > 0.1) & (js.ess_M > 0.1)]
    jc = C[a].join(C[b], lsuffix='_C', rsuffix='_M', how='inner')
    jc = jc[(jc.salt_ess_C > 0.1) & (jc.salt_ess_M > 0.1)]
    lines = [T.round(4).to_string(index=False), '',
             f'site by site ({len(js)} sites usable in both): Spearman {stats.spearmanr(js.phospho_dtotal_C, js.phospho_dtotal_M)[0]:.3f}, '
             f'same sign {(np.sign(js.phospho_dtotal_C) == np.sign(js.phospho_dtotal_M)).mean():.1%}',
             f'salt per chain ({len(jc)} chains usable in both): Spearman {stats.spearmanr(jc.salt_dtotal_C, jc.salt_dtotal_M)[0]:.3f}, '
             f'same sign {(np.sign(jc.salt_dtotal_C) == np.sign(jc.salt_dtotal_M)).mean():.1%}']
    open(os.path.join(d, 'crossmodel_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=f'{TR}/heat_crosses_shared')
    main(ap.parse_args().dir)
