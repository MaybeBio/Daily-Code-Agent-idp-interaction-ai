"""
llps_unevenness_diagnostic.py — DIAGNOSTIC (post hoc, after the pre-registered test of llps_unevenness_plan.md failed):
the exploratory drivers of llps_profiles.py (50-800 residues) through the identical confirmatory pipeline (5 matched
controls each, CRC32 seeds, local shuffles), to see whether matching or the driver set explains the failure.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python llps_unevenness_diagnostic.py  ->  llps_csat/llps_unevenness_diagnostic.txt
"""
import csv
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llps_unevenness as U  # noqa: E402

NAMES = ['P35637_1_285', 'P09651_187_372', 'P22626_194_353', 'Q13148_252_414', 'Q9NQI0_1_252', 'Q92804_1_234',
         'Q01844_1_360', 'P31483_283_386', 'Q9UHD9_98_581', 'Q13283_137_342', 'P06748_122_247', 'Q14444_247_709',
         'P52948_1_166', 'O60885_671_1230', 'P10636_1_758', 'P15502_1_786', 'Q06787_442_632', 'Q96PK6_145_669',
         'P23246_1_285', 'Q13283_410_466']


def main():
    idrome = [r for r in csv.DictReader(open(U.IDROME)) if 50 <= int(r['N']) <= 800]
    byn = {r['seq_name']: r for r in idrome}
    drv = [byn[n] for n in NAMES if n in byn]
    acc_pp = {e['accession'] for e in json.load(open(U.PHASEPRO)).values()}
    pool = [r for r in idrome if r['UniProt_ID'] not in acc_pp and r['UniProt_ID'] not in U.EXPLORATORY]
    Xp = np.array([U.comp_features(r['fasta']) for r in pool])
    mu, sd = Xp.mean(0), Xp.std(0)
    Zp = (Xp - mu) / sd
    used, match = np.zeros(len(pool), bool), {}
    for k in np.random.default_rng(0).permutation(len(drv)):
        dist = np.linalg.norm(Zp - (U.comp_features(drv[k]['fasta']) - mu) / sd, axis=1)
        dist[used] = np.inf
        idx = np.argsort(dist)[:5]
        used[idx] = True
        match[drv[k]['seq_name']] = [pool[i] for i in idx]
    jobs = [(r['seq_name'], r['fasta'], False) for r in drv] + [(c['seq_name'], c['fasta'], False)
                                                              for v in match.values() for c in v]
    with Pool(8, initializer=U.init) as pl:
        res = {o['name']: o for o in pl.imap_unordered(U.one, jobs, chunksize=2)}
    out = [f'DIAGNOSTIC (post hoc): {len(drv)} exploratory drivers (50-800 residues) vs 5 matched controls each']
    for key, lab in (('pct_local', 'local shuffles'), ('pct_full', 'full shuffles')):
        dv = np.array([res[r['seq_name']][key] for r in drv])
        cv_ = np.array([res[c['seq_name']][key] for v in match.values() for c in v])
        d = np.array([res[r['seq_name']][key] - np.mean([res[c['seq_name']][key] for c in match[r['seq_name']]]) for r in drv])
        out.append(f'  {lab}: drivers median {np.median(dv):.0f}, matched controls median {np.median(cv_):.0f}; difference '
                   f'median {np.median(d):+.1f}, positive {int((d > 0).sum())}/{len(d)}, Wilcoxon P = '
                   f'{wilcoxon(d, alternative="greater").pvalue:.2g}')
    print('\n'.join(out))
    open(os.path.join(U.OUT, 'llps_unevenness_diagnostic.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
