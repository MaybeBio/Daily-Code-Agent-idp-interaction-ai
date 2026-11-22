"""
phospho_direct_ceiling.py — noise ceiling of the direct check of phosphorylation (2026-10-25). How large a rank
correlation between the reweighted and the directly simulated change of the total response could the 60 sites give if
the reweighting were exact, given the sampling errors of both?
  true change of each site: the predicted change shrunk towards the mean by the square root of its reliability
  (variance of the predictions minus their mean squared standard error, over the variance of the predictions)
  each draw adds normal noise with the standard error of each site to the true change, once for the prediction and once
  for the direct change, and records the Spearman correlation of the two (20,000 draws, seed 0)
    python phospho_direct_ceiling.py      (openmpipi env)
Input:  heldout_temperature/phospho_heat/capped/sites_direct.csv (phospho_heat_direct_capped.py)
Output: heldout_temperature/phospho_heat/capped/noise_ceiling.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

BM = '/Volumes/User Homes/Nu Project/BENDER-MT'
D = os.path.join(BM, 'heldout_temperature', 'phospho_heat', 'capped')
NDRAW = 20000


def main():
    d = pd.read_csv(os.path.join(D, 'sites_direct.csv'))
    p, ps = d.pred_dtotal.values, d.pred_se.values
    y, ys = d.direct_dtotal.values, d.direct_dtotal_se.values
    obs = stats.spearmanr(p, y).statistic
    rel = (p.var(ddof=1) - (ps ** 2).mean()) / p.var(ddof=1)
    if not 0 < rel <= 1:
        raise SystemExit(f'reliability of the predictions out of range ({rel:.3f})')
    true = p.mean() + (p - p.mean()) * np.sqrt(rel)
    rng = np.random.default_rng(0)
    sims = np.array([stats.spearmanr(true + rng.normal(0, ps), true + rng.normal(0, ys)).statistic
                     for _ in range(NDRAW)])
    lo, med, hi = np.percentile(sims, [5, 50, 95])
    L = ['Noise ceiling of the direct check of phosphorylation (CALVADOS-2, calibrated model, % per K)',
         f'sites {len(d)}; median |predicted| {np.median(np.abs(p)):.3f}, median |direct| {np.median(np.abs(y)):.3f}, '
         f'rms SE predicted {np.sqrt((ps ** 2).mean()):.3f}, direct {np.sqrt((ys ** 2).mean()):.3f}',
         f'reliability of the predictions {rel:.2f}',
         f'observed Spearman {obs:.3f}',
         f'expected Spearman if the reweighting were exact: median {med:.2f}, 90% range {lo:.2f} to {hi:.2f} ({NDRAW} draws)',
         f'share of draws at or below the observed value {(sims <= obs).mean():.3f}']
    open(os.path.join(D, 'noise_ceiling.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
