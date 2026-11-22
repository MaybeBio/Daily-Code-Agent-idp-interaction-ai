"""
kappa_disorder_check.py — does the charge-mixing (kappa) trend with growth temperature differ between proteins that are
disordered over (almost) their whole length and disordered regions of mostly folded proteins? Exploratory, after the fact.

BENDER (fully disordered proteins) gives kappa falling with growth temperature within bacterial phyla; the 44-proteome
replication (MobiDB-lite regions, thermal_predict.py / ogt_checks.check_regions) does not. Here the replication is
repeated on subsets of the same regions by the share of their protein that the region covers, with the same descriptors
(thermal_predict.descriptors), the same organism means, clade fixed effects and organism bootstrap (ogt_checks.fe_slope).
    python kappa_disorder_check.py
Output: Thermal adaptation/predictions/kappa_by_coverage.csv, kappa_by_coverage.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import ogt_checks as OC  # noqa: E402
import thermal_predict as TP  # noqa: E402

MIN_REGIONS = 5   # organisms with fewer regions in a subset are left out of that subset's organism means


def main():
    R = pd.read_csv(f'{TP.OUT}/regions.csv')
    seq = pd.read_csv(f'{TP.OUT}/region_sequences.csv').set_index('region_id').sequence
    R['kappa'] = TP.descriptors(seq.reindex(R.region_id).tolist())['kappa'].to_numpy()
    R['kappa_raw'] = [TP.kappa(s) for s in seq.reindex(R.region_id)]
    R['coverage'] = (R.end - R.start + 1) / R.protein_length
    R = R[R.in_range]
    subsets = {'all regions (as published)': R.coverage >= 0, 'coverage >= 0.9': R.coverage >= 0.9,
               'coverage >= 0.7': R.coverage >= 0.7, 'coverage < 0.5': R.coverage < 0.5}
    rows, lines = [], ['kappa slope with growth temperature, SD of organism means per 10 C, clade fixed effects, organism bootstrap', '']
    for name, m in subsets.items():
        for col in ('kappa', 'kappa_raw'):
            S = R[m].dropna(subset=[col])
            n_reg = S.groupby('upid').size()
            S = S[S.upid.isin(n_reg[n_reg >= (1 if name.startswith('all') else MIN_REGIONS)].index)]
            o = S.groupby(['upid', 'organism', 'clade']).agg(OGT=('OGT_C', 'first'), k=(col, 'mean')).reset_index()
            o = o[~o.clade.str.contains('no mesophile')]
            z = o.assign(zc=o.k / o.k.std())
            est, lo, hi = OC.fe_slope(z, 'zc', 'clade', np.random.default_rng(1))
            r = {'subset': name, 'kappa': 'median-filled (as published)' if col == 'kappa' else 'undefined left out',
                 'regions': len(S), 'organisms': len(o), 'slope': est, 'ci_lo': lo, 'ci_hi': hi}
            rows.append(r)
            lines.append(f"{name:28s} {r['kappa']:28s} {r['regions']:6d} regions, {r['organisms']:2d} organisms: "
                         f"{est:+.2f} ({lo:+.2f}, {hi:+.2f})")
    pd.DataFrame(rows).to_csv(f'{TP.OUT}/kappa_by_coverage.csv', index=False)
    open(f'{TP.OUT}/kappa_by_coverage.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
