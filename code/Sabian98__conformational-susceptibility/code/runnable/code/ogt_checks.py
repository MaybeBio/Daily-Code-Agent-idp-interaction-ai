"""
ogt_checks.py — two checks of the BENDER growth-temperature result (ogt_followup.py: warm-growing bacteria's IDPs are more
expanded at 300 K through better-mixed charges, lower kappa, with fewer D+E+K and more pi/Arg).

1  Independent replication on the 43,104 disordered regions of 44 reference proteomes (thermal_predict.py): sequence
   descriptors (kappa, FCR, |NCPR|, D+E+K, pi, Arg, Ser; same code as thermal_predict.descriptors) per organism against
   OGT, clade fixed effects (organism bootstrap); and within ortholog-matched families (same eggNOG group and position
   within a clade), thermophile minus mesophile, Wilcoxon over families.
2  Finer phylogeny in BENDER: the OGT slopes with genus fixed effects (genera with >= 2 organisms that differ in OGT),
   for nu, swelling alpha at 300 K and projected to the OGT (as ogt_followup.py), kappa and the composition descriptors.
    python ogt_checks.py [--out DIR]
Output: Thermal adaptation/bender_ogt/checks_*.csv, checks_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import ogt_response as OR  # noqa: E402
import thermal_predict as TP  # noqa: E402
import total_response_evidence as E  # noqa: E402

DESC = ['kappa', 'FCR', 'absNCPR', 'f_DEK', 'f_pi', 'f_R', 'f_K', 'f_E', 'f_D', 'f_S', 'f_P', 'f_G']
N_BOOT = 2000


def fe_slope(o, y, fe, rng):
    o = o.dropna(subset=[y])
    def fit(d):
        X = np.column_stack([np.ones(len(d)), d.OGT, pd.get_dummies(d[fe], drop_first=True).astype(float)])
        return 10 * np.linalg.lstsq(X, d[y], rcond=None)[0][1]
    b = [fit(o.iloc[rng.integers(0, len(o), len(o))]) for _ in range(N_BOOT)]
    return fit(o), np.percentile(b, 2.5), np.percentile(b, 97.5)


def check_regions(rng):
    R = pd.read_csv(f'{TP.OUT}/regions.csv')
    seq = pd.read_csv(f'{TP.OUT}/region_sequences.csv').set_index('region_id').sequence
    X = TP.descriptors(seq.reindex(R.region_id).tolist())
    for c in DESC:
        R[c] = X[c].to_numpy()
    R = R[R.in_range]
    sd = R[DESC].std()
    o = R.groupby(['upid', 'organism', 'clade']).agg(OGT=('OGT_C', 'first'), **{c: (c, 'mean') for c in DESC}).reset_index()
    o = o[~o.clade.str.contains('no mesophile')]
    rows = []
    for c in DESC:   # fixed 2026-10-10: SD of organism means (as ogt_followup / check_genus), not the region-level SD
        z = o.assign(zc=o[c] / o[c].std())
        est, lo, hi = fe_slope(z, 'zc', 'clade', rng)
        rows.append({'check': 'proteome regions, organism means, clade fixed effects', 'term': c, 'n': len(o),
                     'slope_SD_per_10C': est, 'ci_lo': lo, 'ci_hi': hi})
    # ortholog-matched
    Ex = R.assign(og=R.eggnog.fillna('').str.split(';')).explode('og')
    Ex = Ex[Ex.og.str.len() > 0]
    diffs = {c: [] for c in DESC}
    for (clade, og, pos), g in Ex.groupby(['clade', 'og', 'position']):
        t, m = g[g['class'] == 'thermophile'], g[g['class'] == 'mesophile']
        if t.upid.nunique() and m.upid.nunique():
            for c in DESC:
                diffs[c].append(t.groupby('upid')[c].mean().mean() - m.groupby('upid')[c].mean().mean())
    for c in DESC:
        d = np.array(diffs[c]) / sd[c]
        rows.append({'check': 'ortholog-matched families, thermophile minus mesophile (region SD)', 'term': c, 'n': len(d),
                     'slope_SD_per_10C': np.median(d), 'ci_lo': np.nan, 'ci_hi': np.nan,
                     'frac_positive': (d > 0).mean(), 'wilcoxon_p': stats.wilcoxon(d).pvalue})
    # bacteria only (fungi dominate the matched families)
    for c in DESC:
        dd = []
        for (clade, og, pos), g in Ex[~Ex.clade.str.startswith('Fungi')].groupby(['clade', 'og', 'position']):
            t, m = g[g['class'] == 'thermophile'], g[g['class'] == 'mesophile']
            if t.upid.nunique() and m.upid.nunique():
                dd.append(t.groupby('upid')[c].mean().mean() - m.groupby('upid')[c].mean().mean())
        d = np.array(dd) / sd[c]
        rows.append({'check': 'ortholog-matched families, bacteria only (SD)', 'term': c, 'n': len(d),
                     'slope_SD_per_10C': np.median(d), 'ci_lo': np.nan, 'ci_hi': np.nan,
                     'frac_positive': (d > 0).mean(), 'wilcoxon_p': stats.wilcoxon(d).pvalue if len(d) > 5 else np.nan})
    return pd.DataFrame(rows)


def check_genus(rng):
    D, _ = E.load()
    bio = E.bio_table()
    bio['OGT'] = OR.ogt_table(bio)
    bio['genus'] = OR.ogt_table(bio, 'genus')
    bio['absNCPR'] = bio.NCPR.abs()
    rows = []
    for k in ('C', 'R'):
        d = D[k][['total', 'nu', 'rg', 'N']].join(bio[['organism', 'kingdom', 'genus', 'OGT', 'kappa', 'FCR', 'absNCPR',
                                                     'f_DEK', 'f_pi', 'f_R', 'f_S']], how='inner').dropna(subset=['OGT', 'genus'])
        d.loc[d.nu <= OR.NU_MIN, 'nu'] = np.nan     # failed scaling fits, as ogt_followup.py (added 2026-10-20)
        d['alpha_300'] = np.log(d.rg / np.sqrt(d.N))
        d['alpha_OGT'] = d.alpha_300 + np.log1p(d.total / 100 * (d.OGT + 273.15 - 300.0))
        cols = ['nu', 'alpha_300', 'alpha_OGT', 'kappa', 'FCR', 'absNCPR', 'f_DEK', 'f_pi', 'f_R', 'f_S']
        o = d.groupby(['organism', 'genus', 'kingdom']).agg(OGT=('OGT', 'first'), **{c: (c, 'mean') for c in cols}).reset_index()
        g = o.groupby('genus').OGT.agg(['size', 'std'])
        keep = g[(g['size'] >= 2) & (g['std'] > 0)].index
        o = o[o.genus.isin(keep)]
        for c in cols:
            z = o.assign(zc=o[c] / (o[c].std() if c not in ('nu', 'alpha_300', 'alpha_OGT') else 1.0))
            est, lo, hi = fe_slope(z, 'zc', 'genus', rng)
            rows.append({'check': f'BENDER {E.MODELS[k]}, genus fixed effects', 'term': c, 'n': len(o),
                         'genera': o.genus.nunique(), 'OGT_range_within_genera_median': float(g.loc[keep, 'std'].median()),
                         'slope_per_10C (SD for descriptors)': est, 'ci_lo': lo, 'ci_hi': hi})
    return pd.DataFrame(rows)


def main(out):
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(4)
    A, B = check_regions(rng), check_genus(rng)
    A.to_csv(os.path.join(out, 'checks_regions.csv'), index=False)
    B.to_csv(os.path.join(out, 'checks_genus.csv'), index=False)
    pd.set_option('display.width', 220)
    lines = ['== 1 independent replication: 44 proteomes, disordered regions', A.round(4).to_string(index=False), '',
             '== 2 BENDER within genera', B.round(4).to_string(index=False)]
    open(os.path.join(out, 'checks_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OR.OUT)
    main(ap.parse_args().out)
