"""
ogt_followup.py — two follow-ups to ogt_response.py (BENDER chains with an optimal growth temperature, OGT).

1  Corresponding states: are IDPs from warm- and cold-growing organisms alike in size where they live?
   Compaction is measured as swelling, alpha = ln(<Rg> / N^0.5) (no power-law fit needed). At 300 K from the simulations;
   at the organism's OGT projected with each chain's own 300 K slope s (% of <Rg> per K):
       alpha(OGT) = alpha(300) + ln(1 + s/100 (T_OGT - 300))
   with s = the calibrated total response, or fixed-water chi alone for comparison. The OGT slope of organism-mean alpha
   (bacteria, phylum fixed effects; all kingdoms, kingdom fixed effects) at 300 K vs at OGT: a slope that shrinks to 0
   at OGT means size is held constant across thermal niches. Linear projection from 300 K: beyond 320 K it extrapolates
   past the solvation calibration (heat-capacity terms are not constrained).
2  What makes warm-growing organisms' IDPs more expanded: OGT slope (organism means, same fixed effects) of sequence
   descriptors, and how much of the OGT slope of nu remains after adding composition covariates.
    python ogt_followup.py [--out DIR]
Output: Thermal adaptation/bender_ogt/followup_*.csv, followup_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import ogt_response as OR  # noqa: E402
import total_response_evidence as E  # noqa: E402

DESC = ['FCR', 'abs_NCPR', 'NCPR', 'f_pi', 'f_DEK', 'f_R', 'f_K', 'f_E', 'f_D', 'f_P', 'f_G', 'f_A', 'f_S', 'f_T', 'f_Q',
        'f_N', 'hydropathy_avg', 'kappa', 'SCD', 'log_length']
COVS = {'+ length': ['log_length'], '+ charge (FCR, |NCPR|, kappa)': ['log_length', 'FCR', 'abs_NCPR', 'kappa'],
        '+ charge + P, G, A': ['log_length', 'FCR', 'abs_NCPR', 'kappa', 'f_P', 'f_G', 'f_A'],
        '+ all composition': ['log_length', 'FCR', 'abs_NCPR', 'kappa', 'f_P', 'f_G', 'f_A', 'f_pi', 'hydropathy_avg',
                              'f_S', 'f_T', 'f_Q', 'f_N']}
N_BOOT = 2000
NU_MIN = 0.30   # added 2026-10-19: nu at or below 0.30 marks a failed fit, as everywhere in the paper


def fe_slope(o, y, fe, extra=(), rng=None):
    o = o.dropna(subset=[y] + list(extra))

    def fit(d):
        X = np.column_stack([np.ones(len(d)), d.OGT] + [d[c] for c in extra] + [pd.get_dummies(d[fe], drop_first=True).astype(float)])
        return 10 * np.linalg.lstsq(X, d[y], rcond=None)[0][1]
    est = fit(o)
    if rng is None:
        return est, np.nan, np.nan
    b = [fit(o.iloc[rng.integers(0, len(o), len(o))]) for _ in range(N_BOOT)]
    return est, np.percentile(b, 2.5), np.percentile(b, 97.5)


def main(out):
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(3)
    D, _ = E.load()
    bio = E.bio_table()
    bio['log_length'] = np.log(bio.length)
    bio['OGT'] = OR.ogt_table(bio.reset_index().set_index('UniProt_ID'))
    bio['phylum'] = OR.ogt_table(bio.reset_index().set_index('UniProt_ID'), 'phylum')
    rows1, rows2 = [], []
    for k in ('C', 'R'):
        d = D[k][['total', 'chi', 'nu', 'rg', 'N']].join(bio[['organism', 'kingdom', 'phylum', 'OGT'] + DESC], how='inner')
        d = d.dropna(subset=['OGT'])
        d.loc[d.nu <= NU_MIN, 'nu'] = np.nan     # failed scaling fits (Supplementary Section S1), left out of nu means
        dT = d.OGT + 273.15 - 300.0
        d['alpha_300'] = np.log(d.rg / np.sqrt(d.N))
        d['alpha_OGT_total'] = d.alpha_300 + np.log1p(d.total / 100 * dT)
        d['alpha_OGT_chi'] = d.alpha_300 + np.log1p(d.chi / 100 * dT)
        d['dRg_pct_at_OGT'] = 100 * (np.exp(d.alpha_OGT_total - d.alpha_300) - 1)
        agg = {c: (c, 'mean') for c in ['alpha_300', 'alpha_OGT_total', 'alpha_OGT_chi', 'dRg_pct_at_OGT', 'nu'] + DESC}
        for label, sub, fe in (('bacteria, phylum fixed effects', d[d.kingdom == 'Bacteria'].dropna(subset=['phylum']), 'phylum'),
                               ('all kingdoms, kingdom fixed effects', d, 'kingdom')):
            o = sub.groupby(['organism', fe]).agg(OGT=('OGT', 'first'), chains=('total', 'size'), **agg).reset_index()
            o = o[o[fe].isin(o[fe].value_counts()[lambda s: s >= 3].index)]
            for y in ('alpha_300', 'alpha_OGT_total', 'alpha_OGT_chi', 'dRg_pct_at_OGT', 'nu'):
                est, lo, hi = fe_slope(o, y, fe, rng=rng)
                rows1.append({'model': E.MODELS[k], 'set': label, 'organisms': len(o), 'quantity': y,
                              'slope_per_10C': est, 'ci_lo': lo, 'ci_hi': hi})
            if fe == 'phylum':
                for c in DESC + ['nu']:
                    z = o.assign(zc=(o[c] - o[c].mean()) / o[c].std())
                    est, lo, hi = fe_slope(z, 'zc', fe, rng=rng)
                    rows2.append({'model': E.MODELS[k], 'analysis': 'descriptor vs OGT (SD per 10 C)', 'term': c,
                                  'organisms': len(o), 'slope_per_10C': est, 'ci_lo': lo, 'ci_hi': hi})
                base = fe_slope(o, 'nu', fe)[0]
                for name, cov in COVS.items():
                    est, lo, hi = fe_slope(o, 'nu', fe, cov, rng=rng)
                    rows2.append({'model': E.MODELS[k], 'analysis': 'nu vs OGT after covariates', 'term': name,
                                  'organisms': len(o), 'slope_per_10C': est, 'ci_lo': lo, 'ci_hi': hi,
                                  'share_of_slope_remaining': est / base})
    A, B = pd.DataFrame(rows1), pd.DataFrame(rows2)
    A.to_csv(os.path.join(out, 'followup_corresponding_states.csv'), index=False)
    B.to_csv(os.path.join(out, 'followup_composition.csv'), index=False)
    pd.set_option('display.width', 220)
    lines = ['== 1 swelling alpha = ln(Rg / N^0.5) of organism means vs OGT, at 300 K and projected to the OGT',
             A.round(4).to_string(index=False), '', '== 2 what changes with OGT (bacteria, phylum fixed effects)',
             B.sort_values(['analysis', 'model', 'slope_per_10C']).round(4).to_string(index=False)]
    open(os.path.join(out, 'followup_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OR.OUT)
    main(ap.parse_args().out)
