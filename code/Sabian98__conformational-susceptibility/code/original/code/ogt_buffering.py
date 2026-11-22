"""
ogt_buffering.py — do the composition shifts with growth temperature offset each other in their effect on compaction?

Within bacterial phyla (phylum fixed effects, phyla with >= 3 organisms; the set of ogt_followup.py), warmer-growing
organisms have more evenly mixed charges (lower kappa), fewer Asp/Glu/Lys and more pi residues, yet their apparent
scaling exponent nu does not change. Here each shift is turned into its own predicted change of nu:

  1  a linear model of nu on sequence descriptors is fitted to all simulated chains of a force field with a valid
     scaling fit (nu > 0.30); descriptors are used in raw units, so the coefficient beta_j is nu per unit descriptor
  2  the within-phylum slope of each organism-mean descriptor against growth temperature, s_j (per 10 C), is
     estimated with phylum fixed effects
  3  contribution_j = beta_j * s_j (nu per 10 C); their sum is the change of nu that composition predicts, compared
     with the observed within-phylum slope of organism-mean nu (failed fits left out)

Uncertainty: 2,000 bootstrap resamples of organisms for the slopes (beta held at the full-data fit, whose standard
errors are reported and are small with ~10,000 chains). Grouped contributions: charge patterning (kappa), Asp+Glu+Lys,
pi residues (Arg + Phe, Trp, Tyr), net charge (|NCPR|), Pro+Gly+Ala, hydropathy, length.

    BENDER_REANALYSIS_DIR=~/Downloads/Scripts/bender_reanalysis python ogt_buffering.py
Output: Thermal adaptation/bender_ogt/followup_nu_filtered/buffering_{coefficients,contributions}.csv and buffering_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import ogt_followup as OF  # noqa: E402
import ogt_response as OR  # noqa: E402
import total_response_evidence as E  # noqa: E402

OUT = os.path.join(OR.OUT, 'followup_nu_filtered')
X = ['kappa', 'f_DEK', 'f_R', 'f_aro', 'abs_NCPR', 'f_P', 'f_G', 'f_A', 'hydropathy_avg', 'log_length']
GROUPS = {'charge patterning (kappa)': ['kappa'], 'Asp + Glu + Lys': ['f_DEK'], 'pi residues (Arg, Phe, Trp, Tyr)': ['f_R', 'f_aro'],
          'net charge (|NCPR|)': ['abs_NCPR'], 'Pro + Gly + Ala': ['f_P', 'f_G', 'f_A'], 'hydropathy': ['hydropathy_avg'],
          'length': ['log_length']}
NU_MIN = OF.NU_MIN


def ols(d, y, cols):
    A = np.column_stack([np.ones(len(d))] + [d[c].values for c in cols])
    b, *_ = np.linalg.lstsq(A, d[y].values, rcond=None)
    res = d[y].values - A @ b
    s2 = res @ res / (len(d) - A.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(A.T @ A)))
    r2 = 1 - res @ res / ((d[y] - d[y].mean()) ** 2).sum()
    return dict(zip(['const'] + cols, b)), dict(zip(['const'] + cols, se)), r2


def main():
    rng = np.random.default_rng(11)
    D, _ = E.load()
    bio = E.bio_table()
    bio['log_length'] = np.log(bio.length)
    bio['OGT'] = OR.ogt_table(bio.reset_index().set_index('UniProt_ID'))
    bio['phylum'] = OR.ogt_table(bio.reset_index().set_index('UniProt_ID'), 'phylum')
    coef_rows, con_rows, lines = [], [], []
    for k in ('C', 'R'):
        model = E.MODELS[k]
        d = D[k][['nu']].join(bio[['organism', 'kingdom', 'phylum', 'OGT'] + X], how='inner')
        d = d[d.nu > NU_MIN].dropna(subset=X)
        beta, se, r2 = ols(d, 'nu', X)                       # step 1, all chains with a valid fit
        for c in X:
            coef_rows.append({'model': model, 'descriptor': c, 'beta': beta[c], 'se': se[c], 'chains': len(d), 'R2': r2})
        # organism means within bacterial phyla (as ogt_followup.py)
        sub = d[d.kingdom == 'Bacteria'].dropna(subset=['phylum', 'OGT'])
        o = sub.groupby(['organism', 'phylum']).agg(OGT=('OGT', 'first'), nu=('nu', 'mean'),
                                                     **{c: (c, 'mean') for c in X}).reset_index()
        o = o[o.phylum.isin(o.phylum.value_counts()[lambda s: s >= 3].index)].reset_index(drop=True)
        o['nu_pred'] = beta['const'] + sum(beta[c] * o[c] for c in X)

        def slopes(oo):
            s = {c: OF.fe_slope(oo, c, 'phylum')[0] for c in X + ['nu', 'nu_pred']}
            return s
        est = slopes(o)
        boots = [slopes(o.iloc[rng.integers(0, len(o), len(o))].reset_index(drop=True)) for _ in range(OF.N_BOOT)]
        B = pd.DataFrame(boots)
        contrib = lambda S: {g: sum(beta[c] * S[c] for c in cs) for g, cs in GROUPS.items()}   # noqa: E731
        ce = contrib(est)
        cb = pd.DataFrame([contrib(B.iloc[i]) for i in range(len(B))])
        tot_e, tot_b = sum(ce.values()), cb.sum(axis=1)
        for g in GROUPS:
            con_rows.append({'model': model, 'term': g, 'nu_per_10C': ce[g], 'ci_lo': np.percentile(cb[g], 2.5),
                             'ci_hi': np.percentile(cb[g], 97.5), 'organisms': len(o), 'phyla': o.phylum.nunique()})
        oth_e = tot_e - ce['pi residues (Arg, Phe, Trp, Tyr)']
        oth_b = tot_b - cb['pi residues (Arg, Phe, Trp, Tyr)']
        chg_e = ce['charge patterning (kappa)'] + ce['Asp + Glu + Lys'] + ce['net charge (|NCPR|)']
        chg_b = cb['charge patterning (kappa)'] + cb['Asp + Glu + Lys'] + cb['net charge (|NCPR|)']
        for term, e_, b_ in (('all charge terms (kappa, Asp+Glu+Lys, |NCPR|)', chg_e, chg_b), ('all shifts except pi', oth_e, oth_b)):
            con_rows.append({'model': model, 'term': term, 'nu_per_10C': e_, 'ci_lo': np.percentile(b_, 2.5),
                             'ci_hi': np.percentile(b_, 97.5), 'organisms': len(o), 'phyla': o.phylum.nunique(),
                             'frac_boot_positive': float((b_ > 0).mean())})
        con_rows.append({'model': model, 'term': 'sum of contributions (predicted)', 'nu_per_10C': tot_e,
                         'ci_lo': np.percentile(tot_b, 2.5), 'ci_hi': np.percentile(tot_b, 97.5), 'organisms': len(o),
                         'phyla': o.phylum.nunique()})
        con_rows.append({'model': model, 'term': 'observed nu', 'nu_per_10C': est['nu'],
                         'ci_lo': np.percentile(B.nu, 2.5), 'ci_hi': np.percentile(B.nu, 97.5), 'organisms': len(o),
                         'phyla': o.phylum.nunique()})
        r_org = np.corrcoef(o.nu, o.nu_pred)[0, 1]
        lines.append(f'{model}: chain model R2 {r2:.3f} on {len(d)} chains; organism-mean nu vs predicted r = {r_org:.3f} '
                     f'({len(o)} organisms, {o.phylum.nunique()} phyla); check sum(beta*s) = slope of nu_pred: '
                     f'{tot_e:.5f} vs {est["nu_pred"]:.5f}')
    C, K = pd.DataFrame(coef_rows), pd.DataFrame(con_rows)
    os.makedirs(OUT, exist_ok=True)
    C.to_csv(os.path.join(OUT, 'buffering_coefficients.csv'), index=False)
    K.to_csv(os.path.join(OUT, 'buffering_contributions.csv'), index=False)
    pd.set_option('display.width', 200)
    txt = '\n'.join(lines + ['', C.round(5).to_string(index=False), '', K.round(5).to_string(index=False)])
    open(os.path.join(OUT, 'buffering_summary.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
