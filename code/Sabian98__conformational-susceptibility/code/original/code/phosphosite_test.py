"""
Are annotated phosphosites the Ser/Thr whose phosphorylation moves the chain most?

Per-site effects: salt_ph_phospho.analyse at 300 K (charge model: dq = -2 in CALVADOS-2, -1.5 in Mpipi-GG; exact
reweighting dRg_pct_reweight). Annotations: UniProt MOD_RES Phosphoserine/Phosphothreonine for the BENDER accessions
(Analysis_FDT/phosphosite_test/phosphosites_annotated.csv, positions mapped onto the simulated sequences; evidence
'experimental' = ECO:0000269 or ECO:0007744, 'similarity' = ECO:0000250).

Test, within each protein (so protein-level differences cannot drive it): percentile rank of a site's |dRg| among
that protein's Ser/Thr (0 = smallest, 1 = largest; ties averaged); mean percentile of annotated sites against 0.5,
with a permutation null that shuffles the labels among Ser/Thr within each protein. Also the signed effect.
Controls: conditional logistic regression (groups = proteins) of annotated ~ within-protein z(|dRg|) + local K/R count
(+-5, kinase basophilic motifs and tryptic-peptide detectability both favour sites near K/R) + local D/E count (+-5)
+ Pro at +1 (proline-directed kinases) + distance to the nearest terminus / N.
Proteins: every protein with at least one mapped annotated site and at least one other Ser/Thr; CALVADOS-2 for all,
Mpipi-GG where the protein is in the common set.
Output: Analysis_FDT/phosphosite_test/{sites.csv, tests.csv}
Sites whose reweighting keeps a Kish effective sample fraction of 0.1 or less are left out (ESS_MIN; 2026-10-20).
Usage: python3 phosphosite_test.py [--workers 8] [--reuse-sites]
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import salt_ph_phospho as sp  # noqa: E402

OUT = os.path.join(sp.L.OUT_ROOT, 'phosphosite_test')
N_PERM = 10000
ESS_MIN = 0.1   # reliable reweighting, as everywhere in the paper (added 2026-10-20; before, 12 of 2,266 Ser/Thr were below)


def site_effects(workers):
    ann = pd.read_csv(os.path.join(OUT, 'phosphosites_annotated.csv'))
    ann = ann[ann.mapped]
    ids = sorted(ann.UniProt_ID.unique())
    have = pd.read_csv(os.path.join(sp.OUT_DIR, 'phospho_sites_common.csv'))
    kingdom = pd.read_csv(sp.bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    todo = [(u, 'CALVADOS-2', kingdom.get(u)) for u in ids if u not in set(have[have.model == 'CALVADOS-2'].UniProt_ID)]
    rows = []
    if todo:
        with mp.Pool(workers) as pool:
            for _, s, err in pool.imap_unordered(sp.analyse, todo):
                if err:
                    print('  failed:', err[:200], flush=True)
                else:
                    rows += s
    S = pd.concat([have[have.UniProt_ID.isin(ids)], pd.DataFrame(rows)], ignore_index=True)
    lab = ann.groupby(['UniProt_ID', 'pos']).evidence.agg(lambda e: 'experimental' if 'experimental' in set(e) else sorted(set(e))[0])
    S = S.join(lab.rename('evidence'), on=['UniProt_ID', 'pos'])
    S['annotated'] = S.evidence.notna()
    return S.sort_values(['model', 'UniProt_ID', 'pos']).reset_index(drop=True)  # same row order whichever worker finishes first


def context(S):
    seq = pd.read_csv(sp.bc.BIO_CSV, usecols=['UniProt_ID', 'sequence']).set_index('UniProt_ID')['sequence']
    kr, de, pro, term = [], [], [], []
    for u, p in zip(S.UniProt_ID, S.pos):
        q = seq[u]
        w = q[max(0, p - 6):p - 1] + q[p:p + 5]
        kr.append(sum(w.count(a) for a in 'KR'))
        de.append(sum(w.count(a) for a in 'DE'))
        pro.append(int(p < len(q) and q[p] == 'P'))
        term.append(min(p - 1, len(q) - p) / len(q))
    return S.assign(local_KR=kr, local_DE=de, pro_plus1=pro, term_dist=term)


def within_protein(S):
    S = S.copy()
    S['abs_d'] = S.dRg_pct_reweight.abs()
    g = S.groupby(['model', 'UniProt_ID'])
    S['pct'] = g.abs_d.rank(pct=False, method='average').sub(1).div(g.abs_d.transform('size').sub(1).clip(lower=1))
    S['z_abs'] = (S.abs_d - g.abs_d.transform('mean')) / g.abs_d.transform('std').replace(0, np.nan)
    keep = g.annotated.transform('any') & ~g.annotated.transform('all')
    return S[keep]


def perm_test(S, mask_ann, rng):
    """Mean within-protein percentile of annotated sites vs labels shuffled within proteins."""
    obs = S.loc[mask_ann, 'pct'].mean()
    # .indices already holds row positions (np.flatnonzero on them gave 0..n-1, i.e. the first rows of the table)
    groups = [(idx, int(mask_ann.iloc[idx].sum())) for idx in S.groupby('UniProt_ID').indices.values()]
    pct = S.pct.to_numpy()
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        tot, n = 0.0, 0
        for idx, na in groups:
            if na:
                tot += pct[rng.choice(idx, na, replace=False)].sum()
                n += na
        null[k] = tot / n
    return obs, float((np.sum(null >= obs) + 1) / (N_PERM + 1)), float(null.mean())


def main(workers, reuse_sites=False):
    rng = np.random.default_rng(0)
    if reuse_sites:                                   # tests only, on the stored per-site effects (no reweighting)
        S = pd.read_csv(os.path.join(OUT, 'sites.csv'))
    else:
        S = context(site_effects(workers))
        S.to_csv(os.path.join(OUT, 'sites.csv'), index=False)
    out = []
    for model, x in S.groupby('model'):
        out.append({'model': model, 'test': f'Ser/Thr with effective sample <= {ESS_MIN} (left out)', 'value': str(int((x.ess_frac <= ESS_MIN).sum()))})
        W = within_protein(x[x.ess_frac > ESS_MIN])
        out.append({'model': model, 'test': 'proteins / sites / annotated', 'value': f'{W.UniProt_ID.nunique()} / {len(W)} / {int(W.annotated.sum())}'})
        for name, m in (('all annotated', W.annotated), ('experimental only', W.evidence.eq('experimental'))):
            if m.sum() < 5:
                continue
            obs, p, null = perm_test(W.reset_index(drop=True), m.reset_index(drop=True), rng)
            out.append({'model': model, 'test': f'{name}: mean within-protein percentile of |dRg| (null {null:.3f})', 'value': f'{obs:.3f}', 'p': p})
        a, b = W.loc[W.annotated, 'dRg_pct_reweight'], W.loc[~W.annotated, 'dRg_pct_reweight']
        out.append({'model': model, 'test': 'median dRg %, annotated vs other Ser/Thr', 'value': f'{a.median():.2f} vs {b.median():.2f}',
                    'p': stats.mannwhitneyu(a, b).pvalue})
        out.append({'model': model, 'test': 'median |dRg| %, annotated vs other', 'value': f'{a.abs().median():.2f} vs {b.abs().median():.2f}'})
        for c in ('local_KR', 'local_DE', 'pro_plus1', 'term_dist'):
            out.append({'model': model, 'test': f'mean {c}, annotated vs other', 'value': f'{W.loc[W.annotated, c].mean():.3f} vs {W.loc[~W.annotated, c].mean():.3f}',
                        'p': stats.mannwhitneyu(W.loc[W.annotated, c], W.loc[~W.annotated, c]).pvalue})
            out.append({'model': model, 'test': f'Spearman |dRg| vs {c} (all Ser/Thr)', 'value': f'{stats.spearmanr(W.abs_d, W[c])[0]:.3f}'})
        try:
            from statsmodels.discrete.conditional_models import ConditionalLogit
            Wz = W.dropna(subset=['z_abs'])
            for cols in (['z_abs'], ['z_abs', 'local_KR', 'local_DE', 'pro_plus1', 'term_dist']):
                X = Wz[cols].astype(float)
                X = (X - X.mean()) / X.std()
                r = ConditionalLogit(Wz.annotated.astype(int).to_numpy(), X.to_numpy(), groups=Wz.UniProt_ID.to_numpy()).fit(disp=False)
                for k, c in enumerate(cols):
                    out.append({'model': model, 'test': f'conditional logit [{" + ".join(cols)}]: {c} (per SD)',
                                'value': f'OR {np.exp(r.params[k]):.2f} [{np.exp(r.conf_int()[k, 0]):.2f}, {np.exp(r.conf_int()[k, 1]):.2f}]',
                                'p': float(r.pvalues[k])})
        except Exception as e:  # noqa: BLE001
            out.append({'model': model, 'test': f'conditional logit failed: {e!r}', 'value': ''})
    T = pd.DataFrame(out)
    T.to_csv(os.path.join(OUT, 'tests.csv'), index=False)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 110)
    print(T.to_string(index=False))
    sp.bc.write_provenance(os.path.join(OUT, 'tests.csv'), __file__, {'n_perm': N_PERM, 'window': 5, 'dq': sp.DQ},
                           [os.path.join(OUT, 'phosphosites_annotated.csv'), os.path.join(sp.OUT_DIR, 'phospho_sites_common.csv')])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--reuse-sites', action='store_true', help='run the tests on the stored sites.csv')
    a = ap.parse_args()
    main(a.workers, a.reuse_sites)
