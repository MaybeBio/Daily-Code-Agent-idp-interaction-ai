"""
phospho_function_test.py — Part F of phospho_charge_function_plan.md (hash 28befc80), Amendment 1 (a78bf482, SCD
control) and Amendment 2 (9739188d, site mapping and rank scale). Written and hashed before it was run; refuses to run
unless its inputs match gaivi/FROZEN.sha256.

Within-protein tests of whether Ochoa et al. (2020) functional scores track network responsiveness (-log ESS):
  F1   pooled within-protein Spearman (scaled, centred ranks) of functional score with -log ESS over scored Ser/Thr
  F2   the same after least-squares residuals on the rank of the sequence-predicted -log ESS
  F2b  the same with |dSCD| as the control (Amendment 1)
  F2c  both controls (secondary)
  Reading (Amendment 1): F1, F2, F2b positive -> beyond local charge; F1 positive, F2 or F2b not -> explained by
  sequence charge context as far as shown; F1 not positive -> no link shown.
Secondary: F3 detected (scored) against other Ser/Thr, with and without both controls; F4 F1 with reweighted R_hub;
F1/F2/F2b in proteins with >= 3 scored sites; pooled within-protein Spearman of the score with K/R and D/E counts and
|dSCD|.

  ~/.venvs/finches/bin/python phospho_function_test.py   ->   phospho_function/{summary.txt, merged_sites.csv}
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
from phospho_direct_test import check_frozen  # noqa: E402

SITES = os.path.join(DYN, 'phospho_function', 'sites_human.csv')
SCORES = os.path.join(DYN, 'external_data', 'ochoa2020', 'csv', 'EMS84831-supplement-Table_S3__functional_score.csv')
BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
OUT = os.path.join(DYN, 'phospho_function')
SEED, B, ALPHA, MISMATCH_MAX = 20261013, 20000, 0.05, 0.20


def scaled_centred_rank(v, g):
    """Within-group ranks (average ties) scaled to [0, 1] by (rank - 1) / (n - 1), then centred within group."""
    s = pd.Series(np.asarray(v, float))
    r = s.groupby(g).rank(method='average')
    n = s.groupby(g).transform('size')
    x = (r - 1) / (n - 1)
    return (x - x.groupby(g).transform('mean')).to_numpy()


def resid(y, C):
    if C is None:
        return y
    C = np.column_stack(C)
    return y - C @ np.linalg.lstsq(C, y, rcond=None)[0]


def group_perm(g, rng):
    """One permutation of positions within groups (an index array)."""
    base = np.lexsort((np.arange(len(g)), g))
    order = np.lexsort((rng.random(len(g)), g))
    p = np.empty(len(g), int)
    p[base] = order
    return p


def corr_test(f, x, g, controls, rng):
    """Pooled correlation of residualised centred ranks; f permuted within group. One-sided positive."""
    fr, xr = resid(f, controls), resid(x, controls)
    stat = lambda a: float(a @ xr / np.sqrt((a @ a) * (xr @ xr)))  # noqa: E731
    obs = stat(fr)
    null = np.array([stat(fr[group_perm(g, rng)]) for _ in range(B)])
    return obs, (1 + int((null >= obs).sum())) / (B + 1)


def label_test(lab, x, g, controls, rng):
    """Mean of the (residualised) centred rank of x among labelled sites; labels permuted within group."""
    xr = resid(x, controls)
    obs = float(xr[lab].mean())
    null = np.array([xr[lab[group_perm(g, rng)]].mean() for _ in range(B)])
    return obs, (1 + int((null >= obs).sum())) / (B + 1)


def main():
    for f in (SITES, SCORES):
        check_frozen(f)
    S = pd.read_csv(SITES)
    sc = pd.read_csv(SCORES)
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'taxonomy_id', 'sequence'], low_memory=False)
    hum = bio[bio.taxonomy_id == 9606].set_index('UniProt_ID')['sequence']
    sc = sc[sc.uniprot.isin(hum.index)].groupby(['uniprot', 'position'], as_index=False).functional_score.mean()

    # site mapping (Amendment 2)
    seqmatch = S.groupby('UniProt_ID').seq_matches_bio.all()
    log, keep = [], []
    for uid, g in sc.groupby('uniprot'):
        seq = hum[uid]
        aa = [seq[p - 1] if 1 <= p <= len(seq) else '-' for p in g.position]
        st = sum(a in 'ST' for a in aa)
        y = sum(a == 'Y' for a in aa)
        mis = len(aa) - st - y
        ok = mis <= MISMATCH_MAX * len(aa) and bool(seqmatch.get(uid, False))
        log.append((uid, len(aa), st, y, mis, ok))
        if ok:
            keep.append(g[[a in 'ST' for a in aa]])
    mp = pd.DataFrame(log, columns=['UniProt_ID', 'rows', 'ST', 'Y', 'mismatch', 'kept'])
    scored = pd.concat(keep).rename(columns={'uniprot': 'UniProt_ID', 'position': 'pos'})
    M = S[S.UniProt_ID.isin(mp[mp.kept].UniProt_ID)].merge(scored, on=['UniProt_ID', 'pos'], how='left')
    M['scored'] = M.functional_score.notna()
    M.to_csv(os.path.join(OUT, 'merged_sites.csv'), index=False)

    rng = np.random.default_rng(SEED)
    L = ['phospho_function Part F (plan 28befc80; Amendments a78bf482, 9739188d)', '',
         f'Score rows on human BENDER entries: {int(mp.rows.sum())} in {len(mp)} proteins; on S/T {int(mp.ST.sum())}, '
         f'on Y {int(mp.Y.sum())}, mismatched {int(mp.mismatch.sum())}; proteins dropped {int((~mp.kept).sum())} '
         f'({list(mp[~mp.kept].UniProt_ID)})']

    def analysis(min_scored):
        A = M[M.scored].copy()
        A = A[A.groupby('UniProt_ID').UniProt_ID.transform('size') >= min_scored].reset_index(drop=True)
        g = pd.factorize(A.UniProt_ID)[0]
        r = {k: scaled_centred_rank(v, g) for k, v in (('f', A.functional_score), ('x', -A.log_ess),
                                                      ('seq', -A.seqpred), ('scd', A.abs_dSCD))}
        return A, g, r

    A, g, r = analysis(2)
    L.append(f'Analysis set (>= 2 scored S/T): {A.UniProt_ID.nunique()} proteins, {len(A)} scored S/T; '
             f'median functional score {A.functional_score.median():.3f}')
    L.append('')
    res = {}
    for k, ctrl in (('F1', None), ('F2', [r['seq']]), ('F2b', [r['scd']]), ('F2c', [r['seq'], r['scd']])):
        res[k] = corr_test(r['f'], r['x'], g, ctrl, rng)
    lab = {'F1': 'F1   score vs -log ESS', 'F2': 'F2   beyond the sequence charge model',
           'F2b': 'F2b  beyond |dSCD|', 'F2c': 'F2c  beyond both (secondary)'}
    L.append(f'DECISION tests (pooled within-protein Spearman of scaled centred ranks, {B} permutations, one-sided)')
    for k in ('F1', 'F2', 'F2b', 'F2c'):
        L.append(f'  {lab[k]}: rho {res[k][0]:+.4f}, P {res[k][1]:.4g}'
                 + ('' if k == 'F2c' else f'  -> {"positive" if res[k][1] < ALPHA else "not positive"}'))
    pos = {k: res[k][1] < ALPHA for k in ('F1', 'F2', 'F2b')}
    reading = ('functional phosphosites sit at network-responsive positions BEYOND LOCAL CHARGE' if all(pos.values())
               else 'functional phosphosites sit at network-responsive positions, EXPLAINED BY SEQUENCE CHARGE CONTEXT as far as shown'
               if pos['F1'] else 'NO LINK to function shown')
    L.append(f'  READING (fixed in Amendment 1): {reading}')
    L.append('')

    L.append('Secondary')
    A3, g3, r3 = analysis(3)
    for k, ctrl in (('F1', None), ('F2', [r3['seq']]), ('F2b', [r3['scd']])):
        o, p = corr_test(r3['f'], r3['x'], g3, ctrl, rng)
        L.append(f'  >= 3 scored ({A3.UniProt_ID.nunique()} proteins, {len(A3)} sites) {lab[k][:4].strip()}: rho {o:+.4f}, P {p:.4g}')
    A4 = A[A.R_hub_rw.notna()].reset_index(drop=True)
    A4 = A4[A4.groupby('UniProt_ID').UniProt_ID.transform('size') >= 2].reset_index(drop=True)
    g4 = pd.factorize(A4.UniProt_ID)[0]
    o, p = corr_test(scaled_centred_rank(A4.functional_score, g4), scaled_centred_rank(A4.R_hub_rw, g4), g4, None, rng)
    L.append(f'  F4 score vs reweighted R_hub (N <= 300; {A4.UniProt_ID.nunique()} proteins, {len(A4)} sites): rho {o:+.4f}, P {p:.4g}')
    for name, col in (('K/R within 5', A.KR5), ('D/E within 5', A.DE5), ('|dSCD|', A.abs_dSCD), ('-sequence-predicted log ESS', -A.seqpred)):
        c = scaled_centred_rank(col, g)
        L.append(f'  score vs {name}: pooled within-protein rho {float(r["f"] @ c / np.sqrt((r["f"] @ r["f"]) * (c @ c))):+.4f}')
    # F3: detected vs other S/T
    T = M[M.groupby('UniProt_ID').scored.transform('any') & ~M.groupby('UniProt_ID').scored.transform('all')].reset_index(drop=True)
    gt = pd.factorize(T.UniProt_ID)[0]
    xt = scaled_centred_rank(-T.log_ess, gt)
    ct = [scaled_centred_rank(-T.seqpred, gt), scaled_centred_rank(T.abs_dSCD, gt)]
    lab3 = T.scored.to_numpy()
    for name, ctrl in (('no control', None), ('both charge controls', ct)):
        o, p = label_test(lab3, xt, gt, ctrl, rng)
        L.append(f'  F3 detected vs other S/T ({T.UniProt_ID.nunique()} proteins, {int(lab3.sum())} detected of {len(T)}), '
                 f'{name}: mean centred rank of -log ESS {o:+.4f}, P {p:.4g}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'summary.txt'), 'w').write(txt)
    mp.to_csv(os.path.join(OUT, 'mapping.csv'), index=False)
    print(txt)


if __name__ == '__main__':
    main()
