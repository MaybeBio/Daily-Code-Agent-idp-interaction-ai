"""
residual_sequence_check.py — do the misses of the calibrated model A0 on the temperature series line up with simple
sequence properties? Exploratory, after the fact, a diagnostic before any change to the hydration forcing; nothing frozen
is refitted.

Series: the four Wuttke et al. 2014 constructs (leave-one-out calibrated slopes, Table S1b) and the seven held-out series
with a measured slope (heldout_temperature/results/primary_<model>/primary.csv, frozen A0). Residual = s_model - s_exp in
% per K (negative: the model contracts too much). Features are computed on the simulated chain (top.pdb for the
held-out runs, input_wuttke/wuttke2014_constructs.fasta for Wuttke):
    N, FCR, |NCPR|, mean hydropathy (Kyte-Doolittle scaled to 0-1, as in the HPS model), SHD (Zheng et al. 2020, with
    that scale), SCD (Sawle & Ghosh 2015), fractions of Gly, Pro, Gly+Pro, aromatic (FWY), aliphatic (AVILM),
    polar (NQST), acidic (DE) and basic (KR).
Statistics: Spearman rho over the 11 series, permutation p (20,000 shuffles), and again with the two flat controls
(sfAFP, Aap) left out (9 series), since both controls are Gly/Pro-rich and are the two largest misses. PNt (two salts)
and the alpha-synuclein loops (one chain) and prothymosin alpha segments are not independent, so the effective number of
proteins is 8 (all series) and 6 (without controls). About 15 features x 2 models are scanned, so a single p < 0.05 is
expected by chance.
    python residual_sequence_check.py      (any env with mdtraj, numpy, pandas, scipy)
Output: heldout_temperature/results/exploratory/residual_sequence.csv, residual_sequence.txt
"""
import sys

sys.dont_write_bytecode = True

import glob  # noqa: E402
import os  # noqa: E402

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

BM = (__import__('os').environ['CS_DATA'] + '/BENDER-MT')
RES = f'{BM}/heldout_temperature/results'
OUT = f'{RES}/exploratory'
WFASTA = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/CALVADOS IDP/input_wuttke/wuttke2014_constructs.fasta'))
HELD = {'PNt-I025': 'HO-PNt-I025-pH74', 'PNt-I175': 'HO-PNt-I175-pH74', 'sfAFP-I175': 'HO-sfAFP-I175-pH75',
        'asyn4W26C-I030': 'HO-asyn4W26C-I030-pH75', 'asyn39W66C-I030': 'HO-asyn39W66C-I030-pH75',
        'AapPGR-I195': 'HO-AapPGR-I195-pH74', 'p53N93-I118': 'HO-p53N93-I118-pH70'}
PROTEIN = {'IN_NTD': 'IN', 'CspM34': 'CspM34', 'ProTa_C2_C56': 'ProTa', 'ProTa_C56_C110': 'ProTa', 'PNt-I025': 'PNt',
           'PNt-I175': 'PNt', 'sfAFP-I175': 'sfAFP', 'asyn4W26C-I030': 'aSyn', 'asyn39W66C-I030': 'aSyn',
           'AapPGR-I195': 'Aap', 'p53N93-I118': 'p53'}
CONTROLS = ['sfAFP-I175', 'AapPGR-I195']
# Table S1b: measured, leave-one-out calibrated C2, MR (% per K)
WUTTKE = {'IN_NTD': (-0.34, -0.30, -0.31), 'CspM34': (-0.29, -0.18, -0.27), 'ProTa_C2_C56': (-0.28, -0.32, -0.28),
          'ProTa_C56_C110': (-0.23, -0.19, -0.24)}
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
      'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
LAM = {a: (v + 4.5) / 9.0 for a, v in KD.items()}
CHG = {'D': -1, 'E': -1, 'K': 1, 'R': 1}
NPERM = 20000


def seq_of_top(name):
    p = sorted(glob.glob(f'{BM}/heldout_CALVADOS2/*K/{name}/top.pdb'))[0]
    return ''.join(r.code for r in md.load_topology(p).residues)


def wuttke_seqs():
    out = {}
    for block in open(WFASTA).read().split('>')[1:]:
        h, *s = block.split('\n')
        out[h.split('|')[0].strip()] = ''.join(s).strip()
    return out


def features(s):
    n = len(s)
    q = np.array([CHG.get(a, 0) for a in s], float)
    lam = np.array([LAM[a] for a in s])
    i, j = np.triu_indices(n, 1)
    sep = (j - i).astype(float)
    f = lambda res: sum(s.count(a) for a in res) / n  # noqa: E731
    return {'N': n, 'FCR': np.abs(q).mean(), 'absNCPR': abs(q.mean()), 'hydropathy': lam.mean(),
            'SHD': ((lam[i] + lam[j]) / sep).sum() / n, 'SCD': (q[i] * q[j] * np.sqrt(sep)).sum() / n,
            'f_G': f('G'), 'f_P': f('P'), 'f_GP': f('GP'), 'f_arom': f('FWY'), 'f_aliph': f('AVILM'),
            'f_polar': f('NQST'), 'f_acid': f('DE'), 'f_basic': f('KR')}


def perm_p(x, y, rng):
    r0 = spearmanr(x, y)[0]
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(NPERM)])
    return r0, float((np.abs(null) >= abs(r0) - 1e-12).mean())


def main():
    rng = np.random.default_rng(20261023)
    ws = wuttke_seqs()
    rows = []
    for k, (se, c2, mr) in WUTTKE.items():
        for model, s in (('CALVADOS-2', c2), ('Mpipi-Recharged', mr)):
            rows.append({'series': k, 'role': 'calibration (LOO)', 'model': model, 's_exp': se, 's_model': s,
                         **features(ws[k])})
    for model, tag in (('CALVADOS-2', 'calvados'), ('Mpipi-Recharged', 'recharged')):
        p = pd.read_csv(f'{RES}/primary_{tag}/primary.csv')
        p = p[(p.variant == 'A0') & (p.model == model)]
        for k, name in HELD.items():
            x = p[p.set == k].iloc[0]
            rows.append({'series': k, 'role': 'held out', 'model': model, 's_exp': x.s_exp, 's_model': x.s_model,
                         **features(seq_of_top(name))})
    df = pd.DataFrame(rows)
    df['protein'] = df.series.map(PROTEIN)
    df['residual'] = df.s_model - df.s_exp
    feats = [c for c in features('ACDE') if c != 'N'] + ['N']
    L = ['Residual of calibrated A0 (s_model - s_exp, % per K) against sequence features; Spearman rho, permutation p',
         f'{NPERM} shuffles; about {len(feats)} features x 2 models scanned, so p < 0.05 once is expected by chance', '']
    stats = []
    for model, d in df.groupby('model'):
        L.append(f'== {model}')
        L.append('  ' + ', '.join(f'{r.series} {r.residual:+.2f}' for r in d.itertuples()))
        for subset, dd in (('all 11', d), ('no controls (9)', d[~d.series.isin(CONTROLS)])):
            for f in feats:
                r, pv = perm_p(dd[f].to_numpy(), dd.residual.to_numpy(), rng)
                g = dd.groupby('protein')[[f, 'residual']].mean()
                rp = spearmanr(g[f], g.residual)[0]
                stats.append({'model': model, 'subset': subset, 'feature': f, 'rho': r, 'perm_p': pv,
                              'rho_by_protein': rp, 'n_series': len(dd), 'n_proteins': len(g)})
        st = pd.DataFrame([s for s in stats if s['model'] == model])
        for subset in ('all 11', 'no controls (9)'):
            t = st[st.subset == subset].sort_values('perm_p')
            L.append(f'  {subset}:')
            for x in t.itertuples():
                L.append(f'    {x.feature:10s} rho {x.rho:+.2f}  p {x.perm_p:.3f}  (by protein {x.rho_by_protein:+.2f}, n={x.n_proteins})')
        L.append('')
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(f'{OUT}/residual_sequence.csv', index=False)
    pd.DataFrame(stats).to_csv(f'{OUT}/residual_sequence_stats.csv', index=False)
    open(f'{OUT}/residual_sequence.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    print(df[['model', 'series', 's_exp', 's_model', 'residual', 'hydropathy', 'SHD', 'FCR', 'f_GP']].round(3).to_string())


if __name__ == '__main__':
    main()
