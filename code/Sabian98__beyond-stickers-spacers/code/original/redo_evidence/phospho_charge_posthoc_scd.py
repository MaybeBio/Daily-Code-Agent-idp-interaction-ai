"""
phospho_charge_posthoc_scd.py — POST HOC (Amendment 1, a78bf482, of phospho_charge_function_plan.md): Part C repeated
with |dSCD| (change in sequence charge decoration when the site goes from 0 to -2) in place of the sequence charge
model. Post hoc because Part C (phospho_charge_test.py, 60a25f47) had already been run.

  ~/.venvs/finches/bin/python phospho_charge_posthoc_scd.py   ->   phospho_charge/posthoc_scd.txt
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
import phospho_charge_test as C  # noqa: E402
from phospho_direct_test import check_frozen  # noqa: E402
from phospho_function_sites import dscd  # noqa: E402

DS = os.path.join(DYN, 'phospho_charge', 'direct_sites.csv')


def main():
    check_frozen(DS)
    d = pd.read_csv(DS)
    d['abs_dSCD'] = [abs(dscd(C.seq_of(u), p)) for u, p in zip(d.UniProt_ID, d.pos)]
    g = d.UniProt_ID.to_numpy()
    c = {k: C.centre(d[k], g) for k in ('log_ess', 'abs_dSCD', 'seqpred', 'y_occ', 'y_deg')}
    rng = np.random.default_rng(C.SEED)
    L = ['POST HOC: Part C with |dSCD| as the sequence control (Amendment 1)', '',
         'Within-protein Spearman with direct rewiring (y = E / R(W, W), centred):',
         f'  |dSCD| {spearmanr(c["abs_dSCD"], c["y_occ"]).correlation:+.3f} (occ), '
         f'{spearmanr(c["abs_dSCD"], c["y_deg"]).correlation:+.3f} (deg);   -log ESS '
         f'{spearmanr(-c["log_ess"], c["y_occ"]).correlation:+.3f} (occ), {spearmanr(-c["log_ess"], c["y_deg"]).correlation:+.3f} (deg)',
         f'  within-protein Spearman of |dSCD| with -log ESS {spearmanr(c["abs_dSCD"], -c["log_ess"]).correlation:+.3f}', '',
         f'Permutation tests ({C.B} permutations, one-sided negative for log ESS, positive for |dSCD|):']
    for name, y in (('E_occ', c['y_occ']), ('E_deg', c['y_deg'])):
        b, p = C.perm_coef(y, c['log_ess'], c['abs_dSCD'], g, rng)
        L.append(f'  {name}: log ESS beyond |dSCD|: b {b:+.4f}, P {p:.4g}')
        b, p = C.perm_coef(-y, c['abs_dSCD'], c['log_ess'], g, rng)       # -y so that "more rewiring" is negative b
        L.append(f'  {name}: |dSCD| beyond log ESS: b {-b:+.4f}, P {p:.4g}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(DYN, 'phospho_charge', 'posthoc_scd.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
