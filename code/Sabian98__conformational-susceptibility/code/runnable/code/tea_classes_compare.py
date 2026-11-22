"""
tea_classes_compare.py — which residues make the TEA hydration term disagree between force fields?

The TEA term (tea_fdt.py / tea_crossff.py --classes) is split exactly into residue classes, aliphatic (A V L I M C),
Pro, aromatic (F Y W), Arg and the rest (G S T N Q H K D E); the classes add up to the total (CHARMM36m Arg; the
CGenFF-Arg total swaps the Arg class). Per model, on the chains all three simulate:
  A  share of each class in the protein-to-protein spread of the TEA term: Cov(dchi_c, dchi)/Var(dchi)
  B  cross-model Spearman of each class term on its own
  C  cross-model Spearman of the TEA term with one class left out, and of the total, both Arg choices
  D  the same for the total response chi + s dchi with the scales matched as in tea_crossff_compare.py
Output: Analysis_FDT/tea_300K/classes/*.csv and summary.txt.
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import tea_crossff_compare as X  # noqa: E402

TEA_DIR = X.TEA_DIR
OUT = os.path.join(TEA_DIR, 'classes')
CLASSES = ['aliph', 'pro', 'arom', 'argC36', 'rest']
FILES = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}


def load(tea_dir=TEA_DIR):
    base = X.load('charmm36', tea_dir)
    out = {}
    for k, f in FILES.items():
        t = pd.read_csv(os.path.join(tea_dir, f)).set_index('UniProt_ID')
        d = pd.DataFrame({c: t[f'dchi_rg_rel_{c}'] for c in CLASSES + ['argCGenFF']})
        d['tot_c36'] = d[CLASSES].sum(axis=1)
        d['tot_cg'] = d.tot_c36 - d.argC36 + d.argCGenFF
        out[k] = d.join(base[k][['chi']], how='inner')
    return out


def main(tea_dir=TEA_DIR, out=OUT):
    os.makedirs(out, exist_ok=True)
    D = load(tea_dir)
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    T3 = {k: d.loc[trio] for k, d in D.items()}
    rows = []
    for k, d in D.items():
        for set_name, dd in (('own full set', d), ('three-model set', T3[k])):
            for tot in ('tot_c36', 'tot_cg'):
                v = dd[tot].var()
                cls = CLASSES[:-2] + (['argC36'] if tot == 'tot_c36' else ['argCGenFF']) + ['rest']
                rows.append({'model': X.NAMES[k], 'set': set_name, 'arg': tot[4:], 'n': len(dd),
                             **{f'share_{c}': np.cov(dd[c], dd[tot])[0, 1] / v for c in cls},
                             'median_tot': dd[tot].median()})
    A = pd.DataFrame(rows)
    rows_B, rows_C, rows_D = [], [], []
    s = {k: 2.0 * T3['C'].tot_c36.median() / T3[k].tot_c36.median() for k in T3}   # matched scales (CHARMM36m Arg)
    for a, b in X.PAIRS:
        x, y = T3[a], T3[b]
        pair = f'{X.NAMES[a]} vs {X.NAMES[b]}'
        rows_B.append({'pair': pair, 'n': len(trio), 'chi': spearmanr(x.chi, y.chi)[0],
                       **{c: spearmanr(x[c], y[c])[0] for c in CLASSES + ['argCGenFF']}})
        for tot, argc in (('tot_c36', 'argC36'), ('tot_cg', 'argCGenFF')):
            r = {'pair': pair, 'arg': tot[4:], 'total': spearmanr(x[tot], y[tot])[0]}
            for c in CLASSES[:-2] + [argc, 'rest']:
                r[f'without_{c}'] = spearmanr(x[tot] - x[c], y[tot] - y[c])[0]
            r['aliph_only_removed_and_arg'] = spearmanr(x[tot] - x.aliph - x[argc], y[tot] - y.aliph - y[argc])[0]
            rows_C.append(r)
            sa, sb = s[a], s[b]
            r = {'pair': pair, 'arg': tot[4:], 'total_response': spearmanr(x.chi + sa * x[tot], y.chi + sb * y[tot])[0]}
            for c in CLASSES[:-2] + [argc, 'rest']:
                r[f'without_{c}'] = spearmanr(x.chi + sa * (x[tot] - x[c]), y.chi + sb * (y[tot] - y[c]))[0]
            r['without_aliph_and_arg'] = spearmanr(x.chi + sa * (x[tot] - x.aliph - x[argc]),
                                                   y.chi + sb * (y[tot] - y.aliph - y[argc]))[0]
            rows_D.append(r)
    T = {'A_class_shares': A, 'B_rank_by_class': pd.DataFrame(rows_B), 'C_rank_leave_one_out': pd.DataFrame(rows_C),
         'D_total_response_leave_one_out': pd.DataFrame(rows_D)}
    lines = [f'{len(trio)} chains in all three models; classes: aliph AVLIMC, pro P, arom FYW, arg R, rest GSTNQHKDE']
    for name, t in T.items():
        t.to_csv(os.path.join(out, f'{name}.csv'), index=False)
        lines += ['', f'== {name}', t.to_string(index=False, float_format=lambda v: f'{v:.3f}')]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(*sys.argv[1:3])
