"""
hub_vs_degree_e49.py — POST-HOC split of Part E of hub_vs_degree_plan.md (hash 9212f51e). Part E as run used all 75 variants
with the same length as a labelled parent, which includes 26 composition-preserving shuffles. Here the 75 are split into
the substitution variants (composition differs from the parent, the 49 of the paper) and the shuffles.

  OMP_NUM_THREADS=4 ~/.venvs/bender/bin/python -I hub_vs_degree_e49.py  ->  hub_vs_degree/part_e_split.txt
"""
import csv, os, pickle, sys
from collections import Counter, defaultdict
import numpy as np
from scipy import stats
HERE = (__import__('os').environ['BSS_REPO'] + '/redo_evidence')
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T
OUT = os.path.join(HERE, 'hub_vs_degree')
deg = pickle.load(open(os.path.join(OUT, 'degree_model_cal.pkl'), 'rb'))
occ = pickle.load(open(os.path.join(HERE, 'hub_predictor/models_v3_2026-10-04.pkl'), 'rb'))[('CALVADOS-2', 'full')]
table = {r['id']: r for r in csv.DictReader(open(os.path.join(HERE, 'gaivi/hub_design_table.csv')))}
lab = defaultdict(lambda: ([], []))
for r in csv.DictReader(open(os.path.join(HERE, 'labels_design_calvados/per_residue.csv'))):
    lab[r['UniProt_ID']][0].append(float(r['occupancy'])); lab[r['UniProt_ID']][1].append(float(r['deg_mean']))
res = {'substitution': {'occupancy': [], 'degree': []}, 'shuffle': {'occupancy': [], 'degree': []}}
for uid, r in table.items():
    if r['category'] != 'variant' or r['parent'] not in table or uid not in lab or r['parent'] not in lab: continue
    sv, sp = r['sequence'], table[r['parent']]['sequence']
    if len(sv) != len(sp) or sv == sp or len(lab[uid][0]) != len(sv) or len(lab[r['parent']][0]) != len(sp): continue
    kind = 'shuffle' if Counter(sv) == Counter(sp) else 'substitution'
    Fv, Fp = T.features(sv, 10), T.features(sp, 10)
    for k, m, j in (('occupancy', occ, 0), ('degree', deg, 1)):
        dp = m.predict(Fv) - m.predict(Fp)
        ds = np.array(lab[uid][j]) - np.array(lab[r['parent']][j])
        res[kind][k].append(stats.spearmanr(dp, ds).correlation if np.std(dp) > 0 and np.std(ds) > 0 else np.nan)
L = ['POST-HOC split of Part E, CALVADOS-2 design labels (not blind)']
for kind in ('substitution', 'shuffle'):
    L.append(f'{kind} variants: n = {len(res[kind]["occupancy"])}')
    for k in ('occupancy', 'degree'):
        v = np.array([x for x in res[kind][k] if x == x])
        L.append(f'   {k:9s} median per-variant Spearman {np.median(v):+.3f}, above zero in {int((v > 0).sum())} of {len(v)}')
    a = np.array(res[kind]['degree']); b = np.array(res[kind]['occupancy'])
    ok = ~np.isnan(a) & ~np.isnan(b)
    L.append(f'   degree better than occupancy in {int((a[ok] > b[ok]).sum())} of {int(ok.sum())} variants; '
             f'paired Wilcoxon two-sided P = {stats.wilcoxon(a[ok], b[ok]).pvalue:.3g}')
open(os.path.join(OUT, 'part_e_split.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
