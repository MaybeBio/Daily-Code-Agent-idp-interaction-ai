"""
fig4_design_degree_export.py — per-variant values behind hub_vs_degree/part_e_split.txt (POST-HOC split of Part E,
hub_vs_degree_plan.md 9212f51e), for plotting degree in Fig. 4d. Code copied from redo_evidence/hub_vs_degree_e49.py, with
the variant id and substitution type kept; no new statistics. Checks the medians (+0.348 occupancy, +0.484 degree).
  OMP_NUM_THREADS=4 ~/.venvs/bender/bin/python -I fig4_design_degree_export.py  ->  cache/design_degree_per_variant.csv
"""
import csv, os, pickle, sys
from collections import Counter, defaultdict
import numpy as np
from scipy import stats
HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T
OUT = os.path.join(HERE, 'hub_vs_degree')
deg = pickle.load(open(os.path.join(OUT, 'degree_model_cal.pkl'), 'rb'))
occ = pickle.load(open(os.path.join(HERE, 'hub_predictor/models_v3_2026-10-04.pkl'), 'rb'))[('CALVADOS-2', 'full')]
table = {r['id']: r for r in csv.DictReader(open(os.path.join(HERE, 'gaivi/hub_design_table.csv')))}
kind_of = {r['id']: r['kind'] for r in csv.DictReader(open(os.path.join(HERE, 'design_response', 'design_response_CALVADOS-2.csv')))}
lab = defaultdict(lambda: ([], []))
for r in csv.DictReader(open(os.path.join(HERE, 'labels_design_calvados/per_residue.csv'))):
    lab[r['UniProt_ID']][0].append(float(r['occupancy'])); lab[r['UniProt_ID']][1].append(float(r['deg_mean']))
rows = []
for uid, r in table.items():
    if r['category'] != 'variant' or r['parent'] not in table or uid not in lab or r['parent'] not in lab: continue
    sv, sp = r['sequence'], table[r['parent']]['sequence']
    if len(sv) != len(sp) or sv == sp or len(lab[uid][0]) != len(sv) or len(lab[r['parent']][0]) != len(sp): continue
    if Counter(sv) == Counter(sp): continue                       # shuffles; the 49 substitution variants remain
    Fv, Fp = T.features(sv, 10), T.features(sp, 10)
    vals = []
    for m, j in ((occ, 0), (deg, 1)):
        dp = m.predict(Fv) - m.predict(Fp)
        ds = np.array(lab[uid][j]) - np.array(lab[r['parent']][j])
        vals.append(stats.spearmanr(dp, ds).correlation if np.std(dp) > 0 and np.std(ds) > 0 else np.nan)
    rows.append((uid, kind_of.get(uid, ''), vals[0], vals[1]))
os.makedirs('cache', exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'design_degree_per_variant.csv'), 'w') as fh:
    w = csv.writer(fh); w.writerow(['id', 'kind', 'rho_occupancy', 'rho_degree']); w.writerows(rows)
o = np.array([x[2] for x in rows if x[2] == x[2]]); d = np.array([x[3] for x in rows if x[3] == x[3]])
print(f'n {len(rows)}  occupancy {np.median(o):+.3f}  degree {np.median(d):+.3f}  kinds {Counter(x[1] for x in rows)}')
