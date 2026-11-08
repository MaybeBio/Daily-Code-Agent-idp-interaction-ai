"""
phospho_network_test.py — step 2 of phospho_network_plan.md (hash 5b113f57). Run only after
phospho_network/sites_network.csv is hashed. Joins BENDER-MT's annotations and context columns, then the positive control,
P1, P1b, S1-S5 and the label check.

  OMP_NUM_THREADS=4 ~/.venvs/bender/bin/python -I phospho_network_test.py  ->  phospho_network/summary.txt
"""
import csv
import os
import pickle
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata, spearmanr

DYN = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, os.path.join(DYN, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

OUT = os.path.join(DYN, 'phospho_network')
MTS = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/phosphosite_test/sites.csv'
NPERM = 10000
L = []


def p(s=''):
    L.append(s)
    print(s, flush=True)


mt = {}
for r in csv.DictReader(open(MTS)):
    if r['model'] == 'CALVADOS-2':
        mt[(r['UniProt_ID'], int(r['pos']))] = r
R = []
for r in csv.DictReader(open(os.path.join(OUT, 'sites_network.csv'))):
    m = mt[(r['UniProt_ID'], int(r['pos']))]
    r.update(annotated=m['annotated'] == 'True', evidence=m['evidence'], dRg_mt=float(m['dRg_pct_reweight']),
             local_KR=float(m['local_KR']), local_DE=float(m['local_DE']), pro_plus1=float(m['pro_plus1']),
             term_dist=float(m['term_dist']))
    R.append(r)
f = lambda k: np.array([float(r[k]) for r in R])
prot = np.array([r['UniProt_ID'] for r in R])
ann = np.array([r['annotated'] for r in R])
expv = np.array([r['evidence'] == 'experimental' for r in R])
ess = f('ess_frac_net')
dRg = f('dRg_pct_rw')
p('Phosphorylation rewiring of the contact network at annotated phosphosites (phospho_network_plan.md)')
p(f'proteins {len(set(prot))}, Ser/Thr sites {len(R)}, annotated {int(ann.sum())} (experimental {int((ann & expv).sum())})')
p(f'check against BENDER-MT: dRg recomputed vs stored, Pearson {np.corrcoef(dRg, f("dRg_mt"))[0, 1]:.4f}, '
  f'max |difference| {np.abs(dRg - f("dRg_mt")).max():.3f} percentage points')
p('')


def groups(mask):
    g = defaultdict(list)
    for i in np.where(mask)[0]:
        g[prot[i]].append(i)
    return [np.array(v) for v in g.values() if len(v) >= 2]


def pct_within(x, G):
    out = np.full(len(x), np.nan)
    for v in G:
        out[v] = (rankdata(x[v]) - 1) / (len(v) - 1)
    return out


def perm_mean(val, lab, G, seed=0):
    """Mean of val over labelled sites, against labels shuffled within groups; two-sided P around the null mean."""
    idx = np.concatenate(G)
    grp = np.concatenate([np.full(len(v), j) for j, v in enumerate(G)])
    v, y = val[idx], lab[idx]
    keep = ~np.isnan(v)
    idx, grp, v, y = idx[keep], grp[keep], v[keep], y[keep]
    obs = v[y].mean()
    n_g = np.bincount(grp)
    k_g = np.bincount(grp, weights=y.astype(float)).astype(int)
    starts = np.concatenate([[0], np.cumsum(n_g)[:-1]])
    rng = np.random.RandomState(seed)
    null = np.empty(NPERM)
    for b in range(NPERM):
        o = np.lexsort((rng.random_sample(len(grp)), grp))
        pos = np.empty(len(grp), int)
        pos[o] = np.arange(len(grp)) - starts[grp[o]]
        lab_b = pos < k_g[grp]
        null[b] = v[lab_b].mean()
    m = null.mean()
    pval = (np.sum(np.abs(null - m) >= abs(obs - m)) + 1) / (NPERM + 1)
    return obs, m, pval


def zwithin(x, G):
    out = np.full(len(x), np.nan)
    for v in G:
        s = x[v].std()
        out[v] = (x[v] - x[v].mean()) / s if s > 0 else 0.0
    return out


def p1b(metric, mask):
    G = groups(mask)
    z = zwithin(metric, G)
    cov = [zwithin(np.abs(dRg), G)] + [zwithin(f(k), G) for k in ('local_KR', 'local_DE', 'pro_plus1', 'term_dist')]
    idx = np.concatenate(G)
    X = np.column_stack([c[idx] for c in cov])
    X = np.nan_to_num(X)
    beta = np.linalg.lstsq(np.column_stack([np.ones(len(idx)), X]), z[idx], rcond=None)[0]
    res = np.full(len(z), np.nan)
    res[idx] = z[idx] - np.column_stack([np.ones(len(idx)), X]) @ beta
    return perm_mean(res, ann, G, seed=1)


allm = np.ones(len(R), bool)
G = groups(allm)
obs, null, pc = perm_mean(pct_within(np.abs(dRg), G), ann, G)
p('POSITIVE CONTROL  annotated sites, mean within-protein percentile of |dRg| (BENDER-MT 0.542)')
p(f'  {obs:.3f} (null {null:.3f}), P {pc:.4f}  -> {"reproduced" if 0.53 <= obs <= 0.55 else "NOT reproduced, check the pipeline"}')
p('')
Rhub = f('R_hub')
o1, n1, q1 = perm_mean(pct_within(Rhub, G), ann, G)
o2, n2, q2 = p1b(Rhub, allm)
p('PRIMARY')
p(f'  P1  annotated sites, mean within-protein percentile of R_hub (hub redistribution): {o1:.3f} (null {n1:.3f}), P {q1:.4f}')
p(f'  P1b mean residual of within-protein z(R_hub) after z(|dRg|), local K/R, local D/E, Pro+1, terminus distance: '
  f'{o2:+.4f} (null {n2:+.4f}), P {q2:.4f}')
ok = q1 < 0.05 and o1 > 0.5 and q2 < 0.05 and o2 > 0
p(f'  decision: {"SUCCESS" if ok else "NOT met"}')
p('')
p('S1  mechanism over all sites (residues at least 3 positions from the site)')
s1 = f('s1_krde')
s1 = s1[~np.isnan(s1)]
p(f'  degree change of K/R minus D/E: median {np.median(s1):+.5f}, above zero at {np.mean(s1 > 0) * 100:.1f}% of {len(s1)} sites')
rc = f('s1_rho_charge')
p(f'  Spearman(degree change, local net charge): median {np.nanmedian(rc):+.3f}, positive at {np.nanmean(rc > 0) * 100:.1f}% of sites')
ncpr = {r['UniProt_ID']: float(r['NCPR']) for r in R}
tot = defaultdict(list)
for r in R:
    tot[r['UniProt_ID']].append(float(r['s1_total']))
agree = [(np.median(v) > 0) == (ncpr[u] > 0) for u, v in tot.items() if ncpr[u] != 0]
p(f'  sign of the median total degree change matches the sign of the protein net charge in {sum(agree)} of {len(agree)} proteins')
p('')
p('S2  hubs as leverage points (within-protein Spearman over Ser/Thr, median over proteins)')
d0 = f('deg0_site')
a = [spearmanr(d0[v], Rhub[v]).correlation for v in G if len(v) >= 5]
b = [spearmanr(d0[v], np.abs(dRg[v])).correlation for v in G if len(v) >= 5]
p(f'  baseline degree of the site vs R_hub {np.nanmedian(a):+.3f};  vs |dRg| {np.nanmedian(b):+.3f}  ({len(a)} proteins)')
p('')
p('S3  glutamate proxy: frozen degree model, Ser/Thr -> Glu, against the reweighted degree change')
model = pickle.load(open(os.path.join(DYN, 'hub_vs_degree/degree_model_cal.pkl'), 'rb'))
rho3, chk = [], []
for fn in sorted(os.listdir(os.path.join(OUT, 'profiles'))):
    if fn.startswith('.') or not fn.endswith('.npz'):
        continue
    z = np.load(os.path.join(OUT, 'profiles', fn))
    seq = str(z['seq'])
    w0 = model.predict(T.features(seq, 10))
    for k, pos in enumerate(z['positions']):
        i = int(pos) - 1
        pm = model.predict(T.features(seq[:i] + 'E' + seq[i + 1:], 10))
        rho3.append(spearmanr(pm - w0, z['ddeg'][k]).correlation)
p(f'  median Spearman over {len(rho3)} sites {np.nanmedian(rho3):+.3f}; positive at {np.nanmean(np.array(rho3) > 0) * 100:.1f}%')
pg = np.load(os.path.join(OUT, 'page4_profiles.npz'))
seq = str(pg['seq'])
w0 = model.predict(T.features(seq, 10))
for form in ('HIPK1T51', 'CLK2full'):
    s = list(seq)
    for t in pg[f'{form}_sites']:
        s[int(t) - 1] = 'E'
    pm = model.predict(T.features(''.join(s), 10))
    p(f'  PAGE4 {form}: proxy vs direct degree change {spearmanr(pm - w0, pg[f"{form}_direct_ddeg"]).correlation:+.3f}; '
      f'proxy vs reweighted {spearmanr(pm - w0, pg[f"{form}_rw_ddeg"]).correlation:+.3f}')
p('')
p('S4  PAGE4 validation (from phospho_network.py)')
for line in open(os.path.join(OUT, 'page4_validation.txt')).read().splitlines()[3:]:
    p('  ' + line)
p('')
p('S5  sensitivity')
for lab, mask in (('effective sample fraction >= 0.2', ess >= 0.2), ('experimental annotations only', ~(ann & ~expv))):
    Gm = groups(mask)
    o, nn, q = perm_mean(pct_within(Rhub, Gm), ann & mask, Gm, seed=2)
    oo, nn2, qq = p1b(Rhub, mask)
    p(f'  {lab:36s} P1 {o:.3f} (null {nn:.3f}, P {q:.4f});  P1b {oo:+.4f} (P {qq:.4f})')
for k in ('R_mag', 'R_shape'):
    o, nn, q = perm_mean(pct_within(f(k), G), ann, G, seed=3)
    oo, _, qq = p1b(f(k), allm)
    p(f'  {k:36s} P1 {o:.3f} (null {nn:.3f}, P {q:.4f});  P1b {oo:+.4f} (P {qq:.4f})')
p('')
lab = defaultdict(list)
for r in csv.DictReader(open(os.path.join(DYN, 'hub_predictor/data_v3/residues.csv'))):
    if r['force_field'] == 'CALVADOS-2' and r['protein_id'] in set(prot):
        lab[r['protein_id']].append(float(r['occupancy']))
cc = []
for fn in os.listdir(os.path.join(OUT, 'profiles')):
    if fn.startswith('.') or not fn.endswith('.npz'):
        continue
    z = np.load(os.path.join(OUT, 'profiles', fn))
    u = fn[:-4]
    if len(lab.get(u, [])) == len(z['occ0']):
        cc.append(spearmanr(z['occ0'], lab[u]).correlation)
p(f'CHECK baseline hub occupancy on these frames vs stored labels: median Spearman {np.median(cc):.3f} ({len(cc)} proteins)')
open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
