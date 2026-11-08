"""
phospho_direct_test.py — statistics of phospho_direct_plan.md (hash 5a69b002) with Amendment 1 (hash 807887ae).
Written and hashed with phospho_direct_labels.py before any run was labelled. Refuses to run unless the selection, the
frozen predictions, the labels and rg.csv match their entries in gaivi/FROZEN.sha256 (the labels are hashed after
labelling and before this script is run).

Definitions (fixed by the plan)
  occupancy, degree   per residue and run from phospho_direct/labels/per_residue.csv ('occupancy', 'deg_mean')
  R(X, Y)             sum over residues of |profile(X) - profile(Y)|
  E(s)                mean of the 9 R(P_i, W_j) minus mean of the 3 R(W_i, W_j), P = runs with site s phosphorylated,
                      W = unmodified runs of the same protein
  Sites per protein   A = annotated, U = matched unannotated partner, L = low-fluctuation unannotated site

Decisions (Amendment 1): one-sided Wilcoxon signed-rank over proteins, alpha 0.05, hub occupancy
  T3  E(A) - E(U) > 0     the original primary
  T2  E(U) - E(L) > 0     co-primary, mechanism (physics against reweighting noise)
  T1  E(A) - E(L) > 0     co-primary
  Holm-adjusted P values across the three are reported, and the reading is repeated with them as a sensitivity.
  Reading: annotation-specific if T3 is positive; otherwise physics if T1 and T2 are positive, reweighting noise if
  neither is, mixed in any other case.
Secondary, reported without a decision
  S1  E(s) > 0 per site class (one-sided sign test), hub occupancy and degree
  S2  T1 to T3 with degree in place of hub occupancy
  S3  Spearman of E(s) with log ESS over all sites (negative expected under physics), hub occupancy and degree
  S4  reweighting against direct over all sites: Spearman of the frozen reweighted R_hub with E(s); the degree analogue
      (sum of |frozen reweighted degree change| with E_degree, declared here); per site, Spearman of the frozen
      reweighted change profile with the direct change profile (mean of the phosphorylated minus mean of the
      unmodified replicates), median over sites and per class, with a one-sided Wilcoxon of the per-site values > 0;
      Spearman of the frozen reweighted Rg change with the direct Rg change, and their sign agreement
  S5  unmodified replicates against the stored BENDER run (phospho_network/profiles, the baseline of the reweighting):
      R(W_i, stored) against R(W_i, W_j), and the profile Spearman correlations, per protein, medians reported
A protein with any run missing from the labels is dropped from every statistic and listed.

  ~/.venvs/finches/bin/python phospho_direct_test.py   ->   phospho_direct/summary.txt, phospho_direct/sites_direct.csv
"""
import hashlib
import itertools
import os
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import binomtest, spearmanr, wilcoxon

DYN = os.path.dirname(os.path.abspath(__file__))
GAIVI = os.path.join(DYN, 'gaivi')
SEL = os.path.join(GAIVI, 'phospho_direct', 'selection_amend1.csv')
PRED = os.path.join(GAIVI, 'phospho_direct', 'frozen_predictions_amend1.npz')
OUT = os.path.join(DYN, 'phospho_direct')
LAB = os.path.join(OUT, 'labels', 'per_residue.csv')
RG = os.path.join(OUT, 'rg.csv')
PROF = os.path.join(DYN, 'phospho_network', 'profiles')
FORM = {'A': 'annotated', 'U': 'partner', 'L': 'low'}
ALPHA = 0.05


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def check_frozen(path):
    """The file must match its most recent entry in gaivi/FROZEN.sha256."""
    entry = None
    for line in open(os.path.join(GAIVI, 'FROZEN.sha256')):
        parts = line.split('#')[0].split()
        if len(parts) == 2 and os.path.normpath(os.path.join(GAIVI, parts[1])) == os.path.normpath(path):
            entry = parts[0]
    if entry is None:
        sys.exit(f'not in FROZEN.sha256: {path}')
    if sha256(path) != entry:
        sys.exit(f'hash mismatch against FROZEN.sha256: {path}')


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj, run = np.empty_like(p), 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - k) * p[i]))
        adj[i] = run
    return adj


def R(x, y):
    return float(np.abs(x - y).sum())


def E(P, W):
    return (np.mean([R(p, w) for p in P for w in W]) - np.mean([R(a, b) for a, b in itertools.combinations(W, 2)]))


def reading(pos):
    if pos['T3']:
        return 'annotation-specific'
    if pos['T1'] and pos['T2']:
        return 'physics'
    if not pos['T1'] and not pos['T2']:
        return 'reweighting noise'
    return 'mixed'


def main():
    for f in (SEL, PRED, LAB, RG):
        check_frozen(f)
    sel = pd.read_csv(SEL)
    pred = np.load(PRED)
    lab = pd.read_csv(LAB, usecols=['UniProt_ID', 'group', 'N', 'pos', 'occupancy', 'deg_mean'])
    rg = pd.read_csv(RG).set_index('run')['rg_mean_nm']
    parts = lab.UniProt_ID.str.extract(r'^PHOS-(\w+)-(WT|[AUL])[ST]?\d*-r(\d)$')
    lab['uid'], lab['form'], lab['rep'] = parts[0], parts[1], parts[2].astype(int)
    assert (lab.uid == lab.group).all()
    prof = {r: g.sort_values('pos') for r, g in lab.groupby('UniProt_ID')}

    uids = sorted(sel.UniProt_ID.unique())
    rows, dropped, s5 = [], [], []
    for uid in uids:
        s = sel[sel.UniProt_ID == uid].set_index('role')
        names = {f: [] for f in ('WT', 'A', 'U', 'L')}
        for f in names:
            tag = 'WT' if f == 'WT' else f'{f}{s.loc[FORM[f], "residue"]}{int(s.loc[FORM[f], "pos"])}'
            names[f] = [f'PHOS-{uid}-{tag}-r{k}' for k in (1, 2, 3)]
        allruns = sum(names.values(), [])
        if any(r not in prof or r not in rg.index for r in allruns):
            dropped.append(uid)
            continue
        n = int(s.N.iloc[0])
        get = {r: (prof[r].occupancy.to_numpy(float), prof[r].deg_mean.to_numpy(float)) for r in allruns}
        assert all(len(v[0]) == n for v in get.values()), uid
        Wo = [get[r][0] for r in names['WT']]
        Wd = [get[r][1] for r in names['WT']]
        rgW = np.mean([rg[r] for r in names['WT']])
        for f in ('A', 'U', 'L'):
            role = FORM[f]
            Po = [get[r][0] for r in names[f]]
            Pd = [get[r][1] for r in names[f]]
            dir_o, dir_d = np.mean(Po, 0) - np.mean(Wo, 0), np.mean(Pd, 0) - np.mean(Wd, 0)
            rw_o, rw_d = pred[f'{uid}_{role}_docc'], pred[f'{uid}_{role}_ddeg']
            rows.append({'UniProt_ID': uid, 'N': n, 'role': role, 'pos': int(s.loc[role, 'pos']),
                         'residue': s.loc[role, 'residue'], 'log_ess': float(np.log(s.loc[role, 'ess'])),
                         'E_occ': E(Po, Wo), 'E_deg': E(Pd, Wd),
                         'R_hub_rw': float(s.loc[role, 'R_hub_rw']), 'R_deg_rw_all': float(np.abs(rw_d).sum()),
                         'dRg_rw_pct': float(s.loc[role, 'dRg_rw']),
                         'dRg_direct_pct': (np.mean([rg[r] for r in names[f]]) / rgW - 1) * 100,
                         'prof_rho_occ': spearmanr(rw_o, dir_o).correlation,
                         'prof_rho_deg': spearmanr(rw_d, dir_d).correlation})
        z = np.load(os.path.join(PROF, f'{uid}.npz'))
        so, sd = z['occ0'], z['deg0']
        assert len(so) == n, uid
        pairs = list(itertools.combinations(range(3), 2))
        s5.append({'UniProt_ID': uid,
                   'R_occ_WW': np.mean([R(Wo[a], Wo[b]) for a, b in pairs]),
                   'R_occ_Wstored': np.mean([R(w, so) for w in Wo]),
                   'R_deg_WW': np.mean([R(Wd[a], Wd[b]) for a, b in pairs]),
                   'R_deg_Wstored': np.mean([R(w, sd) for w in Wd]),
                   'rho_occ_WW': np.median([spearmanr(Wo[a], Wo[b]).correlation for a, b in pairs]),
                   'rho_occ_Wstored': np.median([spearmanr(w, so).correlation for w in Wo]),
                   'rho_deg_WW': np.median([spearmanr(Wd[a], Wd[b]).correlation for a, b in pairs]),
                   'rho_deg_Wstored': np.median([spearmanr(w, sd).correlation for w in Wd])})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'sites_direct.csv'), index=False)
    s5 = pd.DataFrame(s5)
    P = {c: df.pivot(index='UniProt_ID', columns='role', values=c) for c in ('E_occ', 'E_deg')}
    nprot = len(P['E_occ'])

    L = [f'phospho_direct test (plan 5a69b002, Amendment 1 807887ae); scipy {scipy.__version__}',
         f'proteins analysed: {nprot}' + (f'; dropped (missing runs): {dropped}' if dropped else '; none dropped'),
         'E(s) = mean R(phospho, unmodified) - mean R(unmodified, unmodified); R = sum over residues of |difference|', '']

    def tests(c):
        t = {'T3': ('annotated', 'partner'), 'T2': ('partner', 'low'), 'T1': ('annotated', 'low')}
        out = {}
        for k, (a, b) in t.items():
            d = (P[c][a] - P[c][b]).to_numpy()
            out[k] = (float(np.median(d)), int((d > 0).sum()), len(d), float(wilcoxon(d, alternative='greater').pvalue))
        adj = holm([out[k][3] for k in t])
        return {k: (*out[k], float(a)) for k, a in zip(t, adj)}

    names = {'T3': 'T3  E(annotated) - E(partner) > 0   [original primary]',
             'T2': 'T2  E(partner) - E(low) > 0         [co-primary, mechanism]',
             'T1': 'T1  E(annotated) - E(low) > 0       [co-primary]'}
    res = tests('E_occ')
    L.append('DECISION  hub occupancy, one-sided Wilcoxon signed-rank over proteins, alpha 0.05')
    for k in ('T3', 'T2', 'T1'):
        med, npos, nn, p, padj = res[k]
        L.append(f'  {names[k]}  median {med:+.4f}, {npos}/{nn} positive, P {p:.4g} (Holm {padj:.4g})'
                 f'  -> {"positive" if p < ALPHA else "not positive"}')
    read = reading({k: res[k][3] < ALPHA for k in res})
    read_h = reading({k: res[k][4] < ALPHA for k in res})
    L.append(f'  READING (fixed in Amendment 1): {read.upper()}')
    L.append(f'  sensitivity, reading with Holm-adjusted P: {read_h}')
    L.append('')
    L.append('Per class (hub occupancy): median E, and E > 0 one-sided sign test   [S1]')
    for role in ('annotated', 'partner', 'low'):
        for c, lab_ in (('E_occ', 'occupancy'), ('E_deg', 'degree')):
            v = P[c][role].to_numpy()
            L.append(f'  {role:9s} {lab_:9s} median E {np.median(v):+.4f}, {int((v > 0).sum())}/{len(v)} > 0, '
                     f'sign test P {binomtest(int((v > 0).sum()), len(v), 0.5, alternative="greater").pvalue:.4g}')
    L.append('')
    L.append('S2  the same comparisons with degree (no decision)')
    rd = tests('E_deg')
    for k in ('T3', 'T2', 'T1'):
        med, npos, nn, p, padj = rd[k]
        L.append(f'  {names[k].split("[")[0].rstrip()}  median {med:+.4f}, {npos}/{nn} positive, P {p:.4g} (Holm {padj:.4g})')
    L.append(f'  pattern with degree: {reading({k: rd[k][3] < ALPHA for k in rd})}')
    L.append('')
    L.append(f'S3  Spearman of E(s) with log ESS over all {len(df)} sites (negative expected under physics)')
    for c in ('E_occ', 'E_deg'):
        r = spearmanr(df.log_ess, df[c])
        L.append(f'  {c}: rho {r.correlation:+.3f}, P {r.pvalue:.3g}')
    L.append('')
    L.append(f'S4  reweighting (frozen) against direct, {len(df)} sites')
    r = spearmanr(df.R_hub_rw, df.E_occ)
    L.append(f'  Spearman(reweighted R_hub, E_occ) {r.correlation:+.3f}, P {r.pvalue:.3g}')
    r = spearmanr(df.R_deg_rw_all, df.E_deg)
    L.append(f'  Spearman(reweighted sum |degree change|, E_deg) {r.correlation:+.3f}, P {r.pvalue:.3g}   [degree analogue]')
    for c, lab_ in (('prof_rho_occ', 'hub occupancy'), ('prof_rho_deg', 'degree')):
        v = df[c].dropna()
        per = ', '.join(f'{role} {df[df.role == role][c].median():+.3f}' for role in ('annotated', 'partner', 'low'))
        L.append(f'  change profile, reweighted vs direct, {lab_}: median Spearman {v.median():+.3f} over {len(v)} sites '
                 f'({per}); Wilcoxon > 0 P {wilcoxon(v, alternative="greater").pvalue:.3g}')
    r = spearmanr(df.dRg_rw_pct, df.dRg_direct_pct)
    agree = int((np.sign(df.dRg_rw_pct) == np.sign(df.dRg_direct_pct)).sum())
    L.append(f'  Rg change: Spearman(reweighted, direct) {r.correlation:+.3f}, P {r.pvalue:.3g}; same sign {agree}/{len(df)}; '
             f'median |dRg| reweighted {df.dRg_rw_pct.abs().median():.2f}%, direct {df.dRg_direct_pct.abs().median():.2f}%')
    L.append('')
    L.append('S5  unmodified replicates against the stored BENDER run (baseline of the reweighting), medians over proteins')
    L.append(f'  hub occupancy: R(W, W) {s5.R_occ_WW.median():.3f}, R(W, stored) {s5.R_occ_Wstored.median():.3f}; '
             f'Spearman W-W {s5.rho_occ_WW.median():.3f}, W-stored {s5.rho_occ_Wstored.median():.3f}')
    L.append(f'  degree:        R(W, W) {s5.R_deg_WW.median():.3f}, R(W, stored) {s5.R_deg_Wstored.median():.3f}; '
             f'Spearman W-W {s5.rho_deg_WW.median():.3f}, W-stored {s5.rho_deg_Wstored.median():.3f}')
    L.append(f'  proteins with R(W, stored) > 1.5 x R(W, W), hub occupancy: '
             f'{list(s5[s5.R_occ_Wstored > 1.5 * s5.R_occ_WW].UniProt_ID)}')
    s5.to_csv(os.path.join(OUT, 'wt_vs_stored.csv'), index=False)
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
