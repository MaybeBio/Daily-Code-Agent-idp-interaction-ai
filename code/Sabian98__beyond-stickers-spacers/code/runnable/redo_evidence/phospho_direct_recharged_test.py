"""
phospho_direct_recharged_test.py — Mpipi-Recharged replication, phospho_direct_plan.md Amendment 2 (hash d9d1640b).
Written and hashed before any Recharged run of this design was labelled.

  --label            label every finished run (hub_labels_gaivi.py, layout openmm, incremental; unfinished runs are listed
                     in labels/failed.csv and tried again next time) -> phospho_direct_recharged/labels/
  --test             statistics on the labels; refuses to run unless the labels file matches its latest FROZEN entry.
                     Without --interim, all 20 proteins (240 runs) must be labelled and this is the pre-registered test.
  --test --interim   the same statistics over the proteins whose 12 runs are all labelled, written to
                     summary_interim_<n>proteins.txt and marked INTERIM; no decision is taken from it.

Definitions as in phospho_direct_test.py: R(X, Y) = sum over residues of |profile(X) - profile(Y)|, E(s) = mean of the 9
R(P_i, W_j) minus mean of the 3 R(W_i, W_j), hub occupancy for the decisions.
Decisions (Amendment 2): one-sided Wilcoxon signed-rank over proteins, alpha 0.05, Holm across the three also reported.
  R2  E(partner) - E(low) > 0   co-primary (replication of T2)
  R1  E(annotated) - E(low) > 0 co-primary (replication of T1)
  R3  E(annotated) - E(partner) > 0
  Reading: replicated if R1 and R2 positive; partly replicated if exactly one; not replicated if neither. R3 positive
  is reported as annotation-specific in Recharged.
Secondary: E > 0 per class (sign tests); degree; Spearman of Recharged E with CALVADOS-2 E over sites (occupancy and
degree); per site, Spearman of the Recharged and CALVADOS-2 direct change profiles (mean phosphorylated minus mean
unmodified), median over sites; Spearman of Recharged E with CALVADOS-2 log ESS; unmodified replicates against the stored
BENDER Recharged labels (hub_frames_v3_recharged): R(W, W) against R(W, stored).

  ~/.venvs/finches/bin/python phospho_direct_recharged_test.py --label
  ~/.venvs/finches/bin/python phospho_direct_recharged_test.py --test [--interim]
"""
import argparse
import itertools
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr, wilcoxon

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
from phospho_direct_test import E, R, check_frozen, holm  # noqa: E402

ROOT = (__import__('os').environ['BSS_DATA'] + '/BENDER-DYN/phospho_direct_recharged/300K')
FASTA = os.path.join(DYN, 'gaivi', 'phospho_direct', 'OpenMPIPI_snapshot_amend2', 'input_phospho', 'phospho_direct_300K.fasta')
SEL = os.path.join(DYN, 'gaivi', 'phospho_direct', 'selection_amend1.csv')
OUT = os.path.join(DYN, 'phospho_direct_recharged')
LAB = os.path.join(OUT, 'labels', 'per_residue.csv')
CAL_SITES = os.path.join(DYN, 'phospho_direct', 'sites_direct.csv')
CAL_LAB = os.path.join(DYN, 'phospho_direct', 'labels', 'per_residue.csv')
STORED = os.path.join(DYN, 'hub_frames_v3_recharged', 'per_residue.csv')
FORM = {'A': 'annotated', 'U': 'partner', 'L': 'low'}
ALPHA = 0.05


def run_names(sel, uid):
    s = sel[sel.UniProt_ID == uid].set_index('role')
    out = {'WT': [f'PHOS-{uid}-WT-r{k}' for k in (1, 2, 3)]}
    for f, role in FORM.items():
        out[f] = [f'PHOS-{uid}-{f}{s.loc[role, "residue"]}{int(s.loc[role, "pos"])}-r{k}' for k in (1, 2, 3)]
    return out


def profiles(path):
    lab = pd.read_csv(path, usecols=['UniProt_ID', 'pos', 'occupancy', 'deg_mean'])
    return {r: (g.sort_values('pos').occupancy.to_numpy(float), g.sort_values('pos').deg_mean.to_numpy(float))
            for r, g in lab.groupby('UniProt_ID')}


def reading(pos):
    if pos['R1'] and pos['R2']:
        return 'REPLICATED'
    if pos['R1'] or pos['R2']:
        return 'PARTLY REPLICATED'
    return 'NOT REPLICATED'


def label():
    os.makedirs(OUT, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(DYN, 'gaivi', 'hub_labels_gaivi.py'), '--layout', 'openmm', '--root', ROOT,
                    '--fasta', FASTA, '--out', os.path.join(OUT, 'labels'), '--no-require-done', '--workers', '8'],
                   check=True)


def test(interim):
    for f in (SEL, CAL_SITES, CAL_LAB, STORED, LAB):
        check_frozen(f)
    sel = pd.read_csv(SEL)
    prof = profiles(LAB)
    cal = profiles(CAL_LAB)
    cal_sites = pd.read_csv(CAL_SITES).set_index(['UniProt_ID', 'role'])
    uids = sorted(sel.UniProt_ID.unique())
    complete = [u for u in uids if all(r in prof for r in sum(run_names(sel, u).values(), []))]
    if not interim and len(complete) < len(uids):
        sys.exit(f'only {len(complete)}/{len(uids)} proteins fully labelled; use --interim or wait')
    if len(complete) == 0:
        sys.exit('no protein fully labelled yet')
    st = pd.read_csv(STORED, usecols=['UniProt_ID', 'pos', 'occupancy', 'deg_mean'])
    st = st[st.UniProt_ID.isin(complete)]
    rows, s5 = [], []
    for uid in complete:
        names = run_names(sel, uid)
        W = [prof[r] for r in names['WT']]
        Wc = [cal[r] for r in names['WT']]
        for f, role in FORM.items():
            P = [prof[r] for r in names[f]]
            Pc = [cal[r] for r in names[f]]
            d_occ = np.mean([p[0] for p in P], 0) - np.mean([w[0] for w in W], 0)
            dc_occ = np.mean([p[0] for p in Pc], 0) - np.mean([w[0] for w in Wc], 0)
            d_deg = np.mean([p[1] for p in P], 0) - np.mean([w[1] for w in W], 0)
            dc_deg = np.mean([p[1] for p in Pc], 0) - np.mean([w[1] for w in Wc], 0)
            c = cal_sites.loc[(uid, role)]
            rows.append({'UniProt_ID': uid, 'role': role, 'E_occ': E([p[0] for p in P], [w[0] for w in W]),
                         'E_deg': E([p[1] for p in P], [w[1] for w in W]), 'E_occ_cal': c.E_occ, 'E_deg_cal': c.E_deg,
                         'log_ess_cal': c.log_ess, 'prof_rho_occ': spearmanr(d_occ, dc_occ).correlation,
                         'prof_rho_deg': spearmanr(d_deg, dc_deg).correlation})
        g = st[st.UniProt_ID == uid].sort_values('pos')
        so, sd_ = g.occupancy.to_numpy(float), g.deg_mean.to_numpy(float)
        pairs = list(itertools.combinations(range(3), 2))
        s5.append({'UniProt_ID': uid, 'R_occ_WW': np.mean([R(W[a][0], W[b][0]) for a, b in pairs]),
                   'R_occ_Wstored': np.mean([R(w[0], so) for w in W]),
                   'R_deg_WW': np.mean([R(W[a][1], W[b][1]) for a, b in pairs]),
                   'R_deg_Wstored': np.mean([R(w[1], sd_) for w in W])})
    df, s5 = pd.DataFrame(rows), pd.DataFrame(s5)
    tag = f'interim_{len(complete)}proteins' if interim else 'final'
    df.to_csv(os.path.join(OUT, f'sites_recharged_{tag}.csv'), index=False)
    piv = {c: df.pivot(index='UniProt_ID', columns='role', values=c) for c in ('E_occ', 'E_deg')}
    L = [f'phospho_direct Mpipi-Recharged replication (Amendment 2, d9d1640b)'
         + (f'  ***INTERIM, {len(complete)} of {len(uids)} proteins, no decision***' if interim else ''),
         f'proteins: {complete}', '']

    def tests(col):
        t = {'R3': ('annotated', 'partner'), 'R2': ('partner', 'low'), 'R1': ('annotated', 'low')}
        out = {}
        for k, (a, b) in t.items():
            d = (piv[col][a] - piv[col][b]).to_numpy()
            out[k] = (float(np.median(d)), int((d > 0).sum()), len(d), float(wilcoxon(d, alternative='greater').pvalue))
        adj = holm([out[k][3] for k in t])
        return {k: (*out[k], float(x)) for k, x in zip(t, adj)}
    lab = {'R2': 'R2  E(partner) - E(low) > 0     [co-primary]', 'R1': 'R1  E(annotated) - E(low) > 0   [co-primary]',
           'R3': 'R3  E(annotated) - E(partner) > 0'}
    for col, name in (('E_occ', 'hub occupancy (decision)'), ('E_deg', 'degree (secondary)')):
        res = tests(col)
        L.append(f'{name}, one-sided Wilcoxon over proteins')
        for k in ('R2', 'R1', 'R3'):
            med, npos, n, p, padj = res[k]
            L.append(f'  {lab[k]}: median {med:+.4f}, {npos}/{n} positive, P {p:.4g} (Holm {padj:.4g})')
        rd = reading({k: res[k][3] < ALPHA for k in res})
        L.append(f'  reading{" (if this were final)" if interim else ""}: {rd}'
                 + ('; R3 positive: annotation-specific in Recharged' if res['R3'][3] < ALPHA else ''))
        L.append('')
    L.append('Per class: median E and E > 0 sign test (occupancy | degree)')
    for role in ('annotated', 'partner', 'low'):
        vo, vd = piv['E_occ'][role].to_numpy(), piv['E_deg'][role].to_numpy()
        L.append(f'  {role:9s} occ median {np.median(vo):+.4f}, {int((vo > 0).sum())}/{len(vo)} > 0 (P '
                 f'{binomtest(int((vo > 0).sum()), len(vo), 0.5, alternative="greater").pvalue:.3g}) | deg median '
                 f'{np.median(vd):+.4f}, {int((vd > 0).sum())}/{len(vd)} > 0')
    L.append('')
    L.append(f'Agreement with CALVADOS-2 over {len(df)} sites')
    for a, b, n in (('E_occ', 'E_occ_cal', 'E occupancy'), ('E_deg', 'E_deg_cal', 'E degree')):
        r = spearmanr(df[a], df[b])
        L.append(f'  Spearman Recharged vs CALVADOS-2 {n}: {r.correlation:+.3f} (P {r.pvalue:.3g})')
    L.append(f'  per-site change profiles, Recharged vs CALVADOS-2: median Spearman {df.prof_rho_occ.median():+.3f} (occupancy), '
             f'{df.prof_rho_deg.median():+.3f} (degree)')
    for c in ('E_occ', 'E_deg'):
        r = spearmanr(df[c], df.log_ess_cal)
        L.append(f'  Spearman Recharged {c} vs CALVADOS-2 log ESS: {r.correlation:+.3f} (P {r.pvalue:.3g}; negative expected)')
    L.append('')
    L.append('Unmodified replicates vs stored BENDER Recharged labels (medians over proteins)')
    L.append(f'  occupancy R(W, W) {s5.R_occ_WW.median():.3f}, R(W, stored) {s5.R_occ_Wstored.median():.3f}; degree R(W, W) '
             f'{s5.R_deg_WW.median():.3f}, R(W, stored) {s5.R_deg_Wstored.median():.3f}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, f'summary_{tag}.txt'), 'w').write(txt)
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label', action='store_true')
    ap.add_argument('--test', action='store_true')
    ap.add_argument('--interim', action='store_true')
    a = ap.parse_args()
    if a.label:
        label()
    if a.test:
        test(a.interim)


if __name__ == '__main__':
    main()
