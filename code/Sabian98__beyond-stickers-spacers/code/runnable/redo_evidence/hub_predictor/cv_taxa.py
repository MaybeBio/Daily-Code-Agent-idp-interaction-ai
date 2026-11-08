"""
cv_taxa.py — leave-one-taxon-out training (2026-10-11, descriptive supplementary analysis; plan in taxa_plan.md).
Same features, model and labels as cv_full.py (data_v3, CALVADOS-2 occupancy), but the six BENDER groups (Bacteria,
Fungi, Mammals, Plants, Protists, Viruses) replace the hash folds: train on five groups, score every protein of the sixth.
The held-out model is also applied to the Mpipi-Recharged labels of the held-out group (transfer).
Writes cv_taxa_per_protein.csv (protein_id, group, force_field, N, condition, spearman, top10_precision) with conditions
logo (all features) and logo_position (position and length features only), and prints medians per group.
  ~/.venvs/bender/bin/python cv_taxa.py [--data data_v3]
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np

import train_hub_model as T
from cv_full import POS_COLS, hgb, scores

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
GROUPS = ['Bacteria', 'Fungi', 'Mammals', 'Plants', 'Protists', 'Viruses']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.join(HERE, 'data_v3'))
    a = ap.parse_args()
    grp = {r['UniProt_ID']: r['group'] for r in csv.DictReader(open(os.path.join(EV, 'hub_frames_v3_cal', 'per_protein.csv')))}
    prots, y_occ = T.load(a.data, 'occupancy')
    cal = [k for k in prots if k[0] == 'CALVADOS-2' and len(y_occ.get(k, {})) == len(prots[k]['sequence'])
           and grp.get(k[1]) in GROUPS]
    feats = {}

    def F(seq):
        if seq not in feats:
            feats[seq] = T.features(seq, 10).astype(np.float32)
        return feats[seq]
    X = np.vstack([F(prots[k]['sequence']) for k in cal])
    Y = np.concatenate([[y_occ[k][i] for i in range(len(prots[k]['sequence']))] for k in cal]).astype(np.float32)
    P = np.concatenate([np.full(len(prots[k]['sequence']), j) for j, k in enumerate(cal)])
    gi = np.array([GROUPS.index(grp[k[1]]) for k in cal])
    TX = gi[P]
    print(f'{len(cal)} CALVADOS-2 proteins, {len(Y)} residues', flush=True)
    rows = []
    rch = [k for k in prots if k[0] == 'Mpipi-Recharged' and grp.get(k[1]) in GROUPS
           and len(y_occ.get(k, {})) == len(prots[k]['sequence'])]
    for g, name in enumerate(GROUPS):
        tr, te = TX != g, TX == g
        m_all = hgb().fit(X[tr], Y[tr])
        m_pos = hgb().fit(X[tr][:, POS_COLS], Y[tr])
        pa, pp = m_all.predict(X[te]), m_pos.predict(X[te][:, POS_COLS])
        Pte, Yte = P[te], Y[te]
        for j in np.unique(Pte):
            s = Pte == j
            k = cal[j]
            for cond, pr in (('logo', pa[s]), ('logo_position', pp[s])):
                rho, prec = scores(pr, Yte[s])
                rows.append([k[1], name, 'CALVADOS-2', int(s.sum()), cond, rho, prec])
        for k in rch:
            if grp[k[1]] != name:
                continue
            seq = prots[k]['sequence']
            rho, prec = scores(m_all.predict(F(seq)), np.array([y_occ[k][i] for i in range(len(seq))]))
            rows.append([k[1], name, 'Mpipi-Recharged', len(seq), 'logo_transfer', rho, prec])
        print(f'{name} held out: done', flush=True)
    with open(os.path.join(HERE, 'cv_taxa_per_protein.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['protein_id', 'group', 'force_field', 'N', 'condition', 'spearman', 'top10_precision'])
        w.writerows(rows)
    by = defaultdict(list)
    for r in rows:
        by[(r[4], r[1])].append(r[5])
    for cond in ('logo', 'logo_position', 'logo_transfer'):
        print(f'{cond:14s}', '  '.join(f'{g[:4]} {np.nanmedian(by[(cond, g)]):.3f}' for g in GROUPS), flush=True)


if __name__ == '__main__':
    main()
