"""
case_fus_prep.py — data for the FUS case-study figure (FUS N-terminal low-complexity domain, PED00494, 270 residues).

Writes to figures/case_fus/:
  case_fus_residues.csv     pos, aa, occupancy in the experimental ensemble (PED00494), CALVADOS-2 and Mpipi-Recharged
                            simulations of the same construct (design batch), the frozen consensus prediction, and the
                            change simulated (both force fields) and predicted for the aromatic-to-serine variant
  ensemble20.pdb            20 conformers of the experimental ensemble at evenly spaced radius-of-gyration quantiles,
                            all atoms, B-factor = experimental hub occupancy of the residue (for colouring)
  highlight.pdb             the conformer closest to the median radius of gyration, B-factor = experimental occupancy
  highlight_contacts.pb     ChimeraX pseudobond file: the contact edges of that conformer (C-alpha < 8 A, |i-j| >= 3)
  highlight_hubs.txt        residues that are hubs of that single snapshot (top 10% betweenness, as in the paper)
  ~/miniforge3/envs/openmpipi/bin/python case_fus_prep.py
"""
import csv
import os
from collections import defaultdict

import networkx as nx
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, '..', '..', 'redo_evidence')
OUT = os.path.join(HERE, 'case_fus')
PDB = (__import__('os').environ['BSS_DATA'] + '/ZOIDBERG-ICML 2027/experimental_data/raw/pre/ped/PED00494/pdbfile.pdb')
ENTRY, CONSTRUCT, VARIANT = 'PED00494', 'PED_FUSNTD', 'VAR_FUSNTD_aro2S'


def read_models(path):
    """List of models, each a list of ATOM lines (MMMx files may lack ENDMDL)."""
    models, cur = [], None
    for line in open(path):
        if line.startswith('MODEL'):
            if cur:
                models.append(cur)
            cur = []
        elif line.startswith(('ATOM', 'HETATM')):
            if cur is None:
                cur = []
            cur.append(line)
        elif line.startswith('ENDMDL') and cur is not None:
            models.append(cur)
            cur = None
    if cur:
        models.append(cur)
    return models


def ca_coords(atoms):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atoms if l[12:16].strip() == 'CA'])


def labels(d, uid):
    o = {}
    for r in csv.DictReader(open(os.path.join(EV, d, 'per_residue.csv'))):
        if r['UniProt_ID'] == uid:
            o[int(r['pos'])] = float(r['occupancy'])
    return np.array([o[i] for i in range(len(o))])


def main():
    os.makedirs(OUT, exist_ok=True)
    exp, aa = {}, {}
    for r in csv.DictReader(open(os.path.join(EV, 'hub_predictor', 'experimental_occupancy.csv'))):
        if r['entry'] == ENTRY and int(r['pos']) not in exp:
            exp[int(r['pos'])] = float(r['occ_exp'])
            aa[int(r['pos'])] = r['aa']
    n = len(exp)
    exp = np.array([exp[i] for i in range(n)])
    seq = ''.join(aa[i] for i in range(n))
    pred = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, 'gaivi', 'predictions_frozen_2026-10-02.csv'))):
        if r['id'] in (CONSTRUCT, VARIANT):
            pred[r['id']].append(float(r['consensus']))
    cal, rch = labels('labels_design_calvados', CONSTRUCT), labels('labels_design_recharged', CONSTRUCT)
    cal_v, rch_v = labels('labels_design_calvados', VARIANT), labels('labels_design_recharged', VARIANT)
    assert len(cal) == len(rch) == len(pred[CONSTRUCT]) == n, (len(cal), len(rch), len(pred[CONSTRUCT]), n)
    p0, pv = np.array(pred[CONSTRUCT]), np.array(pred[VARIANT])
    with open(os.path.join(OUT, 'case_fus_residues.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['pos', 'aa', 'occ_experiment', 'occ_calvados2', 'occ_recharged', 'pred_consensus',
                    'aro2S_dsim_calvados2', 'aro2S_dsim_recharged', 'aro2S_dpred'])
        for i in range(n):
            w.writerow([i + 1, seq[i], f'{exp[i]:.4f}', f'{cal[i]:.4f}', f'{rch[i]:.4f}', f'{p0[i]:.4f}',
                        f'{cal_v[i] - cal[i]:.4f}', f'{rch_v[i] - rch[i]:.4f}', f'{pv[i] - p0[i]:.4f}'])

    models = read_models(PDB)
    rg = []
    keep = []
    for k, m in enumerate(models):
        X = ca_coords(m)
        if len(X) != n:
            continue
        keep.append(k)
        rg.append(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
    rg = np.array(rg)
    order = np.argsort(rg)
    pick = [keep[order[int(q * (len(order) - 1))]] for q in np.linspace(0.05, 0.95, 20)]
    hi = keep[order[len(order) // 2]]

    def write(path, idxs):
        with open(path, 'w') as fh:
            for j, k in enumerate(idxs, 1):
                fh.write(f'MODEL     {j:4d}\n')
                resmap, r = {}, -1
                for l in models[k]:
                    key = (l[21], l[22:27])
                    if key not in resmap:
                        r += 1
                        resmap[key] = r
                    b = exp[resmap[key]] if resmap[key] < n else 0.0
                    fh.write(l[:60] + f'{b:6.2f}' + l[66:].rstrip('\n') + '\n')
                fh.write('ENDMDL\n')
    write(os.path.join(OUT, 'ensemble20.pdb'), pick)
    write(os.path.join(OUT, 'highlight.pdb'), [hi])

    X = ca_coords(models[hi])
    G = nx.Graph()
    G.add_nodes_from(range(n))
    d = np.linalg.norm(X[:, None] - X[None], axis=-1)
    edges = [(i, j) for i in range(n) for j in range(i + 3, n) if d[i, j] < 8.0]
    G.add_edges_from(edges)
    bc = nx.betweenness_centrality(G, normalized=True)
    b = np.array([bc[i] for i in range(n)])
    hubs = sorted(np.argsort(-b)[:round(0.1 * n)])
    resnums = []
    seen = set()
    for l in models[hi]:
        key = (l[21], l[22:27])
        if l[12:16].strip() == 'CA' and key not in seen:
            seen.add(key)
            resnums.append((l[21], l[22:26].strip()))
    with open(os.path.join(OUT, 'highlight_contacts.pb'), 'w') as fh:
        fh.write('; radius = 0.15\n; color = gray\n; dashes = 0\n')
        for i, j in edges:
            fh.write(f'/{resnums[i][0]}:{resnums[i][1]}@CA /{resnums[j][0]}:{resnums[j][1]}@CA\n')
    open(os.path.join(OUT, 'highlight_hubs.txt'), 'w').write(
        ' '.join(f'{resnums[i][1]}' for i in hubs) + '\n')
    print(f'{len(models)} models ({len(keep)} complete); Rg {rg.min():.1f}-{rg.max():.1f} A; highlight model {hi} '
          f'(Rg {rg[order[len(order) // 2]]:.1f} A), {len(edges)} contacts, {len(hubs)} snapshot hubs; chain {resnums[0][0]}, '
          f'residues {resnums[0][1]}-{resnums[-1][1]}')


if __name__ == '__main__':
    main()
