"""
case_fus_raster.py — per-snapshot hubs of the FUS experimental ensemble (PED00494) for the hub raster panel.
Same definition as the paper (C-alpha < 8 A, |i-j| >= 3, exact betweenness, top round(0.1 N) non-zero, CRC32-seeded
random ties). Writes case_fus/raster_hubs.csv (model index, Rg, hub residues 1-based) for all complete models.
  ~/miniforge3/envs/openmpipi/bin/python case_fus_raster.py
"""
import csv
import os
import zlib

import networkx as nx
import numpy as np

import case_fus_prep as P

OUT = os.path.join(P.OUT, 'raster_hubs.csv')


def main():
    models = P.read_models(P.PDB)
    rows = []
    for k, m in enumerate(models):
        X = P.ca_coords(m)
        n = len(X)
        d = np.linalg.norm(X[:, None] - X[None], axis=-1)
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from((i, j) for i in range(n) for j in range(i + 3, n) if d[i, j] < 8.0)
        b = nx.betweenness_centrality(G, normalized=True)
        bc = np.array([b[i] for i in range(n)])
        rng = np.random.default_rng(zlib.crc32(f'PED00494|{k}'.encode()))
        tb = rng.random(n)
        cand = np.nonzero(bc > 0)[0]
        hubs = cand[np.lexsort((tb[cand], -np.round(bc[cand], 12)))][:round(0.1 * n)]
        rg = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
        rows.append([k, f'{rg:.2f}', ' '.join(str(h + 1) for h in sorted(hubs))])
    with open(OUT, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['model', 'rg_A', 'hubs'])
        w.writerows(rows)
    occ = np.zeros(270)
    for r in rows:
        for h in r[2].split():
            occ[int(h) - 1] += 1
    occ /= len(rows)
    print(f'{len(rows)} snapshots; occupancy from raster: mean {occ.mean():.3f}, max {occ.max():.3f} at residue {occ.argmax() + 1}')


if __name__ == '__main__':
    main()
