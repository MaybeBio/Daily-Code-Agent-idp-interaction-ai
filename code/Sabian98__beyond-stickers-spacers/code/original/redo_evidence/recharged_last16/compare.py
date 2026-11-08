"""
compare.py — do the main Mpipi-Recharged numbers of the RECOMB paper change when the last 16 runs are added? (requested 2026-11-04,
descriptive). The paper's labels (hub_frames_v3_recharged, 2026-10-09) cover 10,266 of 10,282 sequences, since the trajectories of the
other 16 arrived on 2026-10-20. Labels for those 16 come from the same code (gaivi/hub_labels_gaivi.py, same hub_frames hash) in labels/.
Each number is computed on the 10,266, to check that it reproduces the paper, and on all 10,282. Trained models are not refit.
  ~/.venvs/finches/bin/python compare.py   ->   summary.txt
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
COLS = ['UniProt_ID', 'pos', 'aa', 'z_posadj', 'occupancy', 'deg_mean']


def partial(x, y):
    n = len(x)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1]) if ex.std() > 0 and ey.std() > 0 else np.nan


def per_protein(cal, rec):
    rows = []
    cal = {k: g.sort_values('pos') for k, g in cal.groupby('UniProt_ID')}
    for u, b in rec.groupby('UniProt_ID'):
        b = b.sort_values('pos')
        row = {'id': u, 'occ_deg_part': partial(b.occupancy.values, b.deg_mean.values)}
        a = cal.get(u)
        if a is not None and len(a) == len(b):
            for name, col in (('occ', 'occupancy'), ('deg', 'deg_mean')):
                x, y = a[col].values, b[col].values
                if x.std() > 0 and y.std() > 0:
                    row[f'{name}_raw'] = float(spearmanr(x, y).correlation)
                    row[f'{name}_part'] = partial(x, y)
        rows.append(row)
    return pd.DataFrame(rows).set_index('id')


def enrichment(rec, letters):
    """Share of a residue set among hubs (top 10% of occupancy per protein) relative to its share among all residues, pooled."""
    r = rec.copy()
    r['rank'] = r.groupby('UniProt_ID').occupancy.rank(ascending=False, method='first')
    r['n'] = r.groupby('UniProt_ID').pos.transform('size')
    hub = r['rank'] <= np.maximum(1, np.round(0.1 * r.n))
    inset = r.aa.isin(list(letters))
    return float(inset[hub].mean() / inset.mean())


def main():
    pp_old = pd.read_csv(os.path.join(EV, 'hub_frames_v3_recharged/per_protein.csv'))
    pp_new = pd.read_csv(os.path.join(HERE, 'labels/per_protein.csv'))
    old = pd.read_csv(os.path.join(EV, 'hub_frames_v3_recharged/per_residue.csv'), usecols=COLS)
    new = pd.read_csv(os.path.join(HERE, 'labels/per_residue.csv'), usecols=COLS)
    cal = pd.read_csv(os.path.join(EV, 'hub_frames_v3_cal/per_residue.csv'), usecols=COLS)
    assert not set(pp_old.UniProt_ID) & set(pp_new.UniProt_ID)
    allr = pd.concat([old, new])
    allp = pd.concat([pp_old, pp_new])
    t_old, t_new = per_protein(cal, old), per_protein(cal[cal.UniProt_ID.isin(new.UniProt_ID)], new)
    t_all = pd.concat([t_old, t_new])
    sb = lambda s: float((2 * s / (1 + s)).median())  # noqa: E731
    rows = [('Hub set overlap, successive snapshots (Jaccard), median', pp_old.jaccard_successive.median(), allp.jaccard_successive.median(), '0.061 to 0.064 over force fields'),
            ('Split-half reliability of occupancy, Spearman-Brown, median', sb(pp_old.splithalf_r_occupancy), sb(allp.splithalf_r_occupancy), '0.95'),
            ('CALVADOS-2 vs Recharged, occupancy, raw, median', t_old.occ_raw.median(), t_all.occ_raw.median(), '0.69 in the older record'),
            ('CALVADOS-2 vs Recharged, occupancy, position removed, median', t_old.occ_part.median(), t_all.occ_part.median(), '0.62'),
            ('CALVADOS-2 vs Recharged, degree, raw, median', t_old.deg_raw.median(), t_all.deg_raw.median(), ''),
            ('CALVADOS-2 vs Recharged, degree, position removed, median', t_old.deg_part.median(), t_all.deg_part.median(), '0.69'),
            ('Occupancy vs degree within protein, position removed, median', t_old.occ_deg_part.median(), t_all.occ_deg_part.median(), '0.94 to 0.96 (another position adjustment)'),
            ('Enrichment among hubs, F W Y R', enrichment(old, 'FWYR'), enrichment(allr, 'FWYR'), '2.22'),
            ('Enrichment among hubs, F W Y', enrichment(old, 'FWY'), enrichment(allr, 'FWY'), '5.94'),
            ('Enrichment among hubs, G S', enrichment(old, 'GS'), enrichment(allr, 'GS'), '1.37')]
    z_old, z_all = old.groupby('aa').z_posadj.mean(), allr.groupby('aa').z_posadj.mean()
    L = ['Mpipi-Recharged: paper labels (10,266) against all runs (10,282), descriptive, models not refit',
         f'sequences {pp_old.shape[0]} -> {allp.shape[0]}; residues {len(old):,} -> {len(allr):,}; share added {100 * len(new) / len(allr):.3f}% of residues', '']
    w = max(len(r[0]) for r in rows)
    L.append(f'{"":{w}s}  10,266    10,282    change   paper')
    for name, a, b, paper in rows:
        L.append(f'{name:{w}s}  {a:7.3f}   {b:7.3f}   {b - a:+.4f}  {paper}')
    L += ['', f'Mean position-adjusted z per residue type, 20 types: largest change {(z_all - z_old).abs().max():.4f} ({(z_all - z_old).abs().idxmax()}), '
          f'Spearman between the two rankings {spearmanr(z_old, z_all).correlation:.4f}',
          f'Largest change in any median above: {max(abs(r[2] - r[1]) for r in rows):.4f}', '',
          'The 16 sequences themselves: ' + f'median occupancy agreement raw {t_new.occ_raw.median():.3f}, position removed {t_new.occ_part.median():.3f}, length {pp_new.N.min()} to {pp_new.N.max()} residues']
    txt = '\n'.join(L) + '\n'
    open(os.path.join(HERE, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
