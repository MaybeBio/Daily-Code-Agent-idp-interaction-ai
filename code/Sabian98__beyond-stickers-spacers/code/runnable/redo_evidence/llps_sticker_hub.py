"""
llps_sticker_hub.py — EXPLORATORY, predictions only: how do stickers inform the predicted hubs of prion-like and RGG
drivers, and why are the wild-type hubs more uneven than shuffles although aromatic stickers are evenly spaced
(Martin et al. 2020)?

For each IDR (IDRome boundaries, frozen consensus model):
  1. mean predicted occupancy relative to the chain mean, per residue class, and by distance to the nearest aromatic;
  2. aromatic spacing: CV of the gaps between consecutive aromatics, wild type against 300 full shuffles
     (low percentile = more even than random, as Martin et al. reported);
  3. clustering of each residue class: CV of its local density (window +-5) along the chain, wild type percentile among
     the shuffles, and the Spearman correlation across shuffles between that clustering and the unevenness (CV) of the
     predicted occupancy.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python llps_sticker_hub.py  ->  llps_csat/llps_sticker_hub.txt
"""
import csv
import os
import pickle
import sys
import zlib
from multiprocessing import Pool

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

IDROME = os.path.join(HERE, '..', 'external_data', 'idrome_clinvar', 'IDRome_DB.csv')
OUT = os.path.join(HERE, 'llps_csat')
IDRS = [('hnRNPA1 LCD', 'P09651_187_372'), ('FUS N-term LC + RGG1', 'P35637_1_285'), ('hnRNPA2 LCD', 'P22626_194_353'),
        ('TDP-43 CTD', 'Q13148_252_414'), ('DDX4 N-term IDR', 'Q9NQI0_1_252'), ('TAF15 N-term LC', 'Q92804_1_234'),
        ('EWSR1 N-term LC', 'Q01844_1_360'), ('TIA1 PrLD', 'P31483_283_386')]
CLASSES = [('aromatic FWY', 'FWY'), ('R', 'R'), ('K', 'K'), ('D/E', 'DE'), ('G', 'G'), ('S/T', 'ST'), ('N/Q', 'NQ'),
           ('P', 'P'), ('hydrophobic AILMV', 'AILMV')]
K = 300
MODEL = None


def init():
    global MODEL
    MODEL = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl'), 'rb'))['consensus']


def occ_cv(s):
    p = np.clip(MODEL.predict(T.features(s, 10)), 0, None)
    return p, float(p.std() / p.mean())


def local_density(s, aas, w=5):
    x = np.array([c in aas for c in s], float)
    return np.array([x[max(0, i - w):i + w + 1].mean() for i in range(len(s))])


def clustering(s, aas):
    d = local_density(s, aas)
    return float(d.std() / d.mean()) if d.mean() > 0 else np.nan


def aro_spacing_cv(s):
    pos = [i for i, c in enumerate(s) if c in 'FWY']
    g = np.diff(pos)
    return float(g.std() / g.mean()) if len(g) > 2 else np.nan


def one(job):
    label, name, seq = job
    p, cv0 = occ_cv(seq)
    rng = np.random.default_rng(zlib.crc32(name.encode()))
    sh = []
    for _ in range(K):
        s = ''.join(rng.permutation(list(seq)))
        sh.append([occ_cv(s)[1], aro_spacing_cv(s)] + [clustering(s, a) for _, a in CLASSES])
    return label, name, seq, p, cv0, np.array(sh)


def main():
    seqs = {r['seq_name']: r['fasta'] for r in csv.DictReader(open(IDROME))}
    with Pool(4, initializer=init) as pl:
        res = pl.map(one, [(lab, n, seqs[n]) for lab, n in IDRS])
    out = []

    def p_(s=''):
        out.append(s)
        print(s, flush=True)
    p_('EXPLORATORY: stickers and predicted hubs in prion-like and RGG drivers (frozen consensus model)\n')
    p_('1. Mean predicted occupancy relative to the chain mean, by residue class (n residues in brackets)')
    p_(f'  {"IDR":22s} ' + ' '.join(f'{c:>11s}' for c, _ in CLASSES))
    for label, name, seq, p, cv0, sh in res:
        rel = p / p.mean()
        cells = []
        for _, aas in CLASSES:
            m = np.array([c in aas for c in seq])
            cells.append(f'{rel[m].mean():5.2f} ({m.sum():3d})' if m.any() else f'{"-":>11s}')
        p_(f'  {label:22s} ' + ' '.join(cells))
    p_('\n   by distance to the nearest aromatic (0 = the aromatic itself):')
    p_(f'  {"IDR":22s} {"0":>6s} {"1":>6s} {"2":>6s} {"3-5":>6s} {">5":>6s}   top-10% hubs that are aromatic / within 2 of one')
    for label, name, seq, p, cv0, sh in res:
        rel = p / p.mean()
        aro = np.array([i for i, c in enumerate(seq) if c in 'FWY'])
        dist = np.array([np.abs(aro - i).min() for i in range(len(seq))])
        bands = [(0, 0), (1, 1), (2, 2), (3, 5), (6, 10 ** 6)]
        cells = [f'{rel[(dist >= a) & (dist <= b)].mean():6.2f}' if ((dist >= a) & (dist <= b)).any() else f'{"-":>6s}'
                 for a, b in bands]
        top = np.argsort(-p)[:int(np.ceil(0.1 * len(seq)))]
        p_(f'  {label:22s} ' + ' '.join(cells) + f'   {np.mean(dist[top] == 0):.0%} / {np.mean(dist[top] <= 2):.0%}')
    p_('\n2. Aromatic spacing (CV of gaps between consecutive aromatics): wild type and its percentile among shuffles')
    for label, name, seq, p, cv0, sh in res:
        a0 = aro_spacing_cv(seq)
        p_(f'  {label:22s} WT {a0:5.2f}  percentile {100 * np.nanmean(sh[:, 1] < a0):4.0f}   (hub unevenness percentile '
           f'{100 * np.mean(sh[:, 0] < cv0):4.0f})')
    p_('\n3. Clustering of each residue class (CV of local density, +-5): wild-type percentile among shuffles [Spearman across')
    p_('   shuffles between that clustering and hub unevenness]')
    p_(f'  {"IDR":22s} ' + ' '.join(f'{c:>14s}' for c, _ in CLASSES))
    for label, name, seq, p, cv0, sh in res:
        cells = []
        for k, (_, aas) in enumerate(CLASSES):
            col = sh[:, 2 + k]
            c0 = clustering(seq, aas)
            if np.isnan(c0) or np.all(np.isnan(col)):
                cells.append(f'{"-":>14s}')
                continue
            ok = ~np.isnan(col)
            r = spearmanr(col[ok], sh[ok, 0]).correlation
            cells.append(f'{100 * np.mean(col[ok] < c0):4.0f} [{r:+.2f}]   ')
        p_(f'  {label:22s} ' + ' '.join(f'{c:>14s}' for c in cells))
    p_('\n   For each IDR, the top-10% predicted hubs with sequence context (hub residue in lower case):')
    for label, name, seq, p, cv0, sh in res[:2]:
        top = sorted(np.argsort(-p)[:int(np.ceil(0.1 * len(seq)))])
        ctx = [seq[max(0, i - 3):i] + seq[i].lower() + seq[i + 1:i + 4] for i in top]
        p_(f'  {label}: ' + ' '.join(f'{i + 1}:{c}' for i, c in zip(top, ctx)))
    open(os.path.join(OUT, 'llps_sticker_hub.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
