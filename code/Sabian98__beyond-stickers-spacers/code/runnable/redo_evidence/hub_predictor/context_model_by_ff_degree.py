"""
context_model_by_ff_degree.py — degree_models_plan.md (1a1deeb9), analysis 3: copy of context_model_by_ff.py with
target deg_mean, shuffle degree labels, the experimental degree of the 18 development PED ensembles for the neighbourhood
(degree_wetlab/per_residue.csv), full-model numbers from cv_full_per_protein_degree.csv, outputs *_degree. Original:
context_model_by_ff.py — does a sticker-and-spacer description of each residue's surroundings predict hub occupancy in
every force field, is one such description shared between force fields, and what does it learn?

Extends part C of sticker_spacer_test.py (CALVADOS-2 only) to the three label sets of data_v3, with the same folds
(protein-ID hash), booster settings and feature definitions:
  position        position + length (cv_full.py POS_COLS)
  identity        position + one-hot of the residue itself
  context own     position + sticker_spacer features with the force field's own stickiness scale
  context other   the same with the other model family's scale (CALVADOS-2 lambda <-> Mpipi-Recharged mean epsilon)
  full            the 487 sequence features (CALVADOS-2 from cv_full_per_protein.csv, the two Mpipi models trained here)
sticker_spacer features (sticker_spacer_test.sticker_features): the residue's own stickiness and charge; mean stickiness,
aromatic (F/W/Y) fraction, G/S fraction and net charge in windows of +-2, +-5 and +-10; whole-sequence mean stickiness,
aromatic fraction and net charge per residue. Stickiness scales: CALVADOS-2 lambda, and the mean Wang-Frenkel epsilon
over the 20 partner types of Mpipi-GG (run.in files) and Mpipi-Recharged (OpenMpipi recharged_params.txt). Charges
D/E -1, K/R +1, others 0, as before.
Also
  transfer        the context models trained on the CALVADOS-2 folds applied to the same held-out proteins in Mpipi-GG
                  and Mpipi-Recharged (full-model transfer from cv_full_per_protein.csv)
  shuffles        position, identity and context-own models trained on all natural CALVADOS-2 proteins, applied to the
                  1,791 simulated shuffles (hub_frames_v3_shuffle). Post hoc; the frozen full model reached 0.889.
  importance      grouped permutation importance of the context-own model on the held-out folds (drop of the median
                  per-protein Spearman when a feature group is permuted across the residues of the test fold)
  neighbourhood   composition of the +-5 window around each residue (the residue itself excluded) for hubs beyond
                  position (top 10% of position-adjusted occupancy per protein) relative to all residues, for the three
                  force fields and the 27 PED ensembles, with 95% bootstrap intervals over proteins
  ~/.venvs/bender/bin/python context_model_by_ff.py
      ->  context_model_by_ff_summary.txt, context_model_by_ff_per_protein.csv, context_neighbourhood.csv
"""
import csv
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

import sticker_spacer_by_ff as SB
import sticker_spacer_test as SS
import train_hub_model as T

HERE = SS.HERE
EV = os.path.dirname(HERE)
sys.path.insert(0, EV)
import arrangement_test as AT  # noqa: E402

AAS = SS.AAS
POS_COLS = SS.POS_COLS
FFS = ('CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged')
NICE = {'CALVADOS-2': 'CALVADOS-2', 'MPIPI-GG': 'Mpipi-GG', 'Mpipi-Recharged': 'Mpipi-Recharged', 'PED': 'PED'}
OWN = {'CALVADOS-2': 'lambda', 'MPIPI-GG': 'eps GG', 'Mpipi-Recharged': 'eps RC'}
OTHER = {'CALVADOS-2': 'eps RC', 'MPIPI-GG': 'lambda', 'Mpipi-Recharged': 'lambda'}
# column groups of [position | sticker_spacer features] for the permutation importance
GROUPS = {'position and length': [0, 1, 2, 3], 'own stickiness': [4], 'own charge': [5],
          'window stickiness': [6, 10, 14], 'window aromatic fraction': [7, 11, 15], 'window G/S fraction': [8, 12, 16],
          'window net charge': [9, 13, 17], 'whole-sequence features': [18, 19, 20]}
CLASSES = [('aromatic F W Y', 'FWY'), ('R', 'R'), ('K', 'K'), ('acidic D E', 'DE'), ('G', 'G'), ('S', 'S'),
           ('polar N Q T', 'NQT'), ('P', 'P'), ('aliphatic A V I L M', 'AVILM'), ('H', 'H'), ('C', 'C')]
NBOOT = 1000
SMOKE = os.environ.get('SMOKE') == '1'        # 300 proteins per force field, 60 shuffles, outputs *_smoke
SUF = '_degree' + ('_smoke' if SMOKE else '')


def read_deg(dirname):
    """arrangement_test.read with deg_mean in place of occupancy: {id: (sequence, degree array)}."""
    path = os.path.join(EV, dirname, 'per_residue.csv')
    done = {r['UniProt_ID']: int(r['N']) for r in csv.DictReader(open(os.path.join(EV, dirname, 'per_protein.csv')))}
    d, aa = defaultdict(dict), defaultdict(dict)
    with open(path) as fh:
        rd = csv.reader(fh)
        h = next(rd)
        iu, ip, ia, io = h.index('UniProt_ID'), h.index('pos'), h.index('aa'), h.index('deg_mean')
        for row in rd:
            if len(row) <= io or row[iu] not in done:
                continue
            d[row[iu]][int(row[ip])] = float(row[io])
            aa[row[iu]][int(row[ip])] = row[ia]
    res = {}
    for u, v in d.items():
        if len(v) == done[u] and max(v) == len(v) - 1:
            res[u] = (''.join(aa[u][i] for i in range(len(v))), np.array([v[i] for i in range(len(v))]))
    return res


def ped_degree_sets():
    """(sequences, experimental degree arrays) of the 18 development PED ensembles (degree_wetlab/per_residue.csv)."""
    by = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, 'degree_wetlab', 'per_residue.csv'))):
        if r['set'] == 'used18':
            by[r['entry']].append((int(r['pos']), r['aa'], float(r['deg_exp'])))
    seqs, Y = [], []
    for e in sorted(by):
        rows = sorted(by[e])
        seqs.append(''.join(a for _, a, _ in rows))
        Y.append(np.array([v for _, _, v in rows]))
    return seqs, Y
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def scales():
    tab = list(csv.DictReader(open(SS.CAL_TABLE)))
    lam = {r['one']: float(r['lambdas']) for r in tab}
    chg = {a: 0.0 for a in AAS}
    chg.update({'D': -1.0, 'E': -1.0, 'K': 1.0, 'R': 1.0})
    gg = json.load(open(os.path.join(EV, 'mpipi_epsilon_table.json')))
    rc = json.load(open(os.path.join(EV, 'mpipi_recharged_parameters.json')))['epsilon_kJ_per_mol']
    sc = {'lambda': {a: lam[a] for a in AAS},
          'eps GG': {a: float(np.mean([gg[a + b] for b in AAS if a + b in gg])) for a in AAS},
          'eps RC': {a: float(np.mean([rc[a + b] for b in AAS])) for a in AAS}}
    for d in list(sc.values()) + [chg]:                # unusual letters get the mean, so features stay defined
        m = float(np.mean([d[a] for a in AAS]))
        for x in 'BXZUO':
            d.setdefault(x, m if d is not chg else 0.0)
    return sc, chg


def onehot(s):
    return np.eye(20, dtype=np.float32)[[AAS.index(c) if c in AAS else 0 for c in s]]


def per_protein(pred, Y, bounds):
    """Spearman per protein; the residues of protein j are rows bounds[j]:bounds[j + 1]."""
    r = np.full(len(bounds) - 1, np.nan)
    for j in range(len(bounds) - 1):
        a, b = bounds[j], bounds[j + 1]
        if Y[a:b].std() > 0 and pred[a:b].std() > 0:
            r[j] = spearmanr(pred[a:b], Y[a:b]).correlation
    return r


def label_set(prots, y, ff):
    keys = [k for k in prots if k[0] == ff and len(y.get(k, {})) == len(prots[k]['sequence'])]
    keys = keys[:300] if SMOKE else keys
    seqs = [prots[k]['sequence'] for k in keys]
    Y = [np.array([y[k][i] for i in range(len(s))], np.float32) for k, s in zip(keys, seqs)]
    return keys, seqs, Y


def main():
    t0 = time.time()
    sc, chg = scales()
    for f in (f'context_neighbourhood{SUF}.csv', f'context_model_by_ff_per_protein{SUF}.csv', f'context_model_by_ff_summary{SUF}.txt'):
        assert not os.path.exists(os.path.join(HERE, f)), f'{f} exists; refusing to overwrite'
    prots, y = T.load(os.path.join(HERE, 'data_v3'), 'deg_mean')
    feats = {}

    def full(seq):
        if seq not in feats:
            feats[seq] = T.features(seq, 10).astype(np.float32)
        return feats[seq]

    def ctx(seq, scale):
        return SS.sticker_features(seq, sc[scale], chg).astype(np.float32)

    cv_full = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, f'cv_full_per_protein_degree{"_smoke" if SMOKE else ""}.csv'))):
        if r['spearman'] not in ('', 'nan'):
            cv_full[(r['force_field'], r['condition'])][r['protein_id']] = float(r['spearman'])

    rows, med, fold_models, store = [], {}, {}, {}
    p('Sticker-and-spacer context models by force field, TARGET CONTACT DEGREE (data_v3, 5-fold CV by protein, median per-protein Spearman)')
    for ff in FFS:
        keys, seqs, Y = label_set(prots, y, ff)
        keep_full = ff != 'CALVADOS-2'
        P = np.vstack([full(s)[:, POS_COLS] for s in seqs])
        X = {'position': P, 'identity': np.hstack([P, np.vstack([onehot(s) for s in seqs])]),
             'context own': np.hstack([P, np.vstack([ctx(s, OWN[ff]) for s in seqs])]),
             'context other': np.hstack([P, np.vstack([ctx(s, OTHER[ff]) for s in seqs])])}
        if keep_full:
            X['full'] = np.vstack([full(s) for s in seqs])
        else:
            feats.clear()                                # free the CALVADOS-2 full matrices
        Yc = np.concatenate(Y)
        G = np.concatenate([[0], np.cumsum([len(s) for s in seqs])])
        FO = np.concatenate([np.full(len(s), int(prots[k]['fold'])) for k, s in zip(keys, seqs)])
        p(f'\n{NICE[ff]}: {len(keys)} proteins, {len(Yc)} residues, own scale {OWN[ff]}, other scale {OTHER[ff]} '
          f'({time.time() - t0:.0f} s)')
        pred = {c: np.zeros(len(Yc), np.float32) for c in X}
        for f in range(5):
            tr, te = FO != f, FO == f
            for c, M in X.items():
                m = SS.hgb().fit(M[tr], Yc[tr])
                pred[c][te] = m.predict(M[te])
                if c in ('context own', 'context other'):
                    fold_models[(ff, c, f)] = m
            print(f'  fold {f} done ({time.time() - t0:.0f} s)', flush=True)
        res = {c: per_protein(pred[c], Yc, G) for c in X}
        if not keep_full:
            res['full'] = np.array([cv_full[(ff, 'hgb')].get(k[1], np.nan) for k in keys])
        for c, v in res.items():
            med[(ff, c)] = np.nanmedian(v)
            for k, val in zip(keys, v):
                rows.append({'protein_id': k[1], 'force_field': NICE[ff], 'condition': c, 'spearman': val})
        p('  ' + '   '.join(f'{c} {med[(ff, c)]:.3f}' for c in ('position', 'identity', 'context own', 'context other',
                                                               'full')))
        store[ff] = (X['context own'], Yc, G, FO)
        if keep_full:
            del X['full']
        else:
            Xcal = (X, Yc)                               # reused for the shuffles

    # transfer of the CALVADOS-2 context models to the other force fields
    p('\nTransfer of models trained on the CALVADOS-2 folds to the same held-out proteins (median Spearman)')
    sims = {}
    for ff in FFS[1:]:
        keys, seqs, Y = label_set(prots, y, ff)
        line = f'  {NICE[ff]:16s}'
        for c, scale in (('context own', 'lambda'), ('context other', 'eps RC')):
            v = []
            for k, s, yy in zip(keys, seqs, Y):
                m = fold_models[('CALVADOS-2', c, int(prots[k]['fold']))]
                pr = m.predict(np.hstack([full(s)[:, POS_COLS], ctx(s, scale)]))
                v.append(spearmanr(pr, yy).correlation if yy.std() > 0 else np.nan)
                rows.append({'protein_id': k[1], 'force_field': NICE[ff], 'condition': f'transfer {c} ({scale})',
                             'spearman': v[-1]})
            line += f'   {c} ({scale}) {np.nanmedian(v):.3f}'
        cal_y = {k[1]: yy for k, yy in zip(*label_set(prots, y, 'CALVADOS-2')[::2])}
        sv = [spearmanr(cal_y[k[1]], yy).correlation for k, yy in zip(keys, Y) if k[1] in cal_y and yy.std() > 0]
        sims[ff] = np.nanmedian(sv)
        line += f'   full model {np.nanmedian(list(cv_full[(ff, "transfer")].values())):.3f}   simulations agree {sims[ff]:.3f}'
        p(line)
        feats.clear()

    # shuffles
    p('\nSimulated shuffles (1,791; models trained on all natural CALVADOS-2 proteins; post hoc)')
    Xn, Yn = Xcal
    SH = read_deg('hub_frames_v3_shuffle')
    sh_ids = sorted(u for u in SH if '_global_' in u)
    sh_ids = sh_ids[:60] if SMOKE else sh_ids
    sh = {}
    for u in sh_ids:
        s = SH[u][0]
        Pm = full(s)[:, POS_COLS]
        sh[u] = {'position': Pm, 'identity': np.hstack([Pm, onehot(s)]), 'context own': np.hstack([Pm, ctx(s, 'lambda')])}
    feats.clear()
    line = '  '
    for c in ('position', 'identity', 'context own'):
        m = SS.hgb().fit(Xn[c], Yn)
        v = []
        for u in sh_ids:
            occ = SH[u][1]
            v.append(spearmanr(m.predict(sh[u][c]), occ).correlation if occ.std() > 0 else np.nan)
            rows.append({'protein_id': u, 'force_field': 'CALVADOS-2 shuffle', 'condition': c, 'spearman': v[-1]})
        line += f'{c} {np.nanmedian(v):.3f}   '
        print(f'  shuffles {c} done ({time.time() - t0:.0f} s)', flush=True)
    p(line + f'(n = {len(sh_ids)}; occupancy: frozen full model 0.889)')
    del Xn, Xcal, sh

    # permutation importance of the context-own model
    p('\nGrouped permutation importance of the context-own model (drop of the median held-out Spearman)')
    rng = np.random.default_rng(0)
    for ff in FFS:
        M, Yc, G, FO = store[ff]
        base = np.zeros(len(Yc), np.float32)
        perm = {g: np.zeros(len(Yc), np.float32) for g in GROUPS}
        for f in range(5):
            te = FO == f
            m = fold_models[(ff, 'context own', f)]
            Mt = M[te]
            base[te] = m.predict(Mt)
            for g, cols in GROUPS.items():
                Mp = Mt.copy()
                Mp[:, cols] = Mt[rng.permutation(len(Mt))][:, cols]
                perm[g][te] = m.predict(Mp)
        b = np.nanmedian(per_protein(base, Yc, G))
        drops = {g: b - np.nanmedian(per_protein(perm[g], Yc, G)) for g in GROUPS}
        p(f'  {NICE[ff]:16s} baseline {b:.3f}: ' + '  '.join(f'{g} {d:+.3f}' for g, d in
                                                             sorted(drops.items(), key=lambda t: -t[1])))
    store.clear()

    # neighbourhood composition around hubs beyond position
    p('\nComposition of the +-5 window around hubs beyond position (top 10% of position-adjusted degree), relative '
      'to the window around all residues; 95% bootstrap intervals over proteins')
    sets = {ff: label_set(prots, y, ff)[1:] for ff in FFS}
    sets['PED'] = ped_degree_sets()
    nrows = []
    lab_w = max(len(c[0]) for c in CLASSES)
    header = f'  {"neighbour class":{lab_w}s}' + ''.join(f'{NICE[s]:>24s}' for s in sets)
    p(header)
    stats = {}
    for name, (seqs_, Y_) in sets.items():
        Hs, Hn, As, An = [], [], [], []
        for s, occ in zip(seqs_, Y_):
            n = len(s)
            r = AT.posadj(np.asarray(occ, float))
            k = max(1, int(round(0.1 * n)))
            hub = np.zeros(n, bool)
            hub[np.argsort(-(r + np.random.default_rng(len(s)).uniform(0, 1e-9, n)))[:k]] = True
            frac = np.zeros((n, len(CLASSES)))
            for ci, (_, letters) in enumerate(CLASSES):
                x = np.array([a in letters for a in s], float)
                c = np.concatenate([[0.0], np.cumsum(x)])
                lo = np.clip(np.arange(n) - 5, 0, n)
                hi = np.clip(np.arange(n) + 6, 0, n)
                cnt = hi - lo - 1
                frac[:, ci] = np.where(cnt > 0, (c[hi] - c[lo] - x) / np.maximum(cnt, 1), np.nan)
            ok = ~np.isnan(frac).any(1)
            Hs.append(frac[hub & ok].sum(0)), Hn.append((hub & ok).sum())
            As.append(frac[ok].sum(0)), An.append(ok.sum())
        Hs, Hn, As, An = map(np.array, (Hs, Hn, As, An))
        est = (Hs.sum(0) / Hn.sum()) / (As.sum(0) / An.sum())
        W = np.stack([np.bincount(rng.integers(0, len(Hn), len(Hn)), minlength=len(Hn)) for _ in range(NBOOT)]).astype(float)
        bs = ((W @ Hs) / (W @ Hn)[:, None]) / ((W @ As) / (W @ An)[:, None])
        lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        stats[name] = (est, lo, hi)
        for ci, (cname, letters) in enumerate(CLASSES):
            nrows.append({'label_set': NICE[name], 'neighbour_class': cname, 'residues': letters,
                          'enrichment': est[ci], 'lo': lo[ci], 'hi': hi[ci], 'proteins': len(Hn)})
    for ci, (cname, _) in enumerate(CLASSES):
        p(f'  {cname:{lab_w}s}' + ''.join(f'{stats[s][0][ci]:>10.2f} ({stats[s][1][ci]:.2f}-{stats[s][2][ci]:.2f})'
                                         for s in sets))
    with open(os.path.join(HERE, f'context_neighbourhood{SUF}.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(nrows[0]))
        w.writeheader()
        w.writerows(nrows)
    with open(os.path.join(HERE, f'context_model_by_ff_per_protein{SUF}.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    p(f'\n({time.time() - t0:.0f} s)')
    open(os.path.join(HERE, f'context_model_by_ff_summary{SUF}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
