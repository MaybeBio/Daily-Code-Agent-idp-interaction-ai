"""
llps_unevenness_degree — llps_degree_plan.md: copy of llps_unevenness.py with the consensus DEGREE model
(degree_consensus_2026-10-13.pkl) in place of the consensus occupancy model and outputs *_degree; otherwise unchanged.
Original description:
llps_unevenness.py — implements llps_unevenness_plan.md (hashed 2026-10-09 before PhaSePro was downloaded).
Do phase-separation drivers (PhaSePro, human) concentrate predicted hub occupancy by arrangement more than matched
IDRome controls? Predictions only (frozen consensus model).
From PhaSePro only the accession, organism and region boundaries are read.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python llps_unevenness.py  ->  llps_csat/llps_unevenness_summary.txt, .csv
"""
import csv
import json
import math
import os
import pickle
import sys
import zlib
from collections import Counter
from multiprocessing import Pool

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

EXT = os.path.join(HERE, '..', 'external_data')
PHASEPRO = os.path.join(EXT, 'phasepro', 'download_full.json')
IDROME = os.path.join(EXT, 'idrome_clinvar', 'IDRome_DB.csv')
CAL = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')
OUT = os.path.join(HERE, 'llps_csat')
EXPLORATORY = {'P35637', 'P09651', 'P22626', 'Q13148', 'Q9NQI0', 'Q92804', 'Q01844', 'P31483', 'Q9UHD9', 'Q13283',
               'P06748', 'Q14444', 'P52948', 'Q15648', 'O60885', 'P10636', 'P15502', 'Q06787', 'Q99700', 'Q96PK6',
               'P23246'}
K, WINDOW, NCTRL = 100, 25, 5
AA = 'ACDEFGHIKLMNPQRSTVWY'
MODEL = LAM = None


def init():
    global MODEL, LAM
    MODEL = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'degree_consensus_2026-10-13.pkl'), 'rb'))   # degree (llps_degree_plan.md)
    LAM = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(CAL))}


def cv(p):
    p = np.clip(p, 0, None)
    return float(p.std() / p.mean())


def scd(s):
    q = np.array([1 if c in 'KR' else -1 if c in 'DE' else 0 for c in s], float)
    i, j = np.triu_indices(len(s), 1)
    return float((q[i] * q[j] * np.sqrt(j - i)).sum() / len(s))


def shd(s):
    lam = np.array([LAM.get(c, 0.5) for c in s])
    i, j = np.triu_indices(len(s), 1)
    return float(((lam[i] + lam[j]) / (j - i)).sum() / len(s))


def shuffles(seq, rng, local):
    for _ in range(K):
        if local:
            yield ''.join(''.join(rng.permutation(list(seq[i:i + WINDOW]))) for i in range(0, len(seq), WINDOW))
        else:
            yield ''.join(rng.permutation(list(seq)))


def one(job):
    name, seq, patterning = job
    c0 = cv(MODEL.predict(T.features(seq, 10)))
    out = {'name': name, 'cv': c0}
    for local in (True, False):
        rng = np.random.default_rng(zlib.crc32(f'{name}|{int(local)}'.encode()))
        vals, pat = [], []
        for s in shuffles(seq, rng, local):
            vals.append(cv(MODEL.predict(T.features(s, 10))))
            if patterning and local:
                pat.append((scd(s), shd(s)))
        vals = np.array(vals)
        out['pct_local' if local else 'pct_full'] = float(100 * np.mean(vals < c0))
        if patterning and local:
            pat = np.array(pat)
            out['rho_scd'] = float(spearmanr(vals, pat[:, 0]).correlation)
            out['rho_shd'] = float(spearmanr(vals, pat[:, 1]).correlation)
    return out


def comp_features(s):
    n = len(s)
    c = Counter(s)
    f = np.array([c[a] / n for a in AA])
    ent = float(-(f[f > 0] * np.log(f[f > 0])).sum())
    fcr = (c['K'] + c['R'] + c['D'] + c['E']) / n
    ncpr = (c['K'] + c['R'] - c['D'] - c['E']) / n
    return np.array([math.log(n), ent, (c['F'] + c['W'] + c['Y']) / n, (c['G'] + c['S'] + c['Q'] + c['N']) / n, fcr, ncpr])


def regions(b):
    out = []
    for part in b.replace(',', ';').split(';'):
        part = part.strip()
        if '-' in part:
            a, z = part.split('-')[:2]
            out.append((int(a), int(z)))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)
    pp = json.load(open(PHASEPRO))
    ents = [(e['accession'], e['organism'], e['boundaries']) for e in pp.values()]
    all_acc = {a for a, _, _ in ents}
    human = [(a, b) for a, o, b in ents if o == 'Homo sapiens']
    idrome = [r for r in csv.DictReader(open(IDROME)) if 50 <= int(r['N']) <= 800]
    by_acc = {}
    for r in idrome:
        by_acc.setdefault(r['UniProt_ID'], []).append(r)
    drivers, why = [], Counter()
    for acc, b in human:
        if acc in EXPLORATORY:
            why['exploratory'] += 1
            continue
        reg = regions(b)
        best, bo = None, 0
        for r in by_acc.get(acc, []):
            f, l = int(r['first']), int(r['last'])
            ov = sum(max(0, min(l, z) - max(f, a) + 1) for a, z in reg)
            if ov > bo:
                best, bo = r, ov
        if best is None or bo < 30:
            why['no IDRome IDR of 50-800 residues with >= 30 overlap'] += 1
            continue
        drivers.append((acc, best, bo))
    pool = [r for r in idrome if r['UniProt_ID'] not in all_acc and r['UniProt_ID'] not in EXPLORATORY]
    Xp = np.array([comp_features(r['fasta']) for r in pool])
    mu, sd = Xp.mean(0), Xp.std(0)
    Zp = (Xp - mu) / sd
    used = np.zeros(len(pool), bool)
    order = np.random.default_rng(0).permutation(len(drivers))
    match = {}
    for k in order:
        acc, r, _ = drivers[k]
        z = (comp_features(r['fasta']) - mu) / sd
        dist = np.linalg.norm(Zp - z, axis=1)
        dist[used] = np.inf
        idx = np.argsort(dist)[:NCTRL]
        used[idx] = True
        match[acc] = [pool[i] for i in idx]
    p('LLPS drivers vs matched controls: unevenness of predicted hubs (llps_unevenness_plan.md, pre-registered)')
    p(f'PhaSePro: {len(ents)} entries, {len(human)} human; excluded {dict(why)}; {len(drivers)} drivers tested, '
      f'{NCTRL} matched controls each from {len(pool)} candidate IDRs')
    jobs = [(r['seq_name'], r['fasta'], True) for _, r, _ in drivers]
    jobs += [(c['seq_name'], c['fasta'], False) for acc, _, _ in drivers for c in match[acc]]
    with Pool(8, initializer=init) as pl:
        res = {o['name']: o for o in pl.imap_unordered(one, jobs, chunksize=2)}
    ent_med = float(np.median([comp_features(r['fasta'])[1] for r in idrome]))
    rows = []
    for acc, r, ov in drivers:
        ctrl = [res[c['seq_name']] for c in match[acc]]
        o = res[r['seq_name']]
        rows.append(dict(acc=acc, idr=r['seq_name'], n=int(r['N']), ov=ov, ent=comp_features(r['fasta'])[1],
                         loc=o['pct_local'], full=o['pct_full'], cloc=np.mean([c['pct_local'] for c in ctrl]),
                         cfull=np.mean([c['pct_full'] for c in ctrl]), ctrl_loc=[c['pct_local'] for c in ctrl],
                         ctrl_full=[c['pct_full'] for c in ctrl], rscd=o['rho_scd'], rshd=o['rho_shd'],
                         ctrl_names=[c['seq_name'] for c in match[acc]]))
    d = np.array([r['loc'] - r['cloc'] for r in rows])
    w = wilcoxon(d, alternative='greater')
    dl = np.array([r['loc'] for r in rows])
    cl = np.concatenate([r['ctrl_loc'] for r in rows])
    p(f'\nPRIMARY (local shuffles, 25-residue windows): driver minus matched-control percentile, median {np.median(d):+.1f}, '
      f'positive in {int((d > 0).sum())}/{len(d)}; one-sided Wilcoxon signed-rank P = {w.pvalue:.2g}')
    p(f'  decision (P < 0.05): {"confirmed" if w.pvalue < 0.05 else "not confirmed"}')
    p(f'  drivers: median percentile {np.median(dl):.0f}, >= 95 in {np.mean(dl >= 95):.0%}; controls: median {np.median(cl):.0f}, '
      f'>= 95 in {np.mean(cl >= 95):.0%}')
    p('\nSecondary')
    df = np.array([r['full'] - r['cfull'] for r in rows])
    p(f'  full shuffles: median difference {np.median(df):+.1f}, positive in {int((df > 0).sum())}/{len(df)}, '
      f'Wilcoxon P = {wilcoxon(df, alternative="greater").pvalue:.2g}; drivers median {np.median([r["full"] for r in rows]):.0f}, '
      f'controls median {np.median(np.concatenate([r["ctrl_full"] for r in rows])):.0f}')
    lc = [i for i, r in enumerate(rows) if r['ent'] < ent_med]
    if len(lc) >= 5:
        p(f'  low-complexity subset (composition entropy < IDRome median {ent_med:.2f}, n = {len(lc)}): median difference '
          f'{np.median(d[lc]):+.1f}, Wilcoxon P = {wilcoxon(d[lc], alternative="greater").pvalue:.2g}')
    else:
        p(f'  low-complexity subset: only {len(lc)} drivers')
    p(f'  pooled Mann-Whitney drivers > controls (local): P = {mannwhitneyu(dl, cl, alternative="greater").pvalue:.2g}')
    p(f'  across local shuffles of each driver: median Spearman(unevenness, SHD) {np.nanmedian([r["rshd"] for r in rows]):+.2f}, '
      f'(unevenness, SCD) {np.nanmedian([r["rscd"] for r in rows]):+.2f}')
    p('\nPer driver (percentiles among own shuffles; controls = mean of 5 matched)')
    p(f'  {"acc":8s} {"IDR":22s} {"N":>4s} {"ovl":>4s} {"local":>5s} {"ctrl":>5s} {"full":>5s} {"ctrl":>5s}')
    for r in sorted(rows, key=lambda r: -r['loc']):
        p(f'  {r["acc"]:8s} {r["idr"]:22s} {r["n"]:4d} {r["ov"]:4d} {r["loc"]:5.0f} {r["cloc"]:5.0f} {r["full"]:5.0f} {r["cfull"]:5.0f}')
    with open(os.path.join(OUT, 'llps_unevenness_degree.csv'), 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['accession', 'idr', 'N', 'overlap', 'pct_local', 'ctrl_pct_local_mean', 'pct_full', 'ctrl_pct_full_mean',
                     'controls'])
        for r in rows:
            wr.writerow([r['acc'], r['idr'], r['n'], r['ov'], f'{r["loc"]:.1f}', f'{r["cloc"]:.1f}', f'{r["full"]:.1f}',
                         f'{r["cfull"]:.1f}', ';'.join(r['ctrl_names'])])
    open(os.path.join(OUT, 'llps_unevenness_summary_degree.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
