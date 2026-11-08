"""
llps_profiles.py — EXPLORATORY, predictions only (no phase-behaviour outcome is read): how uneven are the predicted hub
profiles of well-known phase-separating IDRs, and how far can composition-preserving shuffles move hub unevenness while
charge and hydropathy patterning stay close to the wild type?

Profiles: frozen consensus model (hub_predictor/frozen_models_2026-10-02.pkl) on the IDRome boundaries
(external_data/idrome_clinvar/IDRome_DB.csv). Reference distribution: the stored IDRome profiles of the ClinVar
analysis (clinvar_hub/wt_profiles.csv, same model), compared with IDRs of similar length (0.67 N to 1.5 N).
Unevenness: CV = std / mean of the predicted occupancy, and top-10% share = fraction of the summed predicted occupancy
carried by the top 10% of residues (0.10 for a flat profile).
Shuffles: K random permutations per IDR; SCD (Sawle & Ghosh) and SHD (CALVADOS-2 lambda) per sequence; "matched"
shuffles have SCD and SHD both within 0.25 SD (of the shuffle distribution) of the wild-type values.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python llps_profiles.py  ->  llps_csat/llps_profiles.txt, llps_profiles.csv
"""
import csv
import math
import os
import pickle
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

IDROME = os.path.join(HERE, '..', 'external_data', 'idrome_clinvar', 'IDRome_DB.csv')
WT = os.path.join(HERE, 'clinvar_hub', 'wt_profiles.csv')
CAL = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')
OUT = os.path.join(HERE, 'llps_csat')
K = 300
LLPS = [  # (label, IDRome seq_name): canonical phase-separating IDRs
    ('FUS N-terminal LC + RGG1', 'P35637_1_285'), ('hnRNPA1 LCD', 'P09651_187_372'),
    ('hnRNPA2 LCD', 'P22626_194_353'), ('TDP-43 CTD', 'Q13148_252_414'), ('DDX4 N-terminal IDR', 'Q9NQI0_1_252'),
    ('TAF15 N-terminal LC', 'Q92804_1_234'), ('EWSR1 N-terminal LC', 'Q01844_1_360'), ('TIA1 PrLD', 'P31483_283_386'),
    ('UBQLN2 IDR', 'Q9UHD9_98_581'), ('G3BP1 acidic IDR', 'Q13283_137_342'), ('G3BP1 RGG', 'Q13283_410_466'),
    ('NPM1 IDR', 'P06748_122_247'), ('CAPRIN1 IDR', 'Q14444_247_709'), ('NUP98 FG (N-term)', 'P52948_1_166'),
    ('MED1 IDR', 'Q15648_525_1581'), ('BRD4 IDR', 'O60885_671_1230'), ('Tau', 'P10636_1_758'), ('Elastin', 'P15502_1_786'),
    ('FMRP C-terminal IDR', 'Q06787_442_632'), ('ATXN2 IDR', 'Q99700_467_1313'), ('RBM14 IDR', 'Q96PK6_145_669'),
    ('SFPQ N-terminal IDR', 'P23246_1_285'),
]
MODEL = None
LAM = None


def init():
    global MODEL, LAM
    MODEL = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl'), 'rb'))['consensus']
    LAM = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(CAL))}


def uneven(p):
    p = np.clip(p, 0, None)
    q = p / p.sum()
    k = math.ceil(0.1 * len(p))
    return float(p.std() / p.mean()), float(np.sort(q)[-k:].sum())


def scd(s):
    q = np.array([1 if c in 'KR' else -1 if c in 'DE' else 0 for c in s], float)
    i, j = np.triu_indices(len(s), 1)
    return float((q[i] * q[j] * np.sqrt(j - i)).sum() / len(s))


def shd(s):
    lam = np.array([LAM.get(c, 0.5) for c in s])
    i, j = np.triu_indices(len(s), 1)
    return float(((lam[i] + lam[j]) / (j - i)).sum() / len(s))


WINDOW = 0   # 0 = full shuffles; w > 0 = shuffle only within consecutive windows of w residues


def shuffle(seq, rng):
    if not WINDOW:
        return ''.join(rng.permutation(list(seq)))
    return ''.join(''.join(rng.permutation(list(seq[i:i + WINDOW]))) for i in range(0, len(seq), WINDOW))


def one(job):
    label, name, seq = job
    rng = np.random.default_rng(abs(hash(name)) % 2 ** 32)
    wt = MODEL.predict(T.features(seq, 10))
    cv0, top0 = uneven(wt)
    rows = []
    for _ in range(K):
        s = shuffle(seq, rng)
        cv, top = uneven(MODEL.predict(T.features(s, 10)))
        rows.append((cv, top, scd(s), shd(s)))
    return label, name, seq, wt, cv0, top0, scd(seq), shd(seq), np.array(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    out = []

    def p(s=''):
        out.append(s)
        print(s, flush=True)
    seqs = {r['seq_name']: r['fasta'] for r in csv.DictReader(open(IDROME))}
    ref = defaultdict(list)
    for r in csv.DictReader(open(WT)):
        ref[r['idr']].append(float(r['pred_consensus']))
    refstats = {k: (len(v),) + uneven(np.array(v)) for k, v in ref.items() if len(v) >= 30}
    jobs = [(lab, n, seqs[n]) for lab, n in LLPS if n in seqs]
    with Pool(8, initializer=init) as pool:
        res = pool.map(one, jobs)
    p('EXPLORATORY (predictions only): unevenness of predicted hub profiles of well-known phase-separating IDRs')
    p(f'reference: {len(refstats)} IDRome IDRs (>= 30 residues); percentile among IDRs of similar length (0.67N-1.5N)')
    p(f'shuffles: {K} per IDR; matched = SCD and SHD both within 0.25 SD of the wild type\n')
    p(f'{"IDR":26s} {"N":>5s} {"CV":>5s} {"pct":>4s} {"top10":>6s} {"pct":>4s} | {"shuffle CV 5-95%":>16s} '
      f'{"WT pct":>6s} | {"matched":>7s} {"matched CV range":>17s} {"rho(CV,SCD)":>11s} {"rho(CV,SHD)":>11s}')
    from scipy.stats import spearmanr
    with open(os.path.join(OUT, 'llps_profiles.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['label', 'idr', 'pos', 'aa', 'pred_consensus'])
        for label, name, seq, wt, cv0, top0, scd0, shd0, R in res:
            n = len(seq)
            sim = [v for k, v in refstats.items() if 0.67 * n <= v[0] <= 1.5 * n]
            pc = 100 * np.mean([v[1] < cv0 for v in sim])
            pt = 100 * np.mean([v[2] < top0 for v in sim])
            sd_scd, sd_shd = R[:, 2].std(), R[:, 3].std()
            m = (np.abs(R[:, 2] - scd0) <= 0.25 * sd_scd) & (np.abs(R[:, 3] - shd0) <= 0.25 * sd_shd)
            lo, hi = np.percentile(R[:, 0], [5, 95])
            wtp = 100 * np.mean(R[:, 0] < cv0)
            mr = f'{R[m, 0].min():.3f}-{R[m, 0].max():.3f}' if m.sum() >= 2 else '-'
            p(f'{label:26s} {n:5d} {cv0:5.3f} {pc:4.0f} {top0:6.3f} {pt:4.0f} | {lo:7.3f}-{hi:<7.3f}  {wtp:6.0f} | '
              f'{int(m.sum()):7d} {mr:>17s} {spearmanr(R[:, 0], R[:, 2]).correlation:+11.2f} '
              f'{spearmanr(R[:, 0], R[:, 3]).correlation:+11.2f}')
            w.writerows([label, name, i, a, f'{v:.5f}'] for i, (a, v) in enumerate(zip(seq, wt)))
    p('\nCV = std/mean of predicted occupancy; top10 = share of summed occupancy in the top 10% of residues (flat = 0.10);')
    p('pct = percentile among similar-length IDRome IDRs; WT pct = percentile of the wild-type CV among its own shuffles.')
    open(os.path.join(OUT, 'llps_profiles.txt'), 'w').write('\n'.join(out) + '\n')


def control(n_idr=300, k=60, seed=7):
    """Same wild-type vs own-shuffles percentile of CV for random IDRome IDRs (100-600 residues, LLPS list excluded)."""
    global K
    K = k
    rows = [r for r in csv.DictReader(open(IDROME)) if 100 <= int(r['N']) <= 600]
    skip = {n for _, n in LLPS}
    rng = np.random.default_rng(seed)
    pick = [rows[i] for i in rng.choice(len(rows), n_idr * 2, replace=False) if rows[i]['seq_name'] not in skip][:n_idr]
    jobs = [('control', r['seq_name'], r['fasta']) for r in pick]
    with Pool(8, initializer=init) as pool:
        res = pool.map(one, jobs)
    pct = np.array([100 * np.mean(R[:, 0] < cv0) for _, _, _, _, cv0, _, _, _, R in res])
    out = [f'CONTROL: {len(res)} random IDRome IDRs (100-600 residues), {k} shuffles each',
           f'  wild-type CV percentile among own shuffles: median {np.median(pct):.0f}; >= 95 in {np.mean(pct >= 95):.0%}; '
           f'>= 85 in {np.mean(pct >= 85):.0%}; <= 15 in {np.mean(pct <= 15):.0%}']
    print('\n'.join(out))
    open(os.path.join(OUT, 'llps_profiles_control.txt'), 'w').write('\n'.join(out) + '\n')
    with open(os.path.join(OUT, 'llps_profiles_control.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['idr', 'N', 'wt_cv', 'wt_pct_among_shuffles'])
        w.writerows([n, len(s), f'{cv0:.4f}', f'{q:.1f}'] for (_, n, s, _, cv0, _, _, _, _), q in zip(res, pct))


def local(k_llps=100, k_ctrl=40, n_idr=300):
    """Wild-type CV percentile among window-25 local shuffles, LLPS list (100-600 residues) and random IDRome IDRs."""
    global K, WINDOW
    WINDOW = 25
    from scipy.stats import mannwhitneyu
    seqs = {r['seq_name']: r['fasta'] for r in csv.DictReader(open(IDROME))}
    K = k_llps
    jobs = [(lab, n, seqs[n]) for lab, n in LLPS if n in seqs and 100 <= len(seqs[n]) <= 600]
    with Pool(8, initializer=init) as pool:
        resL = pool.map(one, jobs)
    rows = [r for r in csv.DictReader(open(IDROME)) if 100 <= int(r['N']) <= 600]
    skip = {n for _, n in LLPS}
    rng = np.random.default_rng(7)
    pick = [rows[i] for i in rng.choice(len(rows), n_idr * 2, replace=False) if rows[i]['seq_name'] not in skip][:n_idr]
    K = k_ctrl
    with Pool(8, initializer=init) as pool:
        resC = pool.map(one, [('control', r['seq_name'], r['fasta']) for r in pick])
    pl = np.array([100 * np.mean(R[:, 0] < cv0) for *_, cv0, _, _, _, R in resL])
    pc = np.array([100 * np.mean(R[:, 0] < cv0) for *_, cv0, _, _, _, R in resC])
    out = ['LOCAL SHUFFLES (within consecutive 25-residue windows; local composition kept)']
    for lab, x in (('LLPS drivers', pl), ('random IDRome', pc)):
        out.append(f'  {lab:14s} n = {len(x):3d}: WT CV percentile median {np.median(x):.0f}; >= 95 in {np.mean(x >= 95):.0%}; '
                   f'>= 85 in {np.mean(x >= 85):.0%}')
    out.append(f'  Mann-Whitney LLPS > random: P = {mannwhitneyu(pl, pc, alternative="greater").pvalue:.2g}')
    out += [f'    {r[0]:26s} {q:5.0f}' for r, q in zip(resL, pl)]
    print('\n'.join(out))
    open(os.path.join(OUT, 'llps_profiles_local.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    if '--control' in sys.argv:
        control()
    elif '--local' in sys.argv:
        local()
    else:
        main()
