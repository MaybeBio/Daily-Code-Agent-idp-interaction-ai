"""
clinvar_hub_conservation.py — EXPLORATORY check 2 (2026-10-09): does the effect-test signal survive control for
evolutionary conservation of the variant residue?

Conservation: for every human IDR that carries a variant, its orthologous IDR fragments from IDRome
(idr_orthologs_and_human_*.csv.gz, OMA orthologs) are aligned to the human IDR one at a time by semi-global alignment
(BLOSUM62, linear gap penalty -4, free end gaps; vectorised row by row). At most 60 distinct orthologs per IDR are used
(seeded random choice). Conservation of a residue = fraction of aligned orthologs carrying the same amino acid at the
aligned position (a gap counts as different); residues of IDRs with fewer than 5 orthologs are left missing.
Then the composition + severity model of clinvar_hub_explore2.py is refitted with conservation added.
  ~/.venvs/bender/bin/python clinvar_hub_conservation.py conservation   ->  clinvar_hub/conservation.csv
  ~/.venvs/bender/bin/python clinvar_hub_conservation.py analyse        ->  clinvar_hub/clinvar_hub_conservation.txt
"""
import csv
import gzip
import os
import random
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import clinvar_hub as C  # noqa: E402
import clinvar_hub_explore2 as E2  # noqa: E402

ORTH = os.path.join(C.EXT, 'orthologs')
MAX_ORTH, MIN_ORTH, GAP = 60, 5, -4
IDX = {a: i for i, a in enumerate(E2.ORDER)}
S = np.vstack([E2.M, np.full((1, 20), -4)])           # row 20: any non-standard residue
S = np.hstack([S, np.full((21, 1), -4)])


def enc(seq):
    return np.array([IDX.get(a, 20) for a in seq], np.int64)


def align_map(h, o):
    """Semi-global alignment of ortholog o to human h (free end gaps on both). Returns, for each human position,
    the aligned ortholog residue index or -1."""
    n, m = len(h), len(o)
    H = np.zeros((n + 1, m + 1))
    P = np.zeros((n + 1, m + 1), np.int8)                 # 0 diag, 1 up (gap in ortholog), 2 left (gap in human)
    P[1:, 0] = 1
    P[0, 1:] = 2
    jg = GAP * np.arange(m + 1)
    for i in range(1, n + 1):
        diag = H[i - 1, :-1] + S[h[i - 1], o]
        up = H[i - 1, 1:] + GAP
        tmp = np.empty(m + 1)
        tmp[0] = H[i - 1, 0]                               # free leading gap in ortholog (human residues unaligned)
        tmp[1:] = np.maximum(diag, up)
        src = np.zeros(m + 1, np.int8)
        src[1:] = np.where(diag >= up, 0, 1)
        src[0] = 1
        # left moves: H[i,j] = max(tmp[j], H[i,j-1] + GAP)  ->  cumulative maximum of tmp - jg, plus jg
        acc = np.maximum.accumulate(tmp - jg) + jg
        left = acc > tmp + 1e-9
        H[i] = acc
        P[i] = np.where(left, 2, src)
    # free trailing gaps: best end in last row or last column
    i, j = n, int(np.argmax(H[n]))
    if H[:, m].max() > H[n, j]:
        i, j = int(np.argmax(H[:, m])), m
    out = np.full(n, -1, np.int64)
    while i > 0 and j > 0:
        d = P[i, j]
        if d == 0:
            out[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif d == 1:
            i -= 1
        else:
            j -= 1
    return out


def one(job):
    name, hseq, orths = job
    h = enc(hseq)
    same = np.zeros(len(hseq))
    for oseq in orths:
        o = enc(oseq)
        mp = align_map(h, o)
        ok = mp >= 0
        same[ok] += (o[mp[ok]] == h[ok])
    return name, same / len(orths), len(orths)


def conservation():
    D = C.idrome()
    need = {r['idr'] for r in csv.DictReader(open(os.path.join(C.OUT, 'variants.csv')))}
    orth = defaultdict(set)
    for k in range(4):
        with gzip.open(os.path.join(ORTH, f'idr_orthologs_and_human_{k}.csv.gz'), 'rt') as fh:
            for r in csv.DictReader(fh):
                if r['human'] in need and r['fasta'] != D[r['human']]['seq']:
                    orth[r['human']].add(r['fasta'])
    rng = random.Random(0)
    jobs = []
    for name in sorted(need):
        o = sorted(orth.get(name, ()))
        if len(o) >= MIN_ORTH:
            rng.shuffle(o)
            jobs.append((name, D[name]['seq'], o[:MAX_ORTH]))
    print(f'{len(need)} IDRs with variants; {len(jobs)} have >= {MIN_ORTH} orthologs', flush=True)
    jobs.sort(key=lambda j: -len(j[1]) * len(j[2]))
    out = os.path.join(C.OUT, 'conservation.csv')
    with Pool(10) as pool, open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['idr', 'pos', 'aa', 'identity', 'n_orthologs'])
        for n_done, (name, cons, k) in enumerate(pool.imap_unordered(one, jobs, chunksize=4), 1):
            w.writerows([name, i, a, f'{c:.4f}', k] for i, (a, c) in enumerate(zip(D[name]['seq'], cons)))
            if n_done % 250 == 0:
                print(f'{n_done}/{len(jobs)}', flush=True)
    print('->', out, flush=True)


def analyse():
    import clinvar_hub_explore as E
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    cons = {}
    for r in csv.DictReader(open(os.path.join(C.OUT, 'conservation.csv'))):
        cons[(r['idr'], int(r['pos']))] = float(r['identity'])
    # rebuild the rows of explore2 (same classes, covariates and severity terms)
    D = C.idrome()
    cls = {}
    rd = csv.reader(gzip.open(os.path.join(C.EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(rd)
    ix = {c: i for i, c in enumerate(h)}
    for row in rd:
        s = row[ix['significance']]
        cls[row[0]] = 1 if 'athogenic' in s else 0 if 'enign' in s else None
    rep_cls = {}
    with gzip.open(os.path.join(C.EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
        r2 = csv.reader(fh, delimiter='\t')
        hh = next(r2)
        jx = {c.lstrip('#'): i for i, c in enumerate(hh)}
        for row in r2:
            if row[jx['Assembly']] == 'GRCh38':
                cs = row[jx['ClinicalSignificance']]
                if cs in C.PATH_RC or cs in C.BEN_RC:
                    rep_cls[row[jx['VariationID']]] = 1 if cs in C.PATH_RC else 0
    wt = {}
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'wt_profiles.csv'))):
        wt.setdefault(rr['idr'], []).append(float(rr['pred_consensus']))
    pct = {k: np.array(v).argsort().argsort() / max(1, len(v) - 1) for k, v in wt.items()}
    rows, missing = [], defaultdict(int)
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'variants.csv'))):
        y = cls.get(rr['variant']) if rr['set'] == 'primary' else rep_cls.get(rr['variation_id'])
        if y is None:
            continue
        k, N = int(rr['k']), int(rr['N'])
        c = cons.get((rr['idr'], k))
        if c is None:
            missing[rr['set']] += 1
            continue
        row = dict(set=rr['set'], idr=rr['idr'], gene=rr['gene'], wt=rr['wt'], mut=rr['mut'], N=N,
                   end=min(k, N - 1 - k), y=y, score=pct[rr['idr']][k], dsite=float(rr['delta_site']),
                   blosum=E2.BLOSUM[(rr['wt'], rr['mut'])], grantham=E2.grantham(rr['wt'], rr['mut']), cons=c)
        row.update(E.local_feats(D[rr['idr']]['seq'], k))
        rows.append(row)
    base = ['fGS', 'fG', 'H', 'idrGS', 'blosum', 'grantham']
    p('EXPLORATORY check 2: evolutionary conservation of the variant residue added to the composition + severity model.')
    p(f'Variants without conservation (IDR with < {MIN_ORTH} orthologs) are left out: {dict(missing)}')
    from scipy.stats import spearmanr
    for s in ('primary', 'replication'):
        Sx = [r for r in rows if r['set'] == s]
        mu, sd = np.mean([r['dsite'] for r in Sx]), np.std([r['dsite'] for r in Sx])
        for r in Sx:
            r['dsite_z'] = (r['dsite'] - mu) / sd
        cv = np.array([r['cons'] for r in Sx])
        for r, v in zip(Sx, (cv - cv.mean()) / cv.std()):
            r['cons_z'] = v
        p(f'\n{s} set with conservation: {len(Sx)} variants ({sum(r["y"] for r in Sx)} pathogenic)')
        p(f'  median identity pathogenic {np.median([r["cons"] for r in Sx if r["y"]]):.2f}, benign '
          f'{np.median([r["cons"] for r in Sx if not r["y"]]):.2f}; Spearman(predicted change, identity) = '
          f'{spearmanr([r["dsite"] for r in Sx], cv).correlation:+.2f}, Spearman(score, identity) = '
          f'{spearmanr([r["score"] for r in Sx], cv).correlation:+.2f}')
        for lab, extra, sub in (('composition + severity (as check 1)', base, Sx),
                                ('+ conservation', base + ['cons'], Sx),
                                ('+ conservation, collagens removed', base + ['cons'],
                                 [r for r in Sx if not r['gene'].startswith('COL')])):
            e = E.fit(sub, 'dsite_z', extra, per=1.0)
            a = E.fit(sub, 'score', extra)
            p(f'  {lab:36s} effect {e[0]:.3f} [{e[1]:.3f}, {e[2]:.3f}], P = {e[3]:.1e}   '
              f'position {a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}], P = {a[3]:.1e}   (n = {len(sub)})')
        o = E.fit(Sx, 'cons_z', base, per=1.0)
        p(f'  for scale: conservation alone (per SD, composition + severity controlled) OR {o[0]:.3f} [{o[1]:.3f}, {o[2]:.3f}]')
    open(os.path.join(C.OUT, 'clinvar_hub_conservation.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    {'conservation': conservation, 'analyse': analyse}[sys.argv[1]]()
