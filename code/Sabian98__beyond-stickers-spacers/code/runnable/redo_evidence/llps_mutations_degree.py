"""
llps_mutations_degree — llps_degree_plan.md: copy of llps_mutations.py with the consensus DEGREE model
(degree_consensus_2026-10-13.pkl) in place of the consensus occupancy model and outputs *_degree; otherwise unchanged.
Original description:
llps_mutations.py — implements llps_mutation_plan.md (hashed 2026-10-09 before any mutant prediction):
are disease mutations that alter condensates (PhaSePro "disease" field) placed at hub-sensitive positions, compared
with the same substitution at the other positions of the same IDR that carry the same wild-type residue?
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python llps_mutations.py  ->  llps_csat/llps_mutations.txt, .csv
"""
import csv
import difflib
import json
import os
import pickle
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from scipy.stats import binomtest, mannwhitneyu, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

EXT = os.path.join(HERE, '..', 'external_data')
PHASEPRO = os.path.join(EXT, 'phasepro', 'download_full.json')
IDROME = os.path.join(EXT, 'idrome_clinvar', 'IDRome_DB.csv')
DESIGN = os.path.join(HERE, 'gaivi', 'hub_design_table.csv')
OUT = os.path.join(HERE, 'llps_csat')
AGG = {('P35637', None), ('P22626-2', 'P298L'), ('Q9UHD9', None), ('P10636-8', None), ('P31483', 'P362L'),
       ('Q13148', 'Q331K')}   # aggregation or loss of dynamics in the notes (plan, secondary)
MODEL = None


def init():
    global MODEL
    MODEL = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'degree_consensus_2026-10-13.pkl'), 'rb'))   # degree (llps_degree_plan.md)


def pred(s):
    return MODEL.predict(T.features(s, 10))


def local(delta, i):
    return float(delta[max(0, i - 3):i + 4].mean())


def cv(p):
    p = np.clip(p, 0, None)
    return float(p.std() / p.mean())


def evaluate(job):
    """job: (key, seq, i, mut). Response of the disease mutation and of the same substitution at the reference positions."""
    key, seq, i, mut = job
    wt = seq[i]
    p0 = pred(seq)
    cv0 = cv(p0)
    out = []
    for j in [i] + [j for j, c in enumerate(seq) if c == wt and j != i]:
        s = seq[:j] + mut + seq[j + 1:]
        p1 = pred(s)
        out.append((j, local(p1 - p0, j), cv(p1) - cv0))
    return key, out


def parse_mutations():
    d = json.load(open(PHASEPRO))
    muts, skipped = [], []
    for e in d.values():
        if e['organism'] != 'Homo sapiens':
            continue
        v = (e.get('disease') or '').strip()
        if not v or v.startswith('Not known'):
            continue
        for line in v.split('\n'):
            f = [x.strip() for x in line.split('|')]
            if not f or not f[0]:
                continue
            m = re.fullmatch(r'([A-Z])(\d+)([A-Z])', f[0])
            if m:
                muts.append(dict(acc=e['accession'], name=f[0], wt=m.group(1), pos=int(m.group(2)), mut=m.group(3),
                                 label=f[4] if len(f) > 4 else ''))
            else:
                skipped.append((e['accession'], f[0], 'not a point substitution'))
    return muts, skipped


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)
    muts, skipped = parse_mutations()
    idrs = defaultdict(list)
    for r in csv.DictReader(open(IDROME)):
        idrs[r['UniProt_ID']].append(r)
    tau441 = [r['sequence'] for r in csv.DictReader(open(DESIGN)) if r['id'] == 'PRE_Tau441'][0]
    tau_can = [r for r in idrs['P10636'] if r['seq_name'] == 'P10636_1_758'][0]['fasta']
    sm = difflib.SequenceMatcher(None, tau441, tau_can, autojunk=False)
    tau_map = {b.a + k + 1: b.b + k + 1 for b in sm.get_matching_blocks() for k in range(b.size)}
    kept = []
    for m in muts:
        base = m['acc'].split('-')[0]
        pos = m['pos']
        if m['acc'] == 'P22626-2':
            pos = pos + 12 if pos > 2 else pos
        elif m['acc'] == 'P10636-8':
            pos = tau_map.get(pos)
            if pos is None:
                skipped.append((m['acc'], m['name'], 'no aligned tau position'))
                continue
        hit = [r for r in idrs.get(base, []) if int(r['first']) <= pos <= int(r['last'])]
        if not hit:
            skipped.append((m['acc'], m['name'], 'not inside an IDRome IDR'))
            continue
        r = hit[0]
        i = pos - int(r['first'])
        if r['fasta'][i] != m['wt']:
            skipped.append((m['acc'], m['name'], f'wild-type mismatch ({r["fasta"][i]} at mapped {pos})'))
            continue
        kept.append(dict(m, idr=r['seq_name'], seq=r['fasta'], i=i, mpos=pos))
    jobs = [(k, m['seq'], m['i'], m['mut']) for k, m in enumerate(kept)]
    with Pool(6, initializer=init) as pl:
        res = dict(pl.imap_unordered(evaluate, jobs))
    rows = []
    for k, m in enumerate(kept):
        out = res[k]
        _, loc, dcv = out[0]
        ref = np.array([x[1] for x in out[1:]])
        refc = np.array([x[2] for x in out[1:]])
        if len(ref) == 0:
            skipped.append((m['acc'], m['name'], 'no reference position'))
            continue

        def pct(x, r):
            return float(100 * (np.mean(r < x) + 0.5 * np.mean(r == x)))
        agg = (m['acc'], None) in AGG or (m['acc'], m['name']) in AGG
        rows.append(dict(m, loc=loc, nref=len(ref), pabs=pct(abs(loc), np.abs(ref)), psign=pct(loc, ref),
                         pcv=pct(abs(dcv), np.abs(refc)), dcv=dcv, agg=agg))
    p('Disease mutations that alter condensates (PhaSePro) at hub-sensitive positions? (llps_mutation_plan.md)')
    p(f'{len(muts)} point substitutions in human entries; {len(rows)} kept')
    for a, n, why in skipped:
        p(f'  excluded {a} {n}: {why}')
    x = np.array([r['pabs'] for r in rows])
    w = wilcoxon(x - 50, alternative='greater')
    p(f'\nPRIMARY: percentile of |local response| among the same substitution at same-residue positions: median '
      f'{np.median(x):.0f}, above 50 in {int((x > 50).sum())}/{len(x)}; one-sided Wilcoxon P = {w.pvalue:.2g}')
    p(f'  decision (P < 0.05): {"confirmed" if w.pvalue < 0.05 else "not confirmed"}')
    p('\nSecondary')
    s = np.array([r['psign'] for r in rows])
    p(f'  signed local response percentile: median {np.median(s):.0f} (> 50 = towards higher predicted occupancy)')
    prot = defaultdict(list)
    for r in rows:
        prot[r['acc'].split('-')[0]].append(r['pabs'])
    med = [np.median(v) for v in prot.values()]
    p(f'  protein level: median percentile above 50 in {sum(v > 50 for v in med)}/{len(med)} proteins, sign test P = '
      f'{binomtest(sum(v > 50 for v in med), len(med), 0.5, "greater").pvalue:.2g}')
    a = [r['pabs'] for r in rows if r['agg']]
    b = [r['pabs'] for r in rows if not r['agg']]
    if a and b:
        p(f'  aggregation / loss of dynamics in notes: median {np.median(a):.0f} (n = {len(a)}); others {np.median(b):.0f} '
          f'(n = {len(b)}); Mann-Whitney P = {mannwhitneyu(a, b, alternative="two-sided").pvalue:.2g}')
    c = np.array([r['pcv'] for r in rows])
    p(f'  |change of hub unevenness| percentile: median {np.median(c):.0f}, Wilcoxon P = '
      f'{wilcoxon(c - 50, alternative="greater").pvalue:.2g}')
    p('\nPer mutation')
    p(f'  {"protein":9s} {"mutation":9s} {"IDR":18s} {"pos":>4s} {"label":10s} {"local":>8s} {"|pct|":>5s} {"sign":>5s} '
      f'{"ref":>4s} {"cv pct":>6s} agg')
    for r in rows:
        p(f'  {r["acc"]:9s} {r["name"]:9s} {r["idr"]:18s} {r["mpos"]:4d} {r["label"][:10]:10s} {r["loc"]:+8.4f} '
          f'{r["pabs"]:5.0f} {r["psign"]:5.0f} {r["nref"]:4d} {r["pcv"]:6.0f} {"yes" if r["agg"] else ""}')
    with open(os.path.join(OUT, 'llps_mutations_degree.csv'), 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['accession', 'mutation', 'idr', 'mapped_pos', 'label', 'local_response', 'pct_abs', 'pct_signed',
                     'n_ref', 'pct_abs_cv', 'aggregation_note'])
        for r in rows:
            wr.writerow([r['acc'], r['name'], r['idr'], r['mpos'], r['label'], f'{r["loc"]:.5f}', f'{r["pabs"]:.1f}',
                         f'{r["psign"]:.1f}', r['nref'], f'{r["pcv"]:.1f}', int(r['agg'])])
    open(os.path.join(OUT, 'llps_mutations_degree.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
