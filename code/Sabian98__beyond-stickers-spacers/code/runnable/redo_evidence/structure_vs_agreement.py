"""
structure_vs_agreement.py — EXPLORATORY (2026-10-09, written after tab:design_exp was known; not pre-registered):
does the agreement of the frozen model (and of the Mpipi-Recharged simulation) with the experimental PED ensembles fall
where the ensembles carry residual secondary structure, which coarse-grained models without backbone hydrogen bonds
cannot represent?

Constructs and ensembles exactly as in design_response/design_vs_experiment.txt (every PED ensemble whose sequence,
phosphoresidues read as S/T/Y, equals the construct; residues matched by difflib; values averaged over ensembles).
Secondary structure per residue, either by DSSP (dssp_numpy.py, Kabsch-Sander from the deposited N, CA, C, O; helix =
H, G or I, strand = E or B; default) or from the Calpha trace (--method calpha, helicity_allatom.one: helical/extended
quadruple criteria, a residue is helical in a model when it lies in two consecutive helical quadruples).
Ensemble-level relative ASA and DSSP entropy from the PED catalogue (only these fields and the entry ID are read).

Reported:
  1. per construct: mean helix fraction, fraction of residues with helix >= 0.3, max helix, extended fraction, rASA,
     DSSP entropy, and the position-partial Spearman of model and simulation with experiment;
  2. across the 20 constructs built from experimental data: Spearman of model / simulation agreement and of the
     model - simulation gap against each structure measure;
  3. within constructs: model and simulation agreement over residues with helix < 0.2 only, against all residues, and
     the Spearman between the rank error (experimental minus predicted rank) and helicity.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python structure_vs_agreement.py [--method dssp|calpha]
      [--labels labels_design_calvados --name _CALVADOS-2]   (default: Mpipi-Recharged design labels)
      ->  design_response/structure_vs_agreement_{dssp,calpha}.txt
"""
import argparse
import csv
import difflib
import json
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, 'hub_predictor')]
import dssp_numpy as DS  # noqa: E402
import helicity_allatom as H  # noqa: E402
import ped_null as pn  # noqa: E402

CAT = (__import__('os').environ['BSS_DATA'] + '/ZOIDBERG-ICML 2027/experimental_data/raw/ped/ped_all_entries.json')
LOCAL = os.path.join(HERE, 'ped_local')
EXTRA = {'PED00001': os.path.join(LOCAL, 'PED00001/PED00001e003.pdb'),
         'PED00158': os.path.join(LOCAL, 'PED00158/PED00158e006.pdb'),
         'PED00017': 'PED00017/PED00017e001.pdb', 'PED00154': 'PED00154/PED00154e006.pdb',
         'PED00216': 'PED00216/PED00216e000.pdb', 'PED00501': 'PED00501/pdbfile.pdb',
         'PED00502': 'PED00502/pdbfile.pdb', 'PED00503': 'PED00503/pdbfile.pdb'}
MD = {'PED_Ab40', 'PED_ACTR'}
PH = {'x': 'S', 'y': 'T', 'z': 'Y'}
_ORIG_FILES_FOR = H.files_for
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def partial(x, y):
    n = len(x)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    return H.partial(x, y, z)


def files_for(e):
    if e in EXTRA:
        return EXTRA[e]
    return _ORIG_FILES_FOR(e)


def helix(e):
    """helicity_allatom.one with this script's file lookup (absolute paths pass through os.path.join unchanged)."""
    H.files_for = files_for
    return H.one(e)


def helix_dssp(e):
    """DSSP fractions over the same models as the labels (models of the first model's length, at most 1,000)."""
    models, names = [], []
    for f in files_for(e).split(';'):
        m, nm = DS.read_backbone(os.path.join(pn.PED, f))
        names = names or nm
        models += m
    models = [m for m in models if len(m) == len(models[0])]
    if len(models) > 1000:
        models = models[::int(np.ceil(len(models) / 1000))]
    S = np.array([list(DS.assign(m, names)) for m in models])
    return e, names, np.isin(S, list('HGI')).mean(0), np.isin(S, list('EB')).mean(0), len(models), None


def catalogue(entries):
    """Only entry_id, relative_asa_mean and entropy_dssp_mean are taken from the catalogue."""
    res = {}
    for x in json.load(open(CAT)):
        if x.get('entry_id') in entries:
            ch = [c for en in x.get('ensembles', []) for c in en.get('chains', [])]
            ra = [c.get('relative_asa_mean') for c in ch if c.get('relative_asa_mean') is not None]
            en = [c.get('entropy_dssp_mean') for c in ch if c.get('entropy_dssp_mean') is not None]
            res[x['entry_id']] = (float(np.mean(ra)) if ra else np.nan, float(np.mean(en)) if en else np.nan)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--method', choices=('dssp', 'calpha'), default='dssp')
    ap.add_argument('--labels', default='labels_design_recharged')
    ap.add_argument('--name', default='')
    args = ap.parse_args()
    method = args.method
    occ, seq = defaultdict(list), defaultdict(str)
    for f in ('experimental_occupancy.csv', 'fresh_ped_occupancy.csv', 'ped_extra_occupancy.csv'):
        for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', f))):
            if r['entry'] in occ and len(occ[r['entry']]) > int(r['pos']):
                continue
            occ[r['entry']].append(float(r['occ_exp']))
            seq[r['entry']] += r['aa']
    plain = {e: ''.join(PH.get(c, c) for c in s) for e, s in seq.items()}
    sim = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, args.labels, 'per_residue.csv'))):
        sim[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
    pred = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02.csv'))):
        pred[r['id']].append(float(r['consensus']))
    import pickle
    import train_hub_model as Tm
    model = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl'), 'rb'))['consensus']
    T = {r['id']: r for r in csv.DictReader(open(os.path.join(HERE, 'gaivi', 'hub_design_table.csv')))
         if r['category'] == 'ped'}
    plan = {}
    for d, t in T.items():
        target = plain[t['parent']]
        plan[d] = [e for e in occ if plain[e] == target]
    entries = sorted({e for v in plan.values() for e in v})
    with Pool(8) as pool:
        hx = {r[0]: r for r in pool.imap_unordered(helix_dssp if method == 'dssp' else helix, entries)}
    cat = catalogue(set(entries))

    rows = []
    for d, t in sorted(T.items()):
        s = t['sequence']
        S = np.array([sim[d][i] for i in range(len(s))])
        P = np.array(pred[d]) if len(pred[d]) == len(s) else model.predict(Tm.features(s, 10))
        hel_sum, ext_sum, cnt = np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s))
        res, ra, en = [], [], []
        for e in plan[d]:
            q = ''.join('E' if (c in 'xyz' and d == 'PED_pSic1E') else PH.get(c, c) for c in seq[e])
            sm = difflib.SequenceMatcher(None, q, s, autojunk=False)
            pr = [(b.a + k, b.b + k) for b in sm.get_matching_blocks() for k in range(b.size)]
            if len(pr) < 0.8 * len(s):
                continue
            ia, ib = np.array([a for a, _ in pr]), np.array([b for _, b in pr])
            X = np.array(occ[e])[ia]
            _, hseq, h, x, nmod, _ = hx[e]
            if len(h) != len(occ[e]):
                p(f'  note: {e} helix length {len(h)} vs occupancy length {len(occ[e])}; skipped for structure')
                continue
            hh, xx = h[ia], x[ia]
            hel_sum[ib] += hh
            ext_sum[ib] += xx
            cnt[ib] += 1
            keep = hh < 0.2
            ok = keep.sum() >= 15
            err = rankdata(X) - rankdata(P[ib])
            res.append([partial(P[ib], X), partial(S[ib], X),
                        partial(P[ib][keep], X[keep]) if ok else np.nan,
                        partial(S[ib][keep], X[keep]) if ok else np.nan,
                        spearmanr(err, hh).correlation if hh.std() > 0 else np.nan])
            if e in cat:
                ra.append(cat[e][0])
                en.append(cat[e][1])
        if not res:
            continue
        m = np.nanmean(np.array(res, float), axis=0)
        c = np.maximum(cnt, 1)
        hel, ext = hel_sum / c, ext_sum / c
        rows.append(dict(d=d, e=t['parent'], n=len(s), ens=len(res), hel=hel.mean(), f30=(hel >= 0.3).mean(),
                         hmax=hel.max(), ext=ext.mean(), rasa=np.nanmean(ra) if ra else np.nan,
                         ent=np.nanmean(en) if en else np.nan, model=m[0], sim=m[1], model_nh=m[2], sim_nh=m[3],
                         err_h=m[4], md=d in MD))

    p(f'EXPLORATORY: residual secondary structure ({method}) of the experimental PED ensembles and agreement with the model')
    p('(not pre-registered; written after tab:design_exp was known)\n')
    p(f'{"construct":11s} {"PED":9s} {"n":>4s} {"ens":>3s} {"helix":>6s} {"h>=.3":>6s} {"hmax":>5s} {"ext":>5s} '
      f'{"rASA":>5s} {"Hdssp":>5s}  {"model":>6s} {"sim":>6s}  {"model nh":>8s} {"sim nh":>7s}  {"err~hel":>7s}')
    for r in sorted(rows, key=lambda r: r['hel']):
        p(f'{r["d"]:11s} {r["e"]:9s} {r["n"]:4d} {r["ens"]:3d} {r["hel"]:6.3f} {r["f30"]:6.2f} {r["hmax"]:5.2f} '
          f'{r["ext"]:5.2f} {r["rasa"]:5.2f} {r["ent"]:5.2f}  {r["model"]:+6.2f} {r["sim"]:+6.2f}  '
          f'{r["model_nh"]:+8.2f} {r["sim_nh"]:+7.2f}  {r["err_h"]:+7.2f}' + ('  (MD)' if r['md'] else ''))
    p('\nhelix = mean helix fraction over residues and ensembles; h>=.3 = fraction of residues helical in >= 30% of models;')
    p('model / sim = position-partial Spearman with experiment (as tab:design_exp); nh = same over residues with helix < 0.2;')
    p('err~hel = Spearman over residues between (experimental rank - predicted rank) and helicity.')

    E = [r for r in rows if not r['md']]
    p(f'\nAcross the {len(E)} constructs built from experimental data (Spearman, n = {len(E)}):')
    for key, lab in (('hel', 'mean helix fraction'), ('f30', 'fraction of residues helix >= 0.3'),
                     ('hmax', 'max helix'), ('ext', 'extended fraction'), ('rasa', 'relative ASA (PED)'),
                     ('ent', 'DSSP entropy (PED)'), ('n', 'length')):
        v = np.array([r[key] for r in E], float)
        ok = ~np.isnan(v)
        cm = spearmanr(v[ok], [r['model'] for r, o in zip(E, ok) if o])
        cs = spearmanr(v[ok], [r['sim'] for r, o in zip(E, ok) if o])
        cg = spearmanr(v[ok], [r['model'] - r['sim'] for r, o in zip(E, ok) if o])
        p(f'  {lab:36s} vs model {cm.correlation:+.2f} (P = {cm.pvalue:.2g})   vs simulation {cs.correlation:+.2f} '
          f'(P = {cs.pvalue:.2g})   vs model - simulation {cg.correlation:+.2f} (P = {cg.pvalue:.2g})')
    hi = [r for r in E if r['f30'] >= 0.1]
    lo = [r for r in E if r['f30'] < 0.1]
    p(f'\n  constructs with >= 10% helical residues (n = {len(hi)}: {", ".join(r["d"][4:] for r in hi)})')
    p(f'    median model {np.median([r["model"] for r in hi]):+.2f}, simulation {np.median([r["sim"] for r in hi]):+.2f}')
    p(f'  constructs with < 10% helical residues (n = {len(lo)})')
    p(f'    median model {np.median([r["model"] for r in lo]):+.2f}, simulation {np.median([r["sim"] for r in lo]):+.2f}')
    w = [r for r in E if not np.isnan(r['model_nh'])]
    p(f'\nWithin constructs, removing residues with helix >= 0.2 ({len(w)} constructs with >= 15 residues left):')
    p(f'  median model all residues {np.median([r["model"] for r in w]):+.2f} -> non-helical {np.median([r["model_nh"] for r in w]):+.2f}; '
      f'simulation {np.median([r["sim"] for r in w]):+.2f} -> {np.median([r["sim_nh"] for r in w]):+.2f}')
    p(f'  median err~helix over constructs {np.nanmedian([r["err_h"] for r in E]):+.2f} '
      f'(positive = helical residues are more central in experiment than predicted)')
    os.makedirs(os.path.join(HERE, 'design_response'), exist_ok=True)
    open(os.path.join(HERE, 'design_response', f'structure_vs_agreement_{method}{args.name}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
