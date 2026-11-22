"""
thermal_predict.py — predictions before any simulation: do IDRs of thermophiles, mesophiles and psychrophiles differ in
their predicted heating response?

1. Disordered regions (UniProt 'Region: Disordered', MobiDB-lite) from the reference proteomes of 44 bacteria and fungi with
   known optimal growth temperatures (OGT, Engqvist 2018 / Zenodo 1175609): overlapping regions merged, 30-600 residues,
   standard residues only.
2. Sequence model trained on BENDER (10,282 disordered proteins): ridge regression on quadratic terms of sequence
   descriptors computed here from the sequence alone (length, residue fractions, charge, hydropathy, CALVADOS-2 stickiness,
   pi fraction, kappa, SCD, SHD), targets at 300 K: the calibrated total heating response (CALVADOS-2, Mpipi-Recharged;
   total_response_evidence.load), fixed-water susceptibility chi (CALVADOS-2) and compaction nu (CALVADOS-2). 5-fold CV R^2.
3. Predictions for every region; flagged when outside BENDER's descriptor range (1st-99th percentile of the main ones).
4. Tests, each organism classed by OGT (thermophile >= 40 C, psychrophile <= 21 C, mesophile otherwise):
   A  ortholog-matched: within a clade, regions of the same eggNOG group and position (N-terminal, C-terminal, internal);
      per group the mean prediction over thermophile (or psychrophile) organisms minus the mean over mesophiles;
      Wilcoxon signed-rank over groups
   B  proteome-wide: median prediction per organism against OGT, with clade fixed effects (and Spearman within clades)
   C  composition of regions by class (the residues that carry the forcing)
    python thermal_predict.py [--out DIR]
Output: Thermal adaptation/predictions/ (regions.csv, A_matched.csv, B_organisms.csv, C_composition.csv,
summary.txt). Predictions, not simulations: regions are excised from their proteins and some lie outside BENDER's range.
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import functools  # noqa: E402
import glob  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.model_selection import KFold, cross_val_predict  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import PolynomialFeatures, StandardScaler  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
TA = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Thermal adaptation')
INP = f'{TA}/inputs'
OUT = f'{TA}/predictions'
BIO_CSV = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/BENDER_BIO.csv')
RES_CSV = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))
AA = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6,
          S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
CHARGE = dict(D=-1, E=-1, K=1, R=1)
MIN_LEN, MAX_LEN = 30, 600
THERMO, PSYCHRO = 40, 21   # fungal thermophily is defined from ~40 C; no organism here lies between 37 and 40 C
CLADE = {**dict.fromkeys([274, 271, 277, 186192, 187137, 1299, 310783, 68909], 'Deinococcus-Thermus'),
         **dict.fromkeys([1462, 33940, 33934, 911092, 496866, 1525, 1423, 1402, 1488, 1639, 1215089], 'Bacillota'),
         **dict.fromkeys([28229, 357794, 84980, 334543, 330922, 56812, 228, 562, 70863, 666], 'Pseudomonadota'),
         **dict.fromkeys([29549, 986, 531, 96345], 'Bacteroidota'),
         **dict.fromkeys([2021, 49319, 1902], 'Actinomycetota'),
         **dict.fromkeys([1097, 1092], 'Chlorobiota'),
         **dict.fromkeys([209285, 78579, 5141, 162425], 'Fungi (Sordariomycetes, Eurotiomycetes)'),
         2336: 'Thermotogota (no mesophile)', 146786: 'Cyanobacteria (no mesophile)'}


# ---------------------------------------------------------------- descriptors
def _delta(q, g):
    """Das-Pappu delta for a charge vector q (+1, -1, 0) and blob size g (as bender_reanalysis/sequence_descriptors)."""
    n = len(q)
    pos, neg = (q > 0).astype(float), (q < 0).astype(float)

    def sig(fp, fn):
        tot = fp + fn
        return np.divide((fp - fn) ** 2, tot, out=np.zeros_like(tot), where=tot > 0)
    s_all = sig(np.array([pos.mean()]), np.array([neg.mean()]))[0]
    if n < g:
        return 0.0
    w = np.ones(g)
    return float(np.mean((sig(np.convolve(pos, w, 'valid') / g, np.convolve(neg, w, 'valid') / g) - s_all) ** 2))


@functools.lru_cache(maxsize=None)
def _delta_max(n_pos, n_neg, n_neu, g):
    """Largest delta over arrangements with each charge sign in one block and the neutrals split between the start, the
    middle and the end (sequence_descriptors._delta_max, cached by composition)."""
    best = 0.0
    for first, second in ((1, -1), (-1, 1)):
        for a in range(g + 2):
            for b in range(g + 2):
                for c_rest in (True, False):
                    if a + b > n_neu:
                        continue
                    rest = n_neu - a - b
                    mid, end = b + (0 if c_rest else rest), (rest if c_rest else 0)
                    q = np.concatenate([np.zeros(a), np.full(n_pos if first > 0 else n_neg, first), np.zeros(mid),
                                        np.full(n_neg if first > 0 else n_pos, second), np.zeros(end)])
                    best = max(best, _delta(q, g))
    return best


def kappa(seq, blobs=(5, 6)):
    """Das & Pappu (2013) kappa: mean over blob sizes 5 and 6 of delta / delta_max; NaN without charged residues.
    Identical to bender_reanalysis/sequence_descriptors.kappa (the BENDER_BIO definition). Fixed 2026-10-10: the earlier
    version used one blocked arrangement as delta_max, divided mean delta by mean delta_max and returned NaN for
    sequences with one charge sign, giving values up to 18 and r = 0.41 with BENDER_BIO kappa."""
    q = np.array([CHARGE.get(a, 0) for a in seq], float)
    q = np.sign(q)
    n_pos, n_neg = int((q > 0).sum()), int((q < 0).sum())
    if n_pos + n_neg == 0:
        return np.nan
    vals = []
    for g in blobs:
        dmax = _delta_max(n_pos, n_neg, len(q) - n_pos - n_neg, g)
        vals.append(_delta(q, g) / dmax if dmax > 0 else 0.0)
    return float(np.mean(vals))


def descriptors(seqs):
    lam = pd.read_csv(RES_CSV).set_index('one').lambdas.to_dict()
    kd_n = {a: (v + 4.5) / 9.0 for a, v in KD.items()}
    rows = []
    for s in seqs:
        s = ''.join(c for c in s if c in AA)   # 3 BENDER sequences carry X
        n = len(s)
        f = {f'f_{a}': s.count(a) / n for a in AA}
        q = np.array([CHARGE.get(a, 0) for a in s], float)
        idx = np.arange(n)
        i, j = np.triu_indices(n, 1)
        sep = np.abs(i - j).astype(float)
        h = np.array([kd_n[a] for a in s])
        r = {'logN': np.log(n), **f, 'FCR': np.abs(q).mean(), 'NCPR': q.mean(), 'absNCPR': abs(q.mean()),
             'hydropathy': np.mean([KD[a] for a in s]), 'lambda_avg': np.mean([lam[a] for a in s]),
             'f_pi': f['f_F'] + f['f_W'] + f['f_Y'] + f['f_R'], 'f_aro': f['f_F'] + f['f_W'] + f['f_Y'],
             'f_DEK': f['f_D'] + f['f_E'] + f['f_K'], 'f_hydres': sum(f[f'f_{a}'] for a in 'DEKRNQSTHG'),
             'SCD': (q[i] * q[j] * np.sqrt(sep)).sum() / n, 'SHD': ((h[i] + h[j]) / sep).sum() / n,
             'kappa': kappa(s)}
        del idx
        rows.append(r)
    d = pd.DataFrame(rows)
    d['kappa'] = d.kappa.fillna(d.kappa.median())
    return d


# ---------------------------------------------------------------- regions
def regions():
    ref = pd.read_csv(f'{INP}/reference_proteomes.tsv', sep='\t', header=None,
                      names=['species_taxid', 'upid', 'organism', 'strain_taxid', 'n_proteins'])
    ref = ref[ref.upid.str.startswith('UP')]
    ogt = pd.read_csv(f'{INP}/ogt/temperature_data.tsv', sep='\t').groupby('taxid').temperature.mean()   # some species listed per strain
    rows = []
    for _, o in ref.iterrows():
        t = pd.read_csv(f'{INP}/proteomes/{o.upid}.tsv.gz', sep='\t').fillna('')
        for _, p in t.iterrows():
            spans = [(int(a), int(b)) for a, b, note in re.findall(r'REGION (\d+)\.\.(\d+); /note="([^"]*)"', p.Region)
                     if note == 'Disordered']
            merged = []
            for a, b in sorted(spans):
                if merged and a <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            for a, b in merged:
                s = p.Sequence[a - 1:b]
                if not (MIN_LEN <= len(s) <= MAX_LEN) or set(s) - set(AA):
                    continue
                pos = 'N' if a <= 5 else ('C' if b >= p.Length - 5 else 'internal')
                rows.append({'region_id': f'{p.Entry}_{a}_{b}', 'accession': p.Entry, 'upid': o.upid,
                             'species_taxid': o.species_taxid, 'organism': o.organism.split(' (')[0],
                             'clade': CLADE[o.species_taxid], 'OGT_C': float(ogt[o.species_taxid]),
                             'gene': p['Gene Names (primary)'], 'protein': p['Protein names'][:80],
                             'protein_length': p.Length, 'start': a, 'end': b, 'position': pos,
                             'eggnog': p.eggNOG.strip(';'), 'sequence': s})
    r = pd.DataFrame(rows)
    r['class'] = np.where(r.OGT_C >= THERMO, 'thermophile', np.where(r.OGT_C <= PSYCHRO, 'psychrophile', 'mesophile'))
    return r


# ---------------------------------------------------------------- model
def bender_targets():
    import total_response_evidence as E
    D, _ = E.load()
    bio = pd.read_csv(BIO_CSV, usecols=['UniProt_ID', 'sequence', 'nu'], low_memory=False).set_index('UniProt_ID')
    t = pd.DataFrame({'total_C': D['C'].total, 'chi_C': D['C'].chi}).join(D['R'].total.rename('total_R'), how='left')
    return bio.join(t, how='inner').join(bio.nu.rename('nu_C'), rsuffix='_').drop(columns=['nu'])


FEATS = ['logN'] + [f'f_{a}' for a in AA] + ['FCR', 'NCPR', 'absNCPR', 'hydropathy', 'lambda_avg', 'f_pi', 'f_aro',
                                             'f_DEK', 'f_hydres', 'SCD', 'SHD', 'kappa']
CORE = ['logN', 'FCR', 'absNCPR', 'hydropathy', 'lambda_avg', 'f_pi', 'SCD', 'kappa']


def model():
    return make_pipeline(StandardScaler(), PolynomialFeatures(2, include_bias=False), StandardScaler(),
                         RidgeCV(alphas=np.logspace(-1, 4, 16)))


def main(out):
    os.makedirs(out, exist_ok=True)
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 30)
    B = bender_targets()
    XB = descriptors(B.sequence.tolist())
    XB.index = B.index
    R = regions()
    XR = descriptors(R.sequence.tolist())
    lo, hi = XB[CORE].quantile(0.01), XB[CORE].quantile(0.99)
    R['in_range'] = ((XR[CORE] >= lo) & (XR[CORE] <= hi)).all(axis=1).to_numpy()
    cv_rows = []
    for tgt in ('total_C', 'total_R', 'chi_C', 'nu_C'):
        ok = B[tgt].notna()
        X, y = XB.loc[ok, FEATS].to_numpy(), B.loc[ok, tgt].to_numpy()
        pred = cross_val_predict(model(), X, y, cv=KFold(5, shuffle=True, random_state=0))
        cv_rows.append({'target': tgt, 'n': len(y), 'cv_R2': 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum(),
                        'spearman': stats.spearmanr(y, pred)[0]})
        m = model().fit(X, y)
        R[f'pred_{tgt}'] = m.predict(XR[FEATS].to_numpy())
    cv = pd.DataFrame(cv_rows)
    for c in ('f_DEK', 'f_hydres', 'f_pi', 'FCR', 'hydropathy', 'f_E', 'f_K', 'f_R', 'f_Q', 'f_N', 'f_S', 'f_T', 'f_G', 'f_P', 'f_A'):
        R[c] = XR[c].to_numpy()
    R.drop(columns=['sequence']).to_csv(os.path.join(out, 'regions.csv'), index=False)
    R[['region_id', 'sequence']].to_csv(os.path.join(out, 'region_sequences.csv'), index=False)
    lines = [f'{len(R)} disordered regions ({MIN_LEN}-{MAX_LEN} residues) from {R.upid.nunique()} proteomes; '
             f'{R.in_range.mean():.0%} inside BENDER\'s descriptor range', '', 'sequence model (BENDER, 5-fold CV):',
             cv.round(3).to_string(index=False), '',
             R.groupby('class').agg(regions=('region_id', 'size'), organisms=('upid', 'nunique')).to_string()]
    P = [c for c in R if c.startswith('pred_')]
    Ri = R[R.in_range]
    # A: ortholog-matched
    rows = []
    E = Ri.assign(og=Ri.eggnog.str.split(';')).explode('og')
    E = E[E.og.str.len() > 0]
    for (clade, og, pos), g in E.groupby(['clade', 'og', 'position']):
        meso = g[g['class'] == 'mesophile']
        if meso.upid.nunique() == 0:
            continue
        for cls in ('thermophile', 'psychrophile'):
            x = g[g['class'] == cls]
            if x.upid.nunique() == 0:
                continue
            xo, mo = x.groupby('upid')[P].mean().mean(), meso.groupby('upid')[P].mean().mean()
            rows.append({'clade': clade, 'eggnog': og, 'position': pos, 'class': cls, 'n_org_class': x.upid.nunique(),
                         'n_org_meso': meso.upid.nunique(), **{f'diff_{p[5:]}': xo[p] - mo[p] for p in P}})
    A = pd.DataFrame(rows)
    # one row per protein family and position: rows with identical region sets and predictions (the same family listed
    # at several eggNOG levels) are merged, keeping the alphabetically first id (comment corrected 2026-10-10)
    A = A.sort_values('eggnog').groupby(['clade', 'position', 'class', 'n_org_class', 'n_org_meso',
                                         *[c for c in A if c.startswith('diff_')]], as_index=False).first()
    summ = []
    for cls, g in A.groupby('class'):
        for p in P:
            d = g[f'diff_{p[5:]}']
            summ.append({'class vs mesophile': cls, 'prediction': p[5:], 'n_matched_groups': len(d), 'median_diff': d.median(),
                         'frac_more_negative': (d < 0).mean(), 'wilcoxon_p': stats.wilcoxon(d).pvalue if len(d) > 5 else np.nan})
    As = pd.DataFrame(summ)
    A.to_csv(os.path.join(out, 'A_matched.csv'), index=False)
    As.to_csv(os.path.join(out, 'A_matched_summary.csv'), index=False)
    # B: organisms
    O = Ri.groupby(['upid', 'organism', 'clade', 'class', 'OGT_C']).agg(regions=('region_id', 'size'),
                                                                       **{p: (p, 'median') for p in P}).reset_index()
    O.to_csv(os.path.join(out, 'B_organisms.csv'), index=False)
    brow = []
    for p in P:
        Oc = O[~O.clade.str.contains('no mesophile')]
        Xd = pd.get_dummies(Oc.clade, drop_first=True).astype(float)
        X1 = np.column_stack([np.ones(len(Oc)), Oc.OGT_C, Xd])
        coef, res, *_ = np.linalg.lstsq(X1, Oc[p], rcond=None)
        resid = Oc[p] - X1 @ coef
        se = np.sqrt((resid ** 2).sum() / (len(Oc) - X1.shape[1]) * np.linalg.inv(X1.T @ X1)[1, 1])
        within = Oc.assign(r_o=Oc.OGT_C - Oc.groupby('clade').OGT_C.transform('mean'),
                           r_p=Oc[p] - Oc.groupby('clade')[p].transform('mean'))
        brow.append({'prediction': p[5:], 'organisms': len(Oc), 'slope_per_10C_clade_fixed': 10 * coef[1],
                     'se_per_10C': 10 * se, 'p': 2 * stats.t.sf(abs(coef[1] / se), len(Oc) - X1.shape[1]),
                     'spearman_within_clade': stats.spearmanr(within.r_o, within.r_p)[0]})
    Bs = pd.DataFrame(brow)
    Bs.to_csv(os.path.join(out, 'B_organisms_summary.csv'), index=False)
    # C: composition
    C = Ri.groupby('class')[['f_DEK', 'f_hydres', 'f_pi', 'FCR', 'hydropathy', 'f_E', 'f_K', 'f_R', 'f_Q', 'f_N', 'f_S',
                             'f_T', 'f_G', 'f_P', 'f_A']].median()
    C.to_csv(os.path.join(out, 'C_composition.csv'))
    lines += ['', '== A ortholog-matched (class minus mesophile, same eggNOG group and position, within clade)',
              As.round(4).to_string(index=False),
              '', '== B organism medians vs OGT (clade fixed effects; clades without a mesophile left out)',
              Bs.round(4).to_string(index=False), '', O.sort_values(['clade', 'OGT_C']).round(3).to_string(index=False),
              '', '== C region composition by class (median)', C.round(3).to_string()]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    main(ap.parse_args().out)
