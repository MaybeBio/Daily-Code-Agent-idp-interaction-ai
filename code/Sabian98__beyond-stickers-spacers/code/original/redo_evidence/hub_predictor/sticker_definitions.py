"""
sticker_definitions.py — does the force-field dependence of sticker enrichment among hubs depend on which residues are
called stickers?

Hubs, enrichment and label sets exactly as sticker_spacer_by_ff.py: the hubs of a protein are the top 10% of its residues
by hub occupancy (ties broken by a tiny jitter, numpy rng seeded 0, proteins in data order), and the enrichment of a group
of residue types is its share of hubs divided by its share of all residues, pooled over proteins. 95% intervals come
from 2,000 bootstrap resamples of proteins.
Groups
  fixed lists        F/W/Y/R (canonical stickers, also the four largest CALVADOS-2 lambda), F/W/Y/H (the four largest
                     mean Wang-Frenkel epsilon in both Mpipi models), F/W/Y, Y/F, Y/R, F/W/Y/R/H, F/W/Y/R/K;
                     spacers G/S and G/S/N/Q
  own parameters     the 4, 5, 6 and 10 residue types with the largest interaction parameter of each force field and the
                     4 with the smallest (CALVADOS-2 lambda; mean Wang-Frenkel epsilon over the 20 partner types, Mpipi-GG
                     from its run.in files, Mpipi-Recharged from OpenMpipi data_files/recharged_params.txt), and the
                     Spearman correlation across the 20 types between enrichment and the own parameter
Beyond geometry: on the 507 sequences simulated with all three force fields, the mean excess (simulation minus the
residue-neutral null, position-adjusted z of time-averaged betweenness, sim_null_v3_per_residue.csv) of each group's
residues, whose average over all residues is about zero.
  ~/.venvs/bender/bin/python sticker_definitions.py [data folder]
      ->  sticker_definitions_summary.txt, sticker_definitions.csv, sticker_definitions_by_type.csv,
          sticker_definitions_table.tex
"""
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

import sticker_spacer_by_ff as S
import train_hub_model as T

HERE, EV, AAS = S.HERE, S.EV, S.AAS
CAL_TABLE = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')
FFS = ('CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged')
PED = 'PED'
SETS_ALL = FFS + (PED,)
NICE = {'CALVADOS-2': 'CALVADOS-2', 'MPIPI-GG': 'Mpipi-GG', 'Mpipi-Recharged': 'Mpipi-Recharged', PED: 'PED'}
NBOOT = 2000
FIXED = [('stickers F, W, Y, R (canonical)', 'FWYR'), ('stickers F, W, Y, H', 'FWYH'), ('aromatic F, W, Y', 'FWY'),
         ('Y, F', 'YF'), ('Y, R', 'YR'), ('F, W, Y, R, H', 'FWYRH'), ('F, W, Y, R, K', 'FWYRK'),
         ('spacers G, S (canonical)', 'GS'), ('spacers G, S, N, Q', 'GSNQ')]
OWN = [('own 4 stickiest', 4), ('own 5 stickiest', 5), ('own 6 stickiest', 6), ('own stickier half (10)', 10),
       ('own 4 least sticky', -4)]
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def scales():
    lam = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(CAL_TABLE))}
    gg = json.load(open(os.path.join(EV, 'mpipi_epsilon_table.json')))
    rc = json.load(open(os.path.join(EV, 'mpipi_recharged_parameters.json')))['epsilon_kJ_per_mol']
    return {'CALVADOS-2': {a: lam[a] for a in AAS},
            'MPIPI-GG': {a: float(np.mean([gg[a + b] for b in AAS if a + b in gg])) for a in AAS},
            'Mpipi-Recharged': {a: float(np.mean([rc[a + b] for b in AAS])) for a in AAS}}


def own_group(scale, k):
    order = sorted(AAS, key=lambda a: -scale[a])
    return ''.join(order[:k]) if k > 0 else ''.join(order[k:])


def label_sets(data):
    prots, y = T.load(os.path.join(HERE, data), 'occupancy')
    sets = {}
    for ff in FFS:
        keys = [k for k in prots if k[0] == ff and len(y.get(k, {})) == len(prots[k]['sequence'])]
        seqs = [prots[k]['sequence'] for k in keys]
        sets[ff] = (seqs, [np.array([y[k][i] for i in range(len(s))]) for k, s in zip(keys, seqs)])
    sets[PED] = S.ped_sets()
    return sets


def counts(seqs, Y):
    """Residues and hubs of each type per protein, with the hubs chosen exactly as in sticker_spacer_by_ff.analyse."""
    rng = np.random.default_rng(0)
    ai = {a: i for i, a in enumerate(AAS)}
    R, H = np.zeros((len(seqs), 20)), np.zeros((len(seqs), 20))
    K, N = np.zeros(len(seqs)), np.zeros(len(seqs))
    for j, (s, occ) in enumerate(zip(seqs, Y)):
        n = len(s)
        k = max(1, int(round(0.1 * n)))
        hubs = set(np.argsort(-(occ + rng.uniform(0, 1e-9, n)))[:k])
        for i, a in enumerate(s):
            if a in ai:
                R[j, ai[a]] += 1
                H[j, ai[a]] += i in hubs
        K[j], N[j] = k, n
    return R, H, K, N


def boot_sums(R, H, K, N, seed=1):
    """Bootstrap over proteins: per-type residue and hub sums and total hubs and residues for every resample."""
    rng = np.random.default_rng(seed)
    P = len(K)
    Rs, Hs, Ks, Ns = [], [], [], []
    for start in range(0, NBOOT, 200):
        W = np.stack([np.bincount(rng.integers(0, P, P), minlength=P) for _ in range(min(200, NBOOT - start))]).astype(float)
        Rs.append(W @ R), Hs.append(W @ H), Ks.append(W @ K), Ns.append(W @ N)
    return np.vstack(Rs), np.vstack(Hs), np.concatenate(Ks), np.concatenate(Ns)


def enrich(group, R, H, K, N):
    """Enrichment of a group among hubs for one sum vector (1D) or many bootstrap sums (2D)."""
    cols = [AAS.index(a) for a in group]
    return (H[..., cols].sum(-1) / K) / (R[..., cols].sum(-1) / N)


def excess_data():
    by = defaultdict(lambda: defaultdict(list))         # force field -> protein -> [(aa, excess)]
    for r in csv.DictReader(open(os.path.join(EV, 'sim_null_v3_per_residue.csv'))):
        by[r['force_field']][r['UniProt_ID']].append((r['aa'], float(r['z_sim_fixed']) - float(r['z_null_fixed'])))
    res = {}
    for ff, prot in by.items():
        ids = sorted(prot)
        S_ = np.zeros((len(ids), 20))                   # summed excess per type
        C_ = np.zeros((len(ids), 20))                   # residue count per type
        for j, u in enumerate(ids):
            for a, e in prot[u]:
                if a in AAS:
                    S_[j, AAS.index(a)] += e
                    C_[j, AAS.index(a)] += 1
        res[ff] = (S_, C_)
    return res


def mean_excess(group, S_, C_, W=None):
    cols = [AAS.index(a) for a in group]
    if W is None:
        return S_[:, cols].sum() / C_[:, cols].sum()
    return (W @ S_[:, cols].sum(1)) / (W @ C_[:, cols].sum(1))


def fmt(est, lo, hi, nd=2):
    return f'{est:.{nd}f} ({lo:.{nd}f}-{hi:.{nd}f})'


def main():
    data = sys.argv[1] if len(sys.argv) > 1 else 'data_v3'
    sc = scales()
    sets = label_sets(data)
    p(f'Sticker and spacer enrichment among hubs under different definitions ({data}; hubs = top 10% occupancy per '
      f'protein; 95% bootstrap intervals over proteins, {NBOOT} resamples)')
    p()
    p('Own parameter rankings (largest first):')
    for ff in FFS:
        p(f'  {NICE[ff]:16s} ' + ' '.join(sorted(AAS, key=lambda a: -sc[ff][a])))
    p()
    rows, by_type = [], []
    C, BS = {}, {}
    for name in SETS_ALL:
        seqs, Y = sets[name]
        C[name] = counts(seqs, Y)
        BS[name] = boot_sums(*C[name])
        R, H, K, N = C[name]
        Rb, Hb, Kb, Nb = BS[name]
        for i, a in enumerate(AAS):
            e = enrich(a, R.sum(0), H.sum(0), K.sum(), N.sum())
            b = enrich(a, Rb, Hb, Kb, Nb)
            by_type.append({'label_set': NICE[name], 'aa': a, 'enrichment': e, 'lo': np.percentile(b, 2.5),
                            'hi': np.percentile(b, 97.5), 'residues': int(R[:, i].sum())})
    groups = [(lab, {s: g for s in SETS_ALL}) for lab, g in FIXED]
    groups += [(lab, {s: own_group(sc[s], k) for s in FFS}) for lab, k in OWN]
    p(f'{"group":34s}' + ''.join(f'{NICE[s]:>30s}' for s in SETS_ALL))
    for lab, gmap in groups:
        line = f'{lab:34s}'
        for s in SETS_ALL:
            g = gmap.get(s)
            if g is None:
                line += f'{"-":>30s}'
                continue
            R, H, K, N = C[s]
            e = enrich(g, R.sum(0), H.sum(0), K.sum(), N.sum())
            b = enrich(g, *BS[s])
            lo, hi = np.percentile(b, [2.5, 97.5])
            rows.append({'measure': 'enrichment among hubs', 'group': lab, 'residues': g, 'label_set': NICE[s],
                         'estimate': e, 'lo': lo, 'hi': hi})
            tag = '' if lab.startswith(('stickers F, W, Y, R (', 'spacers G, S (')) or gmap.get(PED) else f' {g}'
            line += f'{fmt(e, lo, hi) + tag:>30s}'
        p(line)
    p()
    p('List-free: Spearman across the 20 residue types between enrichment among hubs and a force field\'s own parameter')
    for s in SETS_ALL:
        e = {r['aa']: r['enrichment'] for r in by_type if r['label_set'] == NICE[s]}
        parts = []
        for ff in FFS:
            rho = spearmanr([sc[ff][a] for a in AAS], [e[a] for a in AAS])[0]
            parts.append(f'vs {NICE[ff]} {rho:+.2f}')
            rows.append({'measure': 'spearman enrichment vs parameter', 'group': f'parameter of {NICE[ff]}',
                         'residues': '', 'label_set': NICE[s], 'estimate': rho, 'lo': '', 'hi': ''})
        p(f'  {NICE[s]:16s} ' + '   '.join(parts))
    p()
    p('Beyond geometry (507 sequences in all three force fields): mean excess over the residue-neutral null of the '
      'group\'s residues (z units; about 0 over all residues)')
    ex = excess_data()
    rng = np.random.default_rng(2)
    allres = {}
    for ff in FFS:
        S_, C_ = ex[ff]
        P = len(S_)
        W = np.stack([np.bincount(rng.integers(0, P, P), minlength=P) for _ in range(NBOOT)]).astype(float)
        allres[ff] = (S_, C_, W)
        p(f'  {NICE[ff]}: {P} proteins, overall mean excess {S_.sum() / C_.sum():+.3f}')
    p(f'{"group":34s}' + ''.join(f'{NICE[s]:>30s}' for s in FFS))
    for lab, gmap in groups:
        line = f'{lab:34s}'
        for ff in FFS:
            S_, C_, W = allres[ff]
            g = gmap[ff]
            e = mean_excess(g, S_, C_)
            lo, hi = np.percentile(mean_excess(g, S_, C_, W), [2.5, 97.5])
            rows.append({'measure': 'mean excess over null', 'group': lab, 'residues': g, 'label_set': NICE[ff],
                         'estimate': e, 'lo': lo, 'hi': hi})
            line += f'{fmt(e, lo, hi):>30s}'
        p(line)
    with open(os.path.join(HERE, 'sticker_definitions.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'sticker_definitions_by_type.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(by_type[0]))
        w.writeheader()
        w.writerows(by_type)
    open(os.path.join(HERE, 'sticker_definitions_summary.txt'), 'w').write('\n'.join(out) + '\n')
    write_table(rows)


TEX_LABELS = {'stickers F, W, Y, R (canonical)': ('Canonical stickers', 'F W Y R'),
              'stickers F, W, Y, H': ('Stickers', 'F W Y H'), 'aromatic F, W, Y': ('Aromatic residues', 'F W Y'),
              'Y, F': ('Stickers', 'Y F'), 'Y, R': ('Stickers', 'Y R'), 'F, W, Y, R, H': ('Stickers', 'F W Y R H'),
              'F, W, Y, R, K': ('Stickers', 'F W Y R K'), 'spacers G, S (canonical)': ('Canonical spacers', 'G S'),
              'spacers G, S, N, Q': ('Spacers', 'G S N Q'), 'own 4 stickiest': ('Own four stickiest', ''),
              'own 5 stickiest': ('Own five stickiest', ''), 'own 6 stickiest': ('Own six stickiest', ''),
              'own stickier half (10)': ('Own stickier half', ''), 'own 4 least sticky': ('Own four least sticky', '')}


def write_table(rows):
    """Appendix table body: enrichment among hubs and mean excess over the null for every group, with values whose
    95% bootstrap interval includes 1 (enrichment) or 0 (excess) set in italics."""
    get = {(r['measure'], r['group'], r['label_set']): r for r in rows}
    lines = []
    for measure, sets_, null in (('enrichment among hubs', SETS_ALL, 1.0), ('mean excess over null', FFS, 0.0)):
        title = 'Enrichment among hubs' if null else 'Mean excess over the residue-neutral null'
        lines.append(r'\midrule' if lines else r'\midrule')
        lines.append(r'\multicolumn{6}{l}{\emph{' + title + r'}} \\')
        for lab, _ in FIXED + OWN:
            name, res = TEX_LABELS[lab]
            cells = []
            for s in SETS_ALL:
                r = get.get((measure, lab, NICE[s]))
                if r is None:
                    cells.append('--' if null else '')
                    continue
                v = f'{r["estimate"]:.2f}'.replace('-', '$-$')
                cells.append(r'\emph{' + v + '}' if r['lo'] <= null <= r['hi'] else v)
            lines.append(f'{name} & {res} & ' + ' & '.join(cells) + r' \\')
    open(os.path.join(HERE, 'sticker_definitions_table.tex'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
