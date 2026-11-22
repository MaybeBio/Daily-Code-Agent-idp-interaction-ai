"""
ogt_response.py — do BENDER IDPs from organisms that grow at higher temperature respond differently to heat? Simulated
responses (not predictions): the test of Thermal adaptation/PREDICTIONS.md with existing BENDER data.

BENDER chains are matched to optimal growth temperatures (OGT; Engqvist 2018, Zenodo 1175609) by NCBI taxid, else by
species name (genus + species). Per model (CALVADOS-2, Mpipi-Recharged; Mpipi-GG on its subset), at 300 K:
  total     calibrated total heating response (total_response_evidence.load)
  chi       fixed-water susceptibility
  nu        compaction
  L, L_adj  leftover of the total beyond compaction and length, and beyond pi and D+E+K as well (lineage_test.py on
            Analysis_FDT/total_response/lineage_inputs; SD units, negative = more contraction)
Organism means (each organism one vote); slope against OGT with kingdom fixed effects, organism bootstrap CI; Spearman
within kingdoms; OGT classes (<= 20, 21-30, 31-39, >= 40 C) within bacteria; and, because thermophily clusters in
particular phyla, the bacterial slope with phylum fixed effects (phylum from the OGT table; phyla with >= 3 organisms).
    python ogt_response.py [--out DIR]
Output: Thermal adaptation/bender_ogt/ (chains.csv, organisms.csv, slopes.csv, classes.csv, summary.txt)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import total_response_evidence as E  # noqa: E402

TA = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Thermal adaptation')
OUT = f'{TA}/bender_ogt'
LIN = f'{E.TLA}/Analysis_FDT/total_response/lineage'
N_BOOT = 2000


def ogt_table(bio, col='temperature'):
    t = pd.read_csv(f'{TA}/inputs/ogt/temperature_data.tsv', sep='\t')
    agg = 'mean' if col == 'temperature' else 'first'
    by_tax, by_name = t.groupby('taxid')[col].agg(agg), t.groupby('organism')[col].agg(agg)
    sp = bio.organism.str.split(' \\(').str[0].str.split().str[:2].str.join('_').str.lower()
    by_sp = sp.map(by_name)
    if 'kingdom' in bio.columns:    # fixed 2026-10-10: viruses must not inherit their host's OGT through the genus/species name
        by_sp = by_sp.where(bio.kingdom != 'Viruses')
    return bio.taxonomy_id.map(by_tax).fillna(by_sp)


NU_MIN = 0.30   # added 2026-10-20: nu at or below 0.30 marks a failed scaling fit, left out of nu means (as ogt_followup.py)


def slope(o, col, rng, fe='kingdom'):
    """OGT slope per 10 C with fixed effects (kingdom or phylum) on organism means; bootstrap over organisms."""
    o = o.dropna(subset=[col])
    def fit(d):
        X = np.column_stack([np.ones(len(d)), d.OGT, pd.get_dummies(d[fe], drop_first=True).astype(float)])
        return 10 * np.linalg.lstsq(X, d[col], rcond=None)[0][1]
    est = fit(o)
    boots = [fit(o.iloc[rng.integers(0, len(o), len(o))]) for _ in range(N_BOOT)]
    w = o.assign(ro=o.OGT - o.groupby(fe).OGT.transform('mean'), rc=o[col] - o.groupby(fe)[col].transform('mean'))
    return est, np.percentile(boots, 2.5), np.percentile(boots, 97.5), stats.spearmanr(w.ro, w.rc)[0], stats.spearmanr(w.ro, w.rc)[1]


def main(out):
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(1)
    D, _ = E.load()
    bio = pd.read_csv(E.RC.BIO_CSV, usecols=['UniProt_ID', 'organism', 'taxonomy_id', 'kingdom'], low_memory=False).set_index('UniProt_ID')
    bio['OGT'] = ogt_table(bio)
    bio['phylum'] = ogt_table(bio, 'phylum')
    lin = {'C': pd.read_csv(f'{LIN}/per_protein.csv').set_index('UniProt_ID'),
           'R': pd.read_csv(f'{LIN}/recharged/per_protein.csv').set_index('UniProt_ID')}
    rows, cls_rows, chains = [], [], []
    for k in ('C', 'R', 'M'):
        d = pd.DataFrame({'total': D[k].total, 'chi': D[k].chi, 'nu': D[k].nu}).join(bio, how='inner')
        d.loc[d.nu <= NU_MIN, 'nu'] = np.nan     # failed scaling fits
        if k in lin:
            d = d.join(lin[k][['L', 'L_adj']])
        d = d.dropna(subset=['OGT'])
        d['model'] = E.MODELS[k]
        chains.append(d)
        cols = [c for c in ('total', 'chi', 'nu', 'L', 'L_adj') if c in d]
        o = d.groupby(['organism', 'kingdom']).agg(OGT=('OGT', 'first'), chains=('total', 'size'),
                                                   **{c: (c, 'mean') for c in cols}).reset_index()
        o = o[o.kingdom.isin(o.kingdom.value_counts()[lambda s: s >= 5].index)]   # kingdoms with >= 5 organisms
        for c in cols:
            est, lo, hi, rho, p = slope(o, c, rng)
            rows.append({'model': E.MODELS[k], 'quantity': c, 'fixed_effects': 'kingdom', 'chains': int(o.chains.sum()),
                         'organisms': len(o),
                         'slope_per_10C': est, 'ci_lo': lo, 'ci_hi': hi, 'spearman_within_group': rho, 'p_spearman': p,
                         'sd_across_organisms': o[c].std()})
        b = d[d.kingdom == 'Bacteria']
        op = b.dropna(subset=['phylum']).groupby(['organism', 'phylum']).agg(
            OGT=('OGT', 'first'), chains=('total', 'size'), **{c: (c, 'mean') for c in cols}).reset_index()
        op = op[op.phylum.isin(op.phylum.value_counts()[lambda s: s >= 3].index)]
        for c in cols:
            est, lo, hi, rho, p = slope(op, c, rng, fe='phylum')
            rows.append({'model': E.MODELS[k], 'quantity': c, 'fixed_effects': 'phylum (bacteria only)',
                         'chains': int(op.chains.sum()), 'organisms': len(op), 'slope_per_10C': est, 'ci_lo': lo,
                         'ci_hi': hi, 'spearman_within_group': rho, 'p_spearman': p, 'sd_across_organisms': op[c].std()})
        b = b.assign(cls=pd.cut(b.OGT, [-1, 20, 30, 39.9, 100], labels=['<= 20', '21-30', '31-39', '>= 40']))
        ob = b.groupby(['organism', 'cls'], observed=True)[cols].mean().reset_index()
        for cl, g in ob.groupby('cls', observed=True):
            cls_rows.append({'model': E.MODELS[k], 'OGT class (bacteria)': cl, 'organisms': len(g),
                             'chains': int((b.cls == cl).sum()), **{f'{c}_mean': g[c].mean() for c in cols}})
    S, Cl = pd.DataFrame(rows), pd.DataFrame(cls_rows)
    ch = pd.concat(chains)
    ch.to_csv(os.path.join(out, 'chains.csv'))
    S.to_csv(os.path.join(out, 'slopes.csv'), index=False)
    Cl.to_csv(os.path.join(out, 'classes.csv'), index=False)
    lines = [f'BENDER chains with an OGT: {ch[ch.model == "CALVADOS-2"].shape[0]} (CALVADOS-2), '
             f'{ch[ch.model == "Mpipi-Recharged"].shape[0]} (Mpipi-Recharged), {ch[ch.model == "Mpipi-GG"].shape[0]} (Mpipi-GG)',
             '', '== slope against OGT per 10 C, organism means, kingdom fixed effects (95% bootstrap CI over organisms)',
             S.round(4).to_string(index=False), '', '== bacteria by OGT class (organism means)', Cl.round(4).to_string(index=False)]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    main(ap.parse_args().out)
