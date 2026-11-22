"""
make_compact_mpipi_fasta.py — FASTA of compact IDPs to run in Mpipi-GG at 300 K, to test the compact side of the
susceptibility-vs-nu trend (CALVADOS-2: a plateau for nu < 0.50, then a steep fall; the Mpipi-GG set has only 11
proteins below 0.45). The FDT needs one temperature, so 300 K runs are enough.

Selection (from BENDER_BIO.csv, CALVADOS-2 primary nu at 300 K):
  only well-defined nu: power-law fit R^2 >= 0.99 (below nu ~ 0.33 almost no fit reaches it: those chains are not
  power-law polymers, so their nu is not a compaction measure)
  nu < 0.45: every candidate · 0.45 <= nu < 0.50: the rest up to TOTAL, equal numbers per length tercile, fixed seed
  excluded: sequences already simulated in Mpipi-GG at 300 K, non-standard residues, N > 250 (run time ~(N/150)^2)

Output (--out, default Scripts/mpipi_sim/input/compact_mpipi_300K.fasta), header `>ID|kingdom|nu_cal=0.xxx|N=n`
(run_mpipi_from_fasta.py takes the ID before the first '|'), plus <out>.csv with the selection.

Usage: python3 make_compact_mpipi_fasta.py [--out PATH] [--seed 0]
"""
import argparse

import numpy as np
import pandas as pd

BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
MPIPI_300 = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT_MPIPI/mpipi_300K.csv'
OUT = '/Users/taseefrahman/Downloads/Scripts/mpipi_sim/input/compact_mpipi_300K.fasta'
TOTAL, R2_MIN, MAX_N = 250, 0.99, 250


def main(out, seed):
    b = pd.read_csv(BIO, usecols=['UniProt_ID', 'kingdom', 'length', 'nu', 'nu_fit_r2', 'sequence'], low_memory=False)
    done = set(pd.read_csv(MPIPI_300, usecols=['UniProt_ID'])['UniProt_ID'])
    b = b[~b['UniProt_ID'].isin(done) & b['sequence'].str.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+') & (b['length'] <= MAX_N) & (b['nu_fit_r2'] >= R2_MIN)]
    rng = np.random.default_rng(seed)
    compact = b[b['nu'] < 0.45].assign(nu_band='<0.45 (all)')
    picks = [compact]
    n = TOTAL - len(compact)
    band = b[(b['nu'] >= 0.45) & (b['nu'] < 0.50)].copy()
    band['N_tercile'] = pd.qcut(band['length'], 3, labels=['short', 'mid', 'long'])
    per = [n // 3 + (1 if k < n % 3 else 0) for k in range(3)]
    for (t, g), k in zip(band.groupby('N_tercile', observed=True), per):
        idx = rng.choice(len(g), size=min(k, len(g)), replace=False)
        picks.append(g.iloc[np.sort(idx)].assign(nu_band='0.45-0.50'))
    sel = pd.concat(picks).sort_values('nu').reset_index(drop=True)
    with open(out, 'w') as fh:
        for r in sel.itertuples():
            fh.write(f'>{r.UniProt_ID}|{r.kingdom}|nu_cal={r.nu:.3f}|N={r.length}\n{r.sequence}\n')
    sel.drop(columns='sequence').to_csv(out + '.csv', index=False)
    print(f'{len(sel)} sequences -> {out}')
    print(sel.groupby('nu_band').agg(n=('UniProt_ID', 'size'), median_N=('length', 'median'),
                                     max_N=('length', 'max')).to_string())
    print(sel['kingdom'].value_counts().to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    main(a.out, a.seed)
