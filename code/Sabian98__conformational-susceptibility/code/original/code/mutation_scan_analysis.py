"""
mutation_scan_analysis.py — results of mutation_scan.py (CALVADOS-2, 50 chains, every substitution at every position).
Only mutants with Kish effective sample > 0.1 are used (reweighting reliable); 'clear' = |change| > 2 SE.
Sign: the total response is negative (contraction); a negative change = stronger contraction on heating.
  1  reliability: share of mutants usable, by residue introduced
  2  what a substitution does: median change by residue introduced and by residue replaced; the compaction rule for
     mutations (Spearman of the change with the size change)
  3  Arg -> Lys at every Arg (Arg shortlist vs random chains)
  4  per chain: the largest single-mutation effect relative to the chain's total; share of positions with a clear effect
  5  hotspots: per position, mean |change| over substitutions vs that residue's share of the total (heat_crosses.py)
  6  design table: per chain, the strengthening and weakening single mutations with the largest |change| among those
     beyond 2 SE with effective sample > 0.5 ('clearest' below means largest clear effect, not largest z)
    python mutation_scan_analysis.py [--dir DIR] [--hx DIR]
Output: DIR/analysis_*.csv, analysis_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

TLA = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
DIR = f'{TLA}/Analysis_FDT/total_response/mutation_scan'
HX = f'{TLA}/Analysis_FDT/total_response/heat_crosses'


def main(d):
    M = pd.read_csv(f'{d}/mutants.csv').drop_duplicates(['UniProt_ID', 'pos', 'mut'])
    C = pd.read_csv(f'{d}/chains.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    smp = pd.read_csv(f'{d}/sample.csv').set_index('UniProt_ID')
    M = M.join(C[['total_rel', 'N']], on='UniProt_ID').join(smp, on='UniProt_ID')
    M['ok'] = M.ess > 0.1
    M['clear'] = M.ok & (M.dtotal.abs() > 2 * M.dtotal_se)
    lines = [f'{len(M):,} mutants in {M.UniProt_ID.nunique()} chains; usable (effective sample > 0.1): {M.ok.mean():.1%}; '
             f'clear (> 2 SE) among usable: {M.loc[M.ok, "clear"].mean():.1%}']
    U = M[M.ok]
    rel = U.groupby('mut').ok.size() / M.groupby('mut').size()
    lines += ['usable share by residue introduced: ' + ', '.join(f'{a} {v:.0%}' for a, v in rel.sort_values().items())]
    # 2
    to = U.groupby('mut').agg(n=('dtotal', 'size'), median_change=('dtotal', 'median'), frac_stronger=('dtotal', lambda s: (s < 0).mean()),
                              median_dRg=('dRg_pct', 'median')).sort_values('median_change')
    fr = U.groupby('wt').agg(n=('dtotal', 'size'), median_change=('dtotal', 'median')).sort_values('median_change')
    mat = U.pivot_table(index='wt', columns='mut', values='dtotal', aggfunc='median')
    to.to_csv(f'{d}/analysis_by_mut.csv')
    fr.to_csv(f'{d}/analysis_by_wt.csv')
    mat.to_csv(f'{d}/analysis_matrix.csv')
    lines += ['', '== 2 by residue introduced (median change, % per K; negative = stronger contraction on heating)',
              to.round(4).to_string(), '', 'by residue replaced:', fr.round(4).to_string(),
              f'\ncompaction rule for mutations: Spearman(change, size change) = {stats.spearmanr(U.dtotal, U.dRg_pct)[0]:.2f} '
              f'(n = {len(U):,})']
    # 3
    rk = U[(U.wt == 'R') & (U.mut == 'K')]
    g = rk.groupby('arg_shortlist').agg(sites=('dtotal', 'size'), median_change=('dtotal', 'median'),
                                       frac_weaker=('dtotal', lambda s: (s > 0).mean()), clear=('clear', 'mean'))
    per_chain = rk.groupby('UniProt_ID').agg(n_R=('dtotal', 'size'), sum_change=('dtotal', 'sum'), total=('total_rel', 'first'))
    per_chain['sum_rel'] = per_chain.sum_change / per_chain.total.abs()
    lines += ['', '== 3 Arg -> Lys (single sites)', g.round(4).to_string(),
              f'summed over all Arg of a chain (linear sum of single-site effects): median {per_chain.sum_change.median():+.3f} % per K '
              f'({per_chain.sum_rel.median():+.0%} of the chain total)']
    # 4
    def chain_stats(x):
        pos = x.groupby('pos').apply(lambda y: y.clear.any(), include_groups=False)
        return pd.Series({'N': x.N.iloc[0], 'total': x.total_rel.iloc[0], 'max_abs_change': x.dtotal.abs().max(),
                          'max_rel': x.dtotal.abs().max() / abs(x.total_rel.iloc[0]), 'frac_positions_clear': pos.mean()})
    CS = U.groupby('UniProt_ID').apply(chain_stats, include_groups=False).join(smp)
    CS.to_csv(f'{d}/analysis_chains.csv')
    lines += ['', '== 4 per chain', CS.groupby('arg_shortlist')[['max_abs_change', 'max_rel', 'frac_positions_clear']].median().round(3).to_string()]
    # 5
    if os.path.exists(f'{HX}/residues.csv'):
        R = pd.read_csv(f'{HX}/residues.csv')
        P = U.groupby(['UniProt_ID', 'pos']).dtotal.apply(lambda s: s.abs().mean()).rename('sensitivity').reset_index()
        P = P.merge(R[['UniProt_ID', 'pos', 'total_share', 'solv_share']], on=['UniProt_ID', 'pos'])
        rho = P.groupby('UniProt_ID').apply(lambda x: stats.spearmanr(-x.total_share, x.sensitivity)[0], include_groups=False)
        lines += ['', f'== 5 hotspots vs mutational sensitivity ({P.UniProt_ID.nunique()} chains): within-chain Spearman of '
                  f'(contraction share) with (mean |change| over substitutions), median {rho.median():.2f} (IQR {rho.quantile(0.25):.2f}-{rho.quantile(0.75):.2f})']
    # 6
    good = U[(U.ess > 0.5) & U.clear]
    des = []
    for uid, x in good.groupby('UniProt_ID'):
        for lab, row in (('strengthens contraction', x.loc[x.dtotal.idxmin()]), ('weakens contraction', x.loc[x.dtotal.idxmax()])):
            des.append({'UniProt_ID': uid, 'arg_shortlist': bool(row.arg_shortlist), 'effect': lab,
                        'mutation': f'{row.wt}{int(row.pos)}{row.mut}', 'change': row.dtotal, 'se': row.dtotal_se,
                        'chain_total': row.total_rel, 'rel_change': row.dtotal / abs(row.total_rel), 'size_change_pct': row.dRg_pct,
                        'ess': row.ess})
    De = pd.DataFrame(des)
    De.to_csv(f'{d}/analysis_design.csv', index=False)
    lines += ['', '== 6 design: clearest single mutations per chain (effective sample > 0.5, > 2 SE)',
              f'median relative change: strengthening {De.loc[De.effect.str.startswith("strengthens"), "rel_change"].median():+.0%}, '
              f'weakening {De.loc[De.effect.str.startswith("weakens"), "rel_change"].median():+.0%}',
              'most common residues introduced: strengthening ' +
              De[De.effect.str.startswith('strengthens')].mutation.str[-1].value_counts().head(5).to_dict().__repr__() +
              '; weakening ' + De[De.effect.str.startswith('weakens')].mutation.str[-1].value_counts().head(5).to_dict().__repr__()]
    open(f'{d}/analysis_summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=DIR)
    ap.add_argument('--hx', default=HX, help='heat_crosses output with residues.csv (heat_crosses_A0_linear for the linear A0 runs)')
    a = ap.parse_args()
    HX = a.hx
    main(a.dir)
