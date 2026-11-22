"""
mutation_scan_compare.py — do CALVADOS-2 and Mpipi-Recharged agree on how single mutations change the heating response?
(mutation_scan.py, mutation_scan_recharged.py; same 50 chains, same mutants; usable = effective sample > 0.1 in both)
  1  mutant by mutant: Spearman of the change, same sign (all usable; and where both are clear, > 2 SE)
  2  rules: median change by residue introduced and by residue replaced in each model, and their rank agreement
  3  Arg -> Lys per site and summed per chain
  4  design picks: is each chain's clearest strengthening / weakening mutation in one model also in the same direction,
     and in the top 5%, in the other
    python mutation_scan_compare.py [--dir DIR]
Output: DIR/compare_*.csv, compare_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

DIR = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/total_response/mutation_scan'
KEY = ['UniProt_ID', 'pos', 'wt', 'mut']


def main(d):
    C = pd.read_csv(f'{d}/mutants.csv').drop_duplicates(['UniProt_ID', 'pos', 'mut'])
    R = pd.read_csv(f'{d}/mutants_recharged.csv').drop_duplicates(['UniProt_ID', 'pos', 'mut'])
    smp = pd.read_csv(f'{d}/sample.csv').set_index('UniProt_ID')
    J = C.merge(R, on=KEY, suffixes=('_C', '_R'))
    J = J[(J.ess_C > 0.1) & (J.ess_R > 0.1)].join(smp, on='UniProt_ID')
    J['clear_C'] = J.dtotal_C.abs() > 2 * J.dtotal_se_C
    J['clear_R'] = J.dtotal_R.abs() > 2 * J.dtotal_se_R
    both = J[J.clear_C & J.clear_R]
    lines = [f'{len(J):,} mutants usable in both models ({J.UniProt_ID.nunique()} chains)',
             f'mutant by mutant: Spearman {stats.spearmanr(J.dtotal_C, J.dtotal_R)[0]:.2f}, same sign {(np.sign(J.dtotal_C) == np.sign(J.dtotal_R)).mean():.0%}; '
             f'where both are clear ({len(both):,}): Spearman {stats.spearmanr(both.dtotal_C, both.dtotal_R)[0]:.2f}, '
             f'same sign {(np.sign(both.dtotal_C) == np.sign(both.dtotal_R)).mean():.0%}',
             f'size change: Spearman {stats.spearmanr(J.dRg_pct_C, J.dRg_pct_R)[0]:.2f}']
    rows = []
    for g in ('mut', 'wt'):
        t = J.groupby(g)[['dtotal_C', 'dtotal_R']].median()
        t.to_csv(f'{d}/compare_by_{g}.csv')
        rows.append((g, t, stats.spearmanr(t.dtotal_C, t.dtotal_R)[0]))
    for g, t, rho in rows:
        lab = 'residue introduced' if g == 'mut' else 'residue replaced'
        lines += ['', f'== by {lab} (median change, % per K): rank agreement of the 20 medians {rho:.2f}',
                  t.sort_values('dtotal_C').round(4).T.to_string()]
    rk = J[(J.wt == 'R') & (J.mut == 'K')]
    pc = rk.groupby('UniProt_ID')[['dtotal_C', 'dtotal_R']].sum()
    lines += ['', f'== Arg -> Lys: {len(rk)} sites; weaker contraction in CALVADOS-2 {(rk.dtotal_C > 0).mean():.0%}, '
                  f'Mpipi-Recharged {(rk.dtotal_R > 0).mean():.0%}; per-site Spearman {stats.spearmanr(rk.dtotal_C, rk.dtotal_R)[0]:.2f}; '
                  f'summed per chain, median {pc.dtotal_C.median():+.3f} / {pc.dtotal_R.median():+.3f} % per K']
    des = []
    for uid, x in J.groupby('UniProt_ID'):
        for m_a, m_b in (('C', 'R'), ('R', 'C')):
            good = x[(x[f'ess_{m_a}'] > 0.5) & x[f'clear_{m_a}']]
            if not len(good):
                continue
            for lab, row in (('strengthen', good.loc[good[f'dtotal_{m_a}'].idxmin()]), ('weaken', good.loc[good[f'dtotal_{m_a}'].idxmax()])):
                other = x[f'dtotal_{m_b}']
                v = row[f'dtotal_{m_b}']
                pct = (other < v).mean() if lab == 'strengthen' else (other > v).mean()
                des.append({'UniProt_ID': uid, 'picked_in': m_a, 'effect': lab, 'mutation': f'{row.wt}{int(row.pos)}{row.mut}',
                            'same_direction_in_other': (v < 0) if lab == 'strengthen' else (v > 0), 'top5pct_in_other': pct <= 0.05})
    De = pd.DataFrame(des)
    De.to_csv(f'{d}/compare_design.csv', index=False)
    lines += ['', '== design picks transferred to the other model',
              De.groupby(['picked_in', 'effect'])[['same_direction_in_other', 'top5pct_in_other']].mean().round(2).to_string()]
    open(f'{d}/compare_summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=DIR)
    main(ap.parse_args().dir)
