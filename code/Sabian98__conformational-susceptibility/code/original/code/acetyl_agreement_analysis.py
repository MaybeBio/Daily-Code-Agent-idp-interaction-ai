"""
acetyl_agreement_analysis.py — why CALVADOS-2 and Mpipi-Recharged agree on the change of Rg on acetylating the same Lys
but less on the change of the total response (exploratory, after the fact), with the same noise check for single
substitutions (mutation_scan).

  noise ceiling   reliability of each model's changes, 1 - mean(SE^2)/var(change) (10-block jackknife SEs), and the largest
                  Pearson correlation two noise-limited but otherwise identical models could reach, sqrt(rel_C rel_R)
  resolved sites  rank correlation and sign agreement among sites beyond 1, 2 and 3 SE in both models
  parts           acetyl_agreement.py split: ensemble shift, Boltzmann part and electrostatic temperature part of the change;
                  agreement of each between models, its share of the variation (cov with the change / var of the change)
                  and its rank correlation with the change of Rg
Usable: Kish fraction > 0.1 in both models.
    python acetyl_agreement_analysis.py
Output: heat_acetyl/agreement/analysis_summary.txt
"""
import numpy as np
import pandas as pd
from scipy import stats

TR = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/total_response'
AC, MS = f'{TR}/heat_acetyl', f'{TR}/mutation_scan'
sp = lambda a, b: stats.spearmanr(a, b)[0]  # noqa: E731


def noise(M, d, se, lines, label):
    pc = np.corrcoef(M[f'{d}_C'], M[f'{d}_R'])[0, 1]
    rel = [1 - (M[f'{se}_{s}'] ** 2).mean() / M[f'{d}_{s}'].var() for s in 'CR']
    ceil = np.sqrt(rel[0] * rel[1])
    lines.append(f'{label}: {len(M)} usable in both; rank correlation of the change {sp(M[d + "_C"], M[d + "_R"]):.3f}, Pearson {pc:.3f}; '
                 f'reliability {rel[0]:.2f} (C2), {rel[1]:.2f} (MR); noise ceiling {ceil:.3f}; disattenuated {pc / ceil:.2f}')
    for s in 'CR':
        lines.append(f'  {s}: median |change| {M[d + "_" + s].abs().median():.4f}, median SE {M[se + "_" + s].median():.4f}, '
                     f'beyond 2 SE {(M[d + "_" + s].abs() > 2 * M[se + "_" + s]).mean():.0%}')
    for t in (1, 2, 3):
        r = M[(M[f'{d}_C'].abs() > t * M[f'{se}_C']) & (M[f'{d}_R'].abs() > t * M[f'{se}_R'])]
        lines.append(f'  both beyond {t} SE: {len(r)}, rank correlation {sp(r[d + "_C"], r[d + "_R"]):.2f}, '
                     f'same sign {(np.sign(r[d + "_C"]) == np.sign(r[d + "_R"])).mean():.0%}')


def main():
    lines = []
    k = ['UniProt_ID', 'pos']
    C, R = (pd.read_csv(f'{AC}/sites_{m}.csv').drop_duplicates(k) for m in 'CR')
    M = C.merge(R, on=k, suffixes=('_C', '_R'))
    M = M[(M.ess_C > 0.1) & (M.ess_R > 0.1)]
    lines.append(f'Acetylation: rank correlation of the change of Rg {sp(M.acetyl_dRg_pct_C, M.acetyl_dRg_pct_R):.3f}')
    noise(M, 'acetyl_dtotal', 'acetyl_dtotal_se', lines, 'Acetylation, change of the total response')
    P = {m: pd.read_csv(f'{AC}/agreement/parts_{m}.csv').drop_duplicates(k) for m in 'CR'}
    for m in 'CR':
        o = pd.read_csv(f'{AC}/sites_{m}.csv').drop_duplicates(k).merge(P[m], on=k)
        lines.append(f'  parts {m}: sum reproduces the stored change to {np.abs(o.dtotal - o.acetyl_dtotal).max():.1e}')
    Q = P['C'].merge(P['R'], on=k, suffixes=('_C', '_R'))
    Q = Q[(Q.ess_C > 0.1) & (Q.ess_R > 0.1)]
    for c in ('shift', 'boltzmann', 'electro'):
        share = {s: np.cov(Q[f'{c}_{s}'], Q[f'dtotal_{s}'])[0, 1] / Q[f'dtotal_{s}'].var() for s in 'CR'}
        lines.append(f'  {c}: agreement C2-MR {sp(Q[c + "_C"], Q[c + "_R"]):.2f}; share of variation {share["C"]:.2f} (C2), {share["R"]:.2f} (MR); '
                     f'rank correlation with dRg {sp(Q[c + "_C"], Q.dRg_pct_C):.2f} (C2), {sp(Q[c + "_R"], Q.dRg_pct_R):.2f} (MR); '
                     f'median |.| {Q[c + "_C"].abs().median():.4f}, {Q[c + "_R"].abs().median():.4f}')
    km = ['UniProt_ID', 'pos', 'mut']
    S = pd.read_csv(f'{MS}/mutants.csv').drop_duplicates(km).merge(
        pd.read_csv(f'{MS}/mutants_recharged.csv').drop_duplicates(km), on=km + ['wt'], suffixes=('_C', '_R'))
    S = S[(S.ess_C > 0.1) & (S.ess_R > 0.1)]
    lines.append('')
    lines.append(f'Single substitutions: rank correlation of the change of Rg {sp(S.dRg_pct_C, S.dRg_pct_R):.3f}')
    noise(S, 'dtotal', 'dtotal_se', lines, 'Single substitutions, change of the total response')
    open(f'{AC}/agreement/analysis_summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
