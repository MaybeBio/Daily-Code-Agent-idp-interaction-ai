"""
heat_crosses_analysis.py — results of heat_crosses.py (CALVADOS-2, 405 chains: annotated-phosphosite proteins, the Arg
shortlist and 300 random chains).
  Q   checks: residue shares add up to the total minus the bond term; recomputed dV/dT matches the npz
  A   full-resolution dynamics: integrated autocorrelation time of Rg (ns, model time) vs response, nu and log N fixed
  B   thermal hotspots: share of the solvation and total response by residue type (per residue), how concentrated a
      chain's response is (share carried by its top 10% of residues), and position along the chain
  P   phosphosites: within each protein with annotated sites, do annotated Ser/Thr carry more of the heating response
      (more negative total share), or change it more when phosphorylated (|d total|), than the protein's other Ser/Thr?
      Mean within-protein percentile, label-shuffling permutation within proteins (as phosphosite_test.py)
  C   heat x phosphorylation and heat x salt: size of the change in the total response, share beyond 2 SE with adequate
      effective sample (> 0.1), relative to the chain's total, and what sets the sign (net charge)
    python heat_crosses_analysis.py [--dir DIR]
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
DIR = f'{TLA}/Analysis_FDT/total_response/heat_crosses'
BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
N_PERM = 5000


def rn(v):
    return stats.norm.ppf((stats.rankdata(v) - 0.5) / len(v))


def partial(y, x, Z):
    Zr = np.column_stack([np.ones(len(y))] + [rn(z) for z in Z])
    ry = rn(y) - Zr @ np.linalg.lstsq(Zr, rn(y), rcond=None)[0]
    rx = rn(x) - Zr @ np.linalg.lstsq(Zr, rn(x), rcond=None)[0]
    return np.corrcoef(ry, rx)[0, 1]


def within_protein_test(d, score, rng):
    """Mean within-protein percentile of annotated sites by `score` (higher = more), and a permutation p."""
    d = d.copy()
    d['pct'] = d.groupby('UniProt_ID')[score].rank(pct=True)
    obs = d.loc[d.annotated, 'pct'].mean()
    groups = [g.pct.to_numpy() for _, g in d.groupby('UniProt_ID')]
    nann = [int(g.annotated.sum()) for _, g in d.groupby('UniProt_ID')]
    null = np.empty(N_PERM)
    for b in range(N_PERM):
        null[b] = np.mean(np.concatenate([rng.choice(g, k, replace=False) for g, k in zip(groups, nann) if k]))
    return obs, null.mean(), float((np.abs(null - null.mean()) >= abs(obs - null.mean())).mean())


def main(d):
    rng = np.random.default_rng(6)
    C = pd.read_csv(f'{d}/chains.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    R = pd.read_csv(f'{d}/residues.csv').drop_duplicates(['UniProt_ID', 'pos'])   # appended files: guard against resumed runs
    S = pd.read_csv(f'{d}/sites.csv').drop_duplicates(['UniProt_ID', 'pos'])
    smp = pd.read_csv(f'{d}/sample.csv').set_index('UniProt_ID')
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'nu', 'NCPR', 'FCR', 'length'], low_memory=False).set_index('UniProt_ID')
    C = C.join(bio).join(smp)
    lines = []
    # Q
    lines += ['== Q checks', f"max |residue sum - (total - bond)| = {C.residue_sum_minus_total.abs().max():.2e} % per K; "
              f"max |solvation residue sum - solvation| = {C.solv_residue_sum_check.abs().max():.2e}; "
              f"min corr(dV/dT recomputed, npz) = {C.dV_dT_check.min():.6f}"]
    # A
    a = C.dropna(subset=['tau_int_ns', 'nu'])
    lines += ['', '== A full-resolution dynamics (DCD spacing)',
              f"tau_int median {a.tau_int_frames_full.median():.2f} DCD frames ({a.tau_int_ns.median() * 1000:.0f} ps model time), "
              f"5th-95th pct {a.tau_int_ns.quantile(0.05) * 1000:.0f}-{a.tau_int_ns.quantile(0.95) * 1000:.0f} ps",
              f"tau vs total response: Spearman {stats.spearmanr(a.tau_int_ns, a.total_rel)[0]:.2f}, partial | nu, log N "
              f"{partial(a.total_rel.to_numpy(), a.tau_int_ns.to_numpy(), [a.nu.to_numpy(), np.log(a.N).to_numpy()]):.2f}; "
              f"vs chi: {stats.spearmanr(a.tau_int_ns, a.chi_rel)[0]:.2f}, partial "
              f"{partial(a.chi_rel.to_numpy(), a.tau_int_ns.to_numpy(), [a.nu.to_numpy(), np.log(a.N).to_numpy()]):.2f}"]
    # B
    Rr = R[R.UniProt_ID.isin(smp.index[smp.random])].copy()
    by = Rr.groupby('residue').agg(n=('pos', 'size'), solv_share_median=('solv_share', 'median'),
                                   boltz_share_median=('boltz_share', 'median'), diel_share_median=('diel_share', 'median'),
                                   total_share_median=('total_share', 'median')).sort_values('total_share_median')
    def top10(g):
        t = g.total_share.to_numpy()
        k = max(1, int(round(0.1 * len(t))))
        return np.sort(t)[:k].sum() / t.sum() if t.sum() != 0 else np.nan
    conc = Rr.groupby('UniProt_ID').apply(top10, include_groups=False)
    Rr['rel_pos'] = Rr.pos / Rr.groupby('UniProt_ID').pos.transform('max')
    prof = Rr.groupby(pd.cut(Rr.rel_pos, np.linspace(0, 1, 11)), observed=True).total_share.median()
    by.to_csv(f'{d}/analysis_residue_types.csv')
    lines += ['', '== B thermal hotspots (300 random chains; shares in % of <Rg> per K, per residue)', by.round(4).to_string(), '',
              f"share of a chain's total carried by its 10% most-contracting residues: median {conc.median():.2f} "
              f"(IQR {conc.quantile(0.25):.2f}-{conc.quantile(0.75):.2f}; 0.10 if spread evenly)",
              'median total share by position along the chain (tenths): ' + ', '.join(f'{v:.4f}' for v in prof.values)]
    # P
    ann = pd.read_csv(f'{TLA}/Analysis_FDT/phosphosite_test/sites.csv')
    ann = ann[ann.model == 'CALVADOS-2'][['UniProt_ID', 'pos', 'annotated']]
    st = R[R.residue.isin(['S', 'T'])].merge(S, on=['UniProt_ID', 'pos', 'residue']).merge(ann, on=['UniProt_ID', 'pos'], how='inner')
    st = st[st.groupby('UniProt_ID').annotated.transform('any')]
    st['neg_total_share'] = -st.total_share
    st['neg_solv_share'] = -st.solv_share
    st['abs_phospho_dtotal'] = st.phospho_dtotal.abs()
    rows = []
    for score, lab in (('neg_total_share', 'contributes more contraction (total share)'),
                       ('neg_solv_share', 'contributes more contraction (solvation share)'),
                       ('abs_phospho_dtotal', 'phosphorylation changes the heating response more')):
        obs, null, p = within_protein_test(st, score, rng)
        rows.append({'test': lab, 'proteins': st.UniProt_ID.nunique(), 'annotated_sites': int(st.annotated.sum()),
                     'ST_sites': len(st), 'mean_percentile_annotated': obs, 'null': null, 'p_perm': p})
    Pt = pd.DataFrame(rows)
    Pt.to_csv(f'{d}/analysis_phosphosites.csv', index=False)
    lines += ['', '== P phosphosites (within-protein percentile; null from shuffling labels within proteins)', Pt.round(4).to_string(index=False)]
    # C
    s2 = S.merge(C[['total_rel', 'NCPR', 'FCR']], left_on='UniProt_ID', right_index=True)
    ok = s2.ess > 0.1
    sig = ok & (s2.phospho_dtotal.abs() > 2 * s2.phospho_dtotal_se)
    c2 = C.dropna(subset=['salt_dtotal'])
    okc = c2.salt_ess > 0.1
    sigc = okc & (c2.salt_dtotal.abs() > 2 * c2.salt_dtotal_se)
    lines += ['', '== C heat x phosphorylation (per Ser/Thr; change in the total response, % per K)',
              f"{len(s2)} sites in {s2.UniProt_ID.nunique()} chains; effective sample > 0.1 at {ok.mean():.0%}",
              f"|change| median {s2.loc[ok, 'phospho_dtotal'].abs().median():.4f} (5-95%: {s2.loc[ok, 'phospho_dtotal'].quantile(0.05):.3f} to "
              f"{s2.loc[ok, 'phospho_dtotal'].quantile(0.95):.3f}); relative to the chain total: median "
              f"{(s2.loc[ok, 'phospho_dtotal'].abs() / s2.loc[ok, 'total_rel'].abs()).median():.1%}; beyond 2 SE at {sig[ok].mean():.0%} of usable sites",
              f"sign: Spearman of the change with chain NCPR {stats.spearmanr(s2.loc[ok, 'NCPR'], s2.loc[ok, 'phospho_dtotal'])[0]:.2f}, "
              f"with the size change on phosphorylation {stats.spearmanr(s2.loc[ok, 'phospho_dRg_pct'], s2.loc[ok, 'phospho_dtotal'])[0]:.2f}",
              '', '== C heat x salt (doubling the ionic strength; change in the total response, % per K)',
              f"{len(c2)} chains; effective sample > 0.1 at {okc.mean():.0%}; |change| median {c2.loc[okc, 'salt_dtotal'].abs().median():.4f} "
              f"(5-95%: {c2.loc[okc, 'salt_dtotal'].quantile(0.05):.3f} to {c2.loc[okc, 'salt_dtotal'].quantile(0.95):.3f}); relative to the "
              f"chain total: median {(c2.loc[okc, 'salt_dtotal'].abs() / c2.loc[okc, 'total_rel'].abs()).median():.1%}; beyond 2 SE at {sigc[okc].mean():.0%} of usable chains",
              f"Spearman of the change with |NCPR| {stats.spearmanr(c2.loc[okc, 'NCPR'].abs(), c2.loc[okc, 'salt_dtotal'])[0]:.2f}, "
              f"with FCR {stats.spearmanr(c2.loc[okc, 'FCR'], c2.loc[okc, 'salt_dtotal'])[0]:.2f}, with the size change on salt "
              f"{stats.spearmanr(c2.loc[okc, 'salt_dRg_pct'], c2.loc[okc, 'salt_dtotal'])[0]:.2f}"]
    open(f'{d}/analysis_summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=DIR)
    main(ap.parse_args().dir)
