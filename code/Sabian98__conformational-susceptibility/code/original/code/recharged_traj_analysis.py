"""
recharged_traj_analysis.py — results of recharged_traj_pass.py, each next to its CALVADOS-2 / Mpipi-GG counterpart.
  P  pair types (slide 23): median share of the sticker susceptibility and per-pair strength by type; within-protein
     R-aro / K-aro and R-other / K-other per-pair ratios (proteins with >= 20 pairs of both); vs pair_type_split.py
  H  thermal hotspots: per-residue shares by residue type (random chains of heat_crosses' sample); top-10% concentration;
     within-chain agreement of per-residue total shares with CALVADOS-2 (heat_crosses residues)
  S  heat x salt: change of the total vs the size change and |NCPR|; agreement per chain with CALVADOS-2 and Mpipi-GG
  F  phosphomimetics (Ser/Thr -> Glu, -> Asp): the compaction rule; agreement of the size change (leverage) and of the
     heat-response change with CALVADOS-2's charge-model phosphorylation (heat_crosses sites); the phosphosite tests
     (within-protein percentile of annotated sites, label permutation) with Recharged's Glu-mimic leverage and its
     heat-response change
  B  basins (common set): share multi-basin; between-basin share of the plain and total response
    python recharged_traj_analysis.py
Output: Analysis_FDT_RECHARGED/traj_pass/analysis_*.csv, analysis_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

TLA = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
D = f'{TLA}/Analysis_FDT_RECHARGED/traj_pass'
HX = f'{TLA}/Analysis_FDT/total_response/heat_crosses'
SH = f'{TLA}/Analysis_FDT/total_response/heat_crosses_shared'
PT = f'{TLA}/Analysis_FDT_MPIPI/compact/pair_type_split/per_protein.csv'
BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
TYPES = ['aro-aro', 'R-aro', 'K-aro', 'R-other', 'K-other', 'aro-other', 'other-other']
N_PERM = 5000


def within_protein_test(d, score, rng):
    d = d.copy()
    d['pct'] = d.groupby('UniProt_ID')[score].rank(pct=True)
    obs = d.loc[d.annotated, 'pct'].mean()
    groups = [(g.pct.to_numpy(), int(g.annotated.sum())) for _, g in d.groupby('UniProt_ID')]
    null = np.array([np.mean(np.concatenate([rng.choice(g, k, replace=False) for g, k in groups if k])) for _ in range(N_PERM)])
    return obs, null.mean(), float((np.abs(null - null.mean()) >= abs(obs - null.mean())).mean())


def main():
    rng = np.random.default_rng(7)
    C = pd.read_csv(f'{D}/chains.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    P = pd.read_csv(f'{D}/pairs.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    R = pd.read_csv(f'{D}/residues.csv').drop_duplicates(['UniProt_ID', 'pos'])          # appended files (resumable runs)
    F = pd.read_csv(f'{D}/phosphomimetic.csv').drop_duplicates(['UniProt_ID', 'pos', 'mimic'])
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'NCPR'], low_memory=False).set_index('UniProt_ID')
    lines = [f'{len(C)} chains; {len(R):,} residues; {len(F):,} phosphomimetic mutants']
    # P
    pt = pd.read_csv(PT)
    rows = []
    for model, src in (('Mpipi-Recharged', P.loc[P.index.isin(pt.UniProt_ID)].reset_index()), ('CALVADOS-2', pt[pt.model == 'CALVADOS-2']),
                       ('Mpipi-GG', pt[pt.model == 'Mpipi-GG'])):
        r = {'model': model, 'proteins': len(src)}
        for t in TYPES:
            r[f'share_chi {t}'] = src[f'share_chi_{t}'].median()
            r[f'per_pair {t}'] = src.loc[src[f'n_pairs_{t}'] >= 20, f'chi_per_pair_{t}'].median() if f'n_pairs_{t}' in src else np.nan
        for x in ('aro', 'other'):
            ok = (src[f'n_pairs_R-{x}'] >= 20) & (src[f'n_pairs_K-{x}'] >= 20) if f'n_pairs_R-{x}' in src else None
            if ok is not None:
                r[f'R/K {x} per pair'] = (src.loc[ok, f'chi_per_pair_R-{x}'] / src.loc[ok, f'chi_per_pair_K-{x}']).median()
        rows.append(r)
    Pt = pd.DataFrame(rows)
    Pt.to_csv(f'{D}/analysis_pairs.csv', index=False)
    lines += ['', '== P pair types (medians)', Pt.round(3).T.to_string()]
    # H
    smp = pd.read_csv(f'{HX}/sample.csv')
    rnd = set(smp.loc[smp.random, 'UniProt_ID'])
    Rr = R[R.UniProt_ID.isin(rnd)]
    by = Rr.groupby('residue')[['solv_share', 'boltz_share', 'screen_share', 'total_share']].median().sort_values('total_share')
    Rc = pd.read_csv(f'{HX}/residues.csv')
    byc = Rc[Rc.UniProt_ID.isin(rnd)].groupby('residue').total_share.median()
    def top10(g):
        t = np.sort(g.total_share.to_numpy())
        k = max(1, int(round(0.1 * len(t))))
        return t[:k].sum() / t.sum() if t.sum() else np.nan
    conc = Rr.groupby('UniProt_ID').apply(top10, include_groups=False)
    m = Rr.merge(Rc, on=['UniProt_ID', 'pos', 'residue'], suffixes=('_R', '_C'))
    rho = m.groupby('UniProt_ID').apply(lambda x: stats.spearmanr(x.total_share_R, x.total_share_C)[0], include_groups=False)
    by.to_csv(f'{D}/analysis_residue_types.csv')
    lines += ['', '== H thermal hotspots (Mpipi-Recharged, random chains; % per K per residue)', by.round(4).to_string(),
              f'residue-type ranking vs CALVADOS-2: Spearman {stats.spearmanr(by.total_share, byc.reindex(by.index))[0]:.2f}',
              f'top 10% of residues carry {conc.median():.2f} of the total (CALVADOS-2: 0.23)',
              f'within-chain agreement of per-residue shares with CALVADOS-2: median Spearman {rho.median():.2f} (IQR {rho.quantile(0.25):.2f}-{rho.quantile(0.75):.2f})']
    # S
    Cs = C[C.salt_ess > 0.1].join(bio)
    hc = pd.read_csv(f'{HX}/chains.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    sg = pd.read_csv(f'{SH}/chains_gg.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    jc = Cs.join(hc[['salt_dtotal', 'salt_ess']], rsuffix='_C', how='inner')
    jc = jc[jc.salt_ess_C > 0.1]
    jg = Cs.join(sg[['salt_dtotal', 'salt_ess']], rsuffix='_M', how='inner')
    jg = jg[jg.salt_ess_M > 0.1]
    lines += ['', '== S heat x salt (Mpipi-Recharged)',
              f'{len(Cs)} chains: median |change| {Cs.salt_dtotal.abs().median():.4f} % per K; beyond 2 SE {(Cs.salt_dtotal.abs() > 2 * Cs.salt_dtotal_se).mean():.0%}; '
              f'Spearman with the size change {stats.spearmanr(Cs.salt_dtotal, Cs.salt_dRg_pct)[0]:.2f}, with |NCPR| {stats.spearmanr(Cs.salt_dtotal, Cs.NCPR.abs())[0]:.2f}',
              f'per chain vs CALVADOS-2 ({len(jc)}): Spearman {stats.spearmanr(jc.salt_dtotal, jc.salt_dtotal_C)[0]:.2f}; vs Mpipi-GG ({len(jg)}): {stats.spearmanr(jg.salt_dtotal, jg.salt_dtotal_M)[0]:.2f}']
    # F
    Fo = F[F.ess > 0.1]
    lines += ['', '== F phosphomimetics (Mpipi-Recharged)',
              f"usable {(F.ess > 0.1).mean():.0%}; " + '; '.join(
                  f"{b}: median size change {g.dRg_pct.median():+.2f}%, |change of total| {g.dtotal.abs().median():.4f}, beyond 2 SE {(g.dtotal.abs() > 2 * g.dtotal_se).mean():.0%}, "
                  f"Spearman(change, size change) {stats.spearmanr(g.dtotal, g.dRg_pct)[0]:.2f}" for b, g in Fo.groupby('mimic'))]
    sc = pd.concat([pd.read_csv(f'{HX}/sites.csv'), pd.read_csv(f'{SH}/sites.csv')]).drop_duplicates(['UniProt_ID', 'pos'])
    e = Fo[Fo.mimic == 'E'].merge(sc, on=['UniProt_ID', 'pos', 'residue'], suffixes=('_R', '_C'))
    e = e[e.ess_C > 0.1]
    lines += [f'Glu mimic vs CALVADOS-2 charge-model phosphate ({len(e):,} sites): size change Spearman {stats.spearmanr(e.dRg_pct, e.phospho_dRg_pct)[0]:.2f}, '
              f'heat-response change Spearman {stats.spearmanr(e.dtotal, e.phospho_dtotal)[0]:.2f}']
    ann = pd.read_csv(f'{TLA}/Analysis_FDT/phosphosite_test/sites.csv')
    ann = ann[ann.model == 'CALVADOS-2'][['UniProt_ID', 'pos', 'annotated']]
    pe = Fo[Fo.mimic == 'E'].merge(ann, on=['UniProt_ID', 'pos'])
    pe = pe[pe.groupby('UniProt_ID').annotated.transform('any')]
    pe['abs_dRg'], pe['abs_dtotal'] = pe.dRg_pct.abs(), pe.dtotal.abs()
    rows = []
    for score, lab in (('abs_dRg', 'leverage on size (Glu mimic)'), ('abs_dtotal', 'change of the heating response (Glu mimic)')):
        obs, null, p = within_protein_test(pe, score, rng)
        rows.append({'test': lab, 'proteins': pe.UniProt_ID.nunique(), 'annotated': int(pe.annotated.sum()), 'sites': len(pe),
                     'mean_percentile': obs, 'null': null, 'p_perm': p})
    Pf = pd.DataFrame(rows)
    Pf.to_csv(f'{D}/analysis_phosphosites.csv', index=False)
    lines += ['phosphosite tests in Mpipi-Recharged:', Pf.round(4).to_string(index=False)]
    # B
    b = C.dropna(subset=['n_basins'])
    mb = b[b.multi_basin == True]  # noqa: E712
    lines += ['', '== B basins (common set, Mpipi-Recharged)',
              f'{len(b)} chains; multi-basin (both halves) {len(mb)} ({len(mb) / len(b):.1%}); any second basin {(b.n_basins > 1).mean():.1%}',
              f'between-basin share in multi-basin chains: plain {mb.between_share_plain.median():.2f}, total {mb.between_share_total.median():.2f}' if len(mb) else '']
    open(f'{D}/analysis_summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
