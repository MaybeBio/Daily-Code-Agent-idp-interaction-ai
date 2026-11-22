"""
acetyl_heat_analysis.py — results of acetyl_heat.py, set out as part C of heat_crosses_analysis.py for phosphorylation.

Per Lys (sites with Kish effective sample fraction > 0.1, the threshold used throughout):
  size of the change of the total response (% per K), relative to the chain's total, share beyond 2 SE
  compaction rule: Spearman of the change with the change of <Rg> on acetylation (positive = a site that compacts the
  chain on acetylation also strengthens its contraction on heating), and the share of sites whose two signs agree
  sign: Spearman with the chain's net charge per residue (acetylation removes a positive charge, so it should compact
  net-positive chains and expand net-negative ones)
  comparison with phosphorylation in the same chains: rank correlation of chain-median changes
Models: CALVADOS-2 and Mpipi-Recharged (the same 405 chains; Recharged has no phosphorylation run to compare with), and
Mpipi-GG (211 shared chains); site-by-site agreement of CALVADOS-2 and Mpipi-Recharged.
    python acetyl_heat_analysis.py
Output: heat_acetyl/analysis_summary.txt, analysis_by_model.csv (with a section on the fixed-solvent rerun in
        heat_acetyl/fixed_solvent/: compaction rule with and without the hydration forcings)
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

TR = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis/Analysis_FDT/total_response')
DIR = f'{TR}/heat_acetyl'
BIO = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/BENDER_BIO.csv')
# canonical folders hold the calibrated model with saturation since 2026-10-25 (A0 runs in *_A0_linear): the Mpipi-GG
# phosphosites of that model come from forcing_role.py --forcing capped, and the acetylation run carries the fixed-solvent
# change itself (no fixed_solvent subfolder)
SAT = os.path.isdir(f'{TR}/forcing_role_A0_linear')
PHOS = {'C': (f'{TR}/heat_crosses/sites.csv', 'CALVADOS-2'),
        'M': (f'{TR}/forcing_role/phospho_M_sites.csv' if SAT else f'{TR}/heat_crosses_shared/sites_gg.csv', 'Mpipi-GG'),
        'R': (None, 'Mpipi-Recharged')}
ESS_MIN = 0.1


def main():
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'nu', 'NCPR', 'FCR'], low_memory=False).set_index('UniProt_ID')
    lines, rows = [], []
    for k, (pfile, model) in PHOS.items():
        C = pd.read_csv(f'{DIR}/chains_{k}.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID').join(bio)
        S = pd.read_csv(f'{DIR}/sites_{k}.csv').drop_duplicates(['UniProt_ID', 'pos'])
        S = S.merge(C[['total_rel', 'NCPR']], left_on='UniProt_ID', right_index=True)
        ok = S.ess > ESS_MIN
        s = S[ok]
        sig = s.acetyl_dtotal.abs() > 2 * s.acetyl_dtotal_se
        rho_c, p_c = stats.spearmanr(s.acetyl_dRg_pct, s.acetyl_dtotal)
        rho_q, p_q = stats.spearmanr(s.NCPR, s.acetyl_dtotal)
        rho_qr, _ = stats.spearmanr(s.NCPR, s.acetyl_dRg_pct)
        agree = (np.sign(s.acetyl_dRg_pct) == np.sign(s.acetyl_dtotal)).mean()
        if pfile:
            P = pd.read_csv(pfile).drop_duplicates(['UniProt_ID', 'pos'])
            P = P[P.ess > ESS_MIN]
            rho_pc, _ = stats.spearmanr(P.phospho_dRg_pct, P.phospho_dtotal)
            cm = pd.concat([s.groupby('UniProt_ID').acetyl_dtotal.median(), P.groupby('UniProt_ID').phospho_dtotal.median()], axis=1).dropna()
            rho_ap, _ = stats.spearmanr(cm.acetyl_dtotal, cm.phospho_dtotal)
        else:
            rho_pc, rho_ap, cm = np.nan, np.nan, []
        # per chain: sites of one protein that move the response in opposite directions
        both = s.groupby('UniProt_ID').acetyl_dtotal.agg(lambda v: (v.max() > 0) and (v.min() < 0))
        sig_both = s[sig].groupby('UniProt_ID').acetyl_dtotal.agg(lambda v: (v.max() > 0) and (v.min() < 0))
        r = {'model': model, 'chains': S.UniProt_ID.nunique(), 'chains_with_K': s.UniProt_ID.nunique(), 'sites': len(S),
             'frac_ess_ok': ok.mean(), 'sites_ok': int(ok.sum()),
             'abs_change_median': s.acetyl_dtotal.abs().median(), 'q05': s.acetyl_dtotal.quantile(0.05),
             'q95': s.acetyl_dtotal.quantile(0.95), 'rel_to_total_median': (s.acetyl_dtotal.abs() / s.total_rel.abs()).median(),
             'frac_beyond_2se': sig.mean(), 'frac_compacting': (s.acetyl_dRg_pct < 0).mean(),
             'abs_dRg_pct_median': s.acetyl_dRg_pct.abs().median(), 'ess_min': S.ess.min(),
             'frac_compacting_pos': (s[s.NCPR > 0].acetyl_dRg_pct < 0).mean(), 'frac_expanding_neg': (s[s.NCPR < 0].acetyl_dRg_pct > 0).mean(),
             'rho_compaction_rule': rho_c, 'p_compaction_rule': p_c, 'frac_sign_agree': agree,
             'rho_compaction_rule_phospho_same_ess': rho_pc, 'rho_change_vs_NCPR': rho_q, 'p_change_vs_NCPR': p_q,
             'rho_dRg_vs_NCPR': rho_qr, 'rho_chain_median_acetyl_vs_phospho': rho_ap, 'chains_compared': len(cm),
             'chains_sites_both_directions': float(both.mean()), 'chains_sig_sites_both_directions': int(sig_both.sum())}
        rows.append(r)
        lines += [f'== {model}: heat x acetylation (per Lys, change in the total response, % per K)',
                  f"{r['sites']} Lys in {r['chains_with_K']} chains; effective sample > {ESS_MIN} at {r['frac_ess_ok']:.1%} ({r['sites_ok']} sites)",
                  f"|change| median {r['abs_change_median']:.4f} (5-95%: {r['q05']:.3f} to {r['q95']:.3f}); relative to the chain total: "
                  f"median {r['rel_to_total_median']:.1%}; beyond 2 SE at {r['frac_beyond_2se']:.0%} of usable sites",
                  f"size: acetylation compacts the chain at {r['frac_compacting']:.0%} of sites; |dRg| median {r['abs_dRg_pct_median']:.3f}%; "
                  f"compacts net-positive chains at {r['frac_compacting_pos']:.0%}, expands net-negative at {r['frac_expanding_neg']:.0%}; min ESS {r['ess_min']:.2f}",
                  f"compaction rule: Spearman(dRg, dtotal) {rho_c:.2f} (P {p_c:.1e}); signs agree at {agree:.0%} of sites; "
                  f"phosphorylation, same filter: {rho_pc:.2f}",
                  f"sign: Spearman of the change with chain NCPR {rho_q:.2f} (P {p_q:.1e}); of dRg with NCPR {rho_qr:.2f}",
                  f"chain-median change, acetylation vs phosphorylation: Spearman {rho_ap:.2f} over {len(cm)} chains",
                  f"chains whose Lys move the response both ways: {r['chains_sites_both_directions']:.0%}; "
                  f"with sites beyond 2 SE in both directions: {r['chains_sig_sites_both_directions']}", '']
    a = {k: pd.read_csv(f'{DIR}/sites_{k}.csv').drop_duplicates(['UniProt_ID', 'pos']) for k in ('C', 'R')}
    x = a['C'].merge(a['R'], on=['UniProt_ID', 'pos'], suffixes=('_c', '_r'))
    x = x[(x.ess_c > ESS_MIN) & (x.ess_r > ESS_MIN)]
    rho_x = stats.spearmanr(x.acetyl_dtotal_c, x.acetyl_dtotal_r)[0]
    rho_xr = stats.spearmanr(x.acetyl_dRg_pct_c, x.acetyl_dRg_pct_r)[0]
    same = (np.sign(x.acetyl_dtotal_c) == np.sign(x.acetyl_dtotal_r)).mean()
    lines += ['== CALVADOS-2 vs Mpipi-Recharged, same Lys',
              f"{len(x)} sites in {x.UniProt_ID.nunique()} chains: Spearman of the change of the total response {rho_x:.2f} (same sign "
              f"{same:.0%}); of the change of Rg {rho_xr:.2f}", '']
    rows.append({'model': 'C2 vs MR', 'sites': len(x), 'chains': x.UniProt_ID.nunique(), 'rho_dtotal': rho_x, 'rho_dRg': rho_xr,
                 'frac_same_sign': same})
    # with and without the hydration forcings (acetyl_heat.py rerun with the fixed-solvent change, heat_acetyl/fixed_solvent/)
    FS = f'{DIR}/fixed_solvent' if os.path.isdir(f'{DIR}/fixed_solvent') else DIR
    lines += ['== With and without the hydration forcings (fixed solvent: thermal force without the forcings)']
    for k, (_, model) in PHOS.items():
        f = pd.read_csv(f'{FS}/sites_{k}.csv').drop_duplicates(['UniProt_ID', 'pos'])
        f = f.merge(pd.read_csv(f'{FS}/chains_{k}.csv').drop_duplicates('UniProt_ID')[['UniProt_ID', 'chi_rel', 'total_rel']],
                    on='UniProt_ID')
        o = pd.read_csv(f'{DIR}/sites_{k}.csv').drop_duplicates(['UniProt_ID', 'pos']).merge(f, on=['UniProt_ID', 'pos'],
                                                                                              suffixes=('_o', ''))
        f = f[f.ess > ESS_MIN]
        mag_chi = (f.chi_rel + f.acetyl_dchi).abs() - f.chi_rel.abs()
        mag_tot = (f.total_rel + f.acetyl_dtotal).abs() - f.total_rel.abs()
        own = FS == DIR                               # the run carries its own fixed-solvent change: nothing separate to compare
        r = {'model': f'{model}, forcings', 'sites': len(f),
             'max_abs_diff_dtotal_vs_stored': float('nan') if own else float(np.abs(o.acetyl_dtotal_o - o.acetyl_dtotal).max()),
             'frac_chains_chi_positive': float(f.drop_duplicates('UniProt_ID').chi_rel.gt(0).mean()),
             'rho_dRg_dchi': stats.spearmanr(f.acetyl_dRg_pct, f.acetyl_dchi)[0],
             'rho_dRg_dmag_chi': stats.spearmanr(f.acetyl_dRg_pct, mag_chi)[0],
             'rho_dRg_dmag_total': stats.spearmanr(f.acetyl_dRg_pct, mag_tot)[0]}
        rows.append(r)
        chk = 'one run gives both changes' if own else f"stored total changes reproduced to {r['max_abs_diff_dtotal_vs_stored']:.1e}"
        lines += [f"{model}: {len(f)} Lys; {chk}; chains expanding on "
                  f"heating without forcings {r['frac_chains_chi_positive']:.0%}; Spearman with dRg of the change of the fixed-solvent "
                  f"susceptibility {r['rho_dRg_dchi']:.2f}, of its magnitude {r['rho_dRg_dmag_chi']:.2f}, of the magnitude of the total "
                  f"response {r['rho_dRg_dmag_total']:.2f}"]
    pd.DataFrame(rows).to_csv(f'{DIR}/analysis_by_model.csv', index=False)
    open(f'{DIR}/analysis_summary.txt', 'w').write('\n'.join(lines))
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
