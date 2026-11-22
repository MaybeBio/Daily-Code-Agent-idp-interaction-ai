"""
forcing_role_analysis.py — the compaction rule with the calibrated force fields and with the plain force fields (no
hydration forcings), for phosphorylation, Lys acetylation and doubled salt, in one table.

Inputs: forcing_role.py output (phosphorylation, salt) and acetyl_heat.py's fixed-solvent rerun (heat_acetyl/fixed_solvent/).
Usable: Kish effective sample fraction > 0.1, as throughout. Per perturbation and model:
  rho_calibrated        Spearman of the change of the total response with the change of Rg (the rule as in the paper)
  frac_plain_expanding  share of chains that expand on heating in the plain force field (fixed-solvent susceptibility > 0)
  rho_plain             Spearman of the change of the fixed-solvent susceptibility with the change of Rg
  rho_size_calibrated   Spearman of the change of |total response| with the change of Rg
  rho_size_plain        Spearman of the change of |fixed-solvent susceptibility| with the change of Rg
  noise_ceiling         sqrt(1 - mean(SE^2)/var(change of the total response)), the largest correlation of the change with a
                        noise-free quantity that its jackknife sampling error allows
A negative rho_size means that compaction enlarges the response. The calibrated changes are checked against the stored
values of heat_crosses/ (CALVADOS-2). Since 2026-10-25 the canonical folders hold the calibrated model with saturation
(forcing_role/, heat_acetyl/: acetyl_heat.py --forcing capped, whose site tables carry the fixed-solvent change as
well); the Mpipi checks against heat_crosses_shared/ and traj_pass/ apply to the linear A0 runs only and are reported as
not applicable. A missing acetylation table is skipped.
    python forcing_role_analysis.py
Output: total_response/forcing_role/analysis.csv, analysis_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

TLA = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis')
TR = f'{TLA}/Analysis_FDT/total_response'
D, AC = f'{TR}/forcing_role', f'{TR}/heat_acetyl'
ESS_MIN = 0.1
CALIBRATED_SAT = os.path.isdir(f'{TR}/forcing_role_A0_linear')    # canonical folders hold the saturating model
MODELS = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
STORED = {('phospho', 'C'): (f'{TR}/heat_crosses/sites.csv', ['UniProt_ID', 'pos'], 'phospho_dtotal'),
          ('phospho', 'M'): (f'{TR}/heat_crosses_shared/sites_gg.csv', ['UniProt_ID', 'pos'], 'phospho_dtotal'),
          ('salt', 'C'): (f'{TR}/heat_crosses/chains.csv', ['UniProt_ID'], 'salt_dtotal'),
          ('salt', 'M'): (f'{TR}/heat_crosses_shared/chains_gg.csv', ['UniProt_ID'], 'salt_dtotal'),
          ('salt', 'R'): (f'{TLA}/Analysis_FDT_RECHARGED/traj_pass/chains.csv', ['UniProt_ID'], 'salt_dtotal'),
          ('acetyl', 'C'): (f'{AC}/sites_C.csv', ['UniProt_ID', 'pos'], 'acetyl_dtotal'),
          ('acetyl', 'M'): (f'{AC}/sites_M.csv', ['UniProt_ID', 'pos'], 'acetyl_dtotal'),
          ('acetyl', 'R'): (f'{AC}/sites_R.csv', ['UniProt_ID', 'pos'], 'acetyl_dtotal')}


def load(what, k):
    """Rows with columns dtotal, dchi, dRg, ess, total_rel, chi_rel, UniProt_ID (sites, or chains for salt)."""
    if what == 'salt':
        f = pd.read_csv(f'{D}/salt_{k}_chains.csv').drop_duplicates('UniProt_ID')
        f = f.rename(columns={'salt_dtotal': 'dtotal', 'salt_dtotal_se': 'dtotal_se', 'salt_dchi': 'dchi', 'salt_dRg_pct': 'dRg',
                              'salt_ess': 'ess'})
        keys = ['UniProt_ID']
    else:
        base = f'{D}/phospho_{k}' if what == 'phospho' else None
        acd = f'{AC}/fixed_solvent' if os.path.isdir(f'{AC}/fixed_solvent') else AC
        sites = f'{base}_sites.csv' if base else f'{acd}/sites_{k}.csv'
        chains = f'{base}_chains.csv' if base else f'{acd}/chains_{k}.csv'
        f = pd.read_csv(sites).drop_duplicates(['UniProt_ID', 'pos'])
        f = f.rename(columns=lambda c: c.replace(f'{what}_', '')).rename(columns={'dRg_pct': 'dRg'})
        f = f.merge(pd.read_csv(chains).drop_duplicates('UniProt_ID')[['UniProt_ID', 'chi_rel', 'total_rel']], on='UniProt_ID')
        keys = ['UniProt_ID', 'pos']
    if CALIBRATED_SAT:
        # CALVADOS-2 phosphorylation and salt tables are assembled from heat_crosses (forcing_role_assemble_capped.py),
        # so comparing them with heat_crosses would compare a table with itself; their genuine check against a direct
        # capped forcing_role.py run is in forcing_role/assembled_C.txt. The stored Mpipi and acetylation values are A0.
        return f, float('nan'), 0
    path, sk, col = STORED[(what, k)]
    o = pd.read_csv(path).drop_duplicates(sk)[sk + [col]].merge(f[keys + ['dtotal']], on=keys)
    return f, float(np.abs(o[col] - o.dtotal).max()), len(o)


def main():
    rows, lines = [], ['Compaction rule with the calibrated and with the plain force fields (usable: Kish fraction > 0.1)', '']
    for what, label in (('phospho', 'Phosphorylation'), ('acetyl', 'Lys acetylation'), ('salt', 'Doubled salt')):
        for k in (('C', 'M') if what == 'phospho' else ('C', 'M', 'R')):
            try:
                f, dmax, nchk = load(what, k)
            except FileNotFoundError as ex:
                lines.append(f'{label}, {MODELS[k]}: skipped ({ex.filename} missing)')
                continue
            f = f[f.ess > ESS_MIN]
            size_tot = (f.total_rel + f.dtotal).abs() - f.total_rel.abs()
            size_chi = (f.chi_rel + f.dchi).abs() - f.chi_rel.abs()
            r = {'perturbation': label, 'model': MODELS[k], 'n': len(f), 'unit': 'chains' if what == 'salt' else 'sites',
                 'rho_calibrated': stats.spearmanr(f.dRg, f.dtotal)[0],
                 'frac_plain_expanding': float(f.drop_duplicates('UniProt_ID').chi_rel.gt(0).mean()),
                 'rho_plain': stats.spearmanr(f.dRg, f.dchi)[0],
                 'rho_size_calibrated': stats.spearmanr(f.dRg, size_tot)[0],
                 'rho_size_plain': stats.spearmanr(f.dRg, size_chi)[0],
                 'noise_ceiling': float(np.sqrt(1 - (f.dtotal_se ** 2).mean() / f.dtotal.var())),
                 'check_n': nchk, 'check_max_abs_diff': dmax}
            rows.append(r)
            lines.append(f"{label}, {MODELS[k]}: {r['n']} {r['unit']}; calibrated rule {r['rho_calibrated']:.2f}; plain: "
                         f"{r['frac_plain_expanding']:.0%} of chains expand on heating, rule {r['rho_plain']:.2f}; size of the "
                         f"response vs dRg {r['rho_size_calibrated']:.2f} calibrated, {r['rho_size_plain']:.2f} plain "
                         + (f"(calibrated changes match stored values to {dmax:.1e} over {nchk})" if nchk else
                            "(no independent stored values; CALVADOS-2 checked in assembled_C.txt)")
                         + f"; noise ceiling {r['noise_ceiling']:.2f}")
    pd.DataFrame(rows).to_csv(f'{D}/analysis.csv', index=False)
    open(f'{D}/analysis_summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
