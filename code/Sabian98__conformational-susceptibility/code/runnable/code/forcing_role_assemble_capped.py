"""
forcing_role_assemble_capped.py — the CALVADOS-2 tables of forcing_role.py under the capped forcing (the main model from
2026-10-24), assembled from runs that already hold every number instead of a third pass over the same trajectories:
  calibrated change (with the forcings)   heat_crosses.py --forcing capped (total_response/heat_crosses/: sites.csv for
                                           phosphorylation, chains.csv for doubled salt, and each chain's total)
  plain change (without the forcings)     forcing_role.py's A0 run (total_response/forcing_role_A0_linear/ or forcing_role/):
                                           the fixed-solvent change contains no forcing, so it is the same under any forcing
Checks (the script stops otherwise): the same chains and sites in both sources; the calibrated values equal a direct
capped forcing_role.py run on the chains it covers (phosphorylation: the partial run of 2026-10-25, salt: 20 chains run
for this check) to 1e-12, with the same plain values.
    BENDER_REANALYSIS_DIR=... python forcing_role_assemble_capped.py --check-phospho DIR --check-salt DIR
Output: total_response/forcing_role/ (the canonical folder holds the calibrated model since the 2026-10-25 swap; the A0 run
        is in forcing_role_A0_linear) {phospho_C_chains.csv, phospho_C_sites.csv, salt_C_chains.csv}
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

TR = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis/Analysis_FDT/total_response')
OUT = os.path.join(TR, 'forcing_role' if os.path.isdir(os.path.join(TR, 'forcing_role_A0_linear')) else 'forcing_role_capped')


def a0_dir():
    d = os.path.join(TR, 'forcing_role_A0_linear')
    return d if os.path.isdir(d) else os.path.join(TR, 'forcing_role')


def same(a, b, key, cols, what):
    m = a.merge(b, on=key, suffixes=('', '_ref'))
    if len(m) != len(b):
        raise SystemExit(f'{what}: {len(b) - len(m)} reference rows not in the assembled table')
    dev = max(float((m[c] - m[f'{c}_ref']).abs().max()) for c in cols)
    if not dev < 1e-12:
        raise SystemExit(f'{what}: assembled values differ from the direct capped run ({dev:.1e})')
    return len(m), dev


def main(check_phospho, check_salt):
    a0 = a0_dir()
    hs = pd.read_csv(os.path.join(TR, 'heat_crosses', 'sites.csv')).drop_duplicates(['UniProt_ID', 'pos'])
    hc = pd.read_csv(os.path.join(TR, 'heat_crosses', 'chains.csv')).drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    cal_s = ['phospho_dtotal', 'phospho_dtotal_se', 'ess', 'phospho_dRg_pct']
    ps = pd.read_csv(os.path.join(a0, 'phospho_C_sites.csv'))
    m = ps.drop(columns=cal_s).merge(hs[['UniProt_ID', 'pos'] + cal_s], on=['UniProt_ID', 'pos'], how='left')
    if m[cal_s].isna().any().any() or len(m) != len(ps):
        raise SystemExit('phosphosites of the A0 forcing_role run missing from heat_crosses')
    sites = m[ps.columns]
    pc = pd.read_csv(os.path.join(a0, 'phospho_C_chains.csv'))
    if not set(pc.UniProt_ID) <= set(hc.index):
        raise SystemExit('phospho chains missing from heat_crosses')
    pc['total_rel'] = hc.loc[pc.UniProt_ID, 'total_rel'].to_numpy()
    sc = pd.read_csv(os.path.join(a0, 'salt_C_chains.csv'))
    if not set(sc.UniProt_ID) <= set(hc.index):
        raise SystemExit('salt chains missing from heat_crosses')
    for c in ('total_rel', 'salt_dtotal', 'salt_dtotal_se', 'salt_ess', 'salt_dRg_pct'):
        sc[c] = hc.loc[sc.UniProt_ID, c].to_numpy()
    L = []
    ref = pd.read_csv(os.path.join(check_phospho, 'phospho_C_sites.csv'))
    L.append('phospho sites vs direct capped run: %d rows, max |difference| %.1e' %
             same(sites, ref, ['UniProt_ID', 'pos'], cal_s + ['phospho_dchi'], 'phospho sites'))
    ref = pd.read_csv(os.path.join(check_phospho, 'phospho_C_chains.csv'))
    L.append('phospho chains vs direct capped run: %d rows, max |difference| %.1e' %
             same(pc, ref, ['UniProt_ID'], ['total_rel', 'chi_rel'], 'phospho chains'))
    ref = pd.read_csv(os.path.join(check_salt, 'salt_C_chains.csv'))
    L.append('salt chains vs direct capped run: %d rows, max |difference| %.1e' %
             same(sc, ref, ['UniProt_ID'], ['total_rel', 'salt_dtotal', 'salt_dtotal_se', 'salt_dRg_pct', 'salt_dchi'], 'salt'))
    os.makedirs(OUT, exist_ok=True)
    sites.to_csv(os.path.join(OUT, 'phospho_C_sites.csv'), index=False)
    pc.to_csv(os.path.join(OUT, 'phospho_C_chains.csv'), index=False)
    sc.to_csv(os.path.join(OUT, 'salt_C_chains.csv'), index=False)
    open(os.path.join(OUT, 'assembled_C.txt'), 'w').write(
        'CALVADOS-2 tables assembled by forcing_role_assemble_capped.py (calibrated: heat_crosses capped; plain: A0 run)\n'
        + '\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check-phospho', required=True)
    ap.add_argument('--check-salt', required=True)
    a = ap.parse_args()
    main(a.check_phospho, a.check_salt)
