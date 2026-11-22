"""
heat_crosses_ddof_fix.py — one-off correction (2026-10-10) of heat_crosses.py outputs written before its ddof fix.

heat_crosses.perturbed() subtracted a baseline covariance with ddof 1 (fdt._cov) from a reweighted covariance with ddof 0
(wcov), so every phospho_dtotal and salt_dtotal was shifted by -total_rel / n (n = 1000 stored frames per chain; a spurious
weakening of contraction of |total| / 1000 % per K). The point estimates are corrected exactly by adding total_rel / 1000;
the jackknife SEs are unchanged to < 1e-6 and kept. CALVADOS-2 files only (Mpipi-GG used ddof 0 on both sides).
Originals are copied to <dir>/before_ddoffix_2026-10-10/; a marker file prevents a second application.
    python heat_crosses_ddof_fix.py
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402
import shutil  # noqa: E402

import pandas as pd  # noqa: E402

TR = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/total_response'
DIRS = [f'{TR}/heat_crosses', f'{TR}/heat_crosses_shared']
N_FRAMES = 1000


def main():
    for d in DIRS:
        mark = os.path.join(d, 'ddof_fixed_2026-10-10.txt')
        if os.path.exists(mark):
            print(f'{d}: already corrected, skipped')
            continue
        bak = os.path.join(d, 'before_ddoffix_2026-10-10')
        os.makedirs(bak, exist_ok=True)
        for f in ('chains.csv', 'sites.csv'):
            shutil.copy2(os.path.join(d, f), os.path.join(bak, f))
        c = pd.read_csv(os.path.join(d, 'chains.csv'))
        tot = c.drop_duplicates('UniProt_ID').set_index('UniProt_ID').total_rel
        c['salt_dtotal'] = c.salt_dtotal + c.UniProt_ID.map(tot) / N_FRAMES
        s = pd.read_csv(os.path.join(d, 'sites.csv'))
        if s.UniProt_ID.map(tot).isna().any():
            raise ValueError(f'{d}: sites without a chain total')
        s['phospho_dtotal'] = s.phospho_dtotal + s.UniProt_ID.map(tot) / N_FRAMES
        c.to_csv(os.path.join(d, 'chains.csv'), index=False)
        s.to_csv(os.path.join(d, 'sites.csv'), index=False)
        open(mark, 'w').write('salt_dtotal and phospho_dtotal += total_rel / 1000 (heat_crosses_ddof_fix.py); '
                              'originals in before_ddoffix_2026-10-10/\n')
        print(f'{d}: corrected {len(c)} chains, {len(s)} sites; median shift {(tot / N_FRAMES).median():+.5f} % per K')


if __name__ == '__main__':
    main()
