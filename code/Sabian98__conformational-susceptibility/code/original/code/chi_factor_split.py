"""
chi_factor_split.py — why does the Mpipi-GG sticker susceptibility keep rising below the theta point while the
CALVADOS-2 one saturates? Split the sticker term into three factors per protein, at 300 K:

    chi_s / <Rg>  =  r(Rg, E_s) * CV(Rg) * [sd(E_s) / RT] / T

  r        correlation of Rg with the sticker energy over the 1,000 FDT frames (how tightly size and contacts move
           together)
  CV(Rg)   sd(Rg) / <Rg> (how much the chain breathes)
  sd(E_s)/RT  sticker-energy fluctuation in kT (how many contacts come and go, times their strength)
with E_s = Wang-Frenkel energy (Mpipi-GG) or Ashbaugh-Hatch energy (CALVADOS-2), from the FDT npz files.
Also the mean sticker energy per residue in kT (contact density) and sd(E_s)/RT per sqrt(N).
chi_s is recomputed from the frames and checked against chi_WF / chi_AH in the FDT tables.

Per own-nu bin (power-law fit R^2 >= 0.99 in that model): medians, and an exact log split of the rise between the
bin 0.50-0.54 and each more compact bin (geometric means over proteins with r > 0), so the three factor shares
add up to the change in chi_s.

Sets: Mpipi-GG = 992-set + compact set (all with the fit criterion); CALVADOS-2 = up to 600 proteins per nu bin
(random, seed 0) from the 10,282, all below nu 0.50.
Output: Traj Level Analysis/Analysis_FDT_MPIPI/compact/chi_factor_split/ (per_protein.csv, bins.csv, log_split.csv)
Usage: python3 chi_factor_split.py
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import os  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', '/Users/taseefrahman/Downloads/Scripts/bender_reanalysis'))
import bender_common as bc  # noqa: E402

MT = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
MP = f'{MT}/Analysis_FDT_MPIPI'
CAL = f'{MT}/Analysis_FDT'
OUT_DIR = f'{MP}/compact/chi_factor_split'
T = 300.0
R = {'Mpipi-GG': 1.98720425864083e-3, 'CALVADOS-2': 8.314462618e-3}   # kcal / kJ per mol per K
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
REF = '(0.5, 0.54]'
PER_BIN, SEED, R2_MIN = 600, 0, 0.99


def factors(path, key, rgas):
    z = np.load(path)
    rg, e = z['rg'], z[key]
    n = len(str(z['sequence']))
    c = np.cov(rg, e)
    return {'N': n, 'chi_s_rel': 100 * c[0, 1] / (rgas * T ** 2) / rg.mean(), 'r': c[0, 1] / np.sqrt(c[0, 0] * c[1, 1]),
            'cv_rg': rg.std(ddof=1) / rg.mean(), 'sdE_kT': np.sqrt(c[1, 1]) / (rgas * T),
            'meanE_kT_per_res': e.mean() / (rgas * T) / n, 'chi_s_check': c[0, 1] / (rgas * T ** 2)}


def load_mpipi():
    rows = []
    for d, tag, npz in ((MP, 'mpipi_300K', f'{MP}/npz_300K'),
                        (f'{MP}/compact', 'mpipi_compact_300K', f'{MP}/compact/npz_300K')):
        m = pd.read_csv(f'{d}/{tag}.csv').set_index('UniProt_ID')
        f = pd.read_csv(f'{d}/{tag}_fdt.csv')
        f = f[f['metric'] == 'rg'].set_index('UniProt_ID')
        m = m[(m['status'] == 'ok') & (m['nu_r2'] >= R2_MIN)].join(f['chi_WF'], how='inner')
        rows.append(pd.DataFrame({'nu': m['nu'], 'chi_table': m['chi_WF'],
                                  'path': [f'{npz}/{u}.npz' for u in m.index]}, index=m.index))
    return pd.concat(rows)


def load_calvados():
    bio = pd.read_csv(bc.BIO_CSV, usecols=['UniProt_ID', 'nu', 'nu_fit_r2'], low_memory=False).set_index('UniProt_ID')
    bio = bio[bio['nu_fit_r2'] >= R2_MIN]
    bio['b'] = pd.cut(bio['nu'], BINS)
    pick = [g.sample(min(len(g), PER_BIN), random_state=SEED) if iv.left >= 0.5 else g
            for iv, g in bio.groupby('b', observed=True)]
    s = pd.concat(pick)
    f = pd.read_csv(f'{CAL}/fdt_results_300K.csv', usecols=['UniProt_ID', 'metric', 'chi_AH'])
    f = f[f['metric'] == 'rg'].set_index('UniProt_ID')
    s = s.join(f['chi_AH'], how='inner')
    return pd.DataFrame({'nu': s['nu'], 'chi_table': s['chi_AH'], 'path': [f'{CAL}/npz/300K/{u}.npz' for u in s.index]},
                        index=s.index)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    per = []
    for model, df, key in (('Mpipi-GG', load_mpipi(), 'E_WF'), ('CALVADOS-2', load_calvados(), 'V_AH')):
        with ThreadPoolExecutor(8) as ex:
            res = list(ex.map(lambda p: factors(p, key, R[model]), df['path']))
        d = df.drop(columns='path').join(pd.DataFrame(res, index=df.index))
        d.insert(0, 'model', model)
        ratio = d['chi_s_check'] / d['chi_table']
        print(f'{model}: {len(d)} proteins; recomputed / table sticker chi: median {ratio.median():.4f}, '
              f'IQR {ratio.quantile(.25):.4f}-{ratio.quantile(.75):.4f}')
        per.append(d)
    per = pd.concat(per)
    per['b'] = pd.cut(per['nu'], BINS).astype(str)
    per.to_csv(f'{OUT_DIR}/per_protein.csv')

    cols = ['chi_s_rel', 'r', 'cv_rg', 'sdE_kT', 'meanE_kT_per_res', 'N']
    g = per.groupby(['model', 'b'])
    bins = g[cols].median()
    bins.insert(0, 'n', g.size())
    bins['sdE_kT_per_sqrtN'] = g.apply(lambda x: (x['sdE_kT'] / np.sqrt(x['N'])).median())
    bins.to_csv(f'{OUT_DIR}/bins.csv')
    with pd.option_context('display.width', 250):
        print('\nmedians by own-nu bin (chi_s_rel in % Rg per K)')
        print(bins.round(4).to_string())

    # exact log split of the change vs the reference bin (geometric means over r > 0)
    pos = per[per['r'] > 0].copy()
    for c in ('r', 'cv_rg', 'sdE_kT'):
        pos[f'log_{c}'] = np.log(pos[c])
    gm = pos.groupby(['model', 'b'])[['log_r', 'log_cv_rg', 'log_sdE_kT']].mean()
    rows = []
    for model in per['model'].unique():
        ref = gm.loc[(model, REF)]
        for b in sorted(per.loc[per['model'] == model, 'b'].unique()):
            if b == REF or float(b.split(',')[0][1:]) >= 0.5:
                continue
            dlt = gm.loc[(model, b)] - ref
            tot = dlt.sum()
            rows.append({'model': model, 'bin': b, 'vs': REF, 'fold_change_chi_s': np.exp(tot),
                         'fold r': np.exp(dlt['log_r']), 'fold CV(Rg)': np.exp(dlt['log_cv_rg']),
                         'fold sd(E)/kT': np.exp(dlt['log_sdE_kT']),
                         'share r': dlt['log_r'] / tot, 'share CV(Rg)': dlt['log_cv_rg'] / tot,
                         'share sd(E)': dlt['log_sdE_kT'] / tot,
                         'frac_r_pos': (per.loc[(per['model'] == model) & (per['b'] == b), 'r'] > 0).mean()})
    ls = pd.DataFrame(rows)
    ls.to_csv(f'{OUT_DIR}/log_split.csv', index=False)
    with pd.option_context('display.width', 250):
        print('\nchange of the sticker chi below theta vs nu 0.50-0.54, split into factors (geometric means)')
        print(ls.round(3).to_string(index=False))
    bc.write_provenance(f'{OUT_DIR}/bins.csv', __file__, {'T': T, 'bins': BINS, 'ref_bin': REF, 'per_bin_calvados': PER_BIN,
                                                        'seed': SEED, 'nu_fit_r2_min': R2_MIN},
                        [f'{MP}/mpipi_300K_fdt.csv', f'{MP}/compact/mpipi_compact_300K_fdt.csv', f'{CAL}/fdt_results_300K.csv'])
    print(f'\nsaved {OUT_DIR}')


if __name__ == '__main__':
    main()
