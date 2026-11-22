"""
compaction_length_check.py — does "compact chains respond most" survive at fixed chain length and at fixed size of the
energy fluctuations? The FDT response is -Cov(Rg, force) / <Rg> = -r x CV(Rg) x sd(force), and sd(force) grows about as
sqrt(N); compact BENDER chains are also longer. So the compaction dependence could in part be a length or
fluctuation-size effect. Checked on the totals (total_response_evidence.load: chi + TEA + hydrophilic term, frozen
strengths) and on the plain chi, per force field:
  rank correlation with own nu: raw; partial given log N; partial given the fluctuation sizes; partial given both;
  within length bands (40-80, 80-120, 120-200, 200-400 residues), with medians by nu bin.
Fluctuation sizes per chain (only ranks are used, so units do not matter):
  plain force  Mpipi-Recharged sd(H) (recharged_fdt_300K.csv H_sd_kJ); Mpipi-GG sd(potential energy) (pe_sd_kcal);
               CALVADOS-2 sd(Ashbaugh-Hatch energy) from chi_factor_split/per_protein.csv (3,181 chains; no total
               energy SD is stored for CALVADOS-2, the AH term is its main fluctuating part)
  hydration    sd of the hydrophilic forcing energy (D_sd_hyd, as total_response_evidence.factors)
    BENDER_REANALYSIS_DIR=... python compaction_length_check.py       (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/length_check/{partials.csv, bands.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import total_response_evidence as T  # noqa: E402

TLA = T.TLA
OUT = os.path.join(TLA, 'Analysis_FDT', 'total_response', 'length_check')
NAME = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
BANDS = [(40, 80), (80, 120), (120, 200), (200, 400)]
NUB = [0.30, 0.50, 0.54, 0.58, 0.63, 1.0]


def rk(x):
    return stats.rankdata(x)


def partial(y, x, Z):
    """Spearman partial correlation of y and x given the columns of Z (rank-residualized least squares)."""
    if Z is None or Z.shape[1] == 0:
        return float(stats.spearmanr(y, x).statistic)
    Zr = np.column_stack([rk(z) for z in Z.T] + [np.ones(len(y))])
    ry = rk(y) - Zr @ np.linalg.lstsq(Zr, rk(y), rcond=None)[0]
    rx = rk(x) - Zr @ np.linalg.lstsq(Zr, rk(x), rcond=None)[0]
    return float(np.corrcoef(ry, rx)[0, 1])


def plain_sd():
    r = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'recharged_fdt_300K.csv'))
    m = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K_fdt.csv'))
    c = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'compact', 'chi_factor_split', 'per_protein.csv'))
    return {'R': r[r.metric == 'rg'].set_index('UniProt_ID').H_sd_kJ,
            'M': m[m.metric == 'rg'].set_index('UniProt_ID').pe_sd_kcal,
            'C': c[c.model == 'CALVADOS-2'].set_index('UniProt_ID').sdE_kT}


def main():
    os.makedirs(OUT, exist_ok=True)
    D, _ = T.load()
    sdp = plain_sd()
    prow, brow = [], []
    for k in ('C', 'R', 'M'):
        d = D[k].copy()
        d['sd_plain'] = sdp[k].reindex(d.index)
        d = d.dropna(subset=['nu', 'N', 'total', 'chi', 'D_sd_hyd'])
        for resp in ('total', 'chi'):
            y = -d[resp].values          # stronger contraction = larger
            ctrl = {'none': None, 'log N': np.column_stack([np.log(d.N)]),
                    'sd(hydration forcing)': np.column_stack([d.D_sd_hyd])}
            s = d.dropna(subset=['sd_plain'])
            ys = -s[resp].values
            for lab, Z in ctrl.items():
                prow.append({'model': NAME[k], 'response': resp, 'control': lab, 'n': len(d),
                             'rho_nu': partial(y, d.nu.values, Z)})
            for lab, Z in {'sd(plain force)': np.column_stack([s.sd_plain]),
                           'sd(plain) + sd(hydration)': np.column_stack([s.sd_plain, s.D_sd_hyd]),
                           'log N + both sd': np.column_stack([np.log(s.N), s.sd_plain, s.D_sd_hyd])}.items():
                prow.append({'model': NAME[k], 'response': resp, 'control': lab, 'n': len(s),
                             'rho_nu': partial(ys, s.nu.values, Z)})
            prow.append({'model': NAME[k], 'response': resp, 'control': 'length effect at fixed nu (rho with log N)',
                         'n': len(d), 'rho_nu': partial(y, np.log(d.N.values), np.column_stack([d.nu]))})
            for lo, hi in BANDS:
                b = d[(d.N >= lo) & (d.N < hi)]
                if len(b) < 30:
                    continue
                med = b.groupby(pd.cut(b.nu, NUB), observed=False)[resp].median()
                cnt = b.groupby(pd.cut(b.nu, NUB), observed=False)[resp].size()
                brow.append({'model': NAME[k], 'response': resp, 'N_band': f'{lo}-{hi}', 'n': len(b),
                             'rho_nu': float(stats.spearmanr(-b[resp], b.nu).statistic),
                             **{f'med {iv}': (v if cnt[iv] >= 10 else np.nan) for iv, v in med.items()}})
    P, B = pd.DataFrame(prow), pd.DataFrame(brow)
    P.to_csv(os.path.join(OUT, 'partials.csv'), index=False)
    B.to_csv(os.path.join(OUT, 'bands.csv'), index=False)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 20)
    L = ['Compaction dependence of the heating response at fixed length and fluctuation size',
         'rho_nu: Spearman (partial) correlation of contraction (-response) with own nu; negative = compact chains '
         'contract more', '',
         P.pivot_table(index=['model', 'response'], columns='control', values='rho_nu', sort=False).round(2).to_string(),
         '', 'n per model: ' + ', '.join(f"{m} {g.n.max()} (with plain sd {g.n.min()})" for m, g in P.groupby('model')),
         '', 'within length bands: rho and median response (% per K) by own nu bin (bins with >= 10 chains)',
         B.round(3).to_string(index=False)]
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
