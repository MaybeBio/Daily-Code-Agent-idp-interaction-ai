"""
fig_acetyl.py — Supplementary figure on heat x Lys acetylation (acetyl_heat.py), exact reweighting at 300 K, stored
results only; every Lys has an effective sample fraction above 0.4, so no site is filtered.
  a  change of the total response of Rg against the change of chain size on acetylating one Lys, in CALVADOS-2,
     Mpipi-GG and Mpipi-Recharged, with the Spearman correlation (the compaction rule, as Fig. 6a for phosphorylation)
  b  share of Lys whose acetylation compacts the chain, by the net charge per residue of the chain
  c  median change of the total response on acetylation, by the same net-charge bins
     (b, c: points at the median net charge of each bin; 95% intervals from 1,000 bootstrap resamples of chains)
  d  chain-median change on acetylation against that on phosphorylation of one Ser/Thr (CALVADOS-2 and Mpipi-GG,
     phosphosites with effective sample > 0.1 as heat_crosses_analysis.py)
  e  the same Lys in CALVADOS-2 and Mpipi-Recharged (405 shared chains): change of Rg (left) and of the total response
     (right), with the Spearman correlation
Inputs: Analysis_FDT/total_response/heat_acetyl/sites_{C,M,R}.csv, heat_crosses/sites.csv, forcing_role/phospho_M_sites.csv
    python fig_acetyl.py            (openmpipi env)
Output: BENDER-MT/manuscript/figs/fig_acetyl.{pdf,png}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402

TR = FX.TR
AC = f'{TR}/heat_acetyl'
BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
OUT = f'{FX.ROOT}/manuscript/figs/fig_acetyl'
MODELS = [('C', 'CALVADOS-2'), ('M', 'Mpipi-GG'), ('R', 'Mpipi-Recharged')]
PHOS = {'C': f'{TR}/heat_crosses/sites.csv', 'M': f'{TR}/forcing_role/phospho_M_sites.csv'}   # calibrated model (2026-10-25)
EDGES = [-1, -0.15, -0.1, -0.05, -0.02, 0.02, 0.05, 0.1, 0.15, 0.2, 1]
EXT = (-4.2, 7.6, -0.29, 0.27)                     # every site of every model inside the axes (asserted; calibrated model)
N_BOOT = 1000


def load():
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'NCPR'], low_memory=False)
    S = {k: pd.read_csv(f'{AC}/sites_{k}.csv').drop_duplicates(['UniProt_ID', 'pos']).merge(bio, on='UniProt_ID') for k, _ in MODELS}
    for k, s in S.items():
        assert s.ess.min() > 0.4, (k, s.ess.min())
    return S


def tick_minus(ax, axis='y', fmt='{:g}'):
    f = matplotlib.ticker.FuncFormatter(lambda v, _: '0' if v == 0 else FX.minus(fmt.format(v)))
    (ax.yaxis if axis == 'y' else ax.xaxis).set_major_formatter(f)


def rule_panel(ax, s, model, first):
    x, y = s.acetyl_dRg_pct, s.acetyl_dtotal
    assert x.min() > EXT[0] and x.max() < EXT[1] and y.min() > EXT[2] and y.max() < EXT[3], f'{model}: sites outside the axes'
    ax.hexbin(x, y, gridsize=(30, 30), bins='log', cmap=FX.ff_cmap(FX.COL[model]), mincnt=1, linewidths=0, extent=EXT)
    ax.axhline(0, color='0.55', lw=0.5, zorder=0)
    ax.axvline(0, color='0.55', lw=0.5, zorder=0)
    rho = stats.spearmanr(x, y)[0]
    ax.text(0.97, 0.04, f'ρ = {rho:.2f}\n{len(s):,} Lys', transform=ax.transAxes, ha='right', va='bottom')
    ax.set_xlim(EXT[:2])
    ax.set_ylim(EXT[2:])
    ax.set_xticks([-4, -2, 0, 2, 4, 6])
    ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
    tick_minus(ax, 'x')
    tick_minus(ax, 'y')
    ax.set_title(model, loc='left', color=FX.COL[model])
    ax.set_xlabel(f'Change of {FX.RG} (%)')
    if first:
        ax.set_ylabel('Change of total\nresponse (% per K)')
    else:
        ax.tick_params(labelleft=False)
    return rho


def binned(s, col, stat, rng):
    """Per net-charge bin: x at the median NCPR of its sites, the statistic, and a chain-bootstrap 95% interval."""
    s = s.assign(bin=pd.cut(s.NCPR, EDGES))
    rows = []
    for b, g in s.groupby('bin', observed=True):
        ch = {u: v[col].to_numpy() for u, v in g.groupby('UniProt_ID')}
        keys = list(ch)
        boot = [stat(np.concatenate([ch[keys[i]] for i in rng.integers(0, len(keys), len(keys))])) for _ in range(N_BOOT)]
        rows.append((g.NCPR.median(), stat(g[col].to_numpy()), *np.percentile(boot, [2.5, 97.5]), len(g)))
    return np.array(rows)


def charge_panels(axb, axc, S):
    rng = np.random.default_rng(19)
    xs_all = []
    for k, model in MODELS:
        col = FX.COL[model]
        b = binned(S[k], 'acetyl_dRg_pct', lambda v: 100 * (v < 0).mean(), rng)
        c = binned(S[k], 'acetyl_dtotal', np.median, rng)
        for ax, r in ((axb, b), (axc, c)):
            ax.fill_between(r[:, 0], r[:, 2], r[:, 3], color=FX.final_line(col), alpha=0.15, lw=0)   # line shade: a pale fill colour vanishes
            xs_all.extend(r[:, 0])
            ax.plot(r[:, 0], r[:, 1], '-o', color=col, ms=2.6, lw=0.9, label=model)
    for ax in (axb, axc):
        ax.axvline(0, color='0.55', lw=0.5, zorder=0)
        ax.set_xlim(min(-0.33, min(xs_all) - 0.02), max(0.33, max(xs_all) + 0.02))   # every bin inside
        ax.set_xticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3])
        tick_minus(ax, 'x')
        ax.set_xlabel('Net charge per residue of the chain')
    axb.set_ylim(-3, 103)
    axb.set_yticks([0, 25, 50, 75, 100])
    axb.set_ylabel('Lys whose acetylation\ncompacts the chain (%)')
    axb.legend(loc='upper left', handlelength=1.6)
    axc.axhline(0, color='0.55', lw=0.5, zorder=0)
    axc.set_ylabel('Median change of total\nresponse (% per K)')
    tick_minus(axc, 'y', '{:.3f}')


def phospho_panel(ax, S):
    out = {}
    for k, model in MODELS[:2]:
        P = pd.read_csv(PHOS[k]).drop_duplicates(['UniProt_ID', 'pos'])
        P = P[P.ess > 0.1]
        cm = pd.concat([S[k].groupby('UniProt_ID').acetyl_dtotal.median(), P.groupby('UniProt_ID').phospho_dtotal.median()], axis=1).dropna()
        rho = stats.spearmanr(cm.acetyl_dtotal, cm.phospho_dtotal)[0]
        ax.scatter(cm.phospho_dtotal, cm.acetyl_dtotal, s=4, lw=0, color=FX.COL[model], alpha=0.75, zorder=3)
        ax.text(0.03, 0.97 - 0.085 * len(out), f'ρ = {rho:.2f}', color=FX.COL[model], transform=ax.transAxes, ha='left',
                va='top')
        out[model] = (rho, len(cm), cm)
    allc = pd.concat([v[2] for v in out.values()])
    lo = min(allc.min()) * 1.08
    hi = max(allc.max()) * 1.08
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.plot([lo, hi], [lo, hi], color='0.7', lw=0.5, ls='--', zorder=0)
    ax.axhline(0, color='0.55', lw=0.5, zorder=0)
    ax.axvline(0, color='0.55', lw=0.5, zorder=0)
    tick_minus(ax, 'x')
    tick_minus(ax, 'y')
    ax.set_xlabel('Phosphorylation (% per K)')
    ax.set_ylabel('Acetylation (% per K)')
    ax.set_title('Median change per chain', loc='left')
    ax.set_xticks([-0.1, 0, 0.1, 0.2])
    ax.set_yticks([-0.1, 0, 0.1, 0.2])
    return {m: (round(v[0], 2), v[1]) for m, v in out.items()}


def cross_panels(ax1, ax2, S):
    x = S['C'].merge(S['R'], on=['UniProt_ID', 'pos'], suffixes=('_c', '_r'))
    res = {}
    for ax, col, ext, ticks in ((ax1, 'acetyl_dRg_pct', (-4.2, 7.6), [-4, -2, 0, 2, 4, 6]),
                                (ax2, 'acetyl_dtotal', (-0.29, 0.27), [-0.2, -0.1, 0, 0.1, 0.2])):
        a, b = x[f'{col}_c'], x[f'{col}_r']
        assert a.min() > ext[0] and a.max() < ext[1] and b.min() > ext[0] and b.max() < ext[1], col
        ax.hexbin(a, b, gridsize=28, bins='log', cmap=FX.ff_cmap(FX.COL['Mpipi-Recharged']), mincnt=1, linewidths=0,
                  extent=ext + ext)
        ax.plot(ext, ext, color='0.7', lw=0.5, ls='--', zorder=0)
        ax.axhline(0, color='0.55', lw=0.5, zorder=0)
        ax.axvline(0, color='0.55', lw=0.5, zorder=0)
        rho = stats.spearmanr(a, b)[0]
        ax.text(0.97, 0.04, f'ρ = {rho:.2f}\n{len(x):,} Lys', transform=ax.transAxes, ha='right', va='bottom')
        ax.set_xlim(ext)
        ax.set_ylim(ext)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        tick_minus(ax, 'x')
        tick_minus(ax, 'y')
        ax.set_title(f'Same Lys, change of {FX.RG}' if col == 'acetyl_dRg_pct' else 'Same Lys, total response',
                     loc='left')
        res[col] = round(rho, 2)
    ax1.set_xlabel('CALVADOS-2 (%)')
    ax1.set_ylabel('Mpipi-Recharged (%)')
    ax2.set_xlabel('CALVADOS-2 (% per K)')
    ax2.set_ylabel('Mpipi-Recharged (% per K)')
    return res, len(x)


def main():
    FX.house()
    S = load()
    fig = plt.figure(figsize=(FX.WIDTH, 6.6))
    r1 = fig.add_gridspec(1, 3, left=0.105, right=0.985, top=0.95, bottom=0.715, wspace=0.08)
    a = [fig.add_subplot(r1[0, i]) for i in range(3)]
    rho_rule = {m: round(rule_panel(a[i], S[k], m, i == 0), 2) for i, (k, m) in enumerate(MODELS)}
    r2 = fig.add_gridspec(1, 2, left=0.105, right=0.985, top=0.615, bottom=0.405, wspace=0.32)
    axb, axc = fig.add_subplot(r2[0, 0]), fig.add_subplot(r2[0, 1])
    charge_panels(axb, axc, S)
    r3 = fig.add_gridspec(1, 3, left=0.105, right=0.985, top=0.3, bottom=0.07, wspace=0.62)
    axd, axe1, axe2 = (fig.add_subplot(r3[0, i]) for i in range(3))
    rho_ap = phospho_panel(axd, S)
    rho_x, n_x = cross_panels(axe1, axe2, S)
    FX.add_letters(fig, [(a[0], 'a', {'row': 1}), (axb, 'b', {'row': 2}), (axc, 'c', {'row': 2}),
                         (axd, 'd', {'row': 3}), (axe1, 'e', {'row': 3})])
    FX.save(fig, OUT)
    out = {'rule_rho': rho_rule, 'acetyl_vs_phospho': rho_ap, 'C2_vs_MR': rho_x, 'shared_Lys': n_x,
           'Lys': {m: len(S[k]) for k, m in MODELS}}
    print('written', OUT, out)
    return out


if __name__ == '__main__':
    main()
