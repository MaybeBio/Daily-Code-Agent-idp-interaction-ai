"""
fig_heldout.py — manuscript figure for the held-out temperature series (capped main model from 2026-10-25).
  a, b  per data set (CALVADOS-2, Mpipi-Recharged): measured relative slope of size (% per K) against the plain model,
        the linear calibration (A0, the original primary test; open markers) and the capped model (filled; PNt from the
        direct calibrated runs, the other series by reweighting, heldout_temperature/capped_forcing/noise_floor.csv; PNt
        at 525 mM, which has no measured slope, by reweighting here). Faded: effective sample below 0.05 at some
        temperature. PNt is fitted for the capped model (its cap value was chosen on PNt).
  c     mean absolute error over the five series that no model fitted (alpha-synuclein x 2, p53, sfAFP, Aap) for the
        plain, linear and capped models, with the error a perfect model would leave through measurement and sampling
        noise alone (noise_floor.py, mean sqrt(2/pi) sigma over the same series; the dash) and the number of the three
        contracting ones whose contraction is predicted beyond two standard errors
  d-f   example temperature curves in Mpipi-Recharged: sfAFP, alpha-synuclein 39-66 and p53(1-93), measured against the
        plain, linear and capped models, relative to 300 K
Inputs: heldout_temperature/results/primary_<model>/primary.csv (plain and A0), capped_forcing/noise_floor.csv
A force field whose tables do not exist yet is drawn as a blind placeholder.
    python fig_heldout.py
Output: BENDER-MT/manuscript/figs/fig_heldout.{png,pdf}
House style (2026-10-15, fig_extra.house): 6.5 in wide, figure text says "calibrated" and "charged-group form" (was "frozen calibration" until 2026-10-23) for
the variants A0 and S1 of the tables; PNt at 525 mM has no measured slope (dagger, legend).
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import fig_extra as FX  # noqa: E402

ROOT = '/Volumes/User Homes/Nu Project/BENDER-MT'
RES = f'{ROOT}/heldout_temperature/results'
OUT = f'{ROOT}/manuscript/figs/fig_heldout'
MODELS = [('calvados', 'CALVADOS-2', '#2E6FA3'), ('recharged', 'Mpipi-Recharged', '#6B4FA0')]
SETS = [('PNt-I025', 'PNt, 25 mM'), ('PNt-I175', 'PNt, 175 mM'), ('PNt-I525', 'PNt, 525 mM†'),
        ('asyn4W26C-I030', 'α-synuclein 4–26'), ('asyn39W66C-I030', 'α-synuclein 39–66'), ('p53N93-I118', 'p53 (1–93)'),
        ('sfAFP-I175', 'sfAFP (control)'), ('AapPGR-I195', 'Aap PGR (control)')]
CONTRACT = ['PNt-I025', 'PNt-I175', 'asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118']
NOT_FITTED = ['asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118', 'sfAFP-I175', 'AapPGR-I195']
NF = f'{ROOT}/heldout_temperature/capped_forcing/noise_floor.csv'
VARIANTS = ('plain', 'A0', 'capped')
TICK = {'plain': 'plain', 'A0': 'linear', 'capped': 'calibrated'}
XLIM = (-1.55, 0.45)                                  # every model error bar inside (PNt 175 mM CALVADOS-2 to -1.44)


def capped_525(model):
    """PNt at 525 mM (no measured slope) under the capped model, by reweighting as the other reweighted series."""
    import capped_model as CM
    import saturating_forcing as SF
    g = SF.gam_of(model, CM.params()[0]['C' if model == 'calvados' else 'R'][0])
    ser = SF.heldout_series(model, 'PNt-I525', 'S1', CM.CSAT)
    s, se, _ = SF.HA.slope(ser, SF.HA.OBS_OF['PNt-I525'], g)
    return s, se, min(SF.HA.ess(ser[T], g, float(T)) for T in ser)


def load(model, label):
    a = pd.read_csv(os.path.join(RES, f'primary_{model}', 'primary.csv'))
    a = a[a.variant.isin(['plain', 'A0'])][['variant', 'set', 's_model', 'se', 's_exp', 'min_ess']]
    c = pd.read_csv(NF)
    c = c[c.model == label].rename(columns={'series': 'set', 'se_model': 'se'})[['set', 's_model', 'se', 's_exp', 'min_ess']]
    s, se, me = capped_525(model)
    c = pd.concat([c, pd.DataFrame([{'set': 'PNt-I525', 's_model': s, 'se': se, 's_exp': np.nan, 'min_ess': me}])])
    return pd.concat([a, c.assign(variant='capped')], ignore_index=True)


def forest(ax, t, colour, title):
    ax.axvline(0, color='0.75', lw=0.7, zorder=0)
    y = np.arange(len(SETS))[::-1]
    off = {'plain': 0.22, 'A0': 0.0, 'capped': -0.22}
    for k, (st, _) in enumerate(SETS):
        e = t[t.set == st].s_exp.dropna()
        if len(e):
            ax.plot(e.iloc[0], y[k], 'D', ms=4.2, color='#222222', zorder=4)
        for v in VARIANTS:
            r = t[(t.set == st) & (t.variant == v)]
            if not len(r):
                continue
            r = r.iloc[0]
            c = '0.55' if v == 'plain' else colour
            weak = r.min_ess < 0.05
            assert XLIM[0] < r.s_model - r.se and r.s_model + r.se < XLIM[1], (st, v)
            ax.errorbar(r.s_model, y[k] + off[v], xerr=r.se, fmt='o', ms=3.3, color=c,
                        mfc='white' if v == 'A0' else c, mew=0.8, elinewidth=0.7, capsize=0, alpha=0.4 if weak else 1.0,
                        zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in SETS])
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.5, len(SETS) - 0.5)
    ax.set_xlabel('Relative slope of size (% per K)')
    ax.set_title(title, loc='left', color=colour)
    ax.set_xlim(XLIM)
    ax.spines[['top', 'right']].set_visible(False)


def blind(ax, title):
    ax.set_title(title, fontsize=8, loc='left')
    ax.text(0.5, 0.5, 'blind until all 40 runs\nare analysed', ha='center', va='center', transform=ax.transAxes,
            fontsize=8, color='0.45')
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color('0.8')


def stats(t, v):
    x = t[(t.variant == v)].set_index('set').loc[NOT_FITTED]
    mae = float((x.s_model - x.s_exp).abs().mean())
    sign = int(sum(x.loc[s, 's_model'] + 2 * x.loc[s, 'se'] < 0 for s in CONTRACT if s in NOT_FITTED))
    return mae, sign


EXAMPLES = [('sfAFP-I175', 'Rg', 'sfAFP (control)'), ('asyn39W66C-I030', 'mean_loop', 'α-synuclein 39–66'),
            ('p53N93-I118', 'Rh', 'p53 (1–93)')]


def _measured(st):
    """(T in K, relative size to 300 K in %, error or None) for the measured series; p53 only as the published fit."""
    ds = pd.read_csv(os.path.expanduser('~/Downloads/Scripts/GAIVI_HELDOUT/datasets.tsv'), sep='\t')
    if st == 'p53N93-I118':
        T = np.linspace(5, 45, 2) + 273.15
        rh = 32.8 - 0.076 * (T - 298.15)
        r300 = 32.8 - 0.076 * (300 - 298.15)
        return T, 100 * (rh / r300 - 1), None, 'DLS fit'
    d = ds[(ds.set == st) & ds.observable.str.contains('Rg|peak')].dropna(subset=['T_C'])
    T = d.T_C.to_numpy(float) + 273.15
    v = d.value_A.to_numpy(float)
    p = np.polyfit(T - 300, np.log(v), 1)
    r300 = np.exp(np.polyval(p, 0))
    lo, hi = d.low_A.to_numpy(float), d.high_A.to_numpy(float)
    err = None if np.isnan(lo).all() else [100 * (v - lo) / r300, 100 * (hi - v) / r300]
    return T, 100 * (v / r300 - 1), err, 'SAXS end points' if st.startswith('sfAFP') else 'trFRET peaks'


def example_curves(axes):
    import capped_model as CM
    import heldout_analysis as H
    import saturating_forcing as SF
    g_cap = SF.gam_of('recharged', CM.params()[0]['R'][0])
    variants = [('plain', H.frozen_gam('recharged', True), '0.55', 'o', '-', 'series'),
                ('A0', H.frozen_gam('recharged'), '#6B4FA0', 'o', '--', 'series'),
                ('capped', g_cap, '#6B4FA0', 'o', '-', 'capped')]
    for ax, (st, obs, title) in zip(axes, EXAMPLES):
        sers = {'series': H.series('recharged', st), 'capped': SF.heldout_series('recharged', st, 'S1', CM.CSAT)}
        Ts = np.array(sorted(sers['series']), float)
        for lab, gam, col, mk, ls, src in variants:
            ser = sers[src]
            v = np.array([H.ensemble_value(ser[int(T)], obs, gam, float(T)) for T in Ts])
            p = np.polyfit(Ts - 300, np.log(v), 1)
            ax.plot(Ts, 100 * (v / np.exp(np.polyval(p, 0)) - 1), ls=ls, marker=mk, ms=2.6, color=col, lw=0.9,
                    mfc='white' if lab == 'A0' else col)
        T, rel, err, kind = _measured(st)
        if err is None:
            ax.plot(T, rel, '-', color='#222222', lw=1.4)
        else:
            ax.errorbar(T, rel, yerr=err, fmt='D', ms=3.0, color='#222222', elinewidth=0.6, capsize=0, zorder=5)
        if st.startswith('sfAFP'):
            g = pd.read_csv(os.path.join(RES, 'guinier.csv'))
            g = g[g.set == st]
            gv = g.Rg_A.to_numpy(float)
            pg = np.polyfit(g.T_C + 273.15 - 300, np.log(gv), 1)
            ax.plot(g.T_C + 273.15, 100 * (gv / np.exp(np.polyval(pg, 0)) - 1), 's', ms=2.4, color='0.6', mfc='white',
                    mew=0.7, zorder=4)
        ax.axhline(0, color='0.85', lw=0.5, zorder=0)
        ax.set_title(f'{title}, Mpipi-Recharged', loc='left')
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
        ax.set_xlabel('Temperature (K)')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Size change from\n300 K (%)')


def main():
    FX.house()                                        # house style, Arial (figfont), save-time palette
    fig = plt.figure(figsize=(FX.WIDTH, 5.75))
    top = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.92], left=0.175, right=0.99, top=0.95, bottom=0.53, wspace=0.12)
    axes = [fig.add_subplot(top[0, i]) for i in range(2)]
    gc = fig.add_gridspec(1, 1, left=0.765, right=0.99, top=0.95, bottom=0.6)
    axes.append(fig.add_subplot(gc[0, 0]))
    data = {m: load(m, lab) for m, lab, _ in MODELS}
    for i, (m, lab, col) in enumerate(MODELS):
        forest(axes[i], data[m], col, lab)
    axes[1].set_yticklabels([])
    a = axes[2]
    nf = pd.read_csv(NF)
    xs, labels, out = [], [], {}
    for j, (m, lab, col) in enumerate(MODELS):
        for k, v in enumerate(VARIANTS):
            mae, sign = stats(data[m], v)
            out[(lab, v)] = (round(mae, 3), sign)
            x = j * 4 + k
            c = '0.55' if v == 'plain' else col
            a.bar(x, mae, color='white' if v == 'A0' else c, width=0.8, edgecolor=c, lw=0.8 if v == 'A0' else 0)
            a.text(x, mae + 0.008, f'{sign}/3', ha='center', va='bottom', fontsize=6.5)
            xs.append(x)
            labels.append(TICK[v])
        floor = nf[(nf.model == lab) & nf.series.isin(NOT_FITTED)].expected_noise_miss.mean()
        a.plot([j * 4 - 0.45, j * 4 + 2.45], [floor, floor], color='#222222', lw=0.8, ls=(0, (3, 1.5)), zorder=4)
        out[(lab, 'noise floor')] = round(floor, 3)
        a.text(j * 4 + 1, 0.445, lab.replace('Mpipi-', 'Mpipi-\n'), ha='center', va='top', fontsize=6.5, color=col,
               linespacing=1.0)
    print('panel c (mean |error| over the five not-fitted series, contracting predicted of 3):', out)
    a.set_xticks(xs)
    a.set_xticklabels(labels, rotation=90, fontsize=6.5)
    a.tick_params(axis='x', length=0)
    a.set_ylabel('Mean absolute error, series\nnot fitted (% per K)')
    a.set_ylim(0, 0.45)
    a.set_xlim(-0.6, 6.6)
    a.spines[['top', 'right']].set_visible(False)
    hs = [plt.Line2D([], [], marker='D', ls='none', color='#222222', ms=4.2, label='measured'),
          plt.Line2D([], [], marker='o', ls='none', color='0.55', ms=3.3, label='plain'),
          plt.Line2D([], [], marker='o', ls='none', color='0.35', mfc='white', ms=3.3, label='linear form'),
          plt.Line2D([], [], marker='o', ls='none', color='0.35', ms=3.3, label='calibrated model'),
          plt.Line2D([], [], marker='o', ls='none', color='0.35', alpha=0.4, ms=3.3, label='effective sample < 0.05'),
          plt.Line2D([], [], color='#222222', lw=0.8, ls=(0, (3, 1.5)), label='noise floor'),
          plt.Line2D([], [], marker='s', ls='none', color='0.6', mfc='white', mew=0.7, ms=2.6, label='Guinier fits (sfAFP)'),
          plt.Line2D([], [], ls='none', label='† no measured slope')]
    fig.legend(handles=hs, loc='center', ncol=4, bbox_to_anchor=(0.53, 0.395), columnspacing=1.4, handletextpad=0.3)
    g2 = fig.add_gridspec(1, 3, left=0.085, right=0.98, top=0.29, bottom=0.075, wspace=0.3)
    ex = [fig.add_subplot(g2[0, i]) for i in range(3)]
    example_curves(ex)
    FX.add_letters(fig, [(axes[0], 'a', {'row': 1}), (axes[1], 'b', {'row': 1, 'dx': -10}), (axes[2], 'c', {'row': 1})] +
                   [(a_, 'def'[k], {'row': 2}) for k, a_ in enumerate(ex)])
    FX.save(fig, OUT)
    print('written', OUT)


if __name__ == '__main__':
    main()
