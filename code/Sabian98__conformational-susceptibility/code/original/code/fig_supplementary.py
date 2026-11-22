"""
fig_supplementary.py — Supplementary Figs S9-S14 of the manuscript (fig_si1 ... fig_si6), from stored outputs only
(no new simulations).
  si1  validation at every temperature: ratio of the mean local slope (finite difference between neighbouring
       temperatures) to the mean susceptibility for the eight observables, CALVADOS-2 and Mpipi-GG; per-chain Pearson r
       of local slopes and of the 278-320 K integral; mean chi_T(Rg) and mean local slope against T
       (Analysis_FDT/fdt_validation_summary.csv, Analysis_FDT_MPIPI/mpipi_fdt_validation_summary.csv)
  si2  per-protein agreement between force fields: fixed-solvent chi_T(Rg) and total response, CALVADOS-2 against
       Mpipi-Recharged and Mpipi-GG, every chain inside the axes (total_response_evidence.load)
  si3  sequence predictability of the total response and of its residual beyond compaction: cumulative R2 by descriptor
       layer and single-group R2, with the noise ceiling (total_response/F_predictability.csv)
  si4  PAGE4, every directly simulated form: calibrated heating curves simulated directly in both models; reweighted
        against direct change
       of Rg and of the heating response, with the Kish effective sample of the reweighting (page4/, page4_direct)
  si5  PAGE4 at the measurement conditions (295 K, smFRET and SAXS buffers) against the measurements
  si6  mutation landscapes: median change of the total response by residue replaced and introduced, CALVADOS-2 and
       Mpipi-Recharged (reliable mutants, effective sample above 0.1), and every substitution of one example chain
       (total_response/mutation_scan/); undefined or unreliable cells grey
House style (2026-10-15): 6.5 in wide, Arial (figfont), smallest text 6.5 pt, axis labels 7.5 pt, bold 9 pt panel
letters placed left of each panel's y-axis label, force-field colours only for force fields; self-contained (does not
import fig_manuscript2). The helpers (style, letters, save, ff_cmap) are shared with fig_supp_redraw.py.
    BENDER_REANALYSIS_DIR=... python -B fig_supplementary.py [si1 ... si6]      (openmpipi env)
Output: BENDER-MT/manuscript/figs/fig_si{1..6}.{png,pdf}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import figfont  # noqa: E402

ROOT = '/Volumes/User Homes/Nu Project/BENDER-MT'
TLA = f'{ROOT}/Traj Level Analysis'
TR = f'{TLA}/Analysis_FDT/total_response'
OUT = f'{ROOT}/manuscript/figs'
P4 = f'{TR}/page4'
COL = {'CALVADOS-2': '#2E6FA3', 'Mpipi-GG': '#A9502A', 'Mpipi-Recharged': '#6B4FA0'}   # remapped by figfont
MODELS = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
W = 6.5                                                               # print width (in)
RG, REE = r'$R_\mathrm{g}$', r'$R_\mathrm{ee}$'
NU = r'$\mathtt{\nu}$'   # Greek nu from STIX italic (figfont mathtext.tt); Arial's nu reads as v
OBS = [('rg', RG), ('ree', REE), ('delta', 'Asphericity'), ('fragmentation_index', 'Fragmentation'),
       ('global_efficiency', 'Global efficiency'), ('avg_clustering', 'Clustering'), ('transitivity', 'Transitivity'),
       ('degree_assortativity', 'Assortativity')]
TS = ['278', '288', '300', '310', '320']
AA_ORDER = 'DEKRHNQSTGPAVLIMCFWY'
BIN_LAB = ['0.30–0.45', '0.45–0.50', '0.50–0.54', '0.54–0.58', '0.58–0.63', '>0.63']
GREY_BAD = '#B3B3B3'                                                  # undefined / unreliable cells

STYLE = {'font.size': 7, 'axes.labelsize': 7.5, 'axes.titlesize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
         'legend.fontsize': 6.5, 'legend.frameon': False, 'axes.spines.top': False, 'axes.spines.right': False,
         'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.5,
         'ytick.minor.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5,
         'ytick.minor.size': 1.5, 'lines.linewidth': 1.0, 'lines.markersize': 3.5, 'axes.unicode_minus': True,
         'axes.prop_cycle': matplotlib.cycler(color=['#222222']), 'figure.dpi': 100, 'savefig.dpi': 300,
         'axes.titlepad': 4, 'axes.labelpad': 2.5, 'legend.handlelength': 1.6, 'legend.borderaxespad': 0.3}
plt.rcParams.update(STYLE)
figfont.apply()                                                       # Arial, fonttype 42, okabe-ito remapping


def active(c):
    """Colour as printed: the figfont palette counterpart of an original hex colour."""
    role = figfont.ROLE.get(mcolors.to_hex(c).upper())
    return figfont.PALETTES[figfont._active].get(role, c) if role else c


def line_shade(c):
    """The palette's deeper shade for thin marks of a role ('<role>_line', state-light Mpipi-GG), else the colour itself."""
    role = figfont.ROLE.get(mcolors.to_hex(c).upper())
    pal = figfont.PALETTES[figfont._active]
    return pal.get(f'{role}_line', active(c)) if role else c


def mix(c, other, t):
    a, b = np.array(mcolors.to_rgb(c)), np.array(mcolors.to_rgb(other))
    return tuple((1 - t) * a + t * b)


def ff_cmap(model):
    """Single-hue colormap of a force field's printed colour (pale tint to dark shade)."""
    c, cl = active(COL[model]), line_shade(COL[model])
    dark = cl if cl != c else mix(c, 'black', 0.45)                   # a pale fill colour runs to its line shade
    return mcolors.LinearSegmentedColormap.from_list(f'ff_{model}', [mix(c, 'white', 0.82), c, dark])


def letters(fig, axes, names, dy=3.0, dx=0.0):
    """Bold 9 pt panel letters, top-left outside each axes: left edge at the left of the y-axis labels, bottom 3 pt
    above the axes (or level with a title)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pt = fig.dpi / 72
    for ax, s in zip(axes, names):
        yb = ax.yaxis.get_tightbbox(r)
        x0 = min(ax.bbox.x0, yb.x0) if yb is not None else ax.bbox.x0
        y0 = ax.bbox.y1 + dy * pt
        if ax.title.get_text():
            y0 = max(y0, ax.title.get_window_extent(r).y0)
        fig.text((x0 + dx * pt) / fig.bbox.width, y0 / fig.bbox.height, s, fontsize=9, fontweight='bold',
                 ha='left', va='bottom')


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(f'{OUT}/{name}.pdf')
    fig.savefig(f'{OUT}/{name}.png', dpi=300)
    plt.close(fig)
    print('saved', name, f'{fig.get_figwidth():.2f} x {fig.get_figheight():.2f} in')


def si1():
    vc = pd.read_csv(f'{TLA}/Analysis_FDT/fdt_validation_summary.csv')
    vm = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/mpipi_fdt_validation_summary.csv')
    fig, ax = plt.subplots(2, 2, figsize=(W, 4.7))
    ax = ax.ravel()
    cmap = plt.get_cmap('viridis')
    ocol = {k: cmap(0.9 * i / 7) for i, (k, _) in enumerate(OBS)}   # ends at lime, not pale yellow
    allr = []
    for v in (vc, vm):
        for key, _ in OBS:
            g = v[(v.metric == key) & v['T'].isin(TS)].set_index('T').loc[TS]
            allr += list(g.mean_slope / g.mean_chi)
    ylo, yhi = min(0, min(allr) - 0.1), max(allr) + 0.15
    out = {}
    for a, (v, m) in zip(ax[:2], ((vc, 'CALVADOS-2'), (vm, 'Mpipi-GG'))):
        n = int(v[(v.metric == 'rg') & (v['T'] == '300')].n.iloc[0])
        out[m] = n
        for key, name in OBS:
            g = v[(v.metric == key) & v['T'].isin(TS)].set_index('T').loc[TS]
            a.plot([int(t) for t in TS], g.mean_slope / g.mean_chi, 'o-', ms=2.8, lw=0.9, color=ocol[key], label=name)
        a.axhline(1, color='0.45', lw=0.6, ls='--', zorder=0)
        a.set(ylim=(ylo, yhi), xlabel='Temperature (K)', ylabel=r'Mean local slope / mean $\chi_T$',
              xticks=[278, 288, 300, 310, 320], xlim=(275, 323))
        a.set_title(f'{m} ({n:,} chains)', color=COL[m])
    ax[0].legend(ncol=2, loc='lower left', handlelength=1.4, columnspacing=1.0, borderaxespad=0.2)
    a = ax[2]
    xs = np.r_[np.arange(5), 5.35]
    for v, m, off in ((vc, 'CALVADOS-2', -0.17), (vm, 'Mpipi-GG', 0.17)):
        g = v[v.metric == 'rg'].set_index('T').loc[TS + ['integral']]
        a.bar(xs + off, g.pearson_r, 0.32, color=COL[m], label=m)
        out[f'r_{m}'] = g.pearson_r.round(2).tolist()
    a.set(xticks=xs, ylim=(0, 1.05), xlim=(-0.6, 5.95), ylabel=f'Per-chain Pearson $r$, {RG}', xlabel='Temperature (K)')
    a.set_xticklabels(TS + ['278–320\nintegral'])
    a.legend(loc='upper left', handlelength=1.0)
    a = ax[3]
    for v, m in ((vc, 'CALVADOS-2'), (vm, 'Mpipi-GG')):
        g = v[v.metric == 'rg'].set_index('T').loc[TS]
        T = [int(t) for t in TS]
        a.plot(T, 1e3 * g.mean_chi, 'o-', ms=3.2, lw=1, color=COL[m])
        a.plot(T, 1e3 * g.mean_slope, 's--', ms=3.2, lw=0.9, color=COL[m], mfc='white')
    a.legend(handles=[Line2D([], [], color='0.25', marker='o', ms=3.2, label=r'mean $\chi_T$'),
                      Line2D([], [], color='0.25', marker='s', mfc='white', ls='--', ms=3.2, label='mean local slope')],
             loc='upper right')
    a.set(xlabel='Temperature (K)', ylabel=f'd{RG}/d$T$ (10$^{{-3}}$ nm per K)', xticks=[278, 288, 300, 310, 320],
          xlim=(275, 323), ylim=(0.5, 3.45))
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=2.2, w_pad=2.0)
    letters(fig, ax, 'abcd')
    save(fig, 'fig_si1')
    return out


def si2():
    import total_response_evidence as E
    D, _ = E.load()
    fig, ax = plt.subplots(2, 2, figsize=(W, 5.6))
    ax = ax.ravel()
    pairs = [('chi', 'R'), ('chi', 'M'), ('total', 'R'), ('total', 'M')]
    head = {'chi': f'Fixed-solvent susceptibility $\\chi_T$({RG})', 'total': 'Total response, calibrated hydration'}
    out = {}
    for a, (col, other) in zip(ax, pairs):
        mo = MODELS[other]
        j = D['C'][[col]].join(D[other][[col]], rsuffix='_o', how='inner').dropna()
        x, y = j[col].to_numpy(), j[f'{col}_o'].to_numpy()
        px, py = 0.04 * np.ptp(x), 0.04 * np.ptp(y)
        ex = (x.min() - px, x.max() + px, y.min() - py, y.max() + py)
        sx, sy = np.ptp(x), np.ptp(y)
        kde = stats.gaussian_kde(np.vstack([x / sx, y / sy]), bw_method=0.15 if len(x) > 3000 else 0.24)
        dens = kde(np.vstack([x / sx, y / sy]))
        o = np.argsort(dens)                                          # densest points drawn last
        c, cl = active(COL[mo]), line_shade(COL[mo])
        dark = cl if cl != c else mix(c, '#2E2E2E', 0.3)
        start = mix(cl, 'white', 0.75) if cl != c else mix(c, 'white', 0.8)   # pale roles: tint of the line shade stays visible
        soft = mcolors.LinearSegmentedColormap.from_list(f'soft_{other}', [start, c, dark])
        a.scatter(x[o], y[o], c=stats.rankdata(dens[o]) / len(x), cmap=soft, vmin=0, vmax=1.15,
                  s=1.6 if len(x) > 3000 else 4.0, lw=0, alpha=0.9, rasterized=True, zorder=2)
        gx, gy = np.meshgrid(np.linspace(*ex[:2], 160), np.linspace(*ex[2:], 160))
        gz = kde(np.vstack([gx.ravel() / sx, gy.ravel() / sy])).reshape(gx.shape)
        a.contour(gx, gy, gz, levels=np.sort(np.quantile(dens, [0.1, 0.5])), colors=[cl if cl != c else mix(c, '#2E2E2E', 0.45)],
                  linewidths=[0.45, 0.65], zorder=3)                  # regions holding 90% and 50% of chains
        lo, hi = min(ex[0], ex[2]), max(ex[1], ex[3])
        a.plot([lo, hi], [lo, hi], color='0.35', lw=0.6, ls=(0, (4, 3)), zorder=3)
        a.set(xlim=ex[:2], ylim=ex[2:])
        rho = stats.spearmanr(x, y)[0]
        out[f'{col} C vs {other}'] = (round(rho, 2), len(j))
        a.text(0.04, 0.96, f'$\\rho$ = {rho:.2f}\nn = {len(j):,}', transform=a.transAxes, va='top', fontsize=7)
        a.set(xlabel='CALVADOS-2 (% per K)', ylabel=f'{mo} (% per K)')
        a.xaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 5, 10]))
        a.set_title(head[col])
        assert ((x >= ex[0]) & (x <= ex[1]) & (y >= ex[2]) & (y <= ex[3])).all()
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=2.2, w_pad=2.5)
    letters(fig, ax, 'abcd')
    save(fig, 'fig_si2')
    return out


def si3():
    f = pd.read_csv(f'{TR}/F_predictability.csv').set_index('target')
    rows = [('CALVADOS-2 total (own full set)', 'CALVADOS-2, total', COL['CALVADOS-2'], '-', 1.0),
            ('Mpipi-Recharged total (own full set)', 'Mpipi-Recharged, total', COL['Mpipi-Recharged'], '-', 1.0),
            ('CALVADOS-2 leftover (own full set)', 'CALVADOS-2, beyond compaction', COL['CALVADOS-2'], '--', 0.5),
            ('Mpipi-Recharged leftover (own full set)', 'Mpipi-Recharged, beyond compaction', COL['Mpipi-Recharged'],
             '--', 0.5),
            ('shared leftover (three models)', 'three models, beyond compaction', '#333333', ':', 0.35)]
    layers = ['R2 L1 length (cum.)', 'R2 L2 composition (cum.)', 'R2 L3 patterning (cum.)', 'R2 L4 Uversky (cum.)',
              'R2 L5 pi (cum.)', 'R2 quadratic all']
    lab = ['length', '+ compo-\nsition', '+ pattern-\ning', '+ charge,\nhydropathy', '+ $\\pi$\ncontent', '+ quad-\nratic']
    fig, ax = plt.subplots(1, 2, figsize=(W, 3.15), gridspec_kw={'width_ratios': [1.3, 1]})
    a = ax[0]
    for key, name, c, ls, al in rows:
        r = f.loc[key]
        a.plot(range(len(layers)), [max(0, r[x]) for x in layers], marker='o', ms=3, ls=ls, lw=1.1, color=c)
        a.plot([len(layers) - 0.55, len(layers) - 0.15], [r.noise_ceiling] * 2, ls=ls, color=c, lw=1.4)
    a.set(xticks=range(len(layers)), ylim=(0, 1.13), ylabel='Cross-validated $R^2$', xlim=(-0.3, len(layers) - 0.05),
          xlabel='Descriptor groups added in turn')
    a.set_xticklabels(lab, fontsize=6.5)
    a.text(len(layers) - 0.35, 1.12, 'noise\nceiling', fontsize=6.5, ha='center', va='top', color='0.3')
    a = ax[1]
    groups = ['R2 stickers only', 'R2 charge only', 'R2 hydropathy only', 'R2 hydrophilic residues only']
    glab = ['stickers', 'charge', 'hydropathy', 'hydrophilic\nresidues']
    x = np.arange(len(groups))
    for o, (key, name, c, ls, al) in zip((-0.32, -0.16, 0, 0.16, 0.32), rows):
        r = f.loc[key]
        a.bar(x + o, [max(0, r[g]) for g in groups], 0.15, color=c, alpha=al, linewidth=0)
    a.set(xticks=x, ylim=(0, 1.0), ylabel='$R^2$, one descriptor group alone')
    a.set_xticklabels(glab, fontsize=6.5)
    handles = [(Line2D([], [], color=c, ls=ls, lw=1.1), Patch(facecolor=c, alpha=al, linewidth=0))
               for _, _, c, ls, al in rows]
    fig.legend(handles, [r[1] for r in rows], loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=3,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.5)}, handlelength=4.2, columnspacing=1.6)
    fig.tight_layout(rect=(0, 0.13, 1, 0.975), w_pad=2.5)
    letters(fig, ax, 'ab')
    save(fig, 'fig_si3')
    return {r[1]: round(f.loc[r[0], 'noise_ceiling'], 3) for r in rows}


def _page4_tables():
    import page4_direct as PD
    rw = pd.read_csv(f'{P4}/phosphoforms.csv')
    out = {}
    for model, m in PD.MODELS.items():
        runs = pd.read_csv(f'{P4}/direct_{model}/runs.csv')
        # capped main model (page4_capped_fdt.py, 2026-10-25): calibrated total and heating curves of every form; the
        # plain curves and every other column from page4_direct
        cap = f'{ROOT}/heldout_temperature/capped_forcing/page4'
        cc = pd.read_csv(f'{cap}/curves_capped.csv').replace({'curve': {'capped': 'calibrated'}})
        cv = cc[cc.model == m['label']].drop(columns='model')
        # calibrated curves simulated directly with the forcing (page4cal_analysis.py; Recharged 2026-10-25,
        # CALVADOS-2 2026-11-02) replace the reweighted ones
        dc = pd.read_csv(f'{cap}/direct_cal_{model}/curves.csv')
        dc = dc[dc.obs == 'Rg'].rename(columns={f'se_{t}K': f'rel_se_{t}K' for t in (278, 288, 300, 310, 320)})
        dc = dc.assign(curve='calibrated', temps='278,288,300,310,320', min_ess=1.0)[
            ['form', 'curve', 'temps', 'slope300', 'slope300_se', 'min_ess'] + [c for c in dc.columns if c.startswith('rel_')]]
        cv = pd.concat([cv[cv.curve != 'calibrated'], dc], ignore_index=True)
        h = runs[(runs['T'] == 300) & (runs.ionic == 0.15)].set_index('form')
        fd = pd.read_csv(f'{cap}/fdt_300K.csv')
        fd = fd[fd.model == m['label']].set_index('form')
        assert set(h.index) <= set(fd.index)
        h['Rg_total'], h['Rg_total_se'] = fd.loc[h.index, 'total_capped'], fd.loc[h.index, 'total_capped_se']
        r = rw[rw.model == m['label']].set_index('pattern')
        out[m['label']] = (runs, cv, h, r, m['pred'])
    return out


# PAGE4 forms: unmodified black, HIPK1 forms green, CLK2 forms brown, S7 and S9 grey; Glu mimics dashed, open markers
P4_FORMS = ['WT', 'HIPK1occ', 'HIPK1T51', 'HIPK1glu', 'CLK2occ', 'CLK2full', 'CLK2glu', 'Nmotif']
P4_NAME = {'WT': 'unmodified', 'HIPK1occ': 'HIPK1', 'HIPK1T51': 'HIPK1, T51 only', 'HIPK1glu': 'HIPK1, Glu mimic',
           'CLK2occ': 'CLK2', 'CLK2full': 'CLK2, all sites full', 'CLK2glu': 'CLK2, Glu mimic', 'Nmotif': 'S7 and S9'}
P4_STYLE = {'WT': ('#000000', '-', 'o', True), 'HIPK1occ': ('#007A5A', '-', 'o', True),
            'HIPK1T51': ('#5CC29F', '-', '^', True), 'HIPK1glu': ('#007A5A', '--', 'o', False),
            'CLK2occ': ('#8C510A', '-', 'o', True), 'CLK2full': ('#D2A04E', '-', '^', True),
            'CLK2glu': ('#8C510A', '--', 'o', False), 'Nmotif': ('#8A8A8A', '-', 'D', True)}


def si4():
    T = _page4_tables()
    fw, fh = W, 5.7
    fig = plt.figure(figsize=(fw, fh))

    def box(x, y, w, h):                                              # inches -> figure fractions
        return [x / fw, y / fh, w / fw, h / fh]
    ax = [fig.add_axes(box(0.62, 3.62, 2.3, 1.75)), fig.add_axes(box(3.55, 3.62, 2.3, 1.75)),
          fig.add_axes(box(0.62, 0.58, 2.3, 1.9)), fig.add_axes(box(3.55, 0.58, 2.3, 1.9))]
    present, vals = {}, []
    for a, m in zip(ax[:2], ('CALVADOS-2', 'Mpipi-Recharged')):
        _, cv, _, _, _ = T[m]
        c = cv[cv.curve == 'calibrated'].set_index('form')
        present[m] = [f for f in P4_FORMS if f in c.index]
        for f in present[m]:
            col, ls, mk, filled = P4_STYLE[f]
            Ts = [int(t) for t in c.loc[f, 'temps'].split(',')]
            y = [c.loc[f, f'rel_{t}K'] for t in Ts]
            vals += y
            a.plot(Ts, y, marker=mk, ls=ls, ms=3, lw=0.9, color=col, mfc=col if filled else 'white', mew=0.8)
        a.set(xlabel='Temperature (K)', ylabel=f'Mean {RG} relative to 300 K', xticks=[278, 288, 300, 310, 320],
              xlim=(275, 323))
        a.set_title(f'{m} ({len(present[m])} forms)', color=COL[m])
    ylo, yhi = min(vals) - 0.012, max(vals) + 0.012
    for a in ax[:2]:
        a.set_ylim(ylo, yhi)
    # forms legend between the rows; forms simulated in Mpipi-Recharged only are marked
    hs = []
    for f in P4_FORMS:
        col, ls, mk, filled = P4_STYLE[f]
        only = f not in present['CALVADOS-2']
        hs.append(Line2D([], [], color=col, ls=ls, marker=mk, ms=3, lw=0.9, mfc=col if filled else 'white', mew=0.8,
                         label=P4_NAME[f] + (' (b only)' if only else '')))
    fig.legend(handles=hs, loc='center', bbox_to_anchor=(0.5, 2.88 / fh), ncol=4, handlelength=2.4,
               columnspacing=1.5, labelspacing=0.5)
    ess_cmap = mcolors.LinearSegmentedColormap.from_list('ess_grey', ['#F2F2F2', '#111111'])
    sc = None
    for k, (a, (rcol, rse, lab)) in enumerate(zip(ax[2:], (('Rg_pct', 'Rg_pct_se', f'change of {RG} (%)'),
                                                          ('Rg_dtotal', 'Rg_dtotal_se',
                                                           'change of total\nresponse (% per K)')))):
        allv = []
        for m, mk in (('CALVADOS-2', 'o'), ('Mpipi-Recharged', 's')):
            runs, cv, h, r, pred = T[m]
            wt = h.loc['WT']
            for f, pat in pred.items():
                if f not in h.index:
                    continue
                z = h.loc[f]
                if k == 0:
                    d = 100 * (z.Rg_mean / wt.Rg_mean - 1)
                    dse = 100 * z.Rg_mean / wt.Rg_mean * np.hypot(z.Rg_se / z.Rg_mean, wt.Rg_se / wt.Rg_mean)
                else:
                    d, dse = z.Rg_total - wt.Rg_total, np.hypot(z.Rg_total_se, wt.Rg_total_se)
                p = r.loc[pat]
                sc = a.scatter(d, p[rcol], c=[[p.ess]], cmap=ess_cmap, vmin=0, vmax=0.8, marker=mk, s=24,
                               edgecolor=COL[m], linewidth=1.0, zorder=3)
                a.errorbar(d, p[rcol], xerr=dse, yerr=p[rse], fmt='none', ecolor='0.55', elinewidth=0.6, zorder=2)
                allv += [d - dse, d + dse, p[rcol] - p[rse], p[rcol] + p[rse]]
        span = max(allv) - min(allv)
        lo, hi = min(allv) - 0.06 * span, max(allv) + 0.06 * span
        a.plot([lo, hi], [lo, hi], '--', color='0.5', lw=0.7, zorder=1)
        a.set(xlim=(lo, hi), ylim=(lo, hi), xlabel=f'Direct {lab}', ylabel=f'Reweighted {lab}')
    ax[2].legend(handles=[Line2D([], [], marker='o', ls='', ms=4.5, mfc='white', mec=COL['CALVADOS-2'], mew=1.0,
                                label='CALVADOS-2'),
                          Line2D([], [], marker='s', ls='', ms=4.5, mfc='white', mec=COL['Mpipi-Recharged'], mew=1.0,
                                 label='Mpipi-Recharged')], loc='lower right', handletextpad=0.3)
    cax = fig.add_axes(box(5.97, 0.58, 0.08, 1.9))
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label('Effective sample fraction')
    cb.outline.set_linewidth(0.5)
    letters(fig, ax, 'abcd')
    save(fig, 'fig_si4')
    return {m: present[m] for m in present}


def si5():
    T = _page4_tables()
    meas = {'E18_63': [0.55, 0.52, 0.35], 'E63_102': [0.65, 0.65, 0.58], 'Rg': [3.62, 3.47, 4.98]}   # Rg in nm
    meas_se = {'Rg': [0.11, 0.12, 0.19]}
    forms = ['WT', 'HIPK1occ', 'CLK2occ']
    fig, ax = plt.subplots(1, 3, figsize=(W, 2.45))
    x = np.arange(3)
    out = {}
    for a, (k, ionic, lab, ttl) in zip(ax, (('E18_63', 0.128, 'FRET efficiency, residues 18 and 63',
                                             'smFRET buffer, 295 K'),
                                            ('E63_102', 0.128, 'FRET efficiency, residues 63 and 102',
                                             'smFRET buffer, 295 K'),
                                            ('Rg', 0.19, f'{RG} (nm)', 'SAXS buffer, 295 K'))):
        a.errorbar(x, meas[k], yerr=meas_se.get(k), fmt='D-', color='#222222', ms=4, lw=1, elinewidth=0.8,
                   label='measured')
        for m, off in (('CALVADOS-2', -0.14), ('Mpipi-Recharged', 0.14)):
            runs = T[m][0]
            g = runs[(runs['T'] == 295) & (np.isclose(runs.ionic, ionic))].set_index('form').reindex(forms)
            a.errorbar(x + off, g[f'{k}_mean'], yerr=g[f'{k}_se'], fmt='o', ms=4, color=COL[m], elinewidth=0.8,
                       label=m)
            out[f'{k} {m}'] = g[f'{k}_mean'].round(3).tolist()
        a.set(xticks=x, xlim=(-0.45, 2.45), ylabel=lab)
        a.set_xticklabels(['unmodified', 'HIPK1', 'CLK2'])
        a.set_title(ttl)
    ax[0].legend(loc='lower left')
    fig.tight_layout(rect=(0, 0, 1, 0.97), w_pad=2.2)
    letters(fig, ax, 'abc')
    save(fig, 'fig_si5')
    return out


def si6():
    md = f'{TR}/mutation_scan'
    sets = {'CALVADOS-2': pd.read_csv(f'{md}/mutants.csv'), 'Mpipi-Recharged': pd.read_csv(f'{md}/mutants_recharged.csv')}
    fig = plt.figure(figsize=(W, 6.0))
    cm = plt.get_cmap('RdBu').copy()
    cm.set_bad(GREY_BAD)
    vmax = 0.06
    # inches -> figure fractions
    fw, fh = fig.get_figwidth(), fig.get_figheight()
    side, top_y = 2.2, 3.45                                               # heatmap size and bottom edge (in)
    axes = []
    out = {}
    for i, (m, d) in enumerate(sets.items()):
        a = fig.add_axes([(0.42 + i * 2.95) / fw, top_y / fh, side / fw, side / fh])
        d = d[d.ess > 0.1]
        mat = d.groupby(['wt', 'mut']).dtotal.median().unstack().reindex(index=list(AA_ORDER), columns=list(AA_ORDER))
        im = a.imshow(np.ma.masked_invalid(mat.values), cmap=cm, vmin=-vmax, vmax=vmax, aspect='equal',
                      interpolation='nearest')
        a.set_xticks(range(20))
        a.set_xticklabels(list(AA_ORDER), fontsize=6.5)
        a.set_yticks(range(20))
        a.set_yticklabels(list(AA_ORDER), fontsize=6.5)
        a.tick_params(length=1.5, pad=1.5)
        a.set(xlabel='Residue introduced', ylabel='Residue replaced')
        a.set_title(f'{m} ({len(d):,} reliable mutants)', color=COL[m])
        out[m] = len(d)
        axes.append(a)
    cax = fig.add_axes([(0.42 + 2.95 + side + 0.12) / fw, top_y / fh, 0.08 / fw, side / fh])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Median change of total\nresponse (% per K)')
    cb.outline.set_linewidth(0.5)
    # c example chain
    d = sets['CALVADOS-2']
    s = pd.read_csv(f'{md}/sample.csv')
    rand = s[~s.arg_shortlist].UniProt_ID
    L = d[d.UniProt_ID.isin(rand)].groupby('UniProt_ID').pos.max()
    uid = (L - 90).abs().idxmin()
    e = d[d.UniProt_ID == uid].copy()
    n_unrel = int((e.ess <= 0.1).sum())
    e.loc[e.ess <= 0.1, 'dtotal'] = np.nan
    mat = e.pivot_table(index='mut', columns='pos', values='dtotal', dropna=False).reindex(list(AA_ORDER))
    seq = e.drop_duplicates('pos').set_index('pos').wt.sort_index()
    mat = mat.reindex(columns=seq.index)
    a = fig.add_axes([0.42 / fw, 0.45 / fh, 5.2 / fw, 2.2 / fh])
    vm = np.nanpercentile(np.abs(mat.values), 98)
    im2 = a.imshow(np.ma.masked_invalid(mat.values), cmap=cm, vmin=-vm, vmax=vm, aspect='auto',
                   interpolation='nearest')
    a.set_yticks(range(20))
    a.set_yticklabels(list(AA_ORDER), fontsize=6.5)
    step = 5
    a.set_xticks(range(0, len(seq), step))
    a.set_xticklabels([f'{seq.iloc[k]}{seq.index[k]}' for k in range(0, len(seq), step)], fontsize=6.5)
    a.tick_params(length=1.5, pad=1.5)
    a.set(ylabel='Residue introduced', xlabel='Wild-type residue and position')
    a.set_title(f'CALVADOS-2, every substitution of {uid} ({len(seq)} residues)', loc='left')
    a.legend(handles=[Patch(facecolor=GREY_BAD, label='wild type, or effective sample ≤ 0.1')],
             loc='lower right', bbox_to_anchor=(1.0, 1.0), borderaxespad=0.15, handlelength=1.0, handleheight=0.9)
    p = a.get_position()
    cax2 = fig.add_axes([p.x1 + 0.12 / fw, p.y0, 0.08 / fw, p.height])
    cb2 = fig.colorbar(im2, cax=cax2, ticks=[-0.1, -0.05, 0, 0.05, 0.1])
    cb2.set_label('Change of total\nresponse (% per K)')
    cb2.outline.set_linewidth(0.5)
    axes.append(a)
    letters(fig, axes, 'abc')
    save(fig, 'fig_si6')
    out.update({'example chain': uid, 'length': len(seq), 'unreliable cells': n_unrel,
                'colour limit c': round(float(vm), 3)})
    return out


if __name__ == '__main__':
    which = sys.argv[1:] or ['si1', 'si2', 'si3', 'si4', 'si5', 'si6']
    for w in which:
        print(w, globals()[w](), flush=True)
