"""
fig_extra.py — additional manuscript panels (IDRome / ALBATROSS-style schematic, worked examples, maps and sequence
tracks), used by fig_manuscript2.py, fig_page4.py, fig_heldout.py and fig_arial_deckfigs.py. Every panel reads stored
results only (no new simulations).

Fig. 1  schematic(ax)                 the approach, with real conformations and the real Rg-energy covariance of one chain
        example_dist(ax), example_curve(ax)
                                      two named chains: a compact Gly-rich protein (tomato, Q02021) and protamine P1
                                      (marsupial mole, P42143); Rg distributions at 278 and 320 K, and <Rg>(T) against
                                      the tangent given by the fluctuation relation at 300 K (CALVADOS-2)
Fig. 2  fig2(out)                     Wuttke curves (a-e), per-residue forcing (f), left-out slopes (g), residue sets (h)
Fig. 3  fig3(out)                     totals by compaction (a), agreement (b), nu x stickiness map (c), composition (d),
                                      kingdoms (e), residue enrichment in the 5% most/least heat-sensitive chains (f),
                                      CALVADOS-2 against Mpipi-Recharged (g)
Fig. 5  fig5_tracks(out)              sequence tracks of HMGA1 and myelin basic protein (c, d; also drawn inside
                                      fig_apps.py, the one-file Fig. 5)
Fig. S15 si7_s1(out)                  proteome with the frozen calibration against the charged-group form
House style for Figs 2-6 and S15 (2026-10-15): house(), add_letters(), ff_cmap(), save(); 6.5 in wide (the LaTeX text
width), Arial via figfont, smallest text 6.5 pt, axis labels 7.5 pt, bold 9-pt panel letters, every data point inside
the axes, PDF (TrueType) and 300-dpi PNG.
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.expanduser('~/Downloads/Scripts/bender_reanalysis')))
sys.path.insert(0, HERE)
ROOT = '/Volumes/User Homes/Nu Project/BENDER-MT'
TLA = f'{ROOT}/Traj Level Analysis'
FDT = f'{TLA}/Analysis_FDT'
NPZ = f'{FDT}/npz'
TEMPS = (278, 288, 300, 310, 320)
INK = '#222222'
COL = {'CALVADOS-2': '#2E6FA3', 'Mpipi-GG': '#A9502A', 'Mpipi-Recharged': '#6B4FA0'}   # mapped by figfont at save
EX = {'Q02021': {'name': 'Gly-rich protein', 'org': 'tomato', 'col': '#C8792A', 'dir': f'{ROOT}/BENDER_300K/Plants/Q02021'},
      'P42143': {'name': 'Protamine P1', 'org': 'marsupial mole', 'col': '#2E8C6F', 'dir': f'{ROOT}/BENDER_300K/Mammals/P42143'}}


def block_se(a, n=10):
    b = [x.mean() for x in np.array_split(np.asarray(a, float), n)]
    return float(np.std(b, ddof=1) / np.sqrt(n))


def _chi300(uid):
    r = pd.read_csv(f'{FDT}/fdt_results_300K.csv', usecols=['UniProt_ID', 'metric', 'chi_total', 'chi_total_se'])
    r = r[(r.UniProt_ID == uid) & (r.metric == 'rg')].iloc[0]
    return float(r.chi_total), float(r.chi_total_se)


def _rg(uid, T):
    return np.load(f'{NPZ}/{T}K/{uid}.npz')['rg']


def _conformations(uid, n=3, seed=3):
    """n conformations from the 300 K trajectory, projected on their own principal plane (nm), centred."""
    import bender_common as bc
    d = EX[uid]['dir']
    xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=200)
    rng = np.random.default_rng(seed)
    out = []
    for f in rng.choice(len(xyz), n, replace=False):
        x = xyz[f] - xyz[f].mean(0)
        u, s, vt = np.linalg.svd(x, full_matrices=False)
        out.append(x @ vt[:2].T)
    return out


def _covariance(uid):
    z = np.load(f'{NPZ}/300K/{uid}.npz')
    H = z['V_AH'] + z['V_DH'] + z['V_bond']
    return z['rg'], H


# ----------------------------------------------------------------------------- Fig. 1a schematic
def _box(ax, x, y, w, h, fc, ec='none', r=1.2, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw,
                                transform=ax.transData, zorder=1))


def _arrow(ax, x0, y0, x1, y1, col='0.35'):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=9, lw=1.0, color=col, zorder=3,
                                 shrinkA=0, shrinkB=0))


def schematic(ax):
    """Proteome -> three force fields -> one 300 K ensemble -> fluctuation relation -> responses -> tests."""
    ax.set_xlim(0, 100)
    ax.set_ylim(3.2, 26.8)
    ax.axis('off')
    tint = '#F1F3F5'
    # 1 proteome
    _box(ax, 0.3, 4, 15.5, 22, tint)
    ax.text(8.05, 23.6, 'Disordered\nproteomes', ha='center', va='top', fontsize=7.2, fontweight='bold')
    ax.text(8.05, 16.2, '10,282 proteins\n2,156 organisms\nbacteria, plants, fungi,\nmammals, protists,\nviruses', ha='center',
            va='top', fontsize=6.1, linespacing=1.25)
    _arrow(ax, 16.3, 15, 19.2, 15)
    # 2 ensembles
    _box(ax, 19.5, 4, 21, 22, tint)
    ax.text(30, 23.6, 'Coarse-grained ensembles', ha='center', va='top', fontsize=7.2, fontweight='bold')
    confs = _conformations('Q02021')
    for k, (c, cx) in enumerate(zip(confs, (23.7, 30.0, 36.3))):
        sc = 2.6 / max(1e-6, np.abs(c).max())
        xy = c * sc + np.array([cx, 15.6])
        ax.plot(xy[:, 0], xy[:, 1], '-', color='0.55', lw=0.45, zorder=2)
        ax.scatter(xy[:, 0], xy[:, 1], s=2.2, c=np.linspace(0, 1, len(xy)), cmap='viridis', zorder=3, lw=0)
    for k, (m, c) in enumerate(COL.items()):
        ax.text(30, 11.6 - 2.0 * k, m, ha='center', va='top', fontsize=6.1, color=c, fontweight='bold')
    ax.text(30, 5.0, '278–320 K', ha='center', va='center', fontsize=6.1)
    _arrow(ax, 41.0, 15, 43.9, 15)
    # 3 fluctuation relation, with the real covariance of Q02021
    _box(ax, 44.2, 4, 21, 22, tint)
    ax.text(54.7, 23.6, 'Fluctuation–response', ha='center', va='top', fontsize=7.2, fontweight='bold')
    ax.text(54.7, 19.6, r'$\chi_X = -\mathrm{Cov}\left(A,\,\partial(\beta H)/\partial X\right)$', ha='center', va='top',
            fontsize=6.6)
    rg, H = _covariance('Q02021')
    ins = ax.inset_axes([46.2, 5.0, 17.0, 9.6], transform=ax.transData)
    ins.scatter(H, rg, s=0.6, color='#2E6FA3', alpha=0.35, lw=0)
    p = np.polyfit(H, rg, 1)
    hh = np.linspace(np.percentile(H, 1), np.percentile(H, 99), 2)
    ins.plot(hh, np.polyval(p, hh), color=INK, lw=0.9)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_xlabel('energy $H$', fontsize=5.6, labelpad=1)
    ins.set_ylabel(r'$R_g$', fontsize=5.6, labelpad=1)
    for s in ('top', 'right'):
        ins.spines[s].set_visible(False)
    ins.patch.set_alpha(0)
    _arrow(ax, 65.7, 15, 68.6, 15)
    # 4 responses
    _box(ax, 68.9, 4, 14.2, 22, tint)
    ax.text(76, 23.6, 'Responses', ha='center', va='top', fontsize=7.2, fontweight='bold')
    for k, t in enumerate(('temperature', '+ hydration', 'phosphorylation', 'salt', 'substitutions')):
        y = 18.4 - 3.1 * k
        _box(ax, 70.4, y - 1.15, 11.2, 2.3, 'white', ec='0.7', r=1.0, lw=0.5)
        ax.text(76, y, t, ha='center', va='center', fontsize=5.9)
    _arrow(ax, 83.6, 15, 86.5, 15)
    # 5 tests
    _box(ax, 86.8, 4, 13.0, 22, tint)
    ax.text(93.3, 23.6, 'Tests', ha='center', va='top', fontsize=7.2, fontweight='bold')
    ax.text(93.3, 19.2, 'direct runs,\n278–320 K\n\nsmFRET, SAXS,\nDLS: 10 proteins\n\nindependent\ntemperature series',
            ha='center', va='top', fontsize=5.9, linespacing=1.2)


# ----------------------------------------------------------------------------- Fig. 1 worked examples
def example_dist(ax):
    """Rg distributions at 278 K (filled) and 320 K (line) for the two example chains, CALVADOS-2."""
    for uid, e in EX.items():
        lo, hi = _rg(uid, 278), _rg(uid, 320)
        bins = np.linspace(min(lo.min(), hi.min()), max(lo.max(), hi.max()), 40)
        ax.hist(lo, bins=bins, density=True, color=e['col'], alpha=0.35, lw=0)
        ax.hist(hi, bins=bins, density=True, histtype='step', color=e['col'], lw=1.2)
    ymax = ax.get_ylim()[1]
    for (uid, e), (x, ha) in zip(EX.items(), ((0.02, 'left'), (0.98, 'right'))):
        ax.text(x, 0.98, f"{e['name']}\n({e['org']})", transform=ax.transAxes, ha=ha, va='top', fontsize=6.3,
                color=e['col'])
    ax.set_ylim(0, ymax * 1.38)
    ax.set_xlabel(r'$R_g$ (nm)')
    ax.set_ylabel('Probability density')
    h = [matplotlib.patches.Patch(color='0.5', alpha=0.35, lw=0, label='278 K'),
         plt.Line2D([], [], color='0.4', lw=1.2, label='320 K')]
    ax.legend(handles=h, frameon=False, fontsize=6.3, loc='center right', handlelength=1.2, bbox_to_anchor=(1.0, 0.55))


def example_curve(ax):
    """<Rg>(T) relative to 300 K (points, block SE) and the 300 K fluctuation-relation tangent (dashed)."""
    for uid, e in EX.items():
        r300 = _rg(uid, 300).mean()
        m = np.array([_rg(uid, T).mean() for T in TEMPS])
        se = np.array([block_se(_rg(uid, T)) for T in TEMPS])
        ax.errorbar(TEMPS, 100 * (m / r300 - 1), yerr=100 * se / r300, fmt='o', ms=3.5, color=e['col'], elinewidth=0.8,
                    capsize=0, label=e['name'])
        chi, chi_se = _chi300(uid)
        t = np.array([276, 322])
        ax.plot(t, 100 * chi * (t - 300) / r300, '--', color=e['col'], lw=1.0)
    ax.axhline(0, color='0.75', lw=0.6, zorder=0)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'$\langle R_g\rangle$ change from 300 K (%)')
    ax.text(0.03, 0.97, 'dashed: fluctuation relation\nat 300 K', transform=ax.transAxes, va='top', fontsize=6, color='0.35')


# ----------------------------------------------------------------------------- house style (Figs 2-6, S15; 2026-10-15)
WIDTH = 6.5                                   # in, the LaTeX text width (the .tex uses width=\textwidth)
HOUSE = {'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7.5, 'axes.titlesize': 7, 'xtick.labelsize': 6.5,
         'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
         'ytick.major.width': 0.6, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5, 'xtick.major.size': 2.5,
         'ytick.major.size': 2.5, 'xtick.major.pad': 1.8, 'ytick.major.pad': 1.8, 'axes.labelpad': 2.5,
         'axes.titlepad': 3.5, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.0,
         'legend.frameon': False, 'legend.handlelength': 1.4, 'legend.handletextpad': 0.4, 'legend.borderaxespad': 0.3,
         'legend.borderpad': 0.2, 'legend.columnspacing': 1.0, 'legend.labelspacing': 0.3, 'axes.unicode_minus': True,
         'hatch.linewidth': 0.6, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300}
LETTER = dict(fontsize=9, fontweight='bold', color=INK)
BASIC, ACIDIC = '#6F8199', '#B5473A'          # K, R and D, E (figfont: sky blue and red, apart from CALVADOS-2 blue)
NU = r'$\mathtt{\nu}$'   # Greek nu from STIX italic (figfont mathtext.tt); Arial's nu reads as v
RG = r'$R_\mathrm{g}$'


def house():
    """House rcParams (Arial via figfont, 6.5-pt smallest text, 7.5-pt axis labels, no default colour cycle)."""
    import figfont
    from cycler import cycler
    plt.rcParams.update(HOUSE)
    plt.rcParams['axes.prop_cycle'] = cycler(color=[INK])
    figfont.apply()


def minus(s):
    return s.replace('-', '−')


def final_colour(c):
    """The colour a script hex becomes at save time under the active figfont palette (for colormaps, which figfont
    does not recolour)."""
    import figfont
    role = figfont.ROLE.get(c.upper()) if isinstance(c, str) else None
    return figfont.PALETTES[figfont._active].get(role, c) if role else c


def final_line(c):
    """The colour thin marks of a script hex get at save time ('<role>_line' where the palette has one)."""
    import figfont
    role = figfont.ROLE.get(c.upper()) if isinstance(c, str) else None
    pal = figfont.PALETTES[figfont._active]
    return pal.get(f'{role}_line', pal.get(role, c)) if role else c


def ff_cmap(c, name=None):
    """Single-hue colormap in one force field's colour (pale tint to dark shade), for density panels."""
    import figfont
    from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb
    b = np.array(to_rgb(final_colour(c)))
    tint = lambda a: tuple(1 - a * (1 - b))       # noqa: E731
    role = figfont.ROLE.get(to_hex(c).upper())
    line = figfont.PALETTES[figfont._active].get(f'{role}_line') if role else None
    if line:                                                        # a pale fill colour: tints of the line shade, so that
        lb = np.array(to_rgb(line))                                 # single counts stay visible on white
        return LinearSegmentedColormap.from_list(name or f'ff{c}', [tuple(1 - 0.25 * (1 - lb)), tuple(b), tuple(lb)])
    return LinearSegmentedColormap.from_list(name or f'ff{c}', [tint(0.13), tint(0.5), tuple(b), tuple(b * 0.5)])


def add_letters(fig, items, pad=2.5):
    """Bold 9-pt lowercase panel letters, top-left outside the axes. Left edge at the y label (or tick labels) of the
    axes; bottom just above the axes top, or top-aligned with an axes title. items: (ax, letter[, opts]) with opts
    dx, dy (points), row (letters of one row share the highest y), xaxes (axes whose y labels set the left edge),
    x, y (figure fraction overrides)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H, pt = fig.bbox.width, fig.bbox.height, fig.dpi / 72
    pos = []
    for it in items:
        ax, s, o = it[0], it[1], (it[2] if len(it) > 2 else {})
        x0 = ax.bbox.x0
        for a_ in o.get('xaxes', [ax]):
            x0 = min(x0, a_.bbox.x0)
            bb = a_.yaxis.get_tightbbox(r)
            if bb is not None and bb.width > 0:
                x0 = min(x0, bb.x0)
        top = ax.bbox.y1 + (pad + 7.0) * pt
        for t in (ax.title, ax._left_title, ax._right_title):
            if t.get_text() and t.get_visible():
                top = max(top, t.get_window_extent(r).y1 + 0.5 * pt)
        pos.append([x0 - 0.5 * pt + o.get('dx', 0) * pt, top + o.get('dy', 0) * pt, o])
    rows = {}
    for p in pos:
        if 'row' in p[2]:
            rows[p[2]['row']] = max(rows.get(p[2]['row'], -1e9), p[1])
    for (x, y, o), it in zip(pos, items):
        y = rows.get(o.get('row'), y) if 'row' in o else y
        fig.text(o.get('x', x / W), o.get('y', y / H), it[1], ha='left', va='top', **LETTER)


def save(fig, out):
    assert abs(fig.get_figwidth() - WIDTH) < 1e-6, fig.get_figwidth()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=300)
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 2 hydration calibration
WU = f'{TLA}/Wuttke2014'
CONSTRUCTS = [('IN_NTD', 'Integrase NTD'), ('CspM34', 'CspM34'), ('ProTa_C2_C56', 'ProTα, N half'),
              ('ProTa_C56_C110', 'ProTα, C half'), ('lambda_P6C_I84C', 'λ-repressor')]
FF = {'calvados': 'CALVADOS-2', 'recharged': 'Mpipi-Recharged'}
CAPCAL = f'{ROOT}/heldout_temperature/capped_forcing/calibration'
GAMMA_TEA = {'calvados': 2.0, 'recharged': 0.4}


def wuttke_panel(ax, n, title, first):
    import tea_calibrate as TC
    ex = pd.read_csv(f'{WU}/wuttke2014_fig2b_Rg_digitized.csv')
    e = ex[ex.construct == n].sort_values('T_K')
    r0 = np.polyval(np.polyfit(e.T_K - 300, e.Rg_nm, TC.INTERP_DEG), 0.0)
    e = e[e.T_K <= 323]                       # the calibration window 278-320 K (the series continues to 347 K)
    cur = pd.read_csv(CAPCAL + '/curves.csv')                # capped main model (capped_calibration.py, 2026-10-25)
    for model, lab in FF.items():
        c = cur[(cur.model == model) & (cur.construct == n)]
        for kind, kw in (('plain', dict(ls=(0, (3, 2)), lw=0.9, alpha=0.85)), ('capped', dict(lw=1.4))):
            t = c[c.curve == kind].sort_values('T')
            t = t[t['T'] <= 320]
            ax.plot(t['T'], t.rel, color=COL[lab], **kw)
            lowe = t[t.ess < 0.05]                            # reweighting kept < 5% effective sample: flagged
            ax.scatter(lowe['T'], lowe.rel, s=11, facecolors='white', edgecolors=COL[lab], lw=0.8, zorder=4)
    ax.errorbar(e.T_K, e.Rg_nm / r0, yerr=e.err_nm / r0, fmt='o', ms=2.8, color=INK, ecolor='0.45', elinewidth=0.6,
                capsize=0, zorder=5)
    ax.axhline(1, color='0.78', lw=0.5, zorder=0)
    ax.set_title(title)
    ax.set_xlim(275, 323)
    ax.set_xticks([280, 300, 320])
    ax.set_xlabel('Temperature (K)')
    if first:
        ax.set_ylabel(r'$R_\mathrm{g}(T)\,/\,R_\mathrm{g}$(300 K)')
    if n.startswith('lambda'):
        ax.text(0.5, 0.05, 'not fitted', transform=ax.transAxes, ha='center', fontsize=6.5, color='0.4')


def forcing_bars(ax):
    """d lambda_i / dT of each residue at the calibrated CALVADOS-2 strengths of the capped model (below the cap):
    hydrophilic (charged-group, D E K R) and TEA (hydrophobic) parts."""
    import capped_proteome as CP
    import capped_model as CM
    gh, gt = CM.params()
    aa = list('DEKRHNQSTGAVLIMFWYCP')
    hyd = np.array([gh['C'][0] * CP.F_CHG[a] for a in aa]) * 1e3
    tea = np.array([gt['C'] * CP.F_TEA[a] for a in aa]) * 1e3
    x = np.arange(len(aa))
    ax.bar(x - 0.2, hyd, 0.4, color='#2E8C6F', label='hydrophilic (charged groups)')
    ax.bar(x + 0.2, tea, 0.4, color='#C8792A', label='hydrophobic (TEA)')
    ax.axhline(0, color='0.55', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(aa)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.6, len(aa) - 0.4)
    ax.set_ylabel(r'd$\lambda_i$/d$T$ below saturation' + '\n' + r'($10^{-3}$ K$^{-1}$, CALVADOS-2)')
    ax.legend(loc='upper right', handlelength=1)
    ax.text(0.99, 0.66, 'positive values make contacts\nstickier on heating', transform=ax.transAxes, fontsize=6.5,
            color='0.4', ha='right', va='top')


def loo_scatter(ax):
    """Slope at 300 K of each construct predicted with it left out of the fit, against the measurement."""
    lim = (-0.42, 0.16)
    ax.plot(lim, lim, color='0.6', lw=0.6, ls='--')
    ax.axhline(0, color='0.85', lw=0.5)
    ax.axvline(0, color='0.85', lw=0.5)
    loo = pd.read_csv(CAPCAL + '/loo.csv')                   # capped main model (capped_calibration.py)
    for model, lab in FF.items():
        t = loo[loo.model == model]
        ax.scatter(t.slope_exp, t.slope_plain, s=12, facecolors='white', edgecolors=COL[lab], lw=0.8, zorder=3)
        ax.scatter(t.slope_exp, t.slope_fit, s=12, color=COL[lab], lw=0, zorder=4, label=lab)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xticks([-0.4, -0.2, 0.0])
    ax.set_xlabel('Measured slope (% per K)')
    ax.set_ylabel('Predicted slope, construct\nleft out (% per K)')
    h = [plt.Line2D([], [], marker='o', ls='', color='0.3', ms=3.5, label='calibrated'),
         plt.Line2D([], [], marker='o', ls='', mfc='white', color='0.3', ms=3.5, label='plain')]
    ax.legend(handles=h, loc='lower right', handletextpad=0.2)


def residue_sets(ax):
    """Leave-one-out slope error for each residue set of the hydrophilic forcing, with the plain-model error."""
    order = ['all (DEKRNQSTH+G)', 'without G (backbone)', 'without R', 'charged only (DEKR)', 'polar + G only',
             'acidic only (DE)']
    short = ['all', '−Gly', '−Arg', 'DEKR', 'polar+Gly', 'DE']
    x = np.arange(len(order))
    plain = []
    for o, (model, lab) in zip((-0.2, 0.2), FF.items()):
        t = pd.read_csv(f'{WU}/calibration/{model}/sensitivity.csv')
        t = t[(t.interp_deg == 3) & (t.tmax.astype(str) == '320.0')].set_index('hyd residues')
        ax.bar(x + o, t.loc[order, 'loo_slope_err'], 0.38, color=COL[lab], label=lab)
        ax.axhline(t.loo_slope_err_plain.iloc[0], color=COL[lab], lw=0.8, ls=':')
        plain.append(t.loo_slope_err_plain.iloc[0])
    ax.text(len(order) - 0.45, max(plain) + 0.006, 'plain', ha='right', va='bottom', fontsize=6.5, color='0.4')
    ax.set_xticks(x)
    ax.set_xticklabels(short, rotation=45, ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylabel('Left-out slope error\n(% per K)')
    ax.set_ylim(0, 0.38)


def si_residue_sets(out):
    """Fig. S (2026-10-30): the residue-set comparison of the linear form, panel h of the former Fig. 4, moved to the
    supplement when the forcing and calibration figures were merged (fig_forcing.plot_merged)."""
    fig = plt.figure(figsize=(WIDTH, 1.8))      # full width, printed at 100 % like every figure
    ax = fig.add_axes([0.25, 0.25, 0.5, 0.7])
    residue_sets(ax)
    h = [matplotlib.patches.Patch(color=COL[lab], label=lab) for lab in FF.values()]
    ax.legend(handles=h, loc='upper left', bbox_to_anchor=(1.02, 1.0), handlelength=1.0, borderaxespad=0)
    save(fig, out)


def fig2(out):
    fig = plt.figure(figsize=(WIDTH, 4.6))
    top = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1, 0.17, 1], left=0.085, right=0.99, top=0.93, bottom=0.645,
                           wspace=0.11)
    axes = [fig.add_subplot(top[0, i]) for i in (0, 1, 2, 3, 5)]
    for i, (ax, (n, t)) in enumerate(zip(axes, CONSTRUCTS)):
        wuttke_panel(ax, n, t, i == 0)
        ax.set_ylim(0.84, 1.12)
        if 0 < i < 4:
            ax.set_yticklabels([])
    axes[4].set_ylim(0.84, 1.31)              # lambda-repressor: its measured 282 K point (1.27) is above the others
    axes[4].set_yticks([0.9, 1.0, 1.1, 1.2, 1.3])
    axes[4].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
    h = [plt.Line2D([], [], marker='o', ls='', color=INK, ms=2.8, label='measured (smFRET)'),
         plt.Line2D([], [], color=COL['CALVADOS-2'], lw=1.4, label='CALVADOS-2'),
         plt.Line2D([], [], color=COL['Mpipi-Recharged'], lw=1.4, label='Mpipi-Recharged'),
         plt.Line2D([], [], color='0.4', lw=1.4, label='calibrated'),
         plt.Line2D([], [], color='0.4', lw=0.9, ls=(0, (3, 2)), label='plain')]
    fig.legend(handles=h, loc='center', ncol=5, bbox_to_anchor=(0.53, 0.525), handlelength=2.0, columnspacing=1.6)
    bot = fig.add_gridspec(1, 3, width_ratios=[1.65, 1, 1.05], left=0.095, right=0.99, top=0.455, bottom=0.12,
                           wspace=0.62)
    low = []
    for k, fn in enumerate((forcing_bars, loo_scatter, residue_sets)):
        ax = fig.add_subplot(bot[0, k])
        fn(ax)
        low.append(ax)
    add_letters(fig, [(a, 'abcde'[i], {'row': 1}) for i, a in enumerate(axes)] +
                [(a, 'fgh'[i], {'row': 2}) for i, a in enumerate(low)])
    save(fig, out)


# ----------------------------------------------------------------------------- Fig. 3 proteome totals
TR = f'{FDT}/total_response'
NU_BIN_LAB = ['0.30–0.45', '0.45–0.50', '0.50–0.54', '0.54–0.58', '0.58–0.63', '>0.63']
KINGDOMS = ['Bacteria', 'Plants', 'Fungi', 'Mammals', 'Protists', 'Viruses']
AA = list('DEKRHNQSTGPAVLIMCFWY')


def _load_totals():
    import total_response_evidence as T
    D, _ = T.load()
    bio = T.bio_table()
    return D, bio


def total_by_nu(ax):
    a = pd.read_csv(f'{TR}/A_by_nu.csv')
    x = np.arange(6)
    w = 0.26
    for k, m in enumerate(['CALVADOS-2', 'Mpipi-GG', 'Mpipi-Recharged']):
        d = a[(a.model == m) & (a.own_nu_bin != 'all')]
        xx = x + (k - 1) * w
        ax.bar(xx, d.total_median, width=w * 0.92, color=COL[m], label=m)
        ax.errorbar(xx, d.total_median, yerr=[d.total_median - d.total_q25, d.total_q75 - d.total_median], fmt='none',
                    ecolor=INK, lw=0.6, capsize=0)
        ax.scatter(xx, d.chi_median, marker='D', s=8, color='white', edgecolor=COL[m], lw=0.8, zorder=5)
    ax.axhline(0, color=INK, lw=0.6)
    ax.axvline(1.5, color='0.5', lw=0.6, ls=':')
    ax.text(1.56, 0.42, 'θ point', fontsize=6.5, color='0.4', va='top')
    ax.set_xticks(x)
    ax.set_xticklabels(NU_BIN_LAB)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.5, 5.5)
    ax.set_xlabel(f'Compaction, {NU} at 300 K (each model its own)')
    ax.set_ylabel('Total response (% per K)')
    low = a[a.own_nu_bin != 'all'].total_q25.min()           # every bar and interquartile line inside the axes
    ax.set_ylim(np.floor(10 * low) / 10 - 0.35, 0.45)
    h = [matplotlib.patches.Patch(color=COL[m], label=m) for m in COL] + \
        [plt.Line2D([], [], marker='D', ls='', mfc='white', mec='0.4', ms=3.2, label='fixed-solvent susceptibility')]
    ax.legend(handles=h, loc='lower right', ncol=2, handlelength=1, columnspacing=1.0)


def agreement(ax):
    b = pd.read_csv(f'{TR}/B_agreement.csv')
    b = b[b.group == 'all']
    rows = [('CALVADOS-2 vs Mpipi-GG', 'three-model set', 'C2–GG'), ('CALVADOS-2 vs Mpipi-Recharged', 'pair', 'C2–MR'),
            ('Mpipi-GG vs Mpipi-Recharged', 'three-model set', 'GG–MR')]
    for k, (p, sset, lab) in enumerate(rows):
        y = 2 - k
        r = b[(b.pair == p) & (b.set == sset)].iloc[0]
        ax.plot([r.rho_chi, r.rho_total], [y, y], color='0.85', lw=2, zorder=1, solid_capstyle='butt')
        ax.scatter(r.rho_chi, y, s=18, color='white', edgecolor='0.45', lw=0.8, zorder=3)
        ax.scatter(r.rho_total, y, s=20, color=INK, lw=0, zorder=4)
        right = r.rho_total >= r.rho_chi
        ax.text(r.rho_total + (0.016 if right else -0.016), y, f'{r.rho_total:.2f}', va='center',
                ha='left' if right else 'right', fontsize=6.5)
    ax.set_yticks(range(3))
    ax.set_yticklabels([r[2] for r in rows][::-1])
    ax.set_xlim(0.6, 0.95)
    ax.set_xticks([0.6, 0.7, 0.8, 0.9])
    ax.set_ylim(-0.5, 2.5)
    ax.set_xlabel('Rank agreement (Spearman)')
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    h = [plt.Line2D([], [], marker='o', ls='', mfc='white', mec='0.45', ms=3.6, label='fixed-solvent\nsusceptibility'),
         plt.Line2D([], [], marker='o', ls='', color=INK, ms=3.6, label='total\nresponse')]
    ax.legend(handles=h, loc='lower left', bbox_to_anchor=(-0.04, 1.0), ncol=2, handletextpad=0.2, columnspacing=0.8,
              borderaxespad=0)


def nu_lambda_map(ax, D, bio, fig, cax=None):
    """Median total response (CALVADOS-2) on a grid of nu and mean stickiness <lambda>."""
    d = D['C'][['nu', 'total']].join(bio[['lambda_avg']], how='inner')
    d = d[(d.nu > 0.3) & (d.nu < 0.7)]
    xb = np.linspace(0.38, 0.68, 13)
    yb = np.linspace(d.lambda_avg.quantile(0.01), d.lambda_avg.quantile(0.99), 12)
    M = np.full((len(yb) - 1, len(xb) - 1), np.nan)
    for i in range(len(yb) - 1):
        for j in range(len(xb) - 1):
            s = d[(d.nu >= xb[j]) & (d.nu < xb[j + 1]) & (d.lambda_avg >= yb[i]) & (d.lambda_avg < yb[i + 1])]
            if len(s) >= 8:
                M[i, j] = s.total.median()
    im = ax.pcolormesh(xb, yb, M, cmap=ff_cmap(COL['CALVADOS-2']).reversed(), vmin=-1.0, vmax=0.0, shading='flat')
    ax.set_xlabel(f'Compaction, {NU}')
    ax.set_ylabel('Mean stickiness (mean λ)')
    cb = fig.colorbar(im, cax=cax) if cax is not None else fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03)
    cb.set_label('Median total response\n(% per K, CALVADOS-2)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_visible(False)


KCOL = {'Bacteria': '#3B3B3B', 'Plants': '#7A9E6B', 'Fungi': '#E0B44C', 'Mammals': '#C46A5E', 'Protists': '#8E7CB0',
        'Viruses': '#C9C9C9'}                  # taxonomic groups in the composition grid (not figfont roles)


def grid_kingdoms():
    """Share of each taxonomic group among the proteins of each square of the composition grid (E_grid_sharedCR), from
    the same membership as total_response_evidence.grid (the shared CALVADOS-2 + Mpipi-Recharged leftover, all paired
    chains); cached as E_grid_sharedCR_kingdoms.csv, and the counts are checked against E_grid_sharedCR.csv."""
    f = f'{TR}/E_grid_sharedCR_kingdoms.csv'
    if not os.path.exists(f):
        import total_response_evidence as T
        D, _ = T.load()
        bio = T.bio_table()
        y = T.leftovers(D, bio)[-1].L
        b = bio.reindex(y.index)
        tp = pd.cut(b.f_pi, [-1, 0.068, 0.156, 1], labels=['low', 'mid', 'high'])
        tc = pd.cut(b.f_DEK, [-1, 0.18, 0.34, 1], labels=['low', 'mid', 'high'])
        k = pd.crosstab([tp.rename('pi_third'), tc.rename('DEK_third')], b.kingdom).reset_index()
        k.to_csv(f, index=False)
    k = pd.read_csv(f).set_index(['pi_third', 'DEK_third'])
    ref = pd.read_csv(f'{TR}/E_grid_sharedCR.csv').set_index(['pi_third', 'DEK_third'])
    assert (k.sum(axis=1).reindex(ref.index) == ref.n).all(), 'grid counts differ from E_grid_sharedCR.csv'
    return k


def grid_members():
    """The square of the composition grid of every protein (same membership as total_response_evidence.grid, the shared
    CALVADOS-2 + Mpipi-Recharged leftover over all paired chains), cached as E_grid_sharedCR_members.csv and checked
    against the counts and medians of E_grid_sharedCR.csv."""
    f = f'{TR}/E_grid_sharedCR_members.csv'
    if not os.path.exists(f):
        import total_response_evidence as T
        D, _ = T.load()
        bio = T.bio_table()
        y = T.leftovers(D, bio)[-1].L
        b = bio.reindex(y.index)
        m = pd.DataFrame({'pi_third': pd.cut(b.f_pi, [-1, 0.068, 0.156, 1], labels=['low', 'mid', 'high']),
                          'DEK_third': pd.cut(b.f_DEK, [-1, 0.18, 0.34, 1], labels=['low', 'mid', 'high']),
                          'leftover': y})
        m.index.name = 'UniProt_ID'
        m.to_csv(f)
    m = pd.read_csv(f, index_col=0)
    ref = pd.read_csv(f'{TR}/E_grid_sharedCR.csv').set_index(['pi_third', 'DEK_third'])
    g = m.groupby(['pi_third', 'DEK_third']).leftover.agg(['median', 'size']).reindex(ref.index)
    assert (g['size'] == ref.n).all() and np.allclose(g['median'], ref['median']), 'grid membership differs'
    return m


GRID_ROWS = [(p_, c_) for p_ in ('high', 'mid', 'low') for c_ in ('low', 'mid', 'high')]   # top to bottom, as Fig. 4d


def si_grid_composition(out):
    """Fig. S: full composition of the nine squares of Fig. 4d. a extra compaction of each square and the enrichment of
    each residue type (mean fraction in the square over the mean of all squares, log2); b taxonomic groups and number of
    proteins; c-g length, nu (CALVADOS-2), NCPR, kappa and the total response (CALVADOS-2, Mpipi-Recharged)."""
    import total_response_evidence as T
    m = grid_members()
    bio = T.bio_table().reindex(m.index)
    D, _ = T.load()
    lv = {'low': 'low', 'mid': 'mid', 'high': 'high'}
    rows = GRID_ROWS
    lab = [f'π {p_}, D+E+K {c_}' for p_, c_ in rows]
    sq = [m.index[(m.pi_third == p_) & (m.DEK_third == c_)] for p_, c_ in rows]
    n = len(rows)
    W, H = WIDTH, 5.15
    fig = plt.figure(figsize=(W, H))
    ax_in = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])  # noqa: E731
    ytop, yh = 2.95, 1.95
    # a extra compaction strip and residue enrichment
    ref = pd.read_csv(f'{TR}/E_grid_sharedCR.csv').set_index(['pi_third', 'DEK_third'])
    ex = np.array([-ref.loc[r, 'median'] for r in rows])
    a0 = ax_in(1.12, ytop, 0.34, yh)
    a0.imshow(ex[:, None], cmap='PRGn', vmin=-0.75, vmax=0.75, aspect='auto')
    for i, v in enumerate(ex):
        a0.text(0, i, minus(f'{v:+.2f}'), ha='center', va='center', fontsize=6.5, color='white' if abs(v) > 0.45 else INK)
    a0.set_xticks([0], ['extra\ncompaction\n(SD)'], fontsize=6.5)
    a0.set_yticks(range(n), lab, fontsize=6.5)
    a0.tick_params(length=0)
    for sp_ in a0.spines.values():
        sp_.set_visible(False)
    a = ax_in(1.52, ytop, 3.3, yh)
    allm = bio[[f'f_{x}' for x in AA]].mean()
    E = np.array([[np.log2(bio.loc[s_, f'f_{x}'].mean() / allm[f'f_{x}']) for x in AA] for s_ in sq])
    im = a.imshow(np.clip(E, -2, 2), cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    a.set_xticks(range(len(AA)), AA)
    a.xaxis.tick_top()
    a.set_yticks([])
    a.tick_params(length=0)
    for sp_ in a.spines.values():
        sp_.set_visible(False)
    for yb in (2.5, 5.5):
        for ax_ in (a0, a):
            ax_.axhline(yb, color='white', lw=1.5)
    cax = ax_in(4.9, ytop + 0.2, 0.06, yh - 0.4)
    cb = fig.colorbar(im, cax=cax, extend='both')
    cb.set_label('Residue enrichment\n(log$_2$, vs all squares)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_visible(False)
    # b taxonomic groups
    b = ax_in(5.55, ytop, 0.62, yh)
    for i, s_ in enumerate(sq):
        sh = bio.loc[s_, 'kingdom'].value_counts(normalize=True).reindex(KINGDOMS).fillna(0)
        left = 0.0
        for kg in KINGDOMS:
            b.barh(i, sh[kg], left=left, height=0.72, color=KCOL[kg], lw=0)
            left += sh[kg]
        b.text(1.04, i, f'{len(s_):,}', va='center', fontsize=6.5)
    b.set_ylim(n - 0.5, -0.5)
    b.set_xlim(0, 1)
    b.set_xticks([])
    b.set_yticks([])
    for sp_ in b.spines.values():
        sp_.set_visible(False)
    b.text(1.04, -0.75, 'n', fontsize=6.5, va='center', fontstyle='italic')
    b.text(0.5, -0.75, 'groups', fontsize=6.5, va='center', ha='center')
    h = [matplotlib.patches.Patch(color=KCOL[kg], label=kg) for kg in KINGDOMS]
    b.legend(handles=h, loc='upper center', bbox_to_anchor=(0.55, -0.02), ncol=2, handlelength=0.8, columnspacing=0.6,
             handletextpad=0.3, borderaxespad=0, fontsize=6.5)
    # c-g distributions
    dist = [('Length\n(residues)', [bio.loc[s_, 'length'] for s_ in sq], None, True),
            (f'Compaction, {NU}\n(CALVADOS-2)', [D['C'].nu.reindex(s_).dropna() for s_ in sq], None, False),
            ('Net charge\nper residue', [bio.loc[s_, 'NCPR'] for s_ in sq], None, False),
            ('Charge\npatterning, κ', [bio.loc[s_, 'kappa'].dropna() for s_ in sq], None, False),
            ('Total response\n(% per K)', [D['C'].total.reindex(s_).dropna() for s_ in sq],
             [D['R'].total.reindex(s_).dropna() for s_ in sq], False)]
    xs = [1.12, 2.18, 3.24, 4.3, 5.36]
    low = []
    for k, (t, d1, d2, logx) in enumerate(dist):
        ax = ax_in(xs[k], 0.42, 0.98, 1.85)
        pos = np.arange(n)
        sets = [(d1, COL['CALVADOS-2'] if d2 is not None else '0.55', 0.0 if d2 is None else -0.18)]
        if d2 is not None:
            sets.append((d2, COL['Mpipi-Recharged'], 0.18))
        for dd, col, off in sets:
            bp = ax.boxplot([np.asarray(v, float) for v in dd], positions=pos + off, vert=False, showfliers=False,
                            widths=0.3 if d2 is not None else 0.6, patch_artist=True,
                            medianprops=dict(color='white', lw=0.8), whiskerprops=dict(color=col, lw=0.6),
                            capprops=dict(color=col, lw=0.6), boxprops=dict(lw=0))
            for p_ in bp['boxes']:
                p_.set_facecolor(col)
        ax.set_ylim(n - 0.5, -0.5)
        if logx:
            ax.set_xscale('log')
            ax.set_xticks([30, 100, 300], ['30', '100', '300'])
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        if d2 is not None:
            ax.axvline(0, color=INK, lw=0.5)
        for yb in (2.5, 5.5):
            ax.axhline(yb, color='0.85', lw=0.6)
        ax.set_yticks(range(n), lab if k == 0 else [], fontsize=6.5)
        ax.tick_params(axis='y', length=0)
        ax.set_xlabel(t)
        low.append(ax)
    low[-1].legend(handles=[matplotlib.patches.Patch(color=COL[m_], label=m_) for m_ in ('CALVADOS-2', 'Mpipi-Recharged')],
                   loc='lower right', bbox_to_anchor=(1.05, 1.0), handlelength=0.8, borderaxespad=0.1, fontsize=6.5)
    low[0].set_yticks(range(n), lab, fontsize=6.5)
    add_letters(fig, [(a0, 'a', {'row': 0}), (b, 'b', {'row': 0})] +
                [(ax, 'cdefg'[k], {'row': 1}) for k, ax in enumerate(low)])
    save(fig, out)
    return {r: len(s_) for r, s_ in zip(rows, sq)}


def composition_thirds(ax, fig, cax=None, kingdoms=True):
    """Extra compaction on heating by thirds of the pi and D+E+K fractions; with kingdoms, each square also gives its
    number of proteins and a bar of the share of each taxonomic group."""
    g = pd.read_csv(f'{TR}/E_grid_sharedCR.csv')
    lv = ['low', 'mid', 'high']
    M = g.pivot(index='pi_third', columns='DEK_third', values='median').reindex(index=lv[::-1], columns=lv)
    N = g.pivot(index='pi_third', columns='DEK_third', values='n').reindex(index=lv[::-1], columns=lv)
    v = 0.75
    im = ax.imshow(-M.to_numpy(), cmap='PRGn', vmin=-v, vmax=v, aspect='auto')
    K = grid_kingdoms() if kingdoms else None
    for i in range(3):
        for j in range(3):
            val = -M.iloc[i, j]
            col = 'white' if abs(val) > 0.45 else INK
            if not kingdoms:
                ax.text(j, i, minus(f'{val:+.2f}'), ha='center', va='center', fontsize=6.5, color=col)
                continue
            ax.text(j, i - 0.2, minus(f'{val:+.2f}'), ha='center', va='center', fontsize=7, fontweight='bold', color=col)
            ax.text(j, i + 0.06, f'n = {int(N.iloc[i, j]):,}', ha='center', va='center', fontsize=6.5, color=col)
            sh = K.loc[(lv[::-1][i], lv[j])]
            sh = sh / sh.sum()
            x0 = j - 0.42
            ax.add_patch(matplotlib.patches.Rectangle((x0 - 0.02, i + 0.2), 0.88, 0.2, fc='white', ec='none', zorder=3))
            for kg in KINGDOMS:
                ax.add_patch(matplotlib.patches.Rectangle((x0, i + 0.22), 0.84 * sh[kg], 0.16, fc=KCOL[kg], ec='none',
                                                          zorder=4))
                x0 += 0.84 * sh[kg]
    ax.set_xticks(range(3))
    ax.set_xticklabels(lv)
    ax.set_yticks(range(3))
    ax.set_yticklabels(lv[::-1])
    ax.tick_params(length=0)
    ax.set_xlabel('D + E + K fraction')
    ax.set_ylabel('π fraction')
    for s_ in ax.spines.values():
        s_.set_visible(False)
    if kingdoms:
        h = [matplotlib.patches.Patch(color=KCOL[kg], label=kg) for kg in KINGDOMS]
        ax.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=3, handlelength=0.8,
                  columnspacing=0.8, handletextpad=0.3, borderaxespad=0, fontsize=6.5)
    cb = fig.colorbar(im, cax=cax) if cax is not None else fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label('Extra compaction on heating (SD)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_visible(False)


def by_kingdom(ax, D):
    x = np.arange(len(KINGDOMS))
    for o, (k, m) in zip((-0.2, 0.2), (('C', 'CALVADOS-2'), ('R', 'Mpipi-Recharged'))):
        d = D[k]
        if 'kingdom' not in d.columns:
            d = d.join(D['C'][['kingdom']], how='left')
        data = [d.total[d.kingdom == g].dropna().values for g in KINGDOMS]
        bp = ax.boxplot(data, positions=x + o, widths=0.34, showfliers=False, patch_artist=True,
                        medianprops=dict(color='white', lw=0.9), whiskerprops=dict(color=COL[m], lw=0.6),
                        capprops=dict(color=COL[m], lw=0.6), boxprops=dict(lw=0))
        for p in bp['boxes']:
            p.set_facecolor(COL[m])
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(KINGDOMS, rotation=40, ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.6, len(KINGDOMS) - 0.4)
    ax.set_ylabel('Total response (% per K)')
    h = [matplotlib.patches.Patch(color=COL[m], label=m) for m in ('CALVADOS-2', 'Mpipi-Recharged')]
    ax.legend(handles=h, loc='lower right', bbox_to_anchor=(1.03, 1.0), ncol=1, handlelength=1, columnspacing=0.8,
              borderaxespad=0)


def residue_enrichment(ax, D, bio):
    """log2 enrichment of each residue type in the 5% most and 5% least heat-sensitive chains, relative to all."""
    f = bio[[f'f_{a}' for a in AA]]
    x = np.arange(len(AA))
    vals = []
    for o, (k, m) in zip((-0.18, 0.18), (('C', 'CALVADOS-2'), ('R', 'Mpipi-Recharged'))):
        t = D[k][['total']].join(f, how='inner').dropna()
        lo, hi = t.total.quantile(0.05), t.total.quantile(0.95)
        base = t[[f'f_{a}' for a in AA]].mean()
        most = np.log2(t[t.total <= lo][[f'f_{a}' for a in AA]].mean() / base)
        least = np.log2(t[t.total >= hi][[f'f_{a}' for a in AA]].mean() / base)
        ax.scatter(x + o, most, s=10, color=COL[m], lw=0, zorder=3)
        ax.scatter(x + o, least, s=10, facecolors='white', edgecolors=COL[m], lw=0.7, zorder=3)
        vals += list(most) + list(least)
    ax.axhline(0, color='0.6', lw=0.5)
    for xv in x[::2]:
        ax.axvspan(xv - 0.5, xv + 0.5, color='0.95', zorder=0, lw=0)
    ax.set_xticks(x)
    ax.set_xticklabels(AA)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.5, len(AA) - 0.5)
    v = np.array(vals)
    assert np.isfinite(v).all()
    ax.set_ylim(min(-1.6, np.floor(2 * v.min()) / 2 - 0.2), max(1.6, np.ceil(2 * v.max()) / 2 + 0.2))   # every point inside
    ax.set_ylabel(r'log$_2$ enrichment')
    h = [plt.Line2D([], [], marker='o', ls='', color='0.35', ms=3.4, label='5% most heat-sensitive'),
         plt.Line2D([], [], marker='o', ls='', mfc='white', color='0.35', ms=3.4, label='5% least heat-sensitive'),
         matplotlib.patches.Patch(color=COL['CALVADOS-2'], label='CALVADOS-2'),
         matplotlib.patches.Patch(color=COL['Mpipi-Recharged'], label='Mpipi-Recharged')]
    ax.legend(handles=h, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=4, handletextpad=0.3, handlelength=1.0,
              columnspacing=1.3, borderaxespad=0)


def cross_model_scatter(ax, D):
    """Every chain simulated in both force fields (limits include the most contracting ones)."""
    t = D['C'][['total']].join(D['R'][['total']], how='inner', lsuffix='_C', rsuffix='_R')
    lo = np.floor(10 * t.min().min()) / 10 - 0.1               # from the data: every chain inside the axes
    hi = np.ceil(10 * t.max().max()) / 10 + 0.05
    assert t.total_C.min() > lo and t.total_R.min() > lo and t.max().max() < hi, 'chains outside the axes'
    ax.hexbin(t.total_C, t.total_R, gridsize=46, bins='log', cmap=ff_cmap(INK, 'ink'), mincnt=1, linewidths=0,
              extent=(lo, hi, lo, hi))
    ax.plot([lo, hi], [lo, hi], color='0.55', lw=0.6, ls='--', zorder=0)
    from scipy import stats
    rho = stats.spearmanr(t.total_C, t.total_R).statistic
    ax.text(0.04, 0.97, f'ρ = {rho:.2f}\n{len(t):,} chains', transform=ax.transAxes, va='top', fontsize=6.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xticks([-3, -2, -1, 0])
    ax.set_yticks([-3, -2, -1, 0])
    ax.set_aspect('equal')
    ax.set_xlabel('CALVADOS-2 (% per K)')
    ax.set_ylabel('Mpipi-Recharged\n(% per K)')
    return rho, len(t)


def fig3(out):
    D, bio = _load_totals()
    fig = plt.figure(figsize=(WIDTH, 6.15))
    r1 = fig.add_gridspec(1, 2, width_ratios=[2.45, 1], left=0.085, right=0.985, top=0.955, bottom=0.725, wspace=0.42)
    r2 = fig.add_gridspec(1, 8, width_ratios=[0.78, 0.04, 0.8, 1.32, 0.04, 0.8, 0.94, 0.01], left=0.085, right=0.985,
                          top=0.6, bottom=0.39, wspace=0.08)
    r3 = fig.add_gridspec(1, 2, width_ratios=[2.45, 1], left=0.085, right=0.985, top=0.26, bottom=0.065, wspace=0.42)
    cax_c, cax_d = fig.add_subplot(r2[0, 1]), fig.add_subplot(r2[0, 4])
    ax = {k: fig.add_subplot(s) for k, s in (('a', r1[0, 0]), ('b', r1[0, 1]), ('c', r2[0, 0]), ('d', r2[0, 3]),
                                              ('e', r2[0, 6]), ('f', r3[0, 0]), ('g', r3[0, 1]))}
    total_by_nu(ax['a'])
    agreement(ax['b'])
    nu_lambda_map(ax['c'], D, bio, fig, cax_c)
    composition_thirds(ax['d'], fig, cax_d)
    by_kingdom(ax['e'], D)
    residue_enrichment(ax['f'], D, bio)
    rho, n = cross_model_scatter(ax['g'], D)
    add_letters(fig, [(ax['a'], 'a', {'row': 1}), (ax['b'], 'b', {'row': 1}), (ax['c'], 'c', {'row': 2}),
                      (ax['d'], 'd', {'row': 2}), (ax['e'], 'e', {'row': 2}), (ax['f'], 'f', {'row': 3}),
                      (ax['g'], 'g', {'row': 3})])
    save(fig, out)
    return {'rho_g': round(rho, 2), 'n_g': n}


# ----------------------------------------------------------------------------- Fig. 5 sequence tracks
HCX = f'{TR}/heat_crosses'
TRACK_EX = [('P17096', 'HMGA1, human'), ('P81558', 'Myelin basic protein, pig')]
PHOS = '#C8792A'                              # annotated phosphosite (figfont: orange)


def _annotated(uid):
    import phosphosite_test as pt
    S = pd.read_csv(os.path.join(pt.OUT, 'sites.csv'))
    S = S[(S.model == 'CALVADOS-2') & (S.UniProt_ID == uid)]
    return set(S.pos[S.annotated])


def charge_track(a0, pos, seq):
    q = np.array([1 if a in 'KR' else -1 if a in 'DE' else 0 for a in seq])
    a0.bar(pos, q, width=1.0, color=np.where(q > 0, BASIC, np.where(q < 0, ACIDIC, 'none')), lw=0)
    a0.set_ylim(-1.2, 1.2)
    a0.set_yticks([])
    a0.set_ylabel('Charge', rotation=0, ha='right', va='center', labelpad=4)
    plt.setp(a0.get_xticklabels(), visible=False)
    a0.tick_params(axis='x', length=0)
    for sp in ('left', 'bottom'):
        a0.spines[sp].set_visible(False)


def sequence_track(fig, spec, uid, title):
    """Three tracks along the sequence (CALVADOS-2, 300 K): charge; each residue's share of the total response; change
    of the total response on phosphorylating each Ser/Thr (pale bar from zero, marker, 10-block SE with caps; filled:
    annotated phosphosite). Returns the three axes."""
    r = pd.read_csv(f'{HCX}/residues.csv')
    r = r[r.UniProt_ID == uid].sort_values('pos')
    s = pd.read_csv(f'{HCX}/sites.csv')
    s = s[s.UniProt_ID == uid].sort_values('pos')
    ann = _annotated(uid)
    seq = ''.join(r.residue)
    pos = r.pos.to_numpy()
    g = spec.subgridspec(3, 1, height_ratios=[0.3, 1, 1], hspace=0.14)
    a0 = fig.add_subplot(g[0])
    charge_track(a0, pos, seq)
    a0.set_title(title, loc='left')
    a1 = fig.add_subplot(g[1], sharex=a0)
    v = 100 * r.total_share.to_numpy()
    a1.bar(pos, v, width=0.85, color=np.where(v < 0, '0.62', '0.84'), lw=0)
    a1.axhline(0, color='0.5', lw=0.5)
    a1.set_ylabel('Share of total\nresponse (%)', rotation=0, ha='right', va='center', labelpad=4)
    plt.setp(a1.get_xticklabels(), visible=False)
    a2 = fig.add_subplot(g[2], sharex=a0)
    filled = s.pos.isin(ann).to_numpy()
    for f_ in (False, True):
        x = s[filled == f_]
        c = PHOS if f_ else '0.45'
        a2.bar(x.pos, x.phospho_dtotal, width=0.9, color=c, alpha=0.3, lw=0, zorder=1)
        a2.errorbar(x.pos, x.phospho_dtotal, yerr=x.phospho_dtotal_se, fmt='none', ecolor=c, elinewidth=0.6,
                    capsize=1.3, capthick=0.6, zorder=2)
        a2.plot(x.pos, x.phospho_dtotal, 'o', ms=2.6, color=c, mfc=c if f_ else 'white', mew=0.7, ls='', zorder=3)
    a2.axhline(0, color='0.5', lw=0.5, zorder=0)
    a2.set_ylabel('Δ total response\non phosphorylation\n(% per K)', rotation=0, ha='right', va='center', labelpad=4)
    a2.set_xlabel('Residue')
    a2.set_xlim(0.3, len(seq) + 0.7)
    for a in (a1, a2):
        a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    return a0, a1, a2


def track_legend(fig, loc, anchor, ncol=5):
    h = [plt.Line2D([], [], marker='o', color=PHOS, ms=2.6, ls='', label='annotated phosphosite'),
         plt.Line2D([], [], marker='o', color='0.45', mfc='white', ms=2.6, ls='', label='other Ser/Thr'),
         matplotlib.patches.Patch(color=BASIC, label='K, R'), matplotlib.patches.Patch(color=ACIDIC, label='D, E')]
    return fig.legend(handles=h, loc=loc, ncol=ncol, bbox_to_anchor=anchor, handlelength=1.0, columnspacing=1.4)


def fig5_tracks(out):
    fig = plt.figure(figsize=(WIDTH, 4.1))
    gs = fig.add_gridspec(2, 1, hspace=0.42, top=0.94, bottom=0.11, left=0.2, right=0.99)
    tr = [sequence_track(fig, gs[k], uid, t) for k, (uid, t) in enumerate(TRACK_EX)]
    track_legend(fig, 'lower center', (0.56, -0.005))
    add_letters(fig, [(t[0], 'cd'[k], {'xaxes': t}) for k, t in enumerate(tr)])
    save(fig, out)


# ----------------------------------------------------------------------------- Supplementary Fig. S15, charged-group form
S1D = f'{TR}/s1_sensitivity'
CG = 'charged-group'


def si7_s1(out):
    """Proteome results with the frozen hydrophilic forcing against the charged-group form (s1_proteome.py)."""
    from scipy import stats
    fig = plt.figure(figsize=(WIDTH, 4.8))
    r1 = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05], left=0.085, right=0.985, top=0.935, bottom=0.6, wspace=0.5)
    r2 = fig.add_gridspec(1, 2, width_ratios=[1, 1.75], left=0.085, right=0.985, top=0.4, bottom=0.11, wspace=0.42)
    lo, hi = -3.7, 0.4
    out_numbers = {}
    axs = []
    for j, (k, m) in enumerate((('C', 'CALVADOS-2'), ('R', 'Mpipi-Recharged'))):
        ax = fig.add_subplot(r1[0, j])
        d = pd.read_csv(f'{S1D}/chains_{k}.csv')
        assert d[['total', 'total_S1']].min().min() > lo and d[['total', 'total_S1']].max().max() < hi
        ax.hexbin(d.total, d.total_S1, gridsize=44, bins='log', cmap=ff_cmap(COL[m]), mincnt=1, linewidths=0,
                  extent=(lo, hi, lo, hi))
        ax.plot([lo, hi], [lo, hi], color='0.55', lw=0.6, ls='--', zorder=0)
        ax.axhline(0, color='0.75', lw=0.5, zorder=0)
        rho = stats.spearmanr(d.total, d.total_S1).statistic
        pc = 100 * (d.total_S1 < 0).mean()
        out_numbers[m] = (round(rho, 2), round(pc, 1), len(d))
        ax.set_title(m, loc='left', color=COL[m])
        ax.text(0.97, 0.04, f'ρ = {rho:.2f}\n{pc:.1f}% contract\n({CG} form)', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=6.5)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xticks([-3, -2, -1, 0])
        ax.set_yticks([-3, -2, -1, 0])
        ax.set_aspect('equal')
        ax.set_xlabel('Calibrated model (% per K)')
        ax.set_ylabel(f'Charged-group form (% per K)' if j == 0 else '')
        axs.append(ax)
    ax = fig.add_subplot(r1[0, 2])
    b = pd.read_csv(f'{S1D}/by_nu.csv')
    x = np.arange(6)
    for m in ('CALVADOS-2', 'Mpipi-Recharged'):
        t = b[b.model == m]
        ax.plot(x, t.frozen, '-o', ms=2.6, color=COL[m], lw=0.9)
        ax.plot(x, t.S1, '--^', ms=2.8, color=COL[m], lw=0.9)
    ax.axhline(0, color='0.75', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(NU_BIN_LAB, rotation=40, ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.4, 5.4)
    ax.set_ylim(-1.1, 0.05)
    ax.set_xlabel(f'Compaction, {NU}')
    ax.set_ylabel('Median total response (% per K)')
    ax.legend(handles=[plt.Line2D([], [], color='0.4', marker='o', ms=2.6, lw=0.9, label='calibrated'),
                       plt.Line2D([], [], color='0.4', marker='^', ms=2.8, ls='--', lw=0.9, label=CG),
                       matplotlib.patches.Patch(color=COL['CALVADOS-2'], label='CALVADOS-2'),
                       matplotlib.patches.Patch(color=COL['Mpipi-Recharged'], label='Mpipi-Recharged')],
              loc='lower right', handlelength=1.8)
    axs.append(ax)
    ax = fig.add_subplot(r2[0, 0])
    a = pd.read_csv(f'{S1D}/agreement.csv')
    lab = {'CALVADOS-2 vs Mpipi-Recharged': 'C2–MR', 'CALVADOS-2 vs Mpipi-GG': 'C2–GG', 'Mpipi-GG vs Mpipi-Recharged': 'GG–MR'}
    order = ['CALVADOS-2 vs Mpipi-GG', 'CALVADOS-2 vs Mpipi-Recharged', 'Mpipi-GG vs Mpipi-Recharged']
    a = a.set_index('pair').loc[order]
    for i, (p, r) in enumerate(a.iterrows()):
        y = 2 - i
        ax.plot([r.frozen, r.S1], [y, y], color='0.85', lw=2, zorder=1, solid_capstyle='butt')
        ax.scatter(r.frozen, y, s=18, color='white', edgecolor='0.45', lw=0.8, zorder=3)
        ax.scatter(r.S1, y, s=20, color=INK, lw=0, zorder=4)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels([lab[p] for p in order])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Agreement between force fields\n(Spearman)')
    ax.set_xlim(0.75, 0.93)
    ax.set_ylim(-0.5, 2.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', mfc='white', mec='0.45', ms=3.6, label='calibrated'),
                       plt.Line2D([], [], marker='o', ls='', color=INK, ms=3.6, label=CG)],
              loc='lower left', bbox_to_anchor=(-0.04, 1.0), ncol=2, handletextpad=0.2, borderaxespad=0)
    axs.append(ax)
    ax = fig.add_subplot(r2[0, 1])
    pp = pd.read_csv(f'{S1D}/partials.csv')
    labd = {'f_DEK': 'D + E + K', 'f_pi': 'π fraction', 'hydropathy_avg': 'hydropathy', 'f_K': 'Lys'}
    descs = list(labd)
    x = np.arange(len(descs))
    for o, m in zip((-0.2, 0.2), ('CALVADOS-2', 'Mpipi-Recharged')):
        t = pp[pp.model == m].set_index('descriptor').reindex(descs)
        ax.bar(x + o - 0.09, t.frozen, 0.17, color='white', edgecolor=COL[m], lw=0.8)
        ax.bar(x + o + 0.09, t.S1, 0.17, color=COL[m])
    ax.axhline(0, color='0.6', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([labd[d_] for d_ in descs])
    ax.tick_params(axis='x', length=0)
    ax.set_ylim(-0.42, 0.75)
    ax.set_ylabel(f'Partial correlation with\nstronger compaction\non heating (given {NU} and length)')
    ax.legend(handles=[matplotlib.patches.Patch(facecolor='white', edgecolor='0.4', lw=0.8, label='calibrated'),
                       matplotlib.patches.Patch(color='0.4', label=CG),
                       matplotlib.patches.Patch(color=COL['CALVADOS-2'], label='CALVADOS-2'),
                       matplotlib.patches.Patch(color=COL['Mpipi-Recharged'], label='Mpipi-Recharged')],
              ncol=4, loc='upper right', handlelength=1.0, columnspacing=1.2)
    axs.append(ax)
    add_letters(fig, [(axs[0], 'a', {'row': 1}), (axs[1], 'b', {'row': 1}), (axs[2], 'c', {'row': 1}),
                      (axs[3], 'd', {'row': 2}), (axs[4], 'e', {'row': 2})])
    save(fig, out)
    return out_numbers



CVL = f'{FDT}/total_response/capped_vs_linear'


def si7_capped(out):
    """Proteome results with the linear calibration against the capped main model (capped_vs_linear.py, 2026-10-25)."""
    from scipy import stats
    fig = plt.figure(figsize=(WIDTH, 4.8))
    r1 = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05], left=0.085, right=0.985, top=0.935, bottom=0.6, wspace=0.5)
    r2 = fig.add_gridspec(1, 2, width_ratios=[1, 1.75], left=0.085, right=0.985, top=0.4, bottom=0.11, wspace=0.42)
    lo, hi = -3.7, 0.7                                       # every chain inside the axes
    out_numbers = {}
    axs = []
    for j, (k, m) in enumerate((('C', 'CALVADOS-2'), ('R', 'Mpipi-Recharged'))):
        ax = fig.add_subplot(r1[0, j])
        d = pd.read_csv(f'{CVL}/chains_{k}.csv')
        assert d[['total_linear', 'total_capped']].min().min() > lo and d[['total_linear', 'total_capped']].max().max() < hi
        ax.hexbin(d.total_linear, d.total_capped, gridsize=44, bins='log', cmap=ff_cmap(COL[m]), mincnt=1, linewidths=0,
                  extent=(lo, hi, lo, hi))
        ax.plot([lo, hi], [lo, hi], color='0.55', lw=0.6, ls='--', zorder=0)
        ax.axhline(0, color='0.75', lw=0.5, zorder=0)
        rho = stats.spearmanr(d.total_linear, d.total_capped).statistic
        pc = 100 * (d.total_capped < 0).mean()
        out_numbers[m] = (round(rho, 2), round(pc, 1), len(d))
        ax.set_title(m, loc='left', color=COL[m])
        ax.text(0.97, 0.04, f'ρ = {rho:.2f}\n{pc:.1f}% contract\n(calibrated model)', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=6.5)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xticks([-3, -2, -1, 0])
        ax.set_yticks([-3, -2, -1, 0])
        ax.set_aspect('equal')
        ax.set_xlabel('Linear form (% per K)')
        ax.set_ylabel('Calibrated model (% per K)' if j == 0 else '')
        axs.append(ax)
    ax = fig.add_subplot(r1[0, 2])
    b = pd.read_csv(f'{CVL}/by_nu.csv')
    x = np.arange(6)
    for m in ('CALVADOS-2', 'Mpipi-Recharged'):
        t = b[b.model == m]
        ax.plot(x, t.linear, '--o', ms=2.6, color=COL[m], lw=0.9)
        ax.plot(x, t.capped, '-^', ms=2.8, color=COL[m], lw=0.9)
    ax.axhline(0, color='0.75', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(NU_BIN_LAB, rotation=40, ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.4, 5.4)
    ax.set_ylim(-1.1, 0.05)
    ax.set_xlabel(f'Compaction, {NU}')
    ax.set_ylabel('Median total response (% per K)')
    ax.legend(handles=[plt.Line2D([], [], color='0.4', marker='o', ms=2.6, lw=0.9, ls='--', label='linear form'),
                       plt.Line2D([], [], color='0.4', marker='^', ms=2.8, lw=0.9, label='calibrated model'),
                       matplotlib.patches.Patch(color=COL['CALVADOS-2'], label='CALVADOS-2'),
                       matplotlib.patches.Patch(color=COL['Mpipi-Recharged'], label='Mpipi-Recharged')],
              loc='lower right', handlelength=1.8)
    axs.append(ax)
    ax = fig.add_subplot(r2[0, 0])
    a = pd.read_csv(f'{CVL}/agreement.csv')
    lab = {'CALVADOS-2 vs Mpipi-Recharged': 'C2–MR', 'CALVADOS-2 vs Mpipi-GG': 'C2–GG', 'Mpipi-GG vs Mpipi-Recharged': 'GG–MR'}
    order = ['CALVADOS-2 vs Mpipi-GG', 'CALVADOS-2 vs Mpipi-Recharged', 'Mpipi-GG vs Mpipi-Recharged']
    a = a.set_index('pair').loc[order]
    for i, (p, r) in enumerate(a.iterrows()):
        y = 2 - i
        ax.plot([r.linear, r.capped], [y, y], color='0.85', lw=2, zorder=1, solid_capstyle='butt')
        ax.scatter(r.linear, y, s=18, color='white', edgecolor='0.45', lw=0.8, zorder=3)
        ax.scatter(r.capped, y, s=20, color=INK, lw=0, zorder=4)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels([lab[p] for p in order])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Agreement between force fields\n(Spearman)')
    ax.set_xlim(0.75, 0.93)
    ax.set_ylim(-0.5, 2.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', mfc='white', mec='0.45', ms=3.6, label='linear form'),
                       plt.Line2D([], [], marker='o', ls='', color=INK, ms=3.6, label='calibrated model')],
              loc='lower left', bbox_to_anchor=(-0.04, 1.0), ncol=2, handletextpad=0.2, borderaxespad=0)
    axs.append(ax)
    ax = fig.add_subplot(r2[0, 1])
    pp = pd.read_csv(f'{CVL}/partials.csv')
    labd = {'f_DEK': 'D + E + K', 'f_pi': 'π fraction', 'hydropathy_avg': 'hydropathy', 'f_K': 'Lys'}
    descs = list(labd)
    x = np.arange(len(descs))
    for o, m in zip((-0.2, 0.2), ('CALVADOS-2', 'Mpipi-Recharged')):
        t = pp[pp.model == m].set_index('descriptor').reindex(descs)
        ax.bar(x + o - 0.09, t.linear, 0.17, color='white', edgecolor=COL[m], lw=0.8)
        ax.bar(x + o + 0.09, t.capped, 0.17, color=COL[m])
    ax.axhline(0, color='0.6', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([labd[d_] for d_ in descs])
    ax.tick_params(axis='x', length=0)
    v = pp[['linear', 'capped']].to_numpy()                  # every bar inside, legend clear of the bars
    ax.set_ylim(min(-0.42, v.min() - 0.08), max(0.75, v.max() + 0.3))
    ax.set_ylabel(f'Partial correlation with\nstronger compaction\non heating (given {NU} and length)')
    ax.legend(handles=[matplotlib.patches.Patch(facecolor='white', edgecolor='0.4', lw=0.8, label='linear form'),
                       matplotlib.patches.Patch(color='0.4', label='calibrated model'),
                       matplotlib.patches.Patch(color=COL['CALVADOS-2'], label='CALVADOS-2'),
                       matplotlib.patches.Patch(color=COL['Mpipi-Recharged'], label='Mpipi-Recharged')],
              ncol=4, loc='upper right', handlelength=1.0, columnspacing=1.2)
    axs.append(ax)
    add_letters(fig, [(axs[0], 'a', {'row': 1}), (axs[1], 'b', {'row': 1}), (axs[2], 'c', {'row': 1}),
                      (axs[3], 'd', {'row': 2}), (axs[4], 'e', {'row': 2})])
    save(fig, out)
    return out_numbers

if __name__ == '__main__':
    # Figs 2, 3, the Fig. 5 sequence tracks (also drawn inside fig_apps.py) and Fig. S15; Fig. 1 is built by
    # fig_manuscript2.fig1, Fig. 4 by fig_heldout.py, Fig. 5 by fig_apps.py and Fig. 6 by fig_page4.py.
    #   BENDER_REANALYSIS_DIR=... python fig_extra.py       (openmpipi env)
    matplotlib.use('Agg')
    house()
    F = f'{ROOT}/manuscript/figs'
    fig2(f'{F}/fig_hydration')
    print('fig_total', fig3(f'{F}/fig_total'))
    fig5_tracks(f'{F}/fig_tracks')
    print('fig_si7', si7_capped(f'{F}/fig_si7'))
    print('written fig_hydration, fig_total, fig_tracks, fig_si7')
