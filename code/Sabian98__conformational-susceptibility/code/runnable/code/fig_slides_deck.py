"""
fig_slides_deck.py — 16:9 talk slides on the BENDER-MT paper (2026-10-30), not used in the manuscript.
  00 the larger scheme              tree of life to one ensemble to what it reveals (counts from the bio table, rest cartoon)
  01 forcing, how it works           steps 1-2 cartoons, steps 3-5 from the real CALVADOS-2 300 K ensemble of Q02021
  02 forcing on one real chain       Q02021 frames, per-frame forcing against Rg, the response split into its parts
  03 calibrating and testing it      Wuttke 2014 heating curves, left-out slopes, error on the series not fitted
  04 findings overview               four cards, each with real data
  05 one ensemble suffices           direct/predicted change for eight observables, force-field agreement
  06 compaction and charge           heat sensitivity by compaction, composition grid, enrichment of residues
  07 signals and substitutions       phosphorylation, salt and single substitutions (exact reweighting, CALVADOS-2)
  08 PAGE4                           size change against SAXS, heat sensitivity of each phosphoform in direct runs
  09 growth temperature              slopes per 10 C in SD, thermophile minus mesophile orthologues
  10 take-home and limitations
Every plotted value is read from the result files behind the manuscript figures (fig_forcing, fig_extra, fig_heldout,
fig_apps, fig_page4, fig_manuscript2); only the chains of slide 1 steps 1-2 and slide 4 card 3 are cartoons.
Script colours are the original figure colours, which figfont maps to the state-light palette at save time.
    BENDER_REANALYSIS_DIR=$CS_REPO/bender_reanalysis python fig_slides_deck.py      (openmpipi env)
Output: BENDER-MT/slides/slide_NN_<name>.png and BENDER-MT/slides/BENDER-MT_slides.pptx (all slides, one full-bleed\npicture each, written by png_to_pptx.py)
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Patch, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.expanduser((__import__('os').environ['CS_REPO'] + '/bender_reanalysis'))))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402
import figfont  # noqa: E402

OUT = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/slides')
ROOT = FX.ROOT
TA = f'{ROOT}/Thermal adaptation'
# original figure colours (figfont maps them to state-light at save time)
C2, GG, MR = FX.COL['CALVADOS-2'], FX.COL['Mpipi-GG'], FX.COL['Mpipi-Recharged']
HYD, TEA, BASIC, ACID, PHOS, PALE = '#B5473A', '#C8792A', '#7C8BA0', '#B5473A', '#C8792A', '#A8CBE8'
INK, MUTED, CARD, BAND = '#222222', '#6B6B6B', '#F3F4F6', '#EAF2F9'
BEAD, BOND, WATER = '#BDBDBD', '#9E9E9E', '#9CC9E8'
W, H, SX, SY = 13.333, 7.5, 133.33, 75.0


# ----------------------------------------------------------------------------- slide furniture
def new_slide(title, sub=None):
    fig = plt.figure(figsize=(W, H))
    fig.slide_title = title
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, SX)
    ax.set_ylim(0, SY)
    ax.axis('off')
    ax.text(4, 70.2, title, fontsize=22, fontweight='bold', color=INK, va='center')
    if sub:
        ax.text(4, 65.4, sub, fontsize=13, color=MUTED, va='center')
    return fig, ax


def card(ax, x, y, w, h, fc=CARD, r=1.6, z=0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec='none', zorder=z))


def takeaway(ax, text):
    card(ax, 3.5, 2.6, 126.3, 8.0, fc=BAND)
    ax.text(66.7, 6.6, text, ha='center', va='center', fontsize=14.5, fontweight='bold', color='#537A9E')


def pax(fig, x, y, w, h, title=None):
    """Data axes at slide coordinates (x, y, w, h), with a bold panel title above it."""
    a = fig.add_axes([x / SX, y / SY, w / SX, h / SY])
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
    a.tick_params(labelsize=12)
    if title:
        a.set_title(title, loc='left', fontsize=14, fontweight='bold', color=INK, pad=12)
    return a


def arrow(ax, p0, p1, col='#8A8A8A', lw=2.2, ms=16, z=3):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=ms, lw=lw, color=col, zorder=z, shrinkA=0,
                                 shrinkB=0))


def badge(ax, x, y, n):
    ax.add_patch(Circle((x, y), 1.25, fc=INK, ec='none', zorder=5))
    ax.text(x, y - 0.05, str(n), ha='center', va='center', color='white', fontsize=12, fontweight='bold', zorder=6)


def note(ax, x, y, text, size=11.5, ha='center', col=INK):
    ax.text(x, y, text, ha=ha, va='center', fontsize=size, color=col, linespacing=1.3)


WB = dict(fc='white', ec='none', alpha=0.9, pad=1.5)


# ----------------------------------------------------------------------------- cartoon chains (slides 1 and 4)
def walk(n, radius, persist, seed):
    rng = np.random.default_rng(seed)
    for _ in range(500):
        xy, ang, ok = [np.zeros(2)], rng.uniform(0, 2 * np.pi), True
        for i in range(1, n):
            for _ in range(300):
                a = ang + rng.normal(0, persist)
                p = xy[-1] + np.array([np.cos(a), np.sin(a)])
                if np.linalg.norm(p) > radius:
                    continue
                if i > 2 and np.min(np.linalg.norm(np.array(xy[:-2]) - p, axis=1)) < 0.85:
                    continue
                xy.append(p)
                ang = a
                break
            else:
                ok = False
                break
        if ok:
            c = np.array(xy)
            return c - c.mean(0)
    raise RuntimeError('no chain found')


def n_contacts(xy, ch, cut):
    return sum(1 for i in range(len(xy)) for j in range(i + 3, len(xy))
               if ch.get(i) and ch.get(j) and np.linalg.norm(xy[i] - xy[j]) < cut)


def draw_chain(ax, xy, centre, scale, ch, bead=0.45, cut=1.55, phos=()):
    p = xy * scale + np.asarray(centre)
    ax.plot(p[:, 0], p[:, 1], '-', color=BOND, lw=1.6, zorder=3, solid_capstyle='round')
    for i in range(len(p)):
        for j in range(i + 3, len(p)):
            if ch.get(i) and ch.get(j) and np.linalg.norm(xy[i] - xy[j]) < cut:
                ax.plot(*p[[i, j]].T, ':', color=INK, lw=1.3, zorder=4)
    for i, q in enumerate(p):
        ax.add_patch(Circle(q, bead * scale, fc={1: BASIC, -1: ACID}.get(ch.get(i, 0), BEAD), ec='white', lw=0.6,
                            zorder=5))
        if i in phos:
            c = q + np.array([0.75, 0.75]) * scale
            ax.add_patch(Circle(c, 0.62 * scale, fc=PHOS, ec='white', lw=0.6, zorder=6))
            ax.text(*c, 'P', ha='center', va='center', fontsize=8, fontweight='bold', color='white', zorder=7)


def best_walk(n, radius, persist, ch, cut, seeds=60):
    best = None
    for s in range(seeds):
        c = walk(n, radius, persist, s)
        k = n_contacts(c, ch, cut)
        if best is None or k > best[0]:
            best = (k, c)
    return best[1]


def ion_pair(ax, cx, cy, d, released, rng):
    for dx, c, s in ((-d / 2, BASIC, '+'), (d / 2, ACID, '−')):
        ax.add_patch(Circle((cx + dx, cy), 1.5, fc=c, ec='white', lw=1, zorder=4))
        ax.text(cx + dx, cy - 0.1, s, ha='center', va='center', color='white', fontsize=14, fontweight='bold', zorder=5)
        for a in np.linspace(0, 2 * np.pi, 9 if not released else 5, endpoint=False) + rng.uniform(0, 0.4):
            v = np.array([np.cos(a), np.sin(a)])
            if released and np.sign(v[0]) == np.sign(-dx) and abs(v[1]) < 0.75:
                continue
            ax.add_patch(Circle((cx + dx + 2.45 * v[0], cy + 2.45 * v[1]), 0.48, fc=WATER, ec='none', zorder=3))


# ----------------------------------------------------------------------------- data
def forcing_example():
    """Q02021 in CALVADOS-2 at 300 K, exactly as manuscript Fig. 2 (fig_forcing.py, checked against stored values)."""
    import fig_forcing as F
    c = F.compute()
    _, _, _, g = F.check(c)
    val, se = F.response_se(c, g)
    seq = str(c['seq'])
    charged = np.array([r in F.HYD for r in seq])
    return dict(c=c, g=g, val=val, se=se, seq=seq, charged=charged, F=F)


def heldout():
    import fig_heldout as HO
    data = {m: HO.load(m, lab) for m, lab, _ in HO.MODELS}
    nf = pd.read_csv(HO.NF)
    out = {}
    for m, lab, _ in HO.MODELS:
        for v in HO.VARIANTS:
            out[(lab, v)] = HO.stats(data[m], v)
        out[(lab, 'floor')] = nf[(nf.model == lab) & nf.series.isin(HO.NOT_FITTED)].expected_noise_miss.mean()
    return out


def page4():
    import fig_page4 as P
    rw, dr, cv = P.load()
    hs, rg = {}, {}
    for m in ('CALVADOS-2', 'Mpipi-Recharged'):
        d = dr[m]
        if m == 'CALVADOS-2':        # results.tex: fluctuation-response on each direct 300 K run
            hs[m] = [(-d.loc[f, 'Rg_total'], d.loc[f, 'Rg_total_se']) for f in P.DIRECT]
        else:                        # results.tex: slope at 300 K of the heating curves simulated with the forcings
            c = cv[m][cv[m].curve == 'calibrated'].set_index('form')
            hs[m] = [(-c.loc[f, 'slope300'], c.loc[f, 'slope300_se']) for f in P.DIRECT]
        rg[m] = [P.pct(d.loc[f, 'Rg_mean'], d.loc[f, 'Rg_se'], d.loc['WT', 'Rg_mean'], d.loc['WT', 'Rg_se'])
                 for f in P.DIRECT[1:]]
    saxs = [(100 * (P.SAXS[k] / P.SAXS[0] - 1), 100 * P.SAXS_SE[k] / P.SAXS[0]) for k in (1, 2)]
    return dict(hs=hs, rg=rg, saxs=saxs)


def ogt():
    sl = pd.read_csv(f'{TA}/bender_ogt/slopes.csv')
    sl = sl[(sl.quantity == 'total') & (sl.fixed_effects == 'phylum (bacteria only)')].set_index('model')
    fc = pd.read_csv(f'{TA}/bender_ogt/followup_composition.csv')
    fc = fc[(fc.model == 'CALVADOS-2') & fc.analysis.str.startswith('descriptor')].set_index('term')
    rows = []
    for t, lab in (('f_pi', 'π residues'), ('f_R', 'Arg')):
        rows.append((lab, fc.loc[t, 'slope_per_10C'], fc.loc[t, 'ci_lo'], fc.loc[t, 'ci_hi'], INK))
    for m, col in (('CALVADOS-2', C2), ('Mpipi-Recharged', MR)):
        r, sd = sl.loc[m], sl.loc[m, 'sd_across_organisms']
        # heat sensitivity is minus the total response
        rows.append((f'heat sensitivity\n{m}', -r.slope_per_10C / sd, -r.ci_hi / sd, -r.ci_lo / sd, col))
    A = pd.read_csv(f'{TA}/predictions/A_matched.csv')
    diff = A[A['class'] == 'thermophile'].diff_total_C.to_numpy()
    return rows, diff


def load_all():
    import fig_apps as FA
    D, bio = FX._load_totals()
    S, Cs = FA.crosses()
    return dict(fx=forcing_example(), ho=heldout(), p4=page4(), ogt=ogt(), D=D, bio=bio, S=S, Cs=Cs, FA=FA,
                byn=pd.read_csv(f'{FX.TR}/A_by_nu.csv'), agree=pd.read_csv(f'{FX.TR}/B_agreement.csv'),
                grid=pd.read_csv(f'{FX.TR}/E_grid_sharedCR.csv'),
                mut=pd.read_csv(f'{FX.TR}/mutation_scan/analysis_by_mut.csv', index_col=0).sort_values('median_change'))


# ----------------------------------------------------------------------------- shared real-data panels
def saturation_panel(a, X, size=11.5):
    import capped_proteome as CP
    cn = X['c']['c'][:, X['charged']].ravel()
    xc = np.linspace(0, 7, 200)
    hist, edges = np.histogram(cn, bins=np.linspace(0, 7, 36), density=True)
    a.fill_between(edges[:-1], 0, hist / hist.max() * 1.6, step='post', color=PALE, alpha=0.8, lw=0)
    a.plot(xc, xc, '--', color=BOND, lw=2)
    a.plot(xc, CP.h(xc), '-', color=HYD, lw=3)
    a.set(xlim=(0, 7), ylim=(0, 4.6), xticks=[0, 2, 4, 6], yticks=[0, 2, 4])
    a.set_xlabel('Contacts of a residue', fontsize=size)
    a.set_ylabel('Contacts counted', fontsize=size)
    a.text(4.4, 3.7, 'linear form', fontsize=size - 1, color=MUTED)
    a.text(4.9, 1.4, 'saturating', fontsize=size - 0.5, color=HYD, ha='center', fontweight='bold')
    a.text(6.9, 0.25, 'contact numbers\nof charged residues', fontsize=size - 2, color='#537A9E', ha='right')


def frame_scatter(a, X, size=11.5, label=True):
    c, g = X['c'], X['g']
    rg, Dh = c['rg'], c['D_hyd'] * g['C']
    a.scatter(rg, Dh, s=7, color=C2, alpha=0.45, lw=0)
    p = np.polyfit(rg, Dh, 1)
    xx = np.percentile(rg, [1, 99])
    a.plot(xx, np.polyval(p, xx), color=INK, lw=2)
    a.set_xlabel(f'{FX.RG} of a frame (nm)', fontsize=size)
    a.set_ylabel('Hydrophilic forcing\n(kJ mol$^{-1}$ K$^{-1}$)', fontsize=size)
    if label:
        r = np.corrcoef(rg, Dh)[0, 1]
        a.text(0.04, 0.95, f'r = {r:.2f}\n{len(rg):,} frames', transform=a.transAxes, va='top', fontsize=size - 0.5)


def parts_panel(a, X, size=11.5, horizontal=False):
    keys = ['chi', 'TEA', 'H', 'total']
    lab = ['fixed\nsolvent', 'hydro-\nphobic', 'hydro-\nphilic', 'total']
    v = [X['val'][k] for k in keys]
    e = [X['se'][k] for k in keys]
    if horizontal:
        y = np.arange(4)[::-1]
        a.barh(y, v, xerr=e, color=[BEAD, TEA, HYD, INK], height=0.62, error_kw=dict(lw=1.1, capsize=2.5))
        a.axvline(0, color=INK, lw=0.9)
        a.set_yticks(y)
        a.set_yticklabels(['fixed solvent', 'hydrophobic', 'hydrophilic', 'total'], fontsize=size - 1)
        a.tick_params(axis='y', length=0)
        a.spines['left'].set_visible(False)
        a.set_xlabel('Response of $R_\\mathrm{g}$ (% per K)', fontsize=size)
        a.set_xlim(min(v) - 0.12, max(v) + 0.12)
        return
    a.bar(range(4), v, yerr=e, color=[BEAD, TEA, HYD, INK], width=0.66, error_kw=dict(lw=1.2, capsize=3))
    a.axhline(0, color=INK, lw=0.9)
    a.set_xticks(range(4))
    a.set_xticklabels(lab, fontsize=size - 1.5)
    a.tick_params(axis='x', length=0, pad=4)
    a.spines['bottom'].set_visible(False)
    a.set_ylabel('Response of $R_\\mathrm{g}$\n(% per K)', fontsize=size)
    for k, (vv, ee) in enumerate(zip(v, e)):
        a.text(k, vv + (ee + 0.03 if vv > 0 else -ee - 0.03), FX.minus(f'{vv:+.2f}' if abs(vv) >= 0.01 else f'{vv:+.3f}'),
               ha='center',
               va='bottom' if vv > 0 else 'top', fontsize=size - 1.5)
    a.set_ylim(min(v) - 0.17, max(v) + 0.17)


def byn_panel(a, byn, models=('CALVADOS-2', 'Mpipi-GG', 'Mpipi-Recharged'), size=12, legend=True):
    x = np.arange(6)
    w = 0.8 / len(models)
    for k, m in enumerate(models):
        d = byn[(byn.model == m) & (byn.own_nu_bin != 'all')]
        xx = x + (k - (len(models) - 1) / 2) * w
        a.bar(xx, -d.total_median, width=w * 0.92, color=FX.COL[m], label=m)
        a.errorbar(xx, -d.total_median, yerr=[d.total_q75 - d.total_median, d.total_median - d.total_q25], fmt='none',
                   ecolor=INK, lw=0.9, capsize=0)
    a.set_xticks(x)
    a.set_xticklabels(['≤0.45', '0.50', '0.54', '0.58', '0.63', '>0.63'], fontsize=size - 1)
    a.tick_params(axis='x', length=0)
    a.set_xlabel(f'Compaction, {FX.NU} at 300 K (upper bin edge)', fontsize=size)
    a.set_ylabel('Heat sensitivity (% per K)', fontsize=size)
    if legend:
        a.legend(loc='upper right', fontsize=size - 1.5, frameon=False, handlelength=1.0)


def phospho_hex(a, S, size=12, gs=(34, 40)):
    from scipy import stats
    x, y = S.phospho_dRg_pct, S.phospho_dtotal
    ext = (-12.6, 11.0, -0.62, 0.45)
    a.hexbin(x, -y, gridsize=gs, bins='log', cmap=FX.ff_cmap(C2), mincnt=1, linewidths=0,
             extent=(ext[0], ext[1], -ext[3], -ext[2]))
    a.axhline(0, color='0.55', lw=0.6)
    a.axvline(0, color='0.55', lw=0.6)
    rho = stats.spearmanr(x, y)[0]
    a.set(xlim=ext[:2], ylim=(-ext[3], -ext[2]))
    a.set_xlabel(f'Change of {FX.RG} (%)', fontsize=size)
    a.set_ylabel('Change of heat\nsensitivity (% per K)', fontsize=size)
    a.text(0.97, 0.95, FX.minus(f'ρ = {-rho:.2f}') + f'\n{len(S):,} sites', transform=a.transAxes, ha='right', va='top', fontsize=size - 1)
    return rho


# ----------------------------------------------------------------------------- 01 forcing overview
def s01(d):
    X = d['fx']
    fig, ax = new_slide('How the hydration forcing turns one ensemble into a heat response',
                        'Coarse-grained force fields have no temperature-dependent solvent. We add it to the 300 K '
                        'simulations as an energy term that is exactly zero at 300 K.')
    xs, cw, y0, ch = [3.5, 29.8, 56.1, 82.4, 108.7], 23.6, 13.5, 47
    titles = ['One simulated\nensemble at 300 K', 'Heating makes charged\ncontacts stickier',
              'Each residue saturates\nat about two contacts', 'Compact frames are\nfavoured most',
              'Fluctuation–response\ngives the response']
    for k, (x, t) in enumerate(zip(xs, titles)):
        card(ax, x, y0, cw, ch)
        badge(ax, x + 2.6, y0 + ch - 3.2, k + 1)
        ax.text(x + 4.9, y0 + ch - 3.2, t, fontsize=12.5, fontweight='bold', color=INK, va='center', linespacing=1.15)
        if k < 4:
            arrow(ax, (x + cw + 0.3, y0 + ch / 2), (xs[k + 1] - 0.3, y0 + ch / 2))
    CHG = {1: 1, 4: -1, 7: -1, 9: 1, 12: -1, 15: 1, 18: -1, 21: 1, 23: -1, 26: 1, 29: -1, 32: 1, 34: -1}
    x = xs[0]
    draw_chain(ax, best_walk(36, 3.1, 1.3, CHG, 1.55), (x + 6.6, y0 + 27.0), 1.35, CHG, bead=0.42)
    draw_chain(ax, walk(36, 6.5, 0.42, 11), (x + 17.4, y0 + 27.0), 0.82, CHG, bead=0.6)
    note(ax, x + 6.6, y0 + 17.4, 'compact', 11, col=MUTED)
    note(ax, x + 17.4, y0 + 17.4, 'expanded', 11, col=MUTED)
    note(ax, x + cw / 2, y0 + 9.5, 'Thousands of stored frames\nof each IDP, in CALVADOS-2,\nMpipi-GG and Mpipi-Recharged',
         10.5)
    for c, lab, dx in ((BASIC, 'Lys, Arg', 0), (ACID, 'Asp, Glu', 9)):
        ax.add_patch(Circle((x + 5.4 + dx, y0 + 3.2), 0.6, fc=c, ec='none'))
        ax.text(x + 6.5 + dx, y0 + 3.2, lab, va='center', fontsize=10, color=INK)
    x = xs[1]
    rng = np.random.default_rng(4)
    ion_pair(ax, x + cw / 2, y0 + 31.5, 8.0, False, rng)
    arrow(ax, (x + cw / 2, y0 + 26.6), (x + cw / 2, y0 + 22.4), col=HYD, lw=2.4)
    ax.text(x + cw / 2 + 1.2, y0 + 24.5, 'heat', fontsize=11, color=HYD, va='center', fontweight='bold')
    ion_pair(ax, x + cw / 2, y0 + 16.8, 3.0, True, rng)
    for dx, dy in ((-6.5, 20.4), (6.2, 20.6), (7.4, 13.2), (-7.6, 13.0)):
        ax.add_patch(Circle((x + cw / 2 + dx, y0 + dy), 0.48, fc=WATER, ec='none', zorder=3))
    note(ax, x + cw / 2, y0 + 7.2, 'Released water favours the\ncontact more as temperature\nrises (Wuttke et al. 2014)', 10.5)
    a = inset(fig, [xs[2] + 5.6, y0 + 15.5, cw - 7.6, 22.0])
    saturation_panel(a, X, 10.5)
    _clean(a)
    note(ax, xs[2] + cw / 2, y0 + 7.0, 'A buried residue has little\nwater left to release, so it\nadds little beyond two contacts',
         10.5)
    a = inset(fig, [xs[3] + 7.4, y0 + 15.5, cw - 9.2, 22.0])
    frame_scatter(a, X, 10.5, label=False)
    _clean(a)
    note(ax, xs[3] + cw / 2, y0 + 7.0, 'Compact frames have more\ncharged contacts, so heating\nfavours them most',
         10.5)
    a = inset(fig, [xs[4] + 9.6, y0 + 16.5, cw - 11.0, 20.5])
    parts_panel(a, X, 10.5, horizontal=True)
    _clean(a)
    note(ax, xs[4] + cw / 2, y0 + 7.0, 'Negative means compaction\non heating. Its size is the\nheat sensitivity of the IDP',
         10.5)
    ax.text(66.7, y0 - 1.6, 'Steps 3 to 5 show the Gly-rich protein from tomato (UniProt Q02021) in CALVADOS-2 at 300 K, '
            'as in Fig. 3 of the paper', ha='center', va='center', fontsize=10, color=MUTED, style='italic')
    takeaway(ax, 'Every simulated ensemble already contains its response to heat, once hydration is added')
    return fig


def inset(fig, r):
    """Figure-level axes at slide coordinates [x, y, w, h], so that each plot exports as its own picture."""
    return fig.add_axes([r[0] / SX, r[1] / SY, r[2] / SX, r[3] / SY])


def _clean(a):
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
    a.patch.set_alpha(0)
    a.tick_params(labelsize=10)


# ----------------------------------------------------------------------------- 02 forcing on one real chain
def s02(d):
    X = d['fx']
    c, F = X['c'], X['F']
    rg, Dh = c['rg'], c['D_hyd'] * X['g']['C']
    fig, ax = new_slide('The forcing on one real chain',
                        'Gly-rich protein from tomato (UniProt Q02021), 1,000 stored frames of CALVADOS-2 at 300 K')
    k_c, k_e = int(np.argsort(rg)[len(rg) // 20]), int(np.argsort(rg)[-len(rg) // 20])
    a = pax(fig, 4, 20, 52, 33, 'A compact and an expanded frame')
    a.set_title(a.get_title(loc='left'), loc='left', fontsize=14, fontweight='bold', color=INK, pad=30)
    a.set_axis_off()
    off, ymin = 0.0, np.inf
    labels = []
    for k, lab in ((k_c, 'compact frame'), (k_e, 'expanded frame')):
        x = c['xyz'][k].astype(float)
        x = x - x.mean(0)
        _, _, vt = np.linalg.svd(x, full_matrices=False)
        p = x @ vt[:2].T
        p[:, 0] += off - p[:, 0].min()
        a.plot(p[:, 0], p[:, 1], color=BOND, lw=1.0, zorder=1)
        hy = X['charged']
        sign = np.array([1 if r in 'KR' else -1 if r in 'DE' else 0 for r in X['seq']])
        a.scatter(p[~hy, 0], p[~hy, 1], s=26, color=BEAD, edgecolor='white', lw=0.4, zorder=3)
        a.scatter(p[sign > 0, 0], p[sign > 0, 1], s=34, color=BASIC, edgecolor='white', lw=0.4, zorder=4)
        a.scatter(p[sign < 0, 0], p[sign < 0, 1], s=34, color=ACID, edgecolor='white', lw=0.4, zorder=4)
        ci, cj = F.contacts(x, c)
        both = hy[ci] & hy[cj]
        for u, v in zip(ci[both], cj[both]):
            a.plot(p[[u, v], 0], p[[u, v], 1], color=INK, lw=0.9, ls=(0, (1.2, 1.2)), zorder=2)
        labels.append(((p[:, 0].min() + p[:, 0].max()) / 2,
                       f'{lab}\n{FX.RG} = {rg[k]:.2f} nm, {both.sum()} charged contacts\n'
                       f'forcing {FX.minus(f"{Dh[k]:.2f}")} kJ mol$^{{-1}}$ K$^{{-1}}$'))
        ymin = min(ymin, p[:, 1].min())
        off = p[:, 0].max() + 4.0
    for mid, txt in labels:
        a.text(mid, ymin - 0.5, txt, ha='center', va='top', fontsize=11.5, color=INK, linespacing=1.3)
    a.set_aspect('equal', adjustable='datalim')
    a.autoscale_view()
    for k_, (col, lab) in enumerate(((BASIC, 'Lys, Arg'), (ACID, 'Asp, Glu'), (BEAD, 'other residues'))):
        a.scatter([0.02 + 0.2 * k_], [1.02], s=60, color=col, transform=a.transAxes, clip_on=False)
        a.text(0.04 + 0.2 * k_, 1.02, lab, transform=a.transAxes, va='center', fontsize=11.5)
    a2 = pax(fig, 69, 18, 25, 37, 'Forcing of every frame')
    frame_scatter(a2, X, 12.5)
    a2.scatter([rg[k_c], rg[k_e]], [Dh[k_c], Dh[k_e]], s=110, color=[ACID, BASIC], edgecolor='white', lw=1.2, zorder=5)
    a3 = pax(fig, 103, 18, 26, 37, 'Response of this chain')
    parts_panel(a3, X, 12.5)
    takeaway(ax, 'Compact frames carry more charged contacts, so the forcing turns slight expansion into compaction')
    return fig


# ----------------------------------------------------------------------------- 03 calibration and test
def s03(d):
    fig, ax = new_slide('Calibrating and testing the forcing',
                        'One strength per force field, fitted to smFRET heating curves (Wuttke et al. 2014), then '
                        'tested on series to which no parameter was fitted')
    for k, (n, t) in enumerate(FX.CONSTRUCTS[:3]):
        a = pax(fig, 6 + 17.5 * k, 21, 14.5, 32)
        FX.wuttke_panel(a, n, '', k == 0)
        a.set_title(t, loc='left', fontsize=12.5, color=MUTED, pad=6)
        a.set_ylim(0.84, 1.12)
        a.xaxis.label.set_size(12)
        a.yaxis.label.set_size(12)
        if k:
            a.set_yticklabels([])
    ax.text(6, 57.0, 'Heating curves of three calibration constructs', fontsize=14, fontweight='bold', color=INK)
    h = [plt.Line2D([], [], marker='o', ls='', color=INK, ms=6, label='measured'),
         plt.Line2D([], [], color=C2, lw=2.2, label='CALVADOS-2'),
         plt.Line2D([], [], color=MR, lw=2.2, label='Mpipi-Recharged'),
         plt.Line2D([], [], color='0.4', lw=2.2, label='calibrated'),
         plt.Line2D([], [], color='0.4', lw=1.4, ls=(0, (3, 2)), label='plain')]
    fig.legend(handles=h, loc='center', ncol=3, bbox_to_anchor=(31.5 / SX, 13.4 / SY), fontsize=11.5, frameon=False,
               handlelength=2.0, columnspacing=1.2)
    a = pax(fig, 64, 21, 24, 32, 'Construct left out of the fit')
    FX.loo_scatter(a)
    a.xaxis.label.set_size(12)
    a.yaxis.label.set_size(12)
    a.legend(handles=[plt.Line2D([], [], marker='o', ls='', color='0.3', ms=7, label='calibrated'),
                      plt.Line2D([], [], marker='o', ls='', mfc='white', color='0.3', ms=7, label='plain')],
             loc='lower right', fontsize=11.5, frameon=False)
    for coll in a.collections:
        coll.set_sizes([60])
    a = pax(fig, 99, 21, 30, 32, 'Five series not fitted')
    ho = d['ho']
    xs, labs = [], []
    for j, (lab, col) in enumerate((('CALVADOS-2', C2), ('Mpipi-Recharged', MR))):
        for k, v in enumerate(('plain', 'A0', 'capped')):
            mae, sign = ho[(lab, v)]
            x = j * 4 + k
            c_ = '0.6' if v == 'plain' else col
            a.bar(x, mae, color='white' if v == 'A0' else c_, edgecolor=c_, lw=1.6 if v == 'A0' else 0, width=0.8)
            a.text(x, mae + 0.01, f'{mae:.2f}', ha='center', va='bottom', fontsize=11)
            xs.append(x)
            labs.append({'plain': 'plain', 'A0': 'linear', 'capped': 'calibrated'}[v])
        fl = ho[(lab, 'floor')]
        a.plot([j * 4 - 0.45, j * 4 + 2.45], [fl, fl], color=INK, lw=1.4, ls=(0, (3, 1.5)))
        a.text(j * 4 + 1, 0.47, lab.replace('Mpipi-', 'Mpipi-\n'), ha='center', va='top', fontsize=11.5, color=col,
               fontweight='bold')
    a.text(3.5, ho[('CALVADOS-2', 'floor')] + 0.012, 'noise\nfloor', fontsize=10.5, ha='center', va='bottom', color=INK)
    a.set_xticks(xs)
    a.set_xticklabels(labs, rotation=40, ha='right', rotation_mode='anchor', fontsize=11.5)
    a.tick_params(axis='x', length=0)
    a.set(ylim=(0, 0.48), xlim=(-0.6, 7.2))
    a.set_ylabel('Mean absolute error (% per K)', fontsize=12)
    takeaway(ax, 'The calibrated model predicts the series not fitted within measurement noise')
    return fig


# ----------------------------------------------------------------------------- 04 findings overview
def s04(d):
    fig, ax = new_slide('What conformational susceptibility says about IDP biology',
                        '10,282 highly disordered IDPs from 2,156 organisms, three force fields, calibrated hydration. '
                        'Heat sensitivity is the percent drop in Rg per kelvin of heating.')
    xs, cw, y0, ch = [3.5, 36.3, 69.1, 101.9], 28.0, 13.5, 47
    titles = ['Nearly every IDP\ncompacts on heating', 'Signals retune heat\nsensitivity with compaction',
              'One input changes the\nresponse to another', 'Composition adapts to\ngrowth temperature']
    for k, (x, t) in enumerate(zip(xs, titles)):
        card(ax, x, y0, cw, ch)
        badge(ax, x + 2.6, y0 + ch - 3.2, k + 1)
        ax.text(x + 4.9, y0 + ch - 3.2, t, fontsize=12.5, fontweight='bold', color=INK, va='center', linespacing=1.15)
    a = inset(fig, [xs[0] + 6.0, y0 + 17.5, cw - 8.5, 19.5])
    byn_panel(a, d['byn'], ('CALVADOS-2',), 10.5, legend=False)
    a.set_xticks([0, 5])
    a.set_xticklabels(['compact', 'expanded'], fontsize=10.5)
    a.set_xlabel(f'Compaction bins of {FX.NU}', fontsize=10.5)
    _clean(a)
    note(ax, xs[0] + cw / 2, y0 + 8.4, '98.5, 97.1 and 98.2% of IDPs compact\non heating. Compaction and charge\n'
         'content set how much.', 10.5)
    a = inset(fig, [xs[1] + 6.5, y0 + 17.5, cw - 9.0, 19.5])
    phospho_hex(a, d['S'], 10.5, gs=(26, 30))
    _clean(a)
    note(ax, xs[1] + cw / 2, y0 + 8.4, 'One phosphosite each, CALVADOS-2.\nA change that compacts a chain\nraises its heat '
         'sensitivity.', 10.5)
    x = xs[2]
    chg = {0: 1, 1: 1, 2: 1, 3: 1, 14: -1, 15: -1, 16: -1, 17: -1, 18: -1}
    draw_chain(ax, best_walk(22, 3.0, 1.1, chg, 1.7, 80), (x + 6.4, y0 + 30.0), 1.25, chg, cut=1.7)
    draw_chain(ax, walk(22, 6.0, 0.4, 5), (x + 21.4, y0 + 30.0), 0.85, chg, bead=0.6, cut=1.7, phos=(5, 8, 11, 20))
    arrow(ax, (x + 11.8, y0 + 30.0), (x + 15.4, y0 + 30.0), col=PHOS, lw=2.4)
    ax.text(x + 13.6, y0 + 32.3, 'CLK2', ha='center', fontsize=10.5, color=PHOS, fontweight='bold')
    hs = d['p4']['hs']['CALVADOS-2']
    for xx, lab, v in ((x + 6.4, 'PAGE4', hs[0][0]), (x + 21.4, 'CLK2-PAGE4', hs[2][0])):
        note(ax, xx, y0 + 38.0, lab, 10.5, col=MUTED)
        ax.text(xx, y0 + 21.0, f'{v:.2f}', ha='center', fontsize=16, fontweight='bold', color=HYD)
    note(ax, x + cw / 2, y0 + 18.4, 'heat sensitivity (% per K, CALVADOS-2)', 9.5, col=MUTED)
    note(ax, x + cw / 2, y0 + 8.4, 'Predicted from one ensemble, then\nconfirmed by direct simulation.', 10.5)
    a = inset(fig, [xs[3] + 9.0, y0 + 17.5, cw - 11.0, 19.5])
    ogt_bars(a, d['ogt'][0], 10.5)
    _clean(a)
    note(ax, xs[3] + cw / 2, y0 + 8.4, '757 organisms. π residues and Arg\nrise with growth temperature, heat\n'
         'sensitivity barely changes.', 10.5)
    takeaway(ax, 'Compaction and charge link disordered ensembles to temperature sensing across the tree of life')
    return fig


def ogt_bars(a, rows, size=12):
    y = np.arange(len(rows))[::-1]
    for yy, (lab, est, lo, hi, col) in zip(y, rows):
        a.errorbar(est, yy, xerr=[[est - lo], [hi - est]], fmt='o', ms=7 if size > 11 else 5.5, color=col, lw=1.6,
                   capsize=0)
    a.axvline(0, color='0.6', lw=0.8)
    a.set_yticks(y)
    a.set_yticklabels([r[0] for r in rows], fontsize=size - 1)
    a.tick_params(axis='y', length=0)
    a.set_xlabel('Slope per 10 °C (SD)', fontsize=size)
    a.set_ylim(-0.6, len(rows) - 0.4)


# ----------------------------------------------------------------------------- 05 one ensemble suffices
def s05(d):
    fig, ax = new_slide('One ensemble predicts how an IDP responds to temperature',
                        'Integrated over 278 to 320 K, the susceptibility matches the directly simulated change, with 94 '
                        'to 95% of proteins within two standard errors')
    sc = pd.read_csv(f'{FX.FDT}/fdt_validation_summary.csv')
    sc = sc[sc['T'] == 'integral'].set_index('metric')
    sm = pd.read_csv(f'{FX.TLA}/Analysis_FDT_MPIPI/mpipi_fdt_validation_summary.csv')
    sm = sm[sm['T'] == 'integral'].set_index('metric')
    obs = [('rg', FX.RG), ('ree', 'End-to-end'), ('delta', 'Asphericity'), ('fragmentation_index', 'Fragmentation'),
           ('global_efficiency', 'Global efficiency'), ('avg_clustering', 'Clustering'), ('transitivity', 'Transitivity'),
           ('degree_assortativity', 'Assortativity')]
    rc = [sc.loc[m, 'mean_slope'] / sc.loc[m, 'mean_chi'] for m, _ in obs]
    rm = [sm.loc[m, 'ratio_observed_to_integral'] for m, _ in obs]
    a = pax(fig, 7, 20, 62, 35, 'Directly computed over predicted change, eight observables')
    x = np.arange(len(obs))
    a.bar(x - 0.2, rc, 0.38, color=C2, label=f'CALVADOS-2 ({int(sc.loc["rg", "n"]):,} chains)')
    a.bar(x + 0.2, rm, 0.38, color=GG, label=f'Mpipi-GG ({int(sm.loc["rg", "n"]):,} chains)')
    a.axhline(1, color=INK, lw=1.0)
    a.set(ylim=(0.8, 1.25), xticks=x, xlim=(-0.6, len(obs) - 0.4))
    a.set_xticklabels([n for _, n in obs], rotation=30, ha='right', rotation_mode='anchor', fontsize=12)
    a.set_ylabel('Direct / predicted', fontsize=13)
    a.legend(loc='upper left', fontsize=12, frameon=False)
    b = d['agree']
    b = b[b.group == 'all']
    rows = [('CALVADOS-2 vs Mpipi-GG', 'three-model set', 'CALVADOS-2 and\nMpipi-GG'),
            ('CALVADOS-2 vs Mpipi-Recharged', 'pair', 'CALVADOS-2 and\nMpipi-Recharged'),
            ('Mpipi-GG vs Mpipi-Recharged', 'three-model set', 'Mpipi-GG and\nMpipi-Recharged')]
    a = pax(fig, 92, 20, 36, 35, 'Agreement between force fields')
    for k, (p, sset, lab) in enumerate(rows):
        y = 2 - k
        r = b[(b.pair == p) & (b.set == sset)].iloc[0]
        a.plot([r.rho_chi, r.rho_total], [y, y], color='0.85', lw=4, zorder=1, solid_capstyle='butt')
        a.scatter(r.rho_chi, y, s=90, color='white', edgecolor='0.45', lw=1.4, zorder=3)
        a.scatter(r.rho_total, y, s=100, color=INK, lw=0, zorder=4)
        a.text(r.rho_total + 0.012, y + 0.18, f'{r.rho_total:.2f}', fontsize=12)
        a.text(r.rho_chi - 0.012, y + 0.18, f'{r.rho_chi:.2f}', fontsize=12, ha='right', color='0.45')
    a.set_yticks(range(3))
    a.set_yticklabels([r[2] for r in rows][::-1], fontsize=11.5)
    a.set(xlim=(0.6, 0.97), ylim=(-0.6, 2.6))
    a.set_xlabel('Rank correlation (Spearman)', fontsize=13)
    a.tick_params(axis='y', length=0)
    a.spines['left'].set_visible(False)
    a.legend(handles=[plt.Line2D([], [], marker='o', ls='', mfc='white', mec='0.45', ms=9, label='fixed solvent'),
                      plt.Line2D([], [], marker='o', ls='', color=INK, ms=9, label='with hydration')],
             loc='lower right', fontsize=11.5, frameon=False)
    takeaway(ax, 'A single equilibrium ensemble suffices, and hydration makes the three force fields agree more')
    return fig


# ----------------------------------------------------------------------------- 06 compaction and charge
def s06(d):
    from scipy import stats  # noqa: F401
    fig, ax = new_slide('Compaction and charge set heat sensitivity',
                        '98.5, 97.1 and 98.2% of IDPs compact on heating in CALVADOS-2, Mpipi-GG and Mpipi-Recharged')
    a = pax(fig, 7, 18, 42, 37, 'Heat sensitivity falls as chains expand')
    byn_panel(a, d['byn'], size=13)
    a = pax(fig, 59, 18, 21, 37, 'At a given compaction')
    g = d['grid']
    lv = ['low', 'mid', 'high']
    M = -g.pivot(index='pi_third', columns='DEK_third', values='median').reindex(index=lv[::-1], columns=lv).to_numpy()
    im = a.imshow(M, cmap='PRGn', vmin=-0.75, vmax=0.75)
    for i in range(3):
        for j in range(3):
            a.text(j, i, FX.minus(f'{M[i, j]:+.2f}'), ha='center', va='center', fontsize=13,
                   color='white' if abs(M[i, j]) > 0.45 else INK)
    a.set_xticks(range(3))
    a.set_xticklabels(lv, fontsize=12)
    a.set_yticks(range(3))
    a.set_yticklabels(lv[::-1], fontsize=12)
    a.tick_params(length=0)
    a.set_xlabel('Asp + Glu + Lys fraction', fontsize=13)
    a.set_ylabel('π fraction (FWYR)', fontsize=13)
    for s in a.spines.values():
        s.set_visible(False)
    cax = fig.add_axes([81.5 / SX, 22 / SY, 0.8 / SX, 29 / SY])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Extra compaction on heating (SD)', fontsize=11.5)
    cb.ax.tick_params(labelsize=11)
    cb.outline.set_visible(False)
    a = pax(fig, 98, 18, 31, 37, 'Residue enrichment')
    D, bio = d['D'], d['bio']
    AA = FX.AA
    f = bio[[f'f_{x}' for x in AA]]
    t = D['C'][['total']].join(f, how='inner').dropna()
    lo, hi = t.total.quantile(0.05), t.total.quantile(0.95)
    base = t[[f'f_{x}' for x in AA]].mean()
    most = np.log2(t[t.total <= lo][[f'f_{x}' for x in AA]].mean() / base).to_numpy()
    least = np.log2(t[t.total >= hi][[f'f_{x}' for x in AA]].mean() / base).to_numpy()
    xx = np.arange(len(AA))
    a.scatter(xx, most, s=55, color=HYD, lw=0, zorder=3, label='5% most heat-sensitive')
    a.scatter(xx, least, s=55, facecolors='white', edgecolors=HYD, lw=1.3, zorder=3, label='5% least heat-sensitive')
    a.axhline(0, color='0.6', lw=0.8)
    for xv in xx[:4]:
        a.axvspan(xv - 0.5, xv + 0.5, color='#F7E3E1', zorder=0, lw=0)
    a.set_xticks(xx)
    a.set_xticklabels(AA, fontsize=11)
    a.tick_params(axis='x', length=0)
    a.set_xlim(-0.6, len(AA) - 0.4)
    a.set_ylabel('log$_2$ enrichment (CALVADOS-2)', fontsize=12.5)
    a.legend(loc='lower right', fontsize=11, frameon=False, handletextpad=0.2)
    a.text(1.5, a.get_ylim()[1], 'charged', ha='center', va='top', fontsize=11, color=HYD)
    takeaway(ax, 'Compact, charge-rich IDPs compact most on heating, and chains with few charges barely respond')
    return fig


# ----------------------------------------------------------------------------- 07 signals and substitutions
def s07(d):
    from scipy import stats
    fig, ax = new_slide('Signals and substitutions retune heat sensitivity in step with compaction',
                        'Exact reweighting of every 300 K ensemble in CALVADOS-2, without simulating any modified chain')
    a = pax(fig, 10, 18, 28, 37, 'Phosphorylating one Ser or Thr')
    phospho_hex(a, d['S'], 13)
    C = d['Cs']
    a = pax(fig, 49, 18, 27, 37, 'Doubling the ionic strength')
    x, y = C.salt_dRg_pct, -C.salt_dtotal
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list('ncpr', plt.get_cmap('viridis_r')(np.linspace(0.1, 1.0, 256)))
    sc = a.scatter(x, y, c=C.NCPR.abs(), cmap=cm, s=22, lw=0, vmin=0, vmax=0.4, zorder=3)
    a.axhline(0, color='0.55', lw=0.6)
    a.axvline(0, color='0.55', lw=0.6)
    rho = stats.spearmanr(C.salt_dRg_pct, C.salt_dtotal)[0]
    a.text(0.97, 0.95, FX.minus(f'ρ = {-rho:.2f}') + f'\n{len(C):,} chains', transform=a.transAxes, ha='right', va='top', fontsize=12)
    a.set(xlim=(-11.5, 23.5), ylim=(-0.5, 0.56))
    a.set_xlabel(f'Change of {FX.RG} (%)', fontsize=13)
    a.set_ylabel('Change of heat\nsensitivity (% per K)', fontsize=13)
    cax = fig.add_axes([77 / SX, 22 / SY, 0.8 / SX, 29 / SY])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label('|net charge per residue|', fontsize=11.5)
    cb.ax.tick_params(labelsize=11)
    cb.outline.set_visible(False)
    a = pax(fig, 100, 18, 29, 37, 'Introducing each residue')
    t = d['mut'].copy()
    t['hs'] = -t.median_change
    t = t.sort_values('hs', ascending=False)
    xx = np.arange(len(t))
    a.bar(xx, t.hs, color=[d['FA'].CLS[r][1] for r in t.index], width=0.75, lw=0)
    a.axhline(0, color=INK, lw=0.9)
    a.set_xticks(xx)
    a.set_xticklabels(t.index, fontsize=12)
    a.tick_params(axis='x', length=0)
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '0' if v == 0 else FX.minus(f'{v:.3f}')))
    a.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.01))
    a.set_ylabel('Median change of heat\nsensitivity (% per K)', fontsize=13)
    seen = {}
    for r in t.index:
        seen.setdefault(d['FA'].CLS[r][0], d['FA'].CLS[r][1])
    a.legend(handles=[Patch(color=c, label=lab) for lab, c in seen.items()], loc='upper right', fontsize=11.5,
             frameon=False, handlelength=1.0)
    note(ax, 112, 13.6, f'{int(t.n.sum()):,} reliable single mutants of 50 chains', 11, col=MUTED)
    takeaway(ax, 'A modification or a change of salt that compacts a chain raises its heat sensitivity')
    return fig


# ----------------------------------------------------------------------------- 08 PAGE4
def s08(d):
    p4 = d['p4']
    fig, ax = new_slide('In PAGE4, one input changes the response to another',
                        'Two kinases switch PAGE4 between functional states. HIPK1 keeps it compact, CLK2 expands it.')
    a = pax(fig, 7, 18, 30, 37, 'Change of size, direct runs')
    forms = ['HIPK1', 'CLK2']
    x = np.arange(2)
    a.bar(x - 0.27, [v for v, _ in p4['saxs']], 0.25, yerr=[e for _, e in p4['saxs']], color=INK, label='SAXS',
          error_kw=dict(lw=1.2))
    for o, m, col in ((0.0, 'CALVADOS-2', C2), (0.27, 'Mpipi-Recharged', MR)):
        a.bar(x + o, [v for v, _ in p4['rg'][m]], 0.25, yerr=[e for _, e in p4['rg'][m]], color=col,
              label=m, error_kw=dict(lw=1.2))
    a.axhline(0, color=INK, lw=0.9)
    a.set_xticks(x)
    a.set_xticklabels(forms, fontsize=13)
    a.tick_params(axis='x', length=0)
    a.set_ylabel(f'Change of {FX.RG} (%)', fontsize=13)
    a.set_ylim(-9, 62)
    a.legend(loc='upper left', fontsize=11.5, frameon=False)
    a = pax(fig, 47, 18, 34, 37, 'Heat sensitivity of each form')
    x = np.arange(3)
    for o, m, col in ((-0.2, 'CALVADOS-2', C2), (0.2, 'Mpipi-Recharged', MR)):
        v = [h[0] for h in p4['hs'][m]]
        e = [h[1] for h in p4['hs'][m]]
        a.bar(x + o, v, 0.38, yerr=e, color=col, label=m, error_kw=dict(lw=1.2))
        for xx, vv, ee in zip(x + o, v, e):
            a.text(xx, vv + ee + 0.01, f'{vv:.2f}', ha='center', va='bottom', fontsize=11.5)
    a.set_xticks(x)
    a.set_xticklabels(['unmodified', 'HIPK1', 'CLK2'], fontsize=13)
    a.tick_params(axis='x', length=0)
    a.set_ylim(0, 0.8)
    a.set_ylabel('Heat sensitivity (% per K)', fontsize=13)
    a.legend(loc='upper right', fontsize=11.5, frameon=False)
    card(ax, 90, 16, 39.5, 42)
    chg = {0: 1, 1: 1, 2: 1, 3: 1, 14: -1, 15: -1, 16: -1, 17: -1, 18: -1}
    draw_chain(ax, best_walk(22, 3.0, 1.1, chg, 1.7, 80), (98.5, 44.5), 1.3, chg, cut=1.7)
    draw_chain(ax, walk(22, 6.0, 0.4, 5), (120.5, 44.5), 0.9, chg, bead=0.6, cut=1.7, phos=(5, 8, 11, 20))
    arrow(ax, (104.5, 44.5), (112.5, 44.5), col=PHOS, lw=2.6)
    ax.text(108.5, 46.8, 'CLK2', ha='center', fontsize=12, color=PHOS, fontweight='bold')
    ax.text(91.5, 55.3, 'Why it happens', fontsize=14, fontweight='bold', color=INK)
    note(ax, 109.7, 27.5, 'CLK2 cuts contacts between the basic\nN-terminal motif and the acidic region\n2.4- to 4-fold, '
         'as NMR relaxation and\nLin et al. also report. Fewer contacts\nremain to tighten on heating.', 12)
    note(ax, 66.7, 13.6, 'Both force fields capture the direction of the CLK2 expansion but underestimate its size.',
         11, col=MUTED)
    takeaway(ax, 'CLK2 phosphorylation, which expands PAGE4, lowers its heat sensitivity, as predicted from one ensemble')
    return fig


# ----------------------------------------------------------------------------- 09 growth temperature
def s09(d):
    rows, diff = d['ogt']
    fig, ax = new_slide('Composition adapts to growth temperature, heat sensitivity does not',
                        'Slopes in SD of organism means with 95% intervals (left), and 44 independent bacterial and '
                        'fungal proteomes (right)')
    a = pax(fig, 20, 18, 42, 37, 'Within bacterial phyla, 444 organisms')
    ogt_bars(a, rows, 13)
    a = pax(fig, 77, 18, 50, 37, 'Thermophile minus mesophile orthologues')
    a.hist(diff, bins=np.linspace(-1, 1, 61), color=C2, lw=0)
    a.axvline(0, color=INK, lw=1.0)
    med = np.median(diff)
    a.text(0.97, 0.95, f'median {FX.minus(f"{med:+.3f}")} % per K\n{len(diff):,} groups', transform=a.transAxes,
           ha='right', va='top', fontsize=12.5)
    a.set_xlim(-1, 1)
    a.set_xlabel('Predicted total response, difference (% per K)', fontsize=13)
    a.set_ylabel('Orthologous groups', fontsize=13)
    takeaway(ax, 'Heat sensitivity is buffered against growth temperature, whereas composition shifts')
    return fig


# ----------------------------------------------------------------------------- 10 take-home
def s10(d):
    fig, ax = new_slide('Take-home and limitations')
    card(ax, 3.5, 13.5, 61, 49)
    card(ax, 68.8, 13.5, 61, 49)
    ax.text(6.5, 58.5, 'Take-home', fontsize=16, fontweight='bold', color=INK)
    ax.text(71.8, 58.5, 'Limitations', fontsize=16, fontweight='bold', color=INK)
    take = ['One equilibrium ensemble gives the response of an IDP\nto temperature, salt, modification and sequence.',
            'With calibrated hydration, nearly every IDP compacts\non heating, set by compaction and charge.',
            'Phosphorylation, acetylation and salt change heat\nsensitivity in step with the compaction they cause.',
            'In PAGE4, CLK2 phosphorylation lowers heat sensitivity,\npredicted and then confirmed directly.',
            'Composition adapts to growth temperature, whereas\nheat sensitivity stays buffered.']
    lim = ['One bead per residue and implicit solvent.',
           'The forcing rests on four calibration constructs and PNt,\nand is tested without fitted parameters on four '
           'proteins.',
           'It misses the Pro-rich Aap and is least constrained for\nweakly charged chains.',
           'Large phosphorylation-driven expansions are\nunderestimated, CLK2-PAGE4 by three- to fourfold.',
           'Phosphorylation is modelled as a change of charge, and\nIDPs were simulated in isolation.']
    for col0, items, dot in ((6.5, take, '#537A9E'), (71.8, lim, MUTED)):
        for k, t in enumerate(items):
            y = 51.5 - 8.0 * k
            ax.add_patch(Circle((col0 + 0.8, y), 0.55, fc=dot, ec='none'))
            ax.text(col0 + 2.6, y, t, fontsize=13, va='center', color=INK, linespacing=1.3)
    takeaway(ax, 'Compaction and charge link disordered ensembles to temperature sensing across the tree of life')
    return fig


# ----------------------------------------------------------------------------- 00 the larger scheme
def thermometer(ax, cx, cy, s=1.0):
    """Drawn with lines and circles, not rounded boxes, so it stays in the picture layer of the pptx export."""
    for lw, col, top, z in ((6.8 * s, INK, 2.2, 6), (4.6 * s, 'white', 2.2, 7), (2.2 * s, HYD, 1.0, 8)):
        ax.plot([cx, cx], [cy - 1.0 * s, cy + top * s], color=col, lw=lw, solid_capstyle='round', zorder=z)
    ax.add_patch(Circle((cx, cy - 1.3 * s), 1.0 * s, fc=HYD, ec=INK, lw=1.3, zorder=9))


def _scale_band(ax, y, fs):
    """Headline numbers of the paper overview, each as stated in the manuscript text (Methods and Results)."""
    card(ax, 3.5, y, 92.5, 10.5, fc=BAND)
    items = [('10,282', 'highly disordered IDPs'), ('67,415', 'simulated trajectories'),
             ('68 million', 'stored conformations'), ('96,672', 'mutants evaluated\nwithout simulation'),
             ('7,528', 'phosphorylation sites\ntested')]
    for k, (num, lab) in enumerate(items):
        xc = 3.5 + 92.5 * (k + 0.5) / len(items)
        ax.text(xc, y + 7.4, num, ha='center', va='center', fontsize=fs(20), fontweight='bold', color='#537A9E')
        ax.text(xc, y + 2.9, lab, ha='center', va='center', fontsize=fs(10.5), color=INK, linespacing=1.15)


def _trajectory(ax, x, top, w, fs, tr, n_snap=4, dT=20.0):
    """Card f of the paper overview: how a susceptibility comes out of one real stored trajectory (the Gly-rich protein
    from tomato, CALVADOS-2, 300 K, as Figs. 2 and 3), without equations. Step 1: Rg of every frame in time order, with
    four conformations; step 2: how much heating favours each frame (minus its total thermal force G, with the hydration
    terms) against its Rg; step 3: the Rg distribution at 300 K and the same frames reweighted by exp(-dT G) for dT of
    heating, with the mean of each and the total response of this chain (Fig. 3e)."""
    xyz, rg, G, seq = np.asarray(tr['xyz'], float), np.asarray(tr['rg'], float), np.asarray(tr['G'], float), str(tr['seq'])
    # step 1: Rg in time order, with conformations
    xs0, xs1 = x + 6.0, x + 38.0
    idx = np.linspace(80, len(rg) - 80, n_snap).astype(int)
    proj = []
    for k in idx:
        p = xyz[k] - xyz[k].mean(0)
        _, _, vt = np.linalg.svd(p, full_matrices=False)
        proj.append(p @ vt[:2].T)
    sc = min(0.44 * (xs1 - xs0) / n_snap / max(np.abs(q[:, 0]).max() for q in proj),
             2.2 / max(np.abs(q[:, 1]).max() for q in proj))
    chg = np.array([{'K': 1, 'R': 1, 'D': -1, 'E': -1}.get(r, 0) for r in seq])
    fx = lambda i: xs0 + (xs1 - xs0) * i / (len(rg) - 1)  # noqa: E731

    def band(v, y0_, y1_):
        return y0_ + (y1_ - y0_) * (v - v.min()) / (v.max() - v.min())

    yr0, yr1 = top - 14.6, top - 11.0
    ax.plot([fx(i) for i in range(len(rg))], band(rg, yr0, yr1), color=C2, lw=0.6, zorder=3)
    ax.text(x + 2.0, (yr0 + yr1) / 2, r'$R_\mathrm{g}$', fontsize=fs(11), va='center', color=C2)
    for j, (k, q) in enumerate(zip(idx, proj)):
        cx, cy = xs0 + (xs1 - xs0) * (j + 0.5) / n_snap, top - 7.0
        ax.plot(cx + sc * q[:, 0], cy + sc * q[:, 1], color='0.6', lw=0.7, zorder=3)
        for v, col in ((1, BASIC), (-1, ACID)):
            m = chg == v
            ax.scatter(cx + sc * q[m, 0], cy + sc * q[m, 1], s=2.5, color=col, lw=0, zorder=4)
        ax.plot([fx(k), cx], [band(rg, yr0, yr1)[k] + 0.3, cy - 2.3], color='0.75', lw=0.6, ls=(0, (1.5, 1.5)), zorder=2)
        ax.scatter([fx(k)], [band(rg, yr0, yr1)[k]], s=8, color=INK, lw=0, zorder=5)
    ax.text((xs0 + xs1) / 2, top - 16.3, f'{len(rg):,} stored frames at 300 K, in time order', fontsize=fs(9.5),
            ha='center', va='center', color=MUTED)
    # step 2: how much heating favours each frame, against its Rg
    arrow(ax, (xs1 + 1.0, top - 10.0), (xs1 + 3.4, top - 10.0))
    bx0, bx1, by0, by1 = xs1 + 6.6, xs1 + 23.0, top - 14.6, top - 5.4
    ax.plot([bx0, bx0, bx1], [by1, by0, by0], color=INK, lw=0.8, zorder=3)
    fav = -G
    px, py = band(rg, bx0 + 0.5, bx1 - 0.3), band(fav, by0 + 0.4, by1 - 0.3)
    ax.scatter(px, py, s=1.2, color=C2, alpha=0.45, lw=0, zorder=3)
    fit = np.polyfit(px, py, 1)
    ax.plot([bx0 + 0.5, bx1 - 0.3], np.polyval(fit, [bx0 + 0.5, bx1 - 0.3]), color=INK, lw=0.9, ls='--', zorder=4)
    ax.text((bx0 + bx1) / 2, by0 - 1.2, r'$R_\mathrm{g}$ of each frame', fontsize=fs(9.5), ha='center', va='center')
    ax.text(bx0 - 0.6, (by0 + by1) / 2, 'Heating favours', fontsize=fs(9.5), ha='right', va='center',
            rotation=90)
    ax.text(bx1 - 0.2, by1 - 0.2, 'Compact frames\ngain weight', fontsize=fs(9.5), ha='right', va='top', color=C2,
            linespacing=1.1)
    # step 3: the same frames, reweighted for heating
    arrow(ax, (bx1 + 1.0, top - 10.0), (bx1 + 3.4, top - 10.0))
    hx0, hx1, hy0, hy1 = bx1 + 5.5, x + w - 2.5, top - 13.6, top - 5.6
    wgt = np.exp(-dT * (G - G.mean()))
    wgt /= wgt.sum()
    edges = np.linspace(np.percentile(rg, 0.5), np.percentile(rg, 99.5), 26)
    h0, _ = np.histogram(rg, edges, density=True)
    h1, _ = np.histogram(rg, edges, weights=wgt, density=True)
    ex = hx0 + (hx1 - hx0) * (edges - edges[0]) / (edges[-1] - edges[0])
    top_h = max(h0.max(), h1.max())
    yy0 = hy0 + (hy1 - hy0) * h0 / top_h
    yy1 = hy0 + (hy1 - hy0) * h1 / top_h
    ax.fill_between(np.repeat(ex, 2)[1:-1], hy0, np.repeat(yy0, 2), color='0.82', lw=0, zorder=2)
    ax.plot(np.repeat(ex, 2)[1:-1], np.repeat(yy1, 2), color=C2, lw=1.2, zorder=3)
    ax.plot([hx0, hx1], [hy0, hy0], color=INK, lw=0.8, zorder=3)
    m0, m1 = rg.mean(), float((wgt * rg).sum())
    mx = lambda v: hx0 + (hx1 - hx0) * (v - edges[0]) / (edges[-1] - edges[0])  # noqa: E731
    for v, col in ((m0, '0.45'), (m1, C2)):
        ax.plot([mx(v), mx(v)], [hy0, hy1 + 0.2], color=col, lw=0.9, ls=(0, (2, 1.5)), zorder=4)
    ax.add_patch(FancyArrowPatch((mx(m0), hy1 + 0.9), (mx(m1) - 0.2, hy1 + 0.9), arrowstyle='-|>', mutation_scale=8,
                                 lw=1.0, color=C2, zorder=5, shrinkA=0, shrinkB=0))
    ax.text(hx1, hy1 - 0.4, '300 K', fontsize=fs(9.5), ha='right', va='center', color='0.45')
    ax.text(hx1, hy1 - 2.0, f'{dT:.0f} K hotter', fontsize=fs(9.5), ha='right', va='center', color=C2)
    ax.text((hx0 + hx1) / 2, hy0 - 1.1, r'$R_\mathrm{g}$', fontsize=fs(9.5), ha='center', va='center')
    ax.text((hx0 + hx1) / 2, top - 16.3, FX.minus(f'Mean shrinks, {tr["total"]:+.2f}% per K'), fontsize=fs(10),
            fontweight='bold', ha='center', va='center', color=C2)


def _force_fields(ax, x, top, fs):
    """Card c of the paper overview: the three force fields, each with a bead icon and what was simulated."""
    for k, (m, txt) in enumerate((('CALVADOS-2', 'all proteins, five temperatures'),
                                  ('Mpipi-GG', 'subset, five temperatures'),
                                  ('Mpipi-Recharged', 'all proteins, 300 K'))):
        yy = top - 8.2 - 6.2 * k
        col = FX.COL[m]
        for j in range(3):
            ax.add_patch(Circle((x + 2.8 + 1.55 * j, yy), 0.72, fc=col, ec='white', lw=0.6, zorder=4))
        ax.plot([x + 2.8, x + 5.9], [yy, yy], color=BOND, lw=1.4, zorder=3)
        ax.text(x + 8.0, yy, m, fontsize=fs(12), fontweight='bold', color=col, va='center')
        ax.text(x + 2.0, yy - 2.4, txt, fontsize=fs(10.5), color=MUTED, va='center')


def _techniques(ax, x, top, fs):
    """Card d of the paper overview: the four measurement techniques as small icons."""
    from matplotlib.patches import Arc
    pos = {'smFRET': (x + 7.5, top - 10.0), 'trFRET': (x + 21.0, top - 10.0), 'SAXS': (x + 7.5, top - 19.5),
           'DLS': (x + 21.0, top - 19.5)}
    for lab, (cx, cy) in pos.items():
        ax.text(cx, cy - 3.6, lab, ha='center', va='center', fontsize=fs(11), color=INK)
    for lab in ('smFRET', 'trFRET'):                  # chain with donor and acceptor dyes
        cx, cy = pos[lab]
        cx -= 0.8 if lab == 'trFRET' else 0.0
        t = np.linspace(-2.4, 2.4, 50)
        ax.plot(cx + t, cy + 0.55 * np.sin(2.6 * t), color=BOND, lw=1.6, zorder=3)
        ax.add_patch(Circle((cx - 2.4, cy + 0.55 * np.sin(-6.24)), 0.8, fc='#2E8C6F', ec='white', lw=0.6, zorder=4))
        ax.add_patch(Circle((cx + 2.4, cy + 0.55 * np.sin(6.24)), 0.8, fc=HYD, ec='white', lw=0.6, zorder=4))
        ax.add_patch(FancyArrowPatch((cx - 1.9, cy + 1.6), (cx + 1.9, cy + 1.6), arrowstyle='-|>', mutation_scale=9,
                                     lw=1.0, color=MUTED, ls=(0, (2, 1.5)), connectionstyle='arc3,rad=-0.3', zorder=4))
        if lab == 'trFRET':                           # a clock for the time-resolved decay
            ax.add_patch(Circle((cx + 4.3, cy + 1.1), 1.0, fc='white', ec=INK, lw=1.0, zorder=4))
            ax.plot([cx + 4.3, cx + 4.3, cx + 4.9], [cy + 1.8, cy + 1.1, cy + 1.1], color=INK, lw=1.0, zorder=5)
    cx, cy = pos['SAXS']                              # X-ray beam, sample, scattering rings
    ax.add_patch(FancyArrowPatch((cx - 4.2, cy), (cx - 1.3, cy), arrowstyle='-|>', mutation_scale=9, lw=1.4,
                                 color=MUTED, zorder=3))
    ax.add_patch(Circle((cx, cy), 0.95, fc=BEAD, ec='white', lw=0.6, zorder=4))
    for r in (1.8, 2.6, 3.4):
        ax.add_patch(Arc((cx, cy), 2 * r, 2 * r, theta1=-45, theta2=45, color=C2, lw=1.3, zorder=3))
    cx, cy = pos['DLS']                               # laser, sample, scattered light
    t = np.linspace(-4.3, -1.3, 60)
    ax.plot(cx + t, cy + 0.35 * np.sin(9 * t), color=HYD, lw=1.3, zorder=3)
    ax.add_patch(Circle((cx, cy), 0.95, fc=BEAD, ec='white', lw=0.6, zorder=4))
    for ang in (-50, -20, 10, 40):
        a = np.radians(ang)
        ax.plot([cx + 1.4 * np.cos(a), cx + 3.2 * np.cos(a)], [cy + 1.4 * np.sin(a), cy + 3.2 * np.sin(a)],
                color=MUTED, lw=1.1, zorder=3)


def overview(ax, d, f=1.0, floor=0.0, paper=False, rows=False):
    """The larger scheme in three cards and one strip, on slide units. f scales the fonts and floor is the smallest font
    in these units (the paper figure, fig_overview.py, is drawn on the same units and scaled down in LaTeX). rows (the
    paper figure) adds a second row below a and b, with the force fields (c), the measurement techniques (d) and what it
    reveals (e), and the strip under it."""
    fs = lambda s: max(s * f, floor)  # noqa: E731
    y0, ch = 20.0, (41.0 if rows else 42.0)
    yc, hc = (-6.5, 24.5) if rows else (y0, ch)                # second row (rows) or card c
    tb = 'One ensemble answers many perturbations' if paper else 'One ensemble answers many inputs'
    cards = [(3.5, 29.0, 'Sequences across\nthe tree of life', y0, ch), (37.5, 58.5, tb, y0, ch)]
    cards += ([(3.5, 26.5, 'Three force fields', yc, hc), (33.0, 28.5, 'Measurements of 10 proteins', yc, hc),
               (64.5, 31.5, 'What it reveals', yc, hc),
               (3.5, 92.5, 'From one trajectory to a susceptibility', yc - 19.5, 17.5)] if rows else
              [(101.0, 28.8, 'What it reveals', y0, ch)])
    for k, (x, w, t, yb, hh) in enumerate(cards):
        card(ax, x, yb, w, hh)
        badge(ax, x + 2.6, yb + hh - 3.2, 'abcdef'[k] if paper else k + 1)
        ax.text(x + 4.9, yb + hh - 3.2, t, fontsize=fs(12.5), fontweight='bold', color=INK, va='center',
                linespacing=1.15)
    arrow(ax, (33.4, y0 + ch / 2), (36.2, y0 + ch / 2))
    if rows:
        arrow(ax, (80.25, y0 - 0.4), (80.25, yc + hc + 0.4))
    else:
        arrow(ax, (97.0, y0 + ch / 2), (99.8, y0 + ch / 2))

    # 1 cladogram of the six groups of the proteome set (counts from the bio table), viruses set apart
    n = d['bio'].kingdom.value_counts()
    x = 3.5 - (1.0 if paper else 0.0)
    leaf = {'Bacteria': 51.5, 'Protists': 47.3, 'Plants': 43.1, 'Fungi': 38.9, 'Mammals': 34.7, 'Viruses': 28.6}
    lcol = {'Bacteria': C2, 'Protists': '#2E8C6F', 'Plants': '#2E8C6F', 'Fungi': MR, 'Mammals': MR, 'Viruses': MUTED}
    xl = x + 12.5

    def elbow(xn, yn, kids):
        for xk, yk in kids:
            ax.plot([xn, xn, xk], [yn, yk, yk], color=BOND, lw=1.8, zorder=3, solid_capstyle='round')

    fm = (x + 10.0, (leaf['Fungi'] + leaf['Mammals']) / 2)
    pfm = (x + 8.2, (leaf['Plants'] + fm[1]) / 2)
    euk = (x + 6.4, (leaf['Protists'] + pfm[1]) / 2)
    root = (x + 4.2, (leaf['Bacteria'] + euk[1]) / 2)
    elbow(*fm, [(xl, leaf['Fungi']), (xl, leaf['Mammals'])])
    elbow(*pfm, [(xl, leaf['Plants']), fm])
    elbow(*euk, [(xl, leaf['Protists']), pfm])
    elbow(*root, [(xl, leaf['Bacteria']), euk])
    ax.plot([x + 2.6, root[0]], [root[1]] * 2, color=BOND, lw=1.8, zorder=3)
    ax.plot([x + 8.2, xl], [leaf['Viruses']] * 2, color=BOND, lw=1.8, ls=(0, (1.5, 1.5)), zorder=3)
    for g, yy in leaf.items():
        ax.add_patch(Circle((xl, yy), 0.45 + 1.05 * np.sqrt(n[g] / n.max()), fc=lcol[g], ec='white', lw=0.8, zorder=5))
        ax.text(xl + 2.3, yy, g, fontsize=fs(11), color=INK, va='center')
        ax.text(3.5 + 27.6, yy, f'{n[g]:,}', fontsize=fs(11), color=MUTED, va='center', ha='right')
    note(ax, 3.5 + 14.5, y0 + 4.0, f'{d["bio"].taxonomy_id.nunique():,} organisms\nin six groups' if rows else
         f'{len(d["bio"]):,} highly disordered IDPs\nfrom {d["bio"].taxonomy_id.nunique():,} organisms', fs(11))

    # 2 perturbations feeding one ensemble
    x = 37.5
    hub = (x + (45.0 if paper else 42.0), 41.0)
    ax.add_patch(Circle(hub, 9.6, fc='white', ec=PALE, lw=2.5, zorder=2))
    chg = {0: 1, 1: 1, 2: 1, 3: 1, 14: -1, 15: -1, 16: -1, 17: -1, 18: -1}
    draw_chain(ax, best_walk(22, 3.0, 1.1, chg, 1.7, 80), (hub[0] - 4.1, hub[1] + 0.2), 1.3, chg, cut=1.7)
    draw_chain(ax, walk(22, 6.0, 0.4, 5), (hub[0] + 4.6, hub[1] + 0.2), 0.8, chg, bead=0.6, cut=1.7)
    note(ax, hub[0], hub[1] + 12.4, 'One ensemble at 300 K' if paper else 'One simulated ensemble at 300 K', fs(11.5))
    note(ax, hub[0], hub[1] - 12.6, 'Three force fields\nwith calibrated hydration', fs(11), col=MUTED)
    ys = [52.0, 45.0, 38.0, 31.0]
    labs = ['Heat', 'Salt', 'Phosphorylation', 'Substitution']
    cols = [HYD, BASIC, PHOS, '#2E8C6F']
    xi = x + 5.0
    how = ['hydration and electrostatics', 'electrostatic screening', 'charge of the site', 'residue parameters']
    for yy, lab, col, hw in zip(ys, labs, cols, how):
        ax.add_patch(Circle((xi, yy), 2.55, fc='white', ec=col, lw=2.0, zorder=4))
        ax.text(xi + 3.7, yy + (0.9 if rows else 0.0), lab, fontsize=fs(12), color=INK, va='center', fontweight='bold')
        if rows:                                   # what each perturbation changes in the energy function
            ax.text(xi + 3.7, yy - 1.3, hw, fontsize=fs(10), color=MUTED, va='center')
        ang = np.arctan2(yy - hub[1], -9.6 * 2)
        tip = (hub[0] + 10.3 * np.cos(ang), hub[1] + 10.3 * np.sin(ang))
        arrow(ax, (x + (27.6 if paper else 23.2), yy), tip, col=col, lw=2.0, ms=14)
    thermometer(ax, xi, ys[0] + 0.25, 0.62)
    for dx, c, s in ((-0.95, BASIC, '+'), (0.95, ACID, '−')):
        ax.add_patch(Circle((xi + dx, ys[1]), 0.88, fc=c, ec='white', lw=0.6, zorder=6))
        ax.text(xi + dx, ys[1] - 0.05, s, ha='center', va='center', color='white', fontsize=9 * f, fontweight='bold',
                zorder=7)
    ax.add_patch(Circle((xi, ys[2]), 1.3, fc=PHOS, ec='white', lw=0.6, zorder=6))
    ax.text(xi, ys[2] - 0.05, 'P', ha='center', va='center', color='white', fontsize=10 * f, fontweight='bold', zorder=7)
    for dx, c in ((-1.0, BASIC), (1.0, ACID)):
        ax.add_patch(Circle((xi + dx, ys[3]), 0.8, fc=c, ec='white', lw=0.6, zorder=6))
    ax.text(xi, ys[3] + 1.5, '→', ha='center', va='center', fontsize=10 * f, color=INK, zorder=7)
    note(ax, x + 29.25, y0 + 4.0, 'Each perturbation is computed from the same stored conformations' if paper else
         'Each input reweights the same stored frames, so one simulation\nper IDP gives its response to every input',
         fs(11))

    # 3 what it reveals
    x = 64.5 if rows else 101.0
    rows_ = [('Rules', 'Nearly every IDP compacts\non heating. Compaction\nand charge set how much.', C2),
            ('Signals', 'Phosphorylation and salt\nchange heat sensitivity\nin step with compaction.' if paper else
             'Phosphorylation and salt\nretune heat sensitivity\nwith compaction.', PHOS),
            ('Evolution', 'Composition adapts to\ngrowth temperature, heat\nsensitivity is buffered.', '#2E8C6F')]
    if rows:                                       # two-line texts in the narrower second-row card
        rows_ = [('Rules', 'Nearly every IDP compacts on heating,\nset by compaction and charge.', C2),
                 ('Signals', 'Phosphorylation and salt change heat\nsensitivity in step with compaction.', PHOS),
                 ('Evolution', 'Composition adapts to growth\ntemperature, heat sensitivity does not.', '#2E8C6F')]
        xr = [x + 1.7] * 3
        yr = [yc + hc - 8.0 - 6.3 * k for k in range(3)]
        tdy, tx = 1.9, 3.6                         # text to the right of the dot
    else:
        xr = [x + 2.0] * 3
        yr = [51.5 - 11.0 * k for k in range(3)] if paper else [50.8 - 10.6 * k for k in range(3)]
        tdy, tx = (3.5 if paper else 2.6), 0.0
    for xx, yy, (h, t, col) in zip(xr, yr, rows_):
        ax.add_patch(Circle((xx + 1.2, yy + 1.6), 1.45, fc=col, ec='none', zorder=4))
        ax.text(xx + 3.6, yy + 1.6, h, fontsize=fs(12), fontweight='bold', color=INK, va='center')
        ax.text(xx + tx, yy - tdy, t, fontsize=fs(10.5), color=INK, va='center', linespacing=1.25)
    xb = [xx + 1.2 for xx in xr]                   # small white icons inside the row dots
    for k, hgt in enumerate((1.2, 0.8, 0.45)):
        ax.add_patch(Rectangle((xb[0] - 0.85 + 0.6 * k, yr[0] + 0.8), 0.45, hgt, fc='white', ec='none', zorder=5))
    ax.text(xb[1], yr[1] + 1.55, 'P', ha='center', va='center', color='white', fontsize=9.5 * f, fontweight='bold',
            zorder=5)
    ax.plot([xb[2] - 0.8, xb[2], xb[2] + 0.8], [yr[2] + 0.8, yr[2] + 2.0, yr[2] + 0.8], color='white', lw=1.3, zorder=5)
    ax.plot([xb[2], xb[2]], [yr[2] + 2.0, yr[2] + 2.8], color='white', lw=1.3, zorder=5)

    if rows:
        _scale_band(ax, 63.0, fs)
        _force_fields(ax, 3.5, yc + hc, fs)
        _techniques(ax, 33.0, yc + hc, fs)
        _trajectory(ax, 3.5, yc - 2.0, 92.5, fs, d['traj'])

    # the larger scheme, one strip
    ys = -29.5 if rows else 15.4
    ax.text(4.0, ys, 'The larger scheme', fontsize=fs(12), color=MUTED, va='center', fontstyle='italic')
    xs = (26.5, 43.5, 61.0, 84.0) if rows else (31.0, 54.0, 79.5, 111.5) if paper else (27.0, 49.5, 75.5, 108.5)
    steps = list(zip(('Sequence', 'Ensemble', 'Susceptibility', 'Function and adaptation'), xs))
    pill = lambda s: 0.78 * f * len(s) + 5.0  # noqa: E731
    for k, (s, xc) in enumerate(steps):
        hl = s == 'Susceptibility'
        card(ax, xc - pill(s) / 2, ys - 2.2, pill(s), 4.4, fc=BAND if hl else CARD, r=2.2)
        ax.text(xc, ys, s, ha='center', va='center', fontsize=fs(12.5), fontweight='bold',
                color='#537A9E' if hl else INK)
        if k:
            p, xp = steps[k - 1]
            arrow(ax, (xp + pill(p) / 2 + 1.0, ys), (xc - pill(s) / 2 - 1.0, ys), ms=13, lw=1.8)
    ax.text(xs[2], ys - 3.5, 'added by this work', ha='center', va='center', fontsize=fs(10), color='#537A9E',
            fontstyle='italic')


def s00(d):
    fig, ax = new_slide('From sequence to ensemble to susceptibility',
                        'IDP ensembles shift with temperature, salt and modification. We ask what sets that '
                        'susceptibility, protein by protein, across the tree of life.')
    overview(ax, d)
    takeaway(ax, 'One ensemble per IDP shows how disordered proteins across the tree of life respond to their environment')
    return fig


SLIDES = [('00_big_picture', s00), ('01_forcing_steps', s01), ('02_forcing_one_chain', s02), ('03_forcing_calibration', s03),
          ('04_findings_overview', s04), ('05_one_ensemble', s05), ('06_compaction_charge', s06),
          ('07_signals_substitutions', s07), ('08_page4', s08), ('09_growth_temperature', s09),
          ('10_takehome', s10)]


def main():
    FX.house()
    plt.rcParams.update({'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
                         'xtick.major.size': 4, 'ytick.major.size': 4, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
                         'axes.labelsize': 13, 'legend.fontsize': 12, 'axes.titlesize': 13})
    d = load_all()
    import slides_pptx
    os.makedirs(OUT, exist_ok=True)
    deck = slides_pptx.Deck(W, H)
    for name, make in SLIDES:
        fig = make(d)
        fig.savefig(os.path.join(OUT, f'slide_{name}.png'), dpi=150, facecolor='white')     # preview of the whole slide
        deck.add(fig, fig.axes[0], name)
        plt.close(fig)
        print('written', name, flush=True)
    print('written', deck.save(os.path.join(OUT, 'BENDER-MT_slides.pptx')))

if __name__ == '__main__':
    main()
