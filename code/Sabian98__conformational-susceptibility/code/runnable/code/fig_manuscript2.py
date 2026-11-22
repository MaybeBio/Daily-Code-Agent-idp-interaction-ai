"""
fig_manuscript2.py — further manuscript figures (2026-10-10; house style 2026-10-15), from stored outputs only.
  fig1       main Fig. 1: a approach; b chains per taxonomic group in each model; c chain lengths; d apparent nu in each
             model; e, f two worked examples; g-i validation of the fluctuation relation (CALVADOS-2 scatter, share of chains
             within two SE for eight observables in CALVADOS-2 and Mpipi-GG, median chi(Rg) by own-nu bin)
  ed_parts   Fig. S1: median sticker, charge (Boltzmann), electrostatic dH/dT and bond parts of chi_T(Rg) at 300 K
             by own-nu bin, three models (Analysis_FDT_RECHARGED/alongside/A_chi_by_nu.csv)
  pairs      Fig. S3: chi per pair by pair type, three models; CALVADOS-2 and Mpipi-GG from
             Analysis_FDT_MPIPI/compact/pair_type_split/per_protein.csv, Mpipi-Recharged from traj_pass/pairs.csv for the
             same sample (as recharged_traj_analysis.py)
  dynamics   Fig. S4: rank partial correlation (nu and log N fixed) of breathing amplitude, relaxation time,
             compact-expanded switching and bimodality with chi and with the total, from total_response/dynamics/
             per_chain.csv and the current responses (rg_dynamics.partial; Mpipi-GG on its 1,107 main-set chains)
  hotspots   Fig. S5: median contribution of each residue type to the total heating response, CALVADOS-2 and
             Mpipi-Recharged (heat_crosses / traj_pass analysis_residue_types.csv, random chains)
  mutations  Fig. S7: median change of the total by residue introduced and by residue replaced, CALVADOS-2 vs
             Mpipi-Recharged (total_response/mutation_scan/compare_by_{mut,wt}.csv)
  ogt        Fig. S16: growth temperature. a BENDER bacteria, organism means of the total response against OGT after
             removing phylum means, phyla with >= 3 organisms as ogt_response.py (CALVADOS-2, Mpipi-Recharged); b
             predicted total response (CALVADOS-2 predictor), thermophile minus mesophile orthologue groups by clade
             (Thermal adaptation/predictions/A_matched.csv); c composition slopes per 10 C in SD of organism means, BENDER
             bacteria with phylum fixed effects vs the independent proteomes with clade fixed effects
             (bender_ogt/followup_composition.csv, checks_regions.csv)
House style (2026-10-15): 6.5 in wide (printed at 100%), Arial via figfont, text >= 6.5 pt, axis labels 7.5 pt, bold
9-pt panel letters top left outside the axes. Fig. 1 needs nothing from fig_extra / fig_manuscript at run time (its
panels are drawn here).
    BENDER_REANALYSIS_DIR=... python fig_manuscript2.py [names]      (openmpipi env)
Output: BENDER-MT/manuscript/figs/{fig1_dataset_validation,fig_ed_parts,fig_ogt,fig_ed_pairs,fig_ed_hotspots,fig_ed_dynamics,
        fig_ed_mutations}.{png,pdf}
"""
import sys

sys.dont_write_bytecode = True

import functools  # noqa: E402
import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import figfont  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
ROOT = (__import__('os').environ['CS_DATA'] + '/BENDER-MT')
TLA = f'{ROOT}/Traj Level Analysis'
TR = f'{TLA}/Analysis_FDT/total_response'
TA = f'{ROOT}/Thermal adaptation'
OUT = f'{ROOT}/manuscript/figs'
COL = {'CALVADOS-2': '#2E6FA3', 'Mpipi-GG': '#A9502A', 'Mpipi-Recharged': '#6B4FA0'}
MODELS = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
KINGDOMS = ['Bacteria', 'Plants', 'Fungi', 'Mammals', 'Protists', 'Viruses']
BIN_MID = {'(0.3, 0.45]': 0.43, '(0.45, 0.5]': 0.475, '(0.5, 0.54]': 0.52, '(0.54, 0.58]': 0.56, '(0.58, 0.63]': 0.605,
           '(0.63, 1.0]': 0.65}
BIN_LAB = ['0.30–0.45', '0.45–0.50', '0.50–0.54', '0.54–0.58', '0.58–0.63', '>0.63']   # house compaction-bin labels
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.titlesize': 8, 'legend.fontsize': 7})
figfont.apply()                                   # Arial for all figures (fig_supplementary imports this module)

# ----------------------------------------------------------------------------- house style
W = 6.5                                           # LaTeX text width (in); figures print at 100%
INK = '#222222'
MINUS = '−'
RG, REE = r'$R_\mathrm{g}$', r'$R_\mathrm{ee}$'
NU = r'$\mathtt{\nu}$'   # Greek nu from STIX italic (figfont mathtext.tt); Arial's nu reads as v
CHI = r'$\chi_T(R_\mathrm{g})$'
STYLE = {'font.size': 7.5, 'axes.labelsize': 7.5, 'axes.titlesize': 7.5, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
         'legend.fontsize': 6.5, 'legend.title_fontsize': 6.5, 'axes.spines.top': False, 'axes.spines.right': False,
         'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.5,
         'ytick.minor.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5,
         'ytick.minor.size': 1.5, 'xtick.major.pad': 2.0, 'ytick.major.pad': 2.0, 'axes.labelpad': 2.5,
         'axes.titlepad': 3.5, 'legend.frameon': False, 'legend.handlelength': 1.3, 'legend.handletextpad': 0.5,
         'legend.borderaxespad': 0.3, 'legend.borderpad': 0.2, 'legend.labelspacing': 0.3, 'legend.columnspacing': 1.2,
         'lines.linewidth': 1.0, 'patch.linewidth': 0.5, 'axes.unicode_minus': True, 'pdf.fonttype': 42,
         'ps.fonttype': 42}


def house(fn):
    """Run a figure function with the house rcParams (restored afterwards, so importers keep their own)."""
    @functools.wraps(fn)
    def run(*args, **kwargs):
        with plt.rc_context(STYLE):
            figfont.apply()                       # after the rcParams, as figfont requires
            return fn(*args, **kwargs)
    return run


def num(x, fmt):
    """Formatted number with a true minus sign."""
    return format(x, fmt).replace('-', MINUS)


def ax_in(fig, left, bottom, width, height, **kw):
    """Axes placed in inches from the bottom left of the figure."""
    fw, fh = fig.get_size_inches()
    return fig.add_axes([left / fw, bottom / fh, width / fw, height / fh], **kw)


def label(ax, s, x=-0.22, y=1.06):
    """Old-style panel letter (kept for fig_supplementary.py)."""
    ax.text(x, y, s, transform=ax.transAxes, fontweight='bold', fontsize=10, va='bottom')


def letters(fig, items, pad=2.0):
    """Bold 9-pt panel letters, top left outside each axes: left edge at the left of everything drawn for the axes
    (tick and axis labels), bottom just above the axes top, raised above any title, legend or text that would touch it.
    Letters of axes with the same top share one height. items: [(ax, 'a'), ...]."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    px = fig.dpi / 72
    placed = []
    for ax, s in items:
        tb = ax.get_tightbbox(r)
        t = fig.text(0, 0, s, fontsize=9, fontweight='bold', ha='left', va='bottom')
        bb = t.get_window_extent(r)
        w, h = bb.width, bb.height
        x0, y0 = tb.x0, ax.bbox.y1 + pad * px
        blockers = [ax.title, ax._left_title, ax._right_title, *ax.texts]
        if ax.get_legend() is not None:
            blockers.append(ax.get_legend())
        boxes = [a.get_window_extent(r) for a in blockers
                 if a.get_visible() and (not hasattr(a, 'get_text') or a.get_text())]
        yb = ax.yaxis.get_tightbbox(r)            # drawn y tick labels and the y label
        if yb is not None:
            boxes.append(yb)
        for _ in range(5):
            hit = [b for b in boxes if b.x1 > x0 - pad * px and b.x0 < x0 + w + pad * px and b.y1 > y0 - pad * px
                   and b.y0 < y0 + h + pad * px]
            if not hit:
                break
            y0 = max(b.y1 for b in hit) + pad * px
        placed.append([t, ax, x0, y0])
    for p in placed:                              # same top, same letter height
        same = [q[3] for q in placed if abs(q[1].bbox.y1 - p[1].bbox.y1) < 2 * px]
        p[3] = max(same)
    inv = fig.transFigure.inverted()               # store in figure fraction (the saved dpi differs)
    for t, ax, x0, y0 in placed:
        t.set_position(inv.transform((x0, y0)))
    return [p[0] for p in placed]


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=300)
    plt.close(fig)
    print('saved', name)


def load_models():
    import recharged_compare as RC
    D = RC.load()
    bio = pd.read_csv(RC.BIO_CSV, usecols=['UniProt_ID', 'kingdom', 'length'], low_memory=False).set_index('UniProt_ID')
    return {k: d.join(bio, how='left', rsuffix='_bio') for k, d in D.items()}


# ----------------------------------------------------------------------------- Fig. 1a schematic (from fig_extra)
FDT = f'{TLA}/Analysis_FDT'
NPZ = f'{FDT}/npz'
TEMPS = (278, 288, 300, 310, 320)
EX = {'Q02021': {'name': 'Gly-rich protein', 'org': 'tomato', 'col': '#C8792A', 'dir': f'{ROOT}/BENDER_300K/Plants/Q02021'},
      'P42143': {'name': 'Protamine P1', 'org': 'marsupial mole', 'col': '#2E8C6F', 'dir': f'{ROOT}/BENDER_300K/Mammals/P42143'}}


def _block_se(a, n=10):
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
    sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.expanduser((__import__('os').environ['CS_REPO'] + '/bender_reanalysis'))))
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
    return z['rg'], z['V_AH'] + z['V_DH'] + z['V_bond']


def _box(ax, x, y, w, h, fc, ec='none', r=1.2, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw,
                                transform=ax.transData, zorder=1))


def _arrow(ax, x0, y0, x1, y1, col='0.35'):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8, lw=0.9, color=col, zorder=3,
                                 shrinkA=0, shrinkB=0))


def schematic(ax):
    """Proteome, three force fields, one ensemble per temperature, fluctuation relation, responses, tests. Data units
    are square (x 0-100, y 0-18.3 on a 6.28 x 1.15 in axes)."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 18.3)
    ax.axis('off')
    tint, top, mid = '#F1F3F5', 17.5, 9.2
    hd = {'ha': 'center', 'va': 'top', 'fontsize': 7.5, 'fontweight': 'bold', 'linespacing': 1.1}
    tx = {'ha': 'center', 'fontsize': 6.5, 'linespacing': 1.2}
    # 1 proteomes
    _box(ax, 0, 0.2, 15, 17.9, tint)
    ax.text(7.5, top, 'Disordered\nproteomes', **hd)
    ax.text(7.5, 6.9, '10,282 proteins\n2,156 organisms\nbacteria, plants,\nfungi, mammals,\nprotists, viruses', va='center',
            **tx)
    _arrow(ax, 15.4, mid, 17.6, mid)
    # 2 ensembles, with three conformations of the Gly-rich protein
    _box(ax, 18, 0.2, 22, 17.9, tint)
    ax.text(29, top, 'Coarse-grained\nensembles', **hd)
    for c, cx in zip(_conformations('Q02021'), (22.6, 29.0, 35.4)):
        xy = c * (2.3 / max(1e-6, np.abs(c).max())) + np.array([cx, 10.6])
        ax.plot(xy[:, 0], xy[:, 1], '-', color='0.55', lw=0.4, zorder=2)
        ax.scatter(xy[:, 0], xy[:, 1], s=1.6, c=np.linspace(0, 1, len(xy)), cmap='viridis', zorder=3, lw=0)
    for k, (m, c) in enumerate(COL.items()):
        ax.text(29, 6.9 - 1.75 * k, m, va='center', color=c, fontweight='bold', **{**tx, 'linespacing': 1.0})
    ax.text(29, 6.9 - 1.75 * 3, '278–320 K', va='center', **tx)
    _arrow(ax, 40.4, mid, 42.6, mid)
    # 3 fluctuation relation, with the real Rg-energy covariance of the Gly-rich protein at 300 K
    _box(ax, 43, 0.2, 22, 17.9, tint)
    ax.text(54, top, 'Fluctuation–response', **hd)
    ax.text(54, 15.0, r'$\chi_X = -\mathrm{Cov}(A,\,\partial(\beta H)/\partial X)$', ha='center', va='top',
            fontsize=7)
    rg, H = _covariance('Q02021')
    ins = ax.inset_axes([47.5, 2.6, 15.5, 8.6], transform=ax.transData)
    ins.scatter(H, rg, s=0.5, color='#2E6FA3', alpha=0.35, lw=0, rasterized=True)
    p = np.polyfit(H, rg, 1)
    hh = np.linspace(np.percentile(H, 1), np.percentile(H, 99), 2)
    ins.plot(hh, np.polyval(p, hh), color=INK, lw=0.8)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_xlabel('energy $H$', fontsize=6.5, labelpad=1)
    ins.set_ylabel(RG, fontsize=6.5, labelpad=1)
    for s in ('top', 'right'):
        ins.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ins.spines[s].set_linewidth(0.6)
    ins.patch.set_alpha(0)
    _arrow(ax, 65.4, mid, 67.6, mid)
    # 4 responses
    _box(ax, 68, 0.2, 15, 17.9, tint)
    ax.text(75.5, top, 'Responses', **hd)
    for k, t in enumerate(('temperature', '+ hydration', 'phosphorylation', 'salt', 'substitutions')):
        y = 13.9 - 2.95 * k
        _box(ax, 69.0, y - 1.2, 13.0, 2.4, 'white', ec='0.7', r=1.0, lw=0.5)
        ax.text(75.5, y, t, va='center', **tx)
    _arrow(ax, 83.4, mid, 85.6, mid)
    # 5 tests
    _box(ax, 86, 0.2, 14, 17.9, tint)
    ax.text(93, top, 'Tests', **hd)
    ax.text(93, 13.4, 'direct runs\n278–320 K', va='center', **tx)
    ax.text(93, 8.0, '10 proteins by\nsmFRET, trFRET,\nSAXS, DLS', va='center', **tx)
    ax.text(93, 2.6, 'independent\ntemperature series', va='center', **tx)


# ----------------------------------------------------------------------------- Fig. 1 panels
def example_dist(ax):
    """Rg distributions at 278 K (filled) and 320 K (line) for the two example chains, CALVADOS-2."""
    for uid, e in EX.items():
        lo, hi = _rg(uid, 278), _rg(uid, 320)
        bins = np.linspace(min(lo.min(), hi.min()), max(lo.max(), hi.max()), 36)
        ax.hist(lo, bins=bins, density=True, color=e['col'], alpha=0.35, lw=0)
        ax.hist(hi, bins=bins, density=True, histtype='step', color=e['col'], lw=0.9)
    ymax = ax.get_ylim()[1]
    for (uid, e), (x, ha) in zip(EX.items(), ((0.03, 'left'), (1.0, 'right'))):
        ax.text(x, 1.0, f"{e['name']}\n({e['org']})", transform=ax.transAxes, ha=ha, va='top', fontsize=6.5,
                color=e['col'], linespacing=1.1)
    ax.set_ylim(0, ymax * 1.45)
    ax.set_xlabel(f'{RG} (nm)')
    ax.set_ylabel('Probability density')
    h = [Patch(color='0.5', alpha=0.35, lw=0, label='278 K'), Line2D([], [], color='0.4', lw=0.9, label='320 K')]
    ax.legend(handles=h, loc='center right', handlelength=1.1, bbox_to_anchor=(1.03, 0.5))


def example_curve(ax):
    """Mean Rg(T) relative to 300 K (points, 10-block SE) and the 300 K fluctuation-relation tangent (dashed)."""
    for uid, e in EX.items():
        r300 = _rg(uid, 300).mean()
        m = np.array([_rg(uid, T).mean() for T in TEMPS])
        se = np.array([_block_se(_rg(uid, T)) for T in TEMPS])
        ax.errorbar(TEMPS, 100 * (m / r300 - 1), yerr=100 * se / r300, fmt='o', ms=2.8, color=e['col'], elinewidth=0.7,
                    capsize=0)
        chi, _ = _chi300(uid)
        t = np.array([276, 322])
        ax.plot(t, 100 * chi * (t - 300) / r300, '--', color=e['col'], lw=0.9)
    ax.axhline(0, color='0.75', lw=0.5, zorder=0)
    ax.set_xlim(275, 323)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(f'Change of mean {RG}\nfrom 300 K (%)')
    h = [Line2D([], [], color='0.4', marker='o', ms=2.8, ls='none', label='direct simulations'),
         Line2D([], [], color='0.4', ls='--', lw=0.9, label='fluctuation relation')]
    ax.legend(handles=h, loc='lower right', handlelength=1.6, borderaxespad=0.1)


def panel_scatter(ax, lim=(-0.2, 1.2)):
    """Integrated susceptibility (278-320 K) against the directly computed change of mean Rg, CALVADOS-2, all chains;
    chains outside the drawn range are counted on the panel."""
    v = pd.read_csv(f'{FDT}/fdt_validation.csv', usecols=['UniProt_ID', 'metric', 'integral', 'observed_diff'])
    r = v[v.metric == 'rg']
    out = ~(r.integral.between(*lim) & r.observed_diff.between(*lim))
    i = r[~out]
    ax.scatter(i.integral, i.observed_diff, s=1.0, alpha=0.25, color=COL['CALVADOS-2'], rasterized=True, lw=0)
    ax.plot(lim, lim, color='k', lw=0.6)
    ax.set(xlim=lim, ylim=lim, xlabel=f'Predicted change of mean {RG} (nm)', ylabel=f'Direct change of\nmean {RG} (nm)')
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.text(0.03, 0.98, f'CALVADOS-2\n278 to 320 K\n{r.UniProt_ID.nunique():,} chains', transform=ax.transAxes, va='top',
            fontsize=6.5, linespacing=1.15)
    ax.text(1.0, 1.02, f'{int(out.sum())} chains outside the axes', transform=ax.transAxes, ha='right', va='bottom',
            fontsize=6.5, color='0.35')
    return int(out.sum()), len(r)


OBS = [('rg', RG), ('ree', REE), ('delta', 'Asphericity'), ('fragmentation_index', 'Fragmentation'),
       ('global_efficiency', 'Global efficiency'), ('avg_clustering', 'Clustering'), ('transitivity', 'Transitivity'),
       ('degree_assortativity', 'Assortativity')]


def panel_ratio(ax):
    """Directly computed / predicted change (278-320 K) for the eight observables, CALVADOS-2 and Mpipi-GG."""
    sc = pd.read_csv(f'{FDT}/fdt_validation_summary.csv')
    sc = sc[sc['T'] == 'integral'].set_index('metric')
    sm = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/mpipi_fdt_validation_summary.csv')
    sm = sm[sm['T'] == 'integral'].set_index('metric')
    x = np.arange(len(OBS))
    rc = [sc.loc[m, 'mean_slope'] / sc.loc[m, 'mean_chi'] for m, _ in OBS]
    rm = [sm.loc[m, 'ratio_observed_to_integral'] for m, _ in OBS]
    nc, nm = int(sc.loc['rg', 'n']), int(sm.loc['rg', 'n'])
    ax.bar(x - 0.2, rc, 0.38, color=COL['CALVADOS-2'], label=f'CALVADOS-2 ({nc:,} chains)')
    ax.bar(x + 0.2, rm, 0.38, color=COL['Mpipi-GG'], label=f'Mpipi-GG ({nm:,} chains)')
    ax.axhline(1, color='k', lw=0.6)
    top = max(1.25, np.ceil((max(rc + rm) + 0.02) * 20) / 20)
    ax.set(ylim=(0.8, top), ylabel='Direct / predicted\nchange', xticks=x, xlim=(-0.6, len(OBS) - 0.4))
    ax.set_yticks(np.arange(0.8, top + 1e-9, 0.1))
    ax.set_xticklabels([n for _, n in OBS], rotation=30, ha='right', rotation_mode='anchor')
    ax.legend(loc='upper left', handlelength=1.1)
    return rc, rm, nc, nm


def panel_agree(ax):
    """Share of chains whose directly computed change (278-320 K) lies within two standard errors of the integrated
    prediction, eight observables, CALVADOS-2 and Mpipi-GG (frac_abs_z_lt_2 of the validation summaries). The dashed line
    is the 95.4 % expected when the misses are sampling noise alone."""
    sc = pd.read_csv(f'{FDT}/fdt_validation_summary.csv')
    sc = sc[sc['T'] == 'integral'].set_index('metric')
    sm = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/mpipi_fdt_validation_summary.csv')
    sm = sm[sm['T'] == 'integral'].set_index('metric')
    x = np.arange(len(OBS))
    fc = [100 * sc.loc[m, 'frac_abs_z_lt_2'] for m, _ in OBS]
    fm = [100 * sm.loc[m, 'frac_abs_z_lt_2'] for m, _ in OBS]
    nc, nm = int(sc.loc['rg', 'n']), int(sm.loc['rg', 'n'])
    ax.axhline(95.45, color='0.5', lw=0.7, ls='--', zorder=1)
    ax.text(len(OBS) - 0.45, 95.75, 'sampling noise alone', ha='right', va='bottom', fontsize=6.5, color='0.35')
    for xx, fc_, fm_ in zip(x, fc, fm):
        ax.plot([xx - 0.15, xx + 0.15], [fc_, fm_], color='0.8', lw=0.6, zorder=2)
    ax.plot(x - 0.15, fc, 'o', ms=4.2, color=COL['CALVADOS-2'], zorder=3, label=f'CALVADOS-2 ({nc:,} chains)')
    ax.plot(x + 0.15, fm, 'o', ms=4.2, color=COL['Mpipi-GG'], zorder=3,
            label=f'Mpipi-GG ({nm:,} chains)')
    ax.set(ylim=(90, 100), ylabel='Chains within two\nstandard errors (%)', xticks=x, xlim=(-0.6, len(OBS) - 0.4))
    ax.set_yticks(range(90, 101, 2))
    ax.set_xticklabels([n for _, n in OBS], rotation=45, ha='right', rotation_mode='anchor')
    ax.legend(loc='lower left', handlelength=1.0, borderaxespad=0.2, ncol=1)
    return fc, fm, nc, nm


def panel_chi_by_nu(ax):
    """Median chi(Rg) at 300 K of the plain force fields by own-nu bin."""
    a = pd.read_csv(f'{TLA}/Analysis_FDT_RECHARGED/alongside/A_chi_by_nu.csv')
    a = a[a.own_nu_bin != 'all']
    for m in COL:
        g = a[a.model == m].set_index('own_nu_bin').reindex(list(BIN_MID))
        ax.plot(list(BIN_MID.values()), g.chi_median, 'o-', ms=2.6, lw=1.0, color=COL[m], label=m)
    ax.axhline(0, color='0.6', lw=0.5)
    ax.axvline(0.5, color='0.6', lw=0.5, ls=':')
    ax.set(xlabel=f'Compaction bin, {NU}', ylabel=f'Median {CHI} at\n300 K (% per K)')
    ax.set_xticks(list(BIN_MID.values()))
    ax.set_xticklabels(BIN_LAB, rotation=40, ha='right', rotation_mode='anchor')
    ax.set_xlim(0.41, 0.67)
    ax.legend(loc='upper right', handlelength=1.4, borderaxespad=0.0)


KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8,
          P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)       # Kyte-Doolittle hydropathy


def _seq_table(ids):
    """Charge fractions and mean Kyte-Doolittle hydropathy (rescaled to 0-1, as Uversky 2000) of each sequence."""
    import recharged_compare as RC
    b = pd.read_csv(RC.BIO_CSV, usecols=['UniProt_ID', 'sequence'], low_memory=False).set_index('UniProt_ID')
    b = b.loc[b.index.intersection(ids)]
    out = []
    for uid, q in b.sequence.items():
        q = ''.join(c for c in str(q).upper() if c in KD)
        n = len(q)
        fp = sum(q.count(c) for c in 'KR') / n
        fn = sum(q.count(c) for c in 'DE') / n
        out.append((uid, fp, fn, np.mean([(KD[c] + 4.5) / 9.0 for c in q])))
    return pd.DataFrame(out, columns=['UniProt_ID', 'fp', 'fn', 'H']).set_index('UniProt_ID')


def panel_states(ax, t):
    """Das-Pappu diagram of states (fractions of positive and negative residues) with the share in each region."""
    ax.scatter(t.fp, t.fn, s=0.6, alpha=0.3, color=COL['CALVADOS-2'], lw=0, rasterized=True)
    fcr, ncpr = t.fp + t.fn, t.fp - t.fn
    reg = np.select([fcr < 0.25, fcr <= 0.35, ncpr.abs() <= 0.35, ncpr < 0], ['R1', 'R2', 'R3', 'R4'], 'R5')
    share = pd.Series(reg).value_counts(normalize=True) * 100
    g = dict(color='0.35', lw=0.6)
    for f in (0.25, 0.35):
        ax.plot([0, f], [f, 0], **g)
    ax.plot([0.35, 0.675], [0.0, 0.325], **g)            # |NCPR| = 0.35 for FCR above 0.35
    ax.plot([0.0, 0.325], [0.35, 0.675], **g)
    ax.plot([0, 1], [1, 0], color='0.6', lw=0.5, ls=':')
    pos = {'R1': (0.07, 0.07), 'R2': (0.2, 0.2), 'R3': (0.36, 0.36), 'R4': (0.08, 0.62), 'R5': (0.62, 0.08)}
    for r, (x, y) in pos.items():
        ax.text(x, y, f'{r}\n{share.get(r, 0):.0f}%', fontsize=6.5, ha='center', va='center', color=INK, linespacing=1.0,
                bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.3))
    ax.set(xlim=(0, 0.75), ylim=(0, 0.75), xlabel='Fraction positive ($f_+$)', ylabel='Fraction negative ($f_-$)')
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_yticks([0, 0.25, 0.5, 0.75])
    return share


def panel_uversky(ax, t):
    """Charge and hydropathy plot (Uversky 2000): mean net charge against mean hydropathy, with the boundary
    <R> = 2.785 <H> - 1.151 that separates ordered from disordered proteins."""
    R = (t.fp - t.fn).abs()
    ax.scatter(t.H, R, s=0.6, alpha=0.3, color=COL['CALVADOS-2'], lw=0, rasterized=True)
    h = np.array([0.413, 0.65])
    ax.plot(h, 2.785 * h - 1.151, color=INK, lw=0.8)
    dis = float((R > 2.785 * t.H - 1.151).mean() * 100)
    ax.text(0.05, 0.95, f'disordered side\n{dis:.1f}%', transform=ax.transAxes, fontsize=6.5, va='top', color=INK,
            linespacing=1.0)
    ax.text(0.97, 0.06, 'ordered side', transform=ax.transAxes, fontsize=6.5, ha='right', color='0.4')
    ax.set(xlim=(0.15, 0.65), ylim=(0, 0.7), xlabel='Mean hydropathy', ylabel='Mean net charge')
    return dis


@house
def fig1():
    """Fig. 2 of the manuscript (named fig1 for history): data set (a-d), worked example (e, f) and validation of the
    fluctuation relation (g-i). Since 2026-10-31 the approach schematic and the counts per group, which Fig. 1
    (fig_overview.py) now shows, are replaced by the diagram of states (a) and the charge and hydropathy plot (b)."""
    D = load_models()
    fig = plt.figure(figsize=(W, 4.95))
    cols, w = (0.47, 2.66, 4.85), 1.60
    rows = {0: (3.62, 0.95), 1: (2.10, 0.95), 2: (0.62, 0.95)}
    t = _seq_table(D['C'].index)
    assert len(t) == len(D['C']), (len(t), len(D['C']))
    ax_a = ax = ax_in(fig, cols[0], rows[0][0], w, rows[0][1] + 0.12)
    share = panel_states(ax, t)
    ax_b = ax = ax_in(fig, cols[1], rows[0][0], w, rows[0][1] + 0.12)
    dis = panel_uversky(ax, t)
    # c lengths
    ax_c = ax = ax_in(fig, cols[2], rows[0][0], w, rows[0][1] + 0.12)
    L = D['C'].length.dropna()
    ax.hist(L, bins=np.logspace(np.log10(30), np.log10(1200), 40), color=COL['CALVADOS-2'], alpha=0.85, lw=0)
    ax.set_xscale('log')
    ax.axvline(L.median(), color='k', lw=0.7, ls='--')
    ax.set_ylim(0, ax.get_ylim()[1] * 1.08)
    ax.set(xlabel='Chain length (residues)', ylabel='Chains')
    ax.legend(handles=[Line2D([], [], color='k', lw=0.7, ls='--', label=f'median, {L.median():.0f}')], loc='upper right',
              handlelength=1.6, borderaxespad=0.0)
    # d nu
    ax_d = ax = ax_in(fig, cols[0], rows[1][0], w, rows[1][1])
    bins = np.linspace(0.30, 0.75, 46)
    nnu = {}
    for k, m in MODELS.items():
        v = D[k].nu[D[k].nu > 0.30]          # nu <= 0.30 are poor fits, left out as in the binned summaries
        nnu[m] = len(v)
        ax.hist(v, bins=bins, density=True, histtype='step', lw=1.0, color=COL[m], label=m)
    ax.set_ylim(0, ax.get_ylim()[1] * 2.2)                 # room for the legend
    for v, t, ha, dx in ((0.5, 'θ', 'right', -0.006), (0.59, 'good solvent', 'left', 0.006)):
        ax.axvline(v, color='0.5', lw=0.6, ls=':')
        ax.text(v + dx, 0.995, t, transform=ax.get_xaxis_transform(), fontsize=6.5, color='0.35', va='top', ha=ha)
    ax.set(xlabel=f'Apparent scaling exponent, {NU}', ylabel='Density', xlim=(0.3, 0.75))
    ax.legend(loc='upper right', handlelength=1.0, borderaxespad=0.1, bbox_to_anchor=(1.0, 0.88))   # force fields (Fig. 2b used to name them)
    # e, f worked example
    ax_e = ax = ax_in(fig, cols[1], rows[1][0], w, rows[1][1])
    example_dist(ax)
    ax_f = ax = ax_in(fig, cols[2], rows[1][0], w, rows[1][1])
    example_curve(ax)
    # g-i validation
    ax_g = ax = ax_in(fig, cols[0], rows[2][0], w, rows[2][1])
    n_out, n_all = panel_scatter(ax)
    ax_h = ax = ax_in(fig, cols[1], rows[2][0], w, rows[2][1])
    fc, fm, nc, nm = panel_agree(ax)
    ax_i = ax = ax_in(fig, cols[2], rows[2][0], w, rows[2][1])
    panel_chi_by_nu(ax)
    letters(fig, list(zip((ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g, ax_h, ax_i), 'abcdefghi')))
    print('fig1 regions %', share.round(1).to_dict(), 'disordered side %', round(dis, 1))
    print('fig1 nu > 0.30', nnu, 'median length', L.median())
    print('fig1 h % within 2 SE C2', np.round(fc, 1), f'(n {nc}) GG', np.round(fm, 1), f'(n {nm})')
    print(f'fig1 g: {n_out} of {n_all} chains outside the axes')
    save(fig, 'fig1_dataset_validation')


# ----------------------------------------------------------------------------- Fig. S1
@house
def ed_parts():
    """Fig. S1: contributions of each energy term to the susceptibility of Rg by compaction, per force field."""
    parts = pd.read_csv(f'{TLA}/Analysis_FDT_RECHARGED/alongside/A_chi_by_nu.csv')
    parts = parts[parts.own_nu_bin != 'all']
    fig = plt.figure(figsize=(W, 2.6))
    names = {'S_median': 'stickers', 'Q_median': 'charges (Boltzmann)', 'dQ_median': r'charges ($\partial H/\partial T$)',
             'B_median': 'bonds'}
    styles = {'S_median': ('-', 'o'), 'Q_median': ('-', 's'), 'dQ_median': ('--', 'v'), 'B_median': (':', 'D')}
    xs = list(BIN_MID.values())
    axes = []
    for i, (k, m) in enumerate(MODELS.items()):
        ax = ax_in(fig, 0.62 + i * 2.05, 0.66, 1.78, 1.50)
        p = parts[parts.model == m].set_index('own_nu_bin').reindex(list(BIN_MID))
        for c in names:
            ls, mk = styles[c]
            ax.plot(xs, p[c], ls=ls, marker=mk, ms=2.6, lw=0.9, color=COL[m], alpha=1.0 if c == 'S_median' else 0.6)
        ax.plot(xs, p.chi_median, color='k', lw=1.2)
        ax.axhline(0, color='0.6', lw=0.5, zorder=0)
        ax.set_xticks(xs)
        ax.set_xticklabels(BIN_LAB, rotation=40, ha='right', rotation_mode='anchor')
        ax.set_xlim(0.41, 0.67)
        ax.set_title(m, color=COL[m])
        ax.set_xlabel(f'Compaction bin, {NU}')
        if i == 0:
            ax.set_ylabel(f'Median contribution\nto {CHI} (% per K)')
        lo, hi = np.nanmin(p[list(names) + ['chi_median']].to_numpy()), np.nanmax(p[list(names) + ['chi_median']].to_numpy())
        pad = 0.06 * (hi - lo)
        ax.set_ylim(lo - pad, hi + pad)
        axes.append(ax)
    h = [Line2D([], [], color='0.35', ls=styles[c][0], marker=styles[c][1], ms=2.6, lw=0.9, label=names[c]) for c in names]
    h.append(Line2D([], [], color='k', lw=1.2, label=f'{CHI}'))
    fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=5, handlelength=2.2, borderaxespad=0.15)
    letters(fig, list(zip(axes, 'abc')))
    save(fig, 'fig_ed_parts')


# ----------------------------------------------------------------------------- Fig. S16
@house
def ogt():
    """Fig. S16: growth temperature (a BENDER bacteria within phyla, b independent proteomes, c composition)."""
    import ogt_response as OR
    import total_response_evidence as E
    D, _ = E.load()
    bio = pd.read_csv(E.RC.BIO_CSV, usecols=['UniProt_ID', 'organism', 'taxonomy_id', 'kingdom'], low_memory=False).set_index('UniProt_ID')
    bio['OGT'] = OR.ogt_table(bio)
    bio['phylum'] = OR.ogt_table(bio, 'phylum')
    sl = pd.read_csv(f'{TA}/bender_ogt/slopes.csv')
    sl = sl[(sl.quantity == 'total') & (sl.fixed_effects == 'phylum (bacteria only)')].set_index('model')
    fig = plt.figure(figsize=(W, 2.72))
    # a: organism means within phyla, phyla with >= 3 organisms (ogt_response.main)
    ax_a = ax = ax_in(fig, 0.55, 0.42, 1.72, 1.80)
    h, rep = [], {}
    for k in ('R', 'C'):                      # CALVADOS-2 drawn on top
        d = pd.DataFrame({'total': D[k].total}).join(bio, how='inner').dropna(subset=['OGT'])
        d = d[d.kingdom == 'Bacteria'].dropna(subset=['phylum'])
        o = d.groupby(['organism', 'phylum']).agg(OGT=('OGT', 'first'), chains=('total', 'size'),
                                                  total=('total', 'mean')).reset_index()
        o = o[o.phylum.isin(o.phylum.value_counts()[lambda s: s >= 3].index)]
        rx = o.OGT - o.groupby('phylum').OGT.transform('mean')
        ry = o.total - o.groupby('phylum').total.transform('mean')
        b = np.polyfit(rx, ry, 1)
        m = MODELS[k]
        ax.scatter(rx, ry, s=3.5, alpha=0.4 if k == 'C' else 0.3, color=COL[m], lw=0, rasterized=True, zorder=2 + (k == 'C'))
        xx = np.array([rx.min(), rx.max()])
        ax.plot(xx, np.polyval(b, xx), color=COL[m], lw=1.1, zorder=4 + (k == 'C'))
        rep[m] = (len(o), int(o.chains.sum()), 10 * b[0], int(sl.loc[m, 'organisms']), int(sl.loc[m, 'chains']),
                  float(sl.loc[m, 'slope_per_10C']))
        h.append(Line2D([], [], color=COL[m], lw=1.1, label=f'{m}, {len(o)} organisms, {num(10 * b[0], "+.3f")}'))
    print('ogt a (organisms, chains, slope) figure vs slopes.csv:', rep)
    ax.axhline(0, color='0.6', lw=0.5, zorder=1)
    ax.set(xlabel='Growth temperature, within phylum (°C)', ylabel='Total response, within\nphylum (% per K)')
    ax.legend(handles=h[::-1], loc='lower left', bbox_to_anchor=(-0.08, 1.0), title='BENDER bacteria, slope per 10 °C',
              alignment='left', handlelength=1.3, borderaxespad=0.2)
    # b: independent proteomes, CALVADOS-2 predictor
    ax_b = ax = ax_in(fig, 3.32, 0.74, 1.36, 1.48)
    A = pd.read_csv(f'{TA}/predictions/A_matched.csv')
    t = A[A['class'] == 'thermophile']
    order = t.groupby('clade').size().sort_values(ascending=False).index.tolist()
    short = {'Fungi (Sordariomycetes, Eurotiomycetes)': 'Fungi', 'Deinococcus-Thermus': 'Deinococcota'}
    data = [t.loc[t.clade == c, 'diff_total_C'].to_numpy() for c in order]
    ax.boxplot(data[::-1], orientation='horizontal', widths=0.55, showfliers=False,
               medianprops={'color': COL['CALVADOS-2'], 'lw': 1.3}, boxprops={'lw': 0.6}, whiskerprops={'lw': 0.6},
               capprops={'lw': 0.6})
    beyond = 0
    for v in data:
        q1, q3 = np.percentile(v, [25, 75])
        beyond += int(((v < q1 - 1.5 * (q3 - q1)) | (v > q3 + 1.5 * (q3 - q1))).sum())
    ax.axvline(0, color='0.4', lw=0.6, zorder=0)
    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels([f'{short.get(c, c)} ({len(v):,})' for c, v in zip(order, data)][::-1])
    ax.set_xlabel('Predicted total response, thermophile\nminus mesophile orthologue (% per K)')
    ax.set_title(f'Independent proteomes\nCALVADOS-2 predictor\nmedian {num(np.median(t.diff_total_C), "+.3f")}',
                 linespacing=1.1)
    ax.text(0.5, -0.38, f'{beyond:,} of {len(t):,} groups beyond the whiskers not shown', transform=ax.transAxes,
            ha='center', va='top', fontsize=6.5, color='0.35')
    print('ogt b median', np.median(t.diff_total_C), 'groups', len(t), 'beyond whiskers', beyond)
    # c: composition slopes
    ax_c = ax = ax_in(fig, 5.33, 0.42, 1.10, 1.80)
    terms = [('f_R', 'Arg'), ('f_pi', r'$\pi$ (FWYR)'), ('f_P', 'Pro'), ('f_E', 'Glu'), ('f_D', 'Asp'), ('f_S', 'Ser'),
             ('f_K', 'Lys'), ('kappa', r'$\kappa$')]
    fc = pd.read_csv(f'{TA}/bender_ogt/followup_composition.csv')
    fc = fc[(fc.model == 'CALVADOS-2') & fc.analysis.str.startswith('descriptor')].set_index('term')
    cr = pd.read_csv(f'{TA}/bender_ogt/checks_regions.csv')
    cr = cr[cr.check.str.startswith('proteome regions')].set_index('term')
    y = np.arange(len(terms))[::-1]
    for off, (src, lab, col, mk) in zip((0.15, -0.15), ((fc, 'BENDER bacteria', INK, 'o'),
                                                        (cr, 'independent proteomes', '#B07D2B', 's'))):
        e = src.reindex([t_ for t_, _ in terms])
        est = e['slope_per_10C'] if 'slope_per_10C' in e else e['slope_SD_per_10C']
        ax.errorbar(est, y + off, xerr=[est - e.ci_lo, e.ci_hi - est], fmt=mk, ms=2.8, color=col, lw=0.7, capsize=0,
                    label=lab)
    ax.axvline(0, color='0.6', lw=0.5, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in terms])
    ax.set_ylim(-0.6, len(terms) - 0.4)
    ax.set_xlabel('Slope per 10 °C (SD)')
    ax.legend(loc='lower left', bbox_to_anchor=(-0.5, 1.0), borderaxespad=0.2, handlelength=1.0)
    letters(fig, [(ax_a, 'a'), (ax_b, 'b'), (ax_c, 'c')])
    save(fig, 'fig_ogt')


# ----------------------------------------------------------------------------- Fig. S3
PAIR_TYPES = ['aro-aro', 'R-aro', 'K-aro', 'R-other', 'K-other', 'aro-other', 'other-other']


def pair_medians():
    """Median chi per pair by type (pairs >= 20) and median within-protein Arg / Lys ratio, per model: CALVADOS-2 and
    Mpipi-GG from the pair_type_split sample, Mpipi-Recharged for the same proteins (recharged_traj_analysis.main)."""
    pt = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/compact/pair_type_split/per_protein.csv')
    pt = pt[pt.status == 'ok']
    P = pd.read_csv(f'{TLA}/Analysis_FDT_RECHARGED/traj_pass/pairs.csv').drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    out = {}
    for m, src in (('CALVADOS-2', pt[pt.model == 'CALVADOS-2']), ('Mpipi-GG', pt[pt.model == 'Mpipi-GG']),
                   ('Mpipi-Recharged', P.loc[P.index.isin(pt.UniProt_ID)].reset_index())):
        r = {'proteins': len(src)}
        for t in PAIR_TYPES:
            r[t] = src.loc[src[f'n_pairs_{t}'] >= 20, f'chi_per_pair_{t}'].median()
        for x in ('aro', 'other'):
            ok = (src[f'n_pairs_R-{x}'] >= 20) & (src[f'n_pairs_K-{x}'] >= 20)
            r[f'R/K {x}'] = (src.loc[ok, f'chi_per_pair_R-{x}'] / src.loc[ok, f'chi_per_pair_K-{x}']).median()
            r[f'n R/K {x}'] = int(ok.sum())
        out[m] = r
    return out


@house
def pairs():
    """Fig. S3: sticker susceptibility per pair by pair type, and Arg / Lys per pair, three models."""
    Pm = pair_medians()
    s = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/compact/pair_type_split/summary.csv').set_index(['model', 'type'])
    for m in ('CALVADOS-2', 'Mpipi-GG'):          # same numbers as the analysis summary
        assert np.isclose(Pm[m]['R/K aro'], s.loc[(m, 'ratio R-aro / K-aro per pair'), 'median_chi_per_pair'])
        assert np.isclose(Pm[m]['R/K other'], s.loc[(m, 'ratio R-other / K-other per pair'), 'median_chi_per_pair'])
    lab = ['aromatic–\naromatic', 'Arg–\naromatic', 'Lys–\naromatic', 'Arg–\nother', 'Lys–\nother', 'aromatic–\nother',
           'other–\nother']
    fig = plt.figure(figsize=(W, 2.3))
    ax_a = ax = ax_in(fig, 0.62, 0.42, 3.55, 1.65)
    x = np.arange(len(PAIR_TYPES))
    for o, m in zip((-0.27, 0, 0.27), COL):
        ax.bar(x + o, [Pm[m][t] for t in PAIR_TYPES], 0.26, color=COL[m], label=f'{m}, {Pm[m]["proteins"]} proteins')
    ax.set_yscale('log')
    ax.set_ylim(0.08, 60)
    ax.set_yticks([0.1, 1, 10])
    ax.set_yticklabels(['0.1', '1', '10'])
    ax.set_xticks(x)
    ax.set_xticklabels(lab)
    ax.set_xlim(-0.55, len(PAIR_TYPES) - 0.45)
    ax.set_ylabel('Sticker susceptibility per pair\n(relative to the average pair)')
    ax.legend(loc='upper right', handlelength=1.0, borderaxespad=0.0)
    ax_b = ax = ax_in(fig, 4.82, 0.42, 1.62, 1.65)
    xx = np.arange(2)
    for o, m in zip((-0.27, 0, 0.27), COL):
        v = [Pm[m]['R/K aro'], Pm[m]['R/K other']]
        ax.bar(xx + o, v, 0.26, color=COL[m])
        for xi, vi in zip(xx + o, v):
            ax.text(xi, vi * 1.1, f'{vi:.1f}', ha='center', va='bottom', fontsize=6.5)
    ax.axhline(1, color='0.4', lw=0.6, ls='--', zorder=0)
    ax.set_yscale('log')
    ax.set_ylim(0.7, 40)
    ax.set_yticks([1, 10])
    ax.set_yticklabels(['1', '10'])
    ax.set_xticks(xx)
    ax.set_xticklabels(['with\naromatics', 'with other\nresidues'])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel('Arg / Lys per pair\n(same protein)')
    ax.legend(handles=[Line2D([], [], color='0.4', lw=0.6, ls='--', label='equality')], loc='upper right',
              handlelength=1.6, borderaxespad=0.0)
    letters(fig, [(ax_a, 'a'), (ax_b, 'b')])
    for m, r in Pm.items():
        print(f"pairs {m}: proteins {r['proteins']}, R/K aro {r['R/K aro']:.3f} (n {r['n R/K aro']}), "
              f"R/K other {r['R/K other']:.3f} (n {r['n R/K other']}); per pair",
              {t: round(r[t], 3) for t in PAIR_TYPES})
    save(fig, 'fig_ed_pairs')


# ----------------------------------------------------------------------------- Fig. S5
@house
def hotspots():
    """Fig. S5: median contribution of each residue type to the total response, CALVADOS-2 and Mpipi-Recharged."""
    C = pd.read_csv(f'{TR}/heat_crosses/analysis_residue_types.csv').set_index('residue').total_share_median
    R = pd.read_csv(f'{TLA}/Analysis_FDT_RECHARGED/traj_pass_capped/analysis_residue_types.csv').set_index('residue').total_share   # calibrated model (recharged_residues_capped.py)
    order = C.sort_values().index.tolist()
    fig = plt.figure(figsize=(W, 2.1))
    ax = ax_in(fig, 0.62, 0.40, 5.83, 1.62)
    x = np.arange(len(order))
    ax.bar(x - 0.2, 1e3 * C.reindex(order), 0.38, color=COL['CALVADOS-2'], label='CALVADOS-2')
    ax.bar(x + 0.2, 1e3 * R.reindex(order), 0.38, color=COL['Mpipi-Recharged'], label='Mpipi-Recharged')
    ax.axhline(0, color='0.4', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_xlim(-0.6, len(order) - 0.4)
    lo = 1e3 * min(C.min(), R.min())
    hi = 1e3 * max(C.max(), R.max(), 0)                    # positive contributions (e.g. aromatic Boltzmann shares) shown too
    ax.set_ylim(np.floor(lo) - 0.5, max(0, np.ceil(hi) + 0.5))
    ax.set_xlabel('Residue type')
    ax.set_ylabel(f'Contribution per residue\n(10$^{{-3}}$ % of mean {RG} per K)')
    rho = stats.spearmanr(C.reindex(order), R.reindex(order))[0]
    ax.legend(loc='lower right', title=f'rank agreement {rho:.2f}', alignment='left', borderaxespad=0.1)
    print('hotspots rank agreement', rho)
    save(fig, 'fig_ed_hotspots')


# ----------------------------------------------------------------------------- Fig. S4
@house
def dynamics():
    """Fig. S4: rank partial correlations (nu and log N fixed) of the dynamics of Rg with chi and with the total."""
    import rg_dynamics as RD
    import total_response_evidence as E
    D, _ = E.load()
    per = pd.read_csv(f'{TR}/dynamics/per_chain.csv')
    stored = pd.read_csv(f'{TR}/dynamics/correlations.csv').set_index(['model', 'dynamics', 'response'])
    feats = [('rg_cv', 'breathing\namplitude'), ('tau_int_frames', 'relaxation\ntime'),
             ('switches_per_1000', 'compact–expanded\nswitching'), ('bimodality', 'bimodality')]
    c, n = {}, {}
    for k, m in MODELS.items():                   # rg_dynamics.main on the current responses
        x = per[per.model == m].set_index('UniProt_ID').join(D[k][['chi', 'total', 'nu']], how='inner')
        x = x.dropna(subset=['chi', 'total', 'nu', 'N', 'tau_int_frames', 'switches_per_1000', 'bimodality', 'rg_cv'])
        n[m] = len(x)
        for f, _ in feats:
            for resp in ('chi', 'total'):
                c[(m, f, resp)] = RD.partial(x[resp].to_numpy(), x[f].to_numpy(), [x.nu.to_numpy(), np.log(x.N).to_numpy()])
    dif = {key: round(v - stored.loc[key, 'partial | nu, logN'], 4) for key, v in c.items()
           if abs(v - stored.loc[key, 'partial | nu, logN']) > 5e-4}
    print('dynamics chains', n, '| differs from stored correlations.csv (n there',
          stored.groupby('model').n.first().to_dict(), '):', dif)
    fig = plt.figure(figsize=(W, 2.3))
    axes = []
    for i, (resp, title) in enumerate((('chi', 'With the fixed-solvent susceptibility'), ('total', 'With the total response'))):
        ax = ax_in(fig, 0.62 + i * 2.97, 0.42, 2.85, 1.55, sharey=axes[0] if axes else None)
        xs = np.arange(len(feats))
        for o, m in zip((-0.27, 0, 0.27), COL):
            ax.bar(xs + o, [c[(m, f, resp)] for f, _ in feats], 0.26, color=COL[m], label=m)
        ax.axhline(0, color='0.4', lw=0.5)
        ax.set_xticks(xs)
        ax.set_xticklabels([lab for _, lab in feats])
        ax.set_xlim(-0.55, len(feats) - 0.45)
        ax.set_title(title)
        if i:
            ax.tick_params(labelleft=False)
        axes.append(ax)
    axes[0].set_ylim(-0.55, 0.6)
    axes[0].set_ylabel(f'Rank partial correlation\n({NU} and log $N$ fixed)')
    axes[0].legend(loc='lower left', handlelength=1.0, borderaxespad=0.0)
    letters(fig, [(axes[0], 'a'), (axes[1], 'b')])
    save(fig, 'fig_ed_dynamics')


# ----------------------------------------------------------------------------- Fig. S7
def _cross(p1, p2, q1, q2):
    """Whether segments p1-p2 and q1-q2 cross."""
    def orient(a, b, c):
        return np.sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    return orient(p1, p2, q1) != orient(p1, p2, q2) and orient(q1, q2, p1) != orient(q1, q2, p2)


def place_labels(ax, xs, ys, names, fontsize=6.5, avoid_lines=(), manual=None, near=3.2, far=(8.0, 12.0, 16.0, 20.0)):
    """Label scatter points without overlaps (adjustText-like, deterministic): for each label (most crowded first, three
    passes), the candidate position around its point (16 directions; next to the point, or further out with a thin
    leader line) with the least overlap with other labels, points, the listed lines (x and y arrays) and leader lines,
    and inside the axes. manual: {name: (dx_pt, dy_pt)} fixes a label offset."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px = fig.dpi / 72
    P = ax.transData.transform(np.column_stack([xs, ys]))
    rad = 2.6 * px                                  # marker radius (s = 12) plus a margin
    gap = 1.0 * px                                  # minimum clearance between labels
    L = [ax.transData.transform(np.column_stack(seg)) for seg in avoid_lines]
    L = np.vstack([np.vstack([np.linspace(a, b, 400) for a, b in zip(s_[:-1], s_[1:])]) for s_ in L]) if L else np.zeros((0, 2))
    axb = ax.bbox
    ang = np.deg2rad(np.arange(16) * 22.5 + 45)
    cands = [(r * np.cos(a), r * np.sin(a)) for r in (near, *far) for a in ang]
    texts = {nm: ax.annotate(nm, (x, y), xytext=(0, 0), textcoords='offset points', fontsize=fontsize, zorder=4)
             for nm, x, y in zip(names, xs, ys)}
    sizes = {nm: texts[nm].get_window_extent(rend) for nm in names}
    crowd = [-(np.hypot(*(P - P[i]).T) < 30 * px).sum() for i in range(len(names))]
    order = np.argsort(crowd, kind='stable')
    placed, lead = {}, {}

    def align(dx, dy):
        ha = 'center' if abs(dx) < 0.38 * np.hypot(dx, dy) else ('left' if dx > 0 else 'right')
        va = 'center' if abs(dy) < 0.38 * np.hypot(dx, dy) else ('bottom' if dy > 0 else 'top')
        return ha, va

    def box(i, dx, dy):
        w, h = sizes[names[i]].width, sizes[names[i]].height
        ha, va = align(dx, dy)
        x0 = P[i, 0] + dx * px - {'left': 0, 'center': w / 2, 'right': w}[ha]
        y0 = P[i, 1] + dy * px - {'bottom': 0, 'center': h / 2, 'top': h}[va]
        return np.array([x0, y0, x0 + w, y0 + h])

    def inter(a, b, g=0.0):
        return max(0, min(a[2], b[2]) + g - max(a[0], b[0])) * max(0, min(a[3], b[3]) + g - max(a[1], b[1]))

    def inside(pts, b):
        return ((pts[:, 0] > b[0]) & (pts[:, 0] < b[2]) & (pts[:, 1] > b[1]) & (pts[:, 1] < b[3])).sum()

    for _ in range(3):
        for i in order:
            nm = names[i]
            best = None
            opts = [manual[nm]] if manual and nm in manual else cands
            for dx, dy in opts:
                b = box(i, dx, dy)
                r = np.hypot(dx, dy)
                cost = 0.02 * r ** 1.5
                cost += 5.0 * sum(inter(b, placed[o], gap) for o in placed if o != nm)
                cost += 5.0 * sum(inter(b, [p[0] - rad, p[1] - rad, p[0] + rad, p[1] + rad]) for p in P)
                if len(L):
                    cost += 3.0 * inside(L, b)
                if r <= near + 0.5:                 # no leader: the label centre clearly nearer its own point
                    dist = np.hypot(P[:, 0] - (b[0] + b[2]) / 2, P[:, 1] - (b[1] + b[3]) / 2)
                    rival = np.delete(dist, i).min() if len(P) > 1 else np.inf
                    cost += 40.0 * max(0.0, 1.35 * dist[i] - rival) / px
                else:                               # leader line from the point to the label
                    seg = np.linspace(P[i], [(b[0] + b[2]) / 2, (b[1] + b[3]) / 2], 30)[3:]
                    cost += 20.0 * sum(inside(seg, placed[o]) for o in placed if o != nm)
                    others = np.delete(P, i, axis=0)
                    cost += 20.0 * (np.min(np.hypot(*(seg[:, None, :] - others[None]).transpose(2, 0, 1)), axis=1) < rad).sum()
                    c0 = [(b[0] + b[2]) / 2, (b[1] + b[3]) / 2]
                    cost += 50.0 * sum(_cross(P[i], c0, P[o], c) for o, c in lead.items() if o != i)
                for o, c in lead.items():           # other leader lines through this label
                    if o != i:
                        cost += 20.0 * inside(np.linspace(P[o], c, 30)[3:], b)
                out = max(0, axb.x0 - b[0]) + max(0, b[2] - axb.x1) + max(0, axb.y0 - b[1]) + max(0, b[3] - axb.y1)
                cost += 1e4 * out
                if best is None or cost < best[0]:
                    best = (cost, dx, dy, b, r)
            placed[nm] = best[3]
            if best[4] > near + 0.5:
                lead[i] = [(best[3][0] + best[3][2]) / 2, (best[3][1] + best[3][3]) / 2]
            else:
                lead.pop(i, None)
            t = texts[nm]
            t.set_position((best[1], best[2]))
            ha, va = align(best[1], best[2])
            t.set_ha(ha)
            t.set_va(va)
    for i, nm in enumerate(names):                  # leader lines for the labels set apart from their point
        if i in lead:
            t = texts[nm]
            t.remove()
            dx, dy = t.get_position()
            ha, va = align(dx, dy)
            texts[nm] = ax.annotate(nm, (xs[i], ys[i]), xytext=(dx, dy), textcoords='offset points', fontsize=fontsize,
                                    ha=ha, va=va, zorder=4,
                                    arrowprops={'arrowstyle': '-', 'color': '0.55', 'lw': 0.4, 'shrinkA': 0.5,
                                                'shrinkB': 2.4})
    return texts


@house
def mutations():
    """Fig. S7: median change of the total by residue introduced (a) and replaced (b), CALVADOS-2 vs Mpipi-Recharged."""
    d = f'{TR}/mutation_scan'
    fig = plt.figure(figsize=(W, 3.4))
    axes = []
    for i, (f, key, title) in enumerate((('compare_by_mut.csv', 'mut', 'By residue introduced'),
                                         ('compare_by_wt.csv', 'wt', 'By residue replaced'))):
        ax = ax_in(fig, 0.52 + i * 3.26, 0.42, 2.62, 2.62)
        t = pd.read_csv(f'{d}/{f}').set_index(key)
        x, y = 1e3 * t.dtotal_C, 1e3 * t.dtotal_R
        ax.scatter(x, y, s=12, color='0.3', lw=0, zorder=3)
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        span = hi - lo
        lim = [lo - 0.1 * span, hi + 0.1 * span]
        ax.plot(lim, lim, color='0.6', lw=0.6, ls='--', zorder=1)
        ax.axhline(0, color='0.85', lw=0.5, zorder=0)
        ax.axvline(0, color='0.85', lw=0.5, zorder=0)
        ax.set(xlim=lim, ylim=lim, xlabel='CALVADOS-2 (10$^{-3}$ % per K)', ylabel='Mpipi-Recharged (10$^{-3}$ % per K)')
        ax.set_aspect('equal')
        ax.set_title(f'{title} (rank agreement {stats.spearmanr(x, y)[0]:.2f})')
        place_labels(ax, x.to_numpy(), y.to_numpy(), list(t.index),
                     avoid_lines=[(lim, lim), (lim, [0, 0]), ([0, 0], lim)])
        axes.append(ax)
    letters(fig, [(axes[0], 'a'), (axes[1], 'b')])
    save(fig, 'fig_ed_mutations')


if __name__ == '__main__':
    which = sys.argv[1:] or ['fig1', 'ed_parts', 'ogt', 'pairs', 'hotspots', 'dynamics', 'mutations']
    for w in which:
        globals()[w]()
