"""
fig_forcing.py — how the hydration forcing acts on one stored ensemble (Fig. 2), drawn entirely from the real
CALVADOS-2 300 K trajectory of the compact Gly-rich protein from tomato (UniProt Q02021), the example chain of Fig. 1.
Capped main model since 2026-10-25 (capped_proteome.py): hydrophilic term on D, E, K, R, each residue's contact number
c_i = -sum_j [U_AH(r_ij; 1) - U_AH(r_ij; 0)] / eps_AH counted through h(c) = c_sat (1 - exp(-c / c_sat)), c_sat = 2.

Steps shown (Methods Eqs. 5-7):
  a  the stored 300 K ensemble: Rg of the 1,000 analysed frames, with the two frames drawn in b marked
  b  a compact and an expanded frame (principal plane), contacts between charged residues dashed, with each frame's
     per-frame hydrophilic forcing at the calibrated strength
  c  the cap: h(c) against c (the linear forcing would count c itself), with the distribution of the contact numbers of
     the chain's charged residues over all frames
  d  per-frame forcing against Rg over all frames, whose covariance gives the hydration term
  e  the total response of Rg of this chain split into the fixed-solvent susceptibility and the two hydration terms

Every number is recomputed from the stored frames and checked against the stored capped per-protein results
(capped_300K/capped_dchi_calvados_300K.csv and total_response_evidence.load()) before plotting.

    BENDER_REANALYSIS_DIR=$CS_REPO/bender_reanalysis python fig_forcing.py
Output: manuscript/figs/fig_forcing.{pdf,png}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.expanduser((__import__('os').environ['CS_REPO'] + '/bender_reanalysis'))))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402
import tea_fdt as TF  # noqa: E402
import total_response_evidence as TR  # noqa: E402

UID = 'Q02021'
OUT = f'{FX.ROOT}/manuscript/figs/fig_forcing'
OUT_MERGED = f'{FX.ROOT}/manuscript/figs/fig_forcing_calibration'
CACHE = f'{FX.TLA}/Analysis_FDT/forcing_example_capped_{UID}.npz'
HYD = set('DEKR')                                      # the charged-group hydrophilic set of the capped model


def compute():
    """Per-frame Rg, capped D_hyd and D_tea (kJ/mol/K per unit strength), per-residue contact numbers and positions."""
    import capped_proteome as CP
    if os.path.exists(CACHE):
        return dict(np.load(CACHE, allow_pickle=True))
    L, fdt, bc = TF.L, TF.fdt, TF.bc
    d = FX.EX[UID]['dir']
    z = np.load(os.path.join(TF.NPZ, f'{UID}.npz'))
    xyz, seq, info = bc.load_calvados(os.path.join(d, f'{UID}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
    fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
    assert np.array_equal(fi, z['frame_index']), 'frames differ from the stored FDT frames'
    q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
    i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
    dls = {'hyd': TF.DL_CLASSES['hyd']}
    tea = {a: sum(TF.DL_CLASSES[c][a] for c in TR.TEA_CLASSES) for a in TF.DL_CLASSES['hyd']}
    dls['tea'] = tea
    sa = np.array(list(seq))
    dl = {k: 0.5 * (np.array([v[a] for a in seq])[i] + np.array([v[a] for a in seq])[j]) for k, v in dls.items()}
    n = len(xyz)
    Sres = np.zeros((n, len(seq)))
    M = CP.incidence(len(seq), i, j)
    for a in range(n):
        x = xyz[a].astype(float)
        r = np.linalg.norm(x[i] - x[j], axis=1)
        u1 = fdt.ah_pair_energy(r[None], s, np.ones_like(l), shift)[0]
        u0 = fdt.ah_pair_energy(r[None], s, np.zeros_like(l), shift)[0]
        Sres[a] = M.T @ (u1 - u0)
    cres = -Sres / fdt.EPS_AH                                   # contact number of each residue in each frame
    hc = CP.h(cres)
    D = {'hyd': -0.5 * fdt.EPS_AH * (hc @ np.array([CP.F_CHG[a] for a in seq])),
         'tea': -0.5 * fdt.EPS_AH * (hc @ np.array([CP.F_TEA[a] for a in seq]))}
    out = dict(rg=z['rg'].astype(float), D_hyd=D['hyd'], D_tea=D['tea'], xyz=xyz.astype(np.float32), seq=np.array(seq),
               pi=i, pj=j, sig=s, shift=shift, lam=l, c=cres.astype(np.float32))
    np.savez(CACHE, **out)
    return out


def check(c):
    """The recomputed hydration terms must reproduce the stored per-protein values."""
    RT = TF.fdt.R_GAS * TF.T0
    rg = c['rg']
    rel = lambda Dk: -100 * TF.fdt._cov(rg, Dk) / RT / rg.mean()   # noqa: E731  % per K per unit strength
    assert TR.FORCING == 'capped'
    st = pd.read_csv(os.path.join(TR.CAPPED, TR.CAPPED_FILES['C'])).set_index('UniProt_ID').loc[UID]
    assert abs(rel(c['D_hyd']) - st['dchi_rg_rel_chg_cap']) < 1e-6 * max(1, abs(st['dchi_rg_rel_chg_cap'])), \
        (rel(c['D_hyd']), st['dchi_rg_rel_chg_cap'])
    assert abs(rel(c['D_tea']) - st['dchi_rg_rel_tea_cap']) < 1e-6 * max(1, abs(st['dchi_rg_rel_tea_cap'])), \
        (rel(c['D_tea']), st['dchi_rg_rel_tea_cap'])
    Dt, g = TR.load()
    row = Dt['C'].loc[UID]
    parts = {'chi': float(row.chi), 'TEA': g['tea_C'] * rel(c['D_tea']), 'H': g['C'] * rel(c['D_hyd'])}
    assert abs(parts['TEA'] - row.TEA) < 1e-6 and abs(parts['H'] - row.H) < 1e-6, (parts, row)
    parts['total'] = parts['chi'] + parts['TEA'] + parts['H']
    assert abs(parts['total'] - row.total) < 1e-6
    # block standard errors (10 contiguous blocks), as everywhere in the paper
    se = {}
    jk = lambda f: jackknife(f, len(rg))   # noqa: E731
    se['H'] = jk(lambda m: g['C'] * -100 * TF.fdt._cov(rg[m], c['D_hyd'][m]) / RT / rg[m].mean())
    se['TEA'] = jk(lambda m: g['tea_C'] * -100 * TF.fdt._cov(rg[m], c['D_tea'][m]) / RT / rg[m].mean())
    chi, chi_se = FX._chi300(UID)
    print(f'check passed: chi {parts["chi"]:.3f} (stored chi_total {chi:.3f} ± {chi_se:.3f}), '
          f'TEA {parts["TEA"]:.4f} ± {se["TEA"]:.4f}, hyd {parts["H"]:.3f} ± {se["H"]:.3f}, total {parts["total"]:.3f}')
    return parts, se, chi_se, g


def jackknife(f, n, B=10):
    idx = np.array_split(np.arange(n), B)
    vals = np.array([f(np.concatenate([idx[b] for b in range(B) if b != k])) for k in range(B)])
    return float(np.sqrt((B - 1) / B * ((vals - vals.mean()) ** 2).sum()))



def response_se(c, g):
    """10-block jackknife SEs (% per K) of the fixed-solvent susceptibility, both hydration terms and the total."""
    z = np.load(os.path.join(TF.NPZ, f'{UID}.npz'))
    rg, G = c['rg'], z['G_T'].astype(float)
    RT = TF.fdt.R_GAS * TF.T0
    terms = {'chi': G, 'TEA': g['tea_C'] * c['D_tea'] / RT, 'H': g['C'] * c['D_hyd'] / RT}
    terms['total'] = terms['chi'] + terms['TEA'] + terms['H']
    val = {k: -100 * TF.fdt._cov(rg, v) / rg.mean() for k, v in terms.items()}
    se = {k: jackknife(lambda m, v=v: -100 * TF.fdt._cov(rg[m], v[m]) / rg[m].mean(), len(rg)) for k, v in terms.items()}
    return val, se


def contacts(x, c, frame_pairs=None):
    """Contacts |i-j| >= 3 closer than 1.5 times the pair's bead diameter (Supplementary Section S1)."""
    i, j, s = c['pi'], c['pj'], c['sig']
    r = np.linalg.norm(x[i] - x[j], axis=1)
    m = (np.abs(i - j) >= 3) & (r < 1.5 * s)
    return i[m], j[m]


C2_, HY_, OR_, GR_ = FX.COL['CALVADOS-2'], '#2E8C6F', '#C8792A', '0.55'


def frames(c, g):
    """Rg, per-frame forcing at the calibrated strength (kJ/mol/K) and the 5th / 95th percentile frames."""
    rg, Dh = c['rg'], c['D_hyd'] * g['C']
    k_c, k_e = int(np.argsort(rg)[len(rg) // 20]), int(np.argsort(rg)[-len(rg) // 20])
    return rg, Dh, k_c, k_e


def panel_a(a, rg, k_c, k_e, title=None, split=False):
    """The stored 300 K ensemble: Rg of each frame, with the two frames of b marked (split: the two labels set to the
    left and right of their lines, for a narrow panel)."""
    a.hist(rg, bins=30, color=C2_, alpha=0.85, edgecolor='white', linewidth=0.4)
    a.set_ylim(0, a.get_ylim()[1] * (1.15 if split else 1.0))
    for k, col, lab, ha in ((k_c, OR_, 'compact', 'right'), (k_e, GR_, 'expanded', 'left')):
        a.axvline(rg[k], color=col, lw=1.1, ls='--')
        a.text(rg[k] + ((-0.06 if ha == 'right' else 0.06) if split else 0), a.get_ylim()[1] * 0.97, lab, color=col,
               ha=ha if split else 'center', va='top', fontsize=6.5, bbox=dict(fc='white', ec='none', pad=0.6))
    a.set(xlabel=r'$R_\mathrm{g}$ of each frame (nm)', ylabel='Frames')
    a.set_title(title or f'Stored 300 K ensemble, {len(rg):,} frames', fontsize=7.5, pad=3)


def panel_b(fig, b, c, rg, Dh, k_c, k_e, legend=(('hyd', 'other', 'contact'),), y0=1.07, dy=0.085):
    """A compact and an expanded frame on their principal plane, with contacts between charged residues dotted.
    legend: the entries of each legend line above the panel (one tuple per line, from the top)."""
    seq = str(c['seq'])
    b.set_axis_off()
    off, labels, ymin = 0.0, [], np.inf
    for k, col, lab in ((k_c, OR_, 'compact frame'), (k_e, GR_, 'expanded frame')):
        x = c['xyz'][k].astype(float)
        x = x - x.mean(0)
        _, _, vt = np.linalg.svd(x, full_matrices=False)
        p = x @ vt[:2].T
        p[:, 0] += off - p[:, 0].min()
        b.plot(p[:, 0], p[:, 1], color='0.65', lw=0.6, zorder=1)
        hyd = np.array([r in HYD for r in seq])
        b.scatter(p[~hyd, 0], p[~hyd, 1], s=5, color='white', edgecolor='0.45', lw=0.4, zorder=3)
        b.scatter(p[hyd, 0], p[hyd, 1], s=5, color=HY_, edgecolor=HY_, lw=0.3, zorder=3)
        ci, cj = contacts(x, c)
        both = hyd[ci] & hyd[cj]
        for u, v in zip(ci[both], cj[both]):
            b.plot(p[[u, v], 0], p[[u, v], 1], color=HY_, lw=0.4, alpha=0.7, ls=(0, (1.5, 1.5)), zorder=2)
        labels.append(((p[:, 0].min() + p[:, 0].max()) / 2, f'{lab}\n' + r'$R_\mathrm{g}$' + f' = {rg[k]:.2f} nm, '
                       f'{both.sum()} charged contacts\n' + r'$\gamma_\mathrm{hyd}D_\mathrm{hyd}$' +
                       FX.minus(f' = {Dh[k]:.2f} kJ mol') + r'$^{-1}$K$^{-1}$', col))
        ymin = min(ymin, p[:, 1].min())
        off = p[:, 0].max() + 1.6
    for mid, txt, col in labels:
        b.text(mid, ymin - 0.35, txt, ha='center', va='top', fontsize=6.5, color=col)
    b.set_aspect('equal')
    b.autoscale_view()
    # legend above the panel, clear of both frames; entries placed from their measured text widths
    fig.canvas.draw()
    r_ = fig.canvas.get_renderer()
    entry = {'hyd': ('Charged residue (D, E, K, R)', HY_), 'other': ('Other residue', '0.4'),
             'contact': ('Contact between charged residues', HY_)}
    for row, kinds in enumerate(legend):
        y_, x_, gap = y0 - row * dy, 0.005, 0.028
        for kind in kinds:
            txt, col = entry[kind]
            if kind == 'contact':
                b.plot([x_ - 0.008, x_ + 0.025], [y_, y_], color=HY_, lw=0.6, ls=(0, (1.5, 1.5)),
                       transform=b.transAxes, clip_on=False)
                x_ += 0.03
            else:
                b.scatter([x_], [y_], s=6, transform=b.transAxes, clip_on=False,
                          **(dict(color=HY_) if kind == 'hyd' else dict(color='white', edgecolor='0.45', lw=0.4)))
                x_ += 0.015
            t = b.text(x_, y_, txt, transform=b.transAxes, color=col, fontsize=6.5, va='center')
            bb = t.get_window_extent(r_).transformed(b.transAxes.inverted())
            x_ = bb.x1 + gap
            assert bb.x1 <= 1.0, f'legend runs past the panel ({bb.x1:.3f})'


def panel_c(fig, cax, c):
    """The saturation: h(c) against c, with the contact numbers of this chain's charged residues."""
    import capped_proteome as CP
    seq = str(c['seq'])
    cc = np.asarray(c['c'], float)[:, np.array([r in HYD for r in seq])].ravel()
    xc = np.linspace(0, max(6.0, np.percentile(cc, 99.5)), 300)
    hist_ax = cax.twinx()
    hist_ax.hist(cc, bins=np.linspace(0, xc[-1], 40), color=HY_, alpha=0.18, edgecolor='none')
    hist_ax.set_yticks([])
    hist_ax.spines[['top', 'right']].set_visible(False)
    cax.set_zorder(hist_ax.get_zorder() + 1)
    cax.patch.set_visible(False)
    cax.plot(xc, xc, color='0.55', lw=0.9, ls=(0, (3, 2)), label='linear forcing')
    cax.plot(xc, CP.h(xc), color=HY_, lw=1.4, label='saturating forcing')
    cax.axhline(CP.CSAT, color='0.75', lw=0.5, ls=':')
    cax.set(xlabel=r'Contact number $c_i$ of a residue', ylabel=r'Counted contacts, $h(c_i)$', xlim=(0, xc[-1]),
            ylim=(0, xc[-1] * 0.75))
    cax.set_title('Contacts per residue saturate', fontsize=7.5, pad=3)
    cax.legend(loc='upper left', handlelength=1.6, borderaxespad=0.2)
    cax.text(0.97, 0.36, f'{100 * (cc > CP.CSAT).mean():.0f}% of charged residues\n(all frames) beyond ' +
             r'$c_\mathrm{sat}$ = 2', transform=cax.transAxes, ha='right', va='bottom', fontsize=6.5, color=HY_)


def panel_d(d, rg, Dh, k_c, k_e):
    """Per-frame forcing against Rg, whose covariance gives the hydration term."""
    d.scatter(rg, Dh, s=3, color=C2_, alpha=0.45, lw=0)
    for k, col in ((k_c, OR_), (k_e, GR_)):
        d.scatter(rg[k], Dh[k], s=16, color=col, edgecolor='white', lw=0.6, zorder=3)
    fit = np.polyfit(rg, Dh, 1)
    xx = np.array([rg.min(), rg.max()])
    d.plot(xx, np.polyval(fit, xx), color=INK_, lw=0.9, ls='--')
    d.set(xlabel=r'$R_\mathrm{g}$ of each frame (nm)', ylabel=r'$\gamma_\mathrm{hyd}D_\mathrm{hyd}$ (kJ mol$^{-1}$K$^{-1}$)')
    d.set_title('Forcing of each frame', fontsize=7.5, pad=3)
    d.text(0.97, 0.06, r'Cov($R_\mathrm{g}$, $D_\mathrm{hyd}$) > 0', transform=d.transAxes, ha='right', va='bottom',
           fontsize=6.5)


def panel_e(e, val, se, pad=0.42):
    """The total response of Rg split into its parts."""
    lab = ['fixed-solvent\nsusceptibility', 'hydrophobic\n(TEA) term', 'hydrophilic\nterm', 'total response']
    keys = ['chi', 'TEA', 'H', 'total']
    v = [val[k] for k in keys]
    cols = ['0.6', OR_, HY_, C2_]                  # TEA in its role colour (#C8792A), as in the d lambda / dT panel
    yy = np.arange(len(keys))[::-1]
    e.barh(yy, v, xerr=[se[k] for k in keys], color=cols, height=0.62, error_kw=dict(lw=0.7, capsize=1.6))
    e.axvline(0, color='0.3', lw=0.6)
    e.set_yticks(yy, lab, fontsize=6.5)
    e.set_xlabel(r'Response of $R_\mathrm{g}$ at 300 K (% per K)')
    e.set_xlim(min(v) - max(se.values()) - pad, max(v) + max(se.values()) + pad)
    for y_, v_, k in zip(yy, v, keys):
        dx = se[k] + 0.05
        e.text(v_ + (dx if v_ >= 0 else -dx), y_, FX.minus(f'{v_:+.2f}'), va='center',
               ha='left' if v_ >= 0 else 'right', fontsize=6.5)
    e.set_title('Response of this chain', fontsize=7.5, pad=3)


def plot(c, val, se, g):
    """The former stand-alone figure (fig_forcing), kept reproducible."""
    import matplotlib.gridspec as gs
    FX.house()
    rg, Dh, k_c, k_e = frames(c, g)
    fig = plt.figure(figsize=(FX.WIDTH, 4.15))
    G = gs.GridSpec(2, 3, figure=fig, left=0.075, right=0.985, top=0.95, bottom=0.11, hspace=0.62, wspace=0.42,
                    width_ratios=[1, 1, 1])
    a = fig.add_subplot(G[0, 0])
    panel_a(a, rg, k_c, k_e)
    b = fig.add_subplot(G[0, 1:])
    panel_b(fig, b, c, rg, Dh, k_c, k_e)
    cax = fig.add_subplot(G[1, 0])
    panel_c(fig, cax, c)
    d = fig.add_subplot(G[1, 1])
    panel_d(d, rg, Dh, k_c, k_e)
    e = fig.add_subplot(G[1, 2])
    panel_e(e, val, se)
    FX.add_letters(fig, [(a, 'a', {'row': 0}), (b, 'b', {'row': 0}), (cax, 'c', {'row': 1}), (d, 'd', {'row': 1}),
                         (e, 'e', {'row': 1, 'xaxes': [e]})])
    FX.save(fig, OUT)
    return k_c, k_e


def plot_merged(c, val, se, g):
    """Fig. 3 since 2026-10-30, the former Figs. 3 and 4 in one: the forcing on this chain (a-e), the change of stickiness
    per kelvin (f), the five smFRET heating curves (g-k) and the left-out slopes (l). The residue-set comparison went to
    the supplement (fig_extra.si_residue_sets). Axes placed in inches."""
    import matplotlib
    FX.house()
    rg, Dh, k_c, k_e = frames(c, g)
    W, H = FX.WIDTH, 5.18
    fig = plt.figure(figsize=(W, H))
    ax_in = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])  # noqa: E731
    # row 1: a, b, c
    a = ax_in(0.42, 3.83, 1.0, 1.05)
    panel_a(a, rg, k_c, k_e, title='300 K ensemble', split=True)
    b = ax_in(1.68, 3.55, 2.9, 1.36)
    panel_b(fig, b, c, rg, Dh, k_c, k_e, legend=(('hyd', 'other'), ('contact',)), y0=1.16, dy=0.105)
    cc = ax_in(5.08, 3.83, 1.38, 1.05)
    panel_c(fig, cc, c)
    # row 2: d, e, f
    d = ax_in(0.55, 2.2, 1.3, 0.92)
    panel_d(d, rg, Dh, k_c, k_e)
    d.set_ylabel(r'$\gamma_\mathrm{hyd}D_\mathrm{hyd}$' + '\n' + r'(kJ mol$^{-1}$K$^{-1}$)')
    e = ax_in(2.85, 2.2, 1.1, 0.92)
    panel_e(e, val, se, pad=0.5)
    f = ax_in(4.62, 2.2, 1.84, 0.92)
    FX.forcing_bars(f)
    f.set_title('Stickiness change per kelvin', fontsize=7.5, pad=3)
    f.set_ylabel(r'd$\lambda_i$/d$T$ ($10^{-3}$ K$^{-1}$)')
    # row 3: g-k heating curves, l left-out slopes
    xs = [0.5, 1.3, 2.1, 2.9, 4.0]
    cur = []
    for i, (x, (n, t)) in enumerate(zip(xs, FX.CONSTRUCTS)):
        ax = ax_in(x, 0.62, 0.68, 0.93)
        FX.wuttke_panel(ax, n, t, i == 0)
        ax.set_ylim(0.84, 1.12)
        if 0 < i < 4:
            ax.set_yticklabels([])
        ax.set_xlabel('')
        cur.append(ax)
    fig.text((0.5 + 2.9 + 0.68) / 2 / W, 0.33 / H, 'Temperature (K)', ha='center', va='center', fontsize=7.5)
    cur[4].set_xlabel('Temperature (K)')
    cur[4].set_ylim(0.84, 1.31)              # lambda-repressor: its measured 282 K point (1.27) is above the others
    cur[4].set_yticks([0.9, 1.0, 1.1, 1.2, 1.3])
    cur[4].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    lo = ax_in(5.42, 0.62, 1.04, 0.93)
    FX.loo_scatter(lo)
    lo.set_xlabel('Measured (% per K)')
    lo.set_ylabel('Predicted, left out\n(% per K)')
    lo.set_title('Left-out slopes', fontsize=7.5, pad=3)
    h = [plt.Line2D([], [], marker='o', ls='', color=FX.INK, ms=2.8, label='measured (smFRET)'),
         plt.Line2D([], [], color=FX.COL['CALVADOS-2'], lw=1.4, label='CALVADOS-2'),
         plt.Line2D([], [], color=FX.COL['Mpipi-Recharged'], lw=1.4, label='Mpipi-Recharged'),
         plt.Line2D([], [], color='0.4', lw=1.4, label='calibrated'),
         plt.Line2D([], [], color='0.4', lw=0.9, ls=(0, (3, 2)), label='plain')]
    fig.legend(handles=h, loc='center', ncol=5, bbox_to_anchor=(0.5, 0.09 / H), handlelength=2.0, columnspacing=1.6)
    FX.add_letters(fig, [(a, 'a', {'row': 0}), (b, 'b', {'row': 0}), (cc, 'c', {'row': 0}),
                         (d, 'd', {'row': 1}), (e, 'e', {'row': 1, 'xaxes': [e]}), (f, 'f', {'row': 1})] +
                   [(ax, 'ghijk'[i], {'row': 2}) for i, ax in enumerate(cur)] + [(lo, 'l', {'row': 2})])
    FX.save(fig, OUT_MERGED)
    return k_c, k_e


INK_ = FX.INK

if __name__ == '__main__':
    c = compute()
    parts, se0, chi_se, g = check(c)
    val, se = response_se(c, g)
    for k in val:
        assert abs(val[k] - parts[k]) < 1e-6, (k, val[k], parts[k])
    print({k: (round(val[k], 3), round(se[k], 3)) for k in val})
    kc, ke = plot_merged(c, val, se, g) if '--old' not in sys.argv else plot(c, val, se, g)
    print('frames', kc, ke, 'Rg', round(float(c['rg'][kc]), 2), round(float(c['rg'][ke]), 2))
