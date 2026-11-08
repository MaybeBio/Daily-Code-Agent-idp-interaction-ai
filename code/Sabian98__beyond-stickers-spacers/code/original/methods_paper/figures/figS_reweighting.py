"""
figS_reweighting.py — visual companion of the charge-model reweighting (Supplementary Section on phosphorylation),
in the spirit of FINCHES Fig. 1: what the weights are, where they come from in one protein, and the two separate uses
(spread -> fluctuation that ranks sites; weighted averages -> reweighted properties), kept apart from the directly
measured rewiring. Example protein P28667 from the 20-protein direct set (cache from figS_reweighting_prep.py).
  ~/.venvs/finches/bin/python figS_reweighting.py   ->  figS_reweighting.pdf/.png
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.stats import spearmanr

import figstyle as fs

UID = 'P28667'
DQ = -2.0
POS, NEG, SITE = '#2166ac', '#b2182b', '#E6A23C'
CLS = {'low': ('#a6a6a6', 'low fluctuation'), 'partner': ('#7570b3', 'matched'), 'annotated': ('#1b6ca8', 'annotated')}
Q = {'K': 1, 'R': 1, 'D': -1, 'E': -1}


def proj(x):
    """Centre one frame and project onto its first two principal axes (nm)."""
    x = x - x.mean(0)
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    return x @ vt[:2].T


def draw_frame(ax, x3, seq, site, note=None, reach=1.5):
    p = proj(x3)
    ax.plot(p[:, 0], p[:, 1], color='0.6', lw=0.7, zorder=1)
    for i, c in enumerate(seq):
        if c in Q:
            ax.scatter(*p[i], s=7, color=POS if Q[c] > 0 else NEG, lw=0, zorder=3)
    d = np.linalg.norm(x3 - x3[site], axis=1)
    for j in np.where(d < reach)[0]:
        if seq[j] in Q and abs(j - site) >= 2:
            ax.plot(*zip(p[site], p[j]), color=POS if Q[seq[j]] > 0 else NEG, lw=0.6, ls=(0, (2, 1.5)), zorder=2)
    ax.scatter(*p[site], s=55, marker='*', color=SITE, edgecolor='k', lw=0.4, zorder=4)
    ax.set_aspect('equal')
    ax.set_axis_off()
    if note:
        ax.text(0.5, -0.04, note, transform=ax.transAxes, ha='center', va='top', fontsize=5.6, color='#333333')


def ftitle(fig, x, y, letter, title):
    fig.text(x, y, letter, fontsize=9, fontweight='bold', ha='left', va='baseline')
    fig.text(x + 0.022, y, fs.cap(title), fontsize=6.5, ha='left', va='baseline', color='#333333')


def box(fig, x, y, w, h, text, fc, fs_=5.8, ha='center'):
    fig.patches.append(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.012', transform=fig.transFigure,
                                      fc=fc, ec='0.55', lw=0.6, zorder=1))
    fig.text(x + (w / 2 if ha == 'center' else 0.012), y + h / 2, text, ha=ha, va='center', fontsize=fs_, color='#222222',
             linespacing=1.35, zorder=2)


def arrow(fig, a, b, color='0.35'):
    fig.patches.append(FancyArrowPatch(a, b, transform=fig.transFigure, arrowstyle='-|>', mutation_scale=7, lw=0.8,
                                       color=color, zorder=3))


def main():
    z = np.load(os.path.join(fs.HERE, 'cache', f'reweighting_{UID}.npz'))
    seq, xyz, phi, beta = str(z['seq']), z['xyz'], z['phi_sites'], float(z['beta'])
    rg = np.asarray(z['rg'], float)
    pos, roles, E = z['pos'], [str(r) for r in z['roles']], z['E']
    n_f = len(xyz)
    y = -beta * DQ * phi                                   # -beta * dH per frame and site
    w = np.exp(y - y.max(0))
    w /= w.sum(0)
    ess = 1 / (w ** 2).sum(0) / n_f
    ia = roles.index('annotated')
    lab = {r: f'{seq[p - 1]}{p}' for r, p in zip(roles, pos)}

    fig = plt.figure(figsize=(7.2, 3.85))

    # ---------------- row 1: the pipeline, drawn with the example's own numbers
    ia_, il = ia, roles.index('low')
    ax = fig.add_axes([0.01, 0.6, 0.14, 0.3])
    ftitle(fig, 0.01, 0.93, 'a', 'unmodified ensemble')
    rng = np.random.default_rng(3)
    for f in rng.choice(n_f, 30, replace=False):
        p = proj(xyz[f])
        ax.plot(p[:, 0], p[:, 1], color='0.55', lw=0.35, alpha=0.35)
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.text(0.5, -0.05, f'{n_f:,} frames', transform=ax.transAxes, ha='center', va='top', fontsize=5.8, color='#333333')

    ax = fig.add_axes([0.185, 0.6, 0.15, 0.3])
    ftitle(fig, 0.185, 0.93, 'b', 'potential at the site')
    f0 = int(np.argsort(w[:, ia])[n_f // 2])
    draw_frame(ax, xyz[f0], seq, pos[ia] - 1)
    ax.text(0.5, -0.05, r'$\phi_s(\mathbf{x}) = \sum\, q_j\,\varepsilon\,f(r_{sj})$', transform=ax.transAxes, ha='center',
            va='top', fontsize=6.2, color='#333333')

    # c: the weight as a function of the potential, with this site's frames as a rug
    ax = fig.add_axes([0.405, 0.62, 0.15, 0.26])
    ftitle(fig, 0.375, 0.93, 'c', 'phosphorylation as a weight')
    ph = phi[:, ia]
    xs = np.linspace(np.percentile(ph, 0.5), np.percentile(ph, 99.5), 200)
    fac = np.exp(-beta * DQ * (xs - ph.mean()))
    ax.plot(xs, fac, color='k', lw=1.0)
    cut = 0.0 if xs[0] < 0 < xs[-1] else ph.mean()            # sign of the potential: acidic or basic surroundings
    ax.fill_between(xs, 0, fac, where=xs <= cut, color=NEG, alpha=0.25, lw=0)
    ax.fill_between(xs, 0, fac, where=xs >= cut, color=POS, alpha=0.18, lw=0)
    ax.plot(ph, np.full_like(ph, 0.03 * fac.max()), '|', color='0.3', ms=4, mew=0.3, alpha=0.35)
    ax.text(0.04, 0.93, r'$w \propto e^{-\beta\,\Delta q\,\phi_s}$, $\Delta q = -2$', transform=ax.transAxes, fontsize=6,
            va='top')
    ax.set_ylim(0, fac.max() * 1.05)
    ax.set_xlabel(r'$\phi_s$ in each frame (ticks)', labelpad=2)
    ax.set_ylabel('weight', labelpad=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.02, 0.12, 'D, E\nnearby', transform=ax.transAxes, fontsize=5.2, color=NEG, va='bottom')
    ax.text(0.98, 0.12, 'K, R\nnearby', transform=ax.transAxes, fontsize=5.2, color=POS, va='bottom', ha='right')

    # branch 1: spread of the weights (two sites, 60 frames each)
    sel = np.sort(np.random.default_rng(5).choice(n_f, 60, replace=False))
    ymax = max(n_f * w[sel, ia].max(), n_f * w[sel, il].max()) * 1.05
    for k, (yy, r) in enumerate(((0.8, 'low'), (0.715, 'annotated'))):
        axw = fig.add_axes([0.64, yy, 0.17, 0.07])
        kk = roles.index(r)
        axw.bar(np.arange(60), n_f * w[sel, kk], width=0.8, color=CLS[r][0], lw=0)
        axw.axhline(1, color='0.4', lw=0.4, ls=(0, (2, 2)))
        axw.set_ylim(0, ymax)
        axw.set_axis_off()
        axw.text(-0.02, 0.5, lab[r], transform=axw.transAxes, ha='right', va='center', fontsize=5.4, color=CLS[r][0])
    fig.text(0.6, 0.885, 'Spread of the weights', fontsize=6.5, color='#333333', ha='left', va='baseline')
    fig.text(0.825, 0.79, 'ESS fraction\n' + r'$\rightarrow$ ranks sites' + '\n(Fig. 5b, g)', fontsize=5.8, va='center',
             color='#333333', linespacing=1.3)

    # branch 2: weighted averages (radius of gyration of this site)
    axr = fig.add_axes([0.64, 0.545, 0.17, 0.085])
    bins_r = np.linspace(np.percentile(rg, 0.5), np.percentile(rg, 99.5), 26)
    axr.hist(rg, bins=bins_r, histtype='stepfilled', color='0.85', density=True, lw=0)
    axr.hist(rg, bins=bins_r, weights=w[:, ia], histtype='step', color=CLS['annotated'][0], density=True, lw=1.0)
    for v, c in ((rg.mean(), '0.45'), (float(w[:, ia] @ rg), CLS['annotated'][0])):
        axr.axvline(v, color=c, lw=0.8)
    axr.set_yticks([])
    axr.spines['left'].set_visible(False)
    axr.tick_params(axis='x', labelsize=5.2, length=2, pad=1)
    axr.set_xlabel('Rg (nm)', fontsize=5.6, labelpad=1)
    fig.text(0.6, 0.66, 'Weighted averages', fontsize=6.5, color='#333333', ha='left', va='baseline')
    fig.text(0.825, 0.59, 'reweighted Rg\n(Fig. 5d)', fontsize=5.8, va='center', color='#333333', linespacing=1.3)
    fig.text(0.825, 0.505, '$E$: direct runs\n(Fig. 5a, e)', fontsize=5.8, va='center', color='0.45', linespacing=1.3)

    arrow(fig, (0.15, 0.76), (0.185, 0.76))
    arrow(fig, (0.34, 0.76), (0.385, 0.76))
    arrow(fig, (0.565, 0.8), (0.615, 0.83))
    arrow(fig, (0.565, 0.7), (0.615, 0.615))

    # ---------------- row 2: the example protein
    bot, h = 0.11, 0.3
    # d: favoured and disfavoured frames for the annotated site
    fmax, fmin = int(np.argmax(w[:, ia])), int(np.argmin(w[:, ia]))
    ax1 = fig.add_axes([0.0, bot + 0.02, 0.12, h])
    ax2 = fig.add_axes([0.12, bot + 0.02, 0.12, h])
    draw_frame(ax1, xyz[fmax], seq, pos[ia] - 1, note=r'$Nw$' + f' = {n_f * w[fmax, ia]:.1f}')
    draw_frame(ax2, xyz[fmin], seq, pos[ia] - 1, note=r'$Nw$' + f' = {n_f * w[fmin, ia]:.2f}')
    ftitle(fig, 0.01, 0.45, 'd', f'{lab["annotated"]}, high and low weight')
    ax1.scatter([], [], s=7, color=POS, label='K, R')
    ax1.scatter([], [], s=7, color=NEG, label='D, E')
    ax1.legend(frameon=False, fontsize=5.2, loc='upper left', bbox_to_anchor=(0.0, 1.02), handletextpad=0.1,
               borderaxespad=0, ncol=2, columnspacing=0.6)

    # e: spread of -beta dH for the three sites
    ax = fig.add_axes([0.31, bot + 0.04, 0.18, h - 0.02])
    bins = np.linspace(-3, 5, 41)
    for r in ('low', 'partner', 'annotated'):
        k = roles.index(r)
        yc = y[:, k] - y[:, k].mean()
        ax.hist(yc, bins=bins, histtype='step', color=CLS[r][0], lw=1.0, density=True,
                label=f'{lab[r]} {CLS[r][1]}')
    ax.set_xlabel(r'$-\beta\Delta H$ per frame (mean removed)')
    ax.set_ylabel('density')
    ax.legend(frameon=False, fontsize=5.0, loc='upper left', handlelength=1.2, borderaxespad=0)
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 60)
    ax.text(0.97, 0.38, 'upweighted\n' + r'$\rightarrow$', transform=ax.transAxes, ha='right', fontsize=5.2, color='0.35')
    ftitle(fig, 0.27, 0.45, 'e', 'spread of the potential')

    # f: sorted weights
    ax = fig.add_axes([0.575, bot + 0.04, 0.17, h - 0.02])
    rank = np.arange(1, n_f + 1) / n_f
    for r in ('low', 'partner', 'annotated'):
        k = roles.index(r)
        ax.plot(rank, np.sort(n_f * w[:, k])[::-1], color=CLS[r][0], lw=1.0, label=f'{lab[r]}  {ess[k]:.2f}')
    ax.axhline(1, color='0.5', lw=0.5, ls=(0, (3, 2)))
    ax.text(0.98, 1.15, 'equal weights', ha='right', va='bottom', fontsize=5.2, color='0.4', transform=ax.get_yaxis_transform())
    ax.set_yscale('log')
    ax.set_xlabel('frames, ranked by weight (fraction)')
    ax.set_ylabel(r'weight $\times\ N$')
    ax.legend(frameon=False, fontsize=5.0, loc='lower left', handlelength=1.2, borderaxespad=0)
    ftitle(fig, 0.53, 0.45, 'f', 'weights concentrate (ESS fraction)')

    # g: all 60 direct sites, fluctuation against directly measured rewiring
    ax = fig.add_axes([0.82, bot + 0.04, 0.165, h - 0.02])
    sd = pd.read_csv(os.path.join(fs.EV, 'phospho_direct', 'sites_direct.csv'))
    ax.scatter(sd.log_ess, sd.E_occ, s=6, color='#d0d0d0', lw=0)
    for r in ('low', 'partner', 'annotated'):
        k = roles.index(r)
        ax.scatter(z['log_ess'][k], E[k], s=18, color=CLS[r][0], edgecolor='k', lw=0.4, zorder=3)
        ax.annotate(lab[r], (z['log_ess'][k], E[k]), xytext=(4, 2), textcoords='offset points', fontsize=5.2,
                    color=CLS[r][0])
    ax.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)))
    ax.text(0.97, 0.95, f'ρ = {spearmanr(sd.log_ess, sd.E_occ).correlation:.2f}', transform=ax.transAxes, ha='right',
            va='top', fontsize=5.8)
    ax.set_xlim(-1.6, 0.05)
    ax.set_xlabel('log ESS fraction\n(unmodified runs)')
    ax.set_ylabel('$E$ (direct runs)')
    ftitle(fig, 0.78, 0.45, 'g', '60 direct sites')
    fs.save(fig, 'figS_reweighting')
    print('ESS', dict(zip(roles, ess.round(3))), 'E', dict(zip(roles, E.round(3))), 'Nw max/min',
          round(n_f * w[fmax, ia], 2), round(n_f * w[fmin, ia], 3))


if __name__ == '__main__':
    main()
