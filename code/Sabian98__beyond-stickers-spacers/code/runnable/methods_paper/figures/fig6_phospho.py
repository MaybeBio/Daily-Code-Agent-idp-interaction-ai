"""
fig6_phospho.py — Figure 6 (2026-10-13 reframe): the contact network predicts how phosphorylation rewires disordered
chains. All values come from the hashed results in redo_evidence (no new statistics).
  a  Direct rewiring E(s) of hub occupancy for the low-fluctuation, matched unannotated and annotated site of each of the
     20 proteins (240 CALVADOS-2 runs; phospho_direct/sites_direct.csv), as half violins with points and medians.
  b  E(s) against the potential fluctuation of the site in the unmodified ensemble (log effective sample fraction).
  c  Within-protein Spearman correlation of direct rewiring (E / replicate noise, centred) with the ensemble fluctuation,
     the fitted sequence charge model, local charge counts and the SCD change (post hoc); phospho_charge/direct_sites.csv.
  d  Reweighted (BENDER-MT charge model, unmodified ensemble) against direct Rg change, 60 sites.
  e  NMR case series: partial Spearman (predicted change vs amide CSP given distance to the nearest site) per protein for
     contact rewiring E_i (pre-registered), |degree change| and |hub occupancy change| (phospho_nmr/per_protein.csv), with
     the Mpipi-Recharged contact values as open diamonds (phospho_nmr_recharged/per_protein.csv, Amendment 1); since
     2026-10-15 also the four pairs of Amendment 2 (phospho_nmr_amend2/results_{cal,rech}/per_protein.csv).
  f  HMGA1a S64C with S99, S102, S103 phosphorylated: change in contact probability (mean of three phosphorylated minus
     mean of three unmodified runs; case_hmga1a/hmga1a_maps.npz from hmga1a_render_prep.py, descriptive) above the
     diagonal and, below it, the difference between the two most different unmodified replicates (noise), with the
     measured amide CSP and the predicted contact rewiring (each scaled to its maximum) along the top.
  ~/.venvs/finches/bin/python fig6_phospho.py  ->  fig6_phospho.pdf / .png
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import figstyle as fs

EV = fs.EV
CLS = {'low': ('#a6a6a6', 'low fluctuation'), 'partner': ('#7570b3', 'matched'), 'annotated': ('#1b6ca8', 'annotated')}
Q = {'K': 1, 'R': 1, 'D': -1, 'E': -1}


def centre(v, g):
    v = pd.Series(np.asarray(v, float))
    return (v - v.groupby(g).transform('mean')).to_numpy()


def dscd(seq, pos):
    i, n = pos - 1, len(seq)
    return sum(-2.0 * Q.get(c, 0) * abs(i - j) ** 0.5 for j, c in enumerate(seq) if j != i) / n


def main():
    sd = pd.read_csv(os.path.join(EV, 'phospho_direct', 'sites_direct.csv'))
    dc = pd.read_csv(os.path.join(EV, 'phospho_charge', 'direct_sites.csv'))
    npr = pd.read_csv(os.path.join(EV, 'phospho_nmr', 'per_protein.csv')).set_index('protein')
    nres = pd.read_csv(os.path.join(EV, 'phospho_nmr', 'per_residue.csv'))

    fig = plt.figure(figsize=(7.2, 2.78))
    outer = fig.add_gridspec(1, 2, width_ratios=[3, 1.28], wspace=0.22, left=0.062, right=0.935, top=0.93, bottom=0.12)
    gs = outer[0].subgridspec(2, 3, width_ratios=[1.24, 0.92, 1.06], hspace=0.85, wspace=0.7)
    gf = outer[1].subgridspec(3, 2, height_ratios=[0.09, 0.27, 1], width_ratios=[1, 0.055], hspace=0.05, wspace=0.07)

    # a: E by site class, paired by protein
    ax = fig.add_subplot(gs[0, 0])
    order = ['low', 'partner', 'annotated']
    piv = sd.pivot(index='UniProt_ID', columns='role', values='E_occ')[order]
    ax.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=0)
    fs.half_violins(ax, [piv[c] for c in order], [CLS[c][0] for c in order], ['low', 'matched', 'annotated'])
    ax.tick_params(axis='x', labelsize=5.2)
    ax.set_ylabel('rewiring, $E$')
    fs.panel(ax, 'a', x=-0.38, y=1.08, title='240 direct runs, 20 proteins')

    # b: E vs log ESS
    ax = fig.add_subplot(gs[0, 1])
    for c in order:
        s = sd[sd.role == c]
        ax.scatter(s.log_ess, s.E_occ, s=10, color=CLS[c][0], lw=0, label=CLS[c][1])
    r = spearmanr(sd.log_ess, sd.E_occ).correlation
    ax.text(0.97, 0.93, f'ρ = {r:.2f}', transform=ax.transAxes, ha='right', fontsize=6.5)
    ax.axhline(0, color='0.5', lw=0.5, ls=(0, (3, 2)))
    ax.set_xlabel('log effective sample fraction\n(less fluctuation →)')
    ax.set_ylabel('$E$')
    fs.panel(ax, 'b', x=-0.27, y=1.08, title='fluctuation predicts rewiring')

    # c: ensemble vs sequence predictors (within-protein Spearman with y = E / R(W, W))
    ax = fig.add_subplot(gs[0, 2])
    g = dc.UniProt_ID.to_numpy()
    y = centre(dc.y_occ, g)
    seqs = {u: str(np.load(os.path.join(EV, 'phospho_network', 'profiles', f'{u}.npz'))['seq']) for u in dc.UniProt_ID.unique()}
    scd = np.array([abs(dscd(seqs[u], p)) for u, p in zip(dc.UniProt_ID, dc.pos)])
    import sys
    sys.path.insert(0, EV)
    import phospho_charge_posthoc_kappa as KP                      # post hoc kappa control (hashed 2026-10-14)
    dk = np.array([abs(KP.dkappa(seqs[u], p, -2.0)) for u, p in zip(dc.UniProt_ID, dc.pos)])
    vals = [spearmanr(centre(-dc.log_ess, g), y).correlation, spearmanr(centre(-dc.seqpred, g), y).correlation,
            spearmanr(centre(-dc.litpred, g), y).correlation, spearmanr(centre(scd, g), y).correlation,
            spearmanr(centre(dk, g), y).correlation]
    labs = ['ensemble', 'fitted model', 'charges ±5', 'SCD change', 'κ change']
    cols = ['#1b6ca8', '#8fb4d6', '#c9c9c9', '#e0e0e0', '#e0e0e0']
    bars = ax.bar(range(5), vals, color=cols, width=0.7)
    for b_ in bars[3:]:
        b_.set_hatch('////')
        b_.set_edgecolor('0.55')
    ax.axhline(0, color='0.5', lw=0.5)
    ax.set_xticks(range(5), labs, fontsize=5.6, rotation=35, ha='right', rotation_mode='anchor')
    ax.set_ylabel('within-protein\nρ with $E$', linespacing=1.1)
    ax.set_ylim(-0.24, 0.75)
    for k, v in enumerate(vals):
        ax.text(k, v + (0.02 if v > 0 else -0.015), f'{v:.2f}', ha='center', va='bottom' if v > 0 else 'top', fontsize=5.2)
    fs.panel(ax, 'c', x=-0.3, y=1.08, title='ensemble or sequence rules')

    # d: reweighted vs direct Rg change
    ax = fig.add_subplot(gs[1, 0])
    for c in order:
        s = sd[sd.role == c]
        ax.scatter(s.dRg_rw_pct, s.dRg_direct_pct, s=10, color=CLS[c][0], lw=0)
    lim = [min(sd.dRg_rw_pct.min(), sd.dRg_direct_pct.min()) - 0.5, max(sd.dRg_rw_pct.max(), sd.dRg_direct_pct.max()) + 0.5]
    ax.plot(lim, lim, color='0.6', lw=0.6, ls=(0, (3, 2)))
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    r = spearmanr(sd.dRg_rw_pct, sd.dRg_direct_pct).correlation
    ax.text(0.04, 0.9, f'ρ = {r:.2f}', transform=ax.transAxes, fontsize=6.5)
    ax.set_xlabel('reweighted ΔRg (%)')
    ax.set_ylabel('direct ΔRg (%)')
    fs.panel(ax, 'd', x=-0.4, y=1.08, title='one ensemble suffices')

    # e: NMR per protein, the six of the original plan and the four of Amendment 2 (pTyr, alpha-synuclein), both hashed
    ax = fig.add_subplot(gs[1, 1:])
    A2 = os.path.join(EV, 'phospho_nmr_amend2')
    cal = pd.concat([npr, pd.read_csv(os.path.join(A2, 'results_cal', 'per_protein.csv')).set_index('protein')])
    rch = pd.concat([pd.read_csv(os.path.join(EV, 'phospho_nmr_recharged', 'per_protein.csv')).set_index('protein'),
                     pd.read_csv(os.path.join(A2, 'results_rech', 'per_protein.csv')).set_index('protein')])
    prots = ['PAGE4', 'p53TAD', 'Sic1', 'NHE1tail', 'HMGA1a', 'H1-0CTD', 'aSynS129', 'aSynY125', 'CD79b', 'Tir']
    names = ['PAGE4', 'p53', 'Sic1', 'NHE1', 'HMGA1', 'H1.0', 'αSyn\npS129', 'αSyn\npY125', 'CD79b\npY', 'Tir\npY']
    x = np.arange(len(prots)) + np.r_[np.zeros(6), np.full(4, 0.5)]
    w = 0.26
    for k, (col, c, lab) in enumerate((('rho_p', '#1b6ca8', 'contacts, CALVADOS-2'), ('rho_p_ddeg', '#2ca25f', 'degree, CALVADOS-2'),
                                        ('rho_p_docc', '#bdbdbd', 'hub occupancy, CALVADOS-2'))):
        ax.bar(x + (k - 1) * w, cal.loc[prots, col], w, color=c, label=lab)
    ax.scatter(x - w, rch.loc[prots, 'rho_p'], s=12, marker='D', facecolor='white', edgecolor='#2F7DB5',
               lw=0.8, zorder=4, label='contacts, Mpipi-Recharged')
    ax.axhline(0, color='0.5', lw=0.5)
    ax.axvline(5.75, color='0.6', lw=0.5, ls=(0, (3, 2)))
    ax.text(5.85, -0.33, 'added later,\nincl. pTyr', fontsize=4.8, color='0.4', va='bottom')
    ax.set_xticks(x, names, fontsize=4.6)
    ax.set_xlim(x[0] - 0.55, x[-1] + 0.55)
    ax.set_ylabel('partial ρ with NMR')
    ax.legend(frameon=False, fontsize=5.0, loc='upper center', handlelength=1, borderaxespad=0.0, ncol=2, columnspacing=1.0,
              bbox_to_anchor=(0.5, 1.05), handletextpad=0.3)
    ax.set_ylim(-0.35, 0.9)
    ax.set_yticks([-0.25, 0, 0.25, 0.5])
    fs.panel(ax, 'e', x=-0.12, y=1.08, title='ten proteins, NMR conditions')

    # f: HMGA1a, change in contact probability (map, descriptive) with measured CSP and predicted E_i along the top
    mp = np.load(os.path.join(fs.HERE, 'case_hmga1a', 'hmga1a_maps.npz'))
    dp, n = mp['dp'], len(mp['E'])
    hooks, sites = [(21, 31), (53, 63), (78, 89)], [99, 102, 103]       # UniProt P17096 AT-hooks
    ax = fig.add_subplot(gf[2, 0])
    tr = fig.add_subplot(gf[1, 0], sharex=ax)
    cax = fig.add_subplot(gf[2, 1])
    W = mp['pW_reps']
    pairs_ = [(0, 1), (0, 2), (1, 2)]
    a_, b_ = max(pairs_, key=lambda q: np.abs(W[q[0]] - W[q[1]]).max())   # noisiest pair of unmodified replicates
    upper = np.triu(np.ones((n, n), bool), 1).T                             # y > x: phospho change
    half = np.where(upper, dp, W[a_] - W[b_])
    im = ax.imshow(half, origin='lower', cmap='PRGn', vmin=-0.12, vmax=0.12, extent=(0.5, n + 0.5, 0.5, n + 0.5),
                   interpolation='nearest', aspect='auto')
    ax.text(0.03, 0.97, 'phosphorylated\n− unmodified', transform=ax.transAxes, fontsize=5.3, va='top', color='#333333')
    ax.text(0.97, 0.03, 'unmodified replicate\n− replicate (noise)', transform=ax.transAxes, fontsize=5.3, ha='right',
            va='bottom', color='#333333')
    ax.plot([0.5, n + 0.5], [0.5, n + 0.5], color='0.75', lw=0.5)
    for a_, b_ in hooks:
        for span in (ax.axvspan, ax.axhspan):
            span(a_ - 0.5, b_ + 0.5, color='0.55', alpha=0.08, lw=0)
    for s_ in sites:
        ax.plot([s_], [0.5], marker='^', ms=2.6, color='#b2182b', clip_on=False, zorder=5)
    ax.annotate('tail gains contacts\nwith AT-hooks 2 and 3', xy=(72, 100), xytext=(22, 70), fontsize=5.5, color='#333333',
                arrowprops=dict(arrowstyle='-', lw=0.5, color='#555555'))
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(0.5, n + 0.5)
    ax.set_xticks([1, 25, 50, 75, 100])
    ax.set_yticks([1, 25, 50, 75, 100])
    ax.set_xlabel('HMGA1a residue')
    ax.set_ylabel('residue', labelpad=1)
    cb = fig.colorbar(im, cax=cax, ticks=[-0.1, 0, 0.1])
    cb.set_label('change in contact probability', fontsize=5.8, labelpad=2)
    cb.ax.tick_params(labelsize=5.5, length=1.5)
    cb.outline.set_linewidth(0.4)
    d = nres[nres.protein == 'HMGA1a'].sort_values('pos')
    for a_, b_ in hooks:
        tr.axvspan(a_ - 0.5, b_ + 0.5, color='#eeeeee', lw=0, zorder=0)
        tr.text((a_ + b_) / 2, 1.02, f'AT{hooks.index((a_, b_)) + 1}', transform=tr.get_xaxis_transform(), ha='center',
                va='bottom', fontsize=5.2, color='#666666')
    tr.bar(d.pos, d.CSP / d.CSP.max(), width=0.9, color='#d95f02', lw=0)
    pos, e = d.pos.to_numpy(), d.E.to_numpy().astype(float)
    brk = np.where(np.diff(pos) > 6)[0] + 1                   # do not draw across excluded stretches
    tr.plot(np.insert(pos.astype(float), brk, np.nan), np.insert(e / np.nanmax(e), brk, np.nan), color='#1b6ca8', lw=0.9)
    for s_ in sites:
        tr.axvline(s_, color='#b2182b', lw=0.6)
    tr.set_yticks([])
    tr.spines['left'].set_visible(False)
    tr.tick_params(labelbottom=False, bottom=False)
    tr.text(0.01, 1.0, 'measured', color='#d95f02', transform=tr.transAxes, fontsize=5.5, va='top')
    tr.text(0.01, 0.66, 'predicted', color='#1b6ca8', transform=tr.transAxes, fontsize=5.5, va='top')
    fs.panel(tr, 'f', x=-0.3, y=1.36, title='HMGA1a, three phosphoserines')
    fs.save(fig, 'fig6_phospho')
    print('c values', [round(v, 3) for v in vals])


if __name__ == '__main__':
    main()
