"""
figS_phospho_recharged.py — Supplementary figure: the phosphorylation result in Mpipi-Recharged (phospho_direct
Amendments 2 and 3). Values from the hashed results in redo_evidence; no new statistics.
  a  Recharged rewiring E (hub occupancy) by site class, 20 proteins, 240 runs
  b  Recharged E against the Recharged log ESS (full energy change, Amendment 3, M1)
  c  Recharged log ESS by site class (M2)
  d  Recharged against CALVADOS-2 log ESS, 60 sites
  e  Recharged against CALVADOS-2 rewiring (degree), 60 sites
  f  Recharged reweighted against direct Rg change (M3)
  ~/.venvs/finches/bin/python figS_phospho_recharged.py  ->  figS_phospho_recharged.pdf / .png
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import figstyle as fs

EV = fs.EV
CLS = {'low': ('#a6a6a6', 'low fluctuation'), 'partner': ('#7570b3', 'matched'), 'annotated': ('#1b6ca8', 'annotated')}
ORDER = ['low', 'partner', 'annotated']


def paired(ax, piv, ylabel):
    fs.half_violins(ax, [piv[c] for c in ORDER], [CLS[c][0] for c in ORDER], [CLS[c][1] for c in ORDER])
    ax.set_ylabel(ylabel)


def scatter(ax, x, y, roles, xlabel, ylabel, diag=False):
    for c in ORDER:
        m = roles == c
        ax.scatter(x[m], y[m], s=10, color=CLS[c][0], lw=0)
    if diag:
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        pad = 0.05 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color='0.6', lw=0.6, ls=(0, (3, 2)))
    r = spearmanr(x, y).correlation
    ax.text(0.04, 0.9, f'ρ = {r:.2f}', transform=ax.transAxes, fontsize=6.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def main():
    s = pd.read_csv(os.path.join(EV, 'phospho_direct_recharged', 'sites_recharged_final.csv'))
    e = pd.read_csv(os.path.join(EV, 'phospho_direct_recharged', 'ess_recharged.csv'))
    m = s.merge(e, on=['UniProt_ID', 'role'])
    roles = m.role.to_numpy()
    fig = plt.figure(figsize=(7.2, 4.1))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.45, left=0.075, right=0.985, top=0.94, bottom=0.11)
    ax = fig.add_subplot(gs[0, 0])
    paired(ax, m.pivot(index='UniProt_ID', columns='role', values='E_occ')[ORDER], 'rewiring beyond noise, $E$')
    ax.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=0)
    fs.panel(ax, 'a', x=-0.3, y=1.08, title='Mpipi-Recharged, 240 direct runs')
    ax = fig.add_subplot(gs[0, 1])
    scatter(ax, m.log_ess_rech.to_numpy(), m.E_occ.to_numpy(), roles, 'Recharged log effective sample fraction', '$E$')
    ax.axhline(0, color='0.5', lw=0.5, ls=(0, (3, 2)))
    fs.panel(ax, 'b', x=-0.27, y=1.08, title='rewiring follows Recharged fluctuation')
    ax = fig.add_subplot(gs[0, 2])
    paired(ax, m.pivot(index='UniProt_ID', columns='role', values='log_ess_rech')[ORDER], 'Recharged log ESS')
    fs.panel(ax, 'c', x=-0.3, y=1.08, title='classes separate in Recharged')
    ax = fig.add_subplot(gs[1, 0])
    scatter(ax, m.log_ess_cal.to_numpy(), m.log_ess_rech.to_numpy(), roles, 'CALVADOS-2 log ESS', 'Recharged log ESS', diag=True)
    fs.panel(ax, 'd', x=-0.3, y=1.08, title='fluctuation in the two models')
    ax = fig.add_subplot(gs[1, 1])
    scatter(ax, m.E_deg_cal.to_numpy(), m.E_deg.to_numpy(), roles, 'CALVADOS-2 $E$ (degree)', 'Recharged $E$ (degree)')
    fs.panel(ax, 'e', x=-0.27, y=1.08, title='same sites rewire in both models')
    ax = fig.add_subplot(gs[1, 2])
    scatter(ax, m.dRg_rw_rech_pct.to_numpy(), m.dRg_direct_rech_pct.to_numpy(), roles, 'reweighted Rg change (%)',
            'direct Rg change (%)', diag=True)
    fs.panel(ax, 'f', x=-0.3, y=1.08, title='reweighting in Recharged')
    fs.save(fig, 'figS_phospho_recharged')


if __name__ == '__main__':
    main()
