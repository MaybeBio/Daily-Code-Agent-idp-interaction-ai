"""
figS_extra.py — supplementary figures for sections that had none (2026-10-11), from existing outputs only.
  figS_theory   (app:theory)    a  observed time-averaged betweenness of one protein under |i-j| >= 2 and >= 3 with the
                                   closed form j(m-1-j) of the two-path limit; b  per-protein variance explained by the
                                   closed form under both rules (1,198 CALVADOS-2 proteins); c  Spearman-Brown curve of
                                   reliability against snapshots with the observed split-half medians
  figS_ceiling  (app:predictor) per-protein agreement between CALVADOS-2 and Mpipi-Recharged hub occupancy (partial on end
                                   distance) against chain length and Flory exponent, with running medians
  figS_pedv2    (app:pedv2)     partial correlation of the frozen model and the null for every ensemble of the second blind
                                   set (set A) and of the overlapping-construct set (set B)
  figS_llps     (app:llps)      a  wild-type unevenness percentile among local shuffles, canonical drivers against random
                                   IDRome regions; b  PhaSePro drivers against matched controls; c  hnRNPA1 LCD variants,
                                   sticker hub share against the residual of log10 csat; d  disease mutations, percentile of
                                   the local response
  ~/.venvs/bender/bin/python figS_extra.py
"""
import csv
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

import figstyle as fs

EV = fs.EV
L = os.path.join(EV, 'llps_csat')
C3, C2 = fs.RULE_COLORS['sep3'], fs.RULE_COLORS['sep2']
MODEL, NULL = '#E69F00', '#8c8c8c'


def closed_form(n):
    out = np.zeros(n)
    for k in range(n):
        p = k % 2
        m = (n - p + 1) // 2
        j = k // 2
        out[k] = j * (m - 1 - j)
    return out


def load_btw(d):
    by = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, d, 'per_residue.csv'))):
        by[r['UniProt_ID']].append((int(r['pos']), float(r['mean_btw'])))
    return {u: np.array([b for _, b in sorted(v)]) for u, v in by.items()}


def theory():
    d2, d3 = load_btw('hub_frames_sep2'), load_btw('hub_frames')
    ids = sorted(set(d2) & set(d3))
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw={'width_ratios': [1.4, 1, 1]})
    u = 'Q9DUC6' if 'Q9DUC6' in ids else ids[0]
    for d, col, lab in ((d2, C2, '|i−j| ≥ 2'), (d3, C3, '|i−j| ≥ 3')):
        b = d[u]
        ax[0].plot(np.arange(1, len(b) + 1), b / b.max(), color=col, lw=1.0, label=lab)
    c = closed_form(len(d2[u]))
    ax[0].plot(np.arange(1, len(c) + 1), c / c.max(), color='k', lw=0.8, ls=(0, (3, 2)), label='closed form')
    ax[0].set_xlabel('residue')
    ax[0].set_ylabel('betweenness (scaled)')
    ax[0].legend(frameon=False, fontsize=5.8, loc='lower center', ncol=3, handlelength=1.2, columnspacing=0.8)
    ax[0].set_ylim(-0.25, 1.08)
    fs.panel(ax[0], 'a', x=-0.17, title=f'{u}, {len(c)} residues')
    bins = np.linspace(0, 1, 26)
    for d, col, lab in ((d2, C2, '|i−j| ≥ 2'), (d3, C3, '|i−j| ≥ 3')):
        r2 = []
        for v in ids:
            b = d[v]
            if b.std() > 0:
                r2.append(np.corrcoef(b, closed_form(len(b)))[0, 1] ** 2)
        ax[1].hist(r2, bins=bins, color=col, alpha=0.75, label=f'{lab} (median {np.median(r2):.2f})')
    ax[1].set_xlabel('variance explained by closed form')
    ax[1].set_ylabel('proteins')
    ax[1].legend(frameon=False, fontsize=5.8, loc='upper center')
    fs.panel(ax[1], 'b', x=-0.3, title=f'{len(ids):,} CALVADOS-2 proteins')
    txt = open(os.path.join(EV, 'frames_reliability_summary.txt')).read()
    rows = [tuple(map(float, m)) for m in re.findall(r'^\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$', txt, re.M)]
    fr, sh, _ = zip(*rows)
    r1 = float(re.search(r'r1 = ([\d.]+)', open(os.path.join(EV, 'pitfall_theory_summary.txt')).read()).group(1))
    n = np.logspace(0, np.log10(3000), 200)
    ax[2].plot(n, n * r1 / (1 + (n - 1) * r1), color='k', lw=1.0, label=f'Spearman–Brown\n(r$_1$ = {r1:.4f})')
    ax[2].scatter(fr, sh, s=14, color=C3, zorder=3, label='observed (split-half)')
    for t in (0.8, 0.9):
        nt = t * (1 - r1) / (r1 * (1 - t))
        ax[2].plot([nt, nt], [0, t], color='0.6', lw=0.5, ls=(0, (2, 2)))
        ax[2].text(nt, 0.02, f'{nt:.0f}', fontsize=5.5, ha='center', color='0.4')
    ax[2].set_xscale('log')
    ax[2].set_xlim(10, 3000)
    ax[2].set_ylim(0, 1)
    ax[2].set_xlabel('snapshots')
    ax[2].set_ylabel('reliability')
    ax[2].legend(frameon=False, fontsize=5.5, loc='upper left')
    fs.panel(ax[2], 'c', x=-0.3, title='sampling')
    fig.subplots_adjust(left=0.07, right=0.99, top=0.88, bottom=0.2, wspace=0.42)
    fs.save(fig, 'figS_theory')


def ceiling():
    rows = [r for r in csv.DictReader(open(os.path.join(EV, 'agreement_ceiling_per_protein.csv')))
            if r['partial_cal_rech'] not in ('', 'nan')]
    N = np.array([float(r['N']) for r in rows])
    nu = np.array([float(r['nu']) for r in rows])
    p = np.array([float(r['partial_cal_rech']) for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.3), sharey=True)
    for a, x, lab, log in ((ax[0], N, 'chain length (residues)', True), (ax[1], nu, 'Flory exponent ν (CALVADOS-2)', False)):
        a.hexbin(x, p, gridsize=45, cmap='Blues', mincnt=1, xscale='log' if log else 'linear', linewidths=0.1)
        q = np.quantile(x, np.linspace(0, 1, 13))
        mids, meds = [], []
        for lo, hi in zip(q[:-1], q[1:]):
            m = (x >= lo) & (x <= hi)
            mids.append(np.median(x[m]))
            meds.append(np.median(p[m]))
        a.plot(mids, meds, '-o', color='#d6006f', ms=3, lw=1.3, label='running median')
        a.axhline(np.median(p), color='0.3', lw=0.6, ls=(0, (3, 2)), label=f'overall median {np.median(p):.2f}')
        a.set_xlabel(lab)
        a.axhline(0, color='0.7', lw=0.5)
        if not log:
            a.set_xlim(*np.quantile(x, [0.002, 0.998]))
    ax[0].set_ylabel('CALVADOS-2 vs Mpipi-Recharged\n(partial ρ)')
    ax[0].legend(frameon=False, fontsize=5.8, loc='lower right')
    fs.panel(ax[0], 'a', x=-0.2, title=f'{len(p):,} sequences')
    fs.panel(ax[1], 'b', x=-0.1)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.88, bottom=0.2, wspace=0.1)
    fs.save(fig, 'figS_ceiling')


def pedv2():
    rows = list(csv.DictReader(open(os.path.join(EV, 'ped_v2_scores.csv'))))
    rows = [r for r in rows if r['entry'] != 'PED07323']                  # IAPP, excluded by the disorder criterion
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [14, 5]}, sharey=True)
    for ax, st, title in ((axes[0], 'A', 'second test set (14 proteins)'), (axes[1], 'B', 'proteins already tested')):
        R = sorted([r for r in rows if r['set'] == st], key=lambda r: -float(r['partial_model']))
        x = np.arange(len(R))
        for k, r in enumerate(R):
            m, n = float(r['partial_model']), float(r['partial_null'])
            ax.plot([k, k], [n, m], color='#2ca25f' if m > n else '#d6006f', lw=1.2, zorder=1)
        ax.scatter(x, [float(r['partial_null']) for r in R], s=14, facecolor='white', edgecolor=NULL, lw=0.9,
                   zorder=2, label='residue-neutral null')
        ax.scatter(x, [float(r['partial_model']) for r in R], s=16, color=MODEL, zorder=3, label='sequence model')
        ax.axhline(0, color='0.7', lw=0.5)
        ax.set_xticks(x, [r['protein'].replace(' [ortholog]', '*') for r in R], rotation=55, ha='right', fontsize=5.5)
        ax.set_xlim(-0.6, len(R) - 0.4)
        fs.panel(ax, 'ab'[st == 'B'], x=-0.08 if st == 'A' else -0.25, title=title)
    axes[0].set_ylabel('partial correlation\nwith ensemble')
    axes[0].legend(frameon=False, fontsize=5.8, loc='lower left')
    fig.subplots_adjust(left=0.09, right=0.99, top=0.9, bottom=0.36, wspace=0.08)
    fs.save(fig, 'figS_pedv2')


def llps():
    txt = open(os.path.join(L, 'llps_profiles_local.txt')).read()
    drivers = [float(m) for m in re.findall(r'^\s{4}\S.*?\s+(\d+)\s*$', txt, re.M)]
    rnd = [float(r['wt_pct_among_shuffles']) for r in csv.DictReader(open(os.path.join(L, 'llps_profiles_control.csv')))]
    unev = list(csv.DictReader(open(os.path.join(L, 'llps_unevenness.csv'))))
    feat = list(csv.DictReader(open(os.path.join(L, 'llps_features.csv'))))
    muts = list(csv.DictReader(open(os.path.join(L, 'llps_mutations.csv'))))
    fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.2), gridspec_kw={'width_ratios': [1.1, 1, 1, 0.9]})
    bins = np.linspace(0, 100, 11)
    ax[0].hist(rnd, bins=bins, color=NULL, alpha=0.6, density=True, label=f'random IDRs (n = {len(rnd)})')
    ax[0].hist(drivers, bins=bins, color='#d6006f', alpha=0.75, density=True, label=f'canonical drivers (n = {len(drivers)})')
    ax[0].set_xlabel('wild-type percentile\namong local shuffles')
    ax[0].set_ylabel('density')
    ax[0].legend(frameon=False, fontsize=5.2, loc='upper left')
    fs.panel(ax[0], 'a', x=-0.3, title='exploratory')
    d = np.array([float(r['pct_local']) for r in unev])
    c = np.array([float(r['ctrl_pct_local_mean']) for r in unev])
    ax[1].scatter(c, d, s=12, color='#0072b2', edgecolor='white', lw=0.3)
    ax[1].plot([0, 100], [0, 100], color='0.6', lw=0.6, ls=(0, (3, 2)))
    ax[1].set_xlabel('matched controls (mean)')
    ax[1].set_ylabel('PhaSePro driver')
    ax[1].set_xlim(0, 100)
    ax[1].set_ylim(0, 102)
    fs.panel(ax[1], 'b', x=-0.32, title=f'{len(d)} drivers, planned test')
    rows = [r for r in feat if r['csat_M_4C'] not in ('', 'nan')]
    y = np.log10([float(r['csat_M_4C']) for r in rows])
    X = np.column_stack([np.ones(len(rows)), [float(r['n_FWYR']) for r in rows], [float(r['NCPR']) for r in rows]])
    res = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    ax[2].scatter([float(r['F1_sticker_hub_share']) for r in rows], res, s=12, color='#00a651', edgecolor='white',
                  lw=0.3)
    ax[2].axhline(0, color='0.7', lw=0.5)
    ax[2].set_xlabel('sticker share of\npredicted hub occupancy')
    ax[2].set_ylabel('log$_{10}$ c$_{sat}$ residual')
    fs.panel(ax[2], 'c', x=-0.36, title=f'hnRNPA1 LCD ({len(rows)} variants)')
    pa = [float(r['pct_abs']) for r in muts if r['pct_abs'] not in ('', 'nan')]
    ax[3].hist(pa, bins=bins, color='#ff7f00', alpha=0.85)
    ax[3].axvline(50, color='0.4', lw=0.6, ls=(0, (3, 2)))
    ax[3].set_xlabel('percentile of\nlocal response')
    ax[3].set_ylabel('mutations')
    fs.panel(ax[3], 'd', x=-0.36, title=f'{len(pa)} disease mutations')
    fig.subplots_adjust(left=0.07, right=0.99, top=0.86, bottom=0.27, wspace=0.62)
    fs.save(fig, 'figS_llps')
    return len(drivers)


if __name__ == '__main__':
    theory()
    ceiling()
    pedv2()
    print('llps drivers parsed:', llps())
