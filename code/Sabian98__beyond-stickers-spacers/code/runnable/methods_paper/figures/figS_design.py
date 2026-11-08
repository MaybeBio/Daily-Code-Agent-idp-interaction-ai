"""
figS_design.py — the prospective design batch in CALVADOS-2 and Mpipi-Recharged. Panel functions are shared with
main-text Fig. 4 (fig4_validation.py, panels c-e); the supplementary figure itself (2026-10-11) shows
a, b  simulated (Mpipi-Recharged) and predicted change of hub occupancy along the chain for two substitution variants
      (each profile standardised, since predicted changes are smaller than simulated ones, and shown as a 5-residue
      running mean; correlations use the raw profiles); substituted residues marked by triangles.
c     Spearman correlation (chain position removed) with the experimental PED ensembles for the two simulations, the
      frozen model and the residue-neutral null, per construct (20 constructs, MD ensembles excluded); blind = not used
      in development (bold names, above the dashed line).
draw_response (per-variant correlation by type and force field, with the stickiness and charge rule and the noise
ceiling) and draw_sites (mean change at the substituted residues) are drawn in Fig. 4.
  ~/.venvs/bender/bin/python figS_design.py  ->  figS_design.pdf / .png
"""
import csv
import os
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy.stats import spearmanr

import figstyle as fs

sys.path.insert(0, fs.EV)
import design_response as DR  # noqa: E402

TYPES = [('DE2NQ', 'D/E to N/Q'), ('K2R', 'K to R'), ('R2K', 'R to K'), ('aro2S', 'F/W/Y to S')]
FFS = [('CALVADOS-2', 'labels_design_calvados'), ('Mpipi-Recharged', 'labels_design_recharged')]
FFC = {ff: fs.FF_COLORS[ff] for ff, _ in FFS}
SIM, MODEL, NULL = FFC['Mpipi-Recharged'], '#E69F00', '#999999'
BLIND = {'PED00501', 'PED00502', 'PED00503', 'PED00154', 'PED00017', 'PED00216'}
NAMES = {'PED_aSyn': 'α-synuclein', 'PED_bSyn': 'β-synuclein', 'PED_abSyn': 'α/β-synuclein', 'PED_drkN': 'drkN',
         'PED_Sic1': 'Sic1', 'PED_pSic1E': 'pSic1, mimetic', 'PED_eIF4G1': 'eIF4G1', 'PED_V2R': 'V2R tail',
         'PED_GHSR': 'GHSR tail', 'PED_B2AR': 'β2AR tail', 'PED_PDEg': 'PDEγ', 'PED_PaaA2': 'PaaA2',
         'PED_p27KID': 'p27 KID', 'PED_FUSNTD': 'FUS NTD', 'PED_V2Rpm': 'V2R tail, mimetic',
         'PED_GHSRpm': 'GHSR tail, mimetic', 'PED_B2ARpm': 'β2AR tail, mimetic', 'PED_MKK4': 'MKK4 N-term',
         'PED_TauF': 'tau fragment', 'PED_CBPID4': 'CBP ID4'}


def load():
    table = {r['id']: r for r in csv.DictReader(open(os.path.join(DR.GAIVI, 'hub_design_table.csv')))}
    pred = {}
    for r in csv.DictReader(open(os.path.join(DR.GAIVI, 'predictions_frozen_2026-10-02.csv'))):
        pred.setdefault(r['id'], []).append(float(r['consensus']))
    pred = {k: np.array(v) for k, v in pred.items()}
    var = {}
    for ff, d in FFS:
        L = DR.read_labels(d)
        resp = {r['id']: r for r in csv.DictReader(open(os.path.join(DR.OUT, f'design_response_{ff}.csv')))}
        var[ff] = []
        for i in sorted(resp):
            par = table[i]['parent']
            pos = [int(x) for x in table[i]['positions'].split()]
            ds, dp = L[i] - L[par], pred[i] - pred[par]
            var[ff].append(dict(id=i, kind=resp[i]['kind'], rho=float(resp[i]['rho_pred_sim']),
                                ceil=float(resp[i]['ceiling']), base=float(resp[i]['rho_baseline_sim']), ds=ds, dp=dp,
                                pos=pos, site=(float(np.mean(dp[pos])), float(np.mean(ds[pos])))))
    exp = {}
    for ff, _ in FFS:
        for line in open(os.path.join(DR.OUT, f'design_vs_experiment_{ff}.txt')):
            m = re.match(r'(PED_\S+)\s+(PED\d+)\s+\d+\s+\S+/(\S+)\s+\S+/(\S+)\s+\S+/(\S+)\s+\S+\s+\d+\s*(\(MD\))?', line)
            if m and not m.group(6):
                d = exp.setdefault(m.group(2), dict(c=m.group(1), e=m.group(2), model=float(m.group(4)),
                                                    null=float(m.group(5))))
                d[ff] = float(m.group(3))
    exp = list(exp.values())
    return var, exp


def z(x):
    return (x - x.mean()) / x.std()


def smooth(x, w=2):
    """Running mean over 2w + 1 residues, for display only (correlations use the raw profiles)."""
    return np.array([x[max(0, i - w):i + w + 1].mean() for i in range(len(x))])


def pick_examples(var):
    """One variant per type near the median correlation of its type, with a high noise ceiling."""
    picks = []
    for kind in ('DE2NQ', 'aro2S'):
        vv = [v for v in var if v['kind'] == kind and v['ceil'] > 0.7]
        med = np.median([v['rho'] for v in var if v['kind'] == kind])
        picks.append(min(vv, key=lambda v: abs(v['rho'] - med)))
    return picks


def draw_example(ax, v, ylabel=True, legend=False, short=False):
    """Simulated (Mpipi-Recharged) and predicted change along the chain for one variant, standardised and smoothed."""
    x = np.arange(1, len(v['ds']) + 1)
    ys, yp = smooth(z(v['ds'])), smooth(z(v['dp']))
    ax.axhline(0, color='0.6', lw=0.5)
    ax.fill_between(x, 0, ys, color=SIM, alpha=0.22, lw=0, label='simulated change (Mpipi-Recharged)')
    ax.plot(x, ys, color=SIM, lw=0.8)
    ax.plot(x, yp, color=MODEL, lw=1.4, label='predicted change')
    lo = min(ys.min(), yp.min()) - 0.45
    ax.scatter(np.array(v['pos']) + 1, np.full(len(v['pos']), lo + 0.15), marker='^', s=14,
               color='0.3', edgecolor='white', linewidth=0.4, zorder=4, label='substituted residue')
    ax.set_ylim(lo, max(ys.max(), yp.max()) + 0.55)
    name = NAMES.get(v['id'].replace('VAR_', 'PED_').rsplit('_', 1)[0], v['id'])
    lab = dict(TYPES)[v['kind']]
    ax.text(0.01, 0.98, f'{name}, {lab}   ρ = {v["rho"]:.2f}', transform=ax.transAxes, va='top', fontsize=6.5)
    ax.set_xlabel('residue')
    ax.set_xlim(1, len(x))
    if ylabel:
        ax.set_ylabel('standardized change' if short else 'change of occupancy\n(standardized, 5-residue mean)')
    if legend:
        ax.legend(frameon=False, loc=legend if isinstance(legend, str) else 'lower right',
                  bbox_to_anchor=(0.0, 0.9) if legend == 'upper left' else (1.0, 0.1), fontsize=5.5, handlelength=1.2,
                  borderaxespad=0.1)


def draw_response(ax, allvar, legend_fontsize=5.8, short=False):
    """Per-variant correlation of predicted and simulated change, by type and force field, model against the rule."""
    rng = np.random.default_rng(0)
    for t, (kind, lab) in enumerate(TYPES):
        for k, (ff, _) in enumerate(FFS):
            vv = [v for v in allvar[ff] if v['kind'] == kind]
            r = np.array([v['rho'] for v in vv])
            x0 = t + (k - 0.5) * 0.4
            ax.bar(x0, np.median([v['ceil'] for v in vv]), width=0.34, color='0.9', zorder=0)
            ax.scatter(x0 + rng.uniform(-0.08, 0.08, len(r)), r, s=8, color=FFC[ff], edgecolor='white', linewidth=0.3,
                       zorder=3)
            ax.plot([x0 - 0.15, x0 + 0.15], [np.median(r)] * 2, color='k', lw=1.1, zorder=4)
            ax.plot([x0 - 0.15, x0 + 0.15], [np.median([v['base'] for v in vv])] * 2, color='k', lw=0.9,
                    ls=(0, (1.5, 1)), zorder=4)
    ax.axhline(0, color='0.6', lw=0.5)
    ax.set_xticks(range(len(TYPES)), [lab for _, lab in TYPES], rotation=30, ha='right')
    ax.set_ylabel('ρ, predicted\nvs simulated' if short else
                  'Spearman ρ, predicted vs\nsimulated change along chain')
    ax.set_ylim(-0.6 if short else -0.55, 1.5 if short else 1.0)
    ax.set_yticks([-0.5, 0, 0.5, 1.0])
    h = [plt.Line2D([], [], marker='o', ls='', ms=3.5, color=FFC[ff], label=ff) for ff, _ in FFS]
    h += [plt.Line2D([], [], color='k', lw=1.1, label='median, model'),
          plt.Line2D([], [], color='k', lw=0.9, ls=(0, (1.5, 1)),
                     label='median, simple rule' if short else 'median, stickiness + charge rule'),
          Patch(color='0.9', label='noise ceiling')]
    if short:
        ax.legend(handles=h, frameon=False, loc='upper center', fontsize=legend_fontsize, handletextpad=0.3,
                  borderaxespad=0.0, ncol=2, columnspacing=0.6, handlelength=1.3)
    else:
        ax.legend(handles=h, frameon=False, loc='lower left', fontsize=legend_fontsize, handletextpad=0.3,
                  borderaxespad=0.1, ncol=1)


def draw_sites(ax, allvar, short=False):
    """Mean predicted against mean simulated change at the substituted residues, both force fields."""
    ax.axhline(0, color='0.75', lw=0.5)
    ax.axvline(0, color='0.75', lw=0.5)
    txt, stats = [], {}
    for ff, _ in FFS:
        s_ = np.array([v['site'] for v in allvar[ff]])
        ax.scatter(s_[:, 0], s_[:, 1], s=11, color=FFC[ff], edgecolor='white', linewidth=0.4, label=ff, zorder=3)
        rr = spearmanr(s_[:, 0], s_[:, 1]).correlation
        agree = np.mean(np.sign(s_[:, 0]) == np.sign(s_[:, 1]))
        stats[ff] = (rr, agree)
        txt.append((ff, f'ρ = {rr:.2f}, sign {agree:.0%}'))
    for k, (ff, t_) in enumerate(txt):
        ax.text(0.04, 0.97 - 0.09 * k, t_, transform=ax.transAxes, va='top', fontsize=6.3, color=FFC[ff])
    ax.set_xlabel('predicted change' if short else 'predicted change at substituted residues')
    ax.set_ylabel('simulated change' if short else 'simulated change at substituted residues')
    return stats


def draw_experiment(ax, exp):
    """Agreement with the experimental ensembles per construct (position removed)."""
    exp = sorted(exp, key=lambda d: (d['e'] in BLIND, d['model']))
    y = np.arange(len(exp))
    for k, d in enumerate(exp):
        xs = [d['CALVADOS-2'], d['Mpipi-Recharged'], d['model']]
        ax.plot([min(xs), max(xs)], [k, k], color='0.8', lw=0.8, zorder=1)
    ax.scatter([d['null'] for d in exp], y, s=10, facecolor='none', edgecolor=NULL, linewidth=0.6, label='null', zorder=2)
    for ff, _ in FFS:
        ax.scatter([d[ff] for d in exp], y, s=11, color=FFC[ff], label=ff, zorder=3)
    ax.scatter([d['model'] for d in exp], y, s=13, color=MODEL, label='sequence model', zorder=4)
    nb = sum(d['e'] not in BLIND for d in exp)
    ax.axhline(nb - 0.5, color='0.5', lw=0.5, ls=(0, (2, 2)))
    ax.axvline(0, color='0.75', lw=0.5)
    ax.set_yticks(y, [NAMES[d['c']] for d in exp], fontsize=5.8)
    for t, d in zip(ax.get_yticklabels(), exp):
        if d['e'] in BLIND:
            t.set_fontweight('bold')
    ax.set_ylim(-0.7, len(exp) - 0.3)
    ax.set_xlim(-0.45, 1.4)
    ax.set_xticks([0, 0.25, 0.5, 0.75], ['0', '0.25', '0.5', '0.75'])
    ax.set_xlabel('Spearman ρ with experimental\nensemble (position removed)')
    ax.legend(frameon=False, loc='lower right', fontsize=6, handletextpad=0.1, borderaxespad=0.0, markerscale=0.9)


def main():
    """Supplementary figure (2026-10-11): the response along the chain and the change at the substituted residues moved
    to main-text Fig. 4 (fig4_validation.py), so this figure keeps the two example profiles and the per-construct
    agreement with experiment."""
    allvar, exp = load()
    assert all(len(v) == 49 for v in allvar.values()) and len(exp) == 20, (len(exp),)
    var = allvar['Mpipi-Recharged']
    picks = pick_examples(var)
    fig = plt.figure(figsize=(7.2, 4.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.7, 1], hspace=0.55, wspace=0.55, left=0.1, right=0.98, top=0.93,
                          bottom=0.12)
    for k, v in enumerate(picks):
        ax = fig.add_subplot(gs[k, 0])
        draw_example(ax, v, ylabel=True, legend=(k == 1))
        fs.panel(ax, 'ab'[k], x=-0.14)
    ax = fig.add_subplot(gs[:, 1])
    draw_experiment(ax, exp)
    fs.panel(ax, 'c', x=-0.55, y=1.01)
    fs.save(fig, 'figS_design')
    stats = {}
    for ff, _ in FFS:
        s_ = np.array([v['site'] for v in allvar[ff]])
        stats[ff] = (spearmanr(s_[:, 0], s_[:, 1]).correlation, np.mean(np.sign(s_[:, 0]) == np.sign(s_[:, 1])))
    with open(os.path.join(fs.HERE, 'figS_design_values.txt'), 'w') as fh:
        fh.write(f'A examples (Mpipi-Recharged): {[(v["id"], round(v["rho"], 3), round(v["ceil"], 3)) for v in picks]}\n')
        for ff, (rr, agree) in stats.items():
            fh.write(f'C {ff}: site rho {rr:.3f}, sign agreement {agree:.3f}\n')
        fh.write('D medians CALVADOS-2/Mpipi-Recharged/model/null: ' + ' / '.join(
            f'{np.median([d[k] for d in exp]):.3f}' for k in ('CALVADOS-2', 'Mpipi-Recharged', 'model', 'null')) + '\n')


if __name__ == '__main__':
    main()
