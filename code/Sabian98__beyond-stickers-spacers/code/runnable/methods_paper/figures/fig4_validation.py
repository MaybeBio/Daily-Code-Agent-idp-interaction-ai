"""
fig4_validation.py — Figure 4 (v3): predictions made in advance, against experimental ensembles and designed variants.
  a  (2026-10-14, degree beside hub occupancy) Partial Spearman on distance from the chain ends between predicted and
     ensemble labels, for the frozen hub occupancy (consensus) model and the frozen degree model, each next to its
     residue-neutral null: 18 development ensembles, 11 blind ensembles (open symbols: the 3 scored after the partial
     metric was fixed) and the 14 proteins of the second blind set (set A without IAPP, mean over entries per protein),
     with the pre-registered verdicts (ped_scores_all.csv, ped_v2_scores.csv, degree_wetlab/scores.csv).
  b  Raw Spearman of each model over the 29 development and blind ensembles next to the agreement between 16 pairs of
     independent ensembles (ped_ceiling_pairs.csv, hub_vs_degree/part_b.txt); dashed, the square root of that median.
  d  also shows (2026-10-14) the per-variant correlation for degree with the frozen CALVADOS-2 degree model (teal diamonds,
     CALVADOS-2 only, post hoc; cache/design_degree_per_variant.csv from fig4_design_degree_export.py)
  c-e  The prospective design batch (2026-10-11, panel functions from figS_design.py): one variant's simulated
       (Mpipi-Recharged) and predicted change along the chain, the per-variant correlation by substitution type in both
       force fields against the stickiness and charge rule and the noise ceiling, and the mean change at the substituted
       residues.
Input: ../../redo_evidence/hub_predictor/ped_scores_all.csv and ped_ceiling_pairs.csv (ped_blind11_report.py)
  ~/.venvs/bender/bin/python fig4_validation.py  ->  fig4_validation.pdf / .png
"""
import csv
import math
import os
import re

import matplotlib.pyplot as plt
import numpy as np

import figstyle as fs
import figS_design as SD

OCC_C, DEG_C, NULL_C, POS_C, CEIL_C = fs.FF_COLORS['CALVADOS-2'], '#2ca25f', '#999999', '#cccccc', '#555555'
NEW = {'PED00158', 'PED00001', 'PED00489'}


def strip(ax, x, vals, color, open_mask=None, label=None, edge=None):
    vals = np.asarray(vals, float)
    om = np.zeros(len(vals), bool) if open_mask is None else np.asarray(open_mask, bool)
    jit = np.random.default_rng(0).uniform(-0.12, 0.12, len(vals))
    ax.scatter(x + jit[~om], vals[~om], s=12, color=color, edgecolor=edge or 'none', lw=0.5, zorder=3, label=label)
    if om.any():
        ax.scatter(x + jit[om], vals[om], s=13, facecolor='white', edgecolor=color if color != POS_C else '0.45',
                   lw=0.9, zorder=3)
    ax.plot([x - 0.25, x + 0.25], [np.median(vals)] * 2, color='black', lw=1.2, zorder=4)
    ax.text(x, max(vals) + 0.04, f'{np.median(vals):+.2f}', ha='center', va='bottom', fontsize=5.4)


def per_protein(rows, model, null):
    by = {}
    for r in rows:
        by.setdefault(r['protein'], []).append((float(r[model]), float(r[null])))
    return [np.mean([x[0] for x in v]) for v in by.values()], [np.mean([x[1] for x in v]) for v in by.values()]


def main():
    occ = list(csv.DictReader(open(os.path.join(fs.EV, 'hub_predictor', 'ped_scores_all.csv'))))
    deg = list(csv.DictReader(open(os.path.join(fs.EV, 'degree_wetlab', 'scores.csv'))))
    v2 = list(csv.DictReader(open(os.path.join(fs.EV, 'ped_v2_scores.csv'))))
    prot = {r['entry']: r['protein'] for r in v2 if r['set'] == 'A'}
    pairs_o = [float(r['spearman']) for r in csv.DictReader(open(os.path.join(fs.EV, 'hub_predictor', 'ped_ceiling_pairs.csv')))]
    pairs_d = []
    for line in open(os.path.join(fs.EV, 'hub_vs_degree', 'part_b.txt')):
        m = re.match(r'\s+\S+\s+PED\d+-PED\d+\s+\d+\s+\|\s+[+-][0-9.]+\s+[+-][0-9.]+\s+\|\s+([+-][0-9.]+)', line)
        if m:
            pairs_d.append(float(m.group(1)))
    fig = plt.figure(figsize=(7.2, 3.35))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.52, left=0.075, right=0.99, top=0.95, bottom=0.135)
    top = gs[0].subgridspec(1, 2, width_ratios=[1.75, 1], wspace=0.2)
    ax = fig.add_subplot(top[0, 0])
    groups = []
    for sets, name in ((('used18',), 'development'), (('fresh8', 'extra3'), 'test')):
        O = [r for r in occ if r['set'] in sets]
        D = [r for r in deg if r['set'] in sets]
        groups.append((f'{len(O)} {name}', [float(r['partial_consensus']) for r in O], [float(r['partial_null']) for r in O],
                       [r['entry'] in NEW for r in O], [float(r['partial_model']) for r in D],
                       [float(r['partial_null']) for r in D], [r['entry'] in NEW for r in D]))
    vo = [r for r in v2 if r['set'] == 'A' and r['protein'] != 'IAPP']
    vd = [dict(r, protein=prot[r['entry']]) for r in deg if r['set'] == 'v2A' and prot.get(r['entry']) != 'IAPP']
    om_, on_ = per_protein(vo, 'partial_model', 'partial_null')
    dm_, dn_ = per_protein(vd, 'partial_model', 'partial_null')
    groups.append((f'{len(om_)} second test', om_, on_, None, dm_, dn_, None))
    verdict = {'development': ('', ''), 'test': ('met', 'met'), 'second test': ('not met', 'met')}
    for g, (lab, om, on, mask, dm, dn, dmask) in enumerate(groups):
        x0 = 4.6 * g
        first = g == 0
        strip(ax, x0, om, OCC_C, mask, label='hub occupancy model' if first else None)
        strip(ax, x0 + 1, on, NULL_C, mask, label='residue-neutral null' if first else None)
        strip(ax, x0 + 2.2, dm, DEG_C, dmask, label='degree model' if first else None)
        strip(ax, x0 + 3.2, dn, NULL_C, dmask)
        ax.text(x0 + 1.6, -0.66, lab + ' ensembles' if 'second' not in lab else lab.split()[0] + ' proteins, second test set', ha='center', va='top',
                fontsize=6.3)
        key = lab.split(' ', 1)[1]
        vo_, vd_ = verdict[key]
        if key != 'development':
            ax.text(x0 + 0.5, -0.86, vo_, ha='center', va='top', fontsize=5.6, color='#2F7DB5')
            ax.text(x0 + 2.7, -0.86, vd_, ha='center', va='top', fontsize=5.6, color='#2F7F74')
    ax.text(-0.75, -0.86, 'criterion', ha='right', va='top', fontsize=5.6, color='0.4')
    ax.scatter([], [], s=13, facecolor='white', edgecolor='0.3', lw=0.9, label='scored after the metric was fixed')
    ax.axhline(0, color='0.7', lw=0.5)
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.set_xlim(-0.6, 4.6 * 2 + 3.8)
    ax.set_ylim(-0.6, 1.22)
    ax.set_yticks([-0.5, 0, 0.5, 1.0])
    ax.set_ylabel('partial ρ, beyond position')
    ax.legend(frameon=False, loc='upper center', fontsize=5.4, handletextpad=0.1, borderaxespad=0.0, ncol=4,
              columnspacing=0.9, bbox_to_anchor=(0.5, 1.03))
    fs.panel(ax, 'a', x=-0.075, title='predictions recorded in advance against experimental ensembles (PED)')

    ax = fig.add_subplot(top[0, 1])
    allo = [float(r['rho_consensus']) for r in occ]
    alld = [float(r['rho_model']) for r in deg if r['set'] in ('used18', 'fresh8', 'extra3')]
    for x0, mod, prs, col, name in ((0, allo, pairs_o, OCC_C, 'hub occupancy'), (2.6, alld, pairs_d, DEG_C, 'degree')):
        strip(ax, x0, mod, col)
        strip(ax, x0 + 1, prs, CEIL_C)
        c = math.sqrt(max(np.median(prs), 0))
        ax.plot([x0 - 0.4, x0 + 1.4], [c] * 2, color=col, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.text(x0 + 0.5, -0.66, name, ha='center', va='top', fontsize=6.3, color='#2F7DB5' if col == OCC_C else '#2F7F74')
        ax.text(x0, -0.52, 'model', ha='center', va='bottom', fontsize=5.2, color='0.35')
        ax.text(x0 + 1, -0.52, 'ensemble\npairs', ha='center', va='bottom', fontsize=5.2, color='0.35')
    ax.axhline(0, color='0.7', lw=0.5)
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.set_xlim(-0.6, 4.1)
    ax.set_ylim(-0.6, 1.22)
    ax.set_yticks([-0.5, 0, 0.5, 1.0])
    ax.set_ylabel('ρ')
    fs.panel(ax, 'b', x=-0.16, title='against the agreement of experiments')

    # c-e: the prospective design batch (figS_design.py panel functions)
    allvar, _ = SD.load()
    bot = gs[1].subgridspec(1, 3, width_ratios=[1.45, 1, 0.95], wspace=0.42)
    ax = fig.add_subplot(bot[0, 0])
    SD.draw_example(ax, SD.pick_examples(allvar['Mpipi-Recharged'])[0], ylabel=True, legend='upper left', short=True)
    fs.panel(ax, 'c', x=-0.13, title='one designed variant')
    ax = fig.add_subplot(bot[0, 1])
    SD.draw_response(ax, allvar, legend_fontsize=4.9, short=True)
    dd = list(csv.DictReader(open(os.path.join(fs.HERE, 'cache', 'design_degree_per_variant.csv'))))
    rng = np.random.default_rng(1)
    for t, (kind, _) in enumerate(SD.TYPES):
        v = np.array([float(r['rho_degree']) for r in dd if r['kind'] == kind and r['rho_degree'] not in ('', 'nan')])
        ax.scatter(t + rng.uniform(-0.04, 0.04, len(v)), v, s=7, marker='D', color=DEG_C, edgecolor='white', lw=0.3,
                   zorder=5)
        ax.plot([t - 0.07, t + 0.07], [np.median(v)] * 2, color=DEG_C, lw=1.1, zorder=6)
    hd, ld = ax.get_legend_handles_labels() if ax.get_legend() is None else (ax.get_legend().legend_handles,
                                                                             [t.get_text() for t in ax.get_legend().get_texts()])
    hd = list(hd) + [plt.Line2D([], [], marker='D', ls='', ms=3, color=DEG_C)]
    ld = list(ld) + ['degree, CALVADOS-2']
    ax.legend(hd, ld, frameon=False, loc='upper center', fontsize=4.7, handletextpad=0.3, borderaxespad=0.0, ncol=2,
              columnspacing=0.5, handlelength=1.3)
    fs.panel(ax, 'd', x=-0.3, title='49 variants, two force fields')
    ax = fig.add_subplot(bot[0, 2])
    SD.draw_sites(ax, allvar, short=True)
    ax.legend(frameon=False, loc='lower right', fontsize=5.5, handletextpad=0.2, borderaxespad=0.1)
    fs.panel(ax, 'e', x=-0.32, title='at the substituted residues')
    fs.save(fig, 'fig4_validation')


if __name__ == '__main__':
    main()
