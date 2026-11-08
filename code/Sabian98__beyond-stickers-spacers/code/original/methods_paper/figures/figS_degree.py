"""
figS_degree.py — Supplementary figure: contact degree beside hub occupancy. Values from the hashed results in
redo_evidence (hub_vs_degree/, hub_predictor/*_degree*, degree_wetlab/); no new statistics.
  a  held-out accuracy per CALVADOS-2 protein, degree against hub occupancy
  b  model analyses, degree against hub occupancy
  c  grouped permutation importance of the context model: window G/S fraction and window stickiness
  d  PED ensembles: partial correlation of the frozen degree model against that of the occupancy model
  e  glycine around the top 10% of residues by position-adjusted degree or hub occupancy
  ~/.venvs/finches/bin/python figS_degree.py  ->  figS_degree.pdf / .png
"""
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figstyle as fs

EV = fs.EV
HP = os.path.join(EV, 'hub_predictor')
CD, CO = '#2ca25f', '#1b6ca8'          # degree (teal, as in Figs 1 and 2), hub occupancy
FFS = ['CALVADOS-2', 'Mpipi-GG', 'Mpipi-Recharged']


def importance(path):
    out = {}
    for line in open(path):
        m = re.match(r'\s+(CALVADOS-2|Mpipi-GG|Mpipi-Recharged)\s+baseline [0-9.]+: (.*)', line)
        if m:
            d = dict((k.strip(), float(v)) for k, v in re.findall(r'([a-zA-Z/ -]+?) ([+-][0-9.]+)', m.group(2)))
            out[m.group(1)] = d
    return out


def main():
    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(2, 3, hspace=0.66, wspace=0.45, left=0.075, right=0.985, top=0.94, bottom=0.12)

    # a: paired per-protein accuracy
    d = pd.read_csv(os.path.join(HP, 'cv_full_per_protein_degree.csv'))
    o = pd.read_csv(os.path.join(HP, 'cv_full_per_protein.csv'))
    p = d[d.condition == 'hgb'].merge(o[(o.condition == 'hgb') & (o.force_field == 'CALVADOS-2')], on='protein_id',
                                      suffixes=('_d', '_o'))
    ax = fig.add_subplot(gs[0, 0])
    ax.hexbin(p.spearman_o, p.spearman_d, gridsize=40, cmap='Greys', mincnt=1, bins='log', extent=(0.3, 1, 0.3, 1))
    ax.plot([0.3, 1], [0.3, 1], color='0.25', lw=0.7, ls=(0, (3, 2)))
    ax.set_xlim(0.3, 1)
    ax.set_ylim(0.3, 1)
    ax.set_xlabel('hub occupancy model, ρ')
    ax.set_ylabel('degree model, ρ')
    ax.text(0.05, 0.9, f'degree higher in {100 * (p.spearman_d > p.spearman_o).mean():.0f}%', transform=ax.transAxes, fontsize=6)
    fs.panel(ax, 'a', x=-0.3, y=1.08, title='10,281 held-out proteins')

    # b: analyses
    rows = [('full\nmodel', 0.934, 0.894), ('ridge', 0.896, 0.842), ('position\nonly', 0.367, 0.440),
            ('transfer\nRech', 0.724, 0.677), ('transfer\nGG', 0.538, 0.533), ('unseen\nCAL', 0.63, 0.61),
            ('unseen\nGG', 0.56, 0.54), ('unseen\nRech', 0.67, 0.64), ('consensus', 0.92, 0.90)]
    ax = fig.add_subplot(gs[0, 1:])
    x = np.arange(len(rows))
    ax.bar(x - 0.18, [r[1] for r in rows], 0.36, color=CD, label='degree')
    ax.bar(x + 0.18, [r[2] for r in rows], 0.36, color=CO, label='hub occupancy')
    ax.set_xticks(x, [r[0] for r in rows], fontsize=6)
    ax.set_ylabel('median per-protein ρ')
    ax.set_ylim(0, 1.12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.legend(frameon=False, fontsize=6, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.04))
    fs.panel(ax, 'b', x=-0.12, y=1.08, title='sequence models of the two readouts')

    # c: permutation importance
    imp_d = importance(os.path.join(HP, 'context_model_by_ff_summary_degree.txt'))
    imp_o = importance(os.path.join(HP, 'context_model_by_ff_summary.txt'))
    ax = fig.add_subplot(gs[1, 0])
    for k, ff in enumerate(FFS):
        for j, (feat, hatch) in enumerate((('window G/S fraction', None), ('window stickiness', '////'))):
            xo = k + (j - 0.5) * 0.4
            ax.bar(xo - 0.09, imp_d[ff][feat], 0.18, color=CD, hatch=hatch, edgecolor='white', lw=0)
            ax.bar(xo + 0.09, imp_o[ff][feat], 0.18, color=CO, hatch=hatch, edgecolor='white', lw=0)
    ax.set_xticks(range(3), ['CAL', 'GG', 'Rech'])
    ax.set_ylabel('importance (drop in ρ)')
    ax.text(0.02, 0.92, 'solid: G/S fraction   hatched: stickiness', transform=ax.transAxes, fontsize=5.5)
    fs.panel(ax, 'c', x=-0.3, y=1.08, title='spacers matter more for degree')

    # d: PED partial correlations
    sd = pd.read_csv(os.path.join(EV, 'degree_wetlab', 'scores.csv'))
    so = pd.read_csv(os.path.join(HP, 'ped_scores_all.csv'))
    sv = pd.read_csv(os.path.join(EV, 'ped_v2_scores.csv'))
    a = sd[sd.set.isin(['used18', 'fresh8', 'extra3'])].merge(so[['entry', 'partial_consensus']], on='entry')
    b = sd[sd.set == 'v2A'].merge(sv[sv.set == 'A'][['entry', 'protein', 'partial_model']], on='entry', suffixes=('', '_o'))
    b = b[b.protein_o != 'IAPP'].groupby('protein_o')[['partial_model', 'partial_model_o']].mean()   # per protein, as Fig. 4a
    ax = fig.add_subplot(gs[1, 1])
    for name, sel, col in (('development', a.set == 'used18', '#bdbdbd'), ('test', a.set.isin(['fresh8', 'extra3']), '#1b6ca8')):
        ax.scatter(a.partial_consensus[sel], a.partial_model[sel], s=11, color=col, lw=0, label=name)
    ax.scatter(b.partial_model_o, b.partial_model, s=11, color='#7570b3',
               lw=0, label='second test')
    ax.plot([-0.5, 1], [-0.5, 1], color='0.6', lw=0.6, ls=(0, (3, 2)))
    ax.axhline(0, color='0.85', lw=0.5)
    ax.axvline(0, color='0.85', lw=0.5)
    ax.set_xlim(-0.5, 1)
    ax.set_ylim(-0.5, 1)
    ax.set_xlabel('occupancy model, partial ρ')
    ax.set_ylabel('degree model, partial ρ')
    ax.legend(frameon=False, fontsize=5.5, loc='lower right', handletextpad=0.1)
    fs.panel(ax, 'd', x=-0.27, y=1.08, title='ensembles built from experiment')

    # e: glycine around top residues
    nd = pd.read_csv(os.path.join(HP, 'context_neighbourhood_degree.csv'))
    no = pd.read_csv(os.path.join(HP, 'context_neighbourhood.csv'))
    ax = fig.add_subplot(gs[1, 2])
    sets = ['CALVADOS-2', 'Mpipi-GG', 'Mpipi-Recharged', 'PED']
    for k, st in enumerate(sets):
        for df, col, dx in ((nd, CD, -0.18), (no, CO, 0.18)):
            r = df[(df.label_set == st) & (df.neighbour_class == 'G')].iloc[0]
            ax.bar(k + dx, r.enrichment, 0.36, color=col)
            ax.errorbar(k + dx, r.enrichment, yerr=[[r.enrichment - r.lo], [r.hi - r.enrichment]], color='k', lw=0.6, capsize=1.5)
    ax.axhline(1, color='0.5', lw=0.5, ls=(0, (3, 2)))
    ax.set_xticks(range(4), ['CAL', 'GG', 'Rech', 'PED'])
    ax.set_ylabel('glycine ratio')
    fs.panel(ax, 'e', x=-0.3, y=1.08, title='glycine around top residues')
    ax.set_ylim(0, 2.25)
    fs.save(fig, 'figS_degree')
    print('saved; paired n', len(p), 'PED dev/blind', len(a), 'v2A proteins', len(b))


if __name__ == '__main__':
    main()
