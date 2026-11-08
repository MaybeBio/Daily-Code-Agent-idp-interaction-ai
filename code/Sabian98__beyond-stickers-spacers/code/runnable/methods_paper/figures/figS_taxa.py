"""
figS_taxa.py — hub occupancy across the six BENDER groups (taxa_plan.md, hashed 2026-10-11; descriptive).
  a  held-out per-protein Spearman of the full model (5-fold CV) by group, with the position-only baseline and the
     median noise ceiling
  b  leave-one-group-out: model trained on the other five groups against the 5-fold CV score, same proteins
  c  transfer to Mpipi-Recharged against the agreement between CALVADOS-2 and Mpipi-Recharged, one point per group,
     for the 5-fold and the leave-one-group-out models
  d  enrichment of the canonical stickers (F, W, Y, R) and of G, S among hubs, per group and force field
Inputs: hub_predictor/cv_full_per_protein.csv, hub_predictor/cv_taxa_per_protein.csv, agreement_ceiling_per_protein.csv,
        taxa/hub_composition_by_group.csv
  ~/.venvs/bender/bin/python figS_taxa.py  ->  figS_taxa.pdf / .png, figS_taxa_values.txt
"""
import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

import figstyle as fs

EV = fs.EV
GROUPS = ['Bacteria', 'Fungi', 'Mammals', 'Plants', 'Protists', 'Viruses']
GCOL = dict(zip(GROUPS, ['#0072b2', '#e69f00', '#d6006f', '#00a651', '#7b3294', '#5f5f5f']))
GMK = dict(zip(GROUPS, ['o', 's', 'D', '^', 'v', 'P']))


def main():
    grp = {r['UniProt_ID']: r['group'] for r in csv.DictReader(open(os.path.join(EV, 'hub_frames_v3_cal', 'per_protein.csv')))}
    cv, ceil = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, 'hub_predictor', 'cv_full_per_protein.csv'))):
        g = grp.get(r['protein_id'])
        if g not in GROUPS or r['spearman'] in ('', 'nan'):
            continue
        cv[(r['force_field'], r['condition'], g)].append(float(r['spearman']))
        if r['condition'] == 'hgb' and r['ceiling'] not in ('', 'nan'):
            ceil[g].append(float(r['ceiling']))
    lo = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, 'hub_predictor', 'cv_taxa_per_protein.csv'))):
        if r['spearman'] not in ('', 'nan'):
            lo[(r['condition'], r['group'])].append(float(r['spearman']))
    ag = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, 'agreement_ceiling_per_protein.csv'))):
        if r['spearman_cal_rech'] not in ('', 'nan'):
            ag[r['group']].append(float(r['spearman_cal_rech']))
    comp = {(r['force_field'], r['group'], r['set']): float(r['enrichment'])
            for r in csv.DictReader(open(os.path.join(EV, 'taxa', 'hub_composition_by_group.csv')))}
    med = lambda v: float(np.median(v))
    lines = ['group  n  cv_hgb  logo  position  logo_position  ceiling  transfer_cv  transfer_logo  agreement']

    fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.35), gridspec_kw={'width_ratios': [1.35, 1, 1, 1.3]})
    x = np.arange(len(GROUPS))
    # a
    a = ax[0]
    for k, g in enumerate(GROUPS):
        v = cv[('CALVADOS-2', 'hgb', g)]
        b = a.boxplot(v, positions=[k], widths=0.55, showfliers=False, patch_artist=True,
                      medianprops=dict(color='k', lw=1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                      boxprops=dict(lw=0.6))
        b['boxes'][0].set_facecolor(GCOL[g])
    a.plot(x, [med(cv[('CALVADOS-2', 'position', g)]) for g in GROUPS], '_', color='k', ms=10, mew=1.4,
           label='position only')
    a.plot(x, [med(ceil[g]) for g in GROUPS], '_', color='#d6006f', ms=10, mew=1.4, label='noise ceiling')
    a.set_xticks(x, GROUPS, fontsize=5.8, rotation=40, ha='right', rotation_mode='anchor')
    a.set_ylim(-0.2, 1.05)
    a.set_ylabel('Spearman ρ (held out)')
    a.legend(frameon=False, fontsize=5.5, loc='lower right', handlelength=0.8)
    fs.panel(a, 'a', x=-0.3, title='5-fold CV by group')
    # b
    a = ax[1]
    for g in GROUPS:
        a.scatter(med(cv[('CALVADOS-2', 'hgb', g)]), med(lo[('logo', g)]), s=22, color=GCOL[g], marker=GMK[g],
                  edgecolor='white', lw=0.4, label=g, zorder=3)
    lim = (0.84, 0.93)
    a.plot(lim, lim, color='0.6', lw=0.6, ls=(0, (3, 2)))
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel('5-fold CV')
    a.set_ylabel('group held out')
    a.legend(frameon=False, fontsize=5, loc='upper left', handletextpad=0.1, borderaxespad=0.1, labelspacing=0.3)
    fs.panel(a, 'b', x=-0.38, title='leave one group out')
    # c
    a = ax[2]
    for g in GROUPS:
        a.scatter(med(ag[g]), med(cv[('Mpipi-Recharged', 'transfer', g)]), s=22, color=GCOL[g], marker=GMK[g],
                  edgecolor='white', lw=0.4, zorder=3)
        a.scatter(med(ag[g]), med(lo[('logo_transfer', g)]), s=22, facecolor='white', edgecolor=GCOL[g],
                  marker=GMK[g], lw=0.9, zorder=3)
    a.plot([0.55, 0.8], [0.55, 0.8], color='0.6', lw=0.6, ls=(0, (3, 2)))
    a.set_xlim(0.6, 0.78)
    a.set_ylim(0.58, 0.76)
    a.set_xlabel('CALVADOS-2 vs Recharged')
    a.set_ylabel('model on Recharged')
    a.scatter([], [], s=20, color='0.4', label='5-fold')
    a.scatter([], [], s=20, facecolor='white', edgecolor='0.4', label='group held out')
    a.legend(frameon=False, fontsize=5.3, loc='upper left', handletextpad=0.2)
    fs.panel(a, 'c', x=-0.38, title='transfer')
    # d
    a = ax[3]
    w = 0.2
    for j, (ff, st, col, hatch) in enumerate((('CALVADOS-2', 'stickers', fs.FF_COLORS['CALVADOS-2'], None),
                                              ('Mpipi-Recharged', 'stickers', fs.FF_COLORS['Mpipi-Recharged'], None),
                                              ('CALVADOS-2', 'GS', fs.FF_COLORS['CALVADOS-2'], '////'),
                                              ('Mpipi-Recharged', 'GS', fs.FF_COLORS['Mpipi-Recharged'], '////'))):
        vals = [comp[(ff, g, st)] for g in GROUPS]
        a.bar(x + (j - 1.5) * w, vals, width=w, color=col if hatch is None else 'white', edgecolor=col,
              hatch=hatch, lw=0.6, label=f'{"FWYR" if st == "stickers" else "GS"}, {"CAL-2" if ff.startswith("C") else "Rech."}')
    a.axhline(1, color='0.4', lw=0.6)
    a.set_yscale('log')
    a.set_yticks([0.5, 1, 2, 4], ['0.5', '1', '2', '4'])
    a.yaxis.set_minor_locator(plt.NullLocator())
    a.set_xticks(x, GROUPS, fontsize=5.8, rotation=40, ha='right', rotation_mode='anchor')
    a.set_ylabel('enrichment among hubs')
    a.legend(frameon=False, fontsize=5, loc='upper center', ncol=2, handlelength=1.0, columnspacing=0.6, borderaxespad=0.1)
    a.set_ylim(0.3, 12)
    fs.panel(a, 'd', x=-0.3, title='which residues are hubs')
    fig.subplots_adjust(left=0.06, right=0.995, top=0.87, bottom=0.26, wspace=0.55)
    fs.save(fig, 'figS_taxa')
    for g in GROUPS:
        lines.append(f'{g:9s} {len(cv[("CALVADOS-2", "hgb", g)]):5d}  {med(cv[("CALVADOS-2", "hgb", g)]):.3f}  '
                     f'{med(lo[("logo", g)]):.3f}  {med(cv[("CALVADOS-2", "position", g)]):.3f}  '
                     f'{med(lo[("logo_position", g)]):.3f}  {med(ceil[g]):.3f}  '
                     f'{med(cv[("Mpipi-Recharged", "transfer", g)]):.3f}  {med(lo[("logo_transfer", g)]):.3f}  '
                     f'{med(ag[g]):.3f}  stickers CAL {comp[("CALVADOS-2", g, "stickers")]:.2f} '
                     f'RCH {comp[("Mpipi-Recharged", g, "stickers")]:.2f}  GS CAL {comp[("CALVADOS-2", g, "GS")]:.2f} '
                     f'RCH {comp[("Mpipi-Recharged", g, "GS")]:.2f}')
    open(os.path.join(fs.HERE, 'figS_taxa_values.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
