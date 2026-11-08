"""
figS_context.py — three supplementary figures for the sticker and spacer analyses.
  figS_sticker_definitions  enrichment among hubs (top 10% of each protein by occupancy) of seven sticker sets and two
                            spacer sets in each label set, with 95% bootstrap intervals over proteins; rings mark each
                            force field's own four stickiest residue types (sticker_definitions.csv)
  figS_context_ladder       median held-out per-protein Spearman correlation in each force field for chain position,
                            + residue type, + sticker and spacer context, + proline in that context, and the full model
                            (context_model_by_ff_per_protein.csv, context_proline_full_per_protein.csv)
  figS_neighbourhood        A: composition of the +-5 window around hubs relative to residues of the same type, three
                            force fields and the PED ensembles (context_neighbourhood_typed.csv); B: glycine and proline
                            around hubs in the residue-neutral null and in the excess over it, 907 sequences
                            (context_neighbourhood_null_typed.csv)
  ~/.venvs/bender/bin/python figS_context.py  ->  figS_sticker_definitions.pdf/.png, figS_context_ladder.pdf/.png,
                                                   figS_neighbourhood.pdf/.png
"""
import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import figstyle as fs

HP = os.path.join(fs.EV, 'hub_predictor')
FF = [('CALVADOS-2', fs.FF_COLORS['CALVADOS-2']), ('Mpipi-GG', fs.FF_COLORS['MPIPI-GG']),
      ('Mpipi-Recharged', fs.FF_COLORS['Mpipi-Recharged'])]
PED = ('PED', '0.25')
INK = '#1D1D1F'


def rd(name):
    return list(csv.DictReader(open(os.path.join(HP, name))))


def sticker_definitions():
    get = {(r['group'], r['label_set']): r for r in rd('sticker_definitions.csv') if r['measure'] == 'enrichment among hubs'}
    rows = [('stickers F, W, Y, R (canonical)', 'F W Y R (canonical)'), ('stickers F, W, Y, H', 'F W Y H'),
            ('aromatic F, W, Y', 'F W Y'), ('Y, F', 'Y F'), ('Y, R', 'Y R'), ('F, W, Y, R, H', 'F W Y R H'),
            ('F, W, Y, R, K', 'F W Y R K'), None, ('spacers G, S (canonical)', 'G S (canonical)'),
            ('spacers G, S, N, Q', 'G S N Q')]
    own = {('stickers F, W, Y, R (canonical)', 'CALVADOS-2'), ('stickers F, W, Y, H', 'Mpipi-GG'),
           ('stickers F, W, Y, H', 'Mpipi-Recharged')}
    sets = FF + [PED]
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    ys, labs, y = [], [], 0.0
    for row in rows:
        if row is None:
            y -= 0.6
            continue
        g, lab = row
        for j, (name, col) in enumerate(sets):
            r = get[(g, name)]
            yy = y + (1.5 - j) * 0.17
            e, lo, hi = float(r['estimate']), float(r['lo']), float(r['hi'])
            ax.plot([lo, hi], [yy, yy], color=col, lw=0.9, solid_capstyle='butt', zorder=2)
            ax.plot(e, yy, 'D' if name == 'PED' else 'o', ms=3.4, color=col, mec='white', mew=0.3, zorder=3)
            if (g, name) in own:
                ax.plot(e, yy, 'o', ms=7, mfc='none', mec=col, mew=0.8, zorder=3)
        ys.append(y)
        labs.append(lab)
        y -= 1.0
    ax.axvline(1, color='0.4', lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_xscale('log', base=2)
    ax.set_xlim(0.33, 9.5)
    ax.set_xticks([0.5, 1, 2, 4, 8])
    ax.set_xticklabels(['0.5×', '1×', '2×', '4×', '8×'])
    ax.minorticks_off()
    ax.set_yticks(ys)
    ax.set_yticklabels(labs)
    ax.set_ylim(y + 0.45, 1.0)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('share among hubs ÷ share in sequence')
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    handles = [Line2D([], [], marker='D' if n == 'PED' else 'o', color=c, ms=4, lw=0.9, mec='white') for n, c in sets]
    handles.append(Line2D([], [], marker='o', ms=7, mfc='none', mec=INK, mew=0.8, lw=0))
    ax.legend(handles, ['CALVADOS-2', 'Mpipi-GG', 'Mpipi-Recharged', 'PED, built from experiment',
                        "the force field's own four stickiest types"], loc='upper left', bbox_to_anchor=(0.6, 0.6),
              frameon=False, handlelength=1.6)
    fig.tight_layout()
    fs.recolor(fig)
    for ext in ('pdf', 'png'):
        fig.savefig(f'figS_sticker_definitions.{ext}', dpi=300)
    plt.close(fig)


def context_ladder():
    med = defaultdict(list)
    for f in ('context_model_by_ff_per_protein.csv', 'context_proline_full_per_protein.csv'):
        for r in rd(f):
            if r['spearman'] not in ('', 'nan'):
                med[(r['force_field'], r['condition'])].append(float(r['spearman']))
    # an explicit light-to-dark ramp (already in the light palette, so figstyle.recolor leaves it alone)
    conds = [('position', 'chain position only', '#d9d9d9'), ('identity', '+ type of the residue', '#D2E7F5'),
             ('context own', '+ sticker and spacer context', '#9CCBEA'),
             ('context own + proline', '+ proline in that context', '#5A9FD2'), ('full', 'full sequence model', '#2C5F8C')]
    fig, ax = plt.subplots(figsize=(6.5, 2.6))
    w = 0.16
    for i, (ff, _) in enumerate(FF):
        for j, (c, lab, col) in enumerate(conds):
            v = np.median(med[(ff, c)])
            x = i + (j - 2) * w
            ax.bar(x, v, width=w * 0.92, color=col, label=lab if i == 0 else None, zorder=2)
            ax.text(x, v + 0.015, f'{v:.2f}', ha='center', va='bottom', fontsize=6, color=INK)
    ax.set_xticks(range(len(FF)))
    ax.set_xticklabels([f'{ff} ({len(med[(ff, "full")]):,} proteins)' for ff, _ in FF])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel('median per-protein Spearman\n(held-out proteins)')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=5, frameon=False, handlelength=1.0, columnspacing=1.0)
    fig.tight_layout()
    fs.recolor(fig)
    for ext in ('pdf', 'png'):
        fig.savefig(f'figS_context_ladder.{ext}', dpi=300)
    plt.close(fig)


def neighbourhood():
    nb = {(r['label_set'], r['neighbour_class']): r for r in rd('context_neighbourhood_typed.csv')}
    nn = {(r['force_field'], r['hubs_by'], r['neighbour_class']): r for r in rd('context_neighbourhood_null_typed.csv')}
    classes = [('G', 'glycine'), ('S', 'serine'), ('aromatic F W Y', 'aromatic'), ('P', 'proline'), ('K', 'lysine'),
               ('acidic D E', 'D, E')]
    sets = FF + [PED]
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.7), gridspec_kw={'width_ratios': [2.1, 1]})
    w = 0.2
    for j, (name, col) in enumerate(sets):
        for i, (cls, _) in enumerate(classes):
            r = nb[(name, cls)]
            e, lo, hi = float(r['enrichment']), float(r['lo']), float(r['hi'])
            x = i + (j - 1.5) * w
            a.bar(x, e - 1, bottom=1, width=w * 0.9, color=col, zorder=2,
                  label=(name if name != 'PED' else 'PED, built from experiment') if i == 0 else None)
            a.plot([x, x], [lo, hi], color=INK, lw=0.6, zorder=3)
    a.axhline(1, color='0.4', lw=0.7, ls=(0, (3, 2)), zorder=1)
    a.set_yscale('log', base=2)
    a.set_ylim(0.42, 2.6)
    a.set_yticks([0.5, 0.71, 1, 1.41, 2])
    a.set_yticklabels(['0.5×', '0.71×', '1×', '1.41×', '2×'])
    a.minorticks_off()
    a.set_xticks(range(len(classes)))
    a.set_xticklabels([c[1] for c in classes])
    a.set_ylabel('composition of the neighbors of hubs\nrelative to residues of the same type')
    fs.panel(a, 'a', x=-0.11, y=1.05, title='what surrounds hubs')
    a.legend(loc='upper right', frameon=False, ncol=2, handlelength=1.0, columnspacing=0.8)
    xs = np.arange(len(FF))
    for k, (kind, col, lab) in enumerate((('null', '#c7c7cc', 'residue-neutral null'),
                                          ('excess', '#2ca25f', 'excess over the null'))):
        for i, (ff, _) in enumerate(FF):
            for m, (cls, mk) in enumerate((('G', 'o'), ('P', 's'))):
                r = nn[(ff, kind, cls)]
                e, lo, hi = float(r['enrichment']), float(r['lo']), float(r['hi'])
                x = i + (k - 0.5) * 0.34 + (m - 0.5) * 0.12
                b.plot([x, x], [lo, hi], color=col if kind == 'excess' else '0.55', lw=0.7, zorder=2)
                b.plot(x, e, mk, ms=3.6, color=col if kind == 'excess' else '0.55', mec='white', mew=0.3, zorder=3)
    b.axhline(1, color='0.4', lw=0.7, ls=(0, (3, 2)), zorder=1)
    b.set_xticks(xs)
    b.set_xticklabels([ff.replace('Mpipi-', 'Mpipi-\n') for ff, _ in FF])
    b.set_ylabel('around hubs, relative to\nresidues of the same type')
    fs.panel(b, 'b', x=-0.3, y=1.05, title='null and excess, 907 sequences')
    handles = [Line2D([], [], marker='o', color='0.55', lw=0, ms=4), Line2D([], [], marker='o', color='#2ca25f', lw=0, ms=4),
               Line2D([], [], marker='o', color=INK, lw=0, ms=4, mfc='white'), Line2D([], [], marker='s', color=INK, lw=0, ms=4, mfc='white')]
    b.legend(handles, ['null (geometry only)', 'excess over the null', 'glycine', 'proline'], loc='upper center',
             bbox_to_anchor=(0.5, -0.24), frameon=False, ncol=2, handlelength=0.8, columnspacing=0.6, fontsize=5.8)
    for ax in (a, b):
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fs.recolor(fig)
    for ext in ('pdf', 'png'):
        fig.savefig(f'figS_neighbourhood.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    sticker_definitions()
    context_ladder()
    neighbourhood()
    print('wrote figS_sticker_definitions, figS_context_ladder, figS_neighbourhood')
