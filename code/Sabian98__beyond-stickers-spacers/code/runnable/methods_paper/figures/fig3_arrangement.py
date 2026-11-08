"""
fig3_arrangement.py — Figure 3: residue arrangement, not interaction parameters, sets where hubs form, and a sequence model learns it.
  a  Agreement between hub occupancy profiles at matched positions, with the distance from the chain ends removed
     (partial Spearman): a new arrangement at fixed composition and interaction parameters (three composition-preserving shuffles per
     protein, averaged, against the natural CALVADOS-2 profile) and new interaction parameters at fixed arrangement (every tyrosine to
     phenylalanine against the natural CALVADOS-2 profile, and all three pairs of force fields).
  b-d  The residue-type part of hub occupancy (excess of each residue type over the residue-neutral null, the 907 sequences
     simulated with all three force fields) against each force field's own interaction scale (CALVADOS-2 lambda, mean
     Wang-Frenkel epsilon of Mpipi-GG from its run.in files and of Mpipi-Recharged from OpenMpipi's recharged_params.txt),
     both standardised across the 20 types.
  e  5-fold CV by protein in each force field: chain position only, + own residue type, + sticker and spacer context,
     and the full model (context_model_by_ff.py), next to the CALVADOS-2 noise ceiling.
  f  The learned map on arrangements and interaction parameters it never saw: held-out natural CALVADOS-2 proteins (CV), the
     frozen CALVADOS-2 model on the simulated shuffles and on the tyrosine-to-phenylalanine mutants, and transfer to
     Mpipi-GG and Mpipi-Recharged next to the agreement of the simulations themselves.
Inputs: redo_evidence/hub_frames_v3_{cal,shuffle,yf,recharged,mpipi}, sim_null_v3_by_type.csv,
        hub_predictor/cv_full_per_protein.csv, gaivi/predictions_frozen_2026-10-02_arrangement_{simseq_,}2026-10-04.csv,
        mpipi_epsilon_table.json (built from Mpipi-GG run.in files on first use), mpipi_recharged_parameters.json
        (OpenMpipi commit 9c9541b, read 2026-10-04), CALVADOS-2 residue table.
2026-10-14: degree (teal) beside hub occupancy in a (arrangement_test_degree.py, plan a16218bf, all proteins as here),
e (medians and IQR of both readouts, context_model_by_ff_per_protein{,_degree}.csv) and f (cv_full_per_protein_degree.csv
for natural, transfer and sims; the frozen CALVADOS-2 degree model 3a77ea2a on the shuffles, as D5 of that test).
b-d stay hub occupancy with the residue-neutral null (no degree null exists).
  ~/.venvs/bender/bin/python fig3_arrangement.py  ->  fig3_arrangement.pdf / .png, fig3_arrangement_values.txt
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr

import figstyle as fs

sys.path.insert(0, fs.EV)
import arrangement_test as at  # noqa: E402
import arrangement_test_degree as atd  # noqa: E402

HP = os.path.join(fs.EV, 'hub_predictor')
GV = os.path.join(fs.EV, 'gaivi')
CAL_TABLE = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))
EPS_TABLE = os.path.join(fs.EV, 'mpipi_epsilon_table.json')
RC_TABLE = os.path.join(fs.EV, 'mpipi_recharged_parameters.json')
AAS = 'ACDEFGHIKLMNPQRSTVWY'
SHUFFLE_COLOR, YF_COLOR, GR_COLOR = '#8c8c8c', '#6baed6', '#b39ddb'
DEG_C = '#2ca25f'
out = []


def box(ax, data, pos, color, width=0.55, fmt='{:.2f}', fontsize=6, rotation=0):
    b = ax.boxplot(data, positions=[pos], widths=width, showfliers=False, patch_artist=True,
                   medianprops=dict(color='black', lw=1.0), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                   boxprops=dict(lw=0.6))
    b['boxes'][0].set_facecolor(color)
    b['boxes'][0].set_alpha(1.0)
    q1, q3 = np.percentile(data, [25, 75])
    top = min(np.max(data), q3 + 1.5 * (q3 - q1))
    ax.text(pos, top + 0.03, fmt.format(np.median(data)), ha='center', va='bottom', fontsize=fontsize, rotation=rotation)


def mpipi_epsilon():
    if os.path.exists(EPS_TABLE):
        return json.load(open(EPS_TABLE))
    import hub_frames_mpipi as hm
    eps = {}
    for u in sorted(os.listdir(hm.RUNS))[:300]:
        d = os.path.join(hm.RUNS, u)
        try:
            masses = open(os.path.join(d, 'seq.dat')).read().split('Masses')[1].split('Atoms')[0]
            t2aa = {int(m.group(1)): m.group(2) for m in re.finditer(r'^\s*(\d+)\s+[\d.]+\s+#\s*([A-Z])', masses, re.M)}
            for m in re.finditer(r'^pair_coeff\s+(\d+)\s+(\d+)\s+wf/cut\s+([\d.eE+-]+)\s+([\d.eE+-]+)',
                                 open(os.path.join(d, 'run.in')).read(), re.M):
                a, b = t2aa.get(int(m.group(1))), t2aa.get(int(m.group(2)))
                if a and b:
                    eps[a + b] = eps[b + a] = float(m.group(3))
        except Exception:
            continue
        if len(eps) >= 400:
            break
    json.dump(eps, open(EPS_TABLE, 'w'))
    return eps


def partial_agreement(A, B, pairs):
    v = []
    for ua, ub in pairs:
        oa, ob = A[ua][1], B[ub][1]
        if len(oa) == len(ob) and oa.std() > 0 and ob.std() > 0:
            v.append(at.partial(oa, ob))
    return np.array([x for x in v if x == x])


def spearman_agreement(A, B, pairs):
    v = []
    for ua, ub in pairs:
        oa, ob = A[ua][1], B[ub][1]
        if len(oa) == len(ob) and oa.std() > 0 and ob.std() > 0:
            v.append(spearmanr(oa, ob).correlation)
    return np.array([x for x in v if x == x])


def main():
    SH = at.read('hub_frames_v3_shuffle')
    YF = at.read('hub_frames_v3_yf')
    RCH = at.read('hub_frames_v3_recharged')
    GG = at.read('hub_frames_v3_mpipi')
    wt_sh = defaultdict(list)
    for u in SH:
        wt_sh[u.rsplit('_global_', 1)[0]].append(u)
    wt_yf = {u[:-3]: u for u in YF}
    CAL = at.read('hub_frames_v3_cal', set(wt_sh) | set(wt_yf) | set(RCH) | set(GG))
    cohort = sorted(u for u in wt_sh if u in CAL and len(wt_sh[u]) == 3)

    fig = plt.figure(figsize=(7.2, 3.2))
    gs = fig.add_gridspec(2, 2, hspace=0.78, wspace=0.3, left=0.085, right=0.98, top=0.94, bottom=0.15,
                          width_ratios=[1, 1.1])

    # A  arrangement versus interaction parameters
    ax = fig.add_subplot(gs[0, 0])
    sh = np.array([np.mean(partial_agreement(CAL, SH, [(u, v) for v in wt_sh[u]])) for u in cohort])
    yf = partial_agreement(CAL, YF, [(u, wt_yf[u]) for u in wt_yf if u in CAL])
    rc = partial_agreement(CAL, RCH, [(u, u) for u in RCH if u in CAL])
    gg = partial_agreement(CAL, GG, [(u, u) for u in GG if u in CAL])
    gr = partial_agreement(GG, RCH, [(u, u) for u in GG if u in RCH])      # the third pair, without CALVADOS-2
    # degree, the same comparisons (arrangement_test_degree.py, plan a16218bf; all proteins as for occupancy here)
    Dsh, Dyf = atd.read('hub_frames_v3_shuffle', 'deg_mean'), atd.read('hub_frames_v3_yf', 'deg_mean')
    Drc, Dgg = atd.read('hub_frames_v3_recharged', 'deg_mean'), atd.read('hub_frames_v3_mpipi', 'deg_mean')
    Dcal = atd.read('hub_frames_v3_cal', 'deg_mean', set(CAL))
    dsh = np.array([np.mean(partial_agreement(Dcal, Dsh, [(u, v) for v in wt_sh[u]])) for u in cohort])
    dyf = partial_agreement(Dcal, Dyf, [(u, wt_yf[u]) for u in wt_yf if u in Dcal])
    drc = partial_agreement(Dcal, Drc, [(u, u) for u in Drc if u in Dcal])
    dgg = partial_agreement(Dcal, Dgg, [(u, u) for u in Dgg if u in Dcal])
    dgr = partial_agreement(Dgg, Drc, [(u, u) for u in Dgg if u in Drc])
    items = [(0, sh, dsh, 'shuffled', SHUFFLE_COLOR), (1.3, yf, dyf, 'Y→F', YF_COLOR),
             (2.4, gg, dgg, 'CALVADOS-2\nvs Mpipi-GG', fs.FF_COLORS['MPIPI-GG']),
             (3.6, rc, drc, 'CALVADOS-2\nvs Recharged', fs.FF_COLORS['Mpipi-Recharged']),
             (4.8, gr, dgr, 'Mpipi-GG\nvs Recharged', GR_COLOR)]
    for x, d, dd, lab, col in items:
        box(ax, d, x - 0.2, col, width=0.36, fontsize=5)
        box(ax, dd, x + 0.2, DEG_C, width=0.36, fontsize=5)
        out.append(f'A  degree {lab.replace(chr(10), " "):28s} n={len(dd):5d}  median partial {np.median(dd):.2f}')
    ax.text(0.01, 0.99, 'left box hub occupancy, right box degree', transform=ax.transAxes, ha='left', va='top',
            fontsize=5.4, color='#2F7F74')
    ax.set_xticks([i[0] for i in items])
    ax.set_xticklabels([i[3] for i in items], fontsize=5.5)
    tr3 = ax.get_xaxis_transform()
    for xg, g in ((0, 'new arrangement'), (3.6, 'same arrangement, new interaction parameters')):
        ax.annotate(g, (xg, 0), xycoords=tr3, xytext=(0, -22), textcoords='offset points', ha='center', va='top',
                    fontsize=6)
    ax.axhline(0, color='0.7', lw=0.5)
    ax.set_xlim(-0.6, 5.4)
    ax.set_ylim(-0.5, 1.3)
    ax.set_yticks([-0.5, 0, 0.5, 1.0])
    ax.set_ylabel('agreement,\nposition removed')
    fs.panel(ax, 'a', x=-0.2, title='a new arrangement erases both readouts, new parameters do not')
    for name, d in (('shuffled (mean of 3)', sh), ('Y to F', yf), ('Mpipi-GG', gg), ('Mpipi-Recharged', rc),
                    ('Mpipi-GG vs Recharged', gr)):
        out.append(f'A  {name:22s} n={len(d):5d}  median partial {np.median(d):.2f}')

    # B  the residue-type part follows each force field's parameters: one small panel per force field
    lam = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(CAL_TABLE))}
    z = lambda d: {a: (d[a] - np.mean(list(d.values()))) / np.std(list(d.values())) for a in AAS}
    eps = mpipi_epsilon()
    epsm = {a: np.mean([eps[a + b] for b in AAS if a + b in eps]) for a in AAS}
    rc = json.load(open(RC_TABLE))['epsilon_kJ_per_mol']               # Mpipi-Recharged's own table (2026-10-04)
    epsr = {a: np.mean([rc[a + b] for b in AAS]) for a in AAS}
    scales = {'λ': z(lam), 'ε GG': z(epsm), 'ε Recharged': z(epsr)}
    ex = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(fs.EV, 'sim_null_v3_by_type.csv'))):
        if r['rule'] == 'fixed':
            ex[r['force_field']][r['aa']] = float(r['excess'])
    sub = gs[0, 1].subgridspec(1, 3, wspace=0.42)
    # (force field, own scale, its symbol, y-limits, residues to label with offsets in points)
    panels = (('CALVADOS-2', scales['λ'], 'λ', (-0.45, 0.7),
               {'W': (3, -2), 'Y': (3, 0), 'F': (-6, 2), 'R': (3, -5), 'G': (3, -1), 'E': (2, -6), 'K': (2, -6)}),
              ('MPIPI-GG', scales['ε GG'], 'ε', (-0.95, 0.85),
               {'W': (-7, -1), 'Y': (3, -1), 'F': (-6, 2), 'G': (3, -1), 'P': (3, -1), 'K': (-7, 2)}),
              ('Mpipi-Recharged', scales['ε Recharged'], 'ε', (-0.75, 3.0),
               {'W': (-8, -2), 'Y': (3, -1), 'F': (3, -4), 'G': (-2, 4), 'K': (-3, 4)}))
    for j, (ff, par, lab, ylim, labels) in enumerate(panels):
        ax = fig.add_subplot(sub[0, j])
        xs = np.array([par[a] for a in AAS])
        ys = np.array([ex[ff][a] for a in AAS])
        r = pearsonr(xs, ys)[0]
        col = fs.FF_COLORS[ff]
        name = {'MPIPI-GG': 'Mpipi-GG'}.get(ff, ff)
        ax.axhline(0, color='0.8', lw=0.5, zorder=1)
        k, b0 = np.polyfit(xs, ys, 1)
        xx = np.linspace(xs.min(), xs.max(), 10)
        ax.plot(xx, k * xx + b0, color=col, lw=0.8, alpha=0.8, zorder=2)
        ax.scatter(xs, ys, s=8, color=col, edgecolor='white', linewidth=0.3, zorder=3)
        for a, (dx, dy) in labels.items():
            ax.annotate(a, (par[a], ex[ff][a]), xytext=(dx, dy), textcoords='offset points', fontsize=5.5, color=col,
                        fontweight='bold', zorder=4)
        tag = name
        ax.text(0.05, 0.97, tag + '\nr = ' + f'{r:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=6, color=col)
        ax.set_xlim(-2.2, 3.0)
        ax.set_ylim(*ylim)
        ax.set_xticks([-2, 0, 2])
        ax.tick_params(labelsize=6, length=2, pad=1.5)
        fs.panel(ax, 'bcd'[j], x=-0.2 if j == 0 else -0.36, y=1.05)
        if j == 0:
            ax.set_ylabel('hub bonus over null')
        if j == 1:
            ax.set_xlabel('interaction strength in each force field (standardized)')
        out.append(f'B  {ff:16s} excess vs own parameter r = {r:.2f}; vs ' + '; '.join(
            f'{k2} r = {pearsonr([v[a] for a in AAS], ys)[0]:.2f}' for k2, v in scales.items()))

    # C  the context ladder in all three force fields (context_model_by_ff.py, same 5-fold protein CV in each)
    by = defaultdict(list)                     # full-model CV and transfer, used again in panel D
    for r in csv.DictReader(open(os.path.join(HP, 'cv_full_per_protein.csv'))):
        if r['spearman'] not in ('', 'nan'):
            by[(r['force_field'], r['condition'])].append(float(r['spearman']))
    byc = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HP, 'context_model_by_ff_per_protein.csv'))):
        if r['spearman'] not in ('', 'nan'):
            byc[(r['force_field'], r['condition'])].append(float(r['spearman']))
    ceil = [float(r['ceiling']) for r in csv.DictReader(open(os.path.join(HP, 'cv_full_per_protein.csv')))
            if r['condition'] == 'hgb' and r['force_field'] == 'CALVADOS-2' and r['ceiling'] not in ('', 'nan')]
    bycd = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HP, 'context_model_by_ff_per_protein_degree.csv'))):
        if r['spearman'] not in ('', 'nan'):
            bycd[(r['force_field'], r['condition'])].append(float(r['spearman']))
    ax = fig.add_subplot(gs[1, 0])
    conds = [('position', 'position\nonly'), ('identity', '+ own\ntype'), ('context own', '+ sticker and\nspacer context'),
             ('full', 'full\nsequence')]
    ffs = [('CALVADOS-2', 'CALVADOS-2'), ('Mpipi-GG', 'MPIPI-GG'), ('Mpipi-Recharged', 'Mpipi-Recharged')]
    xs = np.arange(len(conds))
    for k, (ff, colkey) in enumerate(ffs):
        col = fs.FF_COLORS[colkey]
        dx = (k - 1) * 0.07
        for data, ls, mfc, lab in ((byc, (0, (3, 2)), 'white', 'hub occupancy'), (bycd, '-', col, 'degree')):
            med = [np.median(data[(ff, c)]) for c, _ in conds]
            q = np.array([np.percentile(data[(ff, c)], [25, 75]) for c, _ in conds])
            ax.vlines(xs + dx, q[:, 0], q[:, 1], color=col, lw=0.6, alpha=0.6, zorder=1)
            ax.plot(xs + dx, med, ls=ls, color=col, lw=0.9, marker='o', ms=3.2, mfc=mfc, mec=col, mew=0.8, zorder=2)
            for c, _ in conds:
                out.append(f'C  {lab:13s} {ff:16s} {c:12s} n={len(data[(ff, c)]):5d}  median {np.median(data[(ff, c)]):.3f}')
    ax.axhline(np.median(ceil), color='0.4', lw=0.6, ls=(0, (1, 1.5)), zorder=0)
    ax.text(-0.35, np.median(ceil) + 0.012, 'noise ceiling, hub occupancy in CALVADOS-2', fontsize=5, va='bottom',
            color='0.35')
    out.append(f'C  noise ceiling CALVADOS-2 n={len(ceil)} median {np.median(ceil):.3f}')
    ax.set_xticks(xs)
    ax.set_xticklabels([c[1] for c in conds], fontsize=6)
    ax.set_xlim(-0.4, len(conds) - 0.6)
    ax.set_ylabel('Spearman ρ')
    ax.set_ylim(0, 1.05)
    hd = [plt.Line2D([], [], color='0.3', lw=0.9, marker='o', ms=3.2, label='degree'),
          plt.Line2D([], [], color='0.3', lw=0.9, ls=(0, (3, 2)), marker='o', ms=3.2, mfc='white', label='hub occupancy')]
    hd += [plt.Line2D([], [], color=fs.FF_COLORS[c], lw=1.5, label=f'{ff} ({len(byc[(ff, "full")]):,})') for ff, c in ffs]
    ax.legend(handles=hd, frameon=False, fontsize=5.2, loc='lower right', ncol=1, handlelength=1.8, borderaxespad=0.1,
              bbox_to_anchor=(1.0, 0.0))
    fs.panel(ax, 'e', x=-0.2, y=1.04, title='context, not type, predicts both readouts')

    # D  generalisation to new arrangements and interaction parameters
    ax = fig.add_subplot(gs[1, 1])
    # frozen CALVADOS-2 model: shuffles from the simulated (top.pdb) sequences, Y-to-F mutants from the first file
    pred = defaultdict(list)
    for f, keep in (('predictions_frozen_2026-10-02_arrangement_simseq_2026-10-04.csv', lambda u: '_global_' in u),
                    ('predictions_frozen_2026-10-02_arrangement_2026-10-04.csv', lambda u: u.endswith('_YF'))):
        for r in csv.DictReader(open(os.path.join(GV, f))):
            if keep(r['id']):
                pred[r['id']].append(float(r['cal']))

    def model_acc(D, keys):
        v = []
        for u in keys:
            if u in pred and len(pred[u]) == len(D[u][1]) and D[u][1].std() > 0:
                v.append(spearmanr(pred[u], D[u][1]).correlation)
        return np.array([x for x in v if x == x])
    acc_sh = model_acc(SH, list(SH))
    import pickle
    mdeg = pickle.load(open(os.path.join(fs.EV, 'hub_vs_degree', 'degree_model_cal.pkl'), 'rb'))
    dacc_sh = np.array([spearmanr(mdeg.predict(atd.T.features(Dsh[u][0], 10)), Dsh[u][1]).correlation
                        for u in Dsh if Dsh[u][1].std() > 0])
    byd = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HP, 'cv_full_per_protein_degree.csv'))):
        if r['spearman'] not in ('', 'nan'):
            byd[(r['force_field'], r['condition'])].append(float(r['spearman']))
    acc_yf = model_acc(YF, list(YF))
    sims_gg = spearman_agreement(CAL, GG, [(u, u) for u in GG if u in CAL])
    sims_rc = spearman_agreement(CAL, RCH, [(u, u) for u in RCH if u in CAL])
    # the Y-to-F mutants are left out here: their wild types were in the training set and they barely differ from them,
    # whereas a shuffle agrees with its wild type at only 0.12 (panel A), so it is a genuine test of new arrangements
    items = [(0, by[('CALVADOS-2', 'hgb')], byd[('CALVADOS-2', 'hgb')], 'natural', fs.FF_COLORS['CALVADOS-2']),
             (1, acc_sh, dacc_sh, 'shuffled', SHUFFLE_COLOR),
             (2.3, by[('MPIPI-GG', 'transfer')], byd[('MPIPI-GG', 'transfer')], 'model', fs.FF_COLORS['MPIPI-GG']),
             (3.1, sims_gg, byd[('MPIPI-GG', 'sims')], 'sims', 'white'),
             (4.4, by[('Mpipi-Recharged', 'transfer')], byd[('Mpipi-Recharged', 'transfer')], 'model',
              fs.FF_COLORS['Mpipi-Recharged']),
             (5.2, sims_rc, byd[('Mpipi-Recharged', 'sims')], 'sims', 'white')]
    for x, d, dd, lab, col in items:
        box(ax, d, x - 0.17, col, width=0.3, fontsize=4.8, rotation=90)
        box(ax, dd, x + 0.17, DEG_C, width=0.3, fontsize=4.8, rotation=90)
        out.append(f'D  degree {lab:10s} at {x}: n={len(dd):5d} median {np.median(dd):.3f}')
    ax.set_xticks([i[0] for i in items])
    ax.set_xticklabels([i[3] for i in items], fontsize=6)
    tr3 = ax.get_xaxis_transform()
    for xg, g in ((0.5, 'CALVADOS-2'), (2.7, 'Mpipi-GG'), (4.8, 'Mpipi-Recharged')):
        ax.annotate(g, (xg, 0), xycoords=tr3, xytext=(0, -22), textcoords='offset points', ha='center', va='top',
                    fontsize=6)
    ax.axhline(0, color='0.7', lw=0.5)
    ax.set_xlim(-0.6, 5.8)
    ax.set_ylim(-0.5, 1.3)
    ax.set_ylabel('Spearman ρ')
    fs.panel(ax, 'f', x=-0.18, title='unseen arrangements and force fields (left occupancy, right degree)')
    for name, d in (('natural CV', by[('CALVADOS-2', 'hgb')]), ('frozen model on shuffles', acc_sh), ('frozen model on Y to F', acc_yf),
                    ('transfer Mpipi-GG', by[('MPIPI-GG', 'transfer')]), ('sims CAL vs GG', sims_gg),
                    ('transfer Recharged', by[('Mpipi-Recharged', 'transfer')]), ('sims CAL vs Recharged', sims_rc)):
        out.append(f'D  {name:26s} n={len(d):5d}  median {np.median(d):.3f}')
    fs.save(fig, 'fig3_arrangement')
    open(os.path.join(fs.HERE, 'fig3_arrangement_values.txt'), 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
