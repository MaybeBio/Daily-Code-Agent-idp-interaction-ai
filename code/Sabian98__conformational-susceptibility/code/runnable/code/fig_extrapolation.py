"""
fig_extrapolation.py — Supplementary figure for the retrospective analysis of the held-out misses
(heldout_extrapolation.py; post hoc, exploratory).
  a  kernel rows: change of the simulated slope per prior SD of each of the ten forcing directions (hydrophilic Asp/Glu,
     Lys, Arg, N/Q/S/T/H, Gly; hydrophobic TEA aliphatic, Pro, aromatic, Arg, other), for the four calibration
     constructs, lambda-repressor (not calibrated) and the eight held-out series, CALVADOS-2 (left) and Mpipi-Recharged
     (right); row labels give each protein's Gly content
  b  miss of the calibrated model against the part of the row outside the span of the four calibration rows, 14 measured
     held-out series (lambda-repressor open), with the pooled Spearman correlation and its permutation P
  c  the part outside the span of each held-out series split into Gly, other hydrophilic and hydrophobic directions
     (bar length = the part outside the span; segments in proportion to their squared components)
Inputs: heldout_temperature/results/exploratory/extrapolation_{rows,series}.csv, extrapolation_summary.txt; sequences
from the Wuttke and held-out FASTA files.
    python fig_extrapolation.py            (openmpipi env)
Output: BENDER-MT/manuscript/figs/fig_extrapolation.{pdf,png}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402
import re  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402

EX = f'{FX.ROOT}/heldout_temperature/results/exploratory'
OUT = f'{FX.ROOT}/manuscript/figs/fig_extrapolation'
FASTA = [os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/CALVADOS IDP/input_wuttke/wuttke2014_constructs.fasta')),
         os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_HELDOUT/calvados/forms.fasta'))]
HYDG = ['DE', 'K', 'R', 'pol', 'G']
TEAC = ['aliph', 'pro', 'arom', 'arg', 'rest']
DIRS = [f'hyd_{g}' for g in HYDG] + [f'tea_{c}' for c in TEAC]
COLLAB = ['Asp, Glu', 'Lys', 'Arg', 'N Q S T H', 'Gly', 'aliphatic', 'Pro', 'aromatic', 'Arg', 'other']
ORDER = ['CspM34', 'IN_NTD', 'ProTa_C2_C56', 'ProTa_C56_C110', 'lambda_P6C_I84C', 'PNt-I025', 'PNt-I175', 'PNt-I525',
         'sfAFP-I175', 'asyn4W26C-I030', 'asyn39W66C-I030', 'AapPGR-I195', 'p53N93-I118']
NAME = {'CspM34': 'CspM34', 'IN_NTD': 'Integrase NTD', 'ProTa_C2_C56': 'ProTα, N half', 'ProTa_C56_C110': 'ProTα, C half',
        'lambda_P6C_I84C': 'λ-repressor', 'PNt-I025': 'PNt, 25 mM', 'PNt-I175': 'PNt, 175 mM', 'PNt-I525': 'PNt, 525 mM',
        'sfAFP-I175': 'sfAFP', 'asyn4W26C-I030': 'α-syn 4–26', 'asyn39W66C-I030': 'α-syn 39–66', 'AapPGR-I195': 'Aap',
        'p53N93-I118': 'p53(1–93)'}
SEQKEY = {'PNt-I025': 'PNt', 'PNt-I175': 'PNt', 'PNt-I525': 'PNt', 'sfAFP-I175': 'sfAFP', 'asyn4W26C-I030': 'asyn4W26C',
          'asyn39W66C-I030': 'asyn39W66C', 'AapPGR-I195': 'AapPGR', 'p53N93-I118': 'p53N93'}
MODELS = ['CALVADOS-2', 'Mpipi-Recharged']
LABEL_OFF = {('CALVADOS-2', 'AapPGR-I195'): (-4, -8), ('CALVADOS-2', 'lambda_P6C_I84C'): (4, 3)}   # keep labels apart
HYDCOL, TEACOL = '#2E8C6F', '#C8792A'            # figfont roles: hydrophilic forcing (teal), hydrophobic (ochre)


def gly_content():
    seq, k = {}, None
    for f in FASTA:
        for line in open(f):
            line = line.strip()
            if line.startswith('>'):
                k = line[1:].split()[0]
                seq[k] = ''
            elif k:
                seq[k] += line
    return {s: seq[SEQKEY.get(s, s)].count('G') / len(seq[SEQKEY.get(s, s)]) for s in ORDER}


def pooled_stat():
    txt = open(f'{EX}/extrapolation_summary.txt').read()
    m = re.search(r'miss vs outside_span: pooled rho ([0-9.]+) \(one-sided permutation P ([0-9.]+)', txt)
    return float(m.group(1)), float(m.group(2))


def heat(ax, R, model, vmax, show_y, gly):
    M = R[R.model == model].set_index('series').loc[ORDER, DIRS].to_numpy()
    im = ax.imshow(M, cmap='RdBu', vmin=-vmax, vmax=vmax, aspect='auto')
    ax.set_xticks(range(len(DIRS)))
    ax.set_xticklabels(COLLAB, rotation=60, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([f'{NAME[s]} ({100 * gly[s]:.0f}% Gly)' for s in ORDER] if show_y else [])
    ax.tick_params(length=0)
    for y in (3.5, 4.5):
        ax.axhline(y, color='white', lw=2.0)
    ax.axvline(4.5, color='white', lw=2.0)
    ax.set_title(model, loc='left', color=FX.COL[model], pad=15)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return im


def main():
    FX.house()
    R = pd.read_csv(f'{EX}/extrapolation_rows.csv')
    S = pd.read_csv(f'{EX}/extrapolation_series.csv')
    gly = gly_content()
    fig = plt.figure(figsize=(FX.WIDTH, 6.6))
    top = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.045], left=0.235, right=0.93, top=0.955, bottom=0.555, wspace=0.08)
    a1, a2, cax = fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1]), fig.add_subplot(top[0, 2])
    vmax = 0.30
    assert R[DIRS].abs().to_numpy().max() < vmax, 'kernel entry beyond the colour range'
    im = heat(a1, R, MODELS[0], vmax, True, gly)
    heat(a2, R, MODELS[1], vmax, False, gly)
    for ax in (a1, a2):
        ax.text(2.0, -0.9, 'hydrophilic', ha='center', va='bottom', fontsize=6.5, color=HYDCOL)
        ax.text(7.0, -0.9, 'hydrophobic (TEA)', ha='center', va='bottom', fontsize=6.5, color=TEACOL)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Change of slope per prior SD\nof the forcing (% per K)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_visible(False)
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '0' if v == 0 else FX.minus(f'{v:.1f}')))
    # b  miss against the part outside the span
    bot = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], left=0.1, right=0.985, top=0.42, bottom=0.085, wspace=0.42)
    b, c = fig.add_subplot(bot[0, 0]), fig.add_subplot(bot[0, 1])
    lab_pts = {'PNt-I025': 'PNt 25', 'PNt-I175': 'PNt 175', 'sfAFP-I175': 'sfAFP', 'asyn39W66C-I030': 'α-syn',
               'AapPGR-I195': 'Aap', 'lambda_P6C_I84C': 'λ'}
    for m in MODELS:
        g = S[(S.model == m) & S.miss.notna()]
        h = g[g.series != 'lambda_P6C_I84C']
        b.scatter(h.outside_span, h.miss, s=18, color=FX.COL[m], zorder=3, label=m, lw=0)
        lam = g[g.series == 'lambda_P6C_I84C']
        b.scatter(lam.outside_span, lam.miss, s=18, facecolor='white', edgecolor=FX.COL[m], lw=0.9, zorder=3)
        for _, r in g.iterrows():
            if r.series in lab_pts and not (r.series == 'asyn39W66C-I030' and m == MODELS[0]):
                off = LABEL_OFF.get((m, r.series), (3, 1))
                b.annotate(lab_pts[r.series], (r.outside_span, r.miss), xytext=off, textcoords='offset points',
                           fontsize=6.0, color=FX.COL[m], ha='right' if off[0] < 0 else 'left')
    rho, p = pooled_stat()
    b.text(0.97, 0.04, f'ρ = {rho:.2f}, P = {p:.3f}\n14 held-out series', transform=b.transAxes, ha='right', va='bottom')
    b.set_xlim(0, 0.25)
    b.set_ylim(0, 0.95)
    b.set_xlabel('Part of the row outside the\ncalibration span (% per K)')
    b.set_ylabel('Miss of the linear form\n|predicted − measured| (% per K)')
    b.legend(loc='upper left', handletextpad=0.2, borderaxespad=0.1)
    # c  part outside the span by direction
    held = ORDER[5:]
    y, yt, yl = 0.0, [], []
    for s in held:
        for k, m in enumerate(MODELS):
            r = S[(S.model == m) & (S.series == s)].iloc[0]
            out = np.array([r[f'out_{d}'] for d in DIRS])
            tot = r.outside_span
            parts = [out[4] ** 2, (out[:4] ** 2).sum(), (out[5:] ** 2).sum()]
            parts = np.array(parts) / sum(parts) * tot
            left = 0.0
            for v, col, a_ in zip(parts, (HYDCOL, HYDCOL, TEACOL), (1.0, 0.35, 1.0)):
                c.barh(y, v, left=left, height=0.36, color=col, alpha=a_, lw=0)
                left += v
            c.plot([-0.0045], [y], 's', ms=3.0, color=FX.COL[m], clip_on=False)
            y -= 0.42
        yt.append(y + 0.63)
        yl.append(NAME[s])
        y -= 0.35
    c.set_yticks(yt)
    c.set_yticklabels(yl)
    c.tick_params(axis='y', length=0, pad=11)
    c.set_xlim(0, 0.25)
    c.set_ylim(y + 0.2, 0.4)
    c.set_xlabel('Part of the row outside the\ncalibration span (% per K)')
    h = [matplotlib.patches.Patch(color=HYDCOL, label='Gly'), matplotlib.patches.Patch(color=HYDCOL, alpha=0.35, label='other hydrophilic'),
         matplotlib.patches.Patch(color=TEACOL, label='hydrophobic')]
    c.legend(handles=h, loc='lower right', handlelength=1.0, handletextpad=0.3, borderaxespad=0.1)
    c.spines['left'].set_visible(False)
    FX.add_letters(fig, [(a1, 'a', {'x': 0.012}), (b, 'b', {'row': 2}), (c, 'c', {'row': 2})])
    FX.save(fig, OUT)
    print('written', OUT, {'rho': rho, 'p': p, 'gly': {k: round(v, 2) for k, v in gly.items()}})


if __name__ == '__main__':
    main()
