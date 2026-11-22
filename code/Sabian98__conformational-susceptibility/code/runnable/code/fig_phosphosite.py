"""
Deck figure for the phosphosite test (CALVADOS-2): where annotated phosphosites fall among their protein's Ser/Thr
ranked by leverage |dRg| (quintiles within each protein), and odds ratios from conditional logistic regression
(groups = proteins) with kinase-context and position controls.
Inputs: Analysis_FDT/phosphosite_test/sites.csv (phosphosite_test.py)
Output: Analysis_FDT/phosphosite_test/fig_phosphosite.{svg,png}
"""
import os
import sys

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from statsmodels.discrete.conditional_models import ConditionalLogit  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import phosphosite_test as pt  # noqa: E402

BG, INK, MUTED, GRID = '#F7F4EE', '#1B2432', '#5B6472', '#E2DCCF'
ACC, CAL, REF = '#8F5E14', '#2E6FA3', '#B9C3CF'
plt.rcParams.update({'font.family': 'Helvetica Neue', 'font.size': 15, 'svg.fonttype': 'path', 'axes.edgecolor': INK,
                     'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})


def odds(W, cols):
    Wz = W.dropna(subset=['z_abs'])
    X = Wz[cols].astype(float)
    X = (X - X.mean()) / X.std()
    r = ConditionalLogit(Wz.annotated.astype(int).to_numpy(), X.to_numpy(), groups=Wz.UniProt_ID.to_numpy()).fit(disp=False)
    ci = np.exp(r.conf_int())
    return {c: (np.exp(r.params[k]), ci[k, 0], ci[k, 1], r.pvalues[k]) for k, c in enumerate(cols)}


def main():
    S = pd.read_csv(os.path.join(pt.OUT, 'sites.csv'))
    W = pt.within_protein(S[S.model == 'CALVADOS-2'])
    W['q'] = np.minimum((W.pct * 5).astype(int), 4)
    fig, axes = plt.subplots(1, 2, figsize=(16.6, 5.4), facecolor=BG, gridspec_kw={'width_ratios': [1, 1.15]})
    ax = axes[0]
    ax.set_facecolor(BG)
    x = np.arange(5)
    for k, (lab, m, c) in enumerate((('all Ser/Thr (reference)', np.ones(len(W), bool), REF),
                                      ('annotated phosphosites', W.annotated.to_numpy(), ACC))):
        f = np.bincount(W.q[m], minlength=5) / m.sum() * 100
        ax.bar(x + (k - 0.5) * 0.38, f, width=0.36, color=c, label=f'{lab} (n = {m.sum()})')
        if k == 1:
            for xx, v in zip(x, f):
                ax.text(xx + 0.19, v + 0.6, f'{v:.0f}%', ha='center', fontsize=13, color=ACC)
    ax.set_xticks(x)
    ax.set_xticklabels(['lowest', '2', '3', '4', 'highest'])
    ax.set_xlabel('Leverage |ΔRg| within the protein (fifths of its Ser/Thr)', fontsize=15)
    ax.set_ylabel('Share of sites (%)', fontsize=15)
    ax.set_title('Where real phosphosites fall', fontsize=17, loc='left', color=INK, pad=10)
    ax.legend(frameon=False, fontsize=13, loc='upper left')
    ax.set_ylim(0, 34)
    for s_ in ('top', 'right'):
        ax.spines[s_].set_visible(False)
    ax.grid(axis='y', color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax = axes[1]
    ax.set_facecolor(BG)
    one = odds(W, ['z_abs'])
    full = odds(W, ['z_abs', 'pro_plus1', 'local_DE', 'local_KR', 'term_dist'])
    rows = [('Leverage |ΔRg| (alone)', one['z_abs'], ACC), ('Leverage |ΔRg| (with controls)', full['z_abs'], ACC),
            ('Pro at +1 (kinase motif)', full['pro_plus1'], MUTED), ('Acidic neighbours (±5)', full['local_DE'], MUTED),
            ('Lys/Arg neighbours (±5)', full['local_KR'], MUTED), ('Distance to chain end', full['term_dist'], MUTED)]
    y = np.arange(len(rows))[::-1]
    ax.axvline(1, color=INK, lw=1.2)
    for yy, (lab, (o, lo, hi, p), c) in zip(y, rows):
        ax.plot([lo, hi], [yy, yy], color=c, lw=3, solid_capstyle='round')
        ax.plot(o, yy, 'o', color=c, ms=10)
        ax.text(3.05, yy, f'{o:.2f}', va='center', fontsize=13, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=14)
    ax.set_xscale('log')
    ax.set_xlim(0.7, 3.4)
    ax.set_xticks([0.75, 1, 1.5, 2, 3])
    ax.set_xticklabels(['0.75', '1', '1.5', '2', '3'])
    ax.set_xlabel('Odds of being a phosphosite, per SD (95% CI)', fontsize=15)
    ax.set_title('Within-protein comparison (conditional logit)', fontsize=17, loc='left', color=INK, pad=10)
    for s_ in ('top', 'right', 'left'):
        ax.spines[s_].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.tight_layout(w_pad=3)
    for ext, kw in (('svg', {}), ('png', {'dpi': 110})):
        fig.savefig(os.path.join(pt.OUT, f'fig_phosphosite.{ext}'), facecolor=BG, **kw)
    print({k: round(v[0], 2) for k, v in full.items()}, 'quintiles annotated %:',
          (np.bincount(W.q[W.annotated], minlength=5) / W.annotated.sum() * 100).round(1).tolist())


if __name__ == '__main__':
    main()
