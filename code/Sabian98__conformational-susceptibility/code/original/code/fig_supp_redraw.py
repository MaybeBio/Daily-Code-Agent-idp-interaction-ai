"""
fig_supp_redraw.py — Supplementary Figs S2, S6 and S8 redrawn in the house style of fig_supplementary.py (6.5 in wide,
vector PDF + 300 dpi PNG), replacing the deck-style rasters drawn by fig_deck_v2.py / fig_total_response.py (through
fig_arial_deckfigs.py). Same quantities as those functions, from the stored outputs only.
  compact_factors  (Fig. S2, was fig_compact_factors_v2.png, fig_deck_v2.compact_factors) median sticker
                   susceptibility of Rg at 300 K by compaction bin and its three factors: coupling r(Rg, E_s), relative
                   breathing CV(Rg) and sticker-energy spread sd(E_s)/kT per N^1/2; CALVADOS-2 and Mpipi-GG
                   (Analysis_FDT_MPIPI/compact/chi_factor_split/bins.csv)
  phosphosite      (Fig. S6, was fig_phosphosite_v2.png, fig_deck_v2.phosphosite) annotated phosphosites per
                   within-protein fifth of leverage |dRg| among their proteins' Ser/Thr (CALVADOS-2), and odds ratios
                   per SD from conditional logistic regression with kinase-context and position controls
                   (Analysis_FDT/phosphosite_test/sites.csv; odds() as in fig_phosphosite.py)
  lineage          (Fig. S8, was fig_lineage_total.png, fig_total_response.lineage_total) extra contraction beyond
                   compaction (minus the leftover of the total response) by taxonomic group, as simulated and adjusted
                   for pi and Asp+Glu+Lys content (CALVADOS-2, calibrated hydration); correlation of the residue-age score
                   with the extra contraction, raw (Spearman) and at fixed composition (rank-based partial), CALVADOS-2
                   and Mpipi-Recharged (Analysis_FDT/total_response/lineage/)
    BENDER_REANALYSIS_DIR=... python -B fig_supp_redraw.py [compact_factors phosphosite lineage]      (openmpipi env)
Output: BENDER-MT/manuscript/figs/{fig_compact_factors, fig_phosphosite, fig_lineage}.{pdf,png}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import fig_supplementary as S  # noqa: E402  (style, palette, letters, save; applies figfont)
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtr  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy import stats  # noqa: E402

COL, RG, NU, W = S.COL, S.RG, S.NU, S.W
MPC = f'{S.TLA}/Analysis_FDT_MPIPI/compact'
LIN = f'{S.TR}/lineage'
PALE_CAL = '#A8CBE8'                                                  # CALVADOS-2 pale tint (remapped by figfont)
DARK = '#4D4D4D'


def compact_factors():
    b = pd.read_csv(f'{MPC}/chi_factor_split/bins.csv')
    panels = [('chi_s_rel', f'Sticker susceptibility of {RG}\n(% per K)', True),
              ('r', f'Coupling, $r$({RG}, $E_\\mathrm{{s}}$)', False),
              ('cv_rg', f'Breathing, CV of {RG}', False),
              ('sdE_kT_per_sqrtN', f'Energy spread, sd($E_\\mathrm{{s}}$)/$k_\\mathrm{{B}}T$\nper $N^{{1/2}}$', False)]
    fig, axes = plt.subplots(2, 2, figsize=(W, 4.6))
    axes = axes.ravel()
    x = np.arange(len(S.BIN_LAB))
    out = {}
    for ax, (col, lab, zero) in zip(axes, panels):
        for m in ('CALVADOS-2', 'Mpipi-GG'):
            g = b[b.model == m]
            ax.plot(x[:len(g)], g[col], '-o', color=COL[m], lw=1.2, ms=3.5,
                    label=f'{m} ({int(g.n.sum()):,} chains)')
            out[f'{col} {m}'] = g[col].round(3).tolist()
        ax.axvline(1.5, color='0.55', lw=0.6, ls=':', zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels(S.BIN_LAB, rotation=35, ha='right', rotation_mode='anchor')
        ax.set_xlabel(f'Compaction bin, {NU} at 300 K')
        ax.set_ylabel(lab)
        v = b[col].to_numpy()
        pad = 0.08 * (v.max() - v.min())
        ax.set_ylim(0 if zero else v.min() - pad, v.max() + pad)
        ax.set_xlim(-0.4, len(x) - 0.6)
    axes[0].text(1.58, 0.004, f'{NU} = 0.5', fontsize=6.5, color='0.4', va='bottom')
    axes[0].legend(loc='upper right')
    out['n'] = b.set_index(['model', 'b']).n.to_dict()
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=2.0, w_pad=2.5)
    S.letters(fig, axes, 'abcd')
    S.save(fig, 'fig_compact_factors')
    return out


def _odds(W_, cols):
    """Conditional-logit odds ratios per SD with 95% CI (copied from fig_phosphosite.odds)."""
    from statsmodels.discrete.conditional_models import ConditionalLogit
    Wz = W_.dropna(subset=['z_abs'])
    X = Wz[cols].astype(float)
    X = (X - X.mean()) / X.std()
    r = ConditionalLogit(Wz.annotated.astype(int).to_numpy(), X.to_numpy(), groups=Wz.UniProt_ID.to_numpy()).fit(disp=False)
    ci = np.exp(r.conf_int())
    return {c: (np.exp(r.params[k]), ci[k, 0], ci[k, 1], r.pvalues[k]) for k, c in enumerate(cols)}


def phosphosite():
    import phosphosite_test as pt
    St = pd.read_csv(os.path.join(pt.OUT, 'sites.csv'))
    Wp = pt.within_protein(St[(St.model == 'CALVADOS-2') & (St.ess_frac > pt.ESS_MIN)])   # reliable reweighting only
    Wp['q'] = np.minimum((Wp.pct * 5).astype(int), 4)
    n_ann = int(Wp.annotated.sum())
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.6), gridspec_kw={'width_ratios': [1, 1.3]})
    ax = axes[0]
    f = np.bincount(Wp.q[Wp.annotated], minlength=5) / n_ann * 100
    x = np.arange(5)
    ax.bar(x, f, width=0.66, color=[PALE_CAL] * 4 + [COL['CALVADOS-2']], linewidth=0)
    for xx, v in zip(x, f):
        ax.text(xx, 0.8, f'{v:.0f}%', ha='center', va='bottom', fontsize=6.5, color='white' if xx == 4 else 'black')
    ax.axhline(20, color='black', lw=0.7, ls='--')
    ax.text(-0.42, 20.5, 'expected without enrichment', fontsize=6.5, va='bottom')
    ax.set_xticks(x)
    ax.set_xticklabels(['lowest', '2', '3', '4', 'highest'])
    ax.set_xlabel(f'Within-protein fifth of leverage |Δ{RG}|')
    ax.set_ylabel(f'Annotated phosphosites (%, n = {n_ann})')
    ax.set_ylim(0, 31)
    ax.set_xlim(-0.5, 4.5)
    ax = axes[1]
    one = _odds(Wp, ['z_abs'])
    full = _odds(Wp, ['z_abs', 'pro_plus1', 'local_DE', 'local_KR', 'term_dist'])
    cal = COL['CALVADOS-2']
    rows = [('Leverage, alone', one['z_abs'], cal), ('Leverage, with controls', full['z_abs'], cal), None,
            ('Pro at +1 (kinase motif)', full['pro_plus1'], DARK), ('Acidic neighbours (±5)', full['local_DE'], DARK),
            ('Lys/Arg neighbours (±5)', full['local_KR'], DARK), ('Distance to chain end', full['term_dist'], DARK)]
    y = np.arange(len(rows))[::-1].astype(float)
    ax.axvline(1, color='black', lw=0.6)
    labels, ticks = [], []
    for yy, row in zip(y, rows):
        if row is None:
            continue
        lab, (o, lo, hi, p), c = row
        ax.plot([lo, hi], [yy, yy], color=c, lw=1.4, solid_capstyle='butt')
        ax.plot(o, yy, 'o', color=c, ms=4.5)
        ax.text(3.3, yy, f'{o:.2f}', va='center', ha='right', fontsize=6.5, color=c)
        labels.append(lab)
        ticks.append(yy)
    tr = mtr.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(-0.03, y[2], 'Controls in the model', transform=tr, ha='right', va='center', fontsize=6.5, color='0.35',
            style='italic')
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xscale('log')
    ax.set_xlim(0.7, 3.4)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xticks([0.75, 1, 1.5, 2, 3])
    ax.set_xticklabels(['0.75', '1', '1.5', '2', '3'])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel('Odds ratio for being a phosphosite, per SD (95% CI)')
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.tight_layout(rect=(0, 0, 1, 0.96), w_pad=2.5)
    S.letters(fig, axes, 'ab')
    S.save(fig, 'fig_phosphosite')
    return {'proteins': Wp.UniProt_ID.nunique(), 'Ser/Thr': len(Wp), 'annotated': n_ann,
            'fifths_pct': f.round(1).tolist(), 'OR_alone': round(one['z_abs'][0], 2),
            'OR_full': {k: round(v[0], 2) for k, v in full.items()},
            'CI_full': {k: (round(v[1], 2), round(v[2], 2)) for k, v in full.items()}}


def lineage():
    order = ['Viruses', 'Bacteria', 'Protists', 'Fungi', 'Plants', 'Mammals']
    cell = {'Viruses': '', 'Bacteria': 'unicellular, ', 'Protists': 'unicellular, ', 'Fungi': 'mixed, ',
            'Plants': 'multicellular, ', 'Mammals': 'multicellular, '}
    K = pd.read_csv(f'{LIN}/kingdoms.csv', index_col=0)
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.9), gridspec_kw={'width_ratios': [1.4, 1]})
    ax = axes[0]
    y = np.arange(len(order))[::-1].astype(float)
    ax.axvline(0, color='black', lw=0.6)
    out, xs = {}, []
    sty = (('L', 'as simulated', 0.15, 'o', 'black', 'black'),
           ('L_adj', 'adjusted for $\\pi$ and\nAsp+Glu+Lys content', -0.15, 's', '0.45', 'white'))
    for col, name, dy, mk, c, fc in sty:
        m = -K.loc[order, f'{col}_mean'].to_numpy()
        lo, hi = -K.loc[order, f'{col}_hi'].to_numpy(), -K.loc[order, f'{col}_lo'].to_numpy()
        for yy, a, b in zip(y + dy, lo, hi):
            ax.plot([a, b], [yy, yy], color=c, lw=1.2, solid_capstyle='butt')
        ax.plot(m, y + dy, mk, color=c, mfc=fc, mew=1.0, ms=4, label=name)
        xs += list(lo) + list(hi)
        out[col] = dict(zip(order, np.round(m, 3)))
    ax.set_yticks(y)
    ax.set_yticklabels([f'{k}\n{cell[k]}{int(K.at[k, "organisms"])} organisms' for k in order], fontsize=6.5)
    pad = 0.05 * (max(xs) - min(xs))
    ax.set_xlim(min(xs) - pad, 0.5)                                    # room at the right for the legend
    ax.set_ylim(-0.6, len(order) - 0.15)
    n_c = int(K.chains.sum())
    ax.set_xlabel('Extra compaction on heating beyond\n$\\nu$ and length, organism mean (SD)')
    ax.set_title('CALVADOS-2, calibrated hydration', color=COL['CALVADOS-2'])
    ax.legend(loc='lower right', handletextpad=0.4, borderaxespad=0.1)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    out['organisms'] = K.loc[order, 'organisms'].astype(int).to_dict()
    out['chains'] = n_c
    ax = axes[1]
    sets = [('CALVADOS-2', f'{LIN}/per_protein.csv'), ('Mpipi-Recharged', f'{LIN}/recharged/per_protein.csv')]
    yy = np.arange(len(sets))[::-1].astype(float)
    ax.axvline(0, color='black', lw=0.6)

    def rn(v):
        return stats.norm.ppf((stats.rankdata(v) - 0.5) / len(v))

    def resid(yv, X):
        X1 = np.column_stack([np.ones(len(yv)), X])
        return yv - X1 @ np.linalg.lstsq(X1, yv, rcond=None)[0]
    xs = []
    for y0, (name, path) in zip(yy, sets):
        c = COL[name]
        p = pd.read_csv(path).dropna(subset=['age', 'L', 'L_adj', 'pi', 'dek'])
        n = len(p)
        r_raw = -stats.spearmanr(p.age, p.L)[0]
        ra = resid(rn(p.age.to_numpy()), np.column_stack([rn(p.pi.to_numpy()), rn(p.dek.to_numpy())]))
        r_par = -stats.pearsonr(ra, p.L_adj)[0]
        out[name] = (round(r_raw, 3), round(r_par, 3), n)
        for r, dy, mk, fill in ((r_raw, 0.15, 'o', c), (r_par, -0.15, 's', 'white')):
            z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
            lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
            ax.plot([lo, hi], [y0 + dy, y0 + dy], color=c, lw=1.2, solid_capstyle='butt')
            ax.plot(r, y0 + dy, mk, color=c, mfc=fill, mew=1.0, ms=4)
            xs += [lo, hi]
    ax.plot([], [], 'o', color=DARK, ms=4, label='raw (Spearman)')
    ax.plot([], [], 's', color=DARK, mfc='white', mew=1.0, ms=4, label='at fixed $\\pi$ and Asp+Glu+Lys\n(partial)')
    ax.set_yticks(yy)
    ax.set_yticklabels([f'{nm}\nn = {out[nm][2]:,}' for nm, _ in sets])
    for t, (nm, _) in zip(ax.get_yticklabels(), sets):
        t.set_color(COL[nm])
    ax.set_ylim(-0.75, 1.5)
    ax.set_xlim(min(-0.02, min(xs) - 0.03), max(0.02, max(xs) + 0.03))    # every interval inside, zero shown
    ax.set_xlabel('Correlation of residue-age score\nwith extra compaction on heating (95% CI)')
    ax.legend(loc='lower right', handletextpad=0.4, borderaxespad=0.1)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.tight_layout(rect=(0, 0, 1, 0.97), w_pad=2.0)
    S.letters(fig, axes, 'ab')
    S.save(fig, 'fig_lineage')
    return out


if __name__ == '__main__':
    which = sys.argv[1:] or ['compact_factors', 'phosphosite', 'lineage']
    for w in which:
        print(w, globals()[w](), flush=True)
