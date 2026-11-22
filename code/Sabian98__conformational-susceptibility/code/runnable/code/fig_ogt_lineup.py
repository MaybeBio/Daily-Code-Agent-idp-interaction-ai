"""
fig_ogt_lineup.py — Supplementary figure: disordered proteins of organisms that grow from 9 to 75 C. Heat sensitivity
follows the compaction of each organism's chains, not its growth temperature (CALVADOS-2 with calibrated hydration).

  a  ten organisms on a growth-temperature scale (PhyloPic silhouettes, all CC0; a relative where the species has none,
     see manuscript/figs/icons_phylopic/CREDITS.csv), with the number of chains
  b  total response of every chain of each organism, with the organism median, the median apparent scaling exponent nu
     and the proteome-wide interquartile band of all chains with a known growth temperature
  c  all 757 organisms: median total response against median nu, coloured by optimal growth temperature
  d  the same organisms: total response at fixed compaction (rank residual on nu) against growth temperature
  e  within bacterial phyla (phylum fixed effects, ogt_followup.py): standardized slopes against growth temperature of
     charge patterning, Asp + Glu + Lys, pi residues and Arg, which shift, and of nu, which does not (95% bootstrap CIs)

Reads Thermal adaptation/bender_ogt/chains.csv (ogt_response.py) only.
    python fig_ogt_lineup.py          Output: manuscript/figs/fig_ogt_lineup.{pdf,png}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402

CHAINS = f'{FX.ROOT}/Thermal adaptation/bender_ogt/chains.csv'
ICONS = f'{FX.ROOT}/manuscript/figs/icons_phylopic'
OUT = f'{FX.ROOT}/manuscript/figs/fig_ogt_lineup'
FOLLOWUP = f'{FX.ROOT}/Thermal adaptation/bender_ogt/followup_nu_filtered/followup_composition.csv'   # ogt_followup.py
LINEUP = [  # prefix of the UniProt organism name, label, icon file
    ('Colwellia psychrerythraea', 'Colwellia', 'colwellia'),
    ('Skeletonema marinoi', 'Skeletonema', 'skeletonema'),
    ('Emiliania huxleyi', 'Emiliania', 'emiliania'),
    ('Chara braunii', 'Chara', 'chara'),
    ('Neurospora crassa', 'Neurospora', 'neurospora'),
    ('Saccharomyces cerevisiae', 'Saccharomyces', 'saccharomyces'),
    ('Plasmodium falciparum', 'Plasmodium', 'plasmodium'),
    ('Geobacillus kaustophilus', 'Geobacillus', 'geobacillus'),
    ('Thermus thermophilus', 'Thermus', 'thermus'),
    ('Thermotoga maritima', 'Thermotoga', 'thermotoga'),
]
COOL, WARM = '#2E6FA3', '#B5452B'
NU_MIN = 0.30
NU = FX.NU


def tcmap():
    return LinearSegmentedColormap.from_list('ogt', [FX.final_colour(COOL), '#D9D4CC', FX.final_colour(WARM)])


def icon(path, ink=FX.INK):
    im = plt.imread(path).astype(float)
    if im.ndim == 2:
        g, a = im, np.ones_like(im)
    elif im.shape[2] == 2:
        g, a = im[..., 0], im[..., 1]
    else:
        g, a = im[..., :3].mean(axis=2), im[..., 3] if im.shape[2] == 4 else np.ones(im.shape[:2])
    from matplotlib.colors import to_rgb
    c = np.array(to_rgb(ink))
    rgb = c[None, None, :] * (1 - g[..., None]) + g[..., None]       # black -> ink, white stays white
    return np.dstack([rgb, a])


def data():
    ch = pd.read_csv(CHAINS)
    c = ch[ch.model == 'CALVADOS-2'].copy()
    ok = c[c.nu > NU_MIN]          # scaling fits at or below 0.30 failed (Supplementary Section S1)
    org = ok.groupby('organism').agg(OGT=('OGT', 'first'), total=('total', 'median'), nu=('nu', 'median'),
                                     n=('total', 'size'), kingdom=('kingdom', 'first'))
    lin = []
    for pref, lab, ic in LINEUP:
        s = c[c.organism.str.startswith(pref)]
        assert len(s), pref
        lin.append(dict(label=lab, icon=ic, OGT=float(s.OGT.median()), n=len(s), total=s.total.values,
                        med=float(s.total.median()), nu=float(s.nu[s.nu > NU_MIN].median()),
                        nuv=s.nu[s.nu > NU_MIN].values))
    return c, org, lin


def resid_rank(y, x):
    ry, rx = rankdata(y), rankdata(x)
    b = np.polyfit(rx, ry, 1)
    return ry - np.polyval(b, rx)


def plot(c, org, lin):
    FX.house()
    cmap, norm = tcmap(), Normalize(0, 80)
    # vertical layout in inches (2026-10-19: 6.1 -> 5.32 in tall so that figure and caption fit one page; text sizes kept)
    H = 5.32
    Y = {'cde': (0.46, 1.25), 'b2': (2.15, 0.55), 'b': (2.79, 0.75), 'a': (3.82, 1.40)}     # bottom, height
    fy = lambda k: (Y[k][0] / H, Y[k][1] / H)                                                 # noqa: E731
    fig = plt.figure(figsize=(FX.WIDTH, H))
    nL = len(lin)
    left, right = 0.085, 0.985
    xs = np.linspace(0, 1, nL)                      # slot centres in axes fraction
    # a  icons and growth-temperature scale (y in inches from the bottom of the axes)
    a = fig.add_axes([left, fy('a')[0], right - left, fy('a')[1]])
    a.set_xlim(-0.06, 1.06)
    a.set_ylim(0, Y['a'][1])
    a.set_axis_off()
    for x, L in zip(xs, lin):
        im = icon(os.path.join(ICONS, L['icon'] + '.png'))
        h, w = im.shape[:2]
        sz = 0.33                                    # icon box height (in)
        ax_w = (right - left) * fig.get_figwidth()
        bw = sz / ax_w * w / h                      # keep the aspect ratio
        bw, bh = (0.075, sz * 0.075 / bw) if bw > 0.075 else (bw, sz)
        a.imshow(im, extent=(x - bw / 2, x + bw / 2, 1.02, 1.02 + bh), aspect='auto', zorder=3, interpolation='antialiased')
        a.text(x, 0.93, L['label'], ha='center', va='top', fontsize=6.3, style='italic')
        a.text(x, 0.80, f"{L['OGT']:.0f} °C\n{L['n']} chains", ha='center', va='top', fontsize=6.0, color='0.35',
               linespacing=1.15)
        tx = -0.04 + 1.08 * L['OGT'] / 80          # position on the temperature bar
        a.plot([x, tx], [0.56, 0.395], color=cmap(norm(L['OGT'])), lw=0.8, zorder=1)
    # temperature bar
    grad = np.linspace(0, 80, 400)[None, :]
    a.imshow(grad, extent=(-0.04, 1.04, 0.198, 0.362), cmap=cmap, norm=norm, aspect='auto', zorder=1)
    for t in range(0, 81, 10):
        tx = -0.04 + 1.08 * t / 80
        a.plot([tx, tx], [0.165, 0.198], color='0.3', lw=0.6)
        a.text(tx, 0.132, f'{t}', ha='center', va='top', fontsize=6.3)
    a.text(0.5, 0.0, 'Optimal growth temperature (°C)', ha='center', va='top', fontsize=7.5)
    # b  medians and interquartile ranges per organism: total response (top) and compaction (bottom)
    b = fig.add_axes([left, fy('b')[0], right - left, fy('b')[1]])
    b2 = fig.add_axes([left, fy('b2')[0], right - left, fy('b2')[1]], sharex=b)
    ok_all = c[c.nu > NU_MIN]
    for ax_, key, allv in ((b, 'total', c.total), (b2, 'nuv', ok_all.nu)):
        q25, q50, q75 = np.percentile(allv, [25, 50, 75])
        ax_.axhspan(q25, q75, color='0.92', zorder=0)
        ax_.axhline(q50, color='0.6', lw=0.7, ls='--', zorder=1)
        for x, L in zip(xs, lin):
            v = L[key]
            lo_, md, hi_ = np.percentile(v, [25, 50, 75])
            ax_.plot([x, x], [lo_, hi_], color=FX.INK, lw=1.0, zorder=2, solid_capstyle='round')
            ax_.scatter([x], [md], s=24, color=cmap(norm(L['OGT'])), edgecolor=FX.INK, lw=0.6, zorder=3)
        ax_.set_xlim(-0.06, 1.06)
        ax_.tick_params(axis='x', bottom=False, labelbottom=False)
        ax_.spines['bottom'].set_visible(False)
    b.set_ylabel('Total response\n(% per K)')
    b2.set_ylabel('Scaling\nexponent ' + NU)
    b.set_ylim(-1.05, 0.0)
    # c  organism medians, total against nu
    cc = fig.add_axes([left, fy('cde')[0], 0.235, fy('cde')[1]])
    sc = cc.scatter(org.nu, org.total, c=org.OGT, cmap=cmap, norm=norm, s=6, lw=0, alpha=0.85)
    for L in lin:
        cc.scatter(L['nu'], L['med'], s=22, facecolor='none', edgecolor=FX.INK, lw=0.7, zorder=3)
    rho_c = spearmanr(org.total, org.nu)[0]
    cc.text(0.03, 0.05, r'$\rho$' + FX.minus(f' = {rho_c:.2f}'), transform=cc.transAxes, fontsize=6.5)
    cc.set(xlabel='Median ' + NU + ' of the organism', ylabel='Median total response\n(% per K)')
    # d  compaction-adjusted response against growth temperature
    d = fig.add_axes([0.415, fy('cde')[0], 0.215, fy('cde')[1]])
    res = resid_rank(org.total.values, org.nu.values)
    d.scatter(org.OGT, res, c=org.OGT, cmap=cmap, norm=norm, s=6, lw=0, alpha=0.85)
    d.axhline(0, color='0.6', lw=0.7, ls='--')
    rho_d = np.corrcoef(res, resid_rank(org.OGT.values, org.nu.values))[0, 1]
    rho_raw = spearmanr(org.total, org.OGT)[0]
    d.text(0.97, 0.04, r'partial $\rho$' + FX.minus(f' = {rho_d:.2f}') + '\n' + r'$\rho$' + FX.minus(f' = {rho_raw:.2f}') +
           ' without\nfixing ' + NU, transform=d.transAxes, ha='right', fontsize=6.0)
    d.set(xlabel='Optimal growth temperature (°C)', ylabel='Total response at fixed\ncompaction (rank residual)')
    # e  within bacterial phyla: composition shifts with growth temperature, nu does not
    e = fig.add_axes([0.79, fy('cde')[0], right - 0.79, fy('cde')[1]])
    F = pd.read_csv(FOLLOWUP)
    F = F[F.analysis.str.startswith('descriptor')]
    rows = [('kappa', 'CALVADOS-2', r'charge patterning $\kappa$', FX.INK), ('f_DEK', 'CALVADOS-2', 'Asp + Glu + Lys', FX.INK),
            ('f_pi', 'CALVADOS-2', r'$\pi$ residues', FX.INK), ('f_R', 'CALVADOS-2', 'Arg', FX.INK),
            ('nu', 'CALVADOS-2', NU + ', CALVADOS-2', FX.COL['CALVADOS-2']),
            ('nu', 'Mpipi-Recharged', NU + ', Mpipi-Rech.', FX.COL['Mpipi-Recharged'])]
    yy = np.arange(len(rows))[::-1]
    for y_, (term, model, lab_, col) in zip(yy, rows):
        r_ = F[(F.term == term) & (F.model == model)].iloc[0]
        e.plot([r_.ci_lo, r_.ci_hi], [y_, y_], color=col, lw=1.1, solid_capstyle='round')
        e.scatter([r_.slope_per_10C], [y_], s=20, color=col, zorder=3)
    e.axvline(0, color='0.6', lw=0.7, ls='--')
    e.axhline(1.5, color='0.85', lw=0.6)
    e.set_yticks(yy, [r[2] for r in rows], fontsize=6.3)
    e.set_xlabel('Slope against growth\ntemperature (SD per 10 °C)')
    e.set_title('Within bacterial phyla', fontsize=7.0, pad=3)
    FX.add_letters(fig, [(a, 'a', {'x': 0.012, 'y': 1 - 0.06 / H}), (b, 'b', {'x': 0.012}), (cc, 'c', {'row': 2}),
                         (d, 'd', {'row': 2}), (e, 'e', {'row': 2, 'xaxes': [e]})])
    FX.save(fig, OUT)
    return dict(rho_total_nu=rho_c, partial=rho_d, rho_total_ogt=rho_raw, n_org=len(org),
                q=(q25, q50, q75), n_chains=len(c))


if __name__ == '__main__':
    c, org, lin = data()
    out = plot(c, org, lin)
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in out.items()})
    for L in lin:
        print(f"{L['label']:14s} OGT {L['OGT']:5.1f}  n {L['n']:3d}  median {L['med']:+.2f}  nu {L['nu']:.3f}")
