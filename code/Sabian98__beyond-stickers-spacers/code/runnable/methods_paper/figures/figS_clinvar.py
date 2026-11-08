"""
figS_clinvar.py — Supplementary figure for app:clinvar, in the style of IDRome Fig. 5e-g: normalised distributions of
predicted hub properties for pathogenic and benign ClinVar missense variants in human IDRs, with Cohen's d (pathogenic
minus benign, bootstrap SE), Brunner-Munzel degrees of freedom and two-sided P for the comparison outside collagens.
Rows: IDRome's mapped variants, newer variants (ClinVar 5 Oct 2026). Columns: within-IDR percentile of the predicted
occupancy of the wild-type residue, predicted change of occupancy at the variant residue.
Benign and pathogenic outside collagens are solid lines; pathogenic collagen variants (94-97% glycine substitutions in
Gly-X-Y repeats) are dashed. Okabe-Ito colours.
  ~/.venvs/bender/bin/python figS_clinvar.py  ->  figS_clinvar.pdf / .png
"""
import csv
import gzip
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import brunnermunzel

import figstyle as fs

sys.path.insert(0, fs.EV)
import clinvar_hub as C  # noqa: E402
import clinvar_hub_idrome_style as I  # noqa: E402

BENIGN, PATH, COLL = '#0072B2', '#D55E00', '#CC79A7'


def load():
    cls = {}
    rd = csv.reader(gzip.open(os.path.join(C.EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(rd)
    ix = {c: i for i, c in enumerate(h)}
    for row in rd:
        s = row[ix['significance']]
        cls[row[0]] = 1 if 'athogenic' in s else 0 if 'enign' in s else None
    rep = {}
    with gzip.open(os.path.join(C.EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
        r2 = csv.reader(fh, delimiter='\t')
        hh = next(r2)
        jx = {c.lstrip('#'): i for i, c in enumerate(hh)}
        for row in r2:
            if row[jx['Assembly']] == 'GRCh38':
                cs = row[jx['ClinicalSignificance']]
                if cs in C.PATH_RC or cs in C.BEN_RC:
                    rep[row[jx['VariationID']]] = 1 if cs in C.PATH_RC else 0
    wt = {}
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'wt_profiles.csv'))):
        wt.setdefault(rr['idr'], []).append(float(rr['pred_consensus']))
    pct = {k: np.array(v).argsort().argsort() / max(1, len(v) - 1) for k, v in wt.items()}
    rows = []
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'variants.csv'))):
        y = cls.get(rr['variant']) if rr['set'] == 'primary' else rep.get(rr['variation_id'])
        if y is None:
            continue
        rows.append(dict(set=rr['set'], coll=rr['gene'].startswith('COL'), y=y, pct=pct[rr['idr']][int(rr['k'])],
                         dsite=float(rr['delta_site'])))
    return rows


def main():
    rows = load()
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.6))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.94, bottom=0.11, hspace=0.75, wspace=0.28)
    # (key, label, histogram bins over the full range, visible window); the change uses fine bins over a wide range so
    # that values outside the window are not piled into the end bins
    cols = [('pct', 'within-IDR percentile of predicted occupancy', np.linspace(0, 1, 21), (0, 1)),
            ('dsite', 'predicted change of occupancy at the variant residue', np.arange(-0.3, 0.3001, 0.002), (-0.04, 0.04))]
    out = []
    for i, (s, sname) in enumerate((('primary', 'IDRome set'), ('replication', 'newer ClinVar variants'))):
        S = [r for r in rows if r['set'] == s]
        for j, (key, xlabel, bins, win) in enumerate(cols):
            ax = axes[i, j]
            ben = np.array([r[key] for r in S if not r['y'] and not r['coll']])
            pat = np.array([r[key] for r in S if r['y'] and not r['coll']])
            col = np.array([r[key] for r in S if r['y'] and r['coll']])
            for v, c, ls, lab in ((ben, BENIGN, '-', f'benign ({len(ben):,})'),
                                  (pat, PATH, '-', f'pathogenic ({len(pat):,})'),
                                  (col, COLL, (0, (3, 1.5)), f'pathogenic, collagens ({len(col):,})')):
                ax.hist(v, bins=bins, density=True, histtype='step', lw=1.2, color=c,
                        linestyle=ls, label=lab)
            bm = brunnermunzel(pat, ben)
            d, se = I.cohen(pat, ben, rng)
            dof = I.bm_df(pat, ben)
            ptxt = f'P = {bm.pvalue:.0e}'.replace('e-', ' × 10$^{-').replace('$^{-0', '$^{-') + '}$'
            ax.text(0.02, 0.97, f"outside collagens\nCohen's d = {d:+.2f} ± {se:.2f}\ndegrees of freedom {dof:,.0f}\n{ptxt}".replace('-', '\u2212', 1) if d < 0 else f"outside collagens\nCohen's d = {d:+.2f} ± {se:.2f}\ndegrees of freedom {dof:,.0f}\n{ptxt}",
                    transform=ax.transAxes, ha='left', va='top', fontsize=6, color='0.15')
            ax.set_xlabel(xlabel)
            ax.set_yticks([])
            if j == 0:
                ax.set_ylabel('normalized distribution')
            ax.set_xlim(*win)
            ax.set_ylim(0, ax.get_ylim()[1] * 1.55)
            if key == 'dsite':
                ax.axvline(0, color='0.75', lw=0.6, zorder=0)
            fs.panel(ax, "ab"[j] if i == 0 else "cd"[j], x=-0.08 if j == 0 else -0.02, y=1.04, title=sname)
            ax.legend(frameon=False, loc='upper right', fontsize=5.8, handlelength=1.8)
            out.append(f'{sname:24s} {key:6s} d {d:+.3f} +- {se:.3f}  DoF {dof:.0f}  P {bm.pvalue:.2e}  '
                       f'n benign {len(ben)}  pathogenic {len(pat)}  collagen pathogenic {len(col)}')
    fs.save(fig, 'figS_clinvar')
    open(os.path.join(fs.HERE, 'figS_clinvar_values.txt'), 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
