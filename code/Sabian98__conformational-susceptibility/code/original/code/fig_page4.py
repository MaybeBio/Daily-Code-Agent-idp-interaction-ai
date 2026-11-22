"""
fig_page4.py — manuscript Fig. 6, PAGE4 phosphoforms: reweighting predictions (page4_phosphoforms.py) against direct
simulations (page4_direct.py) and measurements (Kulkarni et al. PNAS 2017).
  a, b  smFRET efficiency A18:C63 and C63:P102, measured, reweighted (open) and direct (filled), 300 K, 150 mM
  c     change of Rg on phosphorylation: SAXS, smFRET-derived, reweighted and direct
  d     calibrated heating response of Rg at 300 K, reweighted and direct (FDT on each form's own 300 K run)
  e     direct heating curves (calibrated hydration) of WT and CLK2-PAGE4 relative to 300 K, with the FDT tangent at 300 K
  f     FDT response at 300 K against the slope of the direct heating curve, every simulated form, plain and calibrated
  g     PAGE4 sequence map: charge, regions, HIPK1 and CLK2 sites (marker size = occupancy), per-residue share of the
        calibrated heat response of WT (heat_crosses residues.csv, CALVADOS-2)
  h     change of the mean residue-residue distance on CLK2 phosphorylation, direct runs at 300 K (page4_maps.py)
Forms: WT, HIPK1 (T51, S9 at 0.35) and CLK2 (eight sites, occupancy-weighted); Mpipi-Recharged phosphate-like.
House style (2026-10-15, fig_extra.house): 6.5 in wide; figure text says "fluctuation-response relation" for the FDT.
    python fig_page4.py            (openmpipi env)
Output: BENDER-MT/manuscript/figs/fig_page4.{png,pdf}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import fig_extra as FX  # noqa: E402

ROOT = '/Volumes/User Homes/Nu Project/BENDER-MT'
P4 = f'{ROOT}/Traj Level Analysis/Analysis_FDT/total_response/page4'
OUT = f'{ROOT}/manuscript/figs'
COL = {'measured': '#222222', 'CALVADOS-2': '#2E6FA3', 'Mpipi-Recharged': '#6B4FA0'}   # as fig_total_response
DIRS = {'CALVADOS-2': 'direct_calvados', 'Mpipi-Recharged': 'direct_recharged'}
FORMS = ['WT', 'HIPK1', 'CLK2']
DIRECT = ['WT', 'HIPK1occ', 'CLK2occ']
ROWS = {'CALVADOS-2': ['WT', 'HIPK1 (occupancy)', 'CLK2 (occupancy)'],
        'Mpipi-Recharged': ['WT', 'HIPK1 (occupancy) (phosphate-like)', 'CLK2 (occupancy) (phosphate-like)']}
MEAS = {'E18_63': [0.55, 0.52, 0.35], 'E63_102': [0.65, 0.65, 0.58]}
SAXS, SAXS_SE, FRET_RG = [36.2, 34.7, 49.8], [1.1, 1.2, 1.9], [34, 35, 43]
OFF = {('CALVADOS-2', 'rw'): -0.18, ('CALVADOS-2', 'direct'): -0.08,
       ('Mpipi-Recharged', 'rw'): 0.08, ('Mpipi-Recharged', 'direct'): 0.18}


CAP = f'{ROOT}/heldout_temperature/capped_forcing/page4'   # capped main model (page4_capped_fdt.py, 2026-10-25)


def load():
    """Reweighted phosphoforms (capped, page4/phosphoforms.csv); direct runs: FRET, Rg and the plain response from
    page4_direct runs.csv, the calibrated total from page4_capped_fdt.py, and the calibrated heating curves simulated
    directly with the forcing in both force fields (page4cal_analysis.py; Mpipi-Recharged 2026-10-25, CALVADOS-2
    2026-11-02); the plain curves from page4_capped_fdt.py."""
    rw = pd.read_csv(f'{P4}/phosphoforms.csv')
    fd = pd.read_csv(f'{CAP}/fdt_300K.csv')
    cc = pd.read_csv(f'{CAP}/curves_capped.csv').replace({'curve': {'capped': 'calibrated'}})
    dr, cv = {}, {}
    for m, d in DIRS.items():
        r = pd.read_csv(f'{P4}/{d}/runs.csv')
        r = r[(r['T'] == 300) & (r.ionic == 0.15)].set_index('form')
        f = fd[fd.model == m].set_index('form')
        assert set(r.index) <= set(f.index), (m, set(r.index) - set(f.index))
        r['Rg_total'], r['Rg_total_se'] = f.loc[r.index, 'total_capped'], f.loc[r.index, 'total_capped_se']
        dr[m] = r
        cv[m] = cc[cc.model == m].drop(columns='model')
    # calibrated heating curves simulated directly with the forcing (page4cal_analysis.py) replace the reweighted ones,
    # which kept < 5% effective sample in Mpipi-Recharged (2026-10-25); CALVADOS-2 followed on 2026-11-02
    for m, d in (('CALVADOS-2', 'direct_cal_calvados'), ('Mpipi-Recharged', 'direct_cal_recharged')):
        dc = pd.read_csv(f'{CAP}/{d}/curves.csv')
        dc = dc[dc.obs == 'Rg'].rename(columns={f'se_{t}K': f'rel_se_{t}K' for t in (278, 288, 300, 310, 320)})
        dc = dc.assign(curve='calibrated', temps='278,288,300,310,320', min_ess=1.0)[
            ['form', 'curve', 'temps', 'slope300', 'slope300_se', 'min_ess'] + [c for c in dc.columns if c.startswith('rel_')]]
        cv[m] = pd.concat([cv[m][cv[m].curve != 'calibrated'], dc], ignore_index=True)
    return rw, dr, cv


def pct(m1, s1, m0, s0):
    return 100 * (m1 / m0 - 1), 100 * m1 / m0 * np.sqrt((s1 / m1) ** 2 + (s0 / m0) ** 2)


ESS_MIN = 0.1   # reweighted points at or below this Kish effective sample fraction are drawn pale (2026-10-20)


def mark(a, x, y, m, kind, yerr=None, ms=3.3, ess=None):
    filled = kind == 'direct'
    x, y = np.asarray(x, float), np.asarray(y, float)
    yerr = None if yerr is None else np.asarray(yerr, float)
    low = np.zeros(len(x), bool) if ess is None else np.asarray(ess) <= ESS_MIN
    for sel, alpha in ((~low, 1.0), (low, 0.3)):
        if sel.any():
            a.errorbar(x[sel], y[sel], yerr=None if yerr is None else yerr[sel], fmt='o' if filled else 's', ms=ms,
                       color=COL[m], mfc=COL[m] if filled else 'white', mew=0.8, elinewidth=0.7, capsize=0, zorder=3,
                       alpha=alpha)


def main():
    rw, dr, cv = load()
    FX.house()                                        # house style, Arial (figfont), save-time palette
    fig = plt.figure(figsize=(FX.WIDTH, 6.2))
    r1 = fig.add_gridspec(1, 3, left=0.1, right=0.985, top=0.955, bottom=0.735, wspace=0.55)
    r2 = fig.add_gridspec(1, 3, left=0.1, right=0.985, top=0.645, bottom=0.42, wspace=0.55)
    ax = [fig.add_subplot(r1[0, k]) for k in range(3)] + [fig.add_subplot(r2[0, 0])]
    ge = r2[0, 1].subgridspec(1, 2, wspace=0.1)
    ae = [fig.add_subplot(ge[0, 0])]
    ae.append(fig.add_subplot(ge[0, 1], sharey=ae[0]))
    ax += [ae[0], fig.add_subplot(r2[0, 2])]
    x = np.arange(3)
    # a, b  FRET efficiencies
    for p, (k, lab) in enumerate((('E18_63', 'FRET efficiency,\nresidues 18 and 63'),
                                  ('E63_102', 'FRET efficiency,\nresidues 63 and 102'))):
        a = ax[p]
        a.plot(x, MEAS[k], '-', color=COL['measured'], lw=0.9)
        a.plot(x, MEAS[k], 'D', color=COL['measured'], ms=3.6, zorder=4)
        for m in DIRS:
            r = rw[rw.model == m].set_index('pattern').loc[ROWS[m]]
            mark(a, x + OFF[(m, 'rw')], r[f'{k}_mean'].to_numpy(), m, 'rw', ess=r.ess.to_numpy())
            d = dr[m].loc[DIRECT]
            mark(a, x + OFF[(m, 'direct')], d[f'{k}_mean'].to_numpy(), m, 'direct', yerr=d[f'{k}_se'].to_numpy())
        a.set(xticks=x, xlim=(-0.45, 2.45), ylim=(0.3, 0.75), ylabel=lab)
        a.set_xticklabels(FORMS)
    h = [Line2D([], [], color=COL['measured'], marker='D', ms=3.6, lw=0.9, label='measured'),
         Line2D([], [], color='0.3', marker='s', mfc='white', ms=3.3, lw=0, label='reweighted'),
         Line2D([], [], color='0.3', marker='o', ms=3.3, lw=0, label='direct'),
         Line2D([], [], color=COL['CALVADOS-2'], lw=1.6, label='CALVADOS-2'),
         Line2D([], [], color=COL['Mpipi-Recharged'], lw=1.6, label='Mpipi-Recharged')]
    ax[1].legend(handles=h, loc='lower left', ncol=2, handlelength=1.2, borderaxespad=0.1, columnspacing=0.6)
    # c  change of Rg
    a = ax[2]
    rel = lambda v: 100 * (np.array(v) / v[0] - 1)  # noqa: E731
    xs = x[1:]
    a.errorbar(xs - 0.3, rel(SAXS)[1:], yerr=100 * np.array(SAXS_SE[1:]) / SAXS[0], fmt='D', color=COL['measured'],
               ms=3.6, elinewidth=0.7, label='SAXS')
    a.plot(xs - 0.22, rel(FRET_RG)[1:], 'D', color='0.6', ms=3.6, label='smFRET-derived')
    for m in DIRS:
        r = rw[rw.model == m].set_index('pattern').loc[ROWS[m][1:]]
        mark(a, xs + OFF[(m, 'rw')] + 0.04, r.Rg_pct.to_numpy(), m, 'rw', yerr=r.Rg_pct_se.to_numpy(), ess=r.ess.to_numpy())
        d = dr[m]
        v = [pct(d.loc[f, 'Rg_mean'], d.loc[f, 'Rg_se'], d.loc['WT', 'Rg_mean'], d.loc['WT', 'Rg_se']) for f in DIRECT[1:]]
        mark(a, xs + OFF[(m, 'direct')] + 0.04, [q[0] for q in v], m, 'direct', yerr=[q[1] for q in v])
    a.axhline(0, color='0.6', lw=0.5)
    a.set(xticks=xs, xlim=(0.5, 2.5), ylim=(-10, 58), ylabel=f'Change of {FX.RG}\nfrom WT (%)')
    a.set_xticklabels(FORMS[1:])
    a.legend(loc='upper left', handlelength=1, borderaxespad=0.1)
    # d  total response, reweighted and direct
    a = ax[3]
    for m in DIRS:
        r = rw[rw.model == m].set_index('pattern').loc[ROWS[m]]
        tot = r.Rg_total_wt.iloc[0] + r.Rg_dtotal.fillna(0).to_numpy()
        mark(a, x + OFF[(m, 'rw')], tot, m, 'rw', yerr=r.Rg_dtotal_se.fillna(0).to_numpy(), ess=r.ess.to_numpy())
        d = dr[m].loc[DIRECT]
        mark(a, x + OFF[(m, 'direct')], d.Rg_total.to_numpy(), m, 'direct', yerr=d.Rg_total_se.to_numpy())
    a.axhline(0, color='0.6', lw=0.5)
    a.set(xticks=x, xlim=(-0.45, 2.45), ylim=(-0.85, 0), ylabel=f'Total response of {FX.RG}\n(% per K)')
    a.set_xticklabels(FORMS)
    # e  calibrated heating curves of WT and CLK2, one small panel per force field (the fluctuation-response slopes are
    #    compared in f; tangents drawn here overlapped both curves, removed 2026-10-25)
    for a, m in zip(ae, DIRS):
        c = cv[m][cv[m].curve == 'calibrated'].set_index('form')
        for f, mk in (('WT', 'o'), ('CLK2occ', '^')):
            Ts = np.array([int(t) for t in c.loc[f, 'temps'].split(',')])
            y = np.array([c.loc[f, f'rel_{t}K'] for t in Ts])
            e = np.array([c.loc[f, f'rel_se_{t}K'] for t in Ts])
            weak = c.loc[f, 'min_ess'] < 0.05                 # capped reweighting kept < 5% effective sample: faded
            a.errorbar(Ts, y, yerr=e, fmt=mk + '-', ms=2.8, lw=0.9, color=COL[m], mfc=COL[m] if f == 'WT' else 'white',
                       elinewidth=0.7, zorder=3, alpha=0.35 if weak else 1.0)
        top = max(1.15, max(np.nanmax(np.array([c.loc[f, f'rel_{t}K'] + c.loc[f, f'rel_se_{t}K']
                                                for t in c.loc[f, 'temps'].split(',')])) for f in ('WT', 'CLK2occ')) + 0.01)
        a.set(xlim=(269, 331), xticks=[280, 300, 320], ylim=(0.69, top), yticks=[0.8, 0.9, 1.0, 1.1])   # every error bar inside
        a.set_title(m, fontsize=6.5, color=COL[m], pad=2)
    ae[0].set_ylabel(f'Mean {FX.RG} relative\nto 300 K')
    plt.setp(ae[1].get_yticklabels(), visible=False)
    ae[1].spines['left'].set_visible(False)
    ae[1].tick_params(axis='y', length=0)
    for a in ae:
        a.set_xlabel(' ')
    fig.canvas.draw()
    fig.text(0.5 * (ae[0].get_position().x0 + ae[1].get_position().x1), ae[0].xaxis.label.get_window_extent().y0 /
             fig.bbox.height, 'Temperature (K)', ha='center', va='bottom', fontsize=7.5)
    hh = [Line2D([], [], color='0.3', marker='o', ms=2.8, lw=0.9, label='WT'),
          Line2D([], [], color='0.3', marker='^', mfc='white', ms=2.8, lw=0.9, label='CLK2')]
    fig.legend(handles=hh, loc='lower left', bbox_to_anchor=(0, 0), bbox_transform=ae[0].transAxes, ncol=2,
               handlelength=1.4, borderaxespad=0.1, labelspacing=0.2, columnspacing=0.8)   # spans both e panels
    # f  fluctuation-response relation against the direct slope
    a = ax[5]
    for m in DIRS:
        d = dr[m]
        for lab, fcol, scol, fill in (('plain', 'Rg_plain', 'plain', False), ('calibrated', 'Rg_total', 'calibrated', True)):
            c = cv[m][cv[m].curve == scol].set_index('form')
            fs = [f for f in c.index if f in d.index]
            for f in fs:                                      # faded: capped curve kept < 5% effective sample
                a.errorbar(c.loc[f, 'slope300'], d.loc[f, fcol], xerr=c.loc[f, 'slope300_se'], yerr=d.loc[f, f'{fcol}_se'],
                           fmt='o', ms=3.0, color=COL[m], mfc=COL[m] if fill else 'white', mew=0.8, elinewidth=0.5,
                           capsize=0, alpha=0.35 if c.loc[f, 'min_ess'] < 0.05 else 1.0)
    lim = (-0.8, 0.22)
    a.plot(lim, lim, '--', color='0.6', lw=0.6, zorder=0)
    a.axhline(0, color='0.85', lw=0.5, zorder=0)
    a.axvline(0, color='0.85', lw=0.5, zorder=0)
    a.set(xlim=lim, ylim=lim, xlabel='Slope of direct curve (% per K)',
          ylabel='Fluctuation–response relation\nat 300 K (% per K)')
    a.text(-0.1, 0.13, 'plain', fontsize=6.5, color='0.35', va='center', ha='right')
    a.text(-0.42, -0.72, 'calibrated', fontsize=6.5, color='0.35')
    r3 = fig.add_gridspec(1, 2, width_ratios=[3.1, 1], left=0.1, right=0.93, top=0.305, bottom=0.07, wspace=0.25)
    ag = sequence_map(fig, r3[0, 0])
    ah = distance_change(fig, r3[0, 1])
    FX.add_letters(fig, [(ax[0], 'a', {'row': 1}), (ax[1], 'b', {'row': 1}), (ax[2], 'c', {'row': 1}),
                         (ax[3], 'd', {'row': 2}), (ae[0], 'e', {'row': 2}), (ax[5], 'f', {'row': 2}),
                         (ag[0], 'g', {'row': 3, 'xaxes': ag}), (ah, 'h', {'row': 3})])
    os.makedirs(OUT, exist_ok=True)
    FX.save(fig, f'{OUT}/fig_page4')


HIPK1_SITES = {51: 1.0, 9: 0.35}
CLK2_SITES = {7: 1.0, 9: 1.0, 51: 1.0, 85: 1.0, 73: 0.75, 94: 0.6, 71: 0.5, 79: 0.5}
REGIONS = [(4, 12, 'N-motif'), (43, 62, 'central acidic region'), (65, 73, 'helical'), (82, 95, 'C-motif')]
MAPS = f'{P4}/maps'


def sequence_map(fig, spec):
    """g: PAGE4 along its sequence: charge, regions (Lin et al. 2018), kinase sites with occupancy, and each residue's
    share of the calibrated heat response of WT (CALVADOS-2, 300 K). Returns the three axes."""
    r = pd.read_csv(f'{P4}/../heat_crosses/residues.csv')
    r = r[r.UniProt_ID == 'O60829'].sort_values('pos')
    seq = ''.join(r.residue)
    pos = r.pos.to_numpy()
    g = spec.subgridspec(3, 1, height_ratios=[0.75, 0.45, 1], hspace=0.1)
    a0 = fig.add_subplot(g[0])
    for lo, hi, lab in REGIONS:
        a0.axvspan(lo - 0.5, hi + 0.5, color='0.92', lw=0, zorder=0)
        a0.text((lo + hi) / 2, 1.5, lab.replace('central acidic region', 'central acidic\nregion'), ha='center',
                va='bottom', fontsize=6.5, color='0.3', linespacing=0.95)
    FX.charge_track(a0, pos, seq)
    a0.set_ylim(-1.3, 3.45)
    a1 = fig.add_subplot(g[1], sharex=a0)
    for y, (sites, col, lab) in enumerate(((HIPK1_SITES, '#C8792A', 'HIPK1'), (CLK2_SITES, '#2E8C6F', 'CLK2'))):
        for p_, occ in sites.items():
            a1.plot(p_, 1 - y, 'o', ms=2.6 + 2.6 * occ, color=col, mec='white', mew=0.4)
        a1.text(-2.5, 1 - y, lab, ha='right', va='center', fontsize=7.5, color=col)
    a1.set_ylim(-0.6, 1.6)
    a1.set_yticks([])
    plt.setp(a1.get_xticklabels(), visible=False)
    a1.tick_params(axis='x', length=0)
    for sp in ('left', 'bottom'):
        a1.spines[sp].set_visible(False)
    a2 = fig.add_subplot(g[2], sharex=a0)
    v = 100 * r.total_share.to_numpy()
    a2.bar(pos, v, width=0.85, color=np.where(v < 0, '0.62', '0.84'), lw=0)
    a2.axhline(0, color='0.5', lw=0.5)
    a2.set_ylabel('Share of WT\ntotal response\n(%)')
    a2.text(0.995, 0.06, 'CALVADOS-2', transform=a2.transAxes, ha='right', va='bottom', fontsize=6.5,
            color=COL['CALVADOS-2'])
    a2.set_xlabel('Residue')
    a2.set_xlim(0.3, len(seq) + 0.7)
    a2.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    return a0, a1, a2


def distance_change(fig, spec):
    """h: change of the mean residue-residue distance on CLK2 phosphorylation (CLK2 - WT, nm) from the direct runs at
    300 K: CALVADOS-2 below the diagonal (j > i), Mpipi-Recharged above."""
    ax = fig.add_subplot(spec)
    d = {}
    for m in ('calvados', 'recharged'):
        w = np.load(f'{MAPS}/{m}_WT_300K.npz')['mean_r']
        c = np.load(f'{MAPS}/{m}_CLK2occ_300K.npz')['mean_r']
        d[m] = c - w
    n = d['calvados'].shape[0]
    M = np.where(np.triu(np.ones((n, n), bool), 1), d['calvados'], np.where(np.tril(np.ones((n, n), bool), -1),
                                                                             d['recharged'], np.nan))
    v = np.nanpercentile(np.abs(M), 99)
    # softened diverging map (2026-10-24): RdBu_r without its darkest 15% at each end, blended 15% towards white, so the
    # mostly positive (expanded) map reads as soft red rather than maroon; same data and colour limits
    from matplotlib.colors import LinearSegmentedColormap
    base = plt.get_cmap('RdBu_r')(np.linspace(0.15, 0.85, 256))
    base[:, :3] = base[:, :3] + 0.15 * (1 - base[:, :3])
    soft = LinearSegmentedColormap.from_list('RdBu_r_soft', base)
    im = ax.imshow(M, origin='lower', cmap=soft, vmin=-v, vmax=v, extent=(0.5, n + 0.5, 0.5, n + 0.5))
    ax.text(0.03, 0.97, 'Mpipi-Recharged', transform=ax.transAxes, fontsize=6.5, va='top', ha='left', color='0.15',
            bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.0))
    ax.text(0.97, 0.03, 'CALVADOS-2', transform=ax.transAxes, fontsize=6.5, va='bottom', ha='right', color='0.15',
            bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.0))
    ax.set_xlabel('Residue')
    ax.set_ylabel('Residue')
    ax.set_xticks([1, 25, 50, 75, 100])
    ax.set_yticks([1, 25, 50, 75, 100])
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('CLK2 − WT distance (nm)', fontsize=6.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_visible(False)
    return ax


if __name__ == '__main__':
    main()
