"""
page4_fret_series.py — the smFRET temperature series that would test the PAGE4 heat x phosphorylation prediction:
<E> between residues 18-63 and 63-102 (and <Rg>) of WT, HIPK1- and CLK2-PAGE4 at 278-320 K from the direct runs
(page4_direct.py, 150 mM, pH 7.5), without (plain) and with the calibrated hydration term (each T reweighted with
exp(-sum_c gamma_c E_c(T) / RT), as the Wuttke calibration). Uses only the stored per-frame data (no trajectories read).
Per point: mean, 10-block jackknife SE, Kish effective sample fraction. Per form: linear slope over 278-320 K per 10 K
(jackknife leaving the same block out at every T) and the CLK2 - WT and HIPK1 - WT differences.
Figure: rows = force field, columns = E(18-63), E(63-102); calibrated solid, plain dashed; measured values at room
temperature in the smFRET buffer (Kulkarni et al. PNAS 2017; temperature not stated) as diamonds.
    BENDER_REANALYSIS_DIR=... python page4_fret_series.py       (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/page4/fret_series/{series.csv, slopes.csv, summary.txt}
        BENDER-MT/manuscript/figs/fig_page4_fret_series.{png,pdf}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import figfont  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import page4_direct as PD  # noqa: E402

OUT = os.path.join(PD.BM, 'Traj Level Analysis/Analysis_FDT/total_response/page4/fret_series')
FIG = os.path.join(PD.BM, 'manuscript/figs/fig_page4_fret_series')
TAG = {'calvados': '-I150-pH75', 'recharged': '-I150'}
FORMS = ['WT', 'HIPK1occ', 'CLK2occ']
LABEL = {'WT': 'WT', 'HIPK1occ': 'HIPK1-PAGE4', 'CLK2occ': 'CLK2-PAGE4'}
COL = {'WT': '#222222', 'HIPK1occ': '#C8792A', 'CLK2occ': '#2E8C6F'}
MEAS = {'E18_63': {'WT': 0.55, 'HIPK1occ': 0.52, 'CLK2occ': 0.35},
        'E63_102': {'WT': 0.65, 'HIPK1occ': 0.65, 'CLK2occ': 0.58}}
OBS = ['E18_63', 'E63_102', 'Rg']
TEMPS = (278, 288, 300, 310, 320)
NB = 10


def load(model, form):
    d = os.path.join(PD.out_dir(model), 'npz')
    return {T: dict(np.load(os.path.join(d, f'PAGE4-{form}{TAG[model]}_{T}K.npz'))) for T in TEMPS}


def wmean(z, obs, gam, T, keep=None):
    keep = np.ones(len(z['obs_Rg']), bool) if keep is None else keep
    return PD.W.reweight(z[f'obs_{obs}'][keep], {c: z[f'E_{c}'][keep] for c in gam}, gam, float(T))


def keep_mask(n, b):
    k = np.ones(n, bool)
    k[np.array_split(np.arange(n), NB)[b]] = False
    return k


def series(ser, obs, gam):
    """Means at every T and jackknife replicates (same block left out at every T)."""
    n = len(ser[TEMPS[0]]['obs_Rg'])
    full = np.array([wmean(ser[T], obs, gam, T)[0] for T in TEMPS])
    ess = np.array([wmean(ser[T], obs, gam, T)[1] for T in TEMPS])
    jk = np.array([[wmean(ser[T], obs, gam, T, keep_mask(n, b))[0] for T in TEMPS] for b in range(NB)])
    return full, ess, jk


def jk_se(rep):
    rep = np.asarray(rep)
    return np.sqrt((NB - 1) / NB * ((rep - rep.mean(0)) ** 2).sum(0))


def slope10(y):
    return 10 * np.polyfit(np.array(TEMPS, float), y, 1)[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, srows, store = [], [], {}
    for model in ('calvados', 'recharged'):
        gam = PD.MODELS[model]['gam']
        for lab, g in (('plain', {c: 0.0 for c in gam}), ('calibrated', gam)):
            for form in FORMS:
                ser = load(model, form)
                for obs in OBS:
                    full, ess, jk = series(ser, obs, g)
                    se = jk_se(jk)
                    sc = 10.0 if obs == 'Rg' else 1.0                      # Rg in Angstrom
                    for T, m, s, e in zip(TEMPS, full, se, ess):
                        rows.append({'model': PD.MODELS[model]['label'], 'hydration': lab, 'form': form, 'obs': obs,
                                     'T': T, 'mean': sc * m, 'se': sc * s, 'ess': e})
                    sl = slope10(sc * full)
                    rep = np.array([slope10(sc * r) for r in jk])
                    store[(model, lab, form, obs)] = (sl, rep)
                    srows.append({'model': PD.MODELS[model]['label'], 'hydration': lab, 'form': form, 'obs': obs,
                                  'per10K': sl, 'per10K_se': float(jk_se(rep[:, None])[0])})
    df, sl = pd.DataFrame(rows), pd.DataFrame(srows)
    for form in ('HIPK1occ', 'CLK2occ'):                                   # differences from WT (independent runs)
        for (model, lab, f, obs), (s, rep) in list(store.items()):
            if f != form:
                continue
            s0, rep0 = store[(model, lab, 'WT', obs)]
            se = float(np.hypot(jk_se(rep[:, None])[0], jk_se(rep0[:, None])[0]))
            sl.loc[len(sl)] = {'model': PD.MODELS[model]['label'], 'hydration': lab, 'form': f'{form}-WT', 'obs': obs,
                               'per10K': s - s0, 'per10K_se': se}
    df.to_csv(os.path.join(OUT, 'series.csv'), index=False)
    sl.to_csv(os.path.join(OUT, 'slopes.csv'), index=False)

    Lns = ['PAGE4 smFRET temperature series predicted from the direct runs (150 mM, pH 7.5), 278-320 K',
           'mean +- 10-block jackknife SE [Kish effective sample fraction]; slopes per 10 K from a linear fit', '']
    for (model, lab), g in df.groupby(['model', 'hydration'], sort=False):
        Lns.append(f'{model}, {lab}')
        for obs in OBS:
            for form in FORMS:
                h = g[(g.obs == obs) & (g.form == form)]
                pts = '  '.join(f"{int(T)}K {m:.3f}+-{s:.3f}[{e:.2f}]" for T, m, s, e in zip(h['T'], h['mean'], h.se, h.ess))
                Lns.append(f'  {obs:8s} {form:9s} {pts}')
        h = sl[(sl.model == model) & (sl.hydration == lab)]
        for _, r in h.iterrows():
            Lns.append(f"  slope per 10 K {r.obs:8s} {r.form:13s} {r.per10K:+.4f} +- {r.per10K_se:.4f}")
        Lns.append('')
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(Lns) + '\n')
    print('\n'.join(Lns))
    figure(df)


def figure(df):
    plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                         'pdf.fonttype': 42})
    figfont.apply()                                   # Arial for all figures
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True)
    for i, model in enumerate(['CALVADOS-2', 'Mpipi-Recharged']):
        for j, obs in enumerate(['E18_63', 'E63_102']):
            a = ax[i, j]
            for form in FORMS:
                for lab, ls, alpha in (('calibrated', '-', 1.0), ('plain', '--', 0.55)):
                    h = df[(df.model == model) & (df.hydration == lab) & (df.form == form) & (df.obs == obs)]
                    a.errorbar(h['T'], h['mean'], yerr=h.se, fmt='o' if lab == 'calibrated' else 's', ls=ls, ms=3.2,
                               color=COL[form], alpha=alpha, mfc=COL[form] if lab == 'calibrated' else 'white',
                               mew=0.8, lw=1.0, elinewidth=0.7, capsize=0)
                a.plot(295, MEAS[obs][form], 'D', ms=5, color=COL[form], mfc='white', mew=1.1, zorder=4)
            a.set_ylabel(f"smFRET $E$, residues {obs[1:].replace('_', '–')}")
            a.set_title(model, fontsize=8, loc='left')
            a.spines[['top', 'right']].set_visible(False)
            if i == 1:
                a.set_xlabel('Temperature (K)')
            a.text(-0.17, 1.04, 'abcd'[2 * i + j], transform=a.transAxes, fontsize=10, fontweight='bold')
    hs = [plt.Line2D([], [], color=COL[f], marker='o', lw=1.0, ms=3.2, label=LABEL[f]) for f in FORMS]
    hs += [plt.Line2D([], [], color='0.4', ls='-', lw=1.0, label='calibrated hydration'),
           plt.Line2D([], [], color='0.4', ls='--', lw=1.0, alpha=0.55, label='no hydration term'),
           plt.Line2D([], [], color='0.2', marker='D', ls='none', mfc='white', ms=5, label='measured, room temperature')]
    fig.legend(handles=hs, loc='lower center', ncol=3, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}.{ext}', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
