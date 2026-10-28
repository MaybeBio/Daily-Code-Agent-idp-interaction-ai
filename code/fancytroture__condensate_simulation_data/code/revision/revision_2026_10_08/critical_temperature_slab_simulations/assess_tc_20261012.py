"""Read-only density extraction and conditional critical-scaling assessment."""
import csv
import importlib.util
import json
from pathlib import Path
import sys

import check_phase_gpu12 as phase
import numpy as np
import matplotlib.pyplot as plt

# The simulation directory has queue.py, which must not shadow the stdlib.
sys.path = [p for p in sys.path if Path(p or '.').resolve() != Path(__file__).parent.resolve()]
ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'revision_2026_09_30/fus_critical_temperature'
sys.path.insert(0, str(OLD))
spec = importlib.util.spec_from_file_location('old_tc_analysis', OLD / 'analyze.py')
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)
OUT = phase.HERE / 'tc_assessment_20261012'


def readcsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    OUT.mkdir(exist_ok=True)
    phase.OUT = OUT
    cases = ['FUS_335K', 'A1_267K', 'TDP43_315K']
    for case in cases:
        if not (OUT / (case + '.npz')).exists():
            result = phase.analyze(case)
            phase.write_json(OUT / (case + '_snapshot.json'), result[0])
            phase.table(OUT / (case + '_blocks.csv'), result[1])
    original = readcsv(ROOT / 'revision_2026_10_08/critical_temperature_two_point/estimates.csv')
    original = {r['system']: r for r in original
                if r['source'] == 'audited_short_window' and r['core_halfwidth_nm'] == '5'}
    profiles_old = readcsv(ROOT / 'revision_2026_09_30/selection_density/density_profiles.csv')
    fits, windows, points = [], [], []
    plt.rcParams.update({'font.size': 12, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(3, 3, figsize=(16, 11), layout='constrained')
    for row, case in enumerate(cases):
        system, temperature = case.split('_')
        temperature = int(temperature[:-1])
        data = np.load(OUT / (case + '.npz'))
        t, z, rho = data['time_ns'], data['z_centered_nm'], data['centered_g_ml']
        old = original[system]
        temps = [float(old['T_L_K']), float(old['T_H_K']), temperature]
        gaps = [float(old['dense_L_g_ml']) - float(old['vapor_L_g_ml']),
                float(old['dense_H_g_ml']) - float(old['vapor_H_g_ml'])]
        primary = None
        for window, start in [('post5us', 5000), ('post3us', 3000),
                              ('last2us', t[-1] - 2000), ('last1us', t[-1] - 1000)]:
            select = t > start
            if select.sum() < 500:
                continue
            profile = rho[select].mean(0)
            pf = fitmod.fit_profile(z, profile)
            series = rho[select][:, abs(z) < 5].mean(1)
            half = len(series) // 2
            drift = 100 * (series[half:].mean() - series[:half].mean()) / series.mean()
            dense = float(profile[abs(z) < 5].mean())
            vapor = float(profile[abs(z) > 80].mean())
            windows.append(dict(case=case, window=window, start_ns=start, end_ns=float(t[-1]),
                                dense_g_ml=dense, vapor_g_ml=vapor,
                                half_window_drift_pct=float(drift), **pf))
            for core in [3, 5, 7]:
                previous = []
                for temp in temps[:2]:
                    rows = [r for r in profiles_old if r['case'] == f'{system}_{int(temp)}']
                    zz = np.array([float(r['z_nm']) for r in rows])
                    yy = np.array([float(r['density_g_ml']) for r in rows])
                    previous.append(yy[abs(zz) < core].mean() - yy[abs(zz) > 80].mean())
                contrast = previous + [float(profile[abs(z) < core].mean() - vapor)]
                for beta in [.325, .30, .35]:
                    fitted = fitmod.critical_fit(temps, contrast, beta)
                    fits.append(dict(system=system, method='three_point', window=window,
                                     core_nm=core, **fitted,
                                     TH_over_Tc=temps[1]/fitted['Tc_K']))
                    if window == 'post5us' and core == 5 and beta == .325:
                        primary = fitted
                        for temp, gap in zip(temps, contrast):
                            points.append(dict(system=system, temperature_K=temp, gap_g_ml=float(gap),
                                               source='new_post5us' if temp == temperature else 'old_audited_short_window'))
                        axes[row, 2].plot(temps, contrast, 'o', color='#333333', label='Density gaps')
                        xx = np.linspace(min(temps), fitted['Tc_K'], 200)
                        axes[row, 2].plot(xx, fitted['amplitude'] * (fitted['Tc_K']-xx)**.325,
                                          color='#2f7f73', label=f"3 points: {fitted['Tc_K']:.1f} K")
                        prediction = gaps[0] * ((float(old['Tc_initial_K']) - temperature) /
                                                (float(old['Tc_initial_K']) - temps[0]))**.325
                        axes[row, 2].plot(temperature, prediction, 'x', color='#b54b5f',
                                          label='Previous 2-point prediction')
            if window in ['post5us', 'last1us']:
                axes[row, 0].plot(z, profile, label=f'{start/1000:g}-{t[-1]/1000:g} us')
        if system == 'FUS':
            oldrows = readcsv(OLD / 'coexistence_densities.csv')
            for lower in [280, 295]:
                retained = [r for r in oldrows if r['window'] == 'post5us'
                            and lower <= float(r['temperature_K']) <= 325]
                tt = [float(r['temperature_K']) for r in retained] + [temperature]
                yy = [float(r['gap_g_ml']) for r in retained] + [points[-1]['gap_g_ml']]
                fitted = fitmod.critical_fit(tt, yy)
                fits.append(dict(system=system, method=f'multiple_T_from_{lower}', window='post5us',
                                 core_nm=5, **fitted, TH_over_Tc=325/fitted['Tc_K']))
                if lower == 295:
                    xx = np.linspace(min(tt), fitted['Tc_K'], 200)
                    axes[row, 2].plot(tt[:-1], yy[:-1], '+', color='#427cac')
                    axes[row, 2].plot(xx, fitted['amplitude']*(fitted['Tc_K']-xx)**.325,
                                      '--', color='#427cac', label=f"6 points: {fitted['Tc_K']:.1f} K")
        blockrows = readcsv(OUT / (case + '_blocks.csv'))
        for key, label in [('core5_g_ml', 'Core'), ('outer80_g_ml', 'Dilute')]:
            axes[row, 1].plot([float(b['end_ns'])/1000 for b in blockrows],
                              [float(b[key]) for b in blockrows], label=label)
        axes[row, 0].set(xlabel='Centered z (nm)', ylabel='Density (g/mL)', title=case,
                         xlim=(-140, 140), ylim=(0, None))
        axes[row, 1].set(xlabel='Time (us)', ylabel='Density (g/mL)', title='100 ns block means', ylim=(0, None))
        axes[row, 2].set(xlabel='Temperature (K)', ylabel='Dense - dilute density (g/mL)',
                         title='Conditional critical-scaling fit', ylim=(0, None))
        for ax in axes[row]:
            ax.legend(frameon=False, fontsize=9)
        if system == 'TDP43':
            axes[row, 2].set_title('Provisional fit: core-plateau check failed')
        print(case, 'completed ns', t[-1], 'primary fit', primary, flush=True)
    for fitted in fits:
        diagnostic = next(w for w in windows if w['case'].startswith(fitted['system'] + '_')
                          and w['window'] == fitted['window'])
        fitted['new_point_core_plateau'] = diagnostic['core_plateau']
        fitted['new_point_vapor_plateau'] = diagnostic['vapor_plateau80']
        fitted['status'] = ('conditional_finite_box_estimate' if diagnostic['core_plateau']
                            else 'diagnostic_only_core_plateau_check_failed')
    phase.table(OUT / 'window_diagnostics.csv', windows)
    phase.table(OUT / 'fit_sensitivity.csv', fits)
    phase.table(OUT / 'primary_density_gaps.csv', points)
    for suffix in ['png', 'pdf']:
        fig.savefig(OUT / f'critical_temperature_assessment.{suffix}', dpi=180)
    plt.close(fig)
    phase.write_json(OUT / 'assessment_manifest.json', dict(
        script_sha256=phase.sha(__file__), extraction_script_sha256=phase.sha(phase.__file__),
        fit_script_sha256=phase.sha(OLD / 'analyze.py'), simulations_modified=False,
        scope='Conditional finite-box extrapolation; no claim of thermodynamic-limit criticality.',
        caveats=['Original A1/TDP densities use an audited short window.',
                 'TDP43 includes completed segments only; simulation is ongoing.',
                 'TDP43 fails the preset core-plateau geometric check; its Tc fit is diagnostic only.',
                 'A1 final 1 us core density has a -16.8 percent half-window change.',
                 'Window, core and exponent sensitivities are not confidence intervals.',
                 'Slab persistence does not prove equilibrium vapor sampling.']))
    test = fitmod.critical_fit([230, 260, 267], .1 * (275-np.array([230, 260, 267]))**.325)
    assert abs(test['Tc_K']-275) < .001
    print('OUTPUT', OUT, flush=True)


if __name__ == '__main__':
    main()
