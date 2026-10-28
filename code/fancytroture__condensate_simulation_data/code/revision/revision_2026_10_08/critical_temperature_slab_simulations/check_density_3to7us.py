"""Read-only matched-window slab-density check, excluding the first 3 us."""
import csv
from datetime import datetime
import json
import os
from pathlib import Path
from zoneinfo import ZoneInfo

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl-density-3to7us')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MDAnalysis.coordinates.DCD import DCDReader
import numpy as np
from openmm import unit

from common import HERE, load_case, sha, write_json

CACHE = HERE / 'tc_assessment_20261012'
OUT = HERE / 'density_3to7us_20261012'
CASES = ['FUS_335K', 'A1_267K', 'TDP43_315K']
EDGES = np.arange(-140., 141.)
Z = (EDGES[:-1] + EDGES[1:]) / 2
CONVERSION = 1.66053906660e-3
START, END = 3000, 7000


def table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def frame_density(ts, mass):
    np.testing.assert_allclose(ts.dimensions, [200, 200, 2800, 90, 90, 90], atol=.001)
    z = (ts.positions[:, 2].astype(float) * .1) % 280
    mode = np.sum(mass * np.exp(2j * np.pi * z / 280)) / mass.sum()
    center = (np.angle(mode) % (2 * np.pi)) * 280 / (2 * np.pi)
    centered = (z - center + 140) % 280 - 140
    profile = np.histogram(centered, EDGES, weights=mass)[0] / 400 * CONVERSION
    raw = np.histogram(z, EDGES + 140, weights=mass)[0] / 400 * CONVERSION
    return profile, raw, abs(mode)


def read_case(case):
    folder, _, config, system, _ = load_case(case)
    mass = np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                     for i in range(system.getNumParticles())])
    path = CACHE / (case + '.npz')
    with np.load(path) as data:
        keep = (data['time_ns'] > START) & (data['time_ns'] <= END)
        times = list(data['time_ns'][keep])
        profiles = list(data['centered_g_ml'][keep])
        lab = list(data['lab_g_ml'][keep])
        modes = list(data['fourier_mass_amplitude'][keep])
        np.testing.assert_array_equal(data['z_centered_nm'], Z)
    cached_last = int(times[-1])
    assert cached_last % 100 == 0
    old_inputs = {r['completed_step']: r for r in
                  json.loads((CACHE / (case + '_inputs.json')).read_text())}
    inputs = []
    for number in range(START // 100 + 1, END // 100 + 1):
        part = folder / 'runs' / f'part{number:04d}'
        metadata = part / 'complete.json'
        assert json.loads(metadata.read_text())['step'] == number * 10000000
        trajectory = part / 'phase.dcd'
        before = trajectory.stat()
        record = dict(path=str(trajectory), bytes=before.st_size,
                      completed_step=number * 10000000, completion_sha256=sha(metadata))
        if number <= cached_last // 100:
            assert record == old_inputs[record['completed_step']]
        else:
            with DCDReader(str(trajectory)) as reader:
                assert len(reader) == 1000 and reader.n_atoms == len(mass)
                assert np.isclose(reader.dt, 100., rtol=1e-6)
                for index in range(9, 1000, 10):
                    ts = reader[index]
                    time_ns = (number - 1) * 100 + (index + 1) * .1
                    assert np.isclose(ts.time / 1000, time_ns, atol=.001, rtol=0)
                    p, raw, mode = frame_density(ts, mass)
                    times.append(time_ns)
                    profiles.append(p)
                    lab.append(raw)
                    modes.append(mode)
            print(case, 'read through', number * 100, 'ns', flush=True)
        after = trajectory.stat()
        assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
        inputs.append(record)
    times, profiles, lab, modes = map(np.asarray, (times, profiles, lab, modes))
    np.testing.assert_array_equal(times, np.arange(START + 1, END + 1))
    assert profiles.shape == lab.shape == (END - START, len(Z))
    assert np.isfinite(profiles).all() and np.isfinite(lab).all()
    for density in (profiles, lab):
        np.testing.assert_allclose(density.sum(1) * 400 / CONVERSION, mass.sum(), rtol=1e-12)
    # Independently re-read early, middle and final frames of the selected window.
    for time_ns in [START + 1, 5001, END]:
        number = (time_ns - 1) // 100 + 1
        index = (time_ns - (number - 1) * 100) * 10 - 1
        with DCDReader(str(folder / 'runs' / f'part{number:04d}' / 'phase.dcd')) as reader:
            p, raw, mode = frame_density(reader[index], mass)
            i = time_ns - START - 1
            np.testing.assert_allclose(p, profiles[i], atol=1e-12, rtol=0)
            np.testing.assert_allclose(raw, lab[i], atol=1e-12, rtol=0)
            np.testing.assert_allclose(mode, modes[i], atol=1e-12, rtol=0)
    np.savez_compressed(OUT / (case + '.npz'), time_ns=times, z_centered_nm=Z,
                        centered_g_ml=profiles, lab_g_ml=lab, fourier_mass_amplitude=modes)
    write_json(OUT / (case + '_inputs.json'), dict(cache=str(path), cache_sha256=sha(path),
               cached_through_ns=cached_last, source_segments=inputs))
    return config, times, profiles, lab, modes


def describe(case, begin, end, times, rho, modes):
    selected = (times > begin) & (times <= end)
    p = rho[selected].mean(0)
    core = float(p[abs(Z) < 5].mean())
    dilute = float(p[abs(Z) > 80].mean())
    half = np.flatnonzero(p > (core + dilute) / 2)
    block_profiles = rho[selected].reshape(-1, 100, len(Z)).mean(1)
    block_core = block_profiles[:, abs(Z) < 5].mean(1)
    return dict(case=case, start_ns=begin, end_ns=end, n_frames=int(selected.sum()),
                core3_g_ml=float(p[abs(Z) < 3].mean()), core5_g_ml=core,
                core7_g_ml=float(p[abs(Z) < 7].mean()),
                dilute80_g_ml=dilute, dilute100_g_ml=float(p[abs(Z) >= 100].mean()),
                mean_box_g_ml=float(p.mean()),
                core100ns_block_sd=float(block_core.std(ddof=1)) if len(block_core) > 1 else None,
                half_height_width_nm=len(half),
                half_height_contiguous=bool(len(half) and np.all(np.diff(half) == 1)),
                mean_fourier_amplitude=float(modes[selected].mean()))


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(3, 3, figsize=(16, 10.8), layout='constrained')
    colors = ['#427cac', '#2f7f73', '#a46725', '#b54b5f']
    summary, one_us, hundred_ns = [], [], []
    for row, case in enumerate(CASES):
        config, times, rho, lab, modes = read_case(case)
        chunks = [describe(case, begin, begin + 1000, times, rho, modes)
                  for begin in range(START, END, 1000)]
        blocks = [describe(case, begin, begin + 100, times, rho, modes)
                  for begin in range(START, END, 100)]
        total = describe(case, START, END, times, rho, modes)
        first_half = rho[times <= 5000][:, abs(Z) < 5].mean()
        second_half = rho[times > 5000][:, abs(Z) < 5].mean()
        total.update(temperature_K=config['temperature_K'],
                     core_3to5us_g_ml=float(first_half), core_5to7us_g_ml=float(second_half),
                     last_half_vs_first_half_pct=float(100 * (second_half / first_half - 1)),
                     core_5to6us_g_ml=chunks[-2]['core5_g_ml'],
                     core_6to7us_g_ml=chunks[-1]['core5_g_ml'],
                     last1us_vs_previous1us_pct=100 * (chunks[-1]['core5_g_ml'] / chunks[-2]['core5_g_ml'] - 1))
        summary.append(total)
        one_us.extend(chunks)
        hundred_ns.extend(blocks)
        for chunk, color in zip(chunks, colors):
            select = (times > chunk['start_ns']) & (times <= chunk['end_ns'])
            axes[row, 0].plot(Z, rho[select].mean(0), color=color, lw=1.6,
                              label=f"{chunk['start_ns']/1000:g}-{chunk['end_ns']/1000:g}")
        axes[row, 0].set(xlim=(-140, 140), ylim=(0, .43), xlabel=r'Centered $z$ (nm)',
                         ylabel=r'$\rho$ (g/mL)', title=f"{case.split('_')[0]} | {config['temperature_K']} K")
        if row == 0:
            axes[row, 0].legend(frameon=False, fontsize=10, title=r'Time ($\mu$s)', ncol=2)
        grouped = lab.reshape(-1, 10, len(Z)).mean(1)
        im = axes[row, 1].imshow(grouped, origin='lower', aspect='auto',
                                 extent=(0, 280, 3, 7), vmin=0, vmax=.43,
                                 cmap='viridis', interpolation='nearest')
        axes[row, 1].set(xlabel=r'Laboratory $z$ (nm)', ylabel=r'Time ($\mu$s)')
        if row == 0:
            axes[row, 1].set_title('Uncentered density')
        mid = [(b['start_ns'] + b['end_ns']) / 2000 for b in blocks]
        axes[row, 2].plot(mid, [b['core5_g_ml'] for b in blocks], color='#2f7f73',
                          label=r'Core: $|z|<5$ nm')
        axes[row, 2].plot(mid, [b['dilute80_g_ml'] for b in blocks], color='#b54b5f',
                          label=r'Dilute: $|z|>80$ nm')
        axes[row, 2].set(xlim=(3, 7), ylim=(0, .43), xlabel=r'Time ($\mu$s)', ylabel=r'$\rho$ (g/mL)')
        if row == 0:
            axes[row, 2].set_title('100 ns block means')
            axes[row, 2].legend(frameon=False, fontsize=10)
        print(json.dumps(total), flush=True)
    fig.colorbar(im, ax=axes[:, 1].tolist(), label=r'$\rho$ (g/mL)', fraction=.025, pad=.02)
    for suffix in ['pdf', 'png']:
        fig.savefig(OUT / ('Slab_density_3to7us.' + suffix), dpi=180)
    plt.close(fig)
    table(OUT / 'summary.csv', summary)
    table(OUT / 'density_1us_blocks.csv', one_us)
    table(OUT / 'density_100ns_blocks.csv', hundred_ns)
    write_json(OUT / 'validation.json', dict(
        generated_at=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(), script_sha256=sha(__file__),
        selected_window_ns=[START, END], endpoint_convention='3000 < time_ns <= 7000',
        sample_dt_ns=1, bin_width_nm=1, frames_per_case=4000,
        mass_conservation_rtol=1e-12, independent_raw_frames_checked=9,
        frozen_inputs_and_completed_source_metadata_verified=True,
        time_and_box_checks_passed=True, simulations_modified=False,
        scope='Persistence and density stationarity in the requested window, not proof of equilibrium or a critical-temperature fit.',
        variability='100 ns block SD is descriptive variability, not an independent-replica confidence interval.'))
    print('OUTPUT', OUT, flush=True)


if __name__ == '__main__':
    main()
