"""Compact SI density profiles from the validated 3-7 us slab window."""
import csv
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl-density-si')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / 'density_3to7us_20261012'
OUT = DATA / 'si_figure'
STYLE = HERE.parents[1] / 'project_figure/scripts/make_figure2_dynamics_story.py'
CASES = [('FUS_335K', 'FUS-LCD', 335), ('A1_267K', 'A1-LCD', 267),
         ('TDP43_315K', 'TDP-43-LCD', 315)]


def main():
    OUT.mkdir(exist_ok=True)
    spec = importlib.util.spec_from_file_location('maintext_style', STYLE)
    style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(style)
    style.setup_style()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.1))
    source = []
    for col, (case, title, temp) in enumerate(CASES):
        with np.load(DATA / (case + '.npz')) as data:
            np.testing.assert_array_equal(data['time_ns'], np.arange(3001, 7001))
            z, rho = data['z_centered_nm'], data['centered_g_ml']
        assert rho.shape == (4000, len(z))
        blocks = rho.reshape(4, 1000, len(z)).mean(1)
        mean, sd = rho.mean(0), blocks.std(0, ddof=1)
        np.testing.assert_allclose(mean, blocks.mean(0), atol=1e-12)
        ax = axes[col]
        ax.fill_between(z, np.maximum(0, mean - sd), mean + sd,
                        color='#454545', alpha=.16, linewidth=0)
        ax.plot(z, mean, color='#333333', lw=2.3)
        ax.set(xlim=(-80, 80), ylim=(0, .45), xlabel=r'$z$ (nm)')
        ax.set_xticks([-80, -40, 0, 40, 80])
        ax.set_yticks([0, .1, .2, .3, .4])
        ax.set_title(title, pad=12)
        if col == 0:
            ax.set_ylabel(r'$\rho(z)$ (g/mL)')
        ax.text(.95, .93, f'{temp} K', transform=ax.transAxes,
                ha='right', va='top', fontsize=14.8)
        ax.text(-.13, 1.055, chr(65 + col), transform=ax.transAxes,
                fontsize=20, fontweight='bold')
        style.box_axis(ax)
        ax.grid(False)
        for zz, mm, ss in zip(z, mean, sd):
            source.append(dict(system=title, temperature_K=temp, z_nm=float(zz),
                               mean_density_g_ml=float(mm), block_sd_g_ml=float(ss),
                               start_ns=3000, end_ns=7000, n_blocks=4, block_ns=1000))
    fig.subplots_adjust(left=.075, right=.985, bottom=.21, top=.85, wspace=.23)
    for ext in ['png', 'pdf']:
        fig.savefig(OUT / ('Fig_SI_slab_density_3to7us.' + ext), dpi=300)
    plt.close(fig)
    with (OUT / 'density_profiles_source.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(source[0]))
        writer.writeheader()
        writer.writerows(source)
    caption = (
        'Slab density profiles at the additional simulation temperatures. '
        '(A) FUS-LCD at 335 K, (B) A1-LCD at 267 K, and (C) TDP-43-LCD at 315 K. '
        'Solid curves show mean mass-density profiles over 3-7 microseconds, excluding '
        'the first 3 microseconds. Each frame was centered using the mass-weighted '
        'circular mean along the periodic slab-normal direction; profiles were sampled '
        'every 1 ns with 1 nm bins. Shading denotes one standard deviation across '
        'four consecutive 1-microsecond block means, truncated at zero for display; '
        'it describes temporal variability, not a confidence interval. The central '
        'region of the box is shown. All three systems retain a dense slab separated '
        'from low-density surroundings during this window.\n'
    )
    (OUT / 'CAPTION_EN.txt').write_text(caption)
    (OUT / 'figure_manifest.json').write_text(json.dumps(dict(
        input_directory=str(DATA), window_ns=[3000, 7000],
        density_definition='Mass per 1 nm axial bin; full 20 x 20 nm cross-section.',
        style_source=str(STYLE), maintext_style_reused=True,
        smoothing=False, shaded_quantity='SD across four 1 us block means, not CI',
        shown_z_range_nm=[-80, 80], source_csv_includes_full_box=True,
        simulations_modified=False, manuscript_and_SI_not_rebuilt=True), indent=2) + '\n')
    print('OUTPUT', OUT)


if __name__ == '__main__':
    main()
