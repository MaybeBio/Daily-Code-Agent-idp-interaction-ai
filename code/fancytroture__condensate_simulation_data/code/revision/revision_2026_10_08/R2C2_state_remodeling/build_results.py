#!/usr/bin/env python3
"""Read-only C2 reanalysis of saved VAMPNet matrices; no training or MD."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = ROOT / 'revision_2026_09_30/state_audit'
SAMPLING = ROOT / 'revision_2026_09_30/state_sampling'
STYLE = ROOT / 'project_figure/scripts/make_figure2_dynamics_story.py'
SYSTEMS = ('FUS', 'A1', 'TDP43')
FAMILIES = ('free', 'condensate')
ROLES = ('Tl', 'Th')
METRICS = ('common_submatrix', 'common_plus_exit', 'common_conditional')
NAMES = {'FUS': 'FUS-LCD', 'A1': 'A1-LCD', 'TDP43': 'TDP-43-LCD'}
TITLES = {'common_submatrix': 'Common-state transitions', 'common_plus_exit': 'Common states + exit',
          'common_conditional': 'Conditional on staying in C'}
INPUTS = [AUDIT / name for name in ('full_local_transition_matrices.csv', 'distance_comparison.csv',
    'distance_contrasts.csv', 'population_and_kinetics_separated.csv', 'population_distances.csv')]
INPUTS += [SAMPLING / name for name in ('sampling_distances.csv', 'sampling_contrasts.csv',
                                       'absolute_time_exclusions.csv')]


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def write(name, rows):
    with (HERE / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def one(rows, **criteria):
    found = [r for r in rows if all(str(r[k]) == str(v) for k, v in criteria.items())]
    assert len(found) == 1, (criteria, len(found))
    return found[0]


def load_full(rows, system, family, role):
    subset = [r for r in rows if (r['system'], r['family'], r['role']) == (system, family, role)]
    assert len(subset) == 25
    matrix = np.full((5, 5), np.nan)
    mu = np.full(5, np.nan)
    mapping = {}
    for r in subset:
        i, j = int(r['source_local_state']), int(r['target_local_state'])
        matrix[i, j] = float(r['transition_probability'])
        if np.isfinite(mu[i]):
            np.testing.assert_allclose(mu[i], float(r['population']), atol=1e-14)
        mu[i] = float(r['population'])
        if r['source_display_state']:
            mapping[int(r['source_display_state'])] = i
    assert np.isfinite(matrix).all() and np.isfinite(mu).all()
    assert np.all(matrix >= 0) and np.all(mu >= 0)
    np.testing.assert_allclose(matrix.sum(1), 1., atol=1e-12)
    np.testing.assert_allclose(mu.sum(), 1., atol=1e-12)
    assert len(set(mapping.values())) == len(mapping)
    return dict(t=matrix, mu=mu, mapping=mapping)


def analyze():
    full = read(INPUTS[0])
    old = read(AUDIT / 'distance_comparison.csv')
    kinetics = read(AUDIT / 'population_and_kinetics_separated.csv')
    populations, coverage, transitions, distances, population_metrics = [], [], [], [], []
    pairs = {}
    for system in SYSTEMS:
        for family in FAMILIES:
            pair = [load_full(full, system, family, role) for role in ROLES]
            common = sorted(set(pair[0]['mapping']) & set(pair[1]['mapping']))
            union = sorted(set(pair[0]['mapping']) | set(pair[1]['mapping']))
            assert common
            pairs[system, family] = (pair, common, union)
            matrices, mus = [], []
            for role, data in zip(ROLES, pair):
                ids = [data['mapping'][s] for s in common]
                matrices.append(data['t'][np.ix_(ids, ids)])
                mus.append(data['mu'][ids])
                mapped_mass = float(sum(data['mu'][i] for i in data['mapping'].values()))
                common_mass = float(data['mu'][ids].sum())
                outside_mapped = mapped_mass - common_mass
                unmapped = 1. - mapped_mass
                np.testing.assert_allclose(common_mass + outside_mapped + unmapped, 1.)
                coverage.append(dict(system=system, family=family, role=role,
                    common_states_zero_based=' '.join(map(str, common)), common_population=common_mass,
                    mapped_noncommon_population=outside_mapped, unmapped_population=unmapped,
                    total_population=1.))
                for state in union:
                    present = state in data['mapping']
                    populations.append(dict(system=system, family=family, role=role,
                        display_state_zero_based=state, display_label=f's{state+1}', common=state in common,
                        mapped=present, population=float(data['mu'][data['mapping'][state]]) if present else '',
                        missing_policy='NA, not zero' if not present else 'original unrenormalized population'))
                populations.append(dict(system=system, family=family, role=role,
                    display_state_zero_based='', display_label='Unmapped', common=False, mapped=False,
                    population=unmapped, missing_policy='Aggregate only; not a shared physical state'))
            normalized = [mu / mu.sum() for mu in mus]
            weight = .5 * (normalized[0] + normalized[1])
            np.testing.assert_allclose(weight.sum(), 1.)
            exits = [1. - matrix.sum(1) for matrix in matrices]
            assert all(np.all((exit >= -1e-12) & (exit <= 1+1e-12)) for exit in exits)
            dcommon2 = np.sum(weight[:, None] * (matrices[1] - matrices[0])**2)
            dexit2 = np.sum(weight * (exits[1] - exits[0])**2)
            for metric in METRICS:
                if metric == 'common_submatrix':
                    a, b = matrices
                elif metric == 'common_plus_exit':
                    a, b = [np.column_stack((m, e)) for m, e in zip(matrices, exits)]
                    np.testing.assert_allclose(a.sum(1), 1.)
                    np.testing.assert_allclose(b.sum(1), 1.)
                else:
                    assert all(np.all(m.sum(1) > 0) for m in matrices)
                    a, b = [m / m.sum(1)[:, None] for m in matrices]
                distance = float(np.sqrt(np.sum(weight[:, None] * (b - a)**2)))
                reference = one(old, system=system, family=family, definition=metric)
                np.testing.assert_allclose(distance, float(reference['distance']), rtol=1e-10, atol=1e-12)
                if metric == 'common_plus_exit':
                    np.testing.assert_allclose(distance**2, dcommon2 + dexit2, atol=1e-14)
                distances.append(dict(system=system, family=family, metric=metric, distance=distance,
                    model_ci_low=float(reference['model_bootstrap_low']),
                    model_ci_high=float(reference['model_bootstrap_high']),
                    common_squared=float(dcommon2), exit_squared=float(dexit2),
                    exit_fraction_of_plus_exit_squared=float(dexit2/(dcommon2+dexit2))))
            midpoint = .5 * (normalized[0] + normalized[1])
            js = .5 * sum(np.sum(q[q > 0] * np.log2(q[q > 0] / midpoint[q > 0])) for q in normalized)
            population_metrics.append(dict(system=system, family=family,
                common_population_Tl=float(mus[0].sum()), common_population_Th=float(mus[1].sum()),
                common_conditional_JS_distance_base2=float(np.sqrt(max(js, 0.))),
                scope='Population distribution conditional on C, not the full ensemble or transition dynamics'))
            for k, state in enumerate(common):
                row = dict(system=system, family=family, display_state_zero_based=state,
                    display_label=f's{state+1}', lag_ns=2., pooled_common_source_weight=float(weight[k]))
                reference = one(kinetics, system=system, family=family, display_state=state)
                for role, data, mu, matrix, exit in zip(ROLES, pair, mus, matrices, exits):
                    row[f'population_{role}'] = float(mu[k])
                    row[f'self_transition_{role}'] = float(matrix[k,k])
                    row[f'exit_C_{role}'] = float(exit[k])
                    np.testing.assert_allclose(exit[k], float(reference[f'common_set_exit_probability_{role}']), atol=1e-12)
                row['exit_C_delta'] = row['exit_C_Th'] - row['exit_C_Tl']
                transitions.append(row)
    write('state_populations.csv', populations)
    write('state_coverage.csv', coverage)
    write('common_transition_distances_rechecked.csv', distances)
    write('population_only_distances.csv', population_metrics)
    write('per_state_exit_probabilities.csv', transitions)
    return pairs, coverage, transitions, distances


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(HERE / f'{name}.{ext}', dpi=200)
    plt.close(fig)


def population_plot(style, pairs):
    fig, axes = plt.subplots(3, 2, figsize=(13.2, 11.4))
    for row, system in enumerate(SYSTEMS):
        for col, family in enumerate(FAMILIES):
            ax = axes[row, col]
            pair, common, union = pairs[system, family]
            positions = np.arange(len(union)+1)
            for role, data, shift in zip(ROLES, pair, (-.19, .19)):
                for k, state in enumerate(union):
                    if state in data['mapping']:
                        value = 100 * data['mu'][data['mapping'][state]]
                        ax.bar(k+shift, value, width=.36, color=style.TEMP_COLORS[role])
                    else:
                        ax.text(k+shift, 2., 'NA', ha='center', va='bottom', fontsize=10, rotation=90)
                remaining = 100 * (1-sum(data['mu'][v] for v in data['mapping'].values()))
                ax.bar(len(union)+shift, remaining, width=.36, color=style.TEMP_COLORS[role], hatch='//',
                       edgecolor='white', linewidth=.6)
            labels = [rf'$s_{{{s+1}}}$' + ('*' if s not in common else '') for s in union] + ['Unmapped']
            ax.set_xticks(positions, labels, fontsize=13)
            ax.set_ylim(0, 75)
            ax.set_yticks([0,20,40,60])
            ax.set_ylabel('Population (%)')
            label = 'Free chain' if family == 'free' else 'Condensate'
            ax.set_title(f'{NAMES[system]} | {label}', fontsize=15, pad=10)
            masses = [sum(d['mu'][d['mapping'][s]] for s in common) for d in pair]
            ax.text(.98, .96, rf'$\mu_C(T_{{\mathrm{{L}}}})$ = {masses[0]:.1%}'+'\n'+
                    rf'$\mu_C(T_{{\mathrm{{H}}}})$ = {masses[1]:.1%}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=12)
            style.box_axis(ax)
    handles = [Line2D([0],[0], color=style.TEMP_COLORS[r], lw=5, label=style.temp_label(r)) for r in ROLES]
    fig.legend(handles=handles, loc='upper center', ncol=2, bbox_to_anchor=(.5,.963), frameon=False)
    fig.suptitle('State populations and common-state coverage', fontsize=19, y=.994)
    fig.text(.5,.018,'* Mapped at one temperature only; NA is not zero. Unmapped populations are not shared state identities.',
             ha='center', fontsize=11)
    fig.subplots_adjust(left=.08,right=.985,bottom=.075,top=.88,hspace=.43,wspace=.25)
    save(fig, 'C2_01_populations_and_coverage')


def distance_plot(style):
    sample = read(SAMPLING / 'sampling_distances.csv')
    contrast = read(SAMPLING / 'sampling_contrasts.csv')
    seeds = read(AUDIT / 'distance_contrasts.csv')
    summary = []
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 8.7))
    x = np.arange(3)
    for col, metric in enumerate(METRICS):
        top, bottom = axes[:, col]
        for family, shift, face in (('free', -.13, 'white'), ('condensate', .13, '#333333')):
            selected = [one(sample, system=s, family=family, metric=metric, comparison_block_frames=5000) for s in SYSTEMS]
            mean = np.array([float(r['estimate']) for r in selected])
            low = np.array([float(r['ci_low']) for r in selected])
            high = np.array([float(r['ci_high']) for r in selected])
            top.errorbar(x+shift, mean, yerr=np.array([mean-low,high-mean]), fmt='o', color='#333333',
                         markerfacecolor=face, markersize=7, capsize=4, elinewidth=1.3)
        for i, system in enumerate(SYSTEMS):
            a = one(contrast, system=system, metric=metric, block_frames=5000)
            b = one(seeds, system=system, definition=metric)
            for shift, record, keys, marker, color in [(-.10,a,('estimate','ci_low','ci_high'),'o','#333333'),
                    (.10,b,('condensate_minus_free','model_bootstrap_low','model_bootstrap_high'),'s','#427aab')]:
                mean, low, high = [float(record[k]) for k in keys]
                bottom.errorbar(i+shift,mean,yerr=np.array([[mean-low],[high-mean]]),fmt=marker,
                                color=color,markersize=6,capsize=4,elinewidth=1.3)
            summary.append(dict(system=system,metric=metric,
                sampling_contrast=float(a['estimate']),sampling_ci_low=float(a['ci_low']),sampling_ci_high=float(a['ci_high']),
                model_contrast=float(b['condensate_minus_free']),model_ci_low=float(b['model_bootstrap_low']),
                model_ci_high=float(b['model_bootstrap_high'])))
        top.set_title(TITLES[metric], fontsize=15, pad=12)
        top.set_ylim(0,.28)
        bottom.set_ylim(-.115,.155)
        bottom.axhline(0,color='#777777',ls=':',lw=1)
        for ax in (top,bottom):
            ax.set_xticks(x,[NAMES[s].replace('-LCD','') for s in SYSTEMS],fontsize=13)
            ax.set_xlim(-.5,2.5)
            style.box_axis(ax)
        if col == 0:
            top.set_ylabel('Transition distance')
            bottom.set_ylabel(r'$D_{\mathrm{condensate}}-D_{\mathrm{free}}$')
    top_handles = [Line2D([0],[0],marker='o',ls='',color='#333333',markerfacecolor=f,
                         label=l,markersize=7) for f,l in [('white','Free chain'),('#333333','Condensate')]]
    bottom_handles = [Line2D([0],[0],marker=m,ls='',color=c,label=l,markersize=6) for m,c,l in
                      [('o','#333333','Trajectory sampling'),('s','#427aab','Model seeds')]]
    fig.legend(handles=top_handles,loc='upper center',bbox_to_anchor=(.5,.96),ncol=2)
    fig.legend(handles=bottom_handles,loc='upper center',bbox_to_anchor=(.5,.505),ncol=2)
    fig.suptitle('Transition remodeling: definition and uncertainty',fontsize=19,y=.996)
    fig.text(.5,.018,'95% intervals shown separately; sampling uses fixed models and 50 ns condensate time blocks.',
             ha='center',fontsize=11)
    fig.subplots_adjust(left=.07,right=.99,top=.865,bottom=.085,wspace=.28,hspace=.55)
    save(fig,'C2_02_transition_distance_comparison')
    write('contrast_summary.csv',summary)
    sensitivity = [r for r in contrast if r['metric'] in METRICS]
    write('sampling_block_sensitivity.csv',sensitivity)
    return summary


def exit_plot(style, transitions):
    fig, axes = plt.subplots(3,2,figsize=(12.5,10.5))
    for row,system in enumerate(SYSTEMS):
        for col,family in enumerate(FAMILIES):
            ax=axes[row,col]
            data=[r for r in transitions if r['system']==system and r['family']==family]
            x=np.arange(len(data))
            for role,shift in zip(ROLES,(-.13,.13)):
                y=[r[f'exit_C_{role}'] for r in data]
                ax.scatter(x+shift,y,c=style.TEMP_COLORS[role],s=55,zorder=3)
            for i,r in enumerate(data):
                ax.plot([i-.13,i+.13],[r['exit_C_Tl'],r['exit_C_Th']],color='.7',lw=1,zorder=2)
            ax.set_xticks(x,[rf'$s_{{{int(r["display_state_zero_based"])+1}}}$' for r in data])
            ax.set_xlim(-.5,len(data)-.5)
            ax.set_ylim(-.005,.265)
            ax.set_yticks([0,.05,.10,.15,.20,.25])
            ax.set_ylabel(r'$P(i\rightarrow\mathrm{outside}\ C)$',fontsize=15)
            ax.set_title(f'{NAMES[system]} | '+('Free chain' if family=='free' else 'Condensate'),fontsize=15,pad=9)
            style.box_axis(ax)
    handles=[Line2D([0],[0],marker='o',ls='',color=style.TEMP_COLORS[r],label=style.temp_label(r)) for r in ROLES]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.963),ncol=2)
    fig.suptitle('Transitions outside the common-state set',fontsize=19,y=.995)
    fig.text(.5,.018,'Lag = 2 ns. Descriptive ensemble-mean probabilities, not continuous-time rates; no uncertainty intervals shown.',
             ha='center',fontsize=11)
    fig.subplots_adjust(left=.11,right=.985,top=.875,bottom=.075,wspace=.30,hspace=.43)
    save(fig,'C2_03_common_state_exit')


def main():
    spec=importlib.util.spec_from_file_location('maintext_style',STYLE)
    style=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(style)
    style.setup_style()
    pairs,coverage,transitions,distances=analyze()
    population_plot(style,pairs)
    summary=distance_plot(style)
    exit_plot(style,transitions)
    limits=[
        'Only saved VAMPNet results; no traditional MSM, retraining, MD, or changes to submission files.',
        'Existing cross-condition state mappings are held fixed, not newly validated structural equivalences.',
        'State labels are one-based display labels; Free and Condensate labels do not imply identical state definitions.',
        'Unavailable mapped states remain NA, never a measured zero population or transition row.',
        'Common submatrices are substochastic and not renormalized; conditional matrices and plus-exit distributions answer different questions.',
        'Outside C includes both mapped noncommon and unmapped original states; it has no source transition row.',
        'Model and sampling intervals are separate, not a combined confidence interval.',
        'Sampling holds models/mappings fixed; condensate chains are jointly resampled in absolute slab-time blocks, not independent chains.',
        'A1 260K condensate sampling excludes a short trajectory lacking absolute-time metadata (0.388% of frames); model intervals include it.',
        'Different Free/Condensate common sets and model partition uncertainty limit cross-environment distance interpretation.',
        'Probability-weighted common/exit distances depend on the stated source weights; no unique causal separation of population and rates is claimed.',
        'No physical disappearance, activation barrier, or universal condensate amplification is inferred.'
    ]
    metadata=dict(status='independently rechecked from exported full matrices',lag_ns=2.,
        checks=['All 12 full matrices row-normalized; populations sum to one',
                'All 18 common-state distances agree with existing audit to 1e-12 absolute / 1e-10 relative tolerance',
                'Missing mappings remain NA; common + noncommon + unmapped mass = 1',
                'Squared common-plus-exit distance equals common squared distance + weighted exit squared difference',
                'All per-state exit probabilities agree with prior kinetics table'],
        limitations=limits,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in INPUTS+[STYLE,Path(__file__)]})
    (HERE/'validation_and_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    print('Completed: 3 figures (PNG/PDF), 7 CSV tables, provenance. Original files untouched.')


if __name__=='__main__':
    main()
