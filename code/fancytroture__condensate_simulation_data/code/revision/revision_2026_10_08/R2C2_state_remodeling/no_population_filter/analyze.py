#!/usr/bin/env python3
"""Response-only common-state distance: all five states enter matching."""
import csv
import importlib.util
import itertools
import json
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/r2c2_mpl')
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('state_audit', OUT.parent/'state_validation/analyze.py')
a = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(a)
p, s = a.p, a.s
CUTOFF = .90


def write(name, rows):
    with (OUT/name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def assign(corr):
    # Preserve the preceding diagnostic's order: assignment, then acceptance.
    assert corr.shape == (5, 5) and np.isfinite(corr).all()
    high, low = linear_sum_assignment(1 - corr)
    optimum = max(sum(corr[h, perm[h]] for h in range(5))
                  for perm in itertools.permutations(range(5)))
    np.testing.assert_allclose(corr[high, low].sum(), optimum, atol=1e-12)
    all_pairs = [(int(l), int(h)) for h, l in zip(high, low)]
    accepted = [(l, h) for l, h in all_pairs if corr[h, l] >= CUTOFF]
    assert accepted and len({l for l, _ in accepted}) == len(accepted)
    assert len({h for _, h in accepted}) == len(accepted)
    return all_pairs, accepted


def independent_distance(tl, th, ml, mh, pairs):
    il, ih = np.array(pairs).T
    low = tl[:, il][:, :, il]
    high = th[:, ih][:, :, ih]
    pl, ph = ml[:, il], mh[:, ih]
    w = .5*(pl/pl.sum(-1, keepdims=True) + ph/ph.sum(-1, keepdims=True))
    return np.sqrt(np.einsum('bi,bij,bij->b', w, high-low, high-low))


def tests():
    _, pairs = assign(np.eye(5))
    assert pairs == [(i, i) for i in range(5)]
    # A row's retained probability is not renormalized: escape changes matter.
    t0 = np.array([[[.8, .2], [.3, .7]]])
    t1 = np.array([[[.6, .4], [.3, .7]]])
    mu = np.array([[.9, .1]])
    value = independent_distance(t0, t1, mu, mu, [(0, 0)])
    np.testing.assert_allclose(value, .2)
    same = independent_distance(t0, t0, mu, mu, [(0, 0), (1, 1)])
    np.testing.assert_array_equal(same, [0.])
    calc = p.metrics([(t0, mu), (t1, mu)], [{0: 0}, {0: 0}])[0]
    np.testing.assert_allclose(calc['D_common_pooled_weights'], value)


def plot(distances, contrasts):
    plt.rcParams.update({'font.size': 13, 'axes.titlesize': 15, 'axes.labelsize': 14,
                         'xtick.labelsize': 13, 'ytick.labelsize': 12,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8), sharex=True, sharey='col')
    names = ['FUS-LCD', 'A1-LCD', 'TDP-43-LCD']
    for row, source in enumerate(('sampling', 'model_seeds')):
        ax, diff = axes[row]
        for i, system in enumerate(s.SYSTEMS):
            for family, offset, marker, fill in [('free', -.10, 'o', 'white'),
                                                 ('condensate', .10, 's', '0.2')]:
                r = next(r for r in distances if (r['system'], r['family'], r['uncertainty']) ==
                         (system, family, source))
                x = i + offset
                ax.vlines(x, r['ci_low'], r['ci_high'], color='0.2', lw=1.5)
                ax.hlines([r['ci_low'], r['ci_high']], x-.035, x+.035, color='0.2', lw=1.5)
                ax.plot(x, r['estimate'], marker=marker, mfc=fill, mec='0.2', ms=8)
            r = next(r for r in contrasts if (r['system'], r['uncertainty']) == (system, source))
            diff.vlines(i, r['ci_low'], r['ci_high'], color='0.2', lw=1.5)
            diff.hlines([r['ci_low'], r['ci_high']], i-.04, i+.04, color='0.2', lw=1.5)
            diff.plot(i, r['estimate'], 'o', color='0.2', ms=7)
        label = 'Trajectory sampling' if source == 'sampling' else 'Model seeds'
        ax.set_title(label)
        diff.set_title(label)
        ax.set_ylabel(r'$D_T^{(w)}$')
        diff.set_ylabel(r'$D_{\mathrm{condensate}}-D_{\mathrm{free}}$')
        diff.axhline(0, color='.5', ls=':', lw=1.2)
        for target in (ax, diff):
            target.set_xticks(range(3), names)
            target.set_xlim(-.45, 2.45)
            target.grid(axis='y', alpha=.18)
            target.set_axisbelow(True)
            target.tick_params(labelbottom=True)
    axes[0, 0].set_ylim(bottom=0)
    fig.suptitle('Common-state transitions: no population filtering', fontsize=18, y=.98)
    fig.legend(handles=[Line2D([], [], marker='o', mfc='white', mec='0.2', ls='', label='Free chain', ms=8),
                        Line2D([], [], marker='s', color='0.2', ls='', label='Condensate', ms=8)],
               loc='upper center', bbox_to_anchor=(.5, .94), ncol=2, frameon=False)
    fig.subplots_adjust(top=.82, bottom=.08, left=.09, right=.98, hspace=.40, wspace=.36)
    for ext in ('png', 'pdf'):
        fig.savefig(OUT/f'DTw_no_population_filter.{ext}', dpi=220)
    plt.close(fig)


def main():
    tests()
    protected = a.guard('after')
    sources = {str(Path(__file__).resolve()): a.sha(Path(__file__))}
    bootstrap = OUT.parent/'followup/bootstrap_draws.npz'
    sources[str(bootstrap)] = a.sha(bootstrap)
    distances, contrasts, assignment, correlations, coverage, transitions, weights = [], [], [], [], [], [], []
    draws = {}
    diagnostic = a.read(OUT.parent/'state_validation/matching_threshold_sensitivity.csv')
    with np.load(bootstrap) as z:
        for system in s.SYSTEMS:
            for family in s.FAMILIES:
                cm, mu = zip(*(a.contacts(system, family, role) for role in s.ROLES))
                features = [p.contact_features(m) for m in cm]
                corr = features[1] @ features[0].T
                np.testing.assert_allclose(corr, np.corrcoef(np.concatenate(features))[5:, :5], atol=1e-12)
                all_pairs, pairs = assign(corr)
                expected = next(r for r in diagnostic if r['system']==system and r['family']==family
                                and float(r['pop_cut'])==0 and float(r['correlation_cut'])==CUTOFF)
                assert [[l+1, h+1] for l, h in pairs] == json.loads(expected['pairs'])
                maps = [{i: pair[ri] for i, pair in enumerate(pairs)} for ri in (0, 1)]
                for l in range(5):
                    for h in range(5):
                        correlations.append(dict(system=system, family=family, local_TL=l+1,
                            local_TH=h+1, correlation=float(corr[h,l]), assigned=(l,h) in all_pairs,
                            accepted=(l,h) in pairs))
                for l, h in all_pairs:
                    assignment.append(dict(system=system, family=family, local_TL=l+1, local_TH=h+1,
                        population_TL=float(mu[0][l]), population_TH=float(mu[1][h]),
                        correlation=float(corr[h,l]), accepted=(l,h) in pairs))
                for ri, role in enumerate(s.ROLES):
                    d = s.model_dir(system, family, role)
                    source_paths = ([OUT.parent/'tdp43_index_fix'/a.tag(system,family,role)/'fixed_alignment.npz']
                                    if system=='TDP43' else [d/'avg_contacts.npy', d/'avg_mu.npy'])
                    source_paths += [a.SAM/f'{a.tag(system,family,role)}_sufficient_statistics.npz']
                    for source_path in source_paths:
                        sources[str(source_path)] = a.sha(source_path)
                    with np.load(source_paths[-1]) as raw:
                        assert int(raw['lag'])==200 and float(raw['dt'])==.01
                        tcheck = p.transition(raw['counts'].sum((1,2))).mean(0)
                        mcheck = p.pm(raw['masses'].sum((1,2))).mean(0)
                    sampling_t, sampling_mu = tcheck, mcheck
                    if family == 'condensate':
                        absolute_path = a.SAM/f'{a.tag(system,family,role)}_absolute_statistics.npz'
                        sources[str(absolute_path)] = a.sha(absolute_path)
                        with np.load(absolute_path) as raw:
                            sampling_t = p.transition(raw['counts'].sum(1)).mean(0)
                            sampling_mu = p.pm(raw['masses'].sum(1)).mean(0)
                    for source in ('sampling', 'model_seeds'):
                        key = f'{system}_{family}_{role}_{source}'
                        expected_t = sampling_t if source=='sampling' else tcheck
                        expected_mu = sampling_mu if source=='sampling' else mcheck
                        np.testing.assert_allclose(z[key+'_t'][0], expected_t, atol=1e-10)
                        np.testing.assert_allclose(z[key+'_mu'][0], expected_mu, atol=1e-10)
                    np.testing.assert_allclose(mcheck, mu[ri], atol=1e-6)
                for source in ('sampling', 'model_seeds'):
                    arrays = [(z[f'{system}_{family}_{role}_{source}_t'],
                               z[f'{system}_{family}_{role}_{source}_mu']) for role in s.ROLES]
                    for t, m in arrays:
                        assert np.isfinite(t).all() and np.isfinite(m).all() and np.min(t)>=0
                        np.testing.assert_allclose(t.sum(-1), 1., atol=1e-10)
                        np.testing.assert_allclose(m.sum(-1), 1., atol=1e-10)
                    metrics, _, ids, ts, ps = p.metrics(arrays, maps)
                    values = metrics['D_common_pooled_weights']
                    direct = independent_distance(arrays[0][0], arrays[1][0], arrays[0][1], arrays[1][1], pairs)
                    np.testing.assert_allclose(values, direct, atol=1e-14)
                    reverse = independent_distance(arrays[1][0], arrays[0][0], arrays[1][1], arrays[0][1],
                                                   [(h,l) for l,h in pairs])
                    np.testing.assert_allclose(values, reverse, atol=1e-14)
                    draws[f'{system}_{family}_{source}'] = values
                    distances.append(dict(system=system, family=family, uncertainty=source,
                                          n_common=len(pairs), **p.summary(values)))
                    for ri, role in enumerate(s.ROLES):
                        coverage.append(dict(system=system, family=family, role=role, uncertainty=source,
                                             **p.summary(ps[ri].sum(-1))))
                    if source=='model_seeds':
                        w = .5*(ps[0][0]/ps[0][0].sum() + ps[1][0]/ps[1][0].sum())
                        for i, (l,h) in enumerate(pairs):
                            weights.append(dict(system=system, family=family, pair=i+1, local_TL=l+1,
                                                local_TH=h+1, weight=float(w[i])))
                            for j in range(len(pairs)):
                                transitions.append(dict(system=system, family=family, from_pair=i+1,
                                    to_pair=j+1, T_L=float(ts[0][0,i,j]), T_H=float(ts[1][0,i,j]),
                                    weighted_squared_difference=float(w[i]*(ts[1][0,i,j]-ts[0][0,i,j])**2)))
                print(system, family, 'pairs', [(l+1,h+1) for l,h in pairs],
                      'DTw', float(values[0]), flush=True)
            for source in ('sampling', 'model_seeds'):
                value = draws[f'{system}_condensate_{source}'] - draws[f'{system}_free_{source}']
                contrasts.append(dict(system=system, uncertainty=source, **p.summary(value)))
    for name, rows in [('distances.csv',distances), ('condensate_minus_free.csv',contrasts),
                       ('all_assignments.csv',assignment), ('all_correlations.csv',correlations),
                       ('common_population_coverage.csv',coverage), ('common_transition_entries.csv',transitions),
                       ('common_state_weights.csv',weights)]:
        write(name, rows)
    np.savez_compressed(OUT/'distance_bootstrap_draws.npz', **draws)
    plot(distances, contrasts)
    assert protected == a.guard('after')
    metadata = dict(population_filter=None, correlation_threshold=CUTOFF,
        matching='All 5x5 pairs; Hungarian maximum total log-contact correlation; reject assigned pairs below 0.90 afterward',
        contact_features='log10(max(p_contact,1e-4)), upper triangle |i-j|>4; mean-centered Pearson correlation',
        scope='Independent response-only analysis; no original manuscript/SI/figure/mapping changes; no GPU use',
        TDP43='Corrected 148-residue indexing; within-temperature seed permutations unchanged',
        FUS='Unified contact-correlation matcher, not an assertion of reproducing its historical matcher',
        distance='sqrt(sum_i w_i sum_j (T_H[hi,hj]-T_L[li,lj])^2), i,j in accepted common pairs',
        weights='w_i=0.5*(mu_L[li]/sum_C(mu_L)+mu_H[hi]/sum_C(mu_H))',
        lag_ns=2., bootstrap_draws=4000, zero_imputation=False, row_renormalization=False, edge_mask=False,
        uncertainty='Separate sampling and model-seed bootstrap intervals; mapping held fixed',
        sampling='Whole free trajectories; absolute 50 ns slab-time blocks; excludes A1 TH chain139 lacking metadata',
        point_estimates='Each uncertainty row uses its own data scope: model seeds includes all trajectories; sampling excludes unprovenanced A1 TH chain139',
        limitations=['Assignment then rejection is not threshold-constrained maximum-cardinality matching',
                     'Matching is not a shared state function; a contact correlation does not prove identical boundaries',
                     'Rare-state reliability and uncertainty of the matching itself are not established by this calculation',
                     'No population cutoff does not mean all five states are accepted as common'],
        guarded_original_files=protected, original_files_unchanged=True, source_sha256=sources)
    (OUT/'validation_and_settings.json').write_text(json.dumps(metadata, indent=2)+'\n')
    report = ['R2C2: no population filtering; response-only calculation', '',
              'All five states enter Hungarian contact-map matching. Assigned pairs with r >= 0.90 are retained.',
              'No zero rows, no network edge mask, no submatrix row renormalization; lag 2 ns.',
              'TDP43 uses repaired indices. FUS uses the unified diagnostic matcher.', '',
              'system  n_free  n_cond  D_free  D_cond  cond-free']
    for system in s.SYSTEMS:
        f,c = [next(r for r in distances if r['system']==system and r['family']==fam
                    and r['uncertainty']=='model_seeds') for fam in s.FAMILIES]
        report.append(f"{system}  {f['n_common']}  {c['n_common']}  {f['estimate']:.6f}  {c['estimate']:.6f}  {c['estimate']-f['estimate']:+.6f}")
        for source in ('sampling','model_seeds'):
            r = next(r for r in contrasts if r['system']==system and r['uncertainty']==source)
            report.append(f"  {source}: difference point {r['estimate']:+.6f}; 95% interval [{r['ci_low']:+.6f}, {r['ci_high']:+.6f}]")
    report += ['', 'Intervals are conditional on the fixed matching; not total uncertainty.',
               'Main table uses all-trajectory point estimates. Sampling row excludes A1 TH chain139 with missing time metadata, so its A1 point differs slightly.',
               'Do not infer state emergence/absence or uniform condensate enhancement from matching alone.',
               'Assignment changes can alter common-state coverage substantially; consult coverage and assignments tables.',
               f'Original-file hash checks passed for {protected} files.']
    report += ['', 'Interpretation of this particular matching:',
               'All three point contrasts are positive. FUS: both intervals cross zero; A1: both are positive;',
               'TDP43: sampling interval is positive but model-seed interval crosses zero.',
               'This is conditional evidence, not robustness to the state-matching definition.',
               'TDP43 condensate: TL local4 (0.482%) matches TH local3 (25.056%), correlation 0.9152.',
               'TL local1 instead matches TH local2; TL local2 is rejected with TH local4.',
               'Thus the TDP43 sign change is accompanied by changed major-state correspondence, not just adding rare pairs.',
               'A1 free loses the previously accepted TL local2 / TH local3 pair under all-state assignment then rejection.',
               'The algorithm maximizes the total correlation of all five assignments before rejecting weak pairs;',
               'it does not maximize the number or correlation sum of accepted pairs alone.']
    (OUT/'RESULTS.txt').write_text('\n'.join(report)+'\n')
    print('\n'.join(report), flush=True)


if __name__=='__main__':
    main()
