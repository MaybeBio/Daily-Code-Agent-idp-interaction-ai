#!/usr/bin/env python3
"""Independent R2C2 audit with fixed TDP43 indices and no original-file writes."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/r2c2_mpl')
import numpy as np
from scipy.optimize import linear_sum_assignment

OUT = Path(__file__).resolve().parent
PARENT = OUT.parent
ROOT = PARENT.parents[1]
sys.path.insert(0, str(PARENT / 'tdp43_index_fix'))
import compare_fixed as fixed
p = fixed.previous
s = p.s
SAM = p.SAM


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def write(name, rows):
    assert rows
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def guard(stage):
    paths = []
    for directory in ('project_figure', 'project_figure_set_Tl_Th'):
        paths.extend(q for q in (ROOT / directory).rglob('*')
                     if q.is_file() and q.suffix.lower() in ('.png', '.pdf', '.svg', '.eps', '.py', '.csv'))
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            for role in s.ROLES:
                d = s.model_dir(system, family, role)
                paths += [d / 'avg_contacts.npy', d / 'avg_mu.npy']
                paths += list(d.glob('state_map_att*_raw_to_final.npy'))
    manifest = {str(q.relative_to(ROOT)): sha(q) for q in sorted(set(paths))}
    target = OUT / 'original_files_before.json'
    if stage == 'before' and not target.exists():
        target.write_text(json.dumps(manifest, indent=2) + '\n')
    else:
        assert json.loads(target.read_text()) == manifest, 'Original file changed; do not overwrite user edits'
    return len(manifest)


def tag(system, family, role):
    return f'{system}_{family}_{s.TEMPS[system][role]}'


def contacts(system, family, role):
    if system == 'TDP43':
        with np.load(PARENT / 'tdp43_index_fix' / tag(system, family, role) / 'fixed_alignment.npz') as z:
            assert np.array_equal(z['new_to_old'], np.tile(np.arange(5), (20, 1)))
            return z['avg_contacts'].astype(float), z['avg_mu'].astype(float)
    d = s.model_dir(system, family, role)
    return np.load(d / 'avg_contacts.npy').astype(float), np.load(d / 'avg_mu.npy').astype(float)


def mapping(system, family):
    return fixed.match(family, 'fully_fixed')[0] if system == 'TDP43' else p.maps(system, family)


def match(corr, mu_low, mu_high, pop_cut, corr_cut):
    low, high = np.flatnonzero(mu_low >= pop_cut), np.flatnonzero(mu_high >= pop_cut)
    if not len(low) or not len(high):
        return []
    rr, cc = linear_sum_assignment(1 - corr[np.ix_(high, low)])
    return [(int(low[j]), int(high[i])) for i, j in zip(rr, cc)
            if corr[high[i], low[j]] >= corr_cut]


def matching_inventory():
    agreement = read(SAM / 'state_agreement.csv')
    inventory, matrices, scans, time_rows = [], [], [], []
    for system in s.SYSTEMS:
        for family in s.FAMILIES:
            maps = mapping(system, family)
            common = sorted(set(maps[0]) & set(maps[1]))
            pairs = [(maps[0][i], maps[1][i]) for i in common]
            cm, mu = zip(*(contacts(system, family, role) for role in s.ROLES))
            assert cm[0].shape == cm[1].shape
            corr = p.contact_features(cm[1]) @ p.contact_features(cm[0]).T
            pop_cut = .02 if system == 'FUS' else .10
            if system != 'FUS':
                assert set(match(corr, mu[0], mu[1], pop_cut, .9)) == set(pairs)
            for lo in range(5):
                for hi in range(5):
                    matrices.append(dict(system=system, family=family, local_TL=lo+1, local_TH=hi+1,
                                         correlation=float(corr[hi, lo]), accepted=(lo, hi) in pairs))
            for cut in (0., .02, .05, .10):
                for threshold in (.80, .85, .90, .95):
                    accepted = match(corr, mu[0], mu[1], cut, threshold)
                    scans.append(dict(system=system, family=family, pop_cut=cut, correlation_cut=threshold,
                                      pairs=json.dumps([(l+1, h+1) for l, h in accepted]), n_common=len(accepted),
                                      note='diagnostic; FUS primary mapping uses its original matcher'))
            for ri, role in enumerate(s.ROLES):
                with np.load(SAM / f'{tag(system, family, role)}_sufficient_statistics.npz') as z:
                    model_mu = p.pm(z['masses'].sum((1, 2)))
                np.testing.assert_allclose(model_mu.mean(0), mu[ri], atol=1e-6)
                opposite = corr.T if ri == 0 else corr
                kept_other = sorted(maps[1-ri].values())
                inv = {v: k for k, v in maps[ri].items()}
                paper_inv = {v: k for k, v in p.maps(system, family)[ri].items()}
                for local in range(5):
                    state_id = inv.get(local)
                    status = ('common' if state_id in common else
                              'retained_unmatched' if state_id is not None else 'below_population_cutoff')
                    best_all = int(np.argmax(opposite[local]))
                    best_kept = max(kept_other, key=lambda j: opposite[local, j])
                    a = next(r for r in agreement if (r['system'], r['family'], r['role'], int(r['local_state']))
                             == (system, family, role, local))
                    inventory.append(dict(system=system, family=family, role=role, local_state=local+1,
                        paper_state='' if local not in paper_inv else int(paper_inv[local])+1,
                        response_state='' if state_id is None else int(state_id)+1, status=status,
                        population=float(mu[ri][local]), population_model_min=float(model_mu[:, local].min()),
                        population_model_max=float(model_mu[:, local].max()),
                        models_above_population_cut=int((model_mu[:, local] >= pop_cut).sum()), models=len(model_mu),
                        population_cut=pop_cut, soft_dice_mean=float(a['soft_dice_mean']),
                        soft_dice_p05=float(a['soft_dice_p05']), best_other_local=best_all+1,
                        best_other_population=float(mu[1-ri][best_all]), best_other_correlation=float(opposite[local, best_all]),
                        best_retained_other_local=best_kept+1,
                        best_retained_correlation=float(opposite[local, best_kept])))
                if family == 'condensate':
                    with np.load(SAM / f'{tag(system, family, role)}_absolute_statistics.npz') as z:
                        mass = z['masses'].reshape(20, -1, 5, 5).sum(2)
                        origin = int(z['absolute_origin'])
                    coverage = mass.sum(-1).mean(0)
                    use = coverage > 0
                    block_mu = mass[:, use] / mass[:, use].sum(-1, keepdims=True)
                    for j, bi in enumerate(np.flatnonzero(use)):
                        for k in range(5):
                            values = block_mu[:, j, k]
                            time_rows.append(dict(system=system, role=role, local_state=k+1,
                                source_time_ns=(origin+bi*5000+2500)*.01, selected_chain_frames=float(coverage[bi]),
                                population=float(values.mean()), model_p025=float(np.quantile(values,.025)),
                                model_p975=float(np.quantile(values,.975))))
    write('state_inventory.csv', inventory)
    write('all_state_contact_correlations.csv', matrices)
    write('matching_threshold_sensitivity.csv', scans)
    write('absolute_50ns_population_blocks.csv', time_rows)
    return inventory


def separated_metrics():
    rows, contrasts, pops, coverage = [], [], [], []
    values = {}
    with np.load(PARENT / 'followup/bootstrap_draws.npz') as z:
        for system in s.SYSTEMS:
            for family in s.FAMILIES:
                maps = mapping(system, family)
                common = sorted(set(maps[0]) & set(maps[1]))
                for source in ('sampling', 'model_seeds'):
                    arrays = [(z[f'{system}_{family}_{role}_{source}_t'], z[f'{system}_{family}_{role}_{source}_mu'])
                              for role in s.ROLES]
                    metrics, _, ids, ts, ps = p.metrics(arrays, maps)
                    for metric, val in metrics.items():
                        values[system, family, source, metric] = val
                        rows.append(dict(system=system, family=family, uncertainty=source, metric=metric, **p.summary(val)))
                    for ri, role in enumerate(s.ROLES):
                        coverage.append(dict(system=system, family=family, role=role, uncertainty=source,
                                             **p.summary(ps[ri].sum(-1))))
                        for k in range(5):
                            pops.append(dict(system=system, family=family, role=role, local_state=k+1,
                                             uncertainty=source, **p.summary(arrays[ri][1][:, k])))
                    if source == 'model_seeds':
                        reference = (next(r for r in read(PARENT/'tdp43_index_fix/distance_results.csv')
                                     if r['variant']=='fully_fixed' and float(r['cutoff'])==.9 and r['family']==family
                                     and r['uncertainty']==source and r['metric']=='D_common_pooled_weights')
                                     if system=='TDP43' else
                                     next(r for r in read(PARENT/'followup/separated_metrics.csv') if r['system']==system
                                          and r['family']==family and r['uncertainty']==source and r['metric']=='D_common_pooled_weights'))
                        np.testing.assert_allclose(metrics['D_common_pooled_weights'][0], float(reference['estimate']), atol=1e-10)
            for source in ('sampling','model_seeds'):
                for metric in metrics:
                    contrasts.append(dict(system=system, uncertainty=source, metric=metric,
                        **p.summary(values[system,'condensate',source,metric]-values[system,'free',source,metric])))
    write('separated_metrics.csv', rows)
    write('condensate_minus_free.csv', contrasts)
    write('population_intervals.csv', pops)
    write('common_population_coverage.csv', coverage)
    return rows, contrasts


def main():
    n = guard('before')
    inventory = matching_inventory()
    rows, contrasts = separated_metrics()
    (OUT/'settings.json').write_text(json.dumps(dict(
        scope='Independent R2C2 response analysis; originals unchanged; no training or GPU use',
        lag_ns=2., bootstrap_draws=4000, bootstrap_source='followup/bootstrap_draws.npz',
        mapping='Original FUS/A1 mapping; fully repaired TDP43 mapping with unchanged seed permutations',
        weights='Half sum of separately common-normalized populations; no submatrix row renormalization',
        population_metric='Conditional-common-state TV/JS, not total ensemble remodeling',
        sampling='Free whole trajectories; condensate 50 ns absolute time blocks; excludes A1 TH chain139 with missing metadata',
        limitations=['Fixed mappings during bootstrap; seed and sampling intervals are separate, not total uncertainty',
                     '50 ns population traces are conditional on selected chain windows, not whole-slab stationarity tests',
                     'Two temperatures trained independently; no shared state function fitted'],
        structure_stride_frames=20, structure_stride_control=40, structure_contact_cutoff_nm=.66,
        structure_filter='None; all sampled frames weighted by existing soft assignments',
        original_files_guarded=n), indent=2)+'\n')
    print('Unmatched retained states:')
    for row in inventory:
        if row['status']=='retained_unmatched':
            print(json.dumps(row), flush=True)
    print('Weighted-distance contrasts:')
    for row in contrasts:
        if row['metric']=='D_common_pooled_weights': print(json.dumps(row), flush=True)
    print('Original files unchanged:', guard('after'))


if __name__=='__main__':
    main()
