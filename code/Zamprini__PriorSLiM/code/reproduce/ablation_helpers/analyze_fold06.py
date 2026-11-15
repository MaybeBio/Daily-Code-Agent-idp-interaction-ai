from pathlib import Path
import sys
import logging
import gc
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import linear_sum_assignment
from scipy.special import expit
from sklearn.metrics import average_precision_score
import torch
HERE = Path(__file__).resolve().parent
ROOT = Path(__import__('os').environ['PRIORSLIM_DATA_ROOT'])
import group_deletion as dep
RESULTS = HERE / 'results'
OUT = HERE / 'analysis'
SEEDS = tuple(range(2026, 2031))
LABELS = {'A': 'MotifAE-only', 'B': 'Intrinsic + MotifAE', 'C': 'PriorSLiM', 'D': 'Explicit-only', 'B1': 'Intrinsic without conservation', 'B2': 'Intrinsic without motif', 'B3': 'Intrinsic without chemistry', 'CE': 'PriorSLiM + ensemble', 'S': 'Sequential frozen explicit', 'J_SCRATCH': 'Joint from scratch', 'H3': 'PriorSLiM + S/Y/F density', 'AROMATIC_ADD': 'PriorSLiM + F/Y/W descriptors', 'AROMATIC_SPLIT': 'PriorSLiM with split aromatic descriptors'}
DEPENDENCE_ARMS = ['B', 'C', 'B1', 'B2', 'B3', 'CE']

def save(name, rows):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(OUT / name, sep='\t', index=False)
    return frame

def summarize(frame, group, metrics):
    rows = []
    for keys, part in frame.groupby(group, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(group, keys))
        row['n_seeds'] = int(part.seed.nunique())
        for metric in metrics:
            row[metric + '_mean'] = float(part[metric].mean())
            row[metric + '_sd'] = float(part[metric].std(ddof=1))
        rows.append(row)
    return pd.DataFrame(rows)

def verify_all():
    required = [RESULTS / f'{arm}_seed{seed}/complete.json' for arm in LABELS for seed in SEEDS]
    required += [HERE / f'E4/C_seed{seed}/{phase}_complete.json' for seed in SEEDS for phase in ['decoder', 'profile', 'controls', 'refit']]
    required += [RESULTS / f'{arm}_seed{seed}/automatic_k_selection.json' for arm in ['C', 'S', 'J_SCRATCH'] for seed in SEEDS]
    missing = [str(p) for p in required if not p.is_file()]
    dep.require(not missing, 'Required experiment outputs are incomplete: ' + ', '.join(missing[:8]))
    ids = scales = None
    hashes = {}
    for arm in LABELS:
        for seed in SEEDS:
            directory = RESULTS / f'{arm}_seed{seed}'
            cfg = dep.read(directory / 'config.json')
            dep.require(cfg['fold'] == 6 and (not cfg['test_evaluated']), 'Wrong evaluation protocol')
            dep.require(dep.read(directory / 'complete.json')['signature'] == cfg['signature'], 'Completion signature')
            for path in directory.glob('*.json'):
                hashes[str(path)] = dep.digest(path)
            for path in directory.glob('K*/validation_metrics.json'):
                hashes[str(path)] = dep.digest(path)
            for path in directory.glob('K*/final_best.pt'):
                hashes[str(path)] = dep.digest(path)
            if arm == 'D':
                continue
            f = dep.pt(directory / 'frozen_profiles.pt')
            if ids is None:
                ids, scales = (f['feature_ids'].numpy(), f['scales'].numpy())
            dep.require(np.array_equal(f['feature_ids'], ids) and np.array_equal(f['scales'], scales), 'Candidate/P95 drift')
            for folder, pattern in [('points', '*.json'), ('solutions', '*.pt')]:
                for path in (directory / folder).glob(pattern):
                    hashes[str(path)] = dep.digest(path)
            hashes[str(directory / 'frozen_profiles.pt')] = dep.digest(directory / 'frozen_profiles.pt')
            for split in ['train', 'val']:
                for path in [directory / f'explicit_{split}.npy', *(directory / 'evidence' / split).iterdir()]:
                    if path.is_file():
                        hashes[str(path)] = dep.digest(path)
    for seed in SEEDS:
        for path in (RESULTS.parent / f'selection/C_seed{seed}/points').glob('*.json'):
            hashes[str(path)] = dep.digest(path)
    dep.write(OUT / 'source_manifest.json', dict(files=hashes, fold=6, seeds=SEEDS, test_used=False, candidate_count=len(ids), analysis_code_sha256=dep.digest(__file__), protocol_sha256='276c25daad5dd390f24491eb34b82f53972a410adac2943c09310adff36d392d', dependency_analysis_sha256=dep.digest(dep.__file__), experimental_bundles='Stage D checkpoints for ablation analysis, distinct from the released prediction model'))
    return ids

def collect_metrics():
    full = []
    curve = []
    refined = []
    convergence = []
    for arm, label in LABELS.items():
        for seed in SEEDS:
            directory = RESULTS / f'{arm}_seed{seed}'
            if arm == 'D':
                q = dep.read(directory / 'points/explicit.json')
            else:
                q = dep.read(directory / 'points/full_replay.json')
            full.append(dict(arm=arm, model=label, seed=seed, K=q['K'], **q['val']))
            for path in sorted((directory / 'points').glob('*.json')):
                point = dep.read(path)
                row = dict(arm=arm, model=label, seed=seed, point=point['point_id'], K=point['K'], eligible=point['eligible'], lambda_ratio=point['lambda_ratio'], val_ap=point['val']['ap'], val_auc=point['val']['auc'], val_brier=point['val']['brier'])
                convergence.append({**row, **point.get('solver', {}), 'converged': point['converged']})
                if point['eligible']:
                    curve.append(row)
            if arm not in ['C', 'S', 'J_SCRATCH']:
                continue
            selected = dep.read(directory / 'automatic_k_selection.json')['selected']
            for role, point in selected.items():
                path = directory / f"K{point['K']}/validation_metrics.json"
                dep.require(path.is_file(), 'Stage D is unfinished: ' + str(path))
                q = dep.read(path)
                refined.append(dict(arm=arm, model=label, seed=seed, role=role, K=q['K'], stage_c_ap=q['before_refinement']['ap'], **q['metrics']))
    save('full_validation.tsv', full)
    save('stage_c_frontiers.tsv', curve)
    save('convergence_audit.tsv', convergence)
    save('E2_stage_d.tsv', refined)
    full, curve, refined = (pd.DataFrame(full), pd.DataFrame(curve), pd.DataFrame(refined))
    save('full_validation_summary.tsv', summarize(full, ['arm', 'model'], ['K', 'ap', 'auc', 'brier']))
    save('E2_stage_d_summary.tsv', summarize(refined, ['arm', 'model', 'role'], ['K', 'stage_c_ap', 'ap', 'auc', 'brier']))
    differences = []
    for role in ['full', 'performance', 'efficiency']:
        frame = full if role == 'full' else refined[refined.role == role]
        paired = frame[frame.arm.isin(['C', 'S', 'J_SCRATCH'])].pivot(index='seed', columns='arm', values='ap')
        dep.require(list(paired.index) == list(SEEDS) and (not paired.isna().any().any()), 'Incomplete E2 paired seeds')
        for first, second in [('C', 'S'), ('C', 'J_SCRATCH'), ('S', 'J_SCRATCH')]:
            for seed, value in (paired[first] - paired[second]).items():
                differences.append(dict(seed=seed, stage=role, contrast=f'{first}_minus_{second}', ap_difference=value))
    diff = save('E2_paired_differences.tsv', differences)
    save('E2_paired_difference_summary.tsv', summarize(diff, ['stage', 'contrast'], ['ap_difference']))
    return (full, curve, refined)

def dependence(ids):
    directory = OUT / 'dependence'
    directory.mkdir(exist_ok=True)
    dep.E1 = RESULTS
    dep.SEEDS = SEEDS
    dep.MODELS = {'M': ('A', LABELS['A']), 'P': ('C', LABELS['C'])}
    dep.self_test()
    meta = dep.metadata()
    coverage = dep.coverage_stats(meta['train'], ids, directory)
    g, target = dep.discover(meta['train'], ids, coverage, directory)
    controls = dep.controls(g, target, 20, directory)
    dep.MODELS = {'M': ('A', LABELS['A']), 'P': ('C', LABELS['C']), **{arm: (arm, LABELS[arm]) for arm in DEPENDENCE_ARMS if arm != 'C'}}
    for seed in SEEDS:
        for key, (arm, name) in dep.MODELS.items():
            path = directory / f'train_screen_{key}_seed{seed}.tsv'
            if path.exists():
                continue
            logging.info('Additional TRAIN dependence %s seed=%d', name, seed)
            X, alpha, z, _ = dep.full_data(seed, key, 'train', meta['train'], ids)
            mu, cost, importance = dep.screen_columns(X, alpha, meta['train'].label.to_numpy(float), z, dep.group_weights(meta['train']))
            dep.save(path, dict(feature_id=ids, train_mean_contribution=mu, train_centered_bce_increase=cost, train_centered_abs_contribution=importance))
            del X
            gc.collect()
    deleted, endpoints = dep.deletion(meta['val'], ids, target, controls, directory)
    retention, paired = dep.retention(ids, target, controls, directory)
    ci = dep.bootstrap(meta['val'], endpoints, 500, directory)
    n = len(target)
    endpoint = deleted[(deleted.set == 'target') & (deleted['mode'] == 'mean_preserving') & (deleted.deleted_count == n)]
    save('joint_deletion_endpoint.tsv', endpoint)
    save('joint_deletion_endpoint_summary.tsv', endpoint.groupby('model')[['ap_loss', 'bce_change', 'group_bce_change']].agg(['mean', 'std']).reset_index())
    matched = []
    frame = retention[(retention.set == 'target') & (retention.K >= 100)]
    for seed in SEEDS:
        reference = frame[(frame.seed == seed) & (frame.model == LABELS['A'])].reset_index(drop=True)
        for arm in DEPENDENCE_ARMS:
            other = frame[(frame.seed == seed) & (frame.model == LABELS[arm])].reset_index(drop=True)
            delta = abs(reference.K.to_numpy()[:, None] - other.K.to_numpy()[None, :])
            relative = delta / np.maximum(reference.K.to_numpy()[:, None], other.K.to_numpy()[None, :])
            allowed = (delta <= 0.01 * len(ids)) & (relative <= 0.1)
            na, nb = delta.shape
            cost = np.full((na + nb, na + nb), 100000.0)
            cost[:na, :nb] = np.where(allowed, delta, 1000000000.0)
            cost[na:, nb:] = 0.0
            rr, cc = linear_sum_assignment(cost)
            for i, j in zip(rr, cc):
                if i >= na or j >= nb or (not allowed[i, j]):
                    continue
                x, y = (reference.iloc[i], other.iloc[j])
                matched.append(dict(seed=seed, arm=arm, model=LABELS[arm], A_K=int(x.K), K=int(y.K), A_point=x.point, point=y.point, A_retained=int(x.retained_count), retained=int(y.retained_count), retention_difference=y.set_retention - x.set_retention, set_size=n, mean_K=(x.K + y.K) / 2))
    save('near_K_retention_all_variants.tsv', matched)
    m = pd.DataFrame(matched)
    display = []
    for (seed, arm), group in m.groupby(['seed', 'arm']):
        row = group.assign(distance_to_600=abs(group.mean_K - 600)).sort_values(['distance_to_600', 'mean_K']).iloc[0]
        display.append(row.to_dict())
    save('retention_near600_display.tsv', display)
    dep.write(OUT / 'dependence_summary.json', dict(target_count=n, near_zero_80pct=int(target.near_zero_80pct.sum()), ci_A_minus_C_ap_loss=ci, common_target='Training-selected A-versus-C target set shared across branch variants', deletion_unit='latent features jointly removed without refitting', candidate_count=len(ids)))
    return meta

def latent39522(ids, meta):
    index = np.flatnonzero(ids == 39522)
    dep.require(len(index) == 1, 'ID39522 absent from fold candidates')
    index = int(index[0])
    rows = []
    for arm in ['C', 'H3', 'AROMATIC_ADD', 'AROMATIC_SPLIT']:
        dep.MODELS = {'X': (arm, LABELS[arm])}
        for seed in SEEDS:
            X, alpha, z, _ = dep.full_data(seed, 'X', 'train', meta['train'], ids)
            mean_train = float(X[:, index].mean() * alpha[index])
            del X
            X, alpha, z, _ = dep.full_data(seed, 'X', 'val', meta['val'], ids)
            contribution = X[:, index].toarray().ravel() * alpha[index]
            total_abs = np.asarray(abs(X).dot(alpha)).ravel()
            metric = dep.metric(meta['val'].label.to_numpy(float), z, dep.group_weights(meta['val']))
            deleted = dep.metric(meta['val'].label.to_numpy(float), z - contribution + mean_train, dep.group_weights(meta['val']))
            rows.append(dict(arm=arm, model=LABELS[arm], seed=seed, feature_id=39522, mean_abs_contribution=float(abs(contribution).mean()), contribution_share=float(abs(contribution).sum() / total_abs.sum()), val_ap=metric['ap'], deleted_val_ap=deleted['ap'], deletion_ap_loss=metric['ap'] - deleted['ap'], deletion_bce_increase=deleted['bce'] - metric['bce'], group_deletion_bce_increase=deleted['group_bce'] - metric['group_bce']))
            del X
            gc.collect()
    save('ID39522_contribution_and_dependency.tsv', rows)

def collect_e4():
    rows = []
    for seed in SEEDS:
        for path in sorted((RESULTS.parent / f'selection/C_seed{seed}/points').glob('*.json')):
            p = dep.read(path)
            rows.append(dict(seed=seed, point=p['point_id'], method=p['method'], family=p['family'], K=p['K'], eligible=p['eligible'], refitted=p['point_id'].endswith('_refit'), budget_K=p.get('selected_support_K', p['K']), random_repeat=p.get('random_repeat'), val_ap=p['val_ap'], val_auc=p['val_auc'], val_brier=p['val_brier'], reference=p['reference_point']))
    save('E4_all_points.tsv', rows)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(4)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(OUT / 'analysis.log')])
    ids = verify_all()
    collect_metrics()
    meta = dependence(ids)
    latent39522(ids, meta)
    collect_e4()
    dep.write(OUT / 'analysis_complete.json', dict(fold=6, seeds=SEEDS, test_used=False, source_manifest_sha256=dep.digest(OUT / 'source_manifest.json')))
if __name__ == '__main__':
    main()
