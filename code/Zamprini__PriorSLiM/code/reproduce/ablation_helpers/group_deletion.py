from __future__ import annotations
import argparse
import gc
import hashlib
import json
import logging
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.special import expit, logit
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import average_precision_score, roc_auc_score
import torch
HERE = Path(__file__).resolve().parent
ROOT = Path(__import__('os').environ['PRIORSLIM_DATA_ROOT'])
E1 = HERE.parent / 'J_WARM/E1/results'
SEEDS = (2026, 2027, 2028, 2029, 2030)
MODELS = {'M': ('A', 'MotifAE-only'), 'P': ('C', 'PriorSLiM')}

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def read(path):
    return json.loads(Path(path).read_text())

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def native(x):
    if isinstance(x, dict):
        return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [native(v) for v in x]
    if isinstance(x, np.ndarray):
        return native(x.tolist())
    if isinstance(x, np.generic):
        return native(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def write(path, x):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(native(x), ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def table(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, float_precision='round_trip')

def save(path, x):
    pd.DataFrame(x).to_csv(path, sep='\t', index=False)

def pt(path):
    return torch.load(path, map_location='cpu', weights_only=True)

def loss(y, z):
    return np.logaddexp(0.0, z) - y * z

def group_weights(meta):
    g = meta.component_id
    return 1.0 / g.map(g.value_counts()).to_numpy(float) / g.nunique()

def evidence(path, weights):
    info = read(path / 'complete.json')
    wh = hashlib.sha256(np.ascontiguousarray(weights).view(np.uint8)).hexdigest()
    require(info['weight_hash'] == wh, 'Evidence/profile checksum mismatch')
    ptr = np.load(path / 'indptr.npy', mmap_mode='r')
    col = np.memmap(path / 'columns.bin', dtype=np.int32, mode='r')
    val = np.memmap(path / 'values.bin', dtype=np.float64, mode='r')
    require(len(ptr) == info['n'] + 1 and ptr[-1] == len(col) == len(val), 'CSR length mismatch')
    require(np.isfinite(val).all() and col.min() >= 0 and (col.max() < len(weights)), 'Bad evidence')
    return sp.csr_matrix((val, col, ptr), shape=(info['n'], len(weights)), copy=False)

def metadata():
    mapping = table(ROOT / 'new_data/cv10_homology_seed2026/protein_fold_assignments.tsv')
    require(mapping.Primary_Acc.is_unique, 'Ambiguous component mapping')
    lookup = mapping.set_index('Primary_Acc').component_id
    result = {}
    for split in ('train', 'val'):
        x = table(E1 / '_prepared' / f'{split}_metadata.tsv')
        require(x.window_id.is_unique and set(x.split) == {split}, 'Wrong/duplicate split IDs')
        x['component_id'] = x.base_protein_id.map(lookup)
        require(x.component_id.notna().all(), 'Missing homology group')
        result[split] = x
    for key in ('window_id', 'base_protein_id', 'component_id'):
        require(not set(result['train'][key]) & set(result['val'][key]), 'Split overlap: ' + key)
    return result

def coverage_stats(meta, ids, out):
    path = out / 'train_coverage.tsv'
    if path.exists():
        x = table(path)
        require(np.array_equal(x.feature_id, ids), 'Coverage IDs changed')
        return x
    src = E1 / '_prepared/train'
    ptr = np.load(src / 'indptr.npy', mmap_mode='r')
    cols = np.memmap(src / 'columns.bin', dtype=np.int32, mode='r')
    vals = np.memmap(src / 'values.bin', dtype=np.float32, mode='r')
    counts = np.zeros(len(ids), np.int64)
    codes, proteins = pd.factorize(meta.base_protein_id)
    seen = np.zeros((len(proteins), len(ids)), bool)
    for start in range(0, len(meta), 512):
        stop = min(len(meta), start + 512)
        lo, hi = (ptr[start], ptr[stop])
        c = np.asarray(cols[lo:hi]) // 29
        v = np.asarray(vals[lo:hi])
        require(np.isfinite(v).all() and (v > 0).all(), 'Invalid sparse activation')
        owners = np.repeat(np.arange(stop - start), np.diff(ptr[start:stop + 1]))
        pairs = np.unique(owners * len(ids) + c)
        counts += np.bincount(pairs % len(ids), minlength=len(ids))
        seen[codes[start + pairs // len(ids)], pairs % len(ids)] = True
    fraction = counts / len(meta)
    x = pd.DataFrame(dict(feature_id=ids, train_active_windows=counts, train_active_fraction=fraction, train_active_proteins=seen.sum(0), broad=(fraction >= np.quantile(fraction, 0.75)) & (seen.sum(0) >= 10)))
    save(path, x)
    return x

def full_data(seed, model, split, meta, ids):
    src = E1 / f'{MODELS[model][0]}_seed{seed}'
    frozen = pt(src / 'frozen_profiles.pt')
    require(np.array_equal(frozen['feature_ids'].numpy(), ids), 'Feature identity mismatch')
    weights = frozen['weights'].numpy()
    alpha = frozen['alpha_init'].numpy()
    require(np.allclose(np.linalg.norm(weights, axis=1), alpha, atol=1e-12, rtol=0), 'Alpha/profile mismatch')
    X = evidence(src / 'evidence' / split, weights)
    offset = np.load(src / f'explicit_{split}.npy')
    z = offset + X @ alpha
    require(len(z) == len(meta) and np.isfinite(z).all(), 'Bad replay length/value')
    metrics = read(src / 'full_metrics.json')['metrics'][split]['full']
    ap = average_precision_score(meta.label, expit(z))
    require(abs(ap - metrics['ap']) < 1e-09, 'Full AP replay mismatch')
    error = 0.0
    if split == 'val':
        with np.load(src / 'predictions/full_replay.npz') as saved:
            require(np.array_equal(saved['window_ids'].astype(str), meta.window_id), 'Window identity mismatch')
            require(np.array_equal(saved['labels'], meta.label), 'Label mismatch')
            error = float(np.max(abs(z - saved['logits'])))
            require(error <= 1e-08, f'Logit replay failed: {error}')
    return (X, alpha, z, error)

def screen_columns(X, alpha, y, z, group_w):
    X = X.tocsc()
    mu = np.asarray(X.mean(0)).ravel() * alpha
    base = loss(y, z)
    dm = np.empty(len(alpha))
    imp = np.empty(len(alpha))
    for start in range(0, len(alpha), 32):
        end = min(len(alpha), start + 32)
        c = X[:, start:end].toarray() * alpha[None, start:end]
        centered = c - mu[None, start:end]
        change = loss(y[:, None], z[:, None] - centered) - base[:, None]
        dm[start:end] = group_w @ change
        imp[start:end] = group_w @ abs(centered)
    return (mu, dm, imp)

def discover(meta, ids, coverage, out):
    for seed in SEEDS:
        for model in MODELS:
            path = out / f'train_screen_{model}_seed{seed}.tsv'
            if path.exists():
                continue
            logging.info('TRAIN-only dependence screen %s seed=%d', MODELS[model][1], seed)
            X, alpha, z, _ = full_data(seed, model, 'train', meta, ids)
            mu, dep, imp = screen_columns(X, alpha, meta.label.to_numpy(float), z, group_weights(meta))
            save(path, dict(feature_id=ids, train_mean_contribution=mu, train_centered_bce_increase=dep, train_centered_abs_contribution=imp))
            del X
            gc.collect()
    rows = []
    for seed in SEEDS:
        m = table(out / f'train_screen_M_seed{seed}.tsv')
        p = table(out / f'train_screen_P_seed{seed}.tsv')
        for frame in (m, p):
            require(np.array_equal(frame.feature_id, ids), 'Screen ID mismatch')
        rows.append(pd.DataFrame(dict(feature_id=ids, seed=seed, M_dependency=m.train_centered_bce_increase, P_dependency=p.train_centered_bce_increase, M_importance=m.train_centered_abs_contribution, positive_drop=(m.train_centered_bce_increase > 0) & (m.train_centered_bce_increase > p.train_centered_bce_increase))))
    per = pd.concat(rows)
    save(out / 'train_screen_all_seeds.tsv', per)
    g = per.groupby('feature_id').agg(M_dependency=('M_dependency', 'mean'), P_dependency=('P_dependency', 'mean'), M_importance=('M_importance', 'mean'), positive_drop_seeds=('positive_drop', 'sum')).reset_index()
    g = coverage.merge(g, on='feature_id', validate='one_to_one')
    g['dependency_drop'] = g.M_dependency - g.P_dependency
    g['target'] = g.broad & (g.M_dependency >= 1e-05) & (g.positive_drop_seeds >= 4) & (abs(g.P_dependency) <= 0.5 * g.M_dependency)
    g['near_zero_80pct'] = g.target & (abs(g.P_dependency) <= 0.2 * g.M_dependency)
    target = g[g.target].sort_values(['dependency_drop', 'feature_id'], ascending=[False, True])
    require(len(target) > 0, 'No TRAIN-defined target')
    save(out / 'candidate_screen.tsv', g)
    save(out / 'locked_target_features.tsv', target.assign(deletion_rank=np.arange(1, len(target) + 1)))
    logging.info('Locked %d TRAIN-only target IDs (%d meet 80%% band)', len(target), target.near_zero_80pct.sum())
    return (g, target)

def controls(g, target, repeats, out):
    pool = g[~g.target & (g.M_dependency > 0) & (g.M_importance > 0)].copy()
    require(len(pool) >= len(target), 'Insufficient non-target controls')

    def cov(x):
        return np.column_stack([logit(np.clip(x.train_active_fraction, 1e-06, 1 - 1e-06)), np.log10(x.M_dependency.clip(lower=1e-12)), np.log10(x.M_importance.clip(lower=1e-12))])
    both = np.concatenate([cov(target), cov(pool)])
    sd = both.std(0)
    sd[sd < 1e-12] = 1.0
    tx, px = (cov(target) / sd, cov(pool) / sd)
    distance = ((tx[:, None, :] - px[None, :, :]) ** 2).sum(-1)
    rng = np.random.default_rng(20261008)
    rows = []
    diagnostics = []
    sets = []
    for repeat in range(repeats):
        rr, cc = linear_sum_assignment(distance + rng.uniform(0, 0.25, distance.shape))
        require(np.array_equal(rr, np.arange(len(target))), 'Incomplete matching')
        chosen = pool.iloc[cc]
        sets.append(chosen.feature_id.to_numpy(int))
        for rank, (fid, cid, dist) in enumerate(zip(target.feature_id, chosen.feature_id, np.sqrt(distance[rr, cc])), 1):
            rows.append(dict(repeat=repeat, deletion_rank=rank, target_feature_id=int(fid), control_feature_id=int(cid), distance=float(dist)))
        diagnostics.append(dict(repeat=repeat, mean_distance=float(np.sqrt(distance[rr, cc]).mean()), max_distance=float(np.sqrt(distance[rr, cc]).max()), **{name: float(value) for name, value in zip(('smd_logit_coverage', 'smd_log_dependency', 'smd_log_importance'), (tx - px[cc]).mean(0))}))
    save(out / 'matched_controls.tsv', rows)
    save(out / 'matching_diagnostics.tsv', diagnostics)
    return sets

def metric(y, z, w):
    p = expit(z)
    return dict(ap=float(average_precision_score(y, p)), auc=float(roc_auc_score(y, p)), brier=float(np.mean((p - y) ** 2)), bce=float(loss(y, z).mean()), group_ap=float(average_precision_score(y, p, sample_weight=w)), group_bce=float(w @ loss(y, z)))

def deletion(meta, ids, target, controls_, out):
    lookup = {int(v): i for i, v in enumerate(ids)}
    sets = [target.feature_id.to_numpy(int), *controls_]
    positions = [np.array([lookup[int(v)] for v in s]) for s in sets]
    n_target = len(target)
    counts = sorted({v for v in [0, 1, 5, 10, 20, 50, 100, 200, 500, n_target] if v <= n_target})
    y = meta.label.to_numpy(float)
    w = group_weights(meta)
    result = []
    replay = []
    endpoints = {}
    for seed in SEEDS:
        endpoints[seed] = {}
        for model in MODELS:
            logging.info('Validation joint deletion %s seed=%d', MODELS[model][1], seed)
            X, alpha, z, error = full_data(seed, model, 'val', meta, ids)
            baseline = metric(y, z, w)
            replay.append(dict(seed=seed, model=MODELS[model][1], max_logit_error=error, **baseline))
            mu = table(out / f'train_screen_{model}_seed{seed}.tsv').train_mean_contribution.to_numpy()
            X = X.tocsc()
            z_target = None
            for repeat, pos in enumerate(positions):
                ordered = X[:, pos].toarray() * alpha[pos][None, :]
                cumulative = np.cumsum(ordered, axis=1)
                cumulative_mu = np.cumsum(mu[pos])
                for count in counts:
                    raw = cumulative[:, count - 1] if count else np.zeros(len(y))
                    avg = cumulative_mu[count - 1] if count else 0.0
                    for mode in ('mean_preserving', 'zero'):
                        zz = z - raw + (avg if mode == 'mean_preserving' else 0.0)
                        metrics = metric(y, zz, w)
                        row = dict(seed=seed, model=MODELS[model][1], set='target' if repeat == 0 else 'matched_control', repeat=repeat - 1, mode=mode, deleted_count=count, **metrics)
                        row.update({f'{key}_change': metrics[key] - baseline[key] for key in metrics})
                        row['ap_loss'] = baseline['ap'] - metrics['ap']
                        result.append(row)
                        if repeat == 0 and count == n_target and (mode == 'mean_preserving'):
                            z_target = zz.copy()
            endpoints[seed][model] = (z.copy(), z_target)
            del X
            gc.collect()
    save(out / 'joint_deletion_curves.tsv', result)
    save(out / 'full_replay_checks.tsv', replay)
    np.savez_compressed(out / 'validation_target_endpoint_logits.npz', window_ids=meta.window_id.to_numpy(str), labels=y, **{f'{model}_{seed}_{mode}': z for seed, models in endpoints.items() for model, pair in models.items() for mode, z in zip(('full', 'deleted'), pair)})
    return (pd.DataFrame(result), endpoints)

def retention(ids, target, controls_, out):
    lookup = {int(v): i for i, v in enumerate(ids)}
    sets = [target.feature_id.to_numpy(int), *controls_]
    pos = [np.array([lookup[int(v)] for v in s]) for s in sets]
    result = []
    supports = {}
    excluded = []
    for seed in SEEDS:
        supports[seed] = {}
        for model in MODELS:
            src = E1 / f'{MODELS[model][0]}_seed{seed}'
            supports[seed][model] = []
            signature = read(src / 'config.json')['signature']
            for p in sorted((src / 'points').glob('*.json')):
                point = read(p)
                if point['point_id'] == 'full_replay':
                    alpha = pt(src / 'frozen_profiles.pt')['alpha_init'].numpy()
                else:
                    ck = pt(src / 'solutions' / f"{point['point_id']}.pt")
                    require(ck['signature'] == signature, 'Wrong Stage C signature')
                    require(bool(ck['solver']['converged']) == bool(point['eligible']), 'Convergence mismatch')
                    alpha = ck['alpha'].numpy()
                kept = alpha > 1e-10
                require(int(kept.sum()) == point['K'], 'Stage C support K mismatch')
                if not point['eligible']:
                    excluded.append(dict(seed=seed, model=MODELS[model][1], point=p.name, K=point['K']))
                    continue
                info = dict(seed=seed, model=MODELS[model][1], point=point['point_id'], K=int(kept.sum()), retained_fraction=float(kept.mean()), lambda_ratio=point['lambda_ratio'])
                target_retention = float(kept[pos[0]].mean())
                ctrl_retention = float(np.mean([kept[q].mean() for q in pos[1:]]))
                supports[seed][model].append(dict(**info, target_retention=target_retention, control_retention=ctrl_retention))
                for repeat, pidx in enumerate(pos):
                    result.append(dict(**info, set='target' if repeat == 0 else 'matched_control', repeat=repeat - 1, set_size=len(pidx), retained_count=int(kept[pidx].sum()), set_retention=float(kept[pidx].mean())))
    paired = []
    for seed, models in supports.items():
        ma, pb = (models['M'], models['P'])
        distance = np.abs(np.array([r['K'] for r in ma])[:, None] - np.array([r['K'] for r in pb])[None, :])
        n_m, n_p = distance.shape
        cost = np.full((n_m + n_p, n_m + n_p), 100000.0)
        cost[:n_m, :n_p] = np.where(distance <= 0.01 * len(ids), distance, 1000000000.0)
        cost[n_m:, n_p:] = 0.0
        ii, jj = linear_sum_assignment(cost)
        for i, j in zip(ii, jj):
            if i >= n_m or j >= n_p:
                continue
            if distance[i, j] > 0.01 * len(ids):
                continue
            a, b = (ma[i], pb[j])
            paired.append(dict(seed=seed, M_point=a['point'], P_point=b['point'], M_K=a['K'], P_K=b['K'], K_difference=b['K'] - a['K'], mean_retained_fraction=(a['K'] + b['K']) / (2 * len(ids)), M_target_retention=a['target_retention'], P_target_retention=b['target_retention'], P_minus_M_target_retention=b['target_retention'] - a['target_retention'], M_control_retention=a['control_retention'], P_control_retention=b['control_retention'], P_minus_M_control_retention=b['control_retention'] - a['control_retention']))
    save(out / 'stage_c_retention.tsv', result)
    save(out / 'near_matched_K_retention.tsv', paired)
    write(out / 'excluded_stage_c_points.json', excluded)
    return (pd.DataFrame(result), pd.DataFrame(paired))

def bootstrap(meta, endpoints, repeats, out):
    rng = np.random.default_rng(4282026)
    codes, groups = pd.factorize(meta.component_id)
    y = meta.label.to_numpy(float)
    rows = []
    for rep in range(repeats):
        counts = np.bincount(rng.integers(0, len(groups), len(groups)), minlength=len(groups))
        w = counts[codes]
        if len(np.unique(y[w > 0])) < 2:
            continue
        delta = []
        for seed in SEEDS:
            losses = {m: average_precision_score(y, expit(pair[0]), sample_weight=w) - average_precision_score(y, expit(pair[1]), sample_weight=w) for m, pair in endpoints[seed].items()}
            delta.append(losses['M'] - losses['P'])
        rows.append(dict(replicate=rep, mean_M_minus_P_AP_loss=float(np.mean(delta))))
    save(out / 'endpoint_component_bootstrap.tsv', rows)
    x = pd.DataFrame(rows).mean_M_minus_P_AP_loss
    return dict(replicates=len(x), homology_components=len(groups), q025=float(x.quantile(0.025)), q975=float(x.quantile(0.975)), conditional_on_fixed_models_and_train_selected_set=True)

def self_test():
    y = np.array([0.0, 1.0, 0.0, 1.0])
    z = np.array([-0.5, 0.2, 0.8, 1.0])
    c = np.array([[1.0, 0.0, 0.1], [0.0, 1.0, 0.2], [-1.0, 0.0, 0.3], [0.0, -1.0, 0.4]])
    mu = c.mean(0)
    w = np.ones(4) / 4
    _, d, _ = screen_columns(sp.csr_matrix(c), np.ones(3), y, z, w)
    expected = np.array([np.mean(loss(y, z - (c[:, j] - mu[j])) - loss(y, z)) for j in range(3)])
    require(np.allclose(d, expected, atol=1e-14), 'Single deletion implementation')
    all_deleted = z - (c[:, [0, 2]] - mu[[0, 2]]).sum(1)
    require(np.allclose(all_deleted, z - c[:, [0, 2]].sum(1) + mu[[0, 2]].sum()), 'Joint deletion algebra')
    require(not np.isclose(np.mean(loss(y, all_deleted) - loss(y, z)), d[[0, 2]].sum()), 'Nonadditivity test needs nonlinear case')
    return dict(single_deletion_matches_brute_force=True, joint_deletion_algebra=True, nonadditivity_preserved=True)
