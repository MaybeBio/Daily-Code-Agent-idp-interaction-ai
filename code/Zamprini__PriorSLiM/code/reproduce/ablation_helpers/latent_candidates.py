import argparse
import gc
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score
import latent_deletion as a
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_all_candidates'
PLAN = HERE / 'PROTOCOL_ALL_CANDIDATES.md'
PARTITIONS = {'discovery': range(1, 6), 'followup': range(6, 11), 'all10': range(1, 11)}

def safe_ratio(x, y):
    return float(x / y) if y > 0 else float('nan')

def summarize_universe():
    raw = a.join('g_screen.tsv') if hasattr(a, 'join') else pd.concat([a.table(a.OUT / f'fold_{i:02d}/g_screen.tsv') for i in range(1, 11)], ignore_index=True)
    original = a.table(a.OUT / 'discovery_g_candidates.tsv').set_index('feature_id')
    ids = sorted(original.index[original.eligible.astype(bool)].tolist())
    a.t.require(len(ids) == 267, 'Discovery feature count differs from the expected 267')
    records = []
    cols = ['M_center_bce', 'B_center_bce', 'center_dependency_drop', 'M_raw_bce', 'B_raw_bce', 'M_center_abs', 'B_center_abs', 'M_abs_contribution', 'B_abs_contribution', 'train_active_fraction', 'val_active_fraction', 'retained_D']
    for partition, folds in PARTITIONS.items():
        for fid, g in raw[raw.fold.isin(folds)].groupby('feature_id'):
            item = dict(partition=partition, feature_id=int(fid), available_folds=len(g), broad_folds=int(g.broad.sum()), discovery_eligible=int(fid) in ids, discovery_broad=bool(fid in original.index and original.loc[fid, 'broad_folds'] >= 3), positive_M_folds=int((g.M_center_bce > 0).sum()), declining_folds=int(((g.M_center_bce > 0) & (g.center_dependency_drop > 0)).sum()), contribution_declining_folds=int((g.M_center_abs > g.B_center_abs).sum()))
            item.update({c: float(g[c].mean()) for c in cols})
            item['relative_dependency_drop'] = 1 - safe_ratio(item['B_center_bce'], item['M_center_bce'])
            item['absolute_residual_ratio'] = safe_ratio(abs(item['B_center_bce']), item['M_center_bce'])
            item['contribution_ratio'] = safe_ratio(item['B_center_abs'], item['M_center_abs'])
            item['retained_D_count'] = int(g.retained_D.sum())
            for pct in [50, 80, 90]:
                item[f'reduction_ge_{pct}'] = bool(item['M_center_bce'] >= 1e-05 and item['relative_dependency_drop'] >= pct / 100)
            item['near_zero_residual_80'] = bool(item['M_center_bce'] >= 1e-05 and item['absolute_residual_ratio'] <= 0.2)
            item['stringent_descriptive'] = bool(item['discovery_eligible'] and item['available_folds'] == len(folds) and item['near_zero_residual_80'] and (item['contribution_ratio'] <= 0.5) and (item['declining_folds'] >= (4 if len(folds) == 5 else 8)))
            records.append(item)
    summary = pd.DataFrame(records)
    a.save(raw, OUT / 'all_candidate_fold_metrics.tsv')
    a.save(summary, OUT / 'all_candidate_summary.tsv')
    a.save(summary[summary.discovery_eligible], OUT / 'all_267_candidate_summary.tsv')
    a.t.write_json(OUT / 'scope.json', dict(all_union_ids=int(raw.feature_id.nunique()), discovery_union_ids=len(original), discovery_broad_ids=int((original.broad_folds >= 3).sum()), detailed_ids=ids, n_detailed=len(ids), top_n_limit=None, protocol_sha256='6db1d7cb27d5d1cca434224a70c89ed8df11f7d7483cf9528c32e8fd16adc9c8'))
    return ids

def component_aggregator(meta):
    code, names = pd.factorize(meta.component_id, sort=True)
    count = np.bincount(code)
    mat = sp.csr_matrix((1.0 / count[code], (code, np.arange(len(code)))), shape=(len(names), len(code)))
    return (mat, np.asarray(names))

def correlation_groups(tx, targets, ids):
    n = tx.shape[0]
    mean = np.asarray(tx.sum(0)).ravel() / n
    sd = np.sqrt(np.maximum(np.asarray(tx.power(2).sum(0)).ravel() / n - mean ** 2, 0))
    cross = (tx[:, targets].T @ tx).toarray() / n - mean[targets, None] * mean[None, :]
    den = sd[targets, None] * sd[None, :]
    corr = np.divide(cross, den, out=np.zeros_like(cross), where=den > 1e-12)
    groups = []
    for row, j in zip(corr, targets):
        near = sorted((i for i in np.flatnonzero(row >= 0.7) if i != j), key=lambda i: (-row[i], ids[i]))[:19]
        groups.append((np.array([j, *near]), row))
    return groups

def weighted_corr(x, y, w):
    xx = x - w @ x
    yy = y - w @ y
    den = np.sqrt(w @ xx ** 2 * (w @ yy ** 2))
    return float(w @ (xx * yy) / den) if den > 1e-15 else float('nan')

def positive_type_map():
    manifest = a.table(a.ROOT / 'data/new_all_positive_augmentation_manifest.tsv')
    keys = ['window_start', 'window_end', 'slim_start', 'slim_end']
    mapping = {}
    for r in manifest.itertuples(index=False):
        k = (r.uniprot, *(int(float(getattr(r, c))) for c in keys))
        a.t.require(k not in mapping or mapping[k] == r.augmentation_type, 'Conflicting augmentation provenance')
        mapping[k] = r.augmentation_type
    combined = a.table(a.ROOT / 'train/combined_metadata.tsv').set_index('window_id')
    return (mapping, combined)

def detailed_fold(fold, ids, provenance):
    d = OUT / f'fold_{fold:02d}'
    d.mkdir(parents=True, exist_ok=True)
    signature = a.t.identity(dict(ids=ids, code=a.t.digest(__file__), protocol='6db1d7cb27d5d1cca434224a70c89ed8df11f7d7483cf9528c32e8fd16adc9c8'))
    if (d / 'complete.json').exists():
        a.t.require(a.t.read(d / 'complete.json')['signature'] == signature, 'Detail code/plan changed')
        logging.info('Reuse exhaustive fold %02d', fold)
        return
    f = a.load_fold(fold)
    previous = a.t.read(a.OUT / f'fold_{fold:02d}/prepared.json')
    a.t.require(previous['inputs'] == f['fingerprint'], 'Inputs differ from the saved contribution calculation')
    src = a.OUT / f'fold_{fold:02d}'
    lookup = {int(fid): j for j, fid in enumerate(f['ids'])}
    a.t.require(all((fid in lookup for fid in ids)), 'A requested feature ID is missing from the fold feature set')
    targets = np.array([lookup[fid] for fid in ids])
    stats = np.load(src / 'train_stats.npz')
    cm, cb = [sp.load_npz(src / f'{name}.npz').tocsc() for name in ['M', 'B']]
    active = sp.load_npz(src / 'activation.npz').tocsc()
    tx = sp.load_npz(src / 'train_activation_sum.npz').astype(float)
    logging.info('Fold %02d: TRAIN correlated groups for all %d eligible IDs', fold, len(ids))
    groups = correlation_groups(tx, targets, f['ids'])
    del tx
    gc.collect()
    y = f['y']
    sign = 2 * y - 1
    w = a.balanced(f['meta'])
    agg, components = component_aggregator(f['meta'])
    e_mean = float(np.load(f['source'] / 'explicit_train.npy').mean())
    e = f['expb'] - e_mean
    remove_e = a.delta(y, f['zb'], e)
    mtotal = np.asarray(cm.sum(1)).ravel() - stats['mu_m'].sum()
    btotal = np.asarray(cb.sum(1)).ravel() - stats['mu_b'].sum()
    source_type = f['meta'].negative_type.to_numpy(str).copy()
    mapping, combined = provenance
    for n in np.flatnonzero(y == 1):
        r = combined.loc[f['meta'].window_id.iloc[n]]
        key = (r.uniprot, *(int(float(r[k])) for k in ['window_start', 'window_end', 'slim_start', 'slim_end']))
        source_type[n] = 'positive:' + mapping[key]
    source_type = np.array(['positive:' + mapping[combined.loc[wid].uniprot, *(int(float(combined.loc[wid, k])) for k in ['window_start', 'window_end', 'slim_start', 'slim_end'])] if label == 1 else nt for wid, label, nt in zip(f['meta'].window_id, y, f['meta'].negative_type)], dtype=object)
    rows, group_rows, subtypes, same_rows, group_def = ([], [], [], [], [])
    effects = np.zeros((len(components), len(ids), 4), dtype=float)
    ap0 = {name: average_precision_score(y, z) for name, z in [('M', f['zm']), ('B', f['zb'])]}
    auc0 = {name: roc_auc_score(y, z) for name, z in [('M', f['zm']), ('B', f['zb'])]}
    detail_checks = []
    for column, (fid, j, (group, cor)) in enumerate(zip(ids, targets, groups)):
        m = cm[:, j].toarray().ravel() - stats['mu_m'][j]
        b = cb[:, j].toarray().ravel() - stats['mu_b'][j]
        dm = a.delta(y, f['zm'], m)
        db = a.delta(y, f['zb'], b)
        for name, z, c, effect in [('M', f['zm'], m, dm), ('B', f['zb'], b, db)]:
            rows.append(dict(fold=fold, feature_id=fid, model=name, deletion_bce=float(w @ effect), ap_drop=float(ap0[name] - average_precision_score(y, z - c)), auc_drop=float(auc0[name] - roc_auc_score(y, z - c)), mean_abs_centered_contribution=float(w @ abs(c))))
        effects[:, column, 0] = agg @ dm
        effects[:, column, 1] = agg @ db
        effects[:, column, 2] = agg @ abs(m)
        effects[:, column, 3] = agg @ abs(b)
        act = active[:, j].toarray().ravel() > 0
        helpful_m = act & (dm > 0) & (sign * m > 0)
        lost = m - b
        diminished = helpful_m & (sign * lost > 0)
        others = btotal - b - (mtotal - m)
        check = f['zb'] - (e_mean + stats['mu_b'].sum()) - (f['zm'] - (f['intercept'] + stats['mu_m'].sum()))
        err = float(np.max(abs(check - (e - lost + others))))
        a.t.require(err <= 1e-05, 'Centered prediction decomposition failed')
        detail_checks.append(err)
        ww = a.balanced(f['meta'], diminished)
        same_rows.append(dict(fold=fold, feature_id=fid, active_windows=int(act.sum()), M_helpful_windows=int(helpful_m.sum()), diminished_helpful_windows=int(diminished.sum()), diminished_helpful_proteins=f['meta'].loc[diminished, 'base_protein_id'].nunique(), explicit_support_windows=int((diminished & (remove_e > 0)).sum()), explicit_support_fraction=safe_ratio(int((diminished & (remove_e > 0)).sum()), int(diminished.sum())), other_latents_aligned_fraction=safe_ratio(int((diminished & (sign * others > 0)).sum()), int(diminished.sum())), explicit_support_mean_bce=float(ww @ remove_e) if diminished.any() else np.nan, lost_explicit_correlation=weighted_corr(lost, e, w), lost_other_latents_correlation=weighted_corr(lost, others, w), B_latent_deletion_without_explicit_bce=float(w @ a.delta(y, f['zb'] - e, b))))
        for st in sorted(set(source_type)):
            mask = source_type == st
            sw = a.balanced(f['meta'], mask)
            subtypes.append(dict(fold=fold, feature_id=fid, sample_type=st, windows=int(mask.sum()), active_windows=int((act & mask).sum()), M_bce=float(sw @ dm), B_bce=float(sw @ db), dependency_drop=float(sw @ (dm - db)), M_abs=float(sw @ abs(m)), B_abs=float(sw @ abs(b)), diminished_helpful_windows=int((diminished & mask).sum()), explicit_support_windows=int((diminished & mask & (remove_e > 0)).sum())))
        group_def.append(dict(fold=fold, feature_id=fid, group_size=len(group), group_ids=';'.join(map(str, f['ids'][group])), correlations=';'.join((f'{cor[q]:.6f}' for q in group))))
        for name, matrix, z, mus in [('M', cm, f['zm'], stats['mu_m']), ('B', cb, f['zb'], stats['mu_b'])]:
            c = np.asarray(matrix[:, group].sum(1)).ravel() - mus[group].sum()
            group_rows.append(dict(fold=fold, feature_id=fid, model=name, group_size=len(group), deletion_bce=float(w @ a.delta(y, z, c)), ap_drop=float(ap0[name] - average_precision_score(y, z - c)), mean_abs_centered_contribution=float(w @ abs(c))))
        if (column + 1) % 75 == 0:
            logging.info('Fold %02d: %d/%d candidates checked', fold, column + 1, len(ids))
    for name, result in [('single_deletions', rows), ('group_deletions', group_rows), ('sample_types', subtypes), ('same_window_support', same_rows), ('train_groups', group_def)]:
        a.save(result, d / f'{name}.tsv')
    np.savez_compressed(d / 'component_effects.npz', effects=effects, component_ids=components, feature_ids=np.array(ids), metric_names=['M_bce', 'B_bce', 'M_abs', 'B_abs'])
    a.t.write_json(d / 'complete.json', dict(signature=signature, fold=fold, n_ids=len(ids), max_decomposition_error=max(detail_checks), original_replay=previous['replay'], inputs=f['fingerprint'], explicit_train_sha256=a.t.digest(f['source'] / 'explicit_train.npy'), cached_matrices_sha256={name: a.t.digest(src / f'{name}.npz') for name in ['M', 'B', 'activation', 'train_stats', 'train_activation_sum']}))
    logging.info('Fold %02d complete', fold)

def bh(p):
    p = np.asarray(p)
    n = len(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(q, 0, 1)
    return result

def bootstrap(ids):
    rng = np.random.default_rng(20261007)
    draws = 10000
    b = np.zeros((draws, len(ids), 4))
    observed = np.zeros((len(ids), 4))
    total_groups = 0
    for fold in range(6, 11):
        with np.load(OUT / f'fold_{fold:02d}/component_effects.npz') as f:
            a.t.require(np.array_equal(f['feature_ids'], ids), 'Bootstrap ID mismatch')
            e = f['effects'].copy()
        n = len(e)
        total_groups += n
        observed += e.mean(axis=0) / 5
        for start in range(0, draws, 250):
            size = min(250, draws - start)
            counts = rng.multinomial(n, np.full(n, 1 / n), size=size).astype(float) / n
            b[start:start + size] += (counts @ e.reshape(n, -1)).reshape(size, len(ids), 4) / 5
    dm, db = (observed[:, 0], observed[:, 1])
    drop = dm - db
    boot_drop = b[:, :, 0] - b[:, :, 1]
    p_m = (1 + (b[:, :, 0] - dm >= dm).sum(0)) / (draws + 1)
    p_drop = (1 + (boot_drop - drop >= drop).sum(0)) / (draws + 1)
    q = bh(np.maximum(p_m, p_drop))
    rows = []
    for j, fid in enumerate(ids):
        row = dict(feature_id=fid, M_bce=dm[j], B_bce=db[j], dependency_drop=drop[j], p_positive_M=p_m[j], p_positive_drop=p_drop[j], conditional_q_joint=q[j], homology_components=total_groups, bootstrap_draws=draws)
        for name, values in [('M_bce', b[:, j, 0]), ('B_bce', b[:, j, 1]), ('dependency_drop', boot_drop[:, j])]:
            row[name + '_lo'], row[name + '_hi'] = np.quantile(values, [0.025, 0.975])
        valid = b[:, j, 0] > 1e-08
        row['positive_denominator_bootstrap_fraction'] = float(valid.mean())
        if valid.any():
            r = 1 - b[valid, j, 1] / b[valid, j, 0]
            row['relative_drop_lo'], row['relative_drop_hi'] = np.quantile(r, [0.025, 0.975])
        rows.append(row)
    a.save(rows, OUT / 'followup_conditional_bootstrap.tsv')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--phase', choices=['all', 'overview', 'details', 'report'], default='all')
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', stream=sys.stdout)
    a.t.torch.set_num_threads(4)
    OUT.mkdir(exist_ok=True)
    ids = summarize_universe()
    if args.phase in ['all', 'details']:
        provenance = positive_type_map()
        for fold in range(1, 11):
            detailed_fold(fold, ids, provenance)
            gc.collect()
        bootstrap(ids)
    if args.phase in ['all', 'report']:
        pass
if __name__ == '__main__':
    main()
