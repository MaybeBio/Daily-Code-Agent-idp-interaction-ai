from pathlib import Path
import sys, logging, argparse, json, gc
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.special import expit, logit
from sklearn.metrics import average_precision_score, roc_auc_score
HERE = Path(__file__).resolve().parent
ROOT = Path(__import__('os').environ['PRIORSLIM_DATA_ROOT'])
import train_jwarm as t
MODELS = ROOT / 'train/outputs_jwarm_cv10/seed42'
SOLO = ROOT / 'prior_compression_experiments/MOTIFAE_ONLY_COMPARISON/CV10_TEST/results'
CV = ROOT / 'new_data/cv10_homology_seed2026'
OUT = HERE / 'results'

def table(p):
    return pd.read_csv(p, sep='\t', keep_default_na=False, float_precision='round_trip')

def save(x, p):
    pd.DataFrame(x).to_csv(p, sep='\t', index=False)

def loss(y, z):
    return np.logaddexp(0, z) - y * z

def balanced(meta, mask=None):
    mask = np.ones(len(meta), bool) if mask is None else mask
    w = np.zeros(len(meta))
    g = meta.loc[mask, 'component_id']
    counts = g.value_counts()
    if len(g):
        w[np.flatnonzero(mask)] = 1 / g.map(counts).to_numpy(float) / len(counts)
    return w

def delta(y, z, c):
    return loss(y, z - c) - loss(y, z)

def pred(p, meta):
    with np.load(p) as a:
        t.require(np.array_equal(a['window_ids'].astype(str), meta.window_id), 'Prediction window IDs changed')
        t.require(np.array_equal(a['labels'], meta.label.to_numpy(float)), 'Prediction labels changed')
        return (a['logits'].copy(), float(a['threshold']) if 'threshold' in a else None)

def load_fold(fold):
    d = MODELS / f'fold_{fold:02d}'
    s = SOLO / f'fold_{fold:02d}'
    meta = table(CV / f'fold_{fold:02d}/val_metadata.tsv')
    train = table(CV / f'fold_{fold:02d}/train_metadata.tsv')
    t.require(meta.window_id.is_unique and (not set(meta.window_id) & set(train.window_id)), 'Split leakage')
    t.require(set(meta['split']) == {'val'} and set(train['split']) == {'train'}, 'Wrong roles')
    t.require(not set(meta.component_id) & set(train.component_id), 'Component leakage')
    lock = t.read(s / 'motifae_locked.json')
    ck = t.load_pt(s / 'motifae_only/best.pt')
    t.require(t.digest(s / 'motifae_only/best.pt') == lock['best_sha256'], 'Solo lock changed')
    asset = t.load_pt(d / '_prepared/fold_training_assets.pt')
    frozen = t.load_pt(d / 'frozen_profiles.pt')
    t.require(asset['train_window_ids'] == train.window_id.tolist(), 'TRAIN asset identity')
    ids = np.asarray(asset['feature_ids'], int)
    scales = np.asarray(asset['train_p95_scales'])
    t.require(np.array_equal(ids, np.asarray(frozen['feature_ids'])) and np.array_equal(scales, np.asarray(frozen['scales'])), 'Candidate/P95 mismatch')
    config = t.read(s / 'config.json')
    t.require(config['fold_assets_sha256'] == t.digest(d / '_prepared/fold_training_assets.pt'), 'Solo normalization mismatch')
    select = t.read(d / 'automatic_k_selection.json')['selected']['performance']
    bundle = t.load_pt(d / f"K{select['K']}/jwarm_bundle.pt")
    t.require(bundle['K'] == select['K'] and bundle['cv_fold'] == fold and (bundle['seed'] == 42), 'Full lock mismatch')
    wm = ck['state']['profile.weights'].numpy().astype(float)
    wb = frozen['weights'].numpy().astype(float)
    wd = np.zeros_like(wb)
    lookup = {v: i for i, v in enumerate(ids)}
    positions = np.array([lookup[int(v)] for v in bundle['feature_ids']])
    t.require(np.array_equal(scales[positions], np.asarray(bundle['train_p95_scales'])), 'Stage D P95 mismatch')
    wd[positions] = np.asarray(bundle['W_final'])
    kept = np.zeros(len(ids), bool)
    kept[positions] = True
    zm, tm = pred(s / 'MotifAE_only_val.npz', meta)
    za, ta = pred(s / 'Explicit_only_val.npz', meta)
    zd, td = pred(s / 'Full_performance_val.npz', meta)
    t.require(td == bundle['threshold'], 'Full threshold changed')
    zb, _ = pred(d / 'predictions/full_replay.npz', meta)
    expb = np.load(d / 'explicit_val.npy')
    with np.load(d / f"K{select['K']}/val_predictions.npz") as a:
        t.require(np.array_equal(a['window_ids'].astype(str), meta.window_id), 'D replay row identity')
        expd = a['explicit_logits'].copy()
    y = meta.label.to_numpy(float)
    cohorts = dict(FP_corrected=(y == 0) & (expit(za) >= ta) & (expit(zd) < td), FN_rescued=(y == 1) & (expit(za) < ta) & (expit(zd) >= td))
    paths = [s / 'motifae_only/best.pt', s / 'motifae_locked.json', s / 'MotifAE_only_val.npz', s / 'Explicit_only_val.npz', s / 'Full_performance_val.npz', d / 'frozen_profiles.pt', d / '_prepared/fold_training_assets.pt', d / 'automatic_k_selection.json', d / f"K{select['K']}/jwarm_bundle.pt", d / 'explicit_val.npy', d / 'predictions/full_replay.npz', CV / f'fold_{fold:02d}/val_metadata.tsv', CV / f'fold_{fold:02d}/train_metadata.tsv']
    fingerprint = {str(p): t.digest(p) for p in paths}
    for name in ['train', 'val']:
        for f in ['columns.bin', 'values.bin', 'indptr.npy', 'complete.json']:
            p = d / '_prepared' / name / f
            fingerprint[str(p)] = t.digest(p)
    return dict(fold=fold, source=d, meta=meta, train=train, ids=ids, scales=scales, wm=wm, wb=wb, wd=wd, kept=kept, zm=zm, zb=zb, zd=zd, za=za, tm=tm, ta=ta, td=td, expb=expb, expd=expd, y=y, cohorts=cohorts, intercept=float(ck['state']['intercept']), K=int(select['K']), fingerprint=fingerprint)

def aggregate_store(directory, k, weights=None, train_meta=None):
    store = t.CompactStore(directory)
    n = len(store.ptr) - 1
    blocks = []
    contribution_blocks = [[] for _ in weights or []]
    position_sum = np.zeros(k * 29, float)
    for start in range(0, n, 1024):
        end = min(start + 1024, n)
        lo, hi = (store.ptr[start], store.ptr[end])
        cols = np.asarray(store.column[lo:hi])
        v = np.asarray(store.value[lo:hi])
        ptr = np.array(store.ptr[start:end + 1] - lo)
        t.require(np.isfinite(v).all() and (v >= 0).all(), 'Bad normalized values')
        position_sum += np.bincount(cols, weights=v, minlength=k * 29)
        x = sp.csr_matrix((v.copy(), cols // 29, ptr), shape=(end - start, k))
        x.sum_duplicates()
        x.eliminate_zeros()
        blocks.append(x)
        for acc, w in zip(contribution_blocks, weights or []):
            c = sp.csr_matrix((v * w.reshape(-1)[cols], cols // 29, ptr), shape=(end - start, k))
            c.sum_duplicates()
            c.eliminate_zeros()
            acc.append(c)
    x = sp.vstack(blocks, format='csr')
    del blocks
    cs = [sp.vstack(b, format='csr') for b in contribution_blocks]
    return (x, cs, position_sum.reshape(k, 29) / n)

def prepare(fold):
    f = load_fold(fold)
    d = OUT / f'fold_{fold:02d}'
    d.mkdir(parents=True, exist_ok=True)
    sig = t.identity(dict(inputs=f['fingerprint'], code=t.digest(__file__), protocol='3004c3934ef6797bd777ad53ba35e643340b07b2adb1d118156de3c4d8d67d5a'))
    if (d / 'prepared.json').exists():
        t.require(t.read(d / 'prepared.json')['signature'] == sig, 'Analysis inputs differ from the saved run; use a new output directory')
        return (f, d)
    logging.info('Fold %02d sparse TRAIN summaries, then held-out contribution replay', fold)
    tx, _, mean = aggregate_store(f['source'] / '_prepared/train', len(f['ids']))
    trainactive = np.asarray(tx.getnnz(axis=0)).ravel() / len(f['train'])
    codes, proteins = pd.factorize(f['train'].base_protein_id)
    group = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(proteins), len(codes)))
    protein_n = np.asarray((group @ tx).getnnz(axis=0)).ravel()
    broad = (trainactive >= np.quantile(trainactive, 0.75)) & (protein_n >= 10)
    sp.save_npz(d / 'train_activation_sum.npz', tx)
    del tx, group
    gc.collect()
    vx, cs, _ = aggregate_store(f['source'] / '_prepared/val', len(f['ids']), [f['wm'], f['wb'], f['wd']])
    cm, cb, cd = cs
    errors = dict(M=float(np.max(abs(np.asarray(cm.sum(1)).ravel() + f['intercept'] - f['zm']))), B=float(np.max(abs(np.asarray(cb.sum(1)).ravel() + f['expb'] - f['zb']))), D=float(np.max(abs(np.asarray(cd.sum(1)).ravel() + f['expd'] - f['zd']))))
    t.require(max(errors.values()) <= 1e-05, f'Contribution replay failure: {errors}')
    for name, m in zip(['activation', 'M', 'B', 'D'], [vx, cm, cb, cd]):
        sp.save_npz(d / f'{name}.npz', m)
    np.savez_compressed(d / 'train_stats.npz', ids=f['ids'], position_mean=mean, coverage=trainactive, protein_n=protein_n, broad=broad, mu_m=(mean * f['wm']).sum(1), mu_b=(mean * f['wb']).sum(1), mu_d=(mean * f['wd']).sum(1))
    t.write_json(d / 'prepared.json', dict(signature=sig, replay=errors, fold=fold, K=f['K'], inputs=f['fingerprint'], test_read=False, cohort_counts={k: int(v.sum()) for k, v in f['cohorts'].items()}))
    return (f, d)

def screen(f, d):
    if (d / 'screen_complete.json').exists():
        return
    cm, cb, cd = [sp.load_npz(d / f'{n}.npz').tocsc() for n in ['M', 'B', 'D']]
    act = sp.load_npz(d / 'activation.npz').tocsc()
    stats = np.load(d / 'train_stats.npz')
    y = f['y']
    w = balanced(f['meta'])
    weighted_cohort = {k: balanced(f['meta'], m) for k, m in f['cohorts'].items()}
    fs = []
    gs = []
    components = f['meta'].component_id.to_numpy()
    proteins = f['meta'].base_protein_id.to_numpy()
    p = expit(f['zd'])
    baseline_correct = (p >= f['td']) == (y == 1)
    for j, fid in enumerate(f['ids']):
        m = cm[:, j].toarray().ravel()
        b = cb[:, j].toarray().ravel()
        c = cd[:, j].toarray().ravel()
        mu_m = float(stats['mu_m'][j])
        mu_b = float(stats['mu_b'][j])
        dm = delta(y, f['zm'], m)
        db = delta(y, f['zb'], b)
        dmc = delta(y, f['zm'], m - mu_m)
        dbc = delta(y, f['zb'], b - mu_b)
        active = act[:, j].getnnz() / len(y)
        gs.append(dict(fold=f['fold'], feature_id=int(fid), available=True, broad=bool(stats['broad'][j]), train_active_fraction=float(stats['coverage'][j]), train_active_proteins=int(stats['protein_n'][j]), val_active_fraction=active, M_raw_bce=float(w @ dm), B_raw_bce=float(w @ db), M_center_bce=float(w @ dmc), B_center_bce=float(w @ dbc), center_dependency_drop=float(w @ (dmc - dbc)), M_abs_contribution=float(w @ abs(m)), B_abs_contribution=float(w @ abs(b)), M_center_abs=float(w @ abs(m - mu_m)), B_center_abs=float(w @ abs(b - mu_b)), mu_M=mu_m, mu_B=mu_b, retained_D=bool(f['kept'][j])))
        if not f['kept'][j]:
            continue
        delete = delta(y, f['zd'], c)
        pp = expit(f['zd'] - c)
        effect = abs(p - pp)
        new_error = baseline_correct & ((pp >= f['td']) != (y == 1))
        for kind, mask in f['cohorts'].items():
            if not mask.any():
                continue
            strong = mask & (delete > 0) & (effect >= 0.05)
            fs.append(dict(fold=f['fold'], feature_id=int(fid), kind=kind, cohort_n=int(mask.sum()), component_balanced_bce=float(weighted_cohort[kind] @ delete), mean_probability_effect=float(effect[mask].mean()), active_fraction=float((act[:, j].toarray().ravel()[mask] > 0).mean()), strong_windows=int(strong.sum()), strong_protein_n=len(set(proteins[strong])), strong_component_n=len(set(components[strong])), strong_proteins=';'.join(sorted(set(proteins[strong]))), undo_windows=int((mask & new_error).sum()), undo_protein_n=len(set(proteins[mask & new_error])), mean_abs_contribution=float(abs(c[mask]).mean())))
    save(fs, d / 'f_screen.tsv')
    save(gs, d / 'g_screen.tsv')
    cd = cd.tocsr()
    records = []
    for kind, mask in f['cohorts'].items():
        for n in np.flatnonzero(mask):
            row = cd.getrow(n)
            helpful = (2 * y[n] - 1) * row.data
            order = np.argsort(-helpful)
            order = order[helpful[order] > 0]
            value = row.data[order]
            total = abs(value).sum()
            cumulative = np.cumsum(value)
            z0 = f['zd'][n] - value.sum()
            satisfied = expit(z0 + cumulative) >= f['td'] if y[n] else expit(z0 + cumulative) < f['td']
            min_k = 0 if (expit(z0) >= f['td']) == bool(y[n]) else int(np.flatnonzero(satisfied)[0] + 1) if satisfied.any() else -1
            rec = dict(fold=f['fold'], kind=kind, window_id=f['meta'].window_id.iloc[n], base_protein_id=proteins[n], component_id=components[n], favorable_n=len(value), favorable_abs_sum=total, minimum_favorable_k=min_k, full_logit=f['zd'][n], probability=p[n], threshold=f['td'], abs_logit_margin=abs(f['zd'][n] - logit(f['td'])))
            for k in [1, 3, 5]:
                znew = f['zd'][n] - value[:k].sum()
                rec[f'top{k}_fraction'] = float(abs(value[:k]).sum() / total) if total else 0.0
                rec[f'top{k}_deletion_undo'] = bool((expit(znew) >= f['td']) != bool(y[n]))
            records.append(rec)
    save(records, d / 'f_concentration.tsv')
    t.write_json(d / 'screen_complete.json', dict(fold=f['fold']))
    logging.info('Fold %02d screen done; corrected/rescued=%s', f['fold'], {k: int(v.sum()) for k, v in f['cohorts'].items()})

def lock():
    fs = pd.concat([table(OUT / f'fold_{i:02d}/f_screen.tsv') for i in range(1, 6)])
    gs = pd.concat([table(OUT / f'fold_{i:02d}/g_screen.tsv') for i in range(1, 6)])
    frows = []
    for (kind, fid), g in fs.groupby(['kind', 'feature_id']):
        proteins = set(';'.join(g.strong_proteins).split(';')) - {''}
        frows.append(dict(kind=kind, feature_id=int(fid), stable_folds=int((g.strong_protein_n >= 3).sum()), proteins=len(proteins), mean_bce=float(g.component_balanced_bce.sum() / 5), mean_probability_effect=float(g.mean_probability_effect.sum() / 5)))
    fd = pd.DataFrame(frows)
    fd['eligible'] = (fd.stable_folds >= 3) & (fd.proteins >= 10) & (fd.mean_bce > 0)
    chosen_f = fd[fd.eligible].sort_values(['kind', 'mean_bce', 'proteins', 'feature_id'], ascending=[True, False, False, True]).groupby('kind').head(3)
    grows = []
    for fid, g in gs.groupby('feature_id'):
        grows.append(dict(feature_id=int(fid), available_folds=len(g), broad_folds=int(g.broad.sum()), positive_decline_folds=int(((g.M_center_bce > 0) & (g.center_dependency_drop > 0)).sum()), mean_drop=float(g.center_dependency_drop.mean()), mean_M=float(g.M_center_bce.mean()), train_coverage=float(g.train_active_fraction.mean())))
    gd = pd.DataFrame(grows)
    gd['eligible'] = (gd.available_folds >= 3) & (gd.broad_folds >= 3) & (gd.positive_decline_folds >= 3) & (gd.mean_drop > 0) & (gd.mean_M > 0)
    chosen_g = gd[gd.eligible].sort_values(['mean_drop', 'feature_id'], ascending=[False, True]).head(5)
    save(fd, OUT / 'discovery_f_candidates.tsv')
    save(gd, OUT / 'discovery_g_candidates.tsv')
    record = dict(f=chosen_f.to_dict('records'), g=chosen_g.to_dict('records'), discovery_folds=list(range(1, 6)), followup_folds=list(range(6, 11)), protocol_sha256='3004c3934ef6797bd777ad53ba35e643340b07b2adb1d118156de3c4d8d67d5a')
    path = OUT / 'locked_candidates.json'
    if path.exists():
        t.require(t.read(path) == t.jsonable(record), 'Candidate lock changed')
    else:
        t.write_json(path, record)
    logging.info('Locked f=%s g=%s', [(r['kind'], r['feature_id']) for r in record['f']], [r['feature_id'] for r in record['g']])
    return record

def train_group(tx, j, ids):
    n = tx.shape[0]
    mean = np.asarray(tx.sum(0)).ravel() / n
    sd = np.sqrt(np.maximum(np.asarray(tx.power(2).sum(0)).ravel() / n - mean ** 2, 0))
    x = tx[:, j]
    cross = np.asarray((x.T @ tx).toarray()).ravel() / n - mean[j] * mean
    corr = np.divide(cross, sd[j] * sd, out=np.zeros_like(cross), where=sd[j] * sd > 1e-12)
    candidates = np.flatnonzero(corr >= 0.7)
    neighbors = sorted((i for i in candidates if i != j), key=lambda i: (-corr[i], ids[i]))[:19]
    return (np.array([j, *neighbors]), corr)

def detail(f, d, locked):
    signature = t.identity(locked)
    marker = d / 'detail_complete.json'
    if marker.exists():
        t.require(t.read(marker)['lock_signature'] == signature, 'Detail lock differs')
        return
    cm, cb, cd = [sp.load_npz(d / f'{n}.npz').tocsc() for n in ['M', 'B', 'D']]
    stats = np.load(d / 'train_stats.npz')
    lookup = {int(fid): j for j, fid in enumerate(f['ids'])}
    y = f['y']
    w = balanced(f['meta'])
    frows = []
    grows = []
    groups = []
    for entry in locked['f']:
        fid = entry['feature_id']
        kind = entry['kind']
        if fid not in lookup:
            continue
        j = lookup[fid]
        c = cd[:, j].toarray().ravel()
        change = delta(y, f['zd'], c)
        pp = expit(f['zd'] - c)
        p = expit(f['zd'])
        for n in np.flatnonzero(f['cohorts'][kind]):
            r = f['meta'].iloc[n]
            frows.append(dict(fold=f['fold'], kind=kind, feature_id=fid, window_id=r.window_id, base_protein_id=r.base_protein_id, component_id=r.component_id, negative_type=r.negative_type, protein_id=r.protein_id, window_sequence=r.window_sequence, label=int(y[n]), contribution=c[n], bce_deletion=change[n], probability=p[n], without_probability=pp[n], threshold=f['td'], probability_effect=abs(p[n] - pp[n]), undo=bool((pp[n] >= f['td']) != bool(y[n])), retained=bool(f['kept'][j]), logit_margin=f['zd'][n] - logit(f['td']), mean_preserving_bce=float(delta(y[n], f['zd'][n], c[n] - stats['mu_d'][j]))))
    tx = sp.load_npz(d / 'train_activation_sum.npz').astype(float)
    for entry in locked['g']:
        fid = entry['feature_id']
        if fid not in lookup:
            continue
        j = lookup[fid]
        group, corr = train_group(tx, j, f['ids'])
        groups.append(dict(fold=f['fold'], feature_id=fid, group_ids=';'.join(map(str, f['ids'][group])), correlations=';'.join((f'{corr[i]:.6f}' for i in group)), group_size=len(group)))
        for scope, idx in [('single', np.array([j])), ('correlated_group', group)]:
            for name, matrix, z, mu in [('M', cm, f['zm'], stats['mu_m']), ('B', cb, f['zb'], stats['mu_b'])]:
                c = np.asarray(matrix[:, idx].sum(1)).ravel()
                mean = float(mu[idx].sum())
                for centered in [False, True]:
                    removed = c - mean if centered else c
                    dz = delta(y, z, removed)
                    new = z - removed
                    grows.append(dict(fold=f['fold'], feature_id=fid, scope=scope, model=name, centered=centered, mean_bce=float(w @ dz), window_mean_bce=float(dz.mean()), ap_drop=float(average_precision_score(y, expit(z)) - average_precision_score(y, expit(new))), auc_drop=float(roc_auc_score(y, z) - roc_auc_score(y, new)), mean_probability_change=float(w @ abs(expit(z) - expit(new))), group_size=len(idx)))
    save(frows, d / 'f_locked_windows.tsv')
    save(grows, d / 'g_locked_deletions.tsv')
    save(groups, d / 'g_train_correlated_groups.tsv')
    t.write_json(marker, dict(lock_signature=signature))

def cases(locked):
    allrows = pd.concat([table(OUT / f'fold_{i:02d}/f_locked_windows.tsv') for i in range(6, 11)])
    chosen = []
    profiles = []
    selected_rows = []
    for entry in locked['f']:
        fid = entry['feature_id']
        kind = entry['kind']
        rows = allrows[(allrows.feature_id == fid) & (allrows.kind == kind) & (allrows.bce_deletion > 0) & (allrows.probability_effect >= 0.05) & allrows.undo.astype(bool)].copy()
        if not len(rows):
            continue
        rows['distance'] = abs(rows.probability_effect - rows.probability_effect.median())
        rows = rows.sort_values(['distance', 'window_id']).drop_duplicates('component_id').head(8)
        selected_rows.extend(rows.to_dict('records'))
    if not selected_rows:
        save([], OUT / 'f_cases.tsv')
        save([], OUT / 'f_position_profiles.tsv')
        return
    for fold, rows in pd.DataFrame(selected_rows).groupby('fold'):
        f = load_fold(int(fold))
        for _, r in rows.iterrows():
            fid = int(r.feature_id)
            kind = r.kind
            j = int(np.flatnonzero(f['ids'] == fid)[0])
            n = int(np.flatnonzero(f['meta'].window_id == r.window_id)[0])
            store = t.CompactStore(f['source'] / '_prepared/val')
            lo, hi = store.ptr[n:n + 2]
            col = store.column[lo:hi]
            v = store.value[lo:hi]
            keep = col // 29 == j
            track = np.bincount(col[keep] % 29, weights=v[keep], minlength=29)
            con = track * f['wd'][j]
            t.require(abs(con.sum() - r.contribution) < 1e-06, 'Case contribution mismatch')
            chosen.append(r.to_dict())
            for i in range(29):
                profiles.append(dict(feature_id=fid, kind=kind, window_id=r.window_id, fold=int(r.fold), position=i - 14, residue=r.window_sequence[i], normalized_activation=track[i], contribution=con[i]))
    save(chosen, OUT / 'f_cases.tsv')
    save(profiles, OUT / 'f_position_profiles.tsv')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--phase', choices=['all', 'discovery', 'followup', 'report'], default='all')
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', stream=sys.stdout)
    t.torch.set_num_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.phase in ['all', 'discovery']:
        for i in range(1, 6):
            f, d = prepare(i)
            screen(f, d)
            del f
            gc.collect()
        locked = lock()
    else:
        locked = t.read(OUT / 'locked_candidates.json')
    if a.phase in ['all', 'followup']:
        for i in range(1, 11):
            f, d = prepare(i)
            screen(f, d)
            detail(f, d, locked)
            del f
            gc.collect()
        cases(locked)
    if a.phase in ['all', 'report']:
        pass
if __name__ == '__main__':
    main()
