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
from scipy.stats import t as student_t, rankdata
from sklearn.metrics import average_precision_score, roc_auc_score, precision_recall_curve
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sequence_descriptors import construct, CATEGORIES
HERE = Path(__file__).resolve().parent
ROOT = Path(__import__('os').environ['PRIORSLIM_DATA_ROOT'])
CV = ROOT / 'new_data/cv10_homology_seed2026'
MODELS = ROOT / 'train/outputs_jwarm_cv10/seed42'
STORE = ROOT / 'feature_cache/motifae_sparse_windows_idr_coverage'
TRANSITIONS = ['TP_kept', 'FN_rescued', 'TP_lost', 'FN_still', 'TN_kept', 'FP_corrected', 'TN_harmed', 'FP_still']

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def clean(v):
    if isinstance(v, dict):
        return {str(k): clean(x) for (k, x) in v.items()}
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, np.ndarray):
        return clean(v.tolist())
    if isinstance(v, np.generic):
        return clean(v.item())
    if isinstance(v, float) and (not np.isfinite(v)):
        return None
    return v

def js(path, value):
    path.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def read(path):
    return json.loads(path.read_text())

def table(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, low_memory=False, float_precision='round_trip')

def write(frame, path):
    if not len(frame.columns):
        frame = pd.DataFrame(columns=['status'])
    frame.to_csv(path, sep='\t', index=False)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda : f.read(8 * 1024 ** 2), b''):
            h.update(b)
    return h.hexdigest()

def adjust_bh(p):
    p = np.asarray(p, float)
    result = np.ones(len(p))
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    if len(order):
        result[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1])
    return result

def load_bundle(path):
    return torch.load(path, map_location='cpu', weights_only=True)

def transition(y, e, f, threshold):
    ep = expit(e) >= threshold
    fp = expit(f) >= threshold
    y = np.asarray(y) == 1
    return np.select([y & ep & fp, y & ~ep & fp, y & ep & ~fp, y & ~ep & ~fp, ~y & ~ep & ~fp, ~y & ep & ~fp, ~y & ~ep & fp, ~y & ep & fp], TRANSITIONS, default='ERROR')

def metadata():
    combined = table(ROOT / 'train/combined_metadata.tsv')
    frames = []
    for fold in range(1, 11):
        m = table(CV / f'fold_{fold:02d}/val_metadata.tsv')
        m['fold'] = fold
        m['cohort'] = 'development'
        frames.append(m)
    m = pd.concat(frames, ignore_index=True)
    require(len(m) == 87801 and (not m.window_id.duplicated().any()), 'Development OOF IDs')
    test = table(ROOT / 'new_data/test_metadata.tsv')
    base = combined[combined.split == 'test'].sort_values('split_row')
    require(np.array_equal(test.window_sequence, base.sequence), 'Test sequence alignment')
    require(np.array_equal(test.protein_id, base.uniprot), 'Test protein alignment')
    for key in ('window_id', 'tensor_index'):
        test[key] = base[key].values
    test['base_protein_id'] = test.protein_id.str.replace('_(?:C|M|F)\\d+$', '', regex=True)
    test['component_id'] = 'test:' + test.base_protein_id
    test['fold'] = 0
    test['cohort'] = 'test'
    m = pd.concat([m, test], ignore_index=True)
    cache = combined.set_index('window_id').loc[m.window_id]
    for (x, y) in [('sequence', 'window_sequence'), ('uniprot', 'protein_id')]:
        require(np.array_equal(cache[x].values, m[y].values), 'Metadata ' + x)
    for col in ('window_start', 'window_end', 'label', 'tensor_index'):
        require(np.array_equal(cache[col].astype(float), m[col].astype(float)), 'Metadata ' + col)
    for col in ('slim_start', 'slim_end'):
        m[col] = cache[col].values
    require(not m.window_id.duplicated().any(), 'Metadata duplicate window_id')
    for col in ('label', 'window_start', 'window_end', 'tensor_index', 'fold'):
        m[col] = m[col].astype(float).astype(int)
    manifest = table(ROOT / 'data/new_all_positive_augmentation_manifest.tsv')

    def keys(frame, protein):
        return list(zip(frame[protein], frame.window_start.astype(float).astype(int), frame.window_end.astype(float).astype(int), frame.slim_start.astype(float).astype(int), frame.slim_end.astype(float).astype(int)))
    amap = {}
    for (key, kind) in zip(keys(manifest, 'uniprot'), manifest.augmentation_type):
        require(key not in amap or amap[key] == kind, 'Conflicting positive augmentation provenance')
        amap[key] = kind
    m['augmentation_type'] = 'negative'
    pos = m.label == 1
    poskeys = keys(m[pos], 'protein_id')
    require(all((k in amap for k in poskeys)), 'Missing positive augmentation provenance')
    m.loc[pos, 'augmentation_type'] = [amap[k] for k in poskeys]
    m['sample_type'] = np.where(pos, 'positive:' + m.augmentation_type, m.negative_type)
    m['family'] = m.motif_family.replace({'NA': '', 'nan': ''})
    m.loc[pos, 'family'] = m.loc[pos, 'nearest_slim_id'].map(lambda x: str(x).split(':', 2)[-1])
    m['family'] = m.family.replace('', 'unknown')
    m['canonical29'] = m.window_sequence.str.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{29}').fillna(False)
    dev = m[m.cohort == 'development']
    seqcross = set(dev.groupby('window_sequence').fold.nunique().loc[lambda x: x > 1].index)
    bases = dev[['base_protein_id', 'fold']].drop_duplicates().copy()
    bases['isoform_family'] = bases.base_protein_id.str.split('-').str[0]
    isocross = set(bases.groupby('isoform_family').fold.nunique().loc[lambda x: x > 1].index)
    m['cross_fold_exact_sequence'] = m.window_sequence.isin(seqcross)
    m['cross_fold_isoform_family'] = m.base_protein_id.str.split('-').str[0].isin(isocross)
    m['analysis_partition'] = np.where(m.cohort == 'test', 'test_descriptive', np.where(m.fold <= 5, 'discovery', 'replication'))
    m['discovery_eligible'] = (m.cohort == 'development') & m.canonical29 & ~m.cross_fold_exact_sequence & ~m.cross_fold_isoform_family
    original = (m.augmentation_type == 'original') & m.discovery_eligible
    dup = m[original].duplicated('parent_slim_id', keep='first')
    m.loc[m[original].index[dup], 'discovery_eligible'] = False
    require(m[m.cohort == 'development'].groupby('component_id').fold.nunique().max() == 1, 'Component crossed folds')
    return (m, combined)

def prediction_rows(meta, fold, role, cohort):
    d = MODELS / f'fold_{fold:02d}'
    selection = read(d / 'automatic_k_selection.json')
    require(not selection['policy']['test_used'], 'K used Test')
    k = int(selection['selected'][role]['K'])
    folder = d / f'K{k}'
    b = load_bundle(folder / 'jwarm_bundle.pt')
    require(b['K'] == k and b['cv_fold'] == fold and (b['seed'] == 42), 'Bundle identity')
    if cohort == 'development':
        path = folder / 'val_predictions.npz'
        frame = meta[(meta.cohort == cohort) & (meta.fold == fold)]
    else:
        path = MODELS / 'test_evaluation' / f'fold_{fold:02d}' / f'K{k}' / 'test_predictions.npz'
        frame = meta[meta.cohort == cohort]
    with np.load(path) as n:
        ids = n['window_ids'].astype(str)
        require(len(set(ids)) == len(ids), 'Duplicate prediction IDs')
        require(set(frame.window_id) == set(ids), 'Prediction coverage')
        f = frame.set_index('window_id').loc[ids].reset_index()
        require(np.array_equal(f.label.values, n['labels']), 'Prediction labels')
        for (col, key) in [('explicit_logit', 'explicit_logits'), ('final_logit', 'logits'), ('motifae_logit', 'motifae_logits')]:
            f[col] = n[key]
        require(np.allclose(f.explicit_logit + f.motifae_logit, f.final_logit, atol=1e-06), 'Logit decomposition')
    f['model_fold'] = fold
    f['operating_point'] = role
    f['K'] = k
    f['threshold'] = float(b['threshold'])
    f['explicit_probability'] = expit(f.explicit_logit)
    f['final_probability'] = expit(f.final_logit)
    f['explicit_margin'] = f.explicit_logit - logit(float(b['threshold']))
    f['transition'] = transition(f.label, f.explicit_logit, f.final_logit, float(b['threshold']))
    require(np.isfinite(f[['explicit_logit', 'final_logit', 'motifae_logit']]).all().all(), 'Nonfinite score')
    return (f, b, path)

def summarize_transitions(frame):
    rows = []
    for ((cohort, role, fold), d) in frame.groupby(['cohort', 'operating_point', 'model_fold']):
        for (name, s) in [('all', d)] + list(d.groupby('sample_type')):
            counts = s.transition.value_counts().to_dict()
            n = len(s)
            p = int(s.label.sum())
            neg = n - p
            row = dict(cohort=cohort, operating_point=role, fold=fold, sample_type=name, n=n, positive_n=p, negative_n=neg)
            for tr in TRANSITIONS:
                row[tr] = int(counts.get(tr, 0))
            row['positive_rescue_rate_among_explicit_FN'] = row['FN_rescued'] / max(1, row['FN_rescued'] + row['FN_still'])
            row['negative_correction_rate_among_explicit_FP'] = row['FP_corrected'] / max(1, row['FP_corrected'] + row['FP_still'])
            row['positive_harm_rate_among_explicit_TP'] = row['TP_lost'] / max(1, row['TP_lost'] + row['TP_kept'])
            row['negative_harm_rate_among_explicit_TN'] = row['TN_harmed'] / max(1, row['TN_harmed'] + row['TN_kept'])
            row['mean_signed_label_gain'] = float(((2 * s.label - 1) * s.motifae_logit).mean())
            if p and neg:
                row['explicit_ap'] = average_precision_score(s.label, s.explicit_probability)
                row['final_ap'] = average_precision_score(s.label, s.final_probability)
                row['explicit_auc'] = roc_auc_score(s.label, s.explicit_probability)
                row['final_auc'] = roc_auc_score(s.label, s.final_probability)
            rows.append(row)
    return pd.DataFrame(rows)

def sparse_contributions(frame, out):
    cfg = read(STORE / 'motifae_sparse_config.json')
    require(cfg['complete'] and (not cfg['dense_29_by_40960_constructed']), 'Sparse provenance')
    sm = table(STORE / 'motifae_window_metadata.tsv').set_index('window_id')
    sm = sm.loc[frame.window_id]
    require(np.array_equal(sm.base_protein_id.values, frame.base_protein_id.values), 'Sparse parent identity')
    require(np.array_equal(sm.window_start.astype(int), frame.window_start), 'Sparse coordinates')
    require((sm.crop_mode == 'idr_coverage').all() and (sm.context_margin.astype(int) == 128).all(), 'Wrong ESM crop protocol')
    ptr = np.load(STORE / 'motifae_residue_indptr.int64.npy', mmap_mode='r')
    ids = np.memmap(STORE / 'motifae_feature_indices.int32.bin', mode='r', dtype=np.int32)
    values = np.memmap(STORE / 'motifae_feature_values.float32.bin', mode='r', dtype=np.float32)
    (rows, cols, vals) = ([], [], [])
    position = np.zeros((len(frame), 29), np.float32)
    audits = []
    top = []
    selected = np.zeros((11, 40960), bool)
    for (fold, g) in frame.groupby('model_fold', sort=True):
        k = int(g.K.iloc[0])
        b = load_bundle(MODELS / f'fold_{fold:02d}' / f'K{k}' / 'jwarm_bundle.pt')
        fids = b['feature_ids'].numpy()
        scales = b['train_p95_scales'].numpy()
        W = b['W_final'].numpy()
        selected[fold, fids] = True
        lookup = np.full(40960, -1, int)
        lookup[fids] = np.arange(len(fids))
        errs = []
        for (idx, row) in g.iterrows():
            sn = int(sm.iloc[idx].sparse_window_index)
            p = ptr[sn * 29:sn * 29 + 30]
            (a, z) = (int(p[0]), int(p[-1]))
            jj = lookup[ids[a:z]]
            keep = jj >= 0
            pos = np.repeat(np.arange(29), np.diff(p))[keep]
            j = jj[keep]
            v = np.asarray(values[a:z])[keep] / (scales[j] + 1e-08)
            cc = v.astype(float) * W[j, pos].astype(float)
            require(np.isfinite(cc).all(), 'Contribution nonfinite')
            summed = np.bincount(j, weights=cc, minlength=len(fids))
            nz = np.flatnonzero(summed)
            rows.append(np.full(len(nz), idx, np.int32))
            cols.append(fids[nz])
            vals.append(summed[nz].astype(np.float32))
            position[idx] = np.bincount(pos, weights=cc, minlength=29)
            errs.append(abs(float(cc.sum()) - float(row.motifae_logit)))
            largest = np.argsort(np.abs(summed))[::-1][:5]
            for (rank, jj) in enumerate(largest, 1):
                if summed[jj]:
                    top.append(dict(window_id=row.window_id, fold=fold, rank=rank, feature_id=int(fids[jj]), contribution=float(summed[jj]), transition=row.transition, augmentation_type=row.augmentation_type))
        maxerr = max(errs)
        require(maxerr < 1e-05, f'Sparse replay error fold={fold}: {maxerr}')
        audits.append(dict(fold=fold, K=k, n=len(g), max_motifae_logit_error=maxerr))
        logging.info('Contributions fold=%02d rows=%d replay_error=%.3g', fold, len(g), maxerr)
    mat = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(frame), 40960)).tocsr()
    sp.save_npz(out / 'development_contributions_by_feature.npz', mat)
    np.savez_compressed(out / 'development_contributions_by_position.npz', window_ids=np.asarray(frame.window_id.tolist(), dtype=np.str_), contributions=position, selected_features_by_fold=selected)
    write(pd.DataFrame(top), out / 'top5_contributions_per_development_window.tsv')
    js(out / 'contribution_replay_checks.json', audits)
    return (mat, position, selected)

def biological_features(frame):
    rows = frame.tensor_index.values
    manifest = read(ROOT / 'train/combined_feature_manifest.json')
    with np.load(ROOT / 'train/combined_tensor.npz') as a:
        require(np.array_equal(a['window_ids'][rows].astype(str), frame.window_id), 'Tensor ID alignment')
        chem = a['chemical'][rows]
        access = a['accessibility'][rows]
        p = a['pssm_residue'][rows]
        mask = a['pssm_valid_mask'][rows].astype(float)
        peak = a['pssm_window'][rows].reshape(-1)
    existing = pd.DataFrame(chem, columns=manifest['chemical_features'])
    for (i, name) in enumerate(manifest['accessibility_features']):
        existing['access_' + name] = access[:, i]
    for (i, name) in enumerate(manifest['pssm_residue_features']):
        existing['mean_' + name] = (p[:, :, i] * mask).sum(1) / np.maximum(1, mask.sum(1))
    existing['local_conservation_peak'] = peak
    existing['valid_pssm_fraction'] = mask.mean(1)
    existing['log_idr_length'] = np.log1p(frame.idr_end.astype(float).values - frame.idr_start.astype(float).values + 1)
    (new, definitions) = construct(frame.window_sequence.values, p[:, :, 5])
    return (existing, new, definitions)

def design_matrix(existing, frame, mask, common_families):
    x = existing.loc[mask].to_numpy(float)
    columns = [np.ones((len(x), 1)), x, frame.loc[mask, ['explicit_margin']].to_numpy(float)]
    for fold in range(1, 11):
        columns.append((frame.loc[mask, 'fold'].values == fold)[:, None].astype(float))
    for fam in common_families:
        columns.append((frame.loc[mask, 'family'].values == fam)[:, None].astype(float))
    x = np.hstack(columns)
    sd = x.std(0)
    x[:, 1:] = (x[:, 1:] - x[:, 1:].mean(0)) / np.maximum(sd[1:], 1e-10)
    (u, s, _) = np.linalg.svd(x, full_matrices=False)
    rank = int((s > s[0] * 1e-09).sum())
    return (u[:, :rank], rank)

def association_tests(frame, existing, new, mat, selected, signatures, partition, pairs=None, cohort='original'):
    natural_cohort = frame.augmentation_type == 'original' if cohort == 'original' else frame.sample_type == cohort
    base = frame.discovery_eligible & natural_cohort & (frame.analysis_partition == partition)
    counts = frame.loc[base, 'family'].value_counts()
    discovery = frame.discovery_eligible & natural_cohort & (frame.analysis_partition == 'discovery')
    common = frame.loc[discovery, 'family'].value_counts().loc[lambda x: x >= 10].index.tolist()
    results = []
    for fid in signatures.feature_id.astype(int):
        mask = base & selected[frame.fold.values, fid]
        idx = np.flatnonzero(mask)
        n = len(idx)
        if n < 100:
            continue
        cols = list(new.columns) if pairs is None else pairs.loc[pairs.feature_id == fid, 'descriptor'].tolist()
        if not cols:
            continue
        (Q, rank) = design_matrix(existing, frame, mask, common)
        y = mat[idx, fid].toarray().ravel()
        yr = y - Q @ (Q.T @ y)
        X = new.loc[mask, cols].to_numpy(float)
        xr = X - Q @ (Q.T @ X)
        groups = pd.factorize(frame.loc[mask, 'component_id'])[0]
        G = len(np.unique(groups))
        for (k, col) in enumerate(cols):
            r = xr[:, k]
            xx = r @ r
            total = ((X[:, k] - X[:, k].mean()) ** 2).sum()
            if xx < 1e-10 or total < 1e-10 or n <= rank + 2:
                continue
            beta = float(r @ yr / xx)
            err = yr - beta * r
            score = np.bincount(groups, weights=r * err, minlength=G)
            var = score @ score / (xx * xx) * G / max(1, G - 1) * (n - 1) / (n - rank - 1)
            se = float(np.sqrt(max(var, 0)))
            p = float(2 * student_t.sf(abs(beta / max(se, 1e-15)), df=G - 1))
            rho = float(r @ yr / np.sqrt(xx * max(yr @ yr, 1e-20)))
            results.append(dict(partition=partition, feature_id=fid, descriptor=col, n=n, components=G, nuisance_rank=rank, beta=beta, se=se, p=p, partial_r=rho, descriptor_residual_variance_fraction=xx / total, beta_low=beta - student_t.ppf(0.975, G - 1) * se, beta_high=beta + student_t.ppf(0.975, G - 1) * se))
    result = pd.DataFrame(results)
    if len(result):
        result['q'] = adjust_bh(result.p)
    return result
def negative_signatures(dev, mat, selected):
    eligible = dev.discovery_eligible & (dev.sample_type == 'motif_matched') & (dev.analysis_partition == 'discovery') & (dev.transition == 'FP_corrected')
    negative = (-mat[eligible.values]).maximum(0)
    sums = np.asarray(negative.sum(0)).ravel()
    ix = np.flatnonzero(eligible)
    rows = []
    for fid in np.flatnonzero(sums):
        x = negative[:, fid].toarray().ravel()
        g = dev.iloc[ix[x >= 0.05]]
        nf = int(selected[1:6, fid].sum())
        if g.base_protein_id.nunique() < 8 or nf < 3:
            continue
        rows.append(dict(feature_id=int(fid), discovery_negative_correction_magnitude=sums[fid], corrected_window_n=len(g), corrected_protein_n=g.base_protein_id.nunique(), selected_discovery_folds=nf, selected_replication_folds=int(selected[6:11, fid].sum())))
    return pd.DataFrame(rows).sort_values('discovery_negative_correction_magnitude', ascending=False).head(10)
