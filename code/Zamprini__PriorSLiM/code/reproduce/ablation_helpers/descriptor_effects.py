from pathlib import Path
import numpy as np
import pandas as pd
import _training as t
ROOT=t.PROJECT_DIR

EFIELDS = [('local_conservation_peak', 'Local conservation peak'), ('mean_consensus_percentile_IDR', 'Conservation percentile in IDR'), ('mean_pssm_information', 'PSSM information content'), ('motif_min_abs_offset', 'Motif distance from centre'), ('motif_center_enrichment', 'Motif centre enrichment'), ('motif_count_log', 'Motif match count (log)'), ('density_positive_percentile_IDR', 'Positive-charge density percentile'), ('center_enrich_polar_percentile_IDR', 'Polar centre enrichment percentile'), ('density_aromatic_percentile_IDR', 'Aromatic density percentile'), ('access_occlusion_cis', 'Within-IDR occlusion heuristic')]

def table(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, float_precision='round_trip')

def load_prediction(path, meta):
    with np.load(path) as d:
        t.require(np.array_equal(d['window_ids'].astype(str), meta.window_id.astype(str)), 'Prediction ID order differs')
        t.require(np.array_equal(d['labels'], meta.label.to_numpy(float)), 'Prediction labels differ')
        z = d['logits'].copy()
        threshold = float(d['threshold'])
    t.require(np.isfinite(z).all(), 'Nonfinite prediction')
    return (z, threshold)

def phenotype(meta):
    manifest = t.read(ROOT / 'train/combined_feature_manifest.json')
    combined = table(ROOT / 'train/combined_metadata.tsv').set_index('window_id').loc[meta.window_id]
    rows = combined.tensor_index.to_numpy(int)
    for (a, b) in [('sequence', 'window_sequence'), ('uniprot', 'protein_id')]:
        t.require(np.array_equal(combined[a].astype(str), meta[b].astype(str)), 'Test phenotype identity differs')
    t.require(np.array_equal(combined.label.to_numpy(float), meta.label.to_numpy(float)), 'Test phenotype labels differ')
    with np.load(ROOT / 'train/combined_tensor.npz') as f:
        t.require(np.array_equal(f['window_ids'][rows].astype(str), meta.window_id), 'Test tensor IDs differ')
        x = pd.DataFrame(f['chemical'][rows], columns=manifest['chemical_features'])
        for (i, col) in enumerate(manifest['accessibility_features']):
            x['access_' + col] = f['accessibility'][rows, i]
        mask = f['pssm_valid_mask'][rows]
        p = f['pssm_residue'][rows]
        for (i, col) in enumerate(manifest['pssm_residue_features']):
            x['mean_' + col] = (p[:, :, i] * mask).sum(1) / np.maximum(mask.sum(1), 1)
        x['local_conservation_peak'] = f['pssm_window'][rows].reshape(-1)
    t.require(np.isfinite(x.to_numpy()).all(), 'Nonfinite descriptors')
    manifest = table(ROOT / 'data/new_all_positive_augmentation_manifest.tsv')

    def keys(g, protein):
        return list(zip(g[protein], *[g[k].astype(float).astype(int) for k in ['window_start', 'window_end', 'slim_start', 'slim_end']]))
    mapping = {}
    for (key, kind) in zip(keys(manifest, 'uniprot'), manifest.augmentation_type):
        t.require(key not in mapping or mapping[key] == kind, 'Ambiguous augmentation metadata')
        mapping[key] = kind
    meta = meta.copy()
    meta['augmentation_type'] = 'negative'
    pos = meta.label.to_numpy(int) == 1
    meta.loc[pos, 'augmentation_type'] = [mapping[k] for k in keys(combined.iloc[np.flatnonzero(pos)], 'uniprot')]
    return (x, meta)
