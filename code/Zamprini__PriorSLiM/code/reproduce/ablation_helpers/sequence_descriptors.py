import numpy as np
import pandas as pd
CATEGORIES = {'positive': 'RK', 'negative': 'DE', 'aromatic': 'FYW', 'polar': 'STNQCH', 'hydrophobic': 'ILMV', 'proline': 'P', 'glycine': 'G', 'alanine': 'A'}
AA = 'ACDEFGHIKLMNPQRSTVWY'

def construct(sequences, conservation):
    a = np.array([list(s) for s in sequences])
    if a.ndim != 2 or a.shape[1] != 29:
        raise ValueError('Expected exactly 29 characters per window')
    pos = np.arange(-14, 15) / 14.0
    (values, definitions) = ({}, [])

    def add(name, x, formula, kind='new_relationship'):
        values[name] = np.asarray(x, dtype=np.float32)
        definitions.append(dict(feature=name, formula=formula, status=kind))
    tracks = {k: np.isin(a, list(v)).astype(np.float32) for (k, v) in CATEGORIES.items()}
    for (k, b) in tracks.items():
        alphabet = CATEGORIES[k]
        add(f'{k}_CminusN', b[:, 22:].mean(1) - b[:, :7].mean(1), f'mean(1[aa in {alphabet}], positions 23:29) - mean(..., positions 1:7)')
        add(f'{k}_signed_moment', (b * pos).mean(1), f'sum_i 1[aa_i in {alphabet}]*(i-15)/14 / 29; i=1..29')
        patches = (b[:, :-2] + b[:, 1:-1] + b[:, 2:]) / 3.0
        add(f'{k}_center_patch_contrast', patches[:, 7:20].max(1) - np.maximum(patches[:, :5].max(1), patches[:, 22:].max(1)), f'max 3-residue mean of {alphabet} wholly in positions 8:22 minus max wholly in 1:7 or 23:29')
        n = b.sum(1)
        expected = n * (n - 1) / (29 * 28)
        for lag in (1, 3):
            add(f'{k}_lag{lag}_excess', (b[:, :-lag] * b[:, lag:]).mean(1) - expected, f'mean_i 1[aa_i,aa_(i+{lag}) in {alphabet}] - n*(n-1)/(29*28)')
        centered = conservation - conservation.mean(1, keepdims=True)
        add(f'{k}_conservation_covariance', (b * centered).mean(1), f'mean_i (1[aa_i in {alphabet}]-mean indicator)*(consensus_percentile_IDR_i-mean percentile)', 'new_interaction_of_existing_inputs')
    pairs = [('positive', 'negative'), ('proline', 'hydrophobic'), ('proline', 'aromatic'), ('proline', 'negative'), ('proline', 'positive'), ('glycine', 'positive'), ('aromatic', 'negative'), ('hydrophobic', 'negative'), ('polar', 'proline')]
    for (g, h) in pairs:
        for (region, sl) in [('window', slice(None)), ('center', slice(7, 22))]:
            (b, c) = (tracks[g][:, sl], tracks[h][:, sl])
            directed = np.mean([(b[:, :-d] * c[:, d:] - c[:, :-d] * b[:, d:]).mean(1) for d in range(1, 5)], axis=0)
            add(f'{g}_then_{h}_{region}', directed, f'mean over d=1..4 of [mean_i I_{g}(i)*I_{h}(i+d)-I_{h}(i)*I_{g}(i+d)] within {region}')
    grammars = [('SP_TP', 'ST', 'P', 1), ('PxxP', 'P', 'P', 3), ('RG', 'R', 'G', 1), ('STxxDE', 'ST', 'DE', 3), ('RKxxST', 'RK', 'ST', 3)]
    for (name, g, h, d) in grammars:
        for (region, sl) in [('window', slice(None)), ('center', slice(7, 22))]:
            b = np.isin(a[:, sl], list(g))
            c = np.isin(a[:, sl], list(h))
            add(f'{name}_{region}_density', (b[:, :-d] * c[:, d:]).mean(1), f'mean_i 1[aa_i in {g}, aa_(i+{d}) in {h}] within {region}', 'fine_sequence_grammar_may_overlap_existing_motif_regex')
    for aa in AA:
        add(f'aa_{aa}_fraction', (a == aa).mean(1), f'count({aa})/29', 'composition_control')
    if not all((np.isfinite(x).all() for x in values.values())):
        raise ValueError('Nonfinite descriptors')
    return (pd.DataFrame(values), pd.DataFrame(definitions))
