import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb
import fig02_abc as abc
import fig02_detail as detail
ABC = OLD = LATENT = 'data/fig02'
RC = dict(abc.RC, **{'savefig.facecolor': 'white'})
PANELS = 'abcdefg'

def load_abc(ctx):
    wide = ctx.table(ABC + '/aligned_test_predictions.tsv')
    y = wide.label.to_numpy(int)
    assert len(y) == 9765 and y.sum() == 1395
    models = []
    for method in abc.METHODS:
        for fold in range(1, 11) if method in ['PriorSLiM', 'SHARK'] else [0]:
            models.append({'method': method, 'fold': fold, 'score': wide[f'{method}_{fold:02d}_score'].to_numpy(float)})
    return (y, models, ctx.table(ABC + '/a_metrics.tsv'), ctx.table(ABC + '/b_negative_types.tsv'), ctx.table(ABC + '/c_within_motif_auc.tsv'))

def load_detail(ctx, letter):
    read = lambda name: ctx.table(OLD + '/' + name + '.tsv')
    if letter == 'e':
        return {'effects': read('e_feature_effects')}
    if letter == 'g':
        d = {'allg': read('g_all_candidate_dependence'), 'broad': read('g_broad_candidate_evidence'), 'gfolds': ctx.table(LATENT + '/all_candidate_fold_metrics.tsv')}
        assert len(d['allg']) == 5880 and len(d['broad']) == 1459 and (d['broad'].near_zero_stable.sum() == 26)
        return d
    discovery = ctx.table(LATENT + '/discovery_f_candidates.tsv')
    locked = ctx.json(LATENT + '/locked_candidates.json')
    cases = ctx.table(LATENT + '/f_cases.tsv')
    tracks = ctx.table(LATENT + '/f_position_profiles.tsv')
    selected = []
    for kind in ['FP_corrected', 'FN_rescued']:
        fid = int([r for r in locked['f'] if r['kind'] == kind][0]['feature_id'])
        eligible = discovery.query('kind == @kind and eligible').sort_values(['mean_bce', 'proteins', 'feature_id'], ascending=[False, False, True])
        assert fid == int(eligible.iloc[0].feature_id)
        selected.append(cases[(cases.kind == kind) & (cases.feature_id == fid)].head(4))
    import pandas as pd
    return {'fcases': pd.concat(selected, ignore_index=True), 'ftracks': tracks, 'locked': locked, 'discovery': discovery}
