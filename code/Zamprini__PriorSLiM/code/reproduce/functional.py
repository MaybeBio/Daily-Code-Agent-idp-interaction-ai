"""Recompute 43 functional-annotation tests from matched fold-wise hits."""
import json
import numpy as np
import pandas as pd
from scipy import stats
def edge_weights(frame):
    unique = frame[['fg_cluster', 'uniprot', 'unit_id']].drop_duplicates()
    nc = unique.fg_cluster.nunique()
    np_per = unique[['fg_cluster', 'uniprot']].drop_duplicates().groupby('fg_cluster').size()
    na_per = unique.groupby('uniprot').size()
    nk = frame.groupby('unit_id').size()
    w = 1 / (nc * frame.fg_cluster.map(np_per).to_numpy() * frame.uniprot.map(na_per).to_numpy() * frame.unit_id.map(nk).to_numpy())
    assert np.isclose(w.sum(), 1)
    return w

def uncertainty(delta, weights, fg_clusters, bg_clusters):
    estimate = float(np.dot(delta, weights))
    influence = weights * (delta - estimate)
    variances = []
    counts = []
    for labels in [fg_clusters, bg_clusters]:
        unique, indices = np.unique(np.asarray(labels), return_inverse=True)
        n = len(unique)
        counts.append(n)
        sums = np.bincount(indices, weights=influence)
        variances.append(n / (n - 1) * np.dot(sums, sums) if n > 1 else np.nan)
    if min(counts) < 2:
        return dict(se=np.nan, ci_low=np.nan, ci_high=np.nan, p_value=np.nan, df=min(counts) - 1)
    se = float(np.sqrt(sum(variances)))
    df = min(counts) - 1
    critical = stats.t.ppf(0.975, df)
    p = float(2 * stats.t.sf(abs(estimate / se), df)) if se else float(estimate == 0)
    return dict(se=se, ci_low=estimate - critical * se, ci_high=estimate + critical * se, p_value=p, df=df)

def bh(pvalues, planned_tests=None):
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    m = planned_tests if planned_tests is not None else len(p)
    assert m >= len(p)
    q = np.minimum.accumulate((p[order] * m / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1)
    return out
def summarise(pairs, a, b, planned, strata=True):
    rows = []
    subsets = {'all': np.ones(len(pairs), bool)}
    if strata:
        subsets.update(experimental=pairs.experimental.eq(1), no_other_core_overlap=pairs.other_core_overlap.eq(0), length_1_7=pairs.feature_length.le(7), length_8_15=pairs.feature_length.between(8, 15), length_16_29=pairs.feature_length.between(16, 29), length_30_59=pairs.feature_length.between(30, 59), length_ge60=pairs.feature_length.ge(60))
        subsets.update({s: pairs.source_stratum.eq(s) for s in ['B_foreground', 'A_extension']})
        subsets.update({s: pairs.taxon.eq(s) for s in ['Metazoa','Viridiplantae','Fungi','Other_eukaryota']})
    for stratum, mask in subsets.items():
        for (method, sub), g in pairs.loc[mask].groupby(['evaluation_method', 'subclass']):
            w = edge_weights(g)
            ix = g.index.to_numpy()
            fg, bg = (a[ix], b[ix])
            rows.append(dict(evaluation_method=method, subclass=sub, stratum=stratum, annotations=g.unit_id.nunique(), proteins=g.uniprot.nunique(), idrs=len(g[['uniprot', 'idr_start', 'idr_end']].drop_duplicates()), fg_clusters=g.fg_cluster.nunique(), bg_clusters=g.bg_cluster.nunique(), pairs=len(g), observed=float(w @ fg), background=float(w @ bg), delta=float(w @ (fg - bg)), **uncertainty(fg - bg, w, g.fg_cluster, g.bg_cluster)))
    f = pd.DataFrame(rows)
    f['q_joint'] = np.nan
    f['q_within_method'] = np.nan
    mask = f.stratum.eq('all') & f.p_value.notna()
    f.loc[mask, 'q_joint'] = bh(f.loc[mask, 'p_value'], len(planned))
    for method in ['site', 'region']:
        m = mask & f.evaluation_method.eq(method)
        f.loc[m, 'q_within_method'] = bh(f.loc[m, 'p_value'], sum((x['evaluation_method'] == method for x in planned)))
    return f

def run(input_dir, output_dir):
    pairs = pd.read_csv(input_dir / "pairs.tsv.gz", sep="\t")
    planned = json.loads((input_dir / "planned_tests.json").read_text())
    with np.load(input_dir / "hits.npz", allow_pickle=False) as z:
        assert np.array_equal(z["unit_ids"], pairs.unit_id.to_numpy(str))
        a, b = z["observed_hit"], z["background_hit"]
    assert a.shape == b.shape == (len(pairs), 10)
    per = []
    for j in range(10):
        f = summarise(pairs, a[:,j].astype(float), b[:,j].astype(float), planned, strata=False)
        f["fold"] = j + 1
        per.append(f)
    per = pd.concat(per, ignore_index=True)
    main = summarise(pairs, a.mean(axis=1), b.mean(axis=1), planned, strata=False)
    means = per.groupby(["evaluation_method","subclass","stratum"]).agg(positive_folds=("delta",lambda x:int((x>0).sum())),negative_folds=("delta",lambda x:int((x<0).sum())),delta_min=("delta","min"),delta_max=("delta","max")).reset_index()
    main.merge(means,on=["evaluation_method","subclass","stratum"],validate="one_to_one").to_csv(output_dir / "main_results.tsv",sep="\t",index=False)
    per.to_csv(output_dir / "per_fold_summary.tsv",sep="\t",index=False)
