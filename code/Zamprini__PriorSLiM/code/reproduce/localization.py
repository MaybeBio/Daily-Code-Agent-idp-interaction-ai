"""Top-K localization with strict >9-residue center spacing."""
import hashlib
import numpy as np
def require(ok, message):
    if not ok:
        raise ValueError(message)

def nms(centers, order, spacing, max_k):
    chosen = []
    for i in order:
        if all((abs(int(centers[i]) - int(centers[j])) > spacing for j in chosen)):
            chosen.append(int(i))
            if len(chosen) == max_k:
                break
    return np.asarray(chosen, dtype=int)

def hits_for(centers, chosen, starts, ends, first, last, max_k):
    hit = np.zeros((max_k, len(starts)), dtype=bool)
    n = len(chosen)
    if n:
        lo = np.maximum(first, centers[chosen] - 4)
        hi = np.minimum(last, centers[chosen] + 4)
        hit[:n] = np.maximum.accumulate((lo[:, None] <= ends) & (hi[:, None] >= starts), axis=0)
        hit[n:] = hit[n - 1]
    return (hit, np.minimum(np.arange(1, max_k + 1), n))

def evaluate_scan(centers, scores, starts, ends, first, last, token, spacing=9, max_k=10):
    require(np.isfinite(scores).all(), 'Nonfinite score')
    order = np.argsort(-scores, kind='stable')
    selected = nms(centers, order, spacing, max_k)
    if len(selected):
        delta = np.abs(centers[selected, None] - centers[selected][None, :])
        require((delta[np.triu_indices(len(selected), 1)] > spacing).all(), 'Candidate spacing')
    unique, counts = np.unique(scores, return_counts=True)
    tied = bool(len(selected) and (counts[np.searchsorted(unique, scores[selected])] > 1).any())
    if not tied:
        hit, count = hits_for(centers, selected, starts, ends, first, last, max_k)
        return (hit.astype(float), count.astype(float), False)
    seed = int.from_bytes(hashlib.sha256(('20261103:' + token).encode()).digest()[:8], 'little')
    rng = np.random.default_rng(seed)
    hit = np.zeros((max_k, len(starts)))
    count = np.zeros(max_k)
    for _ in range(1000):
        order = np.lexsort((rng.random(len(scores)), -scores))
        selected = nms(centers, order, spacing, max_k)
        one, number = hits_for(centers, selected, starts, ends, first, last, max_k)
        hit += one / 1000
        count += number / 1000
    return (hit, count, True)

def self_test():
    c = np.array([1, 10, 11, 20, 21])
    require(nms(c, np.arange(5), 9, 10).tolist() == [0, 2, 4], 'Strict >9 boundary')
    h, n = hits_for(np.array([1, 100]), np.array([0, 1]), np.array([5, 6, 95, 96]), np.array([5, 6, 95, 96]), 1, 100, 10)
    require(h[0].tolist() == [True, False, False, False], 'N-terminal hit boundary')
    require(h[-1].tolist() == [True, False, False, True], 'C-terminal hit boundary')
    require(n[-1] == 2, 'Fewer-than-K candidate accounting')
    h, n, _ = evaluate_scan(np.array([], int), np.array([]), np.array([5]), np.array([6]), 1, 10, 'empty')
    require(not h.any() and (not n.any()), 'Empty scan')
    rng = np.random.default_rng(20261103)
    for _ in range(200):
        c = np.arange(1, int(rng.integers(2, 600)))
        z = rng.normal(size=len(c))
        chosen = nms(c, np.argsort(-z), 9, 10)
        reference = []
        remaining = list(range(len(c)))
        while remaining and len(reference) < 10:
            best = max(remaining, key=lambda j: z[j])
            reference.append(best)
            remaining = [j for j in remaining if abs(int(c[j]) - int(c[best])) >= 10]
        require(chosen.tolist() == reference, 'Independent greedy-selection replay')
    h, _, tied = evaluate_scan(np.array([1, 10]), np.array([0.5, 0.5]), np.array([1]), np.array([1]), 1, 10, 'test_ties')
    require(tied and 0.4 < h[0, 0] < 0.6, 'Seeded tie averaging')
    return 205

def run(input_dir, output_dir):
    import pandas as pd
    self_test()
    rows = []
    with np.load(input_dir / 'localization.npz', allow_pickle=False) as z:
        recovery = z['recovery']
        for mi, method in enumerate(z['methods']):
            for fi, fold in enumerate(z['folds']):
                for oi, mode in enumerate(z['modes']):
                    for ki, k in enumerate(z['k']):
                        a = recovery[mi,fi,oi,:,ki]
                        rows.append(dict(method=str(method),fold=int(fold),mode=str(mode),k=int(k),recovered=a.sum(),loci=len(a),recall=a.mean()))
    frame = pd.DataFrame(rows)
    frame.to_csv(output_dir / 'pooled_per_fold.tsv',sep='\t',index=False)
    frame.groupby(['method','mode','k']).agg(recall_mean=('recall','mean'),recall_sd=('recall','std'),folds=('fold','nunique'),loci=('loci','first')).reset_index().to_csv(output_dir / 'pooled_summary.tsv',sep='\t',index=False)
