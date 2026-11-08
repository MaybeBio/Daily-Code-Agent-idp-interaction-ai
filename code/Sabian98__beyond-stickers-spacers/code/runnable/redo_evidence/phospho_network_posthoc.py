"""
phospho_network_posthoc.py — POST-HOC robustness of P1b of phospho_network_plan.md (hash 5b113f57), added after the
primary result was known. Concern: sites whose weights vary more (lower effective sample fraction) produce larger apparent
changes partly from reweighting noise, and if annotated sites differ in that respect the result could be an artefact.
Adds within-protein z(log effective sample fraction) to the P1b covariates, overall and for experimental annotations.
  OMP_NUM_THREADS=4 ~/.venvs/bender/bin/python -I phospho_network_posthoc.py  ->  phospho_network/posthoc.txt
"""
import sys, os
sys.argv = [sys.argv[0]]
import numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'phospho_network_test.py')).read().split("allm = np.ones(len(R), bool)")[0])
L.clear()
logess = np.log(ess)
def p1b_ess(metric, mask, seed=4):
    G = groups(mask)
    z = zwithin(metric, G)
    cov = [zwithin(np.abs(dRg), G), zwithin(logess, G)] + [zwithin(f(k), G) for k in ('local_KR', 'local_DE', 'pro_plus1', 'term_dist')]
    idx = np.concatenate(G)
    X = np.nan_to_num(np.column_stack([c[idx] for c in cov]))
    A = np.column_stack([np.ones(len(idx)), X])
    res = np.full(len(z), np.nan)
    res[idx] = z[idx] - A @ np.linalg.lstsq(A, z[idx], rcond=None)[0]
    return perm_mean(res, ann & mask, G, seed=seed)
Rhub = f('R_hub')
G = groups(np.ones(len(R), bool))
import scipy.stats as st
p('POST-HOC: P1b with within-protein z(log effective sample fraction) added as a covariate')
p(f'  Spearman(log ESS, R_hub) over all sites {st.spearmanr(logess, Rhub).correlation:+.3f}; '
  f'annotated vs other median ESS {np.median(ess[ann]):.3f} vs {np.median(ess[~ann]):.3f}')
o, n, q = p1b_ess(Rhub, np.ones(len(R), bool))
p(f'  all annotations:          mean residual {o:+.4f} (null {n:+.4f}), P {q:.4f}')
o, n, q = p1b_ess(Rhub, ~(ann & ~expv))
p(f'  experimental annotations: mean residual {o:+.4f} (null {n:+.4f}), P {q:.4f}')
open(os.path.join(OUT, 'posthoc.txt'), 'w').write('\n'.join(L) + '\n')
