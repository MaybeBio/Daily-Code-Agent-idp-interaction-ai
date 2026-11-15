"""Recompute the 360 primary IDPPI tests from per-record peak features."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import optimize, special, stats
import robust_glm as base
RESULTS, AUDITS, CHECKS = [], [], []
FILE1="IDPPI_1_conditional_folded.tsv"
FILE2="IDPPI_2_remains_disordered.tsv"
def explicit_cov(scores, bread, d):
    n, k = scores.shape
    codes, _ = pd.factorize(pd.concat([d.accession, d.partner_accession], ignore_index=True))
    aa, bb = (codes[:n], codes[n:])
    meat = np.zeros((k, k))
    for i in range(0, n, 256):
        left, right = (aa[i:i + 256, None], bb[i:i + 256, None])
        shared = (left == aa) | (left == bb) | (right == aa) | (right == bb)
        meat += scores[i:i + 256].T @ (shared @ scores)
    return bread @ meat @ bread * n / (n - k)

def collect(fold, rule, analysis, contrast, primary, rows):
    for r in rows:
        z = dict(r, fold=fold, min_run=rule, analysis=analysis, primary_family=primary)
        if contrast is not None:
            z['contrast'] = contrast
        z['effect_type'] = 'IRR' if z.get('family') == 'poisson' else 'OR'
        if z.get('covariance_psd') is False:
            z.update(status='not_estimable', reason='non-PSD dyadic covariance', p_value=np.nan, effect_ratio=np.nan, ci_low=np.nan, ci_high=np.nan)
        RESULTS.append(z)

def glm(fold, rule, name, d, X, y, cs, primary=False, family='binomial', offset=None):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    reason = None
    if len(d) < max(10, 2 * X.shape[1] + 1) or np.linalg.matrix_rank(X) < X.shape[1]:
        reason = 'too few or rank deficient'
    if reason is None and family == 'binomial':
        Z = (2 * y - 1)[:, None] * X
        lp = optimize.linprog(-Z.sum(0), A_ub=-Z, b_ub=np.zeros(len(d)), bounds=[(-1, 1)] * X.shape[1], method='highs')
        if not lp.success or -lp.fun > 1e-07:
            reason = 'complete/quasi separation or LP failure'
    base.RESULTS.clear()
    base.AUDITS.clear()
    base.FAILED.clear()
    if reason:
        rows = [dict(contrast=c, n=len(d), family=family, status='not_estimable', reason=reason, p_value=np.nan) for c in cs]
    else:
        base.fit(name, name, d, X, y, cs, family=family, offset=offset)
        rows = [dict(r) for r in base.RESULTS]
        AUDITS.extend((dict(r, fold=fold, min_run=rule, analysis=name) for r in base.AUDITS))
    collect(fold, rule, name, None, primary, rows)
    if fold == 1 and rule in (3, 5) and (name == 'R_categories_full') and all((r['status'] == 'estimated' for r in rows)):
        beta = np.zeros(X.shape[1])
        for _ in range(100):
            mu = special.expit(X @ beta)
            gradient = X.T @ (mu - y)
            if np.max(abs(gradient)) < 1e-07:
                break
            h = X.T @ ((mu * (1 - mu))[:, None] * X)
            step = np.linalg.solve(h, gradient)
            alpha = 1.0
            old = np.sum(np.logaddexp(0, X @ beta) - y * (X @ beta))
            while alpha > 1e-10:
                trial = beta - alpha * step
                if np.sum(np.logaddexp(0, X @ trial) - y * (X @ trial)) <= old + 1e-08:
                    beta = trial
                    break
                alpha /= 2
        assert np.max(abs(X.T @ (special.expit(X @ beta) - y))) < 1e-05
        mu = special.expit(X @ beta)
        bread = np.linalg.inv(X.T @ ((mu * (1 - mu))[:, None] * X))
        cov = explicit_cov(X * (y - mu)[:, None], bread, d)
        for r in rows:
            c = np.asarray(cs[r['contrast']])
            np.testing.assert_allclose([c @ beta, np.sqrt(c @ cov @ c)], [r['log_effect'], r['robust_se']], atol=1e-07)
        CHECKS.append(dict(fold=fold, min_run=rule, analysis=name, verified='independent Newton and shared-node covariance'))
def categories(fold, rule, name, one, two, count, union=False, primary=False, partner=False):
    d = pd.concat([one, two], ignore_index=True)
    n = d[count]
    cols = [np.ones(len(d)), n.eq(0), n.ge(2), np.log(d.length)]
    if union:
        cols.append(np.log(d.union_length))
    if partner:
        cols.append(d.partner_type.eq('IDR'))
    X = np.column_stack(cols)
    c0 = np.zeros(X.shape[1])
    c2 = c0.copy()
    c0[1] = 1
    c2[2] = 1
    cs = {'file2_zero_vs_one': c0, 'file2_multi_vs_one': c2}
    glm(fold, rule, name, d, X, d.source_file.eq(FILE2), cs, primary)
    if fold == 1 and rule in (3, 5) and (name == 'file2_categories_full'):
        y = d.source_file.eq(FILE2).to_numpy(float)
        beta = np.zeros(X.shape[1])
        for _ in range(100):
            mu = special.expit(X @ beta)
            grad = X.T @ (mu - y)
            if np.max(abs(grad)) < 1e-07:
                break
            step = np.linalg.solve(X.T @ ((mu * (1 - mu))[:, None] * X), grad)
            alpha = 1.0
            old = np.sum(np.logaddexp(0, X @ beta) - y * (X @ beta))
            while alpha > 1e-10:
                trial = beta - alpha * step
                if np.sum(np.logaddexp(0, X @ trial) - y * (X @ trial)) <= old + 1e-08:
                    beta = trial
                    break
                alpha /= 2
        mu = special.expit(X @ beta)
        assert np.max(abs(X.T @ (mu - y))) < 1e-05
        bread = np.linalg.inv(X.T @ ((mu * (1 - mu))[:, None] * X))
        cov = explicit_cov(X * (y - mu)[:, None], bread, d)
        for row in RESULTS[-2:]:
            c = cs[row['contrast']]
            np.testing.assert_allclose([c @ beta, np.sqrt(c @ cov @ c)], [row['log_effect'], row['robust_se']], atol=1e-07)
        CHECKS.append(dict(fold=fold, min_run=rule, analysis=name, verified='independent Newton and explicit shared-node covariance'))

def between(fold, rule, name, one, two, mode, primary=False):
    one = one.copy()
    two = two.copy()
    one['count'] = one[mode + '_in']
    two['count'] = two.full_n
    one['exposure'] = one.union_length
    two['exposure'] = two.length
    d = pd.concat([one, two], ignore_index=True)
    X = np.column_stack([np.ones(len(d)), d.source_file.eq(FILE1), np.log(d.length)])
    glm(fold, rule, name, d, X, d['count'], {'file1_union_vs_file2_full': [0, 1, 0]}, primary, 'poisson', np.log(d.exposure))

def paired(fold, rule, name, one, mode, primary=False):
    eligible = one[one.outside_length.gt(0)]
    tot = eligible[mode + '_in'] + eligible[mode + '_out']
    d = eligible[tot.gt(0)].copy()
    total = tot[tot.gt(0)].to_numpy()
    y = d[mode + '_in'].to_numpy()
    z = dict(n=len(eligible), n_cohort=len(one), n_informative=len(d), n_no_outside=len(one) - len(eligible), family='poisson')
    if len(d) < 10 or y.sum() == 0 or (total - y).sum() == 0:
        z.update(status='not_estimable', reason='too few informative pairs or boundary rate', p_value=np.nan)
    else:
        offset = np.log(d.union_length.to_numpy() / d.outside_length.to_numpy())
        beta = optimize.brentq(lambda b: np.sum(y - total * special.expit(offset + b)), -40, 40, xtol=1e-13)
        p = special.expit(offset + beta)
        scores = (y - total * p)[:, None]
        bread = np.array([[1 / np.sum(total * p * (1 - p))]])
        cov = base.dyadic_cov(scores, bread, d.accession, d.partner_accession)
        if cov[0, 0] <= 0:
            z.update(status='not_estimable', reason='nonpositive dyadic variance', p_value=np.nan)
        else:
            se = np.sqrt(cov[0, 0])
            z.update(status='estimated', log_effect=beta, robust_se=se, effect_ratio=np.exp(beta), ci_low=np.exp(beta - 1.96 * se), ci_high=np.exp(beta + 1.96 * se), p_value=2 * stats.norm.sf(abs(beta / se)))
            if fold in (1, 10) and name == 'file1_union_vs_outside':
                test = optimize.minimize_scalar(lambda b: np.sum(total * np.logaddexp(0, offset + b) - y * (offset + b)), bounds=(-20, 20), method='bounded', options={'xatol': 1e-12})
                np.testing.assert_allclose(test.x, beta, atol=1e-06)
                np.testing.assert_allclose(explicit_cov(scores, bread, d), cov, atol=1e-10)
                CHECKS.append(dict(fold=fold, min_run=rule, analysis=name, verified='independent scalar likelihood and explicit shared-node covariance'))
    collect(fold, rule, name, 'file1_union_vs_file1_outside', primary, [z])

def run(input_dir, output_dir):
    RESULTS.clear()
    all_features = pd.read_csv(input_dir / "structural.tsv.gz", sep="\t")
    for fold in range(1, 11):
        for rule in (1, 3, 5):
            frame = all_features[all_features.fold.eq(fold) & all_features.min_run.eq(rule)].copy()
            for stratum in ("pooled", "Domain", "IDR"):
                d = frame if stratum == "pooled" else frame[frame.partner_type.eq(stratum)]
                one, two = d[d.source_file.eq(FILE1)], d[d.source_file.eq(FILE2)]
                union = one[one.file1_contact_eligible]
                start = len(RESULTS)
                paired(fold, rule, "file1_union_vs_outside", union, "apex", True)
                between(fold, rule, "file1_union_vs_file2_full", union, two, "apex", True)
                categories(fold, rule, "file2_categories_full", one, two, "full_n", primary=True)
                for row in RESULTS[start:]:
                    row["stratum"] = stratum
    result = pd.DataFrame(RESULTS)
    assert len(result) == 360
    p = result.p_value.fillna(1).to_numpy()
    order = np.argsort(p)
    q, h = np.empty(360), np.empty(360)
    q[order] = np.minimum.accumulate((p[order] * 360 / np.arange(1,361))[::-1])[::-1].clip(0,1)
    h[order] = np.maximum.accumulate(p[order] * np.arange(360,0,-1)).clip(0,1)
    q[result.p_value.isna()] = np.nan
    h[result.p_value.isna()] = np.nan
    result["q_bh_primary360"], result["p_holm_primary360"] = q, h
    result.to_csv(output_dir / "primary360_tests.tsv", sep="\t", index=False)
