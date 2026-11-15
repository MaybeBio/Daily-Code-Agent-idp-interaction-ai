"""Dyadic protein-robust generalized linear models."""
import numpy as np
import pandas as pd
from scipy import optimize, special, stats
RESULTS, AUDITS, FAILED = [], [], []
def meat(scores, labels):
    code, unique = pd.factorize(pd.Series(list(labels), dtype=object), sort=False)
    total = np.zeros((len(unique), scores.shape[1]))
    np.add.at(total, code, scores)
    return total.T @ total

def dyadic_cov(scores, bread, a, b):
    """Sandwich allowing arbitrary dependence whenever two rows share a protein.

    Sum node-cluster meats, subtract distinct unordered-pair meats: pairs
    sharing both nodes are then counted once. Self-edges contribute one node.
    HC1 n/(n-k) correction; asymptotic normal inference, not independent windows.
    """
    node_scores, nodes, pair_scores, pairs = ([], [], [], [])
    for s, pa, pb in zip(scores, a, b):
        for p in sorted({pa, pb}):
            nodes.append(p)
            node_scores.append(s)
        if pa != pb:
            pairs.append(tuple(sorted((pa, pb))))
            pair_scores.append(s)
    m = meat(np.array(node_scores), nodes)
    if pairs:
        m -= meat(np.array(pair_scores), pairs)
    n, k = scores.shape
    cov = bread @ m @ bread * n / (n - k)
    return (cov + cov.T) / 2

def fit(view, name, d, X, y, contrasts, *, family='binomial', offset=None):
    """Store dyadic robust contrasts; reject nonidentified/separated designs."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    off = np.zeros(n) if offset is None else np.asarray(offset, float)

    def fail(reason, **extra):
        FAILED.append(dict(view=view, model=name, n=n, reason=reason, **extra))
        for label in contrasts:
            RESULTS.append(dict(view=view, model=name, contrast=label, n=n, family=family, status='not_estimable', reason=reason, p_value=np.nan, covariance_psd=np.nan))
    if n <= k or np.linalg.matrix_rank(X) < k:
        fail('rank deficient or too few observations')
        return

    def calc(beta):
        eta = X @ beta + off
        if family == 'binomial':
            mu = special.expit(eta)
            return (np.sum(np.logaddexp(0, eta) - y * eta), X.T @ (mu - y), mu, mu * (1 - mu))
        if family == 'poisson':
            mu = np.exp(np.clip(eta, -40, 40))
            return (np.sum(mu - y * eta), X.T @ (mu - y), mu, mu)
        return (np.sum((y - eta) ** 2) / 2, X.T @ (eta - y), eta, np.ones(n))
    opt = optimize.minimize(lambda beta: calc(beta)[:2], np.zeros(k), jac=True, method='BFGS', options={'gtol': 1e-08, 'maxiter': 2000})
    beta = opt.x
    root = optimize.root(lambda b: calc(b)[1], beta, jac=lambda b: X.T @ (calc(b)[3][:, None] * X), tol=1e-10)
    if np.max(np.abs(calc(root.x)[1])) < np.max(np.abs(calc(beta)[1])):
        beta = root.x
    _, grad, mu, wt = calc(beta)
    h = X.T @ (wt[:, None] * X)
    if np.max(np.abs(grad)) > 1e-05 or np.linalg.cond(h) > 1000000000000.0 or np.max(np.abs(beta)) > 25:
        fail('nonconvergence, separation or ill-conditioned Hessian', gradient=float(np.max(np.abs(grad))), coefficients=beta.tolist())
        return
    bread = np.linalg.inv(h)
    cov = dyadic_cov(X * (y - mu)[:, None], bread, d.accession, d.partner_accession)
    eigmin = float(np.linalg.eigvalsh(cov).min())
    AUDITS.append(dict(view=view, model=name, family=family, n=n, k=k, n_proteins=len(set(d.accession) | set(d.partner_accession)), gradient=float(np.max(np.abs(grad))), hessian_condition=float(np.linalg.cond(h)), robust_cov_min_eigenvalue=eigmin, coefficients=beta.tolist()))
    for label, c in contrasts.items():
        c = np.asarray(c, float)
        var = float(c @ cov @ c)
        effect = float(c @ beta)
        if var <= 0:
            FAILED.append(dict(view=view, model=name, contrast=label, n=n, reason='nonpositive robust contrast variance'))
            RESULTS.append(dict(view=view, model=name, contrast=label, n=n, family=family, status='not_estimable', reason='nonpositive robust contrast variance', p_value=np.nan, covariance_psd=False))
            continue
        se = np.sqrt(var)
        RESULTS.append(dict(view=view, model=name, contrast=label, n=n, family=family, effect_ratio=float(np.exp(effect)), ci_low=float(np.exp(effect - 1.96 * se)), ci_high=float(np.exp(effect + 1.96 * se)), log_effect=effect, robust_se=float(se), p_value=float(2 * stats.norm.sf(abs(effect / se))), model_se_p=float(2 * stats.norm.sf(abs(effect / np.sqrt(c @ bread @ c)))), covariance='dyadic protein robust HC1', covariance_psd=eigmin >= -1e-10, status='estimated'))
