"""idprg: the distribution of the radius of gyration of an intrinsically disordered protein from its sequence, as
sampled by HPS1 coarse-grained simulations at 298 K and 150 mM salt.

    import idprg
    idprg.predict(["MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"])

returns <Rg>, sigma, skewness, excess kurtosis and the parameters (ln A, S_m, nu) of the predicted two-exponent distribution; the
density itself is idprg.distribution.pdf(rg, lnA, ln S_m, nu). closed_form=True uses the 21-parameter closed form
instead of the network.
"""
from . import distribution  # noqa: F401  (first: enables float64 in JAX)
import numpy as np
import pandas as pd

from . import closed_form as _cf, descriptors as _desc, network as _net

__all__ = ["predict"]
_CACHE = {}


def predict(sequences, closed_form=False):
    """One row per sequence: N, mean_rg, std_rg (nm), skew, exkurt, nu, S_m (nm, the most probable Rg), lnA."""
    if isinstance(sequences, str):
        sequences = [sequences]
    X = _desc.table(sequences)
    if closed_form:
        coef = _CACHE.setdefault("cf", _cf.load())
        lnA, lnSm, u = _cf.parameters(coef, X)
    else:
        net, spec = _CACHE.setdefault("net", _net.load())
        lnA, lnSm, u = _net.parameters(net, spec, X)
    nu = distribution.nu_of(u)
    m, s, k, e = distribution.moments(lnA, lnSm, nu)
    return pd.DataFrame(dict(N=X.N, mean_rg=m, std_rg=s, skew=k, exkurt=e, nu=nu, S_m=np.exp(lnSm), lnA=lnA))
