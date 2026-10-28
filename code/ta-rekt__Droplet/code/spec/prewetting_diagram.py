"""
3D prewetting phase diagram spec: scan (omega_b, R, phi_inf) at fixed chi_b.

Parallelized axis: omega_bs. Inner axes: Rs x phi_infs.
Driver: sim/sweep.py. Pipe builder: runPY.
"""

import numpy as np

from viz_util.bvp_solver import (
    Parameters,
    solve_bvp_robust,
    excess_free_energy,
    compute_mu_inf,
    compute_spinodal,
    compute_binodal,
)


blurb = (
    "3D prewetting phase diagram: scan (omega_b, R, phi_inf) at fixed "
    "(chi_b, kappa, L, nb). Records phi(z), dphi(z), excess grand potential, "
    "scalar summaries, and spinodal/binodal bounds at each grid point."
)

fixed = dict(
    nb=1.0,
    chi_b=4,
    kappa=1.0,
    L=5.0,
    eps=1e-8,
    conserve_volume=False,
    parallel_plates=False,
    z_grid_n=200,
)

alphas = np.linspace(0.01, 0.018, 50)
betas = np.array([5])
gammas = np.linspace(-1.4, -1.6, 50)
# alphas = np.array([0.01])
# gammas = np.array([-1.5])

nparams = len(alphas)


_PARAM_KEYS = {"nb", "chi_b", "kappa", "L", "eps", "conserve_volume", "parallel_plates"}


def _param_kwargs(fixed):
    return {k: v for k, v in fixed.items() if k in _PARAM_KEYS}


def target(res, keya, alpha, keyb, beta, keyg, gamma, fixed, blurb):

    phi_inf = alpha
    R = beta
    omega_b = gamma

    p = Parameters(omega_b=omega_b, R=R, phi_inf=phi_inf, **_param_kwargs(fixed))

    sol = solve_bvp_robust(p, verbose=False)
    z = np.linspace(0.0, p.L, fixed["z_grid_n"])

    success = sol is not None and getattr(sol, "success", True)
    if success:
        y = sol.sol(z)
        phi, dphi = y[0], y[1]
        try:
            dOmega = float(excess_free_energy(phi, dphi, z, p))
        except Exception:
            dOmega = float("nan")
            success = False
    else:
        phi = np.full_like(z, np.nan)
        dphi = np.full_like(z, np.nan)
        dOmega = float("nan")

    phi_spin_low, phi_spin_high = compute_spinodal(p)
    phi_bin_low,  phi_bin_high  = compute_binodal(p)

    res.update(
        blurb=blurb,
        omega_b=float(omega_b),
        R=float(R),
        phi_inf=float(phi_inf),
        success=bool(success),
        mu_inf=float(compute_mu_inf(p)),
        phi_0=float(phi[0]),
        phi_L=float(phi[-1]),
        excess_grand_potential=dOmega,
        # z=z.tolist(),
        # phi=phi.tolist(),
        # dphi=dphi.tolist(),
        phi_spin_low=phi_spin_low,
        phi_spin_high=phi_spin_high,
        phi_bin_low=phi_bin_low,
        phi_bin_high=phi_bin_high,
    )
    return res
