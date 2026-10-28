"""
Fisher information rate from the full 2D two-sphere BVP.

Parallelised axis: alphas = phi_inf grid.  For each phi_inf the target solves
the axisymmetric two-sphere BVP at a grid of separations d, records the
pair potential U_full(d) = F(d) - F(d_ref) and the surface volume fraction
Phi_surf(d) at sphere 1's south pole (outward-facing, on the symmetry axis,
toward the z end-cap), and integrates them against the polymer P_eq(d) to
obtain phi_probe.  The Fisher rate is assembled across phi_inf in the
dashboard (it needs the full array).

Driver: sim/sweep.py.  Pipe builder: runPY.
"""

import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from viz_util.bvp_solver import Parameters
from viz_util.bvp_two_spheres import (
    solve_two_spheres,
    excess_free_energy_two_spheres,
    _phi_on_sphere,
)


blurb = (
    "Fisher information rate from the full 2D two-sphere BVP.  Sweeps "
    "phi_inf (alphas) at fixed (omega_b, R); for each phi_inf solves the "
    "2D axisymmetric BVP on a d-grid, records U_full(d) and Phi_surf(d) "
    "at sphere 1's outward-facing pole, and integrates the polymer P_eq(d) "
    "to get phi_probe.  Fisher rate is computed downstream by the dashboard."
)

fixed = dict(
    # Physical (matches bvp_two_spheres._run_main operating point)
    nb=1.0,
    chi_b=4.0,
    kappa=1.0,
    L=20.0,
    eps=1e-6,
    conserve_volume=False,
    parallel_plates=False,
    # d-grid (surface-to-surface separations)
    d_min_factor=0.05,
    d_max_factor=3.0,
    n_d=20,
    d_ref_factor=8.0,
    # 2D BVP solver
    Nrho=100,
    dz_target=0.05,
    rho_max=10.0,
    z_buf=4.0,
    f_tol=1e-6,
    maxiter=200,
    # Polymer (match fisher_rate_from_profile)
    s=150.0,
    a=1.0,
)

alphas = np.linspace(0.001, 0.01, 20)     # (parallelised)
betas =  np.array([-1.6])
gammas = np.array([5.0])               
# betas  = np.array([-1.5])             

nparams = len(alphas)


_PARAM_KEYS = {"nb", "chi_b", "kappa", "L", "eps", "conserve_volume", "parallel_plates"}


def _param_kwargs(fixed):
    return {k: v for k, v in fixed.items() if k in _PARAM_KEYS}


def target(res, keya, alpha, keyb, beta, keyc, gamma, fixed, blurb):
    t0 = time.time()

    omega_b = beta
    phi_inf = alpha
    R       = gamma

    p = Parameters(omega_b=omega_b, R=R, phi_inf=phi_inf, **_param_kwargs(fixed))

    d_min = fixed["d_min_factor"] * R
    d_max = fixed["d_max_factor"] * R
    n_d   = int(fixed["n_d"])
    d_grid = np.linspace(d_min, d_max, n_d)
    d_ref  = fixed["d_ref_factor"] * R

    Nrho      = int(fixed["Nrho"])
    dz_target = float(fixed["dz_target"])
    rho_max   = float(fixed["rho_max"])
    z_buf     = float(fixed["z_buf"])
    f_tol     = float(fixed["f_tol"])
    maxiter   = int(fixed["maxiter"])

    def _Nz(dd):
        return int((2.0 * R + dd + 2.0 * z_buf) / dz_target) + 1

    solve_kw = dict(rho_max=rho_max, z_buf=z_buf, Nrho=Nrho,
                    f_tol=f_tol, maxiter=maxiter, verbose=False)

    # Reference far-separation free energy.
    res_ref = solve_two_spheres(p, d_ref, Nz=_Nz(d_ref), **solve_kw)
    F_ref = float(excess_free_energy_two_spheres(res_ref)) if res_ref["success"] else float("nan")

    U_full   = np.full(n_d, np.nan)
    Phi_surf = np.full(n_d, np.nan)
    F_arr    = np.full(n_d, np.nan)
    south = np.array([np.pi])

    for k, dd in enumerate(d_grid):
        r = solve_two_spheres(p, float(dd), Nz=_Nz(dd), **solve_kw)
        if not r["success"]:
            continue
        F = float(excess_free_energy_two_spheres(r))
        F_arr[k]    = F
        U_full[k]   = F - F_ref if np.isfinite(F_ref) else np.nan
        Phi_surf[k] = float(_phi_on_sphere(r["phi"], r["geom"],
                                           r["geom"].z_c1, south)[0])

    ok = np.isfinite(U_full) & np.isfinite(Phi_surf)
    if int(ok.sum()) >= 2:
        d_ok   = d_grid[ok]
        U_ok   = U_full[ok]
        Phi_ok = Phi_surf[ok]
        s = float(fixed["s"])
        a = float(fixed["a"])
        U_ent  = 3.0 * (d_ok + R) ** 2 / (2.0 * s * a ** 2)
        P_eq_un = (d_ok + R) ** 2 * np.exp(-U_ent - U_ok)
        denom = float(np.sum(P_eq_un))
        if denom > 0.0 and np.isfinite(denom):
            P_eq = P_eq_un / denom
            phi_probe = float(np.sum(P_eq * Phi_ok))
            success = bool(ok.all() and np.isfinite(F_ref))
        else:
            phi_probe = float("nan")
            success = False
    else:
        phi_probe = float("nan")
        success = False

    res.update(
        blurb=blurb,
        phi_inf=float(phi_inf),
        omega_b=float(omega_b),
        R=float(R),
        chi_b=float(fixed["chi_b"]),
        s=float(fixed["s"]),
        a=float(fixed["a"]),
        d_grid=d_grid.tolist(),
        d_ref=float(d_ref),
        F_ref=F_ref,
        F=F_arr.tolist(),
        U_full=U_full.tolist(),
        Phi_surf=Phi_surf.tolist(),
        phi_probe=phi_probe,
        success=success,
        n_d_failed=int((~ok).sum()),
        runtime_min=float((time.time() - t0) / 60.0),
    )
    return res
