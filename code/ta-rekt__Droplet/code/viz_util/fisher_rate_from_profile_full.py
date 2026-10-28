"""
Polymer prewetting as a concentration sensor

A. linear approximation (analytical)
    1. plotting the equilibrium distance distribution between poly beads
    2. integrating over it for psi_probe
    3. computing the fisher rate from psi_probe

B. nonlinear profile from solver
    4. pair potential U(d): Derjaguin (d < 0.5R) <-> additive (d > 0.5R)
    5. repeat 1-3 with the nonlinear U(d)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import matplotlib
from dataclasses import dataclass, replace
from scipy.interpolate import PchipInterpolator
from scipy.special import expit
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.integrate import solve_bvp
from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, excess_free_energy, compute_mu_inf,
    make_ode, make_bc, generate_initial_guess, clip01,
)

@dataclass
class Probe:
    kb: float = 1.0
    ku: float = 1.0

# ── Parameters ────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=4.0,
    omega_b=-1.0,
    kappa=1,
    phi_inf=0.1,
    L=10.0,
    eps=1e-6,
    R=5.0,
    conserve_volume=True,
    parallel_plates=False,
)

# polymer params
s = 100   # inter-bead distance
a = 1     # Kuhn length

# Geometry / numerics
nsteps   = 1000
L_ISO    = 10.0    # single-sphere / single-plate reference domain
N_ISO    = 200     # discretisation of the single-sphere profile
N_DGRID  = 40      # number of pair separations evaluated for U(d)
N_RHO    = 80      # cylindrical-integration grid (radial)
N_Z      = 120     # cylindrical-integration grid (axial)
CROSS_CENTER = 0.5  # crossover at d = CROSS_CENTER * R
CROSS_WIDTH  = 0.10 # crossover width in units of R


# ── Single-sphere psi_1(s) in open-reservoir mode ────────────────────────
def solve_isolated_sphere(p: Parameters):
    """
    Solve the spherical BVP around a single sphere with open-reservoir BC
    phi(L_iso) = phi_inf.  Returns PCHIP interpolators psi_1(s), dpsi_1(s)
    parameterized by distance s = r - R from the sphere surface, and the
    upper limit s_max beyond which psi_1 is taken to be zero.
    """
    p_iso = replace(p, L=L_ISO, conserve_volume=False, parallel_plates=False)
    sol = solve_auto(p_iso, verbose=False)
    if sol is None or not sol.success:
        raise RuntimeError(f"isolated single-sphere BVP failed at phi_inf={p.phi_inf:.4g}")

    s_grid = np.linspace(0.0, p_iso.L, N_ISO)
    phi  = sol.sol(s_grid)[0]
    dphi = sol.sol(s_grid)[1]
    psi  = phi - p.phi_inf
    psi_interp  = PchipInterpolator(s_grid, psi,  extrapolate=False)
    dpsi_interp = PchipInterpolator(s_grid, dphi, extrapolate=False)
    return psi_interp, dpsi_interp, p_iso.L


def _psi_eval(interp, s):
    """PCHIP evaluation with hard zero outside [0, s_max]."""
    s_arr  = np.asarray(s, dtype=float)
    s_max  = interp.x[-1]
    inside = (s_arr >= 0.0) & (s_arr <= s_max)
    out    = np.zeros_like(s_arr)
    if np.any(inside):
        out[inside] = interp(s_arr[inside])
    return out


# ── Plate-plate W(d) for the Derjaguin branch ───────────────────────────
def _solve_plate_lowest_F(p_plate: Parameters,
                          strategies=('high_enrichment', 'exponential', 'tanh',
                                      'depleted', 'uniform')):
    """
    Solve the open-mode parallel-plates BVP at p_plate from each initial-guess
    strategy in turn, evaluate excess_free_energy on every converged solution,
    and return (sol, F) for the lowest-F branch.  Returns (None, +inf) on
    universal failure.

    This is the right thermodynamic equilibrium for chi_b > 2, where the
    plate-plate system may admit both a dry and a capillary-condensed branch.
    """
    mu_inf = compute_mu_inf(p_plate)
    ode    = make_ode(p_plate, mu_inf)
    bc     = make_bc(p_plate)

    best_sol, best_F = None, float('inf')
    for strat in strategies:
        try:
            z, y_guess = generate_initial_guess(p_plate, strat)
            y_guess[0] = clip01(y_guess[0], p_plate.eps)
            for tol_try in (1e-6, 1e-5, 1e-4):
                try:
                    sol = solve_bvp(ode, bc, z, y_guess,
                                    tol=tol_try, max_nodes=100000, verbose=0)
                except Exception:
                    continue
                if not sol.success:
                    continue
                phi_min = float(np.min(sol.y[0]))
                phi_max = float(np.max(sol.y[0]))
                if phi_min <= 0 or phi_max >= 1:
                    continue
                z_dense    = np.linspace(0.0, p_plate.L, max(50, N_ISO // 2))
                phi_dense  = sol.sol(z_dense)[0]
                dphi_dense = sol.sol(z_dense)[1]
                F = excess_free_energy(phi_dense, dphi_dense, z_dense, p_plate)
                if F < best_F:
                    best_sol, best_F = sol, F
                break  # this strategy converged at this tol; move to next strat
        except Exception:
            continue
    return best_sol, best_F


def _warm_start_plate(p_plate: Parameters, prev_sol):
    """
    Solve the plate-plate BVP warm-started from a previous (larger-d) solution
    by rescaling its profile onto [0, p_plate.L].  Returns (sol, F) or
    (None, +inf) on failure.
    """
    mu_inf = compute_mu_inf(p_plate)
    ode    = make_ode(p_plate, mu_inf)
    bc     = make_bc(p_plate)
    z_new  = np.linspace(0.0, p_plate.L, 200)
    scale  = prev_sol.x[-1] / p_plate.L
    y_new  = prev_sol.sol(np.minimum(z_new * scale, prev_sol.x[-1] * 0.999))
    y_new[0] = clip01(y_new[0], p_plate.eps)
    for tol_try in (1e-6, 1e-5, 1e-4):
        try:
            sol = solve_bvp(ode, bc, z_new, y_new,
                            tol=tol_try, max_nodes=100000, verbose=0)
        except Exception:
            continue
        if (sol.success and float(np.min(sol.y[0])) > 0
                and float(np.max(sol.y[0])) < 1):
            z_dense    = np.linspace(0.0, p_plate.L, max(50, N_ISO // 2))
            phi_dense  = sol.sol(z_dense)[0]
            dphi_dense = sol.sol(z_dense)[1]
            F = excess_free_energy(phi_dense, dphi_dense, z_dense, p_plate)
            return sol, F
    return None, float('inf')


def derjaguin_U(d_grid: np.ndarray, p: Parameters) -> np.ndarray:
    """
    U_Derjaguin(d) = pi * R * [F_plates(d) - 2 F_single] for each d.

    F_plates is the lowest-energy branch (over multiple initial guesses) of
    the open-mode parallel-plates BVP on z in [0, d].  F_single is the
    open-mode single-slab reference on z in [0, L_ISO].

    Iteration order: largest d first, then sweep down using the previous
    larger-d solution as a warm start.  This stabilises the BVP at very
    tight separations where cold initial guesses tend to push phi out of
    (0, 1).  At each d we also try fresh initial guesses and keep whichever
    converged solution has the lowest excess free energy (proper branch
    selection through the capillary-condensation transition).
    """
    R = p.R

    p_single = replace(p, L=L_ISO, R=np.inf,
                       parallel_plates=False, conserve_volume=False)
    sol_s = solve_auto(p_single, verbose=False)
    if sol_s is None or not sol_s.success:
        raise RuntimeError(f"single-plate reference BVP failed at phi_inf={p.phi_inf:.4g}")
    z_s    = np.linspace(0.0, p_single.L, N_ISO)
    phi_s  = sol_s.sol(z_s)[0]
    dphi_s = sol_s.sol(z_s)[1]
    F_single = excess_free_energy(phi_s, dphi_s, z_s, p_single)

    # Sweep d from large to small.
    order  = np.argsort(d_grid)[::-1]
    U      = np.full_like(d_grid, np.nan, dtype=float)
    prev_sol = None
    for idx in order:
        dd = float(d_grid[idx])
        p_j = replace(p, L=dd, R=np.inf,
                      parallel_plates=True, conserve_volume=False)
        sol_cold, F_cold = _solve_plate_lowest_F(p_j)
        if prev_sol is not None:
            sol_warm, F_warm = _warm_start_plate(p_j, prev_sol)
        else:
            sol_warm, F_warm = None, float('inf')
        if F_cold <= F_warm:
            sol, F = sol_cold, F_cold
        else:
            sol, F = sol_warm, F_warm
        if sol is None:
            continue
        prev_sol = sol
        U[idx]   = np.pi * R * (F - 2.0 * F_single)
    return U


# ── Additive U_add(d) via 2D axisymmetric integration ───────────────────
def additive_U(d_grid: np.ndarray, p: Parameters,
               psi_interp, dpsi_interp, s_max: float) -> np.ndarray:
    """
    Cross-term-only additive pair potential:

        U_add(d) = F^(2)[phi_inf + psi(s1) + psi(s2)] - 2 F^(1)[psi]

    Self-pieces cancel; we evaluate the cross volume integrand
        Delta_df + kappa * dpsi(s1) dpsi(s2) (n1.n2)
    on a cylindrical (rho, z) grid, plus the cross surface term
        2 * omega_b * 2pi R^2 * int_0^pi psi(s2(theta)) sin theta dtheta.
    """
    R       = p.R
    phi_inf = p.phi_inf
    mu_ref  = p.nb * dfb_dphi(np.asarray([phi_inf]), p)[0]
    f_inf   = fb(np.asarray([phi_inf]), p)[0]

    def df_of(phi_arr):
        return (fb(phi_arr, p) - f_inf - mu_ref / p.nb * (phi_arr - phi_inf))

    d_max    = float(np.max(d_grid))
    z_grid   = np.linspace(-s_max - R, R + d_max + s_max + R, N_Z)
    rho_grid = np.linspace(0.0, s_max + R, N_RHO)
    Rho, Z   = np.meshgrid(rho_grid, z_grid, indexing='ij')

    theta  = np.linspace(0.0, np.pi, 200)
    sin_t  = np.sin(theta)

    U = np.zeros_like(d_grid, dtype=float)
    for i, dd in enumerate(d_grid):
        c2_z = 2.0 * R + float(dd)         # sphere1 at z=0, sphere2 at z = 2R+d

        r1 = np.sqrt(Rho ** 2 + Z ** 2)
        r2 = np.sqrt(Rho ** 2 + (Z - c2_z) ** 2)
        s1 = r1 - R
        s2 = r2 - R

        outside = (s1 >= 0.0) & (s2 >= 0.0)

        psi1 = _psi_eval(psi_interp, s1)
        psi2 = _psi_eval(psi_interp, s2)
        dps1 = _psi_eval(dpsi_interp, s1)
        dps2 = _psi_eval(dpsi_interp, s2)

        with np.errstate(invalid='ignore', divide='ignore'):
            n_dot = (Rho ** 2 + Z * (Z - c2_z)) / np.maximum(r1 * r2, 1e-30)

        phi_full = phi_inf + psi1 + psi2
        phi_a    = phi_inf + psi1
        phi_b    = phi_inf + psi2
        delta_df = df_of(phi_full) - df_of(phi_a) - df_of(phi_b)
        cross_grad = p.kappa * dps1 * dps2 * n_dot
        integrand  = (delta_df + cross_grad) * (2.0 * np.pi * Rho)
        integrand  = np.where(outside, integrand, 0.0)

        vol_cross = np.trapezoid(np.trapezoid(integrand, z_grid, axis=1),
                                 rho_grid, axis=0)

        rho_s = R * sin_t
        z_s   = R * np.cos(theta)
        s2_on_s1 = np.sqrt(rho_s ** 2 + (z_s - c2_z) ** 2) - R
        psi2_on_s1 = _psi_eval(psi_interp, s2_on_s1)
        surf_cross = 2.0 * p.omega_b * 2.0 * np.pi * R ** 2 * \
                     np.trapezoid(psi2_on_s1 * sin_t, theta)

        U[i] = vol_cross + surf_cross
    return U


# ── Total U(d) with smooth crossover ────────────────────────────────────
def total_U(d_grid: np.ndarray, p: Parameters):
    psi_i, dpsi_i, s_max = solve_isolated_sphere(p)
    U_D = derjaguin_U(d_grid, p)
    U_A = additive_U(d_grid, p, psi_i, dpsi_i, s_max)

    # Where the Derjaguin solve failed (NaN), fall back to additive.  This
    # keeps U(d) well-defined; the worst case is under-attraction in the
    # small-d wet-branch regime if a few plate solves failed.
    U_D_safe = np.where(np.isnan(U_D), U_A, U_D)

    w = expit((CROSS_CENTER * p.R - d_grid) / (CROSS_WIDTH * p.R))
    U = w * U_D_safe + (1.0 - w) * U_A
    return U, U_D, U_A, psi_i


# ── Mean probe concentration as a function of total phi ─────────────────
def compute_P_eq(p: Parameters, d_grid: np.ndarray):
    """Linear-response branch: U_lin = 4 pi R^2 omega_b psi_lin(d)."""
    psi_lin = psi_linear(d_grid, p)[0]
    U_lin   = 4.0 * np.pi * p.R ** 2 * p.omega_b * psi_lin

    P_eq_unnorm = (d_grid + p.R) ** 2 * np.exp(
        -3.0 * (d_grid + p.R) ** 2 / (2.0 * s * a ** 2) - U_lin
    )
    P_eq = P_eq_unnorm / np.sum(P_eq_unnorm)

    Phi_d = p.phi_inf + psi_lin[0] + psi_lin
    return P_eq, Phi_d, U_lin


def compute_P_eq_full(p: Parameters, d_grid: np.ndarray):
    """Nonlinear branch: U(d) blends Derjaguin (small d) and additive (large d)."""
    U, U_D, U_A, psi_interp = total_U(d_grid, p)

    P_eq_unnorm = (d_grid + p.R) ** 2 * np.exp(
        -3.0 * (d_grid + p.R) ** 2 / (2.0 * s * a ** 2) - U
    )
    P_eq = P_eq_unnorm / np.sum(P_eq_unnorm)

    psi_at_surface = float(_psi_eval(psi_interp, 0.0))
    psi_at_d       = _psi_eval(psi_interp, d_grid)
    Phi_d = p.phi_inf + psi_at_surface + psi_at_d
    return P_eq, Phi_d, U, U_D, U_A


# ── phi_probe scan over phi_TF ──────────────────────────────────────────
def fisher_rate(p2: Probe, phi_probe: np.ndarray, phi_tf: np.ndarray):
    fisher_base = phi_tf * p2.kb * p2.ku / (p2.kb + p2.ku)
    dphi_probe  = np.gradient(phi_probe, phi_tf)
    inv_taub    = phi_probe * p2.kb * p2.ku / (p2.kb + p2.ku)
    fisher_probe = phi_tf ** 2 * inv_taub * (dphi_probe / phi_probe) ** 2
    return fisher_base, fisher_probe


if __name__ == "__main__":
    import datetime

    n_params    = 20
    phi_tf_vals = np.linspace(0.001, 0.01, n_params)
    d_grid      = np.linspace(0.02 * p.R, L_ISO, N_DGRID)

    phi_probe_lin  = np.zeros(n_params)
    phi_probe_full = np.zeros(n_params)

    for i, _phi in enumerate(phi_tf_vals):
        p.phi_inf = float(_phi)
        P_eq_l, Phi_d_l, _ = compute_P_eq(p, d_grid)
        phi_probe_lin[i] = float(np.sum(P_eq_l * Phi_d_l))
        try:
            P_eq_f, Phi_d_f, U_f, _, _ = compute_P_eq_full(p, d_grid)
            phi_probe_full[i] = float(np.sum(P_eq_f * Phi_d_f))
        except Exception as e:
            print(f"  phi_inf={_phi:.4f}: nonlinear failed ({e}); falling back to linear")
            phi_probe_full[i] = phi_probe_lin[i]
        print(f"  phi_inf={_phi:.4f}: phi_probe_lin={phi_probe_lin[i]:.4f}, "
              f"phi_probe_full={phi_probe_full[i]:.4f}")

    probe = Probe(ku=1.0, kb=1.0)
    fisher_base, fisher_lin  = fisher_rate(probe, phi_probe_lin,  phi_tf_vals)
    _,           fisher_full = fisher_rate(probe, phi_probe_full, phi_tf_vals)

    fig, axes = plt.subplots(3, 1, figsize=(5.2, 9.0))

    p.phi_inf = 0.01
    _, _, U_lin_ref = compute_P_eq(p, d_grid)
    _, _, U_full_ref, U_D_ref, U_A_ref = compute_P_eq_full(p, d_grid)

    ax = axes[0]
    ax.plot(d_grid, U_lin_ref,  color="black",   lw=1.5, ls="--", label="linear")
    ax.plot(d_grid, U_D_ref,    color="#9333EA", lw=1.0, ls=":",  label="Derjaguin")
    ax.plot(d_grid, U_A_ref,    color="#10B981", lw=1.0, ls=":",  label="additive")
    ax.plot(d_grid, U_full_ref, color="#2563EB", lw=2.0,          label="blended U(d)")
    ax.axvline(CROSS_CENTER * p.R, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("d (surface-to-surface distance)")
    ax.set_ylabel("U(d)")
    ax.set_title(f"Pair potential (phi_inf={p.phi_inf})")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    # ax.plot(phi_tf_vals, phi_probe_lin,  color="C2", lw=1.5, ls="--", label="linear")
    # ax.plot(phi_tf_vals, phi_probe_full, color="C2", lw=2.0,           label="nonlinear")
    # ax.plot(phi_tf_vals, phi_tf_vals,    color="gray", lw=1.0, ls=":",  label="no interaction")
    # ax.set_xlabel(r"$\phi_{TF}$ (total volume fraction)")
    # ax.set_ylabel(r"$\phi_{probe}$")
    # ax.legend(fontsize=9)

    # ax = axes[2]
    # ax.plot(phi_tf_vals, fisher_base, color="red", lw=1.0, ls="--", label="control")
    # ax.plot(phi_tf_vals, fisher_lin,  color="red", lw=1.5, ls=":",  label="linear")
    # ax.plot(phi_tf_vals, fisher_full, color="red", lw=2.0,           label="nonlinear")
    # ax.set_xlabel(r"$\phi_{TF}$ (total volume fraction)")
    # ax.set_ylabel("Fisher information rate")
    # ax.set_yscale("log")
    # ax.legend(fontsize=9)

    fig.tight_layout()
    today = datetime.date.today().strftime("%Y%m%d")
    out = f"media/figures/{today}/fisher_rate_from_profile.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved figure → {out}")