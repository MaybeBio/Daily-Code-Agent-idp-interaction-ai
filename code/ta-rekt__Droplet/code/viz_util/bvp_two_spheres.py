"""
2D axisymmetric BVP solver for two spheres in a Flory-Huggins bulk.

Two identical spheres of radius R sit on the z-axis at z = 0 and z = 2R + d,
so the surface-to-surface separation is d.  The equilibrium volume fraction
phi(rho, z) satisfies

    kappa * Laplacian(phi) = dfb/dphi - mu_inf / nb

on the cylindrical (rho, z) domain with axial symmetry, where the Laplacian
is the cylindrical (rho, z) one:

    Laplacian phi = d2 phi/drho2 + (1/rho) d phi/drho + d2 phi/dz2

Boundary conditions:
  - Robin (adsorption) on each sphere:  -kappa * d phi/dn + omega_b = 0
    with n the outward normal pointing away from the sphere centre.
  - Dirichlet phi = phi_inf on the outer cylinder rho = rho_max and the z
    end-caps z = -z_buf, z = 2R + d + z_buf.  This is the "absorbing"
    truncation of the far-field reservoir.
  - dphi/drho = 0 on the symmetry axis rho = 0.

Sphere interiors are masked out of the FD grid; their surfaces are handled
with a quadratic ghost-cell immersed-boundary method.  Each ghost cell G
is extrapolated from TWO mirror-image probe points P1, P2 along the
outward normal: a quadratic phi(s) = a + b*s + c*s^2 is fit to phi(P1),
phi(P2), and the Robin BC dphi/ds|_{s=delta} = omega_b/kappa exactly at
the sphere surface.  This gives O(dx^2) global accuracy regardless of
whether either probe's bilinear stencil had to be "pushed out" past
inside-sphere cells, which the linear (single-image) variant did not
handle correctly.

Nonlinear solve: scipy.optimize.newton_krylov, warm-started from the linear-
response superposition of two single-sphere psi_linear profiles.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
from scipy.optimize import newton_krylov, NoConvergence

from viz_util.bvp_solver import (
    Parameters, dfb_dphi, d2fb_dphi2, fb, compute_mu_inf, clip01, psi_linear,
)


# ------------------------------------------------------------------
# Geometry and node classification
# ------------------------------------------------------------------

@dataclass
class Geometry:
    p: Parameters
    d: float
    Nrho: int
    Nz: int
    rho_grid: np.ndarray
    z_grid: np.ndarray
    drho: float
    dz: float
    z_c1: float
    z_c2: float
    inside_mask: np.ndarray
    fluid_mask: np.ndarray
    dir_mask: np.ndarray
    unk_mask: np.ndarray
    # ghost-cell data: one row per ghost cell.
    #   ghost_ij[k]            : (i, j) grid index of the ghost cell G_k
    #   ghost_i0j0_1[k]        : base corner of P1's bilinear stencil
    #   ghost_wts_1[k]         : (4,) bilinear weights at P1
    #   ghost_i0j0_2[k]        : base corner of P2's bilinear stencil
    #   ghost_wts_2[k]         : (4,) bilinear weights at P2
    #   ghost_coef_phi1[k]     : alpha1 -- coefficient of phi(P1) in phi(G)
    #   ghost_coef_phi2[k]     : alpha2 -- coefficient of phi(P2) in phi(G)
    #   ghost_coef_bc[k]       : beta   -- multiplies omega_b/kappa in phi(G)
    #   ghost_used_linear[k]   : True iff quadratic was ill-conditioned and
    #                            the ghost falls back to the linear formula
    #                            phi(G) = phi(P1) - d1 * omega_b/kappa
    #                            (i.e. alpha1=1, alpha2=0, beta=-d1)
    #   ghost_delta[k]         : R - |G - nearest sphere centre|  (diagnostic)
    #   ghost_d1[k], ghost_d2[k]: distances G->P1, G->P2          (diagnostic)
    ghost_ij: np.ndarray
    ghost_i0j0_1: np.ndarray
    ghost_wts_1: np.ndarray
    ghost_i0j0_2: np.ndarray
    ghost_wts_2: np.ndarray
    ghost_coef_phi1: np.ndarray
    ghost_coef_phi2: np.ndarray
    ghost_coef_bc: np.ndarray
    ghost_used_linear: np.ndarray
    ghost_delta: np.ndarray
    ghost_d1: np.ndarray
    ghost_d2: np.ndarray


def _decay_length(p: Parameters) -> float:
    d2 = float(d2fb_dphi2(np.asarray([p.phi_inf]), p)[0])
    if d2 > 0.05:
        return float(np.sqrt(p.kappa / d2))
    return float(np.sqrt(p.kappa * p.nb / 0.05))


def _find_clean_stencil(rho_G: float, z_G: float,
                        u_rho: float, u_z: float,
                        d_initial: float, step: float,
                        rho_grid: np.ndarray, z_grid: np.ndarray,
                        drho: float, dz: float,
                        Nrho: int, Nz: int,
                        inside_mask: np.ndarray,
                        max_iter: int = 80):
    """
    Walk a probe point outward along (u_rho, u_z) until its 4-cell bilinear
    stencil contains no inside-sphere cells.  Returns the final probe
    distance and stencil:

        (d_used, i0, j0, wts, clean_bool)

    where wts has shape (4,) ordered (i0,j0), (i0+1,j0), (i0,j0+1),
    (i0+1,j0+1).  `clean_bool` is False if the iteration cap was reached
    without finding a clean stencil (the last attempt is still returned).
    """
    d_use = d_initial
    chosen_i0 = chosen_j0 = 0
    chosen_w = np.zeros(4)
    clean = False
    for _ in range(max_iter):
        rho_I = rho_G + d_use * u_rho
        z_I   = z_G   + d_use * u_z
        rho_I_c = float(np.clip(rho_I, rho_grid[0], rho_grid[-1] - 1e-12 * drho))
        z_I_c   = float(np.clip(z_I,   z_grid[0],   z_grid[-1]   - 1e-12 * dz))
        ti = (rho_I_c - rho_grid[0]) / drho
        tj = (z_I_c   - z_grid[0])   / dz
        i0 = int(np.clip(int(np.floor(ti)), 0, Nrho - 2))
        j0 = int(np.clip(int(np.floor(tj)), 0, Nz   - 2))
        t_rho = ti - i0
        t_z   = tj - j0
        stencil_clean = not (inside_mask[i0, j0] or inside_mask[i0 + 1, j0]
                             or inside_mask[i0, j0 + 1]
                             or inside_mask[i0 + 1, j0 + 1])
        chosen_i0, chosen_j0 = i0, j0
        chosen_w = np.array([
            (1.0 - t_rho) * (1.0 - t_z),
            t_rho         * (1.0 - t_z),
            (1.0 - t_rho) * t_z,
            t_rho         * t_z,
        ])
        if stencil_clean:
            clean = True
            break
        d_use += step
    return d_use, chosen_i0, chosen_j0, chosen_w, clean


def _quadratic_ghost_coefficients(delta: float, d1: float, d2: float):
    """
    Closed-form coefficients (alpha1, alpha2, beta) such that

        phi(G) = alpha1 * phi(P1) + alpha2 * phi(P2) + beta * (omega_b/kappa)

    is the value at G implied by the quadratic phi(s) = a + b*s + c*s^2 that
    interpolates phi(P1) at s=d1, phi(P2) at s=d2, and satisfies the Robin
    BC dphi/ds|_{s=delta} = omega_b/kappa.  s is signed distance from G
    along the outward normal; the sphere surface is at s=delta.
    """
    K = d1 + d2 - 2.0 * delta
    inv_K = 1.0 / K
    inv_dd = 1.0 / (d1 - d2)
    c_phi1 = inv_dd * inv_K
    c_phi2 = -inv_dd * inv_K
    c_bc   = -inv_K
    b_phi1 = -2.0 * delta * c_phi1
    b_phi2 = -2.0 * delta * c_phi2
    b_bc   = 1.0 - 2.0 * delta * c_bc
    alpha1 = 1.0 - d1 * b_phi1 - d1 * d1 * c_phi1
    alpha2 = -d1 * b_phi2 - d1 * d1 * c_phi2
    beta   = -d1 * b_bc   - d1 * d1 * c_bc
    return alpha1, alpha2, beta


def build_geometry(p: Parameters, d: float,
                   rho_max: Optional[float] = None,
                   z_buf:   Optional[float] = None,
                   Nrho: int = 80, Nz: int = 200) -> Geometry:
    R = p.R
    ell = _decay_length(p)
    if rho_max is None:
        rho_max = R + 4.0 * ell
    if z_buf is None:
        z_buf = 4.0 * ell

    z_c1 = 0.0
    z_c2 = 2.0 * R + d

    rho_grid = np.linspace(0.0, rho_max, Nrho)
    # z-box bounds include a buffer beyond each sphere's pole (not just past
    # each centre), so the spheres are fully enclosed with z_buf of clear
    # reservoir on either side.
    z_grid   = np.linspace(z_c1 - R - z_buf, z_c2 + R + z_buf, Nz)
    drho = float(rho_grid[1] - rho_grid[0])
    dz   = float(z_grid[1]   - z_grid[0])

    Rho, Z = np.meshgrid(rho_grid, z_grid, indexing="ij")
    r1_sq = Rho ** 2 + (Z - z_c1) ** 2
    r2_sq = Rho ** 2 + (Z - z_c2) ** 2
    inside_mask = (r1_sq < R ** 2) | (r2_sq < R ** 2)
    fluid_mask  = ~inside_mask

    dir_mask = np.zeros_like(fluid_mask)
    dir_mask[-1, :] = True
    dir_mask[:,  0] = True
    dir_mask[:, -1] = True
    dir_mask &= fluid_mask
    unk_mask = fluid_mask & ~dir_mask

    # Ghost cells: inside-sphere nodes with at least one 4-connected fluid
    # neighbour.  For each ghost G we place TWO probe points P1, P2 along
    # the outward normal:
    #   - P1: natural mirror at d1_init = 2*delta (delta = R - |G-centre|).
    #         Push out along the normal until its bilinear stencil contains
    #         no inside-sphere cells.
    #   - P2: starts at d1_used + step and is pushed out independently.
    # A quadratic phi(s) = a + b*s + c*s^2 is fit to phi(P1), phi(P2) and the
    # Robin BC at s = delta; the resulting linear functional phi(G) =
    # alpha1*phi(P1) + alpha2*phi(P2) + beta*(omega_b/kappa) is precomputed
    # via _quadratic_ghost_coefficients.  Fallback to the linear formula
    # phi(G) = phi(P1) - d1 * omega_b/kappa when the quadratic is ill-
    # conditioned (e.g. the two probes alias the same box-clamped corner).
    step = max(drho, dz)
    max_iter = 80
    sep_min = 0.9 * step
    K_min = 1e-9 * step
    ghost_records = []
    n_fallback = 0
    nbr_offsets = ((1, 0), (-1, 0), (0, 1), (0, -1))
    ii_inside, jj_inside = np.where(inside_mask)
    for i, j in zip(ii_inside, jj_inside):
        has_fluid = False
        for di, dj in nbr_offsets:
            ii, jj = i + di, j + dj
            if 0 <= ii < Nrho and 0 <= jj < Nz and fluid_mask[ii, jj]:
                has_fluid = True
                break
        if not has_fluid:
            continue
        rho_G = rho_grid[i]
        z_G   = z_grid[j]
        dist1 = np.hypot(rho_G, z_G - z_c1)
        dist2 = np.hypot(rho_G, z_G - z_c2)
        if dist1 <= dist2:
            zc, Gc = z_c1, dist1
        else:
            zc, Gc = z_c2, dist2
        if Gc < 1e-12:
            continue
        u_rho = rho_G / Gc
        u_z   = (z_G - zc) / Gc
        delta = R - Gc

        d1, i0_1, j0_1, w1, clean1 = _find_clean_stencil(
            rho_G, z_G, u_rho, u_z, 2.0 * delta, step,
            rho_grid, z_grid, drho, dz, Nrho, Nz, inside_mask, max_iter)
        d2, i0_2, j0_2, w2, clean2 = _find_clean_stencil(
            rho_G, z_G, u_rho, u_z, d1 + step, step,
            rho_grid, z_grid, drho, dz, Nrho, Nz, inside_mask, max_iter)

        # Linear-fallback triggers: either push-out failed, or the two
        # probes collapsed onto each other (e.g. both got clamped against
        # the same outer face), or the quadratic denominator d1+d2-2*delta
        # would underflow.
        use_linear = (not clean1) or (not clean2) \
                     or ((d2 - d1) < sep_min) \
                     or (abs(d1 + d2 - 2.0 * delta) < K_min) \
                     or (abs(d1 - d2) < K_min)
        if use_linear:
            alpha1, alpha2, beta = 1.0, 0.0, -d1
            i0_2, j0_2 = i0_1, j0_1
            w2 = np.zeros(4)
            n_fallback += 1
        else:
            alpha1, alpha2, beta = _quadratic_ghost_coefficients(delta, d1, d2)

        ghost_records.append((i, j,
                              i0_1, j0_1, w1,
                              i0_2, j0_2, w2,
                              alpha1, alpha2, beta,
                              bool(use_linear),
                              delta, d1, d2))

    if ghost_records:
        ghost_ij       = np.array([(r[0], r[1]) for r in ghost_records], dtype=int)
        ghost_i0j0_1   = np.array([(r[2], r[3]) for r in ghost_records], dtype=int)
        ghost_wts_1    = np.array([r[4]         for r in ghost_records], dtype=float)
        ghost_i0j0_2   = np.array([(r[5], r[6]) for r in ghost_records], dtype=int)
        ghost_wts_2    = np.array([r[7]         for r in ghost_records], dtype=float)
        ghost_coef_phi1 = np.array([r[8]  for r in ghost_records], dtype=float)
        ghost_coef_phi2 = np.array([r[9]  for r in ghost_records], dtype=float)
        ghost_coef_bc   = np.array([r[10] for r in ghost_records], dtype=float)
        ghost_used_linear = np.array([r[11] for r in ghost_records], dtype=bool)
        ghost_delta = np.array([r[12] for r in ghost_records], dtype=float)
        ghost_d1    = np.array([r[13] for r in ghost_records], dtype=float)
        ghost_d2    = np.array([r[14] for r in ghost_records], dtype=float)
    else:
        ghost_ij          = np.zeros((0, 2), dtype=int)
        ghost_i0j0_1      = np.zeros((0, 2), dtype=int)
        ghost_wts_1       = np.zeros((0, 4), dtype=float)
        ghost_i0j0_2      = np.zeros((0, 2), dtype=int)
        ghost_wts_2       = np.zeros((0, 4), dtype=float)
        ghost_coef_phi1   = np.zeros((0,), dtype=float)
        ghost_coef_phi2   = np.zeros((0,), dtype=float)
        ghost_coef_bc     = np.zeros((0,), dtype=float)
        ghost_used_linear = np.zeros((0,), dtype=bool)
        ghost_delta       = np.zeros((0,), dtype=float)
        ghost_d1          = np.zeros((0,), dtype=float)
        ghost_d2          = np.zeros((0,), dtype=float)

    return Geometry(
        p=p, d=d, Nrho=Nrho, Nz=Nz,
        rho_grid=rho_grid, z_grid=z_grid,
        drho=drho, dz=dz, z_c1=z_c1, z_c2=z_c2,
        inside_mask=inside_mask, fluid_mask=fluid_mask,
        dir_mask=dir_mask, unk_mask=unk_mask,
        ghost_ij=ghost_ij,
        ghost_i0j0_1=ghost_i0j0_1, ghost_wts_1=ghost_wts_1,
        ghost_i0j0_2=ghost_i0j0_2, ghost_wts_2=ghost_wts_2,
        ghost_coef_phi1=ghost_coef_phi1,
        ghost_coef_phi2=ghost_coef_phi2,
        ghost_coef_bc=ghost_coef_bc,
        ghost_used_linear=ghost_used_linear,
        ghost_delta=ghost_delta,
        ghost_d1=ghost_d1, ghost_d2=ghost_d2,
    )


# ------------------------------------------------------------------
# Scatter (assemble the full phi field including Dirichlet + ghost fills)
# ------------------------------------------------------------------

def _scatter(x_unk: np.ndarray, geom: Geometry) -> np.ndarray:
    p = geom.p
    phi = np.full((geom.Nrho, geom.Nz), p.phi_inf)
    phi[geom.unk_mask] = x_unk
    if geom.ghost_ij.shape[0] > 0:
        gi = geom.ghost_ij[:, 0]
        gj = geom.ghost_ij[:, 1]
        i0_1 = geom.ghost_i0j0_1[:, 0]
        j0_1 = geom.ghost_i0j0_1[:, 1]
        i0_2 = geom.ghost_i0j0_2[:, 0]
        j0_2 = geom.ghost_i0j0_2[:, 1]
        # By construction the bilinear stencils for P1 and P2 either lie in
        # the fluid / on the Dirichlet boundary, or the ghost has fallen
        # back to the linear formula (in which case alpha2=0 and ghost_wts_2
        # is all zeros, so the P2 sample is harmless).  A single pass
        # suffices.
        phi_P1 = (geom.ghost_wts_1[:, 0] * phi[i0_1,     j0_1]
                + geom.ghost_wts_1[:, 1] * phi[i0_1 + 1, j0_1]
                + geom.ghost_wts_1[:, 2] * phi[i0_1,     j0_1 + 1]
                + geom.ghost_wts_1[:, 3] * phi[i0_1 + 1, j0_1 + 1])
        phi_P2 = (geom.ghost_wts_2[:, 0] * phi[i0_2,     j0_2]
                + geom.ghost_wts_2[:, 1] * phi[i0_2 + 1, j0_2]
                + geom.ghost_wts_2[:, 2] * phi[i0_2,     j0_2 + 1]
                + geom.ghost_wts_2[:, 3] * phi[i0_2 + 1, j0_2 + 1])
        bc = p.omega_b / p.kappa
        phi[gi, gj] = (geom.ghost_coef_phi1 * phi_P1
                     + geom.ghost_coef_phi2 * phi_P2
                     + geom.ghost_coef_bc   * bc)
    return phi


# ------------------------------------------------------------------
# Cylindrical Laplacian on the structured grid
# ------------------------------------------------------------------

def _laplacian(phi: np.ndarray, geom: Geometry) -> np.ndarray:
    drho, dz = geom.drho, geom.dz
    lap = np.zeros_like(phi)
    rho_int = geom.rho_grid[1:-1, None]
    lap[1:-1, 1:-1] = (
        (phi[2:, 1:-1] - 2.0 * phi[1:-1, 1:-1] + phi[:-2, 1:-1]) / drho ** 2
        + (phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2.0 * drho * rho_int)
        + (phi[1:-1, 2:] - 2.0 * phi[1:-1, 1:-1] + phi[1:-1, :-2]) / dz ** 2
    )
    # Axis: dphi/drho = 0 by symmetry; L'Hopital gives (1/rho)dphi/drho ->
    # d2phi/drho2, so the axial Laplacian is 2 d2phi/drho2 + d2phi/dz2.
    # With the symmetric ghost phi[-1, j] = phi[1, j], 2 d2phi/drho2 reduces
    # to 4 (phi[1, j] - phi[0, j]) / drho2.
    lap[0, 1:-1] = (
        4.0 * (phi[1, 1:-1] - phi[0, 1:-1]) / drho ** 2
        + (phi[0, 2:] - 2.0 * phi[0, 1:-1] + phi[0, :-2]) / dz ** 2
    )
    return lap


# ------------------------------------------------------------------
# Residual + solve
# ------------------------------------------------------------------

def _make_residual(geom: Geometry, mu_inf: float):
    p = geom.p
    unk = geom.unk_mask

    def residual(x_unk: np.ndarray) -> np.ndarray:
        phi = _scatter(x_unk, geom)
        lap = _laplacian(phi, geom)
        driving = dfb_dphi(phi, p) - mu_inf / p.nb
        return (p.kappa * lap - driving)[unk]

    return residual


def _initial_guess(geom: Geometry) -> np.ndarray:
    p = geom.p
    Rho, Z = np.meshgrid(geom.rho_grid, geom.z_grid, indexing="ij")
    s1 = np.sqrt(Rho ** 2 + (Z - geom.z_c1) ** 2) - p.R
    s2 = np.sqrt(Rho ** 2 + (Z - geom.z_c2) ** 2) - p.R
    s1c = np.clip(s1, 0.0, None)
    s2c = np.clip(s2, 0.0, None)
    psi1, _, _ = psi_linear(s1c, p)
    psi2, _, _ = psi_linear(s2c, p)
    phi0 = p.phi_inf + psi1 + psi2
    phi0 = clip01(phi0, p.eps)
    return phi0[geom.unk_mask]


def solve_two_spheres(p: Parameters, d: float,
                      rho_max: Optional[float] = None,
                      z_buf:   Optional[float] = None,
                      Nrho: int = 80, Nz: int = 200,
                      f_tol: float = 1e-6,
                      maxiter: int = 200,
                      verbose: bool = False,
                      omega_continuation: bool = True) -> dict:
    """
    Solve the two-sphere axisymmetric BVP.

    Returns dict with keys:
      'phi'      : (Nrho, Nz) array including ghost-filled inside cells
                   (and phi_inf in non-ghost inside cells)
      'geom'     : Geometry struct (grids, masks, etc.)
      'p', 'd'   : echo back inputs
      'mu_inf'   : reservoir chemical potential
      'success'  : True if newton_krylov converged
      'x'        : flat solution vector (unknowns only)
    """
    geom   = build_geometry(p, d, rho_max=rho_max, z_buf=z_buf,
                            Nrho=Nrho, Nz=Nz)
    mu_inf = compute_mu_inf(p)

    x0 = _initial_guess(geom)
    success = True
    try:
        x_sol = newton_krylov(_make_residual(geom, mu_inf), x0,
                              f_tol=f_tol, maxiter=maxiter, verbose=verbose)
    except NoConvergence as e:
        x_sol = np.asarray(e.args[0])
        success = False

    # If the cold start failed and continuation is on, ramp omega_b from a
    # small value up to the target, reusing the previous solution each step.
    if not success and omega_continuation and p.omega_b != 0.0:
        if verbose:
            print("Cold start failed; trying omega_b continuation ...")
        n_steps = 6
        x_cur = _initial_guess(build_geometry(replace(p, omega_b=0.0), d,
                                              rho_max=rho_max, z_buf=z_buf,
                                              Nrho=Nrho, Nz=Nz))
        for k in range(1, n_steps + 1):
            frac = k / n_steps
            p_k = replace(p, omega_b=frac * p.omega_b)
            geom_k = build_geometry(p_k, d, rho_max=rho_max, z_buf=z_buf,
                                    Nrho=Nrho, Nz=Nz)
            try:
                x_cur = newton_krylov(
                    _make_residual(geom_k, compute_mu_inf(p_k)),
                    x_cur, f_tol=f_tol, maxiter=maxiter, verbose=False,
                )
            except NoConvergence as e:
                x_cur = np.asarray(e.args[0])
                if verbose:
                    print(f"  continuation step {k}/{n_steps} (frac={frac:.2f}) "
                          f"did not converge")
                break
            if verbose:
                print(f"  continuation step {k}/{n_steps} (frac={frac:.2f}) ok")
        x_sol = x_cur
        geom = geom_k
        mu_inf = compute_mu_inf(p)
        # Final residual check
        res = _make_residual(geom, mu_inf)(x_sol)
        success = float(np.max(np.abs(res))) < 10.0 * f_tol

    phi = _scatter(x_sol, geom)
    return dict(phi=phi, geom=geom, p=p, d=d, mu_inf=mu_inf,
                success=success, x=x_sol)


# ------------------------------------------------------------------
# Bilinear sampling utility (used by the surface integral)
# ------------------------------------------------------------------

def _interp_phi(phi: np.ndarray, geom: Geometry,
                rho_q: np.ndarray, z_q: np.ndarray) -> np.ndarray:
    rho_g, z_g = geom.rho_grid, geom.z_grid
    drho, dz = geom.drho, geom.dz
    rho_c = np.clip(rho_q, rho_g[0], rho_g[-1] - 1e-12 * drho)
    z_c   = np.clip(z_q,   z_g[0],   z_g[-1]   - 1e-12 * dz)
    ti = (rho_c - rho_g[0]) / drho
    tj = (z_c   - z_g[0])   / dz
    i0 = np.clip(np.floor(ti).astype(int), 0, geom.Nrho - 2)
    j0 = np.clip(np.floor(tj).astype(int), 0, geom.Nz   - 2)
    t_rho = ti - i0
    t_z   = tj - j0
    p00 = phi[i0,     j0]
    p10 = phi[i0 + 1, j0]
    p01 = phi[i0,     j0 + 1]
    p11 = phi[i0 + 1, j0 + 1]
    return ((1.0 - t_rho) * (1.0 - t_z) * p00
          + t_rho         * (1.0 - t_z) * p10
          + (1.0 - t_rho) * t_z         * p01
          + t_rho         * t_z         * p11)


# ------------------------------------------------------------------
# Excess free energy
# ------------------------------------------------------------------

def _cell_fluid_fraction(geom: Geometry, nsub: int = 8) -> np.ndarray:
    """
    Per-cell fluid-area fraction in the (rho, z) cross-section.

    Each grid point (rho_i, z_j) owns a drho x dz cell centred on it.  We
    sub-sample the cell with nsub x nsub points and count the fraction
    whose centres lie outside both spheres.  This restores second-order
    spatial accuracy to the volume integral despite the sphere boundary not
    aligning with grid lines: cells fully outside both spheres get weight 1,
    cells fully inside get weight 0, cut cells get a smooth weight in (0, 1).
    Sub-cells with rho < 0 are clipped to rho = 0 (no contribution from
    inside the axis).
    """
    drho, dz = geom.drho, geom.dz
    R2 = geom.p.R ** 2
    offsets = (np.arange(nsub) + 0.5) / nsub - 0.5  # in (-0.5, 0.5)
    counts = np.zeros((geom.Nrho, geom.Nz), dtype=np.int32)
    for a in offsets:
        rho_s = np.maximum(geom.rho_grid + a * drho, 0.0)
        for b in offsets:
            z_s = geom.z_grid + b * dz
            r1sq = rho_s[:, None] ** 2 + (z_s - geom.z_c1)[None, :] ** 2
            r2sq = rho_s[:, None] ** 2 + (z_s - geom.z_c2)[None, :] ** 2
            counts += ((r1sq < R2) | (r2sq < R2)).astype(np.int32)
    return 1.0 - counts / (nsub * nsub)


def _phi_on_sphere(phi: np.ndarray, geom: Geometry,
                   zc: float, theta: np.ndarray,
                   probe_steps: int = 80) -> np.ndarray:
    """
    Sample phi at points on a sphere surface (centre = (0, zc), radius R)
    by quadratic extrapolation from two FLUID probe points along the
    outward normal -- the same closed-form quadratic-fit-with-Robin-BC
    used by the ghost cells, but applied from the fluid side with
    surface offset delta = 0.

    Why this beats bilinear sampling of phi at the surface point: the
    bilinear stencil at a surface point straddles fluid and ghost cells,
    and adjacent ghosts carry values from independently-fitted quadratics
    -- so the sampled phi_surf field has cell-scale jumps that look like
    O(dx) noise in the surface integral.  Probing from the fluid side
    with bilinear interpolation of cells that are guaranteed CLEAN (none
    of the 4 stencil corners is inside a sphere) gives a genuinely smooth
    O(dx^2) phi_surf(theta).

    `probe_steps` is the max push-out iteration count: each probe is
    pushed further along the outward normal until *its bilinear stencil*
    -- not just the probe point itself -- is clear of inside-sphere
    cells.  Checking the stencil (not the point) is critical: a probe
    just outside the surface can still have a stencil corner *inside*
    the sphere if its 4-cell square spans the boundary diagonally.
    """
    p = geom.p
    R = p.R
    drho, dz = geom.drho, geom.dz
    step = max(drho, dz)
    omega_kappa = p.omega_b / p.kappa
    Nrho, Nz = geom.Nrho, geom.Nz
    rho_g = geom.rho_grid
    z_g   = geom.z_grid

    sin_t = np.sin(theta)
    cos_t = np.cos(theta)
    rho_s = R * sin_t
    z_s   = zc + R * cos_t
    u_rho = sin_t
    u_z   = cos_t

    n = theta.shape[0]
    d1 = np.full(n, step)
    d2 = np.full(n, 2.0 * step)

    def _stencil_dirty(d_arr):
        rho_P = rho_s + d_arr * u_rho
        z_P   = z_s   + d_arr * u_z
        rho_Pc = np.clip(rho_P, rho_g[0], rho_g[-1] - 1e-12 * drho)
        z_Pc   = np.clip(z_P,   z_g[0],   z_g[-1]   - 1e-12 * dz)
        ti = (rho_Pc - rho_g[0]) / drho
        tj = (z_Pc   - z_g[0])   / dz
        i0 = np.clip(np.floor(ti).astype(int), 0, Nrho - 2)
        j0 = np.clip(np.floor(tj).astype(int), 0, Nz   - 2)
        return (geom.inside_mask[i0,     j0]
              | geom.inside_mask[i0 + 1, j0]
              | geom.inside_mask[i0,     j0 + 1]
              | geom.inside_mask[i0 + 1, j0 + 1])

    for _ in range(probe_steps):
        bad1 = _stencil_dirty(d1)
        bad2 = _stencil_dirty(d2)
        if not (bad1.any() or bad2.any()):
            break
        d1 = np.where(bad1, d1 + step, d1)
        d2 = np.where(bad2, d2 + step, d2)
        # Maintain at least one cell of separation between P1 and P2.
        d2 = np.maximum(d2, d1 + step)

    rho_P1 = rho_s + d1 * u_rho
    z_P1   = z_s   + d1 * u_z
    rho_P2 = rho_s + d2 * u_rho
    z_P2   = z_s   + d2 * u_z
    phi_P1 = _interp_phi(phi, geom, rho_P1, z_P1)
    phi_P2 = _interp_phi(phi, geom, rho_P2, z_P2)

    # Quadratic phi(s) = a + b*s + c*s^2, s outward distance from S along
    # the normal.  Robin BC: b = omega_b/kappa at s=0.  Closed form (this
    # is the delta=0 specialisation of _quadratic_ghost_coefficients):
    #   c = ((phi_P1 - phi_P2)/(d1-d2) - omega_kappa) / (d1 + d2)
    #   a = phi_P1 - omega_kappa * d1 - c * d1^2 = phi(S)
    c = ((phi_P1 - phi_P2) / (d1 - d2) - omega_kappa) / (d1 + d2)
    return phi_P1 - omega_kappa * d1 - c * d1 * d1


def excess_free_energy_two_spheres(result: dict,
                                   nsub: int = 8,
                                   n_theta: int = 1001) -> float:
    """
    Grand-potential excess over the bulk reservoir at phi_inf:

      Delta Omega
        = int_V [ df(phi) - (1/2) (phi - phi_inf)
                              * (fb'(phi) - fb'(phi_inf)) ] 2 pi rho drho dz
          + (omega_b / 2) sum_{both spheres} int_{surface} (phi - phi_inf) dA

    where df(phi) = fb(phi) - fb(phi_inf) - fb'(phi_inf) (phi - phi_inf) and
    fb' denotes dfb/dphi.  On equilibrium profiles this is mathematically
    identical to the |grad phi|^2 form

      int_V [ df + (kappa/2)|grad phi|^2 ] dV + omega_b * int_dV (phi-phi_inf) dA

    -- the equivalence comes from integration by parts using
    kappa*Laplacian(phi) = fb'(phi) - mu_inf/nb (the Euler-Lagrange equation)
    together with the Robin BC kappa*dphi/dn = omega_b on both spheres.

    This form is much more numerically stable than the |grad phi|^2 form
    because:

      1. It contains no gradient of phi: the FD-gradient term, whose error
         is largest just outside the sphere surfaces where the gradient is
         largest, is gone.
      2. Both df and (1/2) (phi-phi_inf)(fb'(phi)-fb'(phi_inf)) are quadratic
         in psi := phi - phi_inf at leading order with the SAME coefficient
         (1/2) fb''(phi_inf), so their difference is O(psi^3).  The bulk
         integrand decays as psi^3 from the reservoir instead of psi^2, so
         far-field truncation error is much smaller and -- crucially --
         contributions from the box volume cancel cleanly in F(d) - F(d_ref).
      3. The surface coefficient is omega_b/2.  In the original form the
         attractive surface (omega_b<0) and the repulsive bulk-gradient
         term were comparable in magnitude with opposite signs; here that
         cancellation has been done analytically.

    Volume integration uses a cut-cell-aware midpoint Riemann sum: each
    (rho_i, z_j) cell contributes integrand_ij * 2 pi rho_i * drho * dz *
    fluid_fraction_ij, where fluid_fraction is the (nsub x nsub) sub-sampled
    area of the cell outside both spheres.  This avoids the O(drho) error
    of a hard mask at the sphere surfaces.

    Surface integration samples phi at n_theta polar angles via the
    quadratic-from-fluid extrapolation in `_phi_on_sphere` (which avoids
    reading any ghost-cell values), then integrates by Simpson's rule
    (requires n_theta odd).
    """
    from scipy.integrate import simpson

    geom = result["geom"]
    phi  = result["phi"]
    p    = geom.p
    phi_inf = p.phi_inf

    dfb_inf = float(dfb_dphi(np.asarray([phi_inf]), p)[0])
    f_inf   = float(fb(np.asarray([phi_inf]), p)[0])

    psi    = phi - phi_inf
    df_loc = fb(phi, p) - f_inf - dfb_inf * psi
    dfp    = dfb_dphi(phi, p) - dfb_inf
    integrand = df_loc - 0.5 * psi * dfp
    # Inside-sphere cells: phi is either ghost-filled (linear extrapolation)
    # or phi_inf (deep-inside non-ghost cells).  Either way the value is
    # unphysical for the volume integrand; mask to zero.
    integrand = np.where(geom.fluid_mask, integrand, 0.0)

    frac  = _cell_fluid_fraction(geom, nsub=nsub)
    rho_w = 2.0 * np.pi * geom.rho_grid[:, None]
    bulk  = float(np.sum(integrand * rho_w * frac) * geom.drho * geom.dz)

    # Surface integral.  Each sphere parametrised by polar theta in [0, pi];
    # closed form in azimuthal angle gives a factor 2*pi*R^2.
    if n_theta % 2 == 0:
        n_theta += 1
    theta = np.linspace(0.0, np.pi, n_theta)
    sin_t = np.sin(theta)

    surf_total = 0.0
    for zc in (geom.z_c1, geom.z_c2):
        phi_s = _phi_on_sphere(phi, geom, zc, theta)
        surf_total += float(simpson((phi_s - phi_inf) * sin_t, x=theta))
    surf = 0.5 * p.omega_b * 2.0 * np.pi * p.R ** 2 * surf_total
    return bulk + surf


# ------------------------------------------------------------------
# __main__: verification tests + U(d) plot
# ------------------------------------------------------------------

def _run_main():
    import datetime
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from viz_util.fisher_rate_from_profile_full import (
        solve_isolated_sphere, total_U, L_ISO,
    )

    # Parameters at the spinodal threshold chi_b = 2 with a noticeably
    # attractive surface.  Strong enough to produce a clear pair potential
    # but not so strong (chi_b > 2 with large omega_b) that the cold-start
    # initial guess fails to converge to the wetting-layer branch.
    p = Parameters(
        nb=1.0, chi_b=1.5, omega_b=-0.5, kappa=1.0,
        phi_inf=0.1, L=20.0, R=5.0, eps=1e-6,
        conserve_volume=False, parallel_plates=False,
    )

    # ------------------------------------------------------------------
    # Test 1: far-separation matches isolated single-sphere on axis
    # ------------------------------------------------------------------
    print("\n[1] Far-separation comparison (d = 4 R) ...")
    d_far = 4.0 * p.R
    res = solve_two_spheres(p, d_far, Nrho=200, Nz=600, verbose=False)
    print(f"  converged: {res['success']}")
    z_axis   = res["geom"].z_grid
    phi_axis = res["phi"][0, :]
    psi_i, _, s_max = solve_isolated_sphere(p)
    # Region above sphere 1 (z in (R, R + d_far/2))
    region = (z_axis > p.R + 1e-3) & (z_axis < p.R + d_far / 2.0)
    s_ref  = np.clip(z_axis[region] - p.R, 0.0, s_max)
    phi_ref = p.phi_inf + psi_i(s_ref)
    err = float(np.max(np.abs(phi_axis[region] - phi_ref)))
    print(f"  max |two-sphere - single-sphere| on axis, near sphere 1: {err:.3e}")

    # ------------------------------------------------------------------
    # Test 2: linear-response superposition (small omega_b)
    # ------------------------------------------------------------------
    print("\n[2] Linear-response limit (omega_b = -0.05, chi_b = 0) ...")
    p_lin = replace(p, chi_b=0.0, omega_b=-0.05)
    d_lin = 2.0 * p_lin.R
    res_lin = solve_two_spheres(p_lin, d_lin, Nrho=60, Nz=200, verbose=False)
    print(f"  converged: {res_lin['success']}")
    geom = res_lin["geom"]
    Rho, Z = np.meshgrid(geom.rho_grid, geom.z_grid, indexing="ij")
    s1 = np.clip(np.sqrt(Rho ** 2 + (Z - geom.z_c1) ** 2) - p_lin.R, 0.0, None)
    s2 = np.clip(np.sqrt(Rho ** 2 + (Z - geom.z_c2) ** 2) - p_lin.R, 0.0, None)
    psi1 = psi_linear(s1, p_lin)[0]
    psi2 = psi_linear(s2, p_lin)[0]
    phi_lin_pred = p_lin.phi_inf + psi1 + psi2
    diff = np.abs(res_lin["phi"] - phi_lin_pred)
    diff = np.where(geom.fluid_mask, diff, 0.0)
    err_lin = float(np.max(diff))
    print(f"  max |numeric - linear superposition|: {err_lin:.3e}  "
          f"(expected O(omega_b^2) = {p_lin.omega_b ** 2:.3e})")

    # ------------------------------------------------------------------
    # Test 3: U(d) sweep vs Derjaguin / additive
    # ------------------------------------------------------------------
    # Reference F at a large separation: U(d) = F(d) - F(d_ref).  This
    # cancels the systematic discretisation error of the 2D FD method and
    # gives an interaction energy that vanishes at large d by construction.
    print("\n[3] U(d) sweep vs Derjaguin / additive references ...")
    d_grid = np.linspace(0.3 * p.R, 2.5 * p.R, 8)
    Nrho = 500
    dz_target = 0.01
    def _Nz(dd):
        return int((2 * p.R + dd + 8.0) / dz_target) + 1

    d_ref = 8.0 * p.R
    res_ref = solve_two_spheres(p, d_ref, Nrho=Nrho, Nz=_Nz(d_ref),
                                z_buf=4.0, rho_max=10.0, verbose=False)
    F_ref = excess_free_energy_two_spheres(res_ref)
    print(f"  reference F(d={d_ref}) = {F_ref:.4e}")
    U_full = np.full(d_grid.shape, np.nan)
    for k, dd in enumerate(d_grid):
        res = solve_two_spheres(p, float(dd), Nrho=Nrho, Nz=_Nz(dd),
                                z_buf=4.0, rho_max=10.0, verbose=False)
        if not res["success"]:
            print(f"  d={dd:.3f}: NOT converged, skipping")
            continue
        F = excess_free_energy_two_spheres(res)
        U_full[k] = F - F_ref
        print(f"  d={dd:.3f}: F={F:.4e}  U={U_full[k]:+.4e}")

    print("  Computing Derjaguin + additive reference (may take a moment) ...")
    try:
        U_blend, U_D, U_A, _ = total_U(d_grid, p)
    except Exception as exc:
        print(f"  reference failed ({exc}); plotting only U_full")
        U_blend = U_D = U_A = np.full_like(d_grid, np.nan)

    today = datetime.date.today().strftime("%Y%m%d")
    out_dir = f"media/figures/{today}"
    os.makedirs(out_dir, exist_ok=True)
    out = f"{out_dir}/bvp_two_spheres_U_d.png"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.plot(d_grid, U_full, "o-",  color="C0",       label="full 2D BVP",
            lw=2.0)
    ax.plot(d_grid, U_D,    "s--", color="#9333EA",  label="Derjaguin",
            alpha=0.7)
    ax.plot(d_grid, U_A,    "d:",  color="#10B981",  label="additive",
            alpha=0.7)
    ax.plot(d_grid, U_blend, "-",  color="#2563EB",  label="blended (current)",
            alpha=0.8)
    ax.set_xlabel("d (surface-to-surface)")
    ax.set_ylabel("U(d)")
    ax.set_title(f"chi_b={p.chi_b}, omega_b={p.omega_b}, phi_inf={p.phi_inf}")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Right panel: midplane phi(rho) at d = R
    ax = axes[1]
    res_mid = solve_two_spheres(p, p.R, Nrho=Nrho, Nz=_Nz(p.R),
                                z_buf=4.0, rho_max=10.0, verbose=False)
    geom = res_mid["geom"]
    j_mid = np.argmin(np.abs(geom.z_grid - (geom.z_c1 + geom.z_c2) / 2.0))
    ax.plot(geom.rho_grid, res_mid["phi"][:, j_mid], "o-", color="C2",
            label="phi(rho, z=midplane) at d=R")
    ax.axhline(p.phi_inf, ls="--", color="gray", alpha=0.6, label="phi_inf")
    ax.set_xlabel("rho")
    ax.set_ylabel("phi")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    _run_main()
