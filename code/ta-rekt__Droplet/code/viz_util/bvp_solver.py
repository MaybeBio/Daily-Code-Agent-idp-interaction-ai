"""
BVP solver for bulk-surface equilibrium.

Models a Flory-Huggins bulk with a linear surface adsorption energy
    J(phi) = omega_b * phi
so that dJ/dphi = omega_b (a constant).  There is no membrane degree of freedom.

Governing equation (planar or spherical geometry):
    kappa * [phi'' + (2/r) * phi'] = df_b/dphi - mu_inf / nb

Boundary conditions:
    -kappa * phi'(0) + omega_b = 0    (surface adsorption)
     phi(L)          = phi_inf        (far-field reservoir)

Key numerical features:
1. Regularised logarithms with smooth clipping
2. Continuation methods for parameter sweeps
3. Multiple initial-guess strategies
4. Analytical test cases for validation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq
from dataclasses import dataclass
from typing import Tuple, Optional
import warnings


@dataclass
class Parameters:
    """Physical parameters for the bulk-surface system."""

    # Bulk molecule volume relative to solvent
    nb: float = 1.0

    # Bulk Flory-Huggins parameter
    chi_b: float = 0.0

    # Surface adsorption energy: J = omega_b * phi,  dJ/dphi = omega_b
    omega_b: float = -0.0

    # Gradient-squared penalty
    kappa: float = 1.0

    # Far-field (reservoir) volume fraction
    phi_inf: float = 0.06

    # Domain half-length (distance from surface to far-field)
    L: float = 10.0

    # Sphere radius.  R = inf => planar geometry; finite R => spherical.
    # Radial coordinate: r = R + z, z = distance from surface.
    R: float = 3

    # Regularisation floor for log arguments
    eps: float = 1e-10

    # Closed-system mode: conserve the total amount of phi in [0, L]
    # rather than fixing phi(L) = phi_inf.  When True, phi_inf is used
    # only to set the conserved total Phi_tot = phi_inf * integral_of_weight,
    # and the far-field BC is replaced by phi'(L) = 0 plus the integral
    # constraint.  mu_inf becomes an unknown parameter solved for by BVP.
    conserve_volume: bool = False

    # Parallel plate mode: solve the planar problem with two plates
    parallel_plates: bool = False

    @property
    def is_spherical(self) -> bool:
        """True when R is finite (spherical geometry)."""
        return np.isfinite(self.R)


def clip01(x: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """Hard clip to (eps, 1-eps)."""
    return np.clip(x, eps, 1.0 - eps)


def regularized_log(x: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """
    Regularized logarithm: extends log(x) smoothly for x < eps.
    Uses Taylor expansion: log(eps) + (x - eps)/eps for x < eps.
    """
    x = np.asarray(x)
    result = np.where(
        x > eps,
        np.log(np.maximum(x, eps)),
        np.log(eps) + (x - eps) / eps
    )
    return result


def fb(phi: np.ndarray, p: Parameters) -> np.ndarray:
    """Flory-Huggins bulk free energy density."""
    phi_c = clip01(phi, p.eps)
    return (1.0/p.nb) * (phi_c * regularized_log(phi_c, p.eps) +
                         (1-phi_c) * regularized_log(1-phi_c, p.eps)) + \
           p.chi_b * phi_c * (1 - phi_c)


def dfb_dphi(phi: np.ndarray, p: Parameters) -> np.ndarray:
    """First derivative of bulk free energy w.r.t. phi."""
    phi_c = clip01(phi, p.eps)
    return (1.0/p.nb) * (regularized_log(phi_c, p.eps) - regularized_log(1-phi_c, p.eps)) + \
           p.chi_b * (1.0 - 2.0*phi_c)


def d2fb_dphi2(phi: np.ndarray, p: Parameters) -> np.ndarray:
    """Second derivative of bulk free energy (for stability analysis)."""
    phi_c = clip01(phi, p.eps)
    return (1.0/p.nb) * (1.0/phi_c + 1.0/(1-phi_c)) - 2.0*p.chi_b


def psi_linear(z: np.ndarray, p: Parameters):
    """Linearized solution to the equilibrium concentration profile"""
    alpha = d2fb_dphi2(p.phi_inf, p)
    h0 = p.omega_b
    lam = np.sqrt(alpha/p.kappa)
    psi = ( - h0 * p.R * np.exp(- lam * z) 
            / ((z + p.R) * p.kappa * (lam + 1 / p.R)))
    dpsi = psi * (- lam - 1 / (z + p.R))
    return psi, dpsi, lam

def excess_free_energy(phi: np.ndarray, dphi: np.ndarray, z: np.ndarray,
                       p: Parameters, mu_ref: float = None) -> float:
    """
    Grand-potential excess over the bulk reservoir at phi_inf:

        DeltaOmega = int_V [ df(phi) + (kappa/2)(dphi/dz)^2 ] dV
                     + omega_b * int_{dV} (phi - phi_inf) dA

        df(phi) = f_b(phi) - f_b(phi_inf) - mu_ref/n_b * (phi - phi_inf)

    With mu_ref = n_b * dfb/dphi(phi_inf), df vanishes (to linear order) at
    phi_inf so the integral converges as L -> infinity.  This is the single
    energy functional used by both the Derjaguin (plate-plate) and additive
    (cylindrical) pair-potential calculations.

    Geometry is dispatched by `p`:
      - Spherical (p.is_spherical): bulk 4*pi*int (R+z)^2 [..] dz,
        surface 4*pi*R^2 * omega_b * (phi(0) - phi_inf).
      - Planar, single plate (R=inf, parallel_plates=False): bulk per-area
        int [..] dz, surface omega_b * (phi(0) - phi_inf) at z=0.
      - Parallel plates (R=inf, parallel_plates=True): bulk per-area
        int_0^L [..] dz, two surfaces omega_b * (phi(0) + phi(L) - 2 phi_inf).
    """
    phi_inf = p.phi_inf
    if mu_ref is None:
        mu_ref = p.nb * dfb_dphi(np.asarray([phi_inf], dtype=float), p)[0]

    df = (fb(phi, p) - fb(np.asarray([phi_inf], dtype=float), p)[0]
          - mu_ref / p.nb * (phi - phi_inf))
    integrand_per_length = df + 0.5 * p.kappa * dphi ** 2

    if p.is_spherical:
        weight = 4.0 * np.pi * (p.R + z) ** 2
        bulk = np.trapezoid(weight * integrand_per_length, z)
        surface = 4.0 * np.pi * p.R ** 2 * p.omega_b * (phi[0] - phi_inf)
    else:
        bulk = np.trapezoid(integrand_per_length, z)
        if p.parallel_plates:
            surface = p.omega_b * ((phi[0] - phi_inf) + (phi[-1] - phi_inf))
        else:
            surface = p.omega_b * (phi[0] - phi_inf)

    return bulk + surface


# Back-compat shim: old callers used `F(phi, dphi, d, p)`.  Delegate to the
# new functional so they get the corrected energy.
def F(phi: np.ndarray, dphi: np.ndarray, z: np.ndarray, p: Parameters) -> float:
    return excess_free_energy(phi, dphi, z, p)

def compute_spinodal(p: Parameters) -> Tuple[float, float]:
    """
    Compute the spinodal boundaries where d2fb/dphi2 = 0.
    For chi_b > 2 (with nb=1), the system can phase separate.
    Returns (phi_low, phi_high) or (None, None) if no spinodal.
    """
    # Spinodal condition: 1/(nb*phi) + 1/(nb*(1-phi)) = 2*chi_b
    # => 1/(phi(1-phi)) = 2*chi_b*nb
    # => phi(1-phi) = 1/(2*chi_b*nb)
    # => phi^2 - phi + 1/(2*chi_b*nb) = 0
    if p.chi_b <= 0:
        return None, None  # No spinodal for non-positive chi
    discriminant = 1 - 2/(p.chi_b * p.nb)
    if discriminant < 0:
        return None, None  # No spinodal - system is stable
    sqrt_d = np.sqrt(discriminant)
    phi_low = 0.5 * (1 - sqrt_d)
    phi_high = 0.5 * (1 + sqrt_d)
    return phi_low, phi_high


def compute_binodal(p: Parameters) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute binodal (coexistence) compositions via common-tangent construction:
        dfb/dphi(phi_a) = dfb/dphi(phi_b)
        fb(phi_b) - fb(phi_a) = (phi_b - phi_a) * dfb/dphi(phi_a)
    Returns (phi_low, phi_high) or (None, None) if no two-phase region exists.
    """
    from scipy.optimize import fsolve

    phi_spin_low, phi_spin_high = compute_spinodal(p)
    if phi_spin_low is None:
        return None, None

    def equations(x):
        phi_a, phi_b = x
        mu_a = dfb_dphi(phi_a, p)
        mu_b = dfb_dphi(phi_b, p)
        return [mu_a - mu_b,
                (fb(phi_b, p) - fb(phi_a, p)) - (phi_b - phi_a) * mu_a]

    x0 = [0.5 * phi_spin_low, 1.0 - 0.5 * (1.0 - phi_spin_high)]
    try:
        x, _, ier, _ = fsolve(equations, x0, full_output=True, xtol=1e-12)
    except Exception:
        return None, None
    if ier != 1:
        return None, None
    phi_a, phi_b = sorted(float(v) for v in x)
    if phi_a <= p.eps or phi_b >= 1.0 - p.eps or (phi_b - phi_a) < 1e-4:
        return None, None
    return phi_a, phi_b


def is_in_spinodal(phi: float, p: Parameters) -> bool:
    """Check if phi is in the spinodal (unstable) region."""
    phi_low, phi_high = compute_spinodal(p)
    if phi_low is None:
        return False
    return phi_low < phi < phi_high


def compute_mu_inf(p: Parameters) -> float:
    """Reservoir chemical potential: mu_inf = nb * dfb/dphi(phi_inf)."""
    return p.nb * dfb_dphi(p.phi_inf, p)


def _integral_weight(z, p: Parameters):
    """Weight w(z) for the conserved-mode integral: 1 planar, (R+z)^2 spherical."""
    if p.is_spherical:
        return (p.R + z) ** 2
    return np.ones_like(z) if hasattr(z, 'shape') else 1.0


def _domain_volume(p: Parameters) -> float:
    """Integral of the geometric weight over [0, L]."""
    if p.is_spherical:
        return ((p.R + p.L) ** 3 - p.R ** 3) / 3.0
    return p.L


def compute_phi_tot(p: Parameters) -> float:
    """Conserved total: phi_inf * integral of weight over the domain."""
    return p.phi_inf * _domain_volume(p)


def make_ode(p: Parameters, mu_inf: float):
    """
    Return the ODE function for scipy.integrate.solve_bvp.

    Open-reservoir mode (p.conserve_volume = False):
        State y = [phi, phi'].  mu_inf is fixed (passed in).

    Conserved mode (p.conserve_volume = True):
        State y = [phi, phi', I] where I' = w(z) * phi.
        mu_inf is an unknown solve_bvp parameter, read from p_param[0].
        The mu_inf argument is used only as a seed.

    Planar (R = inf):
        kappa * phi'' = dfb/dphi - mu_inf / nb

    Spherical (finite R), with r = R + z:
        kappa * (phi'' + (2/r) * phi') = dfb/dphi - mu_inf / nb
    """
    if p.conserve_volume:
        def ode(z, y, p_param):
            phi  = y[0]
            dphi = y[1]
            mu   = p_param[0]
            driving = dfb_dphi(phi, p) - mu / p.nb
            if p.is_spherical:
                r     = p.R + z
                ddphi = driving / p.kappa - 2.0 * dphi / r
                w     = r * r
            else:
                ddphi = driving / p.kappa
                w     = np.ones_like(z)
            return np.vstack((dphi, ddphi, w * phi))
        return ode

    def ode(z, y):
        phi  = y[0]
        dphi = y[1]
        driving = dfb_dphi(phi, p) - mu_inf / p.nb
        if p.is_spherical:
            r     = p.R + z
            ddphi = driving / p.kappa - 2.0 * dphi / r
        else:
            ddphi = driving / p.kappa
        return np.vstack((dphi, ddphi))

    return ode


def make_bc(p: Parameters):
    """
    Return the boundary-condition function for scipy.integrate.solve_bvp.

    Open-reservoir mode:
        BC[0]: -kappa * phi'(0) + omega_b = 0   (surface)
        BC[1]:  phi(L) = phi_inf                (far-field reservoir)

    Conserved mode (p.conserve_volume = True):
        BC[0]: -kappa * phi'(0) + omega_b = 0   (surface)
        BC[1]:  phi'(L) = 0                     (no flux at outer boundary)
        BC[2]:  I(0) = 0                        (integral starts at zero)
        BC[3]:  I(L) = Phi_tot                  (volume conservation)
    """
    if p.conserve_volume:
        phi_tot = compute_phi_tot(p)

        if not p.is_spherical and p.parallel_plates:
            def bc(ya, yb, p_param):
                return np.array([
                    -p.kappa * ya[1] + p.omega_b,
                    -p.kappa * yb[1] - p.omega_b,
                    ya[2],
                    yb[2] - phi_tot,
                ])
            return bc

        def bc(ya, yb, p_param):
            return np.array([
                -p.kappa * ya[1] + p.omega_b,
                yb[1],
                ya[2],
                yb[2] - phi_tot,
            ])
        return bc
    
    if not p.is_spherical and p.parallel_plates:
        # Open-reservoir parallel plates: surface BC on both plates.
        # mu_inf is fixed (computed from phi_inf and passed via the ODE),
        # so phi(0), phi(L) emerge from the solver.  No free parameter.
        def bc(ya, yb):
            return np.array([
                -p.kappa * ya[1] + p.omega_b,
                +p.kappa * yb[1] + p.omega_b,
            ])
        return bc

    def bc(ya, yb):
        dphi0 = ya[1]
        return np.array([
            -p.kappa * dphi0 + p.omega_b + 2e-2 * regularized_log(ya[0]/(1-ya[0])),
            yb[0] - p.phi_inf,
        ])

    return bc


def generate_initial_guess(
    p: Parameters, strategy: str = 'exponential'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build an initial guess (z, y_guess) where y_guess = [phi, phi'].

    Strategies: 'exponential', 'depleted', 'uniform', 'high_enrichment', 'tanh'.
    """
    z = np.linspace(0.0, p.L, 400)

    phi_spin_low, phi_spin_high = compute_spinodal(p)

    # Characteristic decay length at phi_inf
    d2fb = d2fb_dphi2(p.phi_inf, p)
    ell  = np.sqrt(p.kappa / d2fb) if d2fb > 0.1 else np.sqrt(p.kappa * p.nb / 0.1)
    ell  = min(ell, p.L / 2.0)

    if strategy == 'exponential':
        if phi_spin_low is not None and p.phi_inf < phi_spin_low:
            phi0_g = min(0.8, p.phi_inf + 0.15)
        else:
            phi0_g = min(0.95, max(0.05, p.phi_inf + 0.1))
        phi_g  = p.phi_inf + (phi0_g - p.phi_inf) * np.exp(-z / ell)
        dphi_g = -(phi0_g - p.phi_inf) / ell * np.exp(-z / ell)

    elif strategy == 'depleted':
        phi0_g = max(0.02, p.phi_inf - 0.03)
        phi_g  = p.phi_inf + (phi0_g - p.phi_inf) * np.exp(-z / ell)
        dphi_g = -(phi0_g - p.phi_inf) / ell * np.exp(-z / ell)

    elif strategy == 'uniform':
        phi_g  = np.full_like(z, p.phi_inf)
        dphi_g = np.zeros_like(z)

    elif strategy == 'high_enrichment':
        phi0_g = min(0.95, phi_spin_high + 0.1) if phi_spin_high is not None else 0.7
        phi_g  = p.phi_inf + (phi0_g - p.phi_inf) * np.exp(-z / ell)
        dphi_g = -(phi0_g - p.phi_inf) / ell * np.exp(-z / ell)

    elif strategy == 'tanh':
        z0     = p.L / 4.0
        width  = np.sqrt(p.kappa / max(0.1, p.chi_b))
        phi_hi = min(0.9, p.phi_inf + 0.3)
        phi_g  = p.phi_inf + (phi_hi - p.phi_inf) * 0.5 * (1.0 - np.tanh((z - z0) / width))
        dphi_g = np.gradient(phi_g, z)

    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    if p.conserve_volume:
        from scipy.integrate import cumulative_trapezoid
        w_g = _integral_weight(z, p)
        I_g = cumulative_trapezoid(w_g * phi_g, z, initial=0.0)
        return z, np.vstack((phi_g, dphi_g, I_g))

    return z, np.vstack((phi_g, dphi_g))


def solve_bvp_robust(
    p: Parameters, verbose: bool = True, tol: float = 1e-6
) -> Optional[object]:
    """
    Robust BVP solver.  Tries several initial-guess strategies and tolerances.
    Returns the first physically valid (phi in (0,1)) converged solution.
    """
    mu_inf = compute_mu_inf(p)
    ode    = make_ode(p, mu_inf)
    bc     = make_bc(p)

    phi_spin_low, phi_spin_high = compute_spinodal(p)
    if phi_spin_low is not None and verbose:
        print(f"Note: Spinodal region at phi in ({phi_spin_low:.3f}, {phi_spin_high:.3f})")
        if phi_spin_low < p.phi_inf < phi_spin_high:
            print("WARNING: phi_inf is in the spinodal region — system is unstable!")

    best_sol      = None
    best_residual = float('inf')
    strategies = ['exponential', 'high_enrichment', 'depleted', 'tanh', 'uniform']
    # strategies = ['exponential']

    for strategy in strategies:
        try:
            z, y_guess = generate_initial_guess(p, strategy)
            y_guess[0] = clip01(y_guess[0], p.eps)

            for current_tol in [tol, tol * 10, tol * 100]:
                try:
                    kwargs = dict(tol=current_tol, max_nodes=100000, verbose=0)
                    if p.conserve_volume:
                        kwargs['p'] = [mu_inf]
                    sol = solve_bvp(ode, bc, z, y_guess, **kwargs)
                    if sol.success:
                        phi_min = np.min(sol.y[0])
                        phi_max = np.max(sol.y[0])
                        if phi_min > 0 and phi_max < 1:
                            residual = (
                                np.max(np.abs(sol.rms_residuals))
                                if hasattr(sol, 'rms_residuals') else 0.0
                            )
                            if residual < best_residual:
                                best_sol      = sol
                                best_residual = residual
                            if verbose:
                                print(
                                    f"Converged with strategy '{strategy}'"
                                    f" (tol={current_tol:.0e})"
                                )
                                print(f"  phi(0) = {sol.y[0, 0]:.6f}")
                            return sol
                except Exception:
                    continue

        except Exception as e:
            if verbose:
                print(f"Strategy '{strategy}' failed: {e}")
            continue

    if best_sol is not None:
        return best_sol
    if verbose:
        print("All strategies failed")
    return None


def solve_with_continuation(
    p_start: Parameters,
    p_end: Parameters,
    n_steps: int = 10,
    verbose: bool = True,
    adaptive: bool = True,
) -> Optional[object]:
    """
    Solve via parameter continuation from p_start to p_end.
    Useful for reaching difficult parameter regimes.

    Args:
        adaptive: If True, halve the step size on failure and retry.
    """

    def interp_params(t: float) -> Parameters:
        if np.isinf(p_start.R) and np.isinf(p_end.R):
            R_interp = float('inf')
        elif np.isinf(p_start.R):
            R_interp = p_end.R
        elif np.isinf(p_end.R):
            R_interp = p_start.R
        else:
            R_interp = p_start.R + t * (p_end.R - p_start.R)
        return Parameters(
            nb      = p_start.nb      + t * (p_end.nb      - p_start.nb),
            chi_b   = p_start.chi_b   + t * (p_end.chi_b   - p_start.chi_b),
            omega_b = p_start.omega_b + t * (p_end.omega_b - p_start.omega_b),
            kappa   = p_start.kappa   + t * (p_end.kappa   - p_start.kappa),
            phi_inf = p_start.phi_inf + t * (p_end.phi_inf - p_start.phi_inf),
            L       = p_start.L       + t * (p_end.L       - p_start.L),
            R       = R_interp,
            eps     = p_start.eps,
            conserve_volume = p_start.conserve_volume,
        )

    def try_step(sol, t_current, t_target, tol=1e-5):
        p_c    = interp_params(t_target)
        mu_inf = compute_mu_inf(p_c)
        ode    = make_ode(p_c, mu_inf)
        bc     = make_bc(p_c)

        z_new  = np.linspace(0.0, p_c.L, 500)
        p_prev = interp_params(t_current)
        scale  = p_prev.L / p_c.L if t_current > 0 else 1.0
        y_new  = sol.sol(np.minimum(z_new * scale, p_prev.L * 0.999))
        y_new[0] = clip01(y_new[0], p_c.eps)

        for tol_try in [tol, tol * 10, tol * 100]:
            try:
                kwargs = dict(tol=tol_try, max_nodes=50000, verbose=0)
                if p_c.conserve_volume:
                    mu_prev = sol.p[0] if getattr(sol, 'p', None) is not None else mu_inf
                    kwargs['p'] = [mu_prev]
                sol_new = solve_bvp(ode, bc, z_new, y_new, **kwargs)
                if sol_new.success and np.min(sol_new.y[0]) > 0 and np.max(sol_new.y[0]) < 1:
                    return sol_new
            except Exception:
                continue
        return None

    # Initialise at t = 0
    sol = solve_bvp_robust(interp_params(0), verbose=False, tol=1e-5)
    if sol is None:
        if verbose:
            print("Failed to solve at starting point")
        return None
    if verbose:
        print(f"Step t=0.000: phi(0) = {sol.y[0, 0]:.6f}")

    t_current = 0.0
    dt_base   = 1.0 / n_steps
    dt        = dt_base
    dt_min    = dt_base / 100

    while t_current < 1.0 - 1e-10:
        t_target = min(t_current + dt, 1.0)
        sol_new  = try_step(sol, t_current, t_target)

        if sol_new is not None:
            sol       = sol_new
            t_current = t_target
            if verbose:
                print(f"Step t={t_current:.3f}: phi(0) = {sol.y[0, 0]:.6f}")
            if adaptive and dt < dt_base:
                dt = min(dt * 1.5, dt_base)
        else:
            if adaptive and dt > dt_min:
                dt = dt / 2.0
                if verbose:
                    print(f"  Refining step size to dt={dt:.4f}")
            else:
                if verbose:
                    print(f"Continuation failed at t={t_current:.3f}")
                break

    if t_current >= 1.0 - 1e-6:
        return sol
    if verbose:
        print(f"Continuation incomplete: reached t={t_current:.3f}")
    return None


def solve_shooting(p: Parameters, verbose: bool = True) -> Optional[object]:
    """
    Shooting method.

    Because J = omega_b * phi is linear, the surface BC pins phi'(0) exactly:
        phi'(0) = omega_b / kappa   (independent of phi(0))

    Finding the correct phi(0) is therefore a 1-D root-finding problem:
    shoot from (phi0, omega_b/kappa) and match phi(L) = phi_inf.

    Not implemented in conserved-volume mode (use solve_bvp_robust).
    """
    if p.conserve_volume:
        if verbose:
            print("Shooting not supported with conserve_volume=True")
        return None
    mu_inf = compute_mu_inf(p)
    dphi0  = p.omega_b / p.kappa   # fixed by surface BC

    def phi_at_L(phi0: float) -> float:
        phi0_c = float(clip01(np.array([phi0]), p.eps)[0])

        def ode_ivp(z, y):
            ddphi = (dfb_dphi(y[0], p) - mu_inf / p.nb) / p.kappa
            return [y[1], ddphi]

        try:
            sol = solve_ivp(
                ode_ivp, [0.0, p.L], [phi0_c, dphi0],
                method='RK45', max_step=p.L / 100,
                t_eval=[p.L], rtol=1e-9, atol=1e-11,
            )
            return float(sol.y[0, -1]) if sol.success else np.nan
        except Exception:
            return np.nan

    # Scan phi0 in (0, 1) to locate sign changes in phi(L) - phi_inf
    phi0_grid = np.linspace(0.01, 0.99, 200)
    residuals  = np.array([phi_at_L(v) - p.phi_inf for v in phi0_grid])

    brackets = [
        (phi0_grid[i], phi0_grid[i + 1])
        for i in range(len(residuals) - 1)
        if np.isfinite(residuals[i]) and np.isfinite(residuals[i + 1])
        and residuals[i] * residuals[i + 1] < 0
    ]

    best_phi0     = None
    best_residual = float('inf')

    for a, b in brackets:
        try:
            phi0_sol = brentq(lambda v: phi_at_L(v) - p.phi_inf, a, b, xtol=1e-12)
            res = abs(phi_at_L(phi0_sol) - p.phi_inf)
            if res < best_residual:
                best_phi0     = phi0_sol
                best_residual = res
        except Exception:
            continue

    if best_phi0 is None or best_residual > 1e-6:
        if verbose:
            print("Shooting method failed")
        return None

    phi0 = float(clip01(np.array([best_phi0]), p.eps)[0])
    if verbose:
        print(f"Shooting converged: phi(0) = {phi0:.6f}, |residual| = {best_residual:.2e}")

    # Build a solution object compatible with solve_bvp output
    class ShootingSolution:
        success = True

        def __init__(self):
            def ode_ivp(z, y):
                ddphi = (dfb_dphi(y[0], p) - mu_inf / p.nb) / p.kappa
                return [y[1], ddphi]

            z_eval   = np.linspace(0.0, p.L, 400)
            ivp_sol  = solve_ivp(
                ode_ivp, [0.0, p.L], [phi0, dphi0],
                method='RK45', max_step=p.L / 100,
                t_eval=z_eval, rtol=1e-9, atol=1e-11,
            )
            self.x = ivp_sol.t
            self.y = ivp_sol.y          # shape (2, 400)

        def sol(self, z):               # mimic solve_bvp's .sol interpolant
            return np.array([
                np.interp(z, self.x, self.y[0]),
                np.interp(z, self.x, self.y[1]),
            ])

    return ShootingSolution()


def solve_auto(p: Parameters, verbose: bool = True) -> Optional[object]:
    """
    Try solvers in order: BVP -> shooting -> continuation from easy params.
    Returns the first solution found, or None.
    """
    sol = solve_bvp_robust(p, verbose=verbose)
    if sol is not None:
        return sol

    if verbose:
        print("Trying shooting method...")
    sol = solve_shooting(p, verbose=verbose)
    if sol is not None:
        return sol

    if verbose:
        print("Trying continuation method...")
    p_easy = Parameters(
        nb      = p.nb,
        chi_b   = min(p.chi_b, 1.5),
        omega_b = p.omega_b * 0.1,
        kappa   = p.kappa,
        phi_inf = p.phi_inf,
        L       = p.L,
        R       = p.R,
        eps     = p.eps,
        conserve_volume = p.conserve_volume,
    )
    sol = solve_with_continuation(p_easy, p, n_steps=20, verbose=verbose)
    if sol is not None:
        return sol

    if verbose:
        print("Trying continuation with more steps...")
    return solve_with_continuation(p_easy, p, n_steps=50, verbose=False)


# ============================================================
# TEST CASES
# ============================================================


def test_uniform_solution():
    """
    TEST 1: omega_b = 0.

    The surface BC becomes phi'(0) = 0.  In the stable single-phase regime
    the unique solution satisfying both BCs is phi(z) = phi_inf.
    """
    print("\n" + "=" * 60)
    print("TEST 1: Uniform solution (omega_b = 0)")
    print("=" * 60)
    
    p = Parameters(
        chi_b   = 2.0,   # at the critical value for nb=1, still single-phase
        omega_b = 0.0,
        phi_inf = 0.3,
        kappa   = 1.0,
        L       = 10.0,
    )

    sol = solve_bvp_robust(p, verbose=True)
    if sol is None or not sol.success:
        print("  SOLVER FAILED ✗")
        return False

    z_fine  = np.linspace(0.0, p.L, 200)
    phi     = sol.sol(z_fine)[0]
    phi_std = np.std(phi)
    phi_err = abs(np.mean(phi) - p.phi_inf)

    print(f"  Expected phi = {p.phi_inf:.6f}")
    print(f"  Mean phi     = {np.mean(phi):.6f}")
    print(f"  Std phi      = {phi_std:.2e}")

    passed = phi_std < 1e-4 and phi_err < 1e-3
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_analytical_surface_enrichment():
    """
    TEST 2: Linearised analytical prediction for weak adsorption.

    For chi_b = 0, nb = 1, the ODE linearises around phi_inf to
        u'' = u / ell^2,    ell = sqrt(kappa * phi_inf * (1 - phi_inf))
    with BCs
        u'(0) = omega_b / kappa    (= dphi0, fixed by surface BC)
        u(L)  = 0

    For L >> ell the exact solution collapses to
        u(z) = -omega_b * ell / kappa * exp(-z / ell)
    giving
        phi(0) = phi_inf  -  omega_b * ell / kappa.

    The nonlinear correction is O(omega_b^2).  We use omega_b = -0.05
    (weak, attractive) and verify the numerical phi(0) matches to within
    10 * omega_b^2.
    """
    print("\n" + "=" * 60)
    print("TEST 2: Analytical surface enrichment (chi_b=0, linearised)")
    print("=" * 60)

    phi_inf = 0.3
    omega_b = -0.05   # weak attractive coupling
    kappa   = 1.0
    L       = 30.0    # L >> ell so large-L approximation is accurate

    p = Parameters(
        chi_b   = 0.0,
        omega_b = omega_b,
        phi_inf = phi_inf,
        kappa   = kappa,
        L       = L,
    )

    # Analytical prediction: phi(0) = phi_inf - omega_b * ell / kappa
    ell         = np.sqrt(kappa * phi_inf * (1.0 - phi_inf))
    phi0_theory = phi_inf - omega_b * ell / kappa

    print(f"  ell               = {ell:.4f}")
    print(f"  Analytical phi(0) = {phi0_theory:.6f}")

    sol = solve_bvp_robust(p, verbose=False)
    if sol is None or not sol.success:
        sol = solve_shooting(p, verbose=False)
    if sol is None or not sol.success:
        print("  SOLVER FAILED ✗")
        return False

    phi0_num = float(sol.sol(np.array([0.0]))[0, 0])
    err      = abs(phi0_num - phi0_theory)
    print(f"  Numerical  phi(0) = {phi0_num:.6f}")
    print(f"  Discrepancy       = {err:.2e}  (expected O(omega_b^2) ~ {omega_b**2:.2e})")

    passed = err < 10.0 * omega_b ** 2
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_residual_check():
    """
    TEST 3: Verify the solved profile satisfies the ODE and both BCs.

    Checks:
      - Interior ODE: kappa * phi''(z) = dfb/dphi - mu_inf/nb
      - Surface BC:   -kappa * phi'(0) + omega_b = 0
      - Far-field BC:  phi(L) = phi_inf
    """
    print("\n" + "=" * 60)
    print("TEST 3: Residual check")
    print("=" * 60)

    p = Parameters(
        chi_b   = 2.5,
        omega_b = -0.2,
        phi_inf = 0.1,
        kappa   = 1.0,
        L       = 10.0,
    )

    sol = solve_auto(p, verbose=False)
    if sol is None or not sol.success:
        print("  SOLVER FAILED ✗")
        return False

    mu_inf = compute_mu_inf(p)

    # Interior ODE residual via centred finite differences
    z_test    = np.linspace(0.02 * p.L, 0.98 * p.L, 100)
    h         = 0.05
    phi_c,  _      = sol.sol(z_test)
    _,      dphi_p = sol.sol(z_test + h)
    _,      dphi_m = sol.sol(z_test - h)
    ddphi          = (dphi_p - dphi_m) / (2.0 * h)
    ode_res        = np.max(np.abs(p.kappa * ddphi - (dfb_dphi(phi_c, p) - mu_inf / p.nb)))

    # Boundary conditions
    phi0, dphi0 = sol.sol(np.array([0.0]))[:, 0]
    phiL        = sol.sol(np.array([p.L]))[0, 0]
    bc_surf     = abs(-p.kappa * dphi0 + p.omega_b)
    bc_far      = abs(phiL - p.phi_inf)

    print(f"  Max ODE residual:    {ode_res:.2e}")
    print(f"  Surface BC residual: {bc_surf:.2e}")
    print(f"  Far-field residual:  {bc_far:.2e}")

    tol    = 1e-4
    passed = ode_res < tol and bc_surf < tol and bc_far < tol
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_ideal_solution():
    """
    TEST 4: chi_b = 0, omega_b = 0 (ideal, non-interacting).

    The bulk is ideal and the surface is inert.  The only solution
    satisfying phi'(0) = 0 and phi(L) = phi_inf is phi(z) = phi_inf.
    """
    print("\n" + "=" * 60)
    print("TEST 4: Ideal solution (chi_b=0, omega_b=0)")
    print("=" * 60)

    p = Parameters(
        chi_b   = 0.0,
        omega_b = 0.0,
        phi_inf = 0.2,
        kappa   = 1.0,
        L       = 10.0,
    )

    sol = solve_bvp_robust(p, verbose=True)
    if sol is None or not sol.success:
        print("  SOLVER FAILED ✗")
        return False

    z_fine = np.linspace(0.0, p.L, 200)
    phi    = sol.sol(z_fine)[0]
    err    = np.max(np.abs(phi - p.phi_inf))

    print(f"  Max deviation from phi_inf: {err:.2e}")
    passed = err < 1e-4
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_domain_independence():
    """
    TEST 5: Near-surface profile is independent of L (for large L).

    Solve with L = 10 and L = 20 and compare phi(z) for z in [0, 3].
    The profiles should agree to within 1%.
    """
    print("\n" + "=" * 60)
    print("TEST 5: Domain independence (L=10 vs L=20)")
    print("=" * 60)

    base = dict(chi_b=2.5, omega_b=-0.2, phi_inf=0.1, kappa=1.0)
    p1   = Parameters(L=10.0, **base)
    p2   = Parameters(L=20.0, **base)

    sol1 = solve_auto(p1, verbose=False)
    sol2 = solve_auto(p2, verbose=False)

    if sol1 is None or sol2 is None or not sol1.success or not sol2.success:
        print("  SOLVER FAILED ✗")
        return False

    z_compare = np.linspace(0.0, 3.0, 50)
    phi1      = sol1.sol(z_compare)[0]
    phi2      = sol2.sol(z_compare)[0]
    diff      = np.max(np.abs(phi1 - phi2))

    print(f"  L=10: phi(0) = {sol1.y[0, 0]:.6f}")
    print(f"  L=20: phi(0) = {sol2.y[0, 0]:.6f}")
    print(f"  Max difference near surface: {diff:.2e}")

    passed = diff < 0.01
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_kappa_scaling():
    """
    TEST 6: Interface width scales as sqrt(kappa).

    The characteristic length ell ~ sqrt(kappa / d2fb/dphi2(phi_inf)).
    Compare kappa=1 vs kappa=4; expect width ratio ~ sqrt(4) = 2.
    """
    print("\n" + "=" * 60)
    print("TEST 6: Kappa scaling of interface width")
    print("=" * 60)

    base = dict(chi_b=1.5, omega_b=-0.3, phi_inf=0.3, L=30.0)
    p1   = Parameters(kappa=1.0, **base)
    p2   = Parameters(kappa=4.0, **base)

    sol1 = solve_auto(p1, verbose=False)
    sol2 = solve_auto(p2, verbose=False)

    if sol1 is None or sol2 is None or not sol1.success or not sol2.success:
        print("  SOLVER FAILED ✗")
        return False

    def estimate_width(sol, p):
        z      = np.linspace(0.0, p.L, 1000)
        phi, _ = sol.sol(z)
        delta  = abs(phi[0] - p.phi_inf)
        if delta < 1e-3:
            return 1.0
        return np.trapz(np.abs(phi - p.phi_inf), z) / delta

    w1       = estimate_width(sol1, p1)
    w2       = estimate_width(sol2, p2)
    ratio    = w2 / w1
    expected = np.sqrt(4.0)

    print(f"  kappa=1: width ~ {w1:.3f}")
    print(f"  kappa=4: width ~ {w2:.3f}")
    print(f"  Width ratio: {ratio:.3f}  (expected sqrt(4) = {expected:.3f})")

    passed = abs(ratio - expected) < 0.5
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_spherical_geometry():
    """
    TEST 7: Spherical geometry.

    Checks:
    1. Spherical solution converges.
    2. Very large R approaches the planar solution with O(L/R) error.
    3. Finite R = 5 differs noticeably from planar (curvature effect).
    4. ODE residual for the spherical Laplacian is small.
    """
    print("\n" + "=" * 60)
    print("TEST 7: Spherical geometry")
    print("=" * 60)

    base = dict(chi_b=1.5, omega_b=-0.2, phi_inf=0.2, kappa=1.0, L=15.0)
    p_planar    = Parameters(R=float('inf'), **base)
    p_spherical = Parameters(R=5.0,          **base)
    p_large_R   = Parameters(R=1000.0,       **base)

    sol_planar    = solve_auto(p_planar,    verbose=False)
    sol_spherical = solve_auto(p_spherical, verbose=False)
    sol_large_R   = solve_auto(p_large_R,   verbose=False)

    if any(s is None or not s.success for s in [sol_planar, sol_spherical, sol_large_R]):
        print("  SOLVER FAILED ✗")
        return False

    z = np.linspace(0.0, 15.0, 100)
    phi_planar,    _ = sol_planar.sol(z)
    phi_spherical, _ = sol_spherical.sol(z)
    phi_large_R,   _ = sol_large_R.sol(z)

    print(f"  Planar    (R=inf):  phi(0) = {phi_planar[0]:.4f}")
    print(f"  Spherical (R=5):    phi(0) = {phi_spherical[0]:.4f}")
    print(f"  Large R   (R=1000): phi(0) = {phi_large_R[0]:.4f}")

    diff_large_R   = np.max(np.abs(phi_large_R   - phi_planar))
    diff_spherical = np.max(np.abs(phi_spherical - phi_planar))
    print(f"  |planar - large_R|_max    = {diff_large_R:.6f}"
          f"  (expected O({p_large_R.L / p_large_R.R:.4f}))")
    print(f"  |planar - spherical|_max  = {diff_spherical:.6f}")

    if diff_large_R > 0.05:
        print("  FAILED: Large-R solution does not approach planar ✗")
        return False
    if diff_spherical < 0.005:
        print("  FAILED: Spherical should differ visibly from planar ✗")
        return False

    # ODE residual for the spherical Laplacian (finite-difference check)
    p   = p_spherical
    sol = sol_spherical
    z_t       = np.linspace(0.02 * p.L, 0.98 * p.L, 100)
    phi_t, dphi_t = sol.sol(z_t)
    mu_inf    = compute_mu_inf(p)
    driving   = dfb_dphi(phi_t, p) - mu_inf / p.nb
    dz        = z_t[1] - z_t[0]
    ddphi_fd  = np.gradient(np.gradient(phi_t, dz), dz)
    r         = p.R + z_t
    lap       = ddphi_fd + 2.0 * dphi_t / r
    ode_res   = np.max(np.abs(p.kappa * lap[5:-5] - driving[5:-5]))
    print(f"  Max spherical ODE residual (FD): {ode_res:.2e}")

    if ode_res > 0.1:
        print("  FAILED: Spherical ODE residual too large ✗")
        return False

    print("  TEST PASSED ✓")
    return True


def test_conserve_volume_ideal():
    """
    TEST 8: Conserved mode in the ideal/inert limit.

    With omega_b = 0, chi_b = 0 and conserve_volume = True the unique
    solution is phi(z) = phi_inf, and the recovered mu_inf matches
    nb * dfb/dphi(phi_inf).
    """
    print("\n" + "=" * 60)
    print("TEST 8: Conserved mode, ideal limit")
    print("=" * 60)

    p = Parameters(
        chi_b = 0.0, omega_b = 0.0, phi_inf = 0.2,
        kappa = 1.0, L = 10.0, R = float('inf'),
        conserve_volume = True,
    )
    sol = solve_bvp_robust(p, verbose=False)
    if sol is None or not sol.success:
        print("  SOLVER FAILED ✗")
        return False

    z = np.linspace(0.0, p.L, 200)
    phi = sol.sol(z)[0]
    mu_recovered = sol.p[0]
    mu_expected  = p.nb * dfb_dphi(p.phi_inf, p)

    err_phi = np.max(np.abs(phi - p.phi_inf))
    err_mu  = abs(mu_recovered - mu_expected)
    print(f"  Max |phi - phi_inf|: {err_phi:.2e}")
    print(f"  mu_inf recovered = {mu_recovered:.6f}, expected {mu_expected:.6f}")

    passed = err_phi < 1e-4 and err_mu < 1e-4
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def test_conserve_volume_mass():
    """
    TEST 9: Mass conservation under attractive surface.

    Strong attractive omega_b enriches the surface and depletes the bulk.
    The integral ∫₀ᴸ w(z) phi(z) dz must equal Phi_tot = phi_inf * domain_volume.
    """
    print("\n" + "=" * 60)
    print("TEST 9: Conserved mode, mass conservation (planar)")
    print("=" * 60)

    p = Parameters(
        chi_b = 1.5, omega_b = -0.3, phi_inf = 0.1,
        kappa = 1.0, L = 15.0, R = float('inf'),
        conserve_volume = True,
    )
    sol = solve_bvp_robust(p, verbose=False)
    if sol is None or not sol.success:
        print("  SOLVER FAILED ✗")
        return False

    z = np.linspace(0.0, p.L, 4000)
    phi = sol.sol(z)[0]
    w = _integral_weight(z, p)
    integral = np.trapz(w * phi, z)
    phi_tot  = compute_phi_tot(p)
    phiL     = float(sol.sol(np.array([p.L]))[0, 0])

    err = abs(integral - phi_tot) / abs(phi_tot)
    print(f"  ∫ w·phi dz = {integral:.6f}")
    print(f"  Phi_tot    = {phi_tot:.6f}")
    print(f"  Relative error: {err:.2e}")
    print(f"  phi(0) = {sol.y[0, 0]:.4f}, phi(L) = {phiL:.4f} "
          f"(should be < phi_inf={p.phi_inf})")

    passed = err < 1e-4 and sol.y[0, 0] > p.phi_inf and phiL < p.phi_inf
    print("  TEST PASSED ✓" if passed else "  TEST FAILED ✗")
    return passed


def run_all_tests():
    """Run all test cases and print a summary."""
    print("\n" + "=" * 60)
    print("RUNNING ALL TESTS")
    print("=" * 60)

    tests = [
        ("Uniform solution (omega_b=0)",               test_uniform_solution),
        ("Analytical surface enrichment (linearised)", test_analytical_surface_enrichment),
        ("Residual check",                             test_residual_check),
        ("Ideal solution (chi_b=0, omega_b=0)",        test_ideal_solution),
        ("Domain independence",                        test_domain_independence),
        ("Kappa scaling",                              test_kappa_scaling),
        ("Spherical geometry",                         test_spherical_geometry),
        ("Conserved mode, ideal limit",                test_conserve_volume_ideal),
        ("Conserved mode, mass conservation",          test_conserve_volume_mass),
    ]

    results = [(name, fn()) for name, fn in tests]

    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    passed = 0
    for name, ok in results:
        print(f"  {'✓' if ok else '✗'} {name}")
        passed += ok
    print(f"\nTotal: {passed}/{len(results)} tests passed")
    return passed == len(results)


# ============================================================
# PLOTTING UTILITIES
# ============================================================

def plot_solution(sol, p: Parameters, title: str = ""):
    """Plot phi(z) and phi'(z) profiles."""
    if sol is None or not sol.success:
        print("No valid solution to plot")
        return

    z_fine          = np.linspace(0.0, p.L, 500)
    phi_fine, dphi_fine = sol.sol(z_fine)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(z_fine, phi_fine, 'b-', linewidth=2, label=r'$\phi(z)$')
    axes[0].axhline(p.phi_inf, linestyle='--', color='gray', label=r'$\phi_\infty$')
    axes[0].set_xlabel("z")
    axes[0].set_ylabel(r"$\phi$")
    axes[0].set_ylim((0.0, 1.0))
    axes[0].legend()
    axes[0].set_title(f"Bulk profile {title}")

    axes[1].plot(z_fine, dphi_fine, 'g-', linewidth=2)
    axes[1].axhline(0.0, linestyle='--', color='gray')
    axes[1].set_xlabel("z")
    axes[1].set_ylabel(r"$d\phi/dz$")
    axes[1].set_title("Gradient profile")

    plt.tight_layout()
    plt.show()


def plot_phase_diagram(p: Parameters):
    """Plot bulk free energy and its first two derivatives."""
    phi = np.linspace(0.01, 0.99, 200)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].plot(phi, fb(phi, p), 'b-', linewidth=2)
    axes[0].set_xlabel(r'$\phi$')
    axes[0].set_ylabel(r'$f_b(\phi)$')
    axes[0].set_title('Bulk free energy')

    df = dfb_dphi(phi, p)
    axes[1].plot(phi, df, 'g-', linewidth=2)
    axes[1].axhline(compute_mu_inf(p) / p.nb, linestyle='--', color='red',
                    label=r'$\mu_\infty/\nu_b$')
    axes[1].set_xlabel(r'$\phi$')
    axes[1].set_ylabel(r'$df_b/d\phi$')
    axes[1].set_title('Chemical potential')
    axes[1].legend()

    d2f = d2fb_dphi2(phi, p)
    axes[2].plot(phi, d2f, 'r-', linewidth=2)
    axes[2].axhline(0.0, linestyle='--', color='gray')
    phi_low, phi_high = compute_spinodal(p)
    if phi_low is not None:
        axes[2].axvline(phi_low,  linestyle=':', color='orange', label='Spinodal')
        axes[2].axvline(phi_high, linestyle=':', color='orange')
        axes[2].fill_between([phi_low, phi_high], [-5, -5], [5, 5],
                             alpha=0.2, color='red')
    axes[2].set_xlabel(r'$\phi$')
    axes[2].set_ylabel(r'$d^2f_b/d\phi^2$')
    axes[2].set_title('Stability (negative = unstable)')
    axes[2].set_ylim((-5, 20))
    axes[2].legend()

    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    run_all_tests()

    print("\n" + "=" * 60)
    print("DEMO: Enrichment near attractive surface (chi_b=2.5, omega_b=-0.3)")
    print("=" * 60)

    p = Parameters(
        chi_b   = 2.5,
        omega_b = -0.3,
        phi_inf = 0.1,
        kappa   = 1.0,
        L       = 15.0,
    )
    sol = solve_auto(p, verbose=True)
    if sol is not None and sol.success:
        plot_solution(sol, p)
