"""
Langevin open/close spec: 1-D overdamped Langevin on the effective Rouse +
surface-attraction potential, scanned over concentration phi_inf.

Replaces the expensive full-Rouse sim (spec/rouse_mfpt_attractive.jl) with a fast
1-D surrogate on the end-to-end-distance potential U(r), and extracts empirical
Kramers exit (well -> entropic basin) and re-entry (basin -> well) rates from the
sojourn intervals, alongside the analytical Kramers times.

Parallelized axis: alphas (phi_inf), one value per CPU.  Inner axes: betas (R) x
gammas (omega_b).  Driver: sim/sweep.py.  Pipe builder: runPY.
"""

import numpy as np

from viz_util.bvp_solver import (
    Parameters,
    d2fb_dphi2,
)
from viz_util.overdamped_langevin import OverdampedHarmonicLangevin

blurb = (
    "Simulating overdamped Langevin on the Rouse + surface-attraction potential "
    "at various concentrations, to recover empirical Kramers exit/re-entry rates "
    "and feed the downstream Fisher-rate analysis."
)

# surface model fixed parameters
fixed = dict(
    nb=1.0,
    chi_b=2.5,
    kappa=1.0,
    L=10.0,
    eps=1e-6,
    conserve_volume=False,
    parallel_plates=False,
)

# Rouse polymer descriptors.  kuhnLength (=a) and beadSep (=s) enter the
# potential U(r); kBT and friction set the Langevin bath.  stiffness/D are
# recorded only.
kBT = 1.0
kuhnLength = 1.0          # Kuhn length (a)
friction = 1.0           # friction coefficient
N = 100                  # chain length (beads)
beadSep = N - 1          # number of bonds between tracked beads (s)
stiffness = 2 * kBT / kuhnLength ** 2
D = kBT / friction

# Langevin-simulation controls (notebook exploratory bath: T_bath != kBT).
T_bath = kBT
gamma = 1.0
kB = 1.0
n_walkers = 1000          # independent walkers (sojourn statistics)
nskip = 5000               # subsample factor for the stored trajectory
tau_mult = 200           # total time = tau_mult * tau_K
MAX_STEPS = 5e7          # safety cap on n = tau/dt (clip tau if a barrier is tall)

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.py.
alphas = np.linspace(0.008, 0.1, 32)   # phi_inf
# alphas = [0.1]
betas = [2.5]                            # R
gammas = [-0.5]                          # omega_b

nparams = len(alphas)

# Fixed, shared histogram grids (physical time) -- separate for exit vs re-entry
# so the 16 concentrations are directly comparable (Fisher open/close analysis).
exit_T_grid = np.logspace(-1, 4, 80)
reentry_T_grid = np.logspace(-1, 4, 80)

_PARAM_KEYS = {"nb", "chi_b", "kappa", "L", "eps", "conserve_volume", "parallel_plates"}


def _param_kwargs(fixed):
    return {k: v for k, v in fixed.items() if k in _PARAM_KEYS}


def U(r, *args):
    R, h0, kappa, lam, a, s = args
    profile = np.exp(-lam * (r - R)) * R / r

    return (
        - 2 * np.log(r) + 3 * r ** 2 / (2 * s * a ** 2)
        - 4 * np.pi * R ** 2 * h0 ** 2 * profile / (kappa * (lam + 1 / R))
    )

def d2U(x, *args):
        h = 1e-3
        return (U(x + h, *args) - 2 * U(x, *args) + U(x - h, *args)) / h ** 2


def F(r, *args):
    """Analytic force F = -dU/dr (only valid for r > 0)."""
    R, h0, kappa, lam, a, s = args

    return (
        2 / r - 3 * r / (s * a ** 2)
        - 4 * np.pi * R ** 3 * h0 ** 2 * np.exp(-lam * (r - R)) * (lam / r + 1 / r ** 2)
        / (kappa * (lam + 1 / R))
    )


def target(res, keya, alpha, keyb, beta, keyg, gamma_, fixed, blurb):
    import time
    t0 = time.time()

    phi_inf = alpha
    R = beta
    omega_b = gamma_

    p = Parameters(omega_b=omega_b, R=R, phi_inf=phi_inf, **_param_kwargs(fixed))

    # --- effective potential on the radial (end-to-end) coordinate ----------
    lam = np.sqrt(d2fb_dphi2(p.phi_inf, p) / p.kappa)
    args = (p.R, p.omega_b, p.kappa, lam, kuhnLength, beadSep)

    xs = np.linspace(p.R, p.L, 500)
    U_vals = U(xs, *args)

    # well = global min; barrier = first F sign change (- -> +) past the well
    i_well = int(np.argmin(U_vals))
    x_well, U_well = xs[i_well], U_vals[i_well]

    F_vals = F(xs, *args)
    after_well = xs > x_well
    F_after = F_vals[after_well]
    x_after = xs[after_well]
    sign_change = (np.sign(F_after[:-1]) < 0) & (np.sign(F_after[1:]) > 0)

    surface_params = dict(
        nb=p.nb, chi_b=p.chi_b, kappa=p.kappa, omega_b=p.omega_b,
        R=p.R, phi_inf=p.phi_inf, L=p.L, eps=p.eps,
    )
    rouse_params = dict(
        N=N, kuhn_length=kuhnLength, kBT=kBT, friction=friction,
        stiffness=stiffness, D=D, bead_separation=beadSep,
    )

    if not sign_change.any():
        # No adsorption well/barrier at this (low) concentration: the global
        # minimum is the entropic basin and there is nothing to escape from.
        # Record the potential + params and bail out so the sweep CPU stays
        # alive and the JSON is well-formed.
        res.update(
            blurb=blurb, alpha=float(phi_inf), R=float(R), omega_b=float(omega_b),
            barrier_found=False, lam=float(lam),
            surface_params=surface_params, rouse_params=rouse_params,
            total_potential=U_vals, pot_grid=xs,
            runtime_min=(time.time() - t0) / 60.0,
        )
        return res

    i_cross = int(np.argmax(sign_change))
    x_lo, x_hi = x_after[i_cross], x_after[i_cross + 1]
    F_lo, F_hi = F_after[i_cross], F_after[i_cross + 1]
    x_bar = x_lo - F_lo * (x_hi - x_lo) / (F_hi - F_lo)

    x_bar_a = p.R + 0.05      # region A boundary (inside the adsorption well)
    x_bar_b = x_bar     # region B boundary (fully escaped into the basin)
    U_bar = float(U(x_bar, *args))
    dU = U_bar - U_well
    omega_bar = np.sqrt(-d2U(x_bar, *args))

    # entropic (outer) basin minimum -- replaces the notebook's hardcoded r=8
    outer = xs > x_bar
    x_bot = float(xs[outer][int(np.argmin(U_vals[outer]))])
    omega_bot = np.sqrt(d2U(x_bot, *args))
    dU_bot = float(U(x_bar, *args) - U(x_bot, *args))

    # --- time step, relaxation and analytical Kramers times -----------------
    dV_R = -F(p.R, *args)
    relax_well = gamma * kB * T_bath / dV_R ** 2
    relax_bar = gamma / abs(d2U(x_bar, *args))
    dt_drift = 0.05 * min(relax_well, relax_bar)
    noise_step_target = 0.1 * (x_bar - x_well)
    dt_noise = noise_step_target ** 2 / (2 * kB * T_bath / gamma)
    dt = float(min(dt_drift, dt_noise))

    tau_K = float(
        (gamma / dV_R) * (np.sqrt(2 * np.pi * kB * T_bath) / omega_bar)
        * np.exp(dU / (kB * T_bath))
        * (1 - np.exp(-dV_R * (x_bar - p.R) / (kB * T_bath)))
    )
    tau_bot = float(
        (2 * np.pi * gamma / (omega_bot * omega_bar)) * np.exp(dU_bot / (kB * T_bath))
    )

    # tau = tau_mult * tau_K
    tau = 5_000_000 * dt
    capped = tau / dt > MAX_STEPS
    if capped:
        tau = MAX_STEPS * dt

    # --- simulate -----------------------------------------------------------
    sim = OverdampedHarmonicLangevin(
        tau     = tau, 
        dim     = 1, 
        N       = n_walkers, 
        dt      = dt,
        T       = T_bath, 
        gamma   = gamma, 
        kB      = kB,
        F       = lambda r: F(r, *args),
        r0      = p.R, 
        r_min   = p.R + 0.01,
        seed    = keya + 42, 
        nskip   = nskip,
        ab      = (x_bar_a, x_bar_b), 
        store   = "full",
    )
    sim.generate()

    # --- exit / re-entry sojourn statistics ---------------------------------
    # A = {r <= a} is the adsorption well: an A-sojourn closes on escape, so its
    # duration is an exit (well -> basin) time.  B = {r >= b} is the escaped
    # basin: a B-sojourn duration is a re-entry (basin -> well) time.  The first
    # interval of each walker (entry pinned to step 0) is dropped, as in the nb.
    exit_steps = np.concatenate([a[1:] for a in sim.A_intervals]) if sim.N else np.array([])
    reentry_steps = np.concatenate([b[1:] for b in sim.B_intervals]) if sim.N else np.array([])
    exit_times = exit_steps * dt        # intervals are high-res steps -> *dt
    reentry_times = reentry_steps * dt

    mean_exit = float(exit_times.mean()) if exit_times.size else float("nan")
    mean_reentry = float(reentry_times.mean()) if reentry_times.size else float("nan")
    k_exit = 1.0 / mean_exit if exit_times.size else float("nan")
    k_reentry = 1.0 / mean_reentry if reentry_times.size else float("nan")

    exit_counts, _ = np.histogram(exit_times, bins=exit_T_grid)
    reentry_counts, _ = np.histogram(reentry_times, bins=reentry_T_grid)

    # --- record -------------------------------------------------------------
    res.update(
        blurb               = blurb,
        alpha               = float(phi_inf),
        keya                = keya,
        R                   = float(R),
        omega_b             = float(omega_b),
        barrier_found       = True,
        lam                 = float(lam),
        surface_params      = surface_params,
        rouse_params        = rouse_params,

        langevin_params     = dict(
            T_bath          = T_bath, 
            gamma           = gamma, 
            kB              = kB, 
 	        dt              = dt, 
 	 	    tau             = tau,
            n_walkers       = n_walkers, 
 	 	    nskip           = nskip, 
 	 	    n_steps         = int(tau / dt),
            ns              = sim.ns, 
 	 	    seed            = keya + 42, 
 	 	    capped          = bool(capped),
        ),

        # potential and landmarks
        total_potential     = U_vals,
        pot_grid            = xs,
        x_well              = float(x_well), 
 	 	U_well              = float(U_well),
        x_bar               = float(x_bar), 
 	 	U_bar               = U_bar, 
 	 	dU                  = float(dU),
        x_bar_a             = float(x_bar_a), 
 	 	x_bar_b             = float(x_bar_b), 
 	 	omega_bar           = float(omega_bar),
        x_bot               = x_bot, 
 	 	omega_bot           = float(omega_bot), 
 	 	dU_bot              = dU_bot,
            
        # relaxation / analytical times
        relax_well          = float(relax_well), 
 	 	relax_bar           = float(relax_bar),
        tau_bot             = tau_bot, 
 	 	tau_K               = tau_K, 
 	 	k_kramers           = 1.0 / tau_K, 
 	 	dt                  = dt, 
            
        # empirical rates
        mean_exit_time      = mean_exit, 
 	 	mean_reentry_time   = mean_reentry,
        k_exit              = k_exit, 
 	 	k_reentry           = k_reentry,
        n_exit_events       = int(exit_times.size), 
 	 	n_reentry_events    = int(reentry_times.size),
            
        # distributions
        exit_T_grid         = exit_T_grid, 
 	 	exit_counts         = exit_counts,
        reentry_T_grid      = reentry_T_grid, 
 	 	reentry_counts      = reentry_counts,
        exit_times          = exit_times, 
 	 	reentry_times       = reentry_times,
            
        # walker-0 trajectory (light-storage mode)
        distance_traj       = sim.r[:, 0, :].flatten(),
        sample_skip         = nskip, 
        trace_dt            = dt * nskip,
        runtime_min         = (time.time() - t0) / 60.0,
    )
    return res
