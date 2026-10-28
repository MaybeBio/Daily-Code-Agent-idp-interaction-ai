# First-passage time between the surface-adsorption well and the entropic ridge,
# in the combined attractive-adsorption + Rouse-entropy potential.  Three
# independent estimates of the SAME distribution are compared:
#
#   1. backward-FPE solver (src/fpt_solver.jl) on the reduced 1D coordinate,
#   2. direct 1D Euler–Maruyama Langevin in the SAME reduced potential U_tot,
#   3. the FULL Rouse polymer (src/rouse.jl): a 151-bead chain whose two ends
#      feel the adsorption force, with times read from the sojourn machinery.
#
# `mode` picks which of the two first-passage problems is simulated:
#
#   :exit    x_contact -> x_ridge    climb OUT of the adsorption well.  Reflecting
#                                    at the contact wall, absorbing at the ridge;
#                                    Rouse A-sojourns (bound-state dwells).
#   :entry   x_ridge -> x_contact    fall back INTO the well, i.e. the time to
#                                    return to zero surface gap d = x - R.
#                                    Absorbing at contact, reflecting far out at
#                                    x_max; Rouse B-sojourns.
#
# The reduced coordinate x is the end-to-end (bead-to-bead) distance.  U_tot
# already contains the Rouse entropy (-2 log x is the radial Jacobian, plus the
# Gaussian spring), so the 1D solver with a flat operator is the exact radial
# Smoluchowski equation for that coordinate.  Its diffusion coefficient is the
# relative-coordinate value D = 2 kBT/gamma (two end beads).  The full Rouse
# chain applies only the adsorption force U_ads between the ends — its spring and
# entropy emerge from the chain, so this tests the 1D reduction against the
# genuine many-bead dynamics.
#
# The contact threshold is x_contact = R + contact_tol rather than R itself: the
# tracked Rouse pair has a hard reflecting wall at r_min = R, so the chain never
# attains x = R exactly and the "zero gap" threshold has to sit a hair above the
# wall.  The same threshold is then imposed on the solver and on the 1D Langevin,
# so all three methods solve the identical first-passage problem — and the Rouse
# (a, b) Schmitt-trigger bounds are exactly (x_contact, x_ridge) in both modes.
#
# This script only SIMULATES: it writes every sample and curve to a JSON, and
# viz_util/fpt_solver_validation_rouse_plot.jl draws the figure from those files,
# so the ~10 min of sampling is never repeated for a plotting change.
#
#   julia --project=@. viz_util/fpt_solver_validation_rouse.jl          # :exit
#   julia --project=@. viz_util/fpt_solver_validation_rouse.jl entry    # :entry
#
# Writes media/figures/<today>/fpt_solver_validation_rouse_<mode>.json.

using Droplet
using Random
using Statistics
using Printf
using Dates
import JSON

# --- Which first-passage problem to simulate ---------------------------------
mode = isempty(ARGS) ? :exit : Symbol(ARGS[1])       # :exit or :entry
@assert mode in (:exit, :entry) "mode must be :exit or :entry, got :$mode"

# --- Combined potential ------------------------------------------------------
p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.25, kappa = 1.0,
               phi_inf = 0.1, eps = 1e-6, R = 5.0)     # omega_b weakened to ~7 kBT
d2f_phi_inf = d2fb_dphi2(p.phi_inf, p)
h0  = p.omega_b
lam = sqrt(d2f_phi_inf / p.kappa)
kuhnLength = 1.0
beadSep    = 150

U_ads(x) = -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
            (p.kappa * (lam + 1 / p.R) * x)
U_tot(x) = U_ads(x) - 2 * log(x) + 3 * x^2 / (2 * beadSep * kuhnLength^2)

h_fd = 1e-6
U_tot_prime(x) = -(U_tot(x + h_fd) - U_tot(x - h_fd)) / (2 * h_fd)
U_ads_force(x) = -(U_ads(x + h_fd) - U_ads(x - h_fd)) / (2 * h_fd)  # -dU_ads/dr

# --- Physical constants (shared by all three methods) ------------------------
kBT      = 1.0
friction = 1.0
D_bead   = kBT / friction
D_eff    = 2 * D_bead          # relative coordinate of the two end beads

# --- Landmarks: contact threshold, ridge, outer reflecting wall --------------
contact_tol = 0.05
x_contact   = p.R + contact_tol
x_max       = p.R + 3 * kuhnLength * sqrt(beadSep)   # far end of the re-entry domain
xs = collect(range(p.R, x_max; length = 8000))
Us = U_tot.(xs)
ridge_idx = findfirst(i -> Us[i] > Us[i-1] && Us[i] > Us[i+1], 2:length(xs)-1)
ridge_idx === nothing && error("U_tot has no interior maximum on [R, x_max]: no unstable " *
                               "equilibrium, so neither first-passage problem is defined " *
                               "(the entropic basin has to sit outside the adsorption range)")
x_ridge = xs[ridge_idx + 1]
barrier = U_tot(x_ridge) - U_tot(x_contact)
@printf("mode=%s  lam=%.3f  contact x=%.3f  ridge x=%.3f  barrier=%.2f kBT\n",
        mode, lam, x_contact, x_ridge, barrier)

t_start_run = time()

# --- 1. Backward-FPE solver (delta IC at the start point) --------------------
# exit : domain [R, x_ridge]        — reflecting contact wall, absorbing ridge
# entry: domain [x_contact, x_max]  — absorbing contact, reflecting far wall
Ngrid = 2000
if mode === :exit
    xg = collect(range(p.R, x_ridge; length = Ngrid))
    left_bc, right_bc      = :reflecting, :absorbing
    x_start, x_refl, x_abs = x_contact, p.R, x_ridge
else
    xg = collect(range(x_contact, x_max; length = Ngrid))
    left_bc, right_bc      = :absorbing, :reflecting
    x_start, x_refl, x_abs = x_ridge, x_max, x_contact
end
p0 = zeros(Ngrid); p0[argmin(abs.(xg .- x_start))] = 1.0
prob = FPTProblem(xg, U_tot.(xg); D = D_eff, kBT = kBT,
                  left_bc = left_bc, right_bc = right_bc, p0 = p0)
res  = solve_fpt(prob)

# --- 2. Direct 1D Langevin in U_tot ------------------------------------------
# `force` is -dU/dx (same convention as U_tot_prime and the Rouse pair force);
# the overdamped drift is D*β*force*dt.  The walker starts at x0, reflects off
# xrefl and is absorbed on crossing xabs; the direction is inferred from the
# geometry, and the absorbing threshold carries the standard 0.5826*sqrt(2 D dt)
# Mannella/Milstein shift toward the reflecting side.  Samples that exceed
# `maxstep` are censored (counted, not recorded) — the re-entry distribution has
# a long tail and a hard cap keeps one stuck walker from stalling the run.
function langevin_fpt(force, x0, xrefl, xabs; D, kBT, dt, nsamples, maxstep, seed)
    rng = MersenneTwister(seed)
    c1  = D * dt / kBT
    c2  = sqrt(2 * D * dt)
    up  = xabs > xrefl
    xa  = up ? xabs - 0.5826 * c2 : xabs + 0.5826 * c2
    times     = Float64[]
    ncensored = 0
    for _ in 1:nsamples
        x   = x0
        hit = false
        for step in 1:maxstep
            x += c1 * force(x) + c2 * randn(rng)
            if up
                x <= xrefl && (x = 2 * xrefl - x)
                x >= xa && (push!(times, step * dt); hit = true; break)
            else
                x >= xrefl && (x = 2 * xrefl - x)
                x <= xa && (push!(times, step * dt); hit = true; break)
            end
        end
        hit || (ncensored += 1)
    end
    return times, ncensored
end

dt_bd      = 2e-4       # small dt: the well is steep, EM biases the MFPT otherwise
n_langevin = 10_000
t_bd, bd_censored = langevin_fpt(U_tot_prime, x_start, x_refl, x_abs;
                                 D = D_eff, kBT = kBT, dt = dt_bd, nsamples = n_langevin,
                                 maxstep = ceil(Int, 100 * res.mfpt / dt_bd), seed = 7)

# --- 3. Full Rouse polymer (151 beads); times from the sojourn machinery -----
# a = x_contact and b = x_ridge are the same two thresholds used above: an
# A-sojourn (entered at the down-crossing of a, closed at the up-crossing of b)
# is an exit time, a B-sojourn is a re-entry time — one trajectory yields both,
# so `mode` only selects which list is read out.  The number of sojourns is set
# by t_total/(MFPT_exit + MFPT_entry), hence the long run.
dt_r    = 5e-4
nsteps  = 10_000_000                       # ~3 min at 151 beads; ~130 sojourns per region
r = RousePolymer(; boltzmannConst = 1.0, temperature = kBT, friction = friction,
                 kuhnLength = kuhnLength, r_min = p.R,
                 nbeads = Int32(beadSep + 1), nsteps = Int32(nsteps), dt = dt_r,
                 trackedBead1 = 1, trackedBead2 = 1 + beadSep,
                 attractiveForce = U_ads_force, saveMode = :sojourn_only,
                 a = x_contact, b = x_ridge, seed = 20260701)
init!(r); simulate!(r)
# The first interval of each region is dropped: whichever region the chain starts
# in has its clock pinned to step 0 rather than to a genuine threshold crossing.
intervals = mode === :exit ? r.A_intervals : r.B_intervals
t_rouse   = Float64.(length(intervals) > 1 ? intervals[2:end] : Int[]) .* dt_r

# The re-entry distribution is heavy-tailed (a walker that wanders up the
# entropic basin takes a long time back), so the median is the honest summary at
# these sample sizes — the sample mean is dominated by the largest few draws.
_mean(v)   = isempty(v) ? NaN : mean(v)
_median(v) = isempty(v) ? NaN : median(v)
runtime_min = (time() - t_start_run) / 60
@printf("MFPT    solver=%.3f   1D-Langevin=%.3f (n=%d, %d censored)   full-Rouse=%.3f (n=%d)\n",
        res.mfpt, _mean(t_bd), length(t_bd), bd_censored, _mean(t_rouse), length(t_rouse))
@printf("median  1D-Langevin=%.3f   full-Rouse=%.3f   [%.1f min]\n",
        _median(t_bd), _median(t_rouse), runtime_min)

# --- Dump everything the plotter needs ---------------------------------------
# Includes the potential curve itself, so the plotting script never re-derives
# the physics; summary statistics are left to the plotter (they are one line off
# the raw samples, and duplicating them here would let the two drift apart).
pot_x = collect(range(p.R, x_max; length = 1200))
outdir = joinpath("media", "figures", Dates.format(today(), "yyyymmdd"))
mkpath(outdir)
out = joinpath(outdir, "fpt_solver_validation_rouse_$(mode).json")
open(out, "w") do io
    JSON.print(io, Dict(
        "mode" => String(mode), "omega_b" => p.omega_b, "R" => p.R, "lam" => lam,
        "barrier" => barrier, "contact_tol" => contact_tol,
        "x_contact" => x_contact, "x_ridge" => x_ridge, "x_max" => x_max,
        "x_start" => x_start, "x_abs" => x_abs,
        "U_contact" => U_tot(x_contact), "U_ridge" => U_tot(x_ridge),
        "pot_x" => pot_x, "pot_U" => U_tot.(pot_x),
        "D_eff" => D_eff, "beadSep" => beadSep, "kuhn_length" => kuhnLength,
        "dt_bd" => dt_bd, "n_langevin" => n_langevin, "langevin_censored" => bd_censored,
        "dt_rouse" => dt_r, "nsteps_rouse" => nsteps,
        "mfpt_solver" => res.mfpt, "kramers_rate" => res.rates[1],
        "solver_t" => res.t, "solver_pdf" => res.pdf, "solver_survival" => res.survival,
        "t_bd" => t_bd, "t_rouse" => t_rouse, "runtime_min" => runtime_min))
end
println("saved -> ", out)
