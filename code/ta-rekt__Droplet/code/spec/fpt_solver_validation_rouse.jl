# Stress-test of the 1D Langevin/Smoluchowski reduction of Rouse cyclization dynamics.
# Sweep-spec recast of viz_util/fpt_solver_validation_rouse.jl (sim/sweep.jl driver +
# target! contract), so the comparison can be run in parallel over many parameters:
#
#   ./runJL sim/sweep.jl spec/fpt_solver_validation_rouse.jl <ncpus> data/YYYYMMDD/<name>
#
# Three independent estimates of the SAME two first-passage-time distributions on the
# end-to-end (bead-to-bead) distance x are compared at every swept parameter:
#
#   1. backward-FPE solver (src/fpt_solver.jl) on the reduced 1D coordinate,
#   2. direct 1D Euler-Maruyama Langevin in the same reduced potential U_tot,
#   3. the FULL Rouse chain (src/rouse.jl), whose two tracked ends feel only the
#      adsorption force -- its spring and entropy emerge from the chain itself.
#
# U_tot already contains the Rouse entropy (-2 log x is the radial Jacobian, plus the
# Gaussian spring 3x^2/2sa^2), so 1 and 2 are the exact radial Smoluchowski dynamics of
# that coordinate WITH the chain slaved to it.  Method 3 keeps the internal modes, so the
# gap between 3 and 1/2 is exactly the error of the Markovian 1D reduction.
#
# The two distributions are
#   exit:     x_contact -> x_ridge   (climb OUT of the adsorption well)
#   re-entry: x_ridge   -> x_contact (fall back INTO the well -- "time to reach x = 0"
#                                     in the surface-gap coordinate d = x - R)
# where x_ridge is the unstable equilibrium of U_tot nearest contact (its first interior
# local maximum) and x_contact = R + contact_tol is the bound-state threshold.  All three
# methods use the SAME pair of thresholds, so the comparison is definitional rather than
# merely qualitative: for the Rouse chain they are exactly the (a, b) Schmitt-trigger
# bounds of the sojourn machinery, whose A-intervals are then the exit times and whose
# B-intervals are the re-entry times.
#
# contact_tol > 0 is required, not cosmetic: the tracked Rouse pair has a hard reflecting
# wall at r_min = R, so the chain never attains x = R exactly and the "zero gap" threshold
# has to sit a hair above the wall.  The same threshold is then imposed on the solver and
# the 1D Langevin.
#
# Goodness of fit: two-sided one-sample Kolmogorov-Smirnov of each sampled set (Langevin,
# Rouse) against the solver CDF, evaluated exactly at any t from the modal decomposition
# S(t) = sum_k a_k exp(-lambda_k t) returned by solve_fpt (no histogram binning, no
# interpolation).  D is the headline number; the p-value uses the asymptotic Kolmogorov
# distribution and should be read as indicative only, since sojourns cut from one long
# trajectory are not perfectly independent.

using Droplet
using Random
using Statistics
using Printf

# ── Sweep contract ───────────────────────────────────────────────────────────────────
slist           = []        # off-lattice study: the driver's lattice species list is unused
checkpoint_path = ""

# ── Physics held fixed across the sweep ──────────────────────────────────────────────
kBT        = 1.0
friction   = 1.0
kuhnLength = 1.0                    # a
D_bead     = kBT / friction
D_eff      = 2 * D_bead             # relative coordinate of the two tracked beads

# Surface-adsorption model (linear-response condensation profile on a sphere).
nb, chi_b, kappa, phi_inf, eps, R = 1.0, 2.5, 1.0, 0.1, 1e-6, 5.0

# ── Sweep axes ───────────────────────────────────────────────────────────────────────
# alphas = beadSep s, the number of bonds between the tracked beads.  This is the axis the
#          stress test lives on: the subchain relaxation time grows like s^2 while the
#          barrier crossing time does not, so the Markovian reduction must fail at large s.
#          It is also the Rouse cost axis (~s work per step), hence the parallelised one.
# betas  = omega_b, the adsorption strength (sets the barrier height).
alphas  = [100]
betas   = [-0.25]
nparams = length(alphas)

# ── Numerics ─────────────────────────────────────────────────────────────────────────
contact_tol  = 0.05         # x_contact = R + contact_tol (bound-state / "zero gap" threshold)
Ngrid        = 2000         # FPE grid points per problem
ridge_scan_n = 8000         # fine scan used to locate the ridge
x_max_mult   = 3.0          # outer (reflecting) end of the re-entry domain: R + mult*a*sqrt(s)
dt_bd        = 2e-4         # 1D Langevin step: the well is steep, EM biases the MFPT otherwise
n_langevin   = 2_000        # independent 1D Langevin first-passage samples per problem
maxstep_mult = 100          # per-sample step cap, in units of the solver MFPT
dt_rouse     = 5e-4
nsteps_rouse = 2_000_000
trace_skip   = 200          # subsampling of the stored Rouse distance trace
seed0        = 20260728

blurb = "Stress-testing the reduction of Rouse cyclization dynamics to a 1D Langevin
         equation in the effective end-to-end potential.  At each (beadSep, omega_b) the
         exit time (contact -> ridge) and the re-entry time (ridge -> contact) are
         estimated three ways -- backward-FPE solver, 1D Langevin in the reduced
         potential, and the full Rouse chain -- and the sampled CDFs are compared to the
         solver CDF with a Kolmogorov-Smirnov test."

# ── Helpers ──────────────────────────────────────────────────────────────────────────

# 1D overdamped Euler-Maruyama first-passage sampler.  `force(x) = -dU/dx`; the walker
# starts at x0, reflects off xrefl, and is absorbed on crossing xabs.  The direction is
# inferred from the geometry (xabs above xrefl => escaping upward).  The absorbing
# threshold carries the standard 0.5826*sqrt(2 D dt) Mannella/Milstein shift toward the
# reflecting side, which removes the leading O(sqrt(dt)) bias of a hard discrete cutoff.
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

# Survival S(t) = sum_k a_k exp(-lambda_k t) from a solve_fpt result: exact at any t, so
# the reference CDF needs no interpolation off the solver's output time grid.
function fpt_survival(rates, weights, t)
    s = 0.0
    @inbounds for k in eachindex(rates)
        s += weights[k] * exp(-rates[k] * t)
    end
    return s
end

# Two-sided one-sample Kolmogorov-Smirnov statistic of `samples` against the solver CDF,
# D = sup_t |F_n(t) - F(t)|, with the asymptotic p-value Q_KS(lambda) = 2 sum (-1)^(j-1)
# exp(-2 j^2 lambda^2) and the Stephens finite-n correction on lambda.
function ks_vs_solver(samples, rates, weights)
    n = length(samples)
    n == 0 && return Dict{String,Any}("D" => NaN, "p" => NaN, "n" => 0)
    xs = sort(samples)
    D  = 0.0
    @inbounds for i in 1:n
        F = clamp(1 - fpt_survival(rates, weights, xs[i]), 0.0, 1.0)
        D = max(D, max(i / n - F, F - (i - 1) / n))
    end
    en = sqrt(n)
    λ  = (en + 0.12 + 0.11 / en) * D
    q  = 0.0
    for j in 1:100
        q += 2 * (-1)^(j - 1) * exp(-2 * j^2 * λ^2)
    end
    return Dict{String,Any}("D" => D, "p" => clamp(q, 0.0, 1.0), "n" => n)
end

# ── Target ───────────────────────────────────────────────────────────────────────────
function target!(res, args...)
    slist, keya, alpha, beta, blurb = args
    t0 = time()

    beadSep = Int(alpha)
    omega_b = Float64(beta)
    seed    = seed0 + 7919 * keya + round(Int, 1e4 * abs(omega_b))

    p   = Parameters(nb = nb, chi_b = chi_b, omega_b = omega_b, kappa = kappa,
                     phi_inf = phi_inf, eps = eps, R = R)
    lam = sqrt(d2fb_dphi2(p.phi_inf, p) / p.kappa)
    h0  = p.omega_b

    U_ads(x) = -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
                (p.kappa * (lam + 1 / p.R) * x)
    U_tot(x) = U_ads(x) - 2 * log(x) + 3 * x^2 / (2 * beadSep * kuhnLength^2)

    # Forces by central difference off the single potential definition, so the reduced
    # dynamics and the Rouse pair force can never drift apart from U_tot.
    h_fd = 1e-6
    U_tot_force(x) = -(U_tot(x + h_fd) - U_tot(x - h_fd)) / (2 * h_fd)   # -dU_tot/dx
    U_ads_force(x) = -(U_ads(x + h_fd) - U_ads(x - h_fd)) / (2 * h_fd)   # -dU_ads/dx

    surface_params = Dict("nb" => p.nb, "chi_b" => p.chi_b, "kappa" => p.kappa,
                          "omega_b" => p.omega_b, "R" => p.R, "L" => p.L,
                          "phi_inf" => p.phi_inf, "lam" => lam)
    rouse_params   = Dict("kuhn_length" => kuhnLength, "kBT" => kBT, "friction" => friction,
                          "D_bead" => D_bead, "D_eff" => D_eff, "bead_separation" => beadSep,
                          "nbeads" => beadSep + 1, "dt" => dt_rouse, "nsteps" => nsteps_rouse,
                          "t_total" => nsteps_rouse * dt_rouse)

    # ── Landmarks: contact threshold, ridge, outer end of the re-entry domain ─────────
    x_contact = p.R + contact_tol
    x_max     = p.R + x_max_mult * kuhnLength * sqrt(beadSep)
    x_scan    = collect(range(p.R, x_max; length = ridge_scan_n))
    U_scan    = U_tot.(x_scan)
    ir        = findfirst(i -> U_scan[i] > U_scan[i-1] && U_scan[i] > U_scan[i+1],
                          2:length(x_scan)-1)

    res["blurb"]          = blurb
    res["alpha"]          = alpha
    res["beta"]           = beta
    res["keya"]           = keya
    res["seed"]           = seed
    res["surface_params"] = surface_params
    res["rouse_params"]   = rouse_params
    res["pot_grid"]       = x_scan[1:8:end]
    res["total_potential"] = U_scan[1:8:end]
    res["ads_potential"]  = U_ads.(x_scan[1:8:end])
    res["x_contact"]      = x_contact
    res["x_max"]          = x_max
    res["contact_tol"]    = contact_tol

    if ir === nothing
        # No interior local maximum: U_tot rises monotonically out of the adsorption well,
        # so there is no unstable equilibrium and neither first-passage problem is defined.
        # Record the landscape and bail out so the CPU stays alive and the JSON is well formed.
        @printf("s=%d omega_b=%+.3f | NO RIDGE (U_tot monotonic) -- skipped\n", beadSep, omega_b)
        res["barrier_found"] = false
        res["runtime_min"]   = (time() - t0) / 60.0
        return res
    end

    x_ridge = x_scan[ir + 1]
    barrier = U_tot(x_ridge) - U_tot(x_contact)
    res["barrier_found"] = true
    res["x_ridge"]       = x_ridge
    res["barrier"]       = barrier
    res["well_depth"]    = U_tot(x_contact)

    # ── 1. Backward-FPE solver ───────────────────────────────────────────────────────
    t_solver = time()
    # exit: [R, x_ridge], reflecting at the contact wall, absorbing at the ridge, started
    # from a delta at x_contact (the same point the Rouse A-sojourn clock starts from).
    xg_exit = collect(range(p.R, x_ridge; length = Ngrid))
    p0_exit = zeros(Ngrid); p0_exit[argmin(abs.(xg_exit .- x_contact))] = 1.0
    sol_exit = solve_fpt(FPTProblem(xg_exit, U_tot.(xg_exit); D = D_eff, kBT = kBT,
                                    left_bc = :reflecting, right_bc = :absorbing, p0 = p0_exit))
    # re-entry: [x_contact, x_max], absorbing at contact, reflecting far out, started from
    # a delta at the ridge.
    xg_ret = collect(range(x_contact, x_max; length = Ngrid))
    p0_ret = zeros(Ngrid); p0_ret[argmin(abs.(xg_ret .- x_ridge))] = 1.0
    sol_ret = solve_fpt(FPTProblem(xg_ret, U_tot.(xg_ret); D = D_eff, kBT = kBT,
                                   left_bc = :absorbing, right_bc = :reflecting, p0 = p0_ret))
    runtime_solver = time() - t_solver

    # ── 2. 1D Langevin in the reduced potential ──────────────────────────────────────
    t_lang = time()
    bd_exit, bd_exit_cens = langevin_fpt(U_tot_force, x_contact, p.R, x_ridge;
                                         D = D_eff, kBT = kBT, dt = dt_bd,
                                         nsamples = n_langevin,
                                         maxstep = ceil(Int, maxstep_mult * sol_exit.mfpt / dt_bd),
                                         seed = seed + 1)
    bd_ret, bd_ret_cens   = langevin_fpt(U_tot_force, x_ridge, x_max, x_contact;
                                         D = D_eff, kBT = kBT, dt = dt_bd,
                                         nsamples = n_langevin,
                                         maxstep = ceil(Int, maxstep_mult * sol_ret.mfpt / dt_bd),
                                         seed = seed + 2)
    runtime_langevin = time() - t_lang

    # ── 3. Full Rouse chain ──────────────────────────────────────────────────────────
    # a = x_contact and b = x_ridge are the same two thresholds used above: an A-sojourn
    # (entered at the down-crossing of a, closed at the up-crossing of b) is an exit time,
    # a B-sojourn is a re-entry time.  Only the adsorption force acts between the ends --
    # the spring and the -2 log x entropy are generated by the chain.
    t_rouse = time()
    r = RousePolymer(; boltzmannConst = 1.0, temperature = kBT, friction = friction,
                     kuhnLength = kuhnLength, r_min = p.R,
                     nbeads = Int32(beadSep + 1), nsteps = Int32(nsteps_rouse), dt = dt_rouse,
                     trackedBead1 = 1, trackedBead2 = 1 + beadSep,
                     attractiveForce = U_ads_force, saveMode = :distance_only,
                     a = x_contact, b = x_ridge, seed = seed + 3)
    init!(r)
    simulate!(r)
    runtime_rouse = time() - t_rouse

    # The first interval of each region is dropped: whichever region the chain starts in
    # has its clock pinned to step 0 rather than to a genuine threshold crossing.
    rouse_exit = Float64.(length(r.A_intervals) > 1 ? r.A_intervals[2:end] : Int[]) .* dt_rouse
    rouse_ret  = Float64.(length(r.B_intervals) > 1 ? r.B_intervals[2:end] : Int[]) .* dt_rouse

    # ── Goodness of fit vs the solver CDF ────────────────────────────────────────────
    ks = Dict{String,Any}(
        "langevin_exit"    => ks_vs_solver(bd_exit,    sol_exit.rates, sol_exit.weights),
        "langevin_reentry" => ks_vs_solver(bd_ret,     sol_ret.rates,  sol_ret.weights),
        "rouse_exit"       => ks_vs_solver(rouse_exit, sol_exit.rates, sol_exit.weights),
        "rouse_reentry"    => ks_vs_solver(rouse_ret,  sol_ret.rates,  sol_ret.weights))

    _mean(v) = isempty(v) ? NaN : mean(v)

    res["ks"] = ks
    res["mfpt"] = Dict(
        "solver_exit"      => sol_exit.mfpt,   "solver_reentry"   => sol_ret.mfpt,
        "langevin_exit"    => _mean(bd_exit),  "langevin_reentry" => _mean(bd_ret),
        "rouse_exit"       => _mean(rouse_exit), "rouse_reentry"  => _mean(rouse_ret))
    res["kramers_rate_exit"]    = sol_exit.rates[1]
    res["kramers_rate_reentry"] = sol_ret.rates[1]

    # solver curves (log-spaced auto grids, one per problem)
    res["solver_exit_t"]        = sol_exit.t
    res["solver_exit_pdf"]      = sol_exit.pdf
    res["solver_exit_survival"] = sol_exit.survival
    res["solver_ret_t"]         = sol_ret.t
    res["solver_ret_pdf"]       = sol_ret.pdf
    res["solver_ret_survival"]  = sol_ret.survival
    res["x_grid_exit"]          = xg_exit[1:8:end]
    res["U_tot_exit"]           = U_tot.(xg_exit[1:8:end])
    res["x_grid_ret"]           = xg_ret[1:8:end]
    res["U_tot_ret"]            = U_tot.(xg_ret[1:8:end])

    # raw first-passage samples (the notebook builds ECDFs from these)
    res["langevin_exit_times"]    = bd_exit
    res["langevin_reentry_times"] = bd_ret
    res["rouse_exit_times"]       = rouse_exit
    res["rouse_reentry_times"]    = rouse_ret
    res["langevin_censored"]      = Dict("exit" => bd_exit_cens, "reentry" => bd_ret_cens)
    res["langevin_params"]        = Dict("dt" => dt_bd, "nsamples" => n_langevin,
                                         "maxstep_mult" => maxstep_mult)

    res["distance_traj"]   = r.distanceTraj === nothing ? Float64[] :
                             r.distanceTraj[1:trace_skip:end]
    res["sample_skip"]     = trace_skip
    res["trace_dt"]        = dt_rouse * trace_skip
    res["bead_relax_time"] = r.timeScale
    res["bead_relax_scale"] = r.posScale
    # Rouse time of the tracked subchain, tau_R = timeScale * s^2 / pi^2: the reduction is
    # expected to hold only while this is small compared to the first-passage times.
    res["rouse_time"]      = r.timeScale * beadSep^2 / π^2

    res["runtime_solver_s"]   = runtime_solver
    res["runtime_langevin_s"] = runtime_langevin
    res["runtime_rouse_s"]    = runtime_rouse
    res["runtime_total_s"]    = time() - t0
    res["runtime_min"]        = (time() - t0) / 60.0

    @printf("s=%d omega_b=%+.3f | ridge=%.3f barrier=%.2f kBT | MFPT exit %.3g/%.3g/%.3g  ret %.3g/%.3g/%.3g (solver/BD/Rouse)\n",
            beadSep, omega_b, x_ridge, barrier,
            sol_exit.mfpt, _mean(bd_exit), _mean(rouse_exit),
            sol_ret.mfpt, _mean(bd_ret), _mean(rouse_ret))
    @printf("            KS-D  exit: BD %.3f (p=%.3g, n=%d)  Rouse %.3f (p=%.3g, n=%d)\n",
            ks["langevin_exit"]["D"], ks["langevin_exit"]["p"], ks["langevin_exit"]["n"],
            ks["rouse_exit"]["D"], ks["rouse_exit"]["p"], ks["rouse_exit"]["n"])
    @printf("            KS-D  ret : BD %.3f (p=%.3g, n=%d)  Rouse %.3f (p=%.3g, n=%d)  [%.1f min]\n",
            ks["langevin_reentry"]["D"], ks["langevin_reentry"]["p"], ks["langevin_reentry"]["n"],
            ks["rouse_reentry"]["D"], ks["rouse_reentry"]["p"], ks["rouse_reentry"]["n"],
            res["runtime_min"])

    return res
end
