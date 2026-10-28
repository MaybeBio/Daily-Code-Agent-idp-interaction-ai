# test_bvp_solver.jl — solve_isolated_sphere (src/bvp_solver.jl) checked against
# the Cahn "graphical construction" for surface adsorption / prewetting, in the
# zero-curvature limit R -> infinity (here R = 1e4).
#
# THE CONSTRUCTION.  In planar geometry the isolated-surface BVP is
#
#     kappa * phi'' = dfb/dphi(phi) - mu_inf/nb,      mu_inf = nb*dfb_dphi(phi_inf)
#
# with phi -> phi_inf, phi' -> 0 in the far field.  Multiplying by phi' and
# integrating in from infinity gives the first integral
#
#     (kappa/2) * phi'(z)^2 = Delta_f(phi(z)),
#     Delta_f(phi) = fb(phi) - fb(phi_inf) - dfb_dphi(phi_inf)*(phi - phi_inf)
#
# (Delta_f is fb measured from the tangent at phi_inf, so it vanishes quadratically
# there and is >= 0 for any phi_inf outside the two-phase region).  The Robin
# adsorption BC is -kappa*phi'(0) + omega_b = 0, i.e. phi'(0) = h0/kappa, so the
# surface composition is fixed by a single scalar equation -- the horizontal line
# cutting the Delta_f curve that the notebook draws:
#
#     Delta_f(phi_s) = h0^2 / (2*kappa).                                     (*)
#
# phi_s is the SMALLEST root of (*) above phi_inf: that is the branch reached by
# integrating in from the dilute far field, and the one the solver's continuation
# tracks.  If (*) has no root at all then Delta_f < h0^2/(2 kappa) everywhere,
# phi' never vanishes, the profile runs into phi = 1 at finite z, and the surface
# is completely wet -- the solver must then either saturate or fail, never report
# a converged finite film.
#
# bc_reg BREAKS (*).  With bc_reg > 0 the surface BC carries the regularizing
# log term (see _surf_slope), phi'(0) = (h0 + bc_reg*log(phi_s/(1-phi_s)))/kappa,
# and the construction becomes
#
#     Delta_f(phi_s) = (h0 + bc_reg*log(phi_s/(1-phi_s)))^2 / (2*kappa).     (**)
#
# Below phi_s = 1/2 the log term is negative and so reinforces an attractive
# h0 < 0: a bc_reg > 0 solve sits at a systematically HIGHER phi_s than the bare
# h0^2/(2 kappa) construction predicts.  Comparing the two is a mismatch of
# boundary conditions, not a solver error.  "BcRegShiftsTheConstruction" below
# pins both forms down.
#
# FREE ENERGY.  fb here is the ASYMMETRIC (polymer-in-solvent) Flory-Huggins form
# of src/bvp_solver.jl, whose critical point is chi_c = (1 + 1/sqrt(nb))^2/2 and
# whose binodal must be found by common tangent.  It is NOT the symmetric blend
# that spec/.../param.jl's `dilute_binodal` assumes (that one is the binodal of
# an effective chi = nb*chi_b and agrees only at nb = 1), so this file computes
# its own binodal.

using Test
using Random
using Printf
using LinearAlgebra
using Droplet

# ---------------------------------------------------------------------------
# The graphical construction, evaluated exactly
# ---------------------------------------------------------------------------

# Bisection on a bracketed sign change; 80 halvings exhaust Float64.
function _bisect(f, a, b; n = 80)
    fa = f(a)
    for _ in 1:n
        m = 0.5 * (a + b)
        (fa * f(m) <= 0) ? (b = m) : (a = m; fa = f(m))
    end
    return 0.5 * (a + b)
end

# fb measured from its tangent at phi_inf -- the curve the construction reads off.
deltaf(phi, p::Parameters) =
    fb(phi, p) - fb(p.phi_inf, p) - dfb_dphi(p.phi_inf, p) * (phi - p.phi_inf)

chi_critical(nb) = (1 + 1 / sqrt(nb))^2 / 2

# Common-tangent binodal (phi_dilute, phi_dense) of the asymmetric fb.  dfb_dphi
# rises, falls between the two spinodal points, then rises again, so each branch
# holds exactly one root of dfb_dphi = mu; bisect on mu until the two branches
# have equal grand potential fb - mu*phi.
function binodal(chi, nb)
    p = Parameters(nb = nb, chi_b = chi, omega_b = 0.0, kappa = 1.0,
                   phi_inf = 0.1, L = 1.0, R = Inf, eps = 1e-16)
    lo_sp, hi_sp = compute_spinodal(p)
    lo_sp === nothing && return (nothing, nothing)
    dilute(mu) = _bisect(x -> dfb_dphi(x, p) - mu, 1e-300, lo_sp)
    dense(mu)  = _bisect(x -> dfb_dphi(x, p) - mu, hi_sp, 1 - 1e-16)
    gap(mu)    = (fb(dense(mu), p) - mu * dense(mu)) -
                 (fb(dilute(mu), p) - mu * dilute(mu))
    mu = _bisect(gap, dfb_dphi(hi_sp, p), dfb_dphi(lo_sp, p))
    return dilute(mu), dense(mu)
end

# Delta_f has at most two interior turning points -- the other two solutions of
# dfb_dphi(phi) = dfb_dphi(phi_inf), one in each spinodal bracket -- so it is
# monotone on three segments: rising to the (barrier) local max, falling to the
# (well) local min on the dense branch, rising again into phi = 1.
function turning_points(p::Parameters)
    lo_sp, hi_sp = compute_spinodal(p)
    slope = dfb_dphi(p.phi_inf, p)
    return (_bisect(x -> dfb_dphi(x, p) - slope, lo_sp, hi_sp),
            _bisect(x -> dfb_dphi(x, p) - slope, hi_sp, 1 - 1e-16))
end

# Smallest root of Delta_f(phi) = E(phi) above phi_inf, or NaN when there is none
# (complete wetting).  E is the level line: the constant h0^2/(2 kappa) of (*),
# or the phi_s-dependent right-hand side of (**) when bc_reg > 0.
function surface_root(p::Parameters, E)
    barrier, well = turning_points(p)
    g(x) = deltaf(x, p) - (E isa Function ? E(x) : E)
    for (a, b) in ((p.phi_inf * (1 + 1e-12), barrier), (barrier, well),
                   (well, 1 - 1e-16))
        g(a) * g(b) <= 0 && return _bisect(g, a, b)
    end
    return NaN
end

# ---------------------------------------------------------------------------
# Both branches of the construction, and their free energies (bc_reg = 0)
# ---------------------------------------------------------------------------
# Inside the prewetting window the h0^2/(2 kappa) line cuts Delta_f once on each
# monotone segment: THIN (rising to the barrier), MIDDLE (falling into the
# well), THICK (rising out of the well).  The window is therefore
#
#     Delta_f(well) < h0^2/(2 kappa) < Delta_f(barrier),
#
# whose two edges are the folds where the thick resp. thin branch ends.  Only
# meaningful where Delta_f has its two turning points (phi_inf close enough to
# the binodal); every case below is chosen so.
function planar_roots(p::Parameters)
    barrier, well = turning_points(p)
    g(x) = deltaf(x, p) - p.omega_b^2 / (2 * p.kappa)
    root(a, b) = g(a) * g(b) <= 0 ? _bisect(g, a, b) : NaN
    return root(p.phi_inf * (1 + 1e-12), barrier), root(barrier, well),
           root(well, 1 - 1e-16)
end

function planar_window(p::Parameters)
    barrier, well = turning_points(p)
    return -sqrt(2 * p.kappa * deltaf(well, p)), -sqrt(2 * p.kappa * deltaf(barrier, p))
end

# Cahn's excess free energy per unit area of the planar film with surface
# composition phi_s.  On the profile (kappa/2) phi'^2 = Delta_f, so the bulk and
# gradient terms are equal and Int (Delta_f + kappa/2 phi'^2) dz becomes an
# integral over phi:
#
#     gamma(phi_s) = Int_{phi_inf}^{phi_s} [ sqrt(2 kappa Delta_f(phi)) + h0 ] dphi.
#
# d gamma/d phi_s vanishes exactly on (*); the thin and thick roots are minima,
# the middle one a maximum, and gamma_thick - gamma_thin is Cahn's equal-area
# integral between the curve sqrt(2 kappa Delta_f) and the line |h0|.  Composite
# Simpson on the closed interval: the integrand is smooth there (sqrt(Delta_f) is
# linear in phi - phi_inf at the lower end, and Delta_f > 0 in the well below
# the binodal).  This is the planar limit of excess_free_energy / (4 pi R^2).
function planar_gamma(p::Parameters, phi_s::Real; n::Int = 20000)
    f(x) = sqrt(2 * p.kappa * max(deltaf(x, p), 0.0)) + p.omega_b
    h = (phi_s - p.phi_inf) / n
    s = f(p.phi_inf) + f(phi_s)
    for i in 1:n - 1
        s += (isodd(i) ? 4 : 2) * f(p.phi_inf + i * h)
    end
    return s * h / 3
end

# Planar prewetting field: the h0 in the window where gamma_thin = gamma_thick.
function prewetting_field(p::Parameters)
    lo, hi = planar_window(p)
    function dgamma(h0)
        q = with(p; omega_b = h0)
        thin, _, thick = planar_roots(q)
        return planar_gamma(q, thick) - planar_gamma(q, thin)
    end
    return _bisect(dgamma, lo + 1e-6 * (hi - lo), hi - 1e-6 * (hi - lo); n = 60)
end

# ---------------------------------------------------------------------------
# Solver domain
# ---------------------------------------------------------------------------
# The truncation z in [0, L] must be long enough that the Dirichlet wall at z = L
# does no work: too short a box props the profile up and the solver happily
# "converges" to a finite film that the construction says cannot exist.  The
# first integral gives the exact profile width by quadrature,
# dz = dphi / sqrt(2*Delta_f(phi)/kappa) (cubic-clustered near the logarithmically
# divergent phi_inf end, cut off at phi_inf*(1+1e-3)), and the far tail decays over
# ell = sqrt(kappa/fb''(phi_inf)); L = 2*width + 20*ell covers both with margin.
# Complete-wetting cases have no phi_s, so their width runs all the way to phi = 1.
function _film_width(p::Parameters, phi_t::Real; n = 1200)
    lo = p.phi_inf * (1 + 1e-3)
    w, x_prev = 0.0, lo
    for k in 1:n
        x = lo + (phi_t - lo) * (k / n)^3
        d = max(deltaf(0.5 * (x + x_prev), p), 0.0)
        d > 0 && (w += (x - x_prev) / sqrt(2 * d / p.kappa))
        x_prev = x
    end
    return w
end

solver_domain(p::Parameters, phi_s::Real) =
    2 * _film_width(p, isnan(phi_s) ? 1 - 1e-6 : phi_s) +
    12 * sqrt(p.kappa / d2fb_dphi2(p.phi_inf, p))

# ---------------------------------------------------------------------------
# Branch tracking (solve_isolated_sphere, track_branches = true), planar limit
# ---------------------------------------------------------------------------
# Prewetting cases as (nb, chi_b, phi_inf); phi_inf = 0 means 0.9 x the dilute
# binodal.  Each has a planar bistable window (planar_window).
const PREWETTING_CASES = ((1.0, 2.5, 0.13), (1.0, 2.5, 0.14), (1.0, 2.8, 0.08),
                          (5.0, 1.5, 0.0))

# Planar-limit sphere at field h0 (R = 1e8: curvature corrections ~1e-8, below
# everything tested).
function prewetting_case(nb, chi, phi_inf, h0)
    phi_b, _ = binodal(chi, nb)
    return Parameters(nb = nb, chi_b = chi, omega_b = h0, kappa = 1.0,
                      phi_inf = phi_inf == 0 ? 0.9 * phi_b : phi_inf,
                      L = 1.0, R = 1e8, eps = 1e-16)
end

# Branch-tracked solve, the box sized on the thickest film the construction
# allows so both branches fit.  The box depends on p only, not on N, so solves
# at N and 2N differ by dz alone (needed for the Richardson checks).  The thick
# candidate is seeded by the solver itself (_dense_well) unless phi_dense is given.
function solve_tracked(p::Parameters, N::Int; phi_dense = nothing)
    thin, _, thick = planar_roots(p)
    L = solver_domain(p, isnan(thick) ? thin : thick)
    return solve_isolated_sphere(with(p; L = L); L_iso = L, N = N, bc_reg = 0.0,
                                 track_branches = true, phi_dense = phi_dense)
end

# Surface composition of one branch of a bistable solve, whichever was returned.
branch_phis(s, b::Symbol) = s.branch == b ? s.phi[1] : s.phi_other[1]

# ---------------------------------------------------------------------------

@testset "bvp_solver" begin

    # -- 1. The first-integral identity on a few hand-picked points -----------
    # Weak, moderate and near-prewetting adsorption at the pipeline's own
    # chi_b = 2.5, nb = 1: phi_s must sit on Delta_f = h0^2/(2 kappa) to the
    # 1/R accuracy of the R = 1e4 sphere.
    @testset "FirstIntegral" begin
        for (h0, phi_inf) in ((-0.2, 0.05), (-0.5, 0.05), (-0.5, 0.09),
                              (-0.5, 0.12), (-0.7, 0.03))
            p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = h0, kappa = 1.0,
                           phi_inf = phi_inf, L = 1.0, R = 1e4, eps = 1e-16)
            E = h0^2 / (2 * p.kappa)
            r = surface_root(p, E)
            L = solver_domain(p, r)
            s = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 1500,
                                      bc_reg = 0.0)
            @test s.success
            @test isapprox(s.phi[1], r; atol = 1e-3)
            @test isapprox(deltaf(s.phi[1], p), E; atol = 1e-3)
            @test s.phi[end] ≈ phi_inf                 # Dirichlet far field
            @test all(diff(s.phi) .<= 1e-9)            # monotone decay
        end
    end

    # -- 2. bc_reg moves the construction, it does not break the solver -------
    # The regression this file was written for: spec/jobs/20260804/run1 solves
    # with bc_reg = 5e-2 and was compared against the bare h0^2/(2 kappa) line.
    # At bc_reg = 0 the solver lands on (*); at bc_reg > 0 it lands on (**), and
    # the two differ by far more than any numerical tolerance.
    @testset "BcRegShiftsTheConstruction" begin
        bc_reg = 5e-2
        for phi_inf in (0.05, 0.065, 0.0786)
            p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.5, kappa = 1.0,
                           phi_inf = phi_inf, L = 1.0, R = 1e4, eps = 1e-16)
            E = p.omega_b^2 / (2 * p.kappa)
            r = surface_root(p, E)
            L = solver_domain(p, r)

            bare = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 1500,
                                         bc_reg = 0.0)
            @test bare.success
            @test isapprox(bare.phi[1], r; atol = 1e-3)

            # (**): Delta_f(phi_s) = (h0 + bc_reg*log(phi_s/(1-phi_s)))^2/(2 kappa).
            reg = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 1500,
                                        bc_reg = bc_reg)
            r_reg = surface_root(p, x -> (p.omega_b +
                                          bc_reg * log(x / (1 - x)))^2 / (2 * p.kappa))
            @test reg.success
            @test isapprox(reg.phi[1], r_reg; atol = 1e-3)

            # ... and the shift is large.  These phi_s all sit below 1/2, where the
            # log term is negative and so reinforces the attractive h0: the
            # regularized surface is markedly WETTER than the bare construction.
            @test reg.phi[1] > bare.phi[1] + 0.02
        end
    end

    # -- 3. Deep robustness scan against the construction ---------------------
    # 1000 random (omega_b, chi_b, nb) x 10 far-field compositions spanning
    # 0.1-0.8 of the dilute binodal (equivalently 10 values of mu_inf, which is
    # monotone in phi_inf below the binodal).
    #
    # chi_b is drawn from [1.5, 2.5] EXCEPT that a binodal only exists above
    # chi_critical(nb); chi_critical = 2.0 at nb = 1 (and 1.05 / 0.87 at nb =
    # 5 / 10), so at nb = 1 the draw is from [2.04, 2.5] -- below that there is
    # no two-phase region and "a fraction of the binodal" is undefined.
    #
    # Every sample is checked.  Where (*) has a root the surface composition must
    # match it; where it has none the construction predicts complete wetting, and
    # the solver must not answer with a converged finite film.  Note that most of
    # the requested omega_b range is past the wetting transition: at 1000 trials
    # 4185 of the 10000 samples have a finite phi_s, 5769 are completely wet.
    #
    # ~7 min at the default 1000 trials; BVP_SCAN_TRIALS shortens it.  The seed is
    # fixed and smaller runs are a prefix of the same stream, so the tolerances
    # below (set with 2-3x headroom over the measured 1000-trial distribution)
    # stay valid at any trial count.
    @testset "ConstructionScan" begin
        ntrial = parse(Int, get(ENV, "BVP_SCAN_TRIALS", "1000"))
        nphi   = 10
        Random.seed!(20260804)

        dphi_err = Float64[]      # |phi_s(solver) - phi_s(construction)|
        energy_err = Float64[]    # |Delta_f(phi_s) - h0^2/(2 kappa)|
        n_root = n_stall = n_wet = n_wet_ok = n_false_film = 0

        for _ in 1:ntrial
            nb   = rand((1.0, 5.0, 10.0))
            cmin = max(1.5, 1.02 * chi_critical(nb))
            chi  = cmin + rand() * (2.5 - cmin)
            h0   = -2.0 * rand()
            phi_b, _ = binodal(chi, nb)

            for frac in range(0.1, 0.8; length = nphi)
                p = Parameters(nb = nb, chi_b = chi, omega_b = h0, kappa = 1.0,
                               phi_inf = frac * phi_b, L = 1.0, R = 1e4,
                               eps = 1e-16)
                E = h0^2 / (2 * p.kappa)
                r = surface_root(p, E)
                L = solver_domain(p, r)
                s = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 600,
                                          maxiter = 100, bc_reg = 0.0)

                if isnan(r)                       # construction: complete wetting
                    n_wet += 1
                    s.phi[1] > 0.99 && (n_wet_ok += 1)
                    (s.success && s.phi[1] < 0.99) && (n_false_film += 1)
                elseif !s.success                 # Newton stalled; counted below
                    n_stall += 1
                else
                    n_root += 1
                    push!(dphi_err, abs(s.phi[1] - r))
                    push!(energy_err, abs(deltaf(s.phi[1], p) - E))
                end
            end
        end

        sort!(dphi_err); sort!(energy_err)
        q(v, f) = v[clamp(round(Int, f * length(v)), 1, length(v))]
        @printf("\n  ConstructionScan: %d samples | finite phi_s %d (stalled %d) | complete wetting %d\n",
                ntrial * nphi, n_root, n_stall, n_wet)
        @printf("    |dphi_s|          median %.2e  p90 %.2e  p99 %.2e  max %.2e\n",
                q(dphi_err, 0.5), q(dphi_err, 0.9), q(dphi_err, 0.99), dphi_err[end])
        @printf("    |Delta_f - h0^2/2k| median %.2e  p99 %.2e  max %.2e\n",
                q(energy_err, 0.5), q(energy_err, 0.99), energy_err[end])

        # Surface composition matches the construction on every sample (measured:
        # median 1.5e-5, p99 1.4e-3, max 3.0e-3).  The residual is the O(1/R)
        # curvature the R = 1e4 sphere still carries -- it falls off as 1/R, and at
        # R = 1e8 or in the exact planar slab the worst case drops to 3e-6.  It is
        # amplified where the h0^2/(2 kappa) line grazes the top of the Delta_f
        # barrier and the root is ill-conditioned, hence the loose ceiling against
        # a median two orders of magnitude tighter.  For scale, mistaking (**) for
        # (*) -- the bc_reg confusion of testset 2 -- costs 0.03 to 0.2.
        @test dphi_err[end] <= 1e-2
        @test q(dphi_err, 0.99) <= 3e-3
        @test q(dphi_err, 0.5) <= 1e-4

        # The same statement read off the well-conditioned axis: whatever phi_s
        # the solver returns, it sits on the Delta_f = h0^2/(2 kappa) level line.
        @test energy_err[end] <= 3e-2
        @test q(energy_err, 0.99) <= 1e-2

        # Complete wetting: the solver may saturate or stall, but it must never
        # report a converged finite adsorbed film where the construction has no
        # intersection at all.
        @test n_false_film == 0
        @test n_wet_ok >= 0.5 * n_wet

        # Newton is allowed to stall on the stiffest near-prewetting samples, but
        # must carry the overwhelming majority.
        @test n_stall <= 0.05 * (n_root + n_stall)
    end

    # -- 4. Both branches against the construction, to discretization order ---
    # Inside the planar window (20/50/80 % of the way from the thick fold to the
    # thin fold) the tracked solve must return BOTH films, each sitting on its
    # own root of (*) with Cahn's free energy gamma.  Stringency comes from
    # resolution: at fixed box, the error must fall as dz^2 between N = 1500 and
    # 3000 (ratio ~4), and the Richardson extrapolation (4 e_2N - e_N)/3 must
    # vanish to the 1e-8 level -- i.e. the solver converges to the graphical
    # construction itself, not merely near it.
    @testset "PlanarBranches" begin
        A = 4 * pi * 1e8^2
        raw_phi = Float64[]; raw_gam = Float64[]       # |error| at N = 3000
        ratios  = Float64[]                            # e(1500) / e(3000)
        ext_phi = Float64[]; ext_gam = Float64[]       # |Richardson error|
        for (nb, chi, phi_inf) in PREWETTING_CASES, frac in (0.2, 0.5, 0.8)
            p0 = prewetting_case(nb, chi, phi_inf, 0.0)
            lo, hi    = planar_window(p0)
            p         = with(p0; omega_b = lo + frac * (hi - lo))
            thin, _, thick = planar_roots(p)
            g_thin, g_thick = planar_gamma(p, thin), planar_gamma(p, thick)

            s1 = solve_tracked(p, 1500)
            s2 = solve_tracked(p, 3000)
            @test s1.success && s2.success
            @test s1.bistable && s2.bistable
            (s1.bistable && s2.bistable) || continue

            err(s) = (branch_phis(s, :thin) - thin, branch_phis(s, :thick) - thick,
                      s.F_thin / A - g_thin, s.F_thick / A - g_thick)
            e1, e2 = err(s1), err(s2)
            append!(raw_phi, abs.(e2[1:2])); append!(raw_gam, abs.(e2[3:4]))
            # The order test needs the error well clear of the ~1e-8 floor of the
            # comparison itself (quadrature of gamma, 1/R): below 1e-7 the ratio
            # measures that floor (e.g. 3.4 at nb = 5), not the dz^2 error.
            append!(ratios, [e1[k] / e2[k] for k in 1:4 if abs(e2[k]) > 1e-7])
            ext = (4 .* e2 .- e1) ./ 3
            append!(ext_phi, abs.(ext[1:2])); append!(ext_gam, abs.(ext[3:4]))

            # Selection is the equal-area rule: the lower-gamma film is returned,
            # the other one is kept alongside it.
            @test s2.branch == (g_thin < g_thick ? :thin : :thick)
            @test s2.F == min(s2.F_thin, s2.F_thick)
            @test s2.phi[1] == branch_phis(s2, s2.branch)
            @test length(s2.phi_other) == length(s2.phi)
        end
        @printf("\n  PlanarBranches: |dphi_s| %.1e  |dgamma| %.1e at N=3000 | e(1500)/e(3000) in [%.2f, %.2f] | Richardson |dphi_s| %.1e |dgamma| %.1e\n",
                maximum(raw_phi), maximum(raw_gam), minimum(ratios), maximum(ratios),
                maximum(ext_phi), maximum(ext_gam))
        @test maximum(raw_phi) <= 3e-5
        @test maximum(raw_gam) <= 1e-5
        @test all(3.5 .<= ratios .<= 4.5)              # second-order convergence
        @test maximum(ext_phi) <= 1e-7
        @test maximum(ext_gam) <= 5e-8
    end

    # -- 5. The prewetting point --------------------------------------------
    # The planar prewetting field h_pw solves gamma_thin = gamma_thick (Cahn's
    # equal area).  Just on either side of it the tracked solve must switch from
    # the thin to the thick film, and linear interpolation of its own
    # F_thick - F_thin across that +-1e-4 bracket must put the crossing on h_pw:
    # second-order convergent, and to ~1e-7 after Richardson extrapolation.  (The
    # bracket is that narrow because linear interpolation of the curved
    # F_thick - F_thin is biased by O(bracket^2) independently of dz, which
    # Richardson cannot remove: +-1e-3 leaves a 5e-7 floor.)
    @testset "PlanarPrewettingPoint" begin
        for (nb, chi, phi_inf) in PREWETTING_CASES
            p0 = prewetting_case(nb, chi, phi_inf, 0.0)
            h_pw  = prewetting_field(p0)
            h_w, h_s = h_pw * (1 - 1e-4), h_pw * (1 + 1e-4)   # weaker / stronger
            rel = map((1500, 3000)) do N
                s_w = solve_tracked(with(p0; omega_b = h_w), N)
                s_s = solve_tracked(with(p0; omega_b = h_s), N)
                @test s_w.bistable && s_s.bistable
                @test s_w.branch == :thin
                @test s_s.branch == :thick
                d_w, d_s = s_w.F_thick - s_w.F_thin, s_s.F_thick - s_s.F_thin
                (h_w - d_w * (h_s - h_w) / (d_s - d_w)) / h_pw - 1
            end
            ext = (4 * rel[2] - rel[1]) / 3
            @printf("  PlanarPrewettingPoint (nb=%g chi=%g phi_inf=%.4f): h_pw = %.6f | rel. error N=1500 %.1e  N=3000 %.1e  Richardson %.1e\n",
                    nb, chi, p0.phi_inf, h_pw, rel[1], rel[2], ext)
            @test abs(rel[2]) <= 2e-5
            @test abs(rel[2]) < abs(rel[1]) / 3
            @test abs(ext) <= 2e-7
        end
    end

    # -- 6. The whole window, and just outside it -----------------------------
    # Hysteresis must be detected across the ENTIRE planar window (39 interior
    # fields, the outermost 1/40 of the window from each fold), F_thick - F_thin
    # must change sign exactly once -- at h_pw -- and the returned phi_s must
    # jump exactly there.  Just past each fold only one branch survives: no
    # bistability, and the lone film sits on the surviving root.
    @testset "PlanarWindowCoverage" begin
        A = 4 * pi * 1e8^2
        for (nb, chi, phi_inf) in PREWETTING_CASES
            p0 = prewetting_case(nb, chi, phi_inf, 0.0)
            lo, hi = planar_window(p0)
            h_pw   = prewetting_field(p0)
            hs     = collect(range(lo, hi; length = 41))[2:end - 1]   # weak -> strong
            sols   = [solve_tracked(with(p0; omega_b = h), 1500) for h in hs]
            @test all(s -> s.success && s.bistable, sols)

            flips = findall(diff(sign.([s.F_thick - s.F_thin for s in sols])) .!= 0)
            @test length(flips) == 1
            @test hs[flips[1]] > h_pw > hs[flips[1] + 1]
            @test findall(diff([s.phi[1] for s in sols]) .> 0.1) == flips

            for (h, survivor) in ((lo - 0.05 * (hi - lo), 1), (hi + 0.05 * (hi - lo), 3))
                p = with(p0; omega_b = h)
                r = planar_roots(p)
                isnan(r[survivor]) && continue      # past the thick fold: complete wetting
                @test count(isnan, r) == 2          # the construction has one root
                s = solve_tracked(p, 1500)
                @test s.success && !s.bistable && s.branch == :unique
                @test isapprox(s.phi[1], r[survivor]; atol = 1e-4)
                @test isapprox(s.F / A, planar_gamma(p, r[survivor]); atol = 5e-5)
            end
        end
    end

    # -- 7. Envelope identity, per branch, with and without bc_reg ------------
    # F is stationary in phi, so dF/d(omega_b) is the explicit derivative of the
    # surface term: 4 pi R^2 (phi_s - phi_inf), on EACH branch separately.  It
    # holds only if F carries the surface free energy the BC was derived from --
    # at bc_reg > 0 that includes bc_reg*[h(phi_s) - h(phi_inf)], h(x) = x log x
    # + (1-x) log(1-x).  The negative control drops that term and must fail.
    # Finite R = 50, where the window is narrower but survives bc_reg = 0.05.
    @testset "EnvelopeIdentity" begin
        hs(x) = x * log(x) + (1 - x) * log(1 - x)
        dev = Float64[]; dev_noent = Float64[]
        for (bc_reg, h0s) in ((0.0, (-0.23, -0.27)), (0.05, (-0.28, -0.30)))
            for h0 in h0s
                p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = h0, kappa = 1.0,
                               phi_inf = 0.14, L = 20.0, R = 50.0, eps = 1e-16)
                A = 4 * pi * p.R^2
                δ = 1e-4
                s, sp, sm = (solve_isolated_sphere(with(p; omega_b = h0 + d); L_iso = 20.0,
                                                   N = 3000, bc_reg = bc_reg,
                                                   track_branches = true)
                             for d in (0.0, δ, -δ))
                @test s.bistable && sp.bistable && sm.bistable
                (s.bistable && sp.bistable && sm.bistable) || continue
                for b in (:thin, :thick)
                    F(x)  = b == :thin ? x.F_thin : x.F_thick
                    phi_s = branch_phis(s, b)
                    rhs   = A * (phi_s - p.phi_inf)
                    push!(dev, (F(sp) - F(sm)) / (2δ) / rhs - 1)
                    noent(x) = F(x) - A * bc_reg * (hs(branch_phis(x, b)) - hs(p.phi_inf))
                    bc_reg > 0 && push!(dev_noent, (noent(sp) - noent(sm)) / (2δ) / rhs - 1)
                end
            end
        end
        @printf("  EnvelopeIdentity: max |ratio - 1| = %.1e (with surface entropy), %.2f without\n",
                maximum(abs, dev), maximum(abs, dev_noent))
        @test maximum(abs, dev) <= 2e-4
        @test maximum(abs, dev_noent) >= 0.05
    end

    # -- 8. Default path (track_branches = false) ----------------------------
    @testset "DefaultPath" begin
        for bc_reg in (0.0, 0.05)
            p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.3, kappa = 1.0,
                           phi_inf = 0.14, L = 10.0, R = 5.0, eps = 1e-10)
            s = solve_isolated_sphere(p; L_iso = 10.0, N = 400, bc_reg = bc_reg)
            @test s.success
            @test !s.bistable && s.branch == :unique
            @test isnan(s.F_thin) && isnan(s.F_thick) && isempty(s.phi_other)
            @test s.F == excess_free_energy(s.phi, s.dphi, s.z, p; bc_reg = bc_reg)
        end
    end

    # -- 8b. The thick-film seed --------------------------------------------
    # The solver seeds the thick candidate itself with the well of Delta_f (the
    # dense root of fb'(phi) = fb'(phi_inf)), which must coincide with the
    # construction's well.  Where Delta_f has no well the construction has a
    # single root, so tracking must reduce exactly to the default solve.  An
    # explicit phi_dense (here the dense binodal) only changes the starting
    # guess, never the branches found.
    @testset "ThickSeed" begin
        for (nb, chi, phi_inf) in PREWETTING_CASES
            p = prewetting_case(nb, chi, phi_inf, -0.3)
            @test isapprox(Droplet._dense_well(p), turning_points(p)[2]; atol = 1e-12)
        end

        # No well (phi_inf far below the binodal): one root, no thick candidate.
        p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.3, kappa = 1.0,
                       phi_inf = 0.05, L = 10.0, R = 5.0, eps = 1e-10)
        @test Droplet._dense_well(p) === nothing
        s_def = solve_isolated_sphere(p; L_iso = 10.0, N = 400)
        s_trk = solve_isolated_sphere(p; L_iso = 10.0, N = 400, track_branches = true)
        @test s_trk.phi == s_def.phi
        @test !s_trk.bistable && s_trk.branch == :unique

        # Explicit seed: the same two branches as the automatic one.
        p0     = prewetting_case(1.0, 2.8, 0.08, 0.0)
        lo, hi = planar_window(p0)
        p      = with(p0; omega_b = 0.5 * (lo + hi))
        s_auto = solve_tracked(p, 1500)
        s_bin  = solve_tracked(p, 1500; phi_dense = binodal(2.8, 1.0)[2])
        @test s_auto.bistable && s_bin.bistable
        @test s_auto.branch == s_bin.branch
        @test isapprox(s_auto.F_thin, s_bin.F_thin; rtol = 1e-8)
        @test isapprox(s_auto.F_thick, s_bin.F_thick; rtol = 1e-8)
        @test maximum(abs, s_auto.phi .- s_bin.phi) < 1e-6

        # Scale-free seed: the production regime of spec/local/20261004/large_R
        # (kappa = 10, nb = 2230, R = 1e4), whose thick film is ~15 wide with an
        # interface ~25 wide.  The default reach = 2 alone found bistability at 1
        # of these 5 consecutive phi_inf (the finite-R window, shifted from the
        # planar one by curvature); the doubling reach ladder finds all 5, and
        # none just outside.
        q = Parameters(nb = 2230.0, chi_b = 0.544821635760637,
                       omega_b = -0.03576257306287501, kappa = 10.12228759025409,
                       phi_inf = 2e-5, L = 200.0, R = 1e4, eps = 1e-8)
        window = (2.0568646565168776e-5, 2.0726246194131617e-5, 2.0883845823094452e-5,
                  2.104144545205729e-5, 2.1199045081020126e-5)
        outside = (2.041104693620594e-5, 2.1356644709982964e-5)
        bi(x) = solve_isolated_sphere(with(q; phi_inf = x); L_iso = 200.0, N = 400,
                                      f_tol = 1e-6, track_branches = true).bistable
        @test all(bi, window)
        @test !any(bi, outside)
    end

    # -- 9. Morse index = number of negative eigenvalues of -J ----------------
    # _morse_index counts negative pivots of the unpivoted elimination of -J,
    # relying on Sylvester's law of inertia (J is a row-weighted symmetric
    # matrix).  Checked against a dense eigensolver on random tridiagonals with
    # positive off-diagonal products (the structure of every Jacobian here), and
    # on the actual isolated-sphere stencil with an f'' dip deep enough to make
    # 0, 1, 2, ... directions unstable.
    @testset "MorseIndex" begin
        Random.seed!(20261004)
        for _ in 1:300
            n  = rand(2:60)
            lo = rand(n - 1) .+ 0.05                 # J[i+1, i]
            up = rand(n - 1) .+ 0.05                 # J[i, i+1]
            dd = 3.0 .* randn(n)
            λ  = eigvals(Matrix(Tridiagonal(-lo, -dd, -up)))
            @test maximum(abs, imag.(λ)) < 1e-8
            @test Droplet._morse_index(vcat(0.0, lo), dd, vcat(up, 0.0)) ==
                  count(<(0), real.(λ))
        end

        # Sphere stencil of solve_isolated_sphere's jacobian! (bc_reg = 0).
        R, kappa, n, dz = 5.0, 1.0, 400, 0.02
        z  = (0:n - 1) .* dz
        counts = Int[]
        for depth in range(0.0, 60.0; length = 31)
            f2 = 0.5 .- depth .* exp.(-((z .- 2.0) ./ 0.7) .^ 2)   # fb'' with a dip
            dl = [i == 1 ? 0.0 : kappa * (1 / dz^2 - (1 / (R + z[i])) / dz) for i in 1:n]
            du = [i == 1 ? 2kappa / dz^2 :
                  i == n ? 0.0 : kappa * (1 / dz^2 + (1 / (R + z[i])) / dz) for i in 1:n]
            dd = -2kappa / dz^2 .- f2
            λ  = real.(eigvals(Matrix(Tridiagonal(-dl[2:n], -dd, -du[1:n - 1]))))
            minimum(abs, λ) > 1e-6 || continue      # skip a dip sitting on a fold
            k = Droplet._morse_index(dl, dd, du)
            @test k == count(<(0), λ)
            push!(counts, k)
        end
        @test maximum(counts) >= 3                   # the dip did reach unstable modes
    end
end
