# test_bvp_arclength.jl — pseudo-arclength continuation of the two-sphere BVP
# in phi_inf (src/bvp_arclength.jl).
#
# Four things are worth pinning down, in increasing cost:
#   1. the bordered linear algebra, including the singular-J case that IS a fold;
#   2. agreement with plain phi_inf continuation where the branch is a graph
#      over phi_inf (no fold), so the machinery cannot be silently wrong;
#   3. traversal of a real S-curve: two limit points, three branches, and the
#      free energy continuous with dF/ds -> 0 at each fold;
#   4. the envelope identity dF/d(omega_b) == Int (phi_s - phi_inf) dA, which
#      holds for any correct free energy at a stationary state and is the
#      cheapest guard on excess_free_energy_two_spheres.

using Test
using Random
using LinearAlgebra
using SparseArrays
using Droplet
using Droplet: _bordered_solve, _bordered_solve_full, _arc_dot, n_unknowns

# Int (phi_s - phi_inf) dA over both spheres; Simpson in theta, as in
# excess_free_energy_two_spheres.
function _surface_adsorption(g::Geometry, phi::AbstractMatrix)
    theta = collect(range(0.0, pi; length = 1001))
    h = theta[2] - theta[1]
    total = 0.0
    for zc in (g.z_c1, g.z_c2)
        phi_s = phi_on_sphere(phi, g, zc, theta)
        s = 0.0
        for k in eachindex(theta)
            w = (k == 1 || k == length(theta)) ? 1.0 : (iseven(k) ? 4.0 : 2.0)
            s += w * (phi_s[k] - g.p.phi_inf) * sin(theta[k])
        end
        total += s * h / 3
    end
    return 2pi * g.p.R^2 * total
end

# A cheap fold-free two-sphere setup (n ~ 2.5e3).
function _benign_setup()
    p = Parameters(chi_b = 1.0, omega_b = -0.2, kappa = 1.0, phi_inf = 0.10,
                   R = 3.0, eps = 1e-6)
    d, Nrho, dz, z_buf, rho_max = 1.0, 40, 0.2, 2.0, 8.0
    Nz = Int(round((4p.R + d + 2z_buf) / dz)) + 1
    return (p = p, d = d, Nrho = Nrho, Nz = Nz, z_buf = z_buf, rho_max = rho_max)
end

# A cheap setup that DOES fold: near-saturated bulk, gap narrow enough that the
# two adsorption films cannot coexist with a dilute neck (n ~ 7e3, ~5 s).
function _folding_setup()
    p = Parameters(nb = 1.0, chi_b = 2.8, omega_b = -0.5, kappa = 1.0,
                   phi_inf = 0.055, L = 5.0, R = 5.0, eps = 1e-8)
    d, Nrho, dz, rho_max = 1.95, 50, 0.15, 12.0
    z_buf = ceil(Int, (p.R + 3.0) / dz) * dz - p.R
    Nz = Int(round((4p.R + d + 2z_buf) / dz)) + 1
    return (p = p, d = d, Nrho = Nrho, Nz = Nz, z_buf = z_buf, rho_max = rho_max)
end

@testset "bvp_arclength" begin

@testset "Bordered solve" begin
    # Well-conditioned: Keller bordering and the full assembly must agree.
    rng = MersenneTwister(12345)
    n = 60
    J = sparse(5.0I, n, n) .+ sprandn(rng, n, n, 0.05)
    Gl = collect(range(-1.0, 1.0; length = n))
    tx = fill(1.0 / n, n)
    tl = 0.7
    r  = collect(range(0.3, -0.2; length = n))
    rn = 0.05

    dx1, dl1 = _bordered_solve(J, Gl, tx, tl, r, rn)
    dx2, dl2 = _bordered_solve_full(J, Gl, tx, tl, r, rn)
    @test dl1 ≈ dl2 atol = 1e-9
    @test dx1 ≈ dx2 atol = 1e-9
    # and it really solves the bordered system
    @test maximum(abs, J * dx1 .+ Gl .* dl1 .- r) < 1e-9
    @test abs(dot(tx, dx1) + tl * dl1 - rn) < 1e-9

    # AT a fold: G_x is exactly singular.  The toy G(x,lambda) = x^2 - lambda has
    # G_x = 2x = 0 and G_lam = -1 at the limit point x = 0, lambda = 0.  Keller
    # bordering must detect the collapsing pivot and fall through to the full
    # assembly, which is nonsingular there.
    Js  = sparse(reshape([0.0], 1, 1))
    Gls = [-1.0]
    txs = [1.0]
    tls = 0.0
    rs  = [0.25]
    rns = 0.5
    dxs, dls = _bordered_solve(Js, Gls, txs, tls, rs, rns)
    @test all(isfinite, dxs) && isfinite(dls)
    @test abs(Js[1, 1] * dxs[1] + Gls[1] * dls - rs[1]) < 1e-10
    @test abs(txs[1] * dxs[1] + tls * dls - rns) < 1e-10
end

@testset "Scaled inner product" begin
    a = [3.0, 4.0]; b = [1.0, 0.0]
    @test _arc_dot(a, 2.0, b, 5.0) ≈ 3.0 / 2 + 10.0
    @test _arc_dot(a, 0.0, a, 0.0) ≈ dot(a, a) / 2
end

@testset "polish_at_phi_inf" begin
    s = _benign_setup()
    geom = build_geometry(s.p, s.d; rho_max = s.rho_max, z_buf = s.z_buf,
                          Nrho = s.Nrho, Nz = s.Nz)
    res = solve_two_spheres(s.p, s.d; rho_max = s.rho_max, z_buf = s.z_buf,
                            Nrho = s.Nrho, Nz = s.Nz, omega_continuation = false)
    @test res[:success]

    lam = 0.13
    x, ok = polish_at_phi_inf(geom, res[:x], lam)
    @test ok
    # it must leave the geometry AT lambda, or a following free-energy call
    # would use the wrong bulk reference
    @test geom.p.phi_inf == lam
    # ... and land on the same solution a fresh solve finds
    ref = solve_two_spheres(Droplet.with(s.p; phi_inf = lam), s.d;
                            rho_max = s.rho_max, z_buf = s.z_buf,
                            Nrho = s.Nrho, Nz = s.Nz, omega_continuation = false)
    @test ref[:success]
    @test maximum(abs, x .- ref[:x]) < 1e-6
end

@testset "Agrees with natural continuation (no fold)" begin
    s = _benign_setup()
    geom = build_geometry(s.p, s.d; rho_max = s.rho_max, z_buf = s.z_buf,
                          Nrho = s.Nrho, Nz = s.Nz)
    seed = solve_two_spheres(s.p, s.d; rho_max = s.rho_max, z_buf = s.z_buf,
                             Nrho = s.Nrho, Nz = s.Nz, omega_continuation = false)
    targets = collect(range(0.10, 0.20; length = 6))

    work = (Matrix{Float64}(undef, geom.Nrho, geom.Nz),
            Matrix{Float64}(undef, geom.Nrho, geom.Nz))
    got = Dict{Float64, Float64}()
    function on_hit(x, lam, seg, stable)
        xp, ok = polish_at_phi_inf(geom, x, lam; work = work)
        ok && (got[lam] = gap_center_phi(geom, scatter(xp, geom)))
    end
    trace = arclength_phi_inf(geom, seed[:x], 0.10, 0.20;
                              targets = targets, on_hit = on_hit)
    @test trace.status == :target
    @test isempty(trace.folds)              # this branch is a graph over phi_inf
    @test all(trace.segment .== 0)
    @test all(>=(0.0), trace.tangent_lambda)   # phi_inf increases monotonically

    # Every target reached, and matching a plain warm-started phi_inf march.
    x = seed[:x]
    for t in targets
        @test haskey(got, t)
        ref = solve_two_spheres(Droplet.with(s.p; phi_inf = t), s.d;
                                rho_max = s.rho_max, z_buf = s.z_buf,
                                Nrho = s.Nrho, Nz = s.Nz, x0_init = x,
                                omega_continuation = false)
        @test ref[:success]
        x = ref[:x]
        @test got[t] ≈ gap_center_phi(ref[:geom], ref[:phi]) atol = 1e-6
    end

    # The corrector should not be limping: a smooth branch is 1-2 Newton steps
    # per continuation point.  (A mismatch between the bordering row and the
    # scaled arclength residual silently costs ~10x here.)
    @test sum(trace.newton_iters) / max(length(trace.lambda) - 1, 1) < 4.0
end

@testset "Traverses the capillary-bridging S-curve" begin
    s = _folding_setup()
    geom = build_geometry(s.p, s.d; rho_max = s.rho_max, z_buf = s.z_buf,
                          Nrho = s.Nrho, Nz = s.Nz, bc_reg = 0.0)
    seed = solve_two_spheres(s.p, s.d; rho_max = s.rho_max, z_buf = s.z_buf,
                             Nrho = s.Nrho, Nz = s.Nz, bc_reg = 0.0,
                             omega_continuation = false)
    @test seed[:success]

    frac = cell_fluid_fraction(geom)
    lam = Float64[]; Fv = Float64[]; sv = Float64[]; gc = Float64[]; seg = Int[]
    function on_point(x, l, sarc, sg, stable)
        phi = scatter(x, geom)
        push!(lam, l); push!(sv, sarc); push!(seg, sg)
        push!(gc, gap_center_phi(geom, phi))
        push!(Fv, excess_free_energy_two_spheres(
                    Dict{Symbol, Any}(:geom => geom, :phi => phi); frac = frac))
    end
    trace = arclength_phi_inf(geom, seed[:x], 0.055, 0.093;
                              on_point = on_point, dlambda0 = 6e-4,
                              max_steps = 500)

    # Two limit points, and the sweep comes out the far side on branch 2.
    @test length(trace.folds) >= 2
    @test maximum(trace.segment) >= 2
    @test trace.status == :target
    # The thin branch folds ABOVE the bridged one in phi_inf: that overlap is
    # the hysteresis window where all three states coexist.
    @test trace.folds[1] > trace.folds[2]
    # Dilute gap through to a fully condensed bridge.
    @test minimum(gc) < 0.2
    @test maximum(gc) > 0.85

    # on_point is index-aligned with the trace, so a fold sits exactly at the
    # index where the segment counter increments (equivalently where
    # dlambda/ds changes sign).  Do NOT look it up by nearest lambda: lambda
    # turns around there, so several points share nearly the same value.
    @test length(lam) == length(trace.lambda)
    fold_idx = [k for k in 3:(length(trace.segment) - 2)
                if trace.segment[k] > trace.segment[k - 1]]
    @test length(fold_idx) >= 2
    span = maximum(Fv) - minimum(Fv)
    secant(j) = (Fv[j] - Fv[j - 1]) / (sv[j] - sv[j - 1])
    for k in fold_idx
        # The two branches meeting at a limit point ARE the same solution, so F
        # is continuous through it.
        @test abs(Fv[k] - Fv[k - 1]) < 0.05 * span
        # And dF/ds = (dF/dlambda)(dlambda/ds) passes through zero there,
        # because dlambda/ds does and dF/dlambda stays finite (envelope
        # theorem).  The turning point falls BETWEEN samples, so the secant
        # flips sign somewhere in the neighbouring window rather than at a
        # fixed offset from the segment increment.
        window = [secant(j) for j in (k - 1):(k + 1)]
        @test minimum(window) < 0.0 < maximum(window)
    end
end

@testset "Envelope identity for the free energy" begin
    # dF/d(omega_b) = Int phi_s dA at a stationary point (envelope theorem: the
    # dphi/d(omega_b) term dies against delta(Omega)/delta(phi) = 0).  The code's
    # F drops the constant omega_b * phi_inf * A, so the identity it must satisfy
    # is dF/d(omega_b) = Int (phi_s - phi_inf) dA.
    #
    # This holds for the DISCRETE problem only if the discretisation is
    # variational.  It converges to 1 under refinement at bc_reg = 0 (ratio
    # 0.992 / 0.997 / 0.998 at n = 2.6e3 / 6.5e3 / 1.5e4).
    #
    # It must hold at bc_reg > 0 as well.  That BC is
    #     kappa dphi/ds = omega_b + bc_reg * log(phi_s / (1 - phi_s)),
    # whose surface free energy carries an extra
    #     bc_reg * Int [phi_s log phi_s + (1 - phi_s) log(1 - phi_s)] dA.
    # bc_reg still drops out of d/d(omega_b), so the identity is unchanged --
    # which makes it the guard on that term being present in
    # excess_free_energy_two_spheres.  Before it was added the ratio fell to
    # ~0.87 (thin film) / ~0.81 (bridged) at bc_reg = 0.05, branch-dependently,
    # and free energies at bc_reg > 0 were not comparable across branches.
    s = _benign_setup()
    sol(w, reg) = solve_two_spheres(Droplet.with(s.p; omega_b = w), s.d;
                                    rho_max = s.rho_max, z_buf = s.z_buf,
                                    Nrho = s.Nrho, Nz = s.Nz, bc_reg = reg,
                                    omega_continuation = false)
    h = 0.005
    for reg in (0.0, 0.05)
        r_lo, r_0, r_hi = sol(s.p.omega_b - h, reg), sol(s.p.omega_b, reg),
                          sol(s.p.omega_b + h, reg)
        @test r_lo[:success] && r_0[:success] && r_hi[:success]

        dF = (excess_free_energy_two_spheres(r_hi) -
              excess_free_energy_two_spheres(r_lo)) / (2h)
        ads = _surface_adsorption(r_0[:geom], r_0[:phi])
        @test dF / ads ≈ 1.0 atol = 0.02
    end
end

end  # @testset "bvp_arclength"
