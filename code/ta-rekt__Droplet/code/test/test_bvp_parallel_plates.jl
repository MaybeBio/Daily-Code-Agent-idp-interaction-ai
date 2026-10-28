# test_bvp_parallel_plates.jl — Tests for the planar slab BVP in src/bvp_solver.jl
# (solve_parallel_plates, psi_linear_plates, excess_free_energy,
# plates_pair_potential), the port of the planar branches of
# viz_util/bvp_solver.py and of viz_util/bvp_parallel_plates.py.
#
# Coverage:
#   1. Single plate, weak adsorption: phi(0) vs the linearised prediction
#      phi_inf - omega_b*ell/kappa (the Julia analogue of the Python TEST 2).
#   2. Two plates, weak adsorption: profile vs the cosh linearisation, and the
#      mirror symmetry the two identical plates must produce.
#   3. Interior ODE residual and both Robin surface BCs, with bc_reg on.
#   4. Wide gap: the plate pair decouples into two single plates, so the
#      Derjaguin pair potential vanishes; and U(d) is attractive as d shrinks.
#   5. Geometry dispatch: parallel_plates overrides a finite R (planar, per-area
#      energy), and the sphere solver is untouched by the flag.

using Test
using Droplet

@testset "bvp_parallel_plates" begin

    # -- 1. Single plate, linearised surface enrichment ------------------------
    @testset "SinglePlateLinear" begin
        # chi_b = 0 => alpha = 1/(phi(1-phi)), ell = sqrt(kappa/alpha);
        # the nonlinear correction is O(omega_b^2).
        p = Parameters(chi_b = 0.0, omega_b = -0.05, phi_inf = 0.3, kappa = 1.0,
                       L = 30.0, R = 5.0, eps = 1e-10, parallel_plates = false)
        s = solve_parallel_plates(p; plates = false, N = 3000)
        ell = sqrt(p.kappa * p.phi_inf * (1 - p.phi_inf))

        @test s.success
        @test isapprox(s.phi[1], p.phi_inf - p.omega_b * ell / p.kappa;
                       atol = 10 * p.omega_b^2)
        @test s.phi[end] ≈ p.phi_inf            # Dirichlet far field
        @test all(0 .< s.phi .< 1)
    end

    # -- 2. Two plates: cosh linearisation and mirror symmetry ----------------
    @testset "PlatePairLinear" begin
        p = Parameters(chi_b = 0.0, omega_b = -0.05, phi_inf = 0.3, kappa = 1.0,
                       L = 4.0, R = 5.0, eps = 1e-10, parallel_plates = true)
        s = solve_parallel_plates(p; N = 4000)
        psi, _, _ = psi_linear_plates(s.z, p)

        @test s.success
        @test maximum(abs, s.phi .- (p.phi_inf .+ psi)) < 10 * p.omega_b^2
        @test maximum(abs, s.phi .- reverse(s.phi)) < 1e-12   # identical plates
        @test s.phi[1] > p.phi_inf                            # attractive => enriched
        # gap centre is the least enriched point, and the tail never dips below phi_inf
        @test argmin(s.phi) == cld(length(s.phi), 2)
        @test minimum(s.phi) >= p.phi_inf
    end

    # -- 3. ODE residual and both Robin BCs (regularised) ---------------------
    @testset "Residuals" begin
        bc_reg = 5e-2
        p = Parameters(chi_b = 2.5, omega_b = -0.5, phi_inf = 0.05, kappa = 1.0,
                       L = 6.0, R = 5.0, eps = 1e-8, parallel_plates = true)
        s = solve_parallel_plates(p; N = 2000, bc_reg = bc_reg)
        @test s.success

        dz = s.z[2] - s.z[1]
        mu = compute_mu_inf(p)
        ode_res = maximum(abs(p.kappa * (s.phi[i+1] - 2s.phi[i] + s.phi[i-1]) / dz^2 -
                              (dfb_dphi(s.phi[i], p) - mu / p.nb))
                          for i in 2:length(s.z) - 1)
        @test ode_res < 1e-6

        # -kappa*phi'(0) + omega_b + bc_reg*log(phi_s/(1-phi_s)) = 0 at z = 0,
        # and the mirrored condition at z = L; one-sided second-order stencils.
        surf(phi_s) = (p.omega_b + bc_reg * log(phi_s / (1 - phi_s))) / p.kappa
        d0 = (-3s.phi[1] + 4s.phi[2] - s.phi[3]) / (2dz)
        dL = (3s.phi[end] - 4s.phi[end-1] + s.phi[end-2]) / (2dz)
        @test isapprox(d0,  surf(s.phi[1]);   atol = 1e-5)
        @test isapprox(dL, -surf(s.phi[end]); atol = 1e-5)
    end

    # -- 4. Derjaguin pair potential ------------------------------------------
    @testset "PairPotential" begin
        p = Parameters(chi_b = 2.5, omega_b = -0.5, phi_inf = 0.05, kappa = 1.0,
                       L = 40.0, R = 5.0, eps = 1e-8, parallel_plates = true)

        # A wide gap is two independent plates: W_plates = 2*W_single exactly.
        wide   = solve_parallel_plates(p; N = 4000, bc_reg = 5e-2)
        single = solve_parallel_plates(p; plates = false, N = 4000, bc_reg = 5e-2)
        @test wide.success && single.success
        @test isapprox(wide.phi[1], single.phi[1]; atol = 1e-9)   # films don't overlap
        @test isapprox(excess_free_energy(wide), 2 * excess_free_energy(single);
                       atol = 1e-12)

        # Adsorbed films overlap as the plates approach: U < 0 and deepening.
        ds = [1.0, 2.0, 4.0, 8.0]
        U, ok, W, W_single = plates_pair_potential(p, ds; N = 4000, bc_reg = 5e-2,
                                                   L_ref = 40.0)
        @test all(ok)
        @test isapprox(W_single, excess_free_energy(single); atol = 1e-12)
        @test all(U .< 0)
        @test issorted(U)                        # deepest at the smallest gap
        @test abs(U[end]) < abs(U[1]) / 100      # and negligible at the widest
    end

    # -- 5. Geometry dispatch --------------------------------------------------
    @testset "GeometryDispatch" begin
        # parallel_plates wins over a finite R: the slab energy is per unit area,
        # not the 4*pi*R^2-weighted spherical one, so it must not scale with R.
        base = (chi_b = 2.5, omega_b = -0.4, phi_inf = 0.05, kappa = 1.0,
                L = 6.0, eps = 1e-8, parallel_plates = true)
        s5  = solve_parallel_plates(Parameters(; R = 5.0,  base...); N = 1500, bc_reg = 5e-2)
        s50 = solve_parallel_plates(Parameters(; R = 50.0, base...); N = 1500, bc_reg = 5e-2)
        @test s5.phi ≈ s50.phi                                   # R unused by the plate BVP
        @test isapprox(excess_free_energy(s5), excess_free_energy(s50); atol = 1e-12)

        # ... and the sphere solver ignores the flag entirely.
        ps = Parameters(chi_b = 2.5, omega_b = -0.5, phi_inf = 0.09, kappa = 1.0,
                        L = 20.0, R = 5.0, eps = 1e-6)
        a = solve_isolated_sphere(ps; L_iso = ps.L, N = 400, bc_reg = 5e-2)
        b = solve_isolated_sphere(with(ps; parallel_plates = true);
                                  L_iso = ps.L, N = 400, bc_reg = 5e-2)
        @test a.phi == b.phi

        # A dense seed pins the capillary-condensed branch near contact.  Which
        # branch the default (linearised, dilute) guess drifts onto there depends
        # on the grid, so seeding is how a caller asks for a definite one.
        pc = Parameters(chi_b = 2.5, omega_b = -0.5, phi_inf = 0.05, kappa = 1.0,
                        L = 1.0, R = 5.0, eps = 1e-8, parallel_plates = true)
        dense = solve_parallel_plates(pc; N = 4000, bc_reg = 5e-2, phi_seed = 0.85)
        @test dense.success
        @test dense.phi[cld(length(dense.phi), 2)] > 0.5   # dense gap, not a thin film
    end
end
