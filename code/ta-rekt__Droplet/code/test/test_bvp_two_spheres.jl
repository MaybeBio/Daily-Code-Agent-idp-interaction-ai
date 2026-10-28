# test_bvp_two_spheres.jl — Comprehensive tests for the 2D axisymmetric
# two-sphere BVP solver in src/bvp_two_spheres.jl.  Ported from
# test/test_bvp_two_spheres.py (14 test classes).
#
# Within each nested testset, expensive BVP solves are computed once via a
# `let` block and reused across assertions (the Julia analogue of Python's
# class-scoped fixtures).

using Random
using LinearAlgebra
using Droplet
using Droplet: make_residual_inplace, _run_newton

# ---------------------------------------------------------------------------
# Shared parameter presets
# ---------------------------------------------------------------------------

P_benign() = Parameters(chi_b = 1.0, omega_b = -0.2, kappa = 1.0,
                        phi_inf = 0.1, R = 5.0, eps = 1e-6)

P_linear(omega_b = -0.05) =
    Parameters(chi_b = 0.0, omega_b = omega_b, kappa = 1.0,
               phi_inf = 0.2, R = 5.0, eps = 1e-6)

P_strong() = Parameters(chi_b = 2.0, omega_b = -0.5, kappa = 1.0,
                        phi_inf = 0.1, R = 5.0, eps = 1e-6)

P_strong2() = Parameters(chi_b = 4.0, omega_b = -1.35, kappa = 1.0,
                        phi_inf = 0.008, R = 10.0, eps = 1e-6)  
                        
P_strong3() = Parameters(chi_b = 4.0, omega_b = -1.35, kappa = 1.0,
                        phi_inf = 0.01, R = 10.0, eps = 1e-6)                          

function _Nz_for(p, d; dz_target = 0.15, z_buf = 4.0)
    span = 2.0 * (p.R + z_buf) + d
    return Int(floor(span / dz_target)) + 1
end

function _trapezoid(y::AbstractVector, x::AbstractVector)
    s = 0.0
    @inbounds for i in 1:length(x) - 1
        s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return s
end

function _avg_surface_enrichment(result; sphere::Int = 1)
    geom = result[:geom]
    phi  = result[:phi]
    p    = geom.p
    R    = p.R
    zc   = sphere == 1 ? geom.z_c1 : geom.z_c2
    theta = collect(range(0.0, π; length = 400))
    sin_t = sin.(theta)
    cos_t = cos.(theta)
    phi_s = interp_phi(phi, geom, R .* sin_t, zc .+ R .* cos_t)
    return _trapezoid((phi_s .- p.phi_inf) .* sin_t, theta) / 2.0
end

# Map a 1D isolated-sphere profile onto the 2D unknown vector, using each cell's
# distance to the nearer sphere surface.  For well-separated spheres this is an
# almost-exact solution, so it is a robust warm start for the stiff prewetting
# cases where the cold-start Newton stalls on the sharp adsorption front.
function _x0_from_1d(geom, sol1d)
    p = geom.p
    x0 = Vector{Float64}(undef, count(geom.unk_mask))
    k = 0
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        geom.unk_mask[i, j] || continue
        k += 1
        rho = geom.rho_grid[i]; z = geom.z_grid[j]
        s1 = sqrt(rho^2 + (z - geom.z_c1)^2) - p.R
        s2 = sqrt(rho^2 + (z - geom.z_c2)^2) - p.R
        x0[k] = clamp(sol1d(max(0.0, min(s1, s2))), p.eps, 1.0 - p.eps)
    end
    return x0
end

@testset "bvp_two_spheres" begin

# ===========================================================================
# 1. TestGeometry
# ===========================================================================

@testset "Geometry" begin
    geom = build_geometry(P_benign(), 2.0; Nrho = 40, Nz = 120)

    # test_grid_shapes
    @test size(geom.rho_grid) == (geom.Nrho,)
    @test size(geom.z_grid)   == (geom.Nz,)
    for m in (geom.inside_mask, geom.fluid_mask, geom.dir_mask, geom.unk_mask)
        @test size(m) == (geom.Nrho, geom.Nz)
        @test eltype(m) == Bool
    end

    # test_grid_endpoints
    p   = geom.p
    ell = decay_length(p)
    z_buf = 4.0 * ell
    @test geom.rho_grid[1]   ≈ 0.0 atol = 1e-12
    @test geom.rho_grid[end] ≈ p.R + 4.0 * ell
    @test geom.z_grid[1]     ≈ -p.R - z_buf
    @test geom.z_grid[end]   ≈ geom.z_c2 + p.R + z_buf

    # test_spacing_positive
    @test geom.drho ≈ geom.rho_grid[2] - geom.rho_grid[1]
    @test geom.dz   ≈ geom.z_grid[2]   - geom.z_grid[1]
    @test geom.drho > 0 && geom.dz > 0

    # test_sphere_centres
    @test geom.z_c1 ≈ 0.0
    @test geom.z_c2 ≈ 2.0 * geom.p.R + geom.d

    # test_mask_partition: inside XOR fluid covers every cell exactly once
    @test all(geom.inside_mask .⊻ geom.fluid_mask)
    @test !any(geom.inside_mask .& geom.fluid_mask)

    # test_dirichlet_unknown_split
    @test !any(geom.dir_mask .& geom.unk_mask)
    @test all((geom.dir_mask .| geom.unk_mask) .== geom.fluid_mask)

    # test_dirichlet_faces: outer cylinder + z end-caps Dirichlet wherever fluid
    @test all(geom.dir_mask[end, :] .== geom.fluid_mask[end, :])
    @test all(geom.dir_mask[:, 1]   .== geom.fluid_mask[:, 1])
    @test all(geom.dir_mask[:, end] .== geom.fluid_mask[:, end])

    # test_sphere_classification: cells near each centre are inside, far cells fluid
    i0   = 1     # rho = 0
    j_c1 = argmin(abs.(geom.z_grid .- geom.z_c1))
    j_c2 = argmin(abs.(geom.z_grid .- geom.z_c2))
    @test geom.inside_mask[i0, j_c1]
    @test geom.inside_mask[i0, j_c2]
    i_far = geom.Nrho
    j_mid = argmin(abs.(geom.z_grid .- 0.5 * (geom.z_c1 + geom.z_c2)))
    @test geom.fluid_mask[i_far, j_mid]

    # test_no_sphere_overlap
    n_inside = count(geom.inside_mask)
    cell_area = geom.drho * geom.dz
    upper = 2.5 * π * geom.p.R^2 / cell_area
    @test 0 < n_inside < upper

    # test_ghost_cells_are_inside
    if length(geom.ghost_i) > 0
        for k in 1:length(geom.ghost_i)
            @test geom.inside_mask[geom.ghost_i[k], geom.ghost_j[k]]
        end
    end

    # test_ghost_cells_border_fluid
    for k in 1:length(geom.ghost_i)
        i = geom.ghost_i[k]; j = geom.ghost_j[k]
        has_fluid = false
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            if 1 <= ii <= geom.Nrho && 1 <= jj <= geom.Nz && geom.fluid_mask[ii, jj]
                has_fluid = true
                break
            end
        end
        @test has_fluid
    end

    # test_ghost_weights_sum_to_one: P1 sums to 1; P2 sums to 1 (quadratic) or 0 (linear)
    if length(geom.ghost_i) > 0
        sum1 = vec(sum(geom.ghost_wts_1; dims = 1))
        @test all(isapprox.(sum1, 1.0; atol = 1e-12))
        sum2 = vec(sum(geom.ghost_wts_2; dims = 1))
        quad = .!geom.ghost_used_linear
        @test all(isapprox.(sum2[quad], 1.0; atol = 1e-12))
        @test all(isapprox.(sum2[.!quad], 0.0; atol = 1e-12))
    end

    # test_ghost_image_stencil_avoids_inside
    for k in 1:length(geom.ghost_i)
        i0 = geom.ghost_i0_1[k]; j0 = geom.ghost_j0_1[k]
        for (di, dj) in ((0, 0), (1, 0), (0, 1), (1, 1))
            @test !geom.inside_mask[i0 + di, j0 + dj]
        end
    end
    for k in 1:length(geom.ghost_i)
        geom.ghost_used_linear[k] && continue
        i0 = geom.ghost_i0_2[k]; j0 = geom.ghost_j0_2[k]
        for (di, dj) in ((0, 0), (1, 0), (0, 1), (1, 1))
            @test !geom.inside_mask[i0 + di, j0 + dj]
        end
    end

    # test_ghost_d1_d2_positive
    if length(geom.ghost_d1) > 0
        @test all(geom.ghost_d1 .> 0)
        quad = .!geom.ghost_used_linear
        sep = geom.ghost_d2[quad] .- geom.ghost_d1[quad]
        min_sep = 0.9 * max(geom.drho, geom.dz)
        @test all(sep .>= min_sep - 1e-12)
    end

    # test_ghost_constant_partition: alpha1 + alpha2 = 1 for quadratic ghosts
    if length(geom.ghost_coef_phi1) > 0
        s = geom.ghost_coef_phi1 .+ geom.ghost_coef_phi2
        @test all(isapprox.(s, 1.0; atol = 1e-12))
    end
end


# ===========================================================================
# 2. TestScatter
# ===========================================================================

@testset "Scatter" begin
    geom = build_geometry(P_benign(), 2.0; Nrho = 40, Nz = 120)

    # test_output_shape
    let x = fill(geom.p.phi_inf, count(geom.unk_mask))
        phi = scatter(x, geom)
        @test size(phi) == (geom.Nrho, geom.Nz)
    end

    # test_unknowns_roundtrip
    let rng = MersenneTwister(0)
        n = count(geom.unk_mask)
        x = 0.5 .+ 0.1 .* randn(rng, n)
        phi = scatter(x, geom)
        @test isapprox(phi[geom.unk_mask], x; atol = 1e-12)
    end

    # test_dirichlet_cells_phi_inf
    let x = fill(0.7, count(geom.unk_mask))
        phi = scatter(x, geom)
        @test all(isapprox.(phi[geom.dir_mask], geom.p.phi_inf; atol = 1e-12))
    end

    # test_uniform_omega_zero
    let p = with(P_benign(); omega_b = 0.0)
        g2 = build_geometry(p, 2.0; Nrho = 40, Nz = 120)
        x = fill(p.phi_inf, count(g2.unk_mask))
        phi = scatter(x, g2)
        @test all(isapprox.(phi, p.phi_inf; atol = 1e-12))
    end

    # test_ghost_fill_formula
    let p = P_benign()
        g3 = build_geometry(p, 2.0; Nrho = 40, Nz = 120)
        rng = MersenneTwister(7)
        n = count(g3.unk_mask)
        x = p.phi_inf .+ 0.05 .* randn(rng, n)
        phi = scatter(x, g3)
        if length(g3.ghost_i) > 0
            for k in 1:length(g3.ghost_i)
                i01 = g3.ghost_i0_1[k]; j01 = g3.ghost_j0_1[k]
                i02 = g3.ghost_i0_2[k]; j02 = g3.ghost_j0_2[k]
                phi_P1 = g3.ghost_wts_1[1, k] * phi[i01,     j01]     +
                         g3.ghost_wts_1[2, k] * phi[i01 + 1, j01]     +
                         g3.ghost_wts_1[3, k] * phi[i01,     j01 + 1] +
                         g3.ghost_wts_1[4, k] * phi[i01 + 1, j01 + 1]
                phi_P2 = g3.ghost_wts_2[1, k] * phi[i02,     j02]     +
                         g3.ghost_wts_2[2, k] * phi[i02 + 1, j02]     +
                         g3.ghost_wts_2[3, k] * phi[i02,     j02 + 1] +
                         g3.ghost_wts_2[4, k] * phi[i02 + 1, j02 + 1]
                expected = g3.ghost_coef_phi1[k] * phi_P1 +
                           g3.ghost_coef_phi2[k] * phi_P2 +
                           g3.ghost_coef_bc[k] * (p.omega_b / p.kappa)
                actual = phi[g3.ghost_i[k], g3.ghost_j[k]]
                @test isapprox(actual, expected; atol = 1e-12)
            end
        end
    end

    # test_ghost_robin_bc_quadratic
    let p = P_benign()
        g3 = build_geometry(p, 2.0; Nrho = 50, Nz = 150)
        rng = MersenneTwister(3)
        x = p.phi_inf .+ 0.1 .* randn(rng, count(g3.unk_mask))
        phi = scatter(x, g3)
        if length(g3.ghost_i) > 0
            for k in 1:length(g3.ghost_i)
                g3.ghost_used_linear[k] && continue
                i01 = g3.ghost_i0_1[k]; j01 = g3.ghost_j0_1[k]
                i02 = g3.ghost_i0_2[k]; j02 = g3.ghost_j0_2[k]
                phi_P1 = g3.ghost_wts_1[1, k] * phi[i01,     j01]     +
                         g3.ghost_wts_1[2, k] * phi[i01 + 1, j01]     +
                         g3.ghost_wts_1[3, k] * phi[i01,     j01 + 1] +
                         g3.ghost_wts_1[4, k] * phi[i01 + 1, j01 + 1]
                phi_P2 = g3.ghost_wts_2[1, k] * phi[i02,     j02]     +
                         g3.ghost_wts_2[2, k] * phi[i02 + 1, j02]     +
                         g3.ghost_wts_2[3, k] * phi[i02,     j02 + 1] +
                         g3.ghost_wts_2[4, k] * phi[i02 + 1, j02 + 1]
                d1, d2, delta = g3.ghost_d1[k], g3.ghost_d2[k], g3.ghost_delta[k]
                # Inverted closed form: solve c, b from the (P1, P2, Robin BC)
                slope_avg = (phi_P1 - phi_P2) / (d1 - d2)
                c = (slope_avg - p.omega_b / p.kappa) / (d1 + d2 - 2.0 * delta)
                b = p.omega_b / p.kappa - 2.0 * c * delta
                bc_check = b + 2.0 * c * delta
                @test isapprox(bc_check, p.omega_b / p.kappa; atol = 1e-12)
            end
        end
    end
end


# ===========================================================================
# 3. TestLaplacian
# ===========================================================================

@testset "Laplacian" begin
    geom = build_geometry(P_benign(), 2.0; Nrho = 50, Nz = 150)
    interior_i = 2:geom.Nrho - 1
    interior_j = 2:geom.Nz - 1

    # test_constant
    let phi = fill(0.3, geom.Nrho, geom.Nz)
        lap = laplacian(phi, geom)
        @test all(isapprox.(lap[interior_i, interior_j], 0.0; atol = 1e-10))
    end

    # test_linear_in_z
    let Z = repeat(geom.z_grid', geom.Nrho, 1)
        lap = laplacian(Float64.(Z), geom)
        @test all(isapprox.(lap[interior_i, interior_j], 0.0; atol = 1e-10))
    end

    # test_quadratic_in_z
    let Z = repeat(geom.z_grid', geom.Nrho, 1)
        lap = laplacian(Float64.(Z) .^ 2, geom)
        @test all(isapprox.(lap[interior_i, interior_j], 2.0; atol = 1e-10))
    end

    # test_quadratic_in_rho: phi = rho^2 ⇒ Lap = 4 in cylindrical coords
    let Rho = repeat(geom.rho_grid, 1, geom.Nz)
        lap = laplacian(Float64.(Rho) .^ 2, geom)
        @test all(isapprox.(lap[interior_i, interior_j], 4.0; atol = 1e-10))
    end

    # test_axis_quadratic_rho: axis row of phi = rho^2 gives Lap = 4
    let Rho = repeat(geom.rho_grid, 1, geom.Nz)
        lap = laplacian(Float64.(Rho) .^ 2, geom)
        @test all(isapprox.(lap[1, interior_j], 4.0; atol = 1e-10))
    end

    # test_axis_linear_in_z: axis row of phi = z gives Lap = 0
    let Z = repeat(geom.z_grid', geom.Nrho, 1)
        lap = laplacian(Float64.(Z), geom)
        @test all(isapprox.(lap[1, interior_j], 0.0; atol = 1e-10))
    end
end


# ===========================================================================
# 4. TestInterp
# ===========================================================================

@testset "Interp" begin
    geom = build_geometry(P_benign(), 2.0; Nrho = 40, Nz = 120)

    # test_grid_point_recovery
    let rng = MersenneTwister(0)
        phi = randn(rng, geom.Nrho, geom.Nz)
        ii = 1:5:geom.Nrho
        jj = 1:5:geom.Nz
        rho_q = Float64[geom.rho_grid[i] for i in ii for _ in jj]
        z_q   = Float64[geom.z_grid[j]   for _ in ii for j in jj]
        expected = Float64[phi[i, j] for i in ii for j in jj]
        out = interp_phi(phi, geom, rho_q, z_q)
        keep = (rho_q .< geom.rho_grid[end] - 1e-9) .&
               (z_q   .< geom.z_grid[end]   - 1e-9)
        @test all(isapprox.(out[keep], expected[keep]; atol = 1e-10))
    end

    # test_exact_on_linear_field
    let a = 0.7, b = -1.3
        Rho = repeat(geom.rho_grid, 1, geom.Nz)
        Z   = repeat(geom.z_grid', geom.Nrho, 1)
        phi = a .* Rho .+ b .* Z
        rng = MersenneTwister(1)
        rho_q = 0.1 .+ (geom.rho_grid[end] - 0.2) .* rand(rng, 20)
        z_q   = geom.z_grid[1] + 0.1 .+ (geom.z_grid[end] - geom.z_grid[1] - 0.2) .* rand(rng, 20)
        out = interp_phi(phi, geom, rho_q, z_q)
        expected = a .* rho_q .+ b .* z_q
        @test all(isapprox.(out, expected; atol = 1e-10))
    end

    # test_out_of_range_clamped
    let phi = fill(0.42, geom.Nrho, geom.Nz)
        out = interp_phi(phi, geom,
                         [-1.0, geom.rho_grid[end] + 5.0, 1.0],
                         [geom.z_grid[1] - 5.0, 0.0, geom.z_grid[end] + 5.0])
        @test all(isapprox.(out, 0.42))
    end
end


# ===========================================================================
# 5. TestInitialGuess
# ===========================================================================

@testset "InitialGuess" begin
    # test_length_matches_unknowns
    let geom = build_geometry(P_benign(), 2.0; Nrho = 40, Nz = 120)
        x0 = initial_guess(geom)
        @test length(x0) == count(geom.unk_mask)
    end

    # test_values_clipped
    let geom = build_geometry(P_benign(), 2.0; Nrho = 40, Nz = 120)
        x0 = initial_guess(geom)
        @test all(x0 .> geom.p.eps)
        @test all(x0 .< 1.0 - geom.p.eps)
    end

    # test_omega_zero_uniform
    let p = with(P_benign(); omega_b = 0.0)
        geom = build_geometry(p, 2.0; Nrho = 40, Nz = 120)
        x0 = initial_guess(geom)
        @test all(isapprox.(x0, p.phi_inf; atol = 1e-12))
    end
end


# ===========================================================================
# 6. TestSolverBasic
# ===========================================================================

@testset "SolverBasic" begin
    result = solve_two_spheres(P_benign(), 2.0;
                               Nrho = 60, Nz = 180, verbose = false)

    # test_keys
    for k in (:phi, :geom, :p, :d, :mu_inf, :success, :x)
        @test haskey(result, k)
    end

    # test_success
    @test result[:success]

    # test_phi_in_unit_interval
    let phi_fluid = result[:phi][result[:geom].fluid_mask]
        @test all(isfinite, phi_fluid)
        @test all(phi_fluid .> 0)
        @test all(phi_fluid .< 1)
    end

    # test_residual_below_tol
    let geom = result[:geom]
        res_fn = make_residual(geom, result[:mu_inf])
        F = res_fn(result[:x])
        @test maximum(abs, F) < 1e-5
    end

    # test_midplane_reflection_symmetry: phi(rho, z) == phi(rho, z_c1+z_c2-z)
    let geom = result[:geom], phi = result[:phi]
        z_mirror = geom.z_c1 + geom.z_c2 .- geom.z_grid
        rho_full = repeat(geom.rho_grid, 1, geom.Nz)
        z_full   = repeat(reshape(z_mirror, 1, :), geom.Nrho, 1)
        phi_mirror = reshape(
            interp_phi(phi, geom, vec(rho_full), vec(z_full)),
            geom.Nrho, geom.Nz)
        diff = abs.(phi .- phi_mirror)
        diff[.!geom.fluid_mask] .= 0.0
        @test maximum(diff) < 1e-5
    end

    # test_far_field_phi_inf
    let geom = result[:geom], phi = result[:phi], p = result[:p]
        i_far = geom.Nrho - 2
        j_mid = argmin(abs.(geom.z_grid .- 0.5 * (geom.z_c1 + geom.z_c2)))
        @test abs(phi[i_far, j_mid] - p.phi_inf) < 0.01
    end
end


# ===========================================================================
# 7. TestOmegaZero
# ===========================================================================

@testset "OmegaZero" begin
    p = Parameters(chi_b = 0.0, omega_b = 0.0, kappa = 1.0,
                   phi_inf = 0.2, R = 5.0, eps = 1e-6)
    result = solve_two_spheres(p, 2.0; Nrho = 50, Nz = 140, verbose = false)

    @test result[:success]
    @test all(isapprox.(result[:phi][result[:geom].fluid_mask],
                        result[:geom].p.phi_inf; atol = 1e-8))
    @test abs(excess_free_energy_two_spheres(result)) < 1e-8
end


# ===========================================================================
# 8. TestSingleSphereLimit
# ===========================================================================

@testset "SingleSphereLimit" begin
    # At large separation each sphere is effectively isolated, so the 2D radial
    # profile must match the 1D spherical BVP.  The cases sweep from a benign
    # regime up to P_strong2 / P_strong3, which sit just before and just after a
    # surface (prewetting) phase transition (only phi_inf differs).  On the
    # post-transition branch the bare adsorption BC is bistable — a thin film and
    # a runaway phi_s -> 1 — so those cases share a regularized surface BC
    # (bc_reg > 0, see _surf_slope / scatter!) between the 1D and 2D solvers,
    # which pins the same finite branch in both.  The 1D solution also warm-starts
    # the 2D solve where the sharp front stalls a cold Newton start.
    cases = [
        (name = "benign",  p = P_benign(),  bc_reg = 0.0,  d = 8.0 * P_benign().R,
         Nrho = 100, dz_target = 0.10, rho_buf = 8.0, z_buf = 4.0,
         N1d = 800,  L_iso = 10.0, warm = false, tol_prof = 0.005, tol_iso = 0.01),
        (name = "strong",  p = P_strong(),  bc_reg = 0.02, d = 6.0 * P_strong().R,
         Nrho = 140, dz_target = 0.06, rho_buf = 8.0, z_buf = 4.0,
         N1d = 800,  L_iso = 10.0, warm = false, tol_prof = 0.005, tol_iso = 0.01),
        (name = "strong2", p = P_strong2(), bc_reg = 0.02, d = 10.0,
         Nrho = 351, dz_target = 0.04, rho_buf = 4.0, z_buf = 3.0,
         N1d = 1400, L_iso = 12.0, warm = true,  tol_prof = 0.005, tol_iso = 0.01),
        (name = "strong3", p = P_strong3(), bc_reg = 0.05, d = 10.0,
         Nrho = 467, dz_target = 0.03, rho_buf = 4.0, z_buf = 3.0,
         N1d = 1400, L_iso = 12.0, warm = true,  tol_prof = 0.006, tol_iso = 0.01),
    ]

    for c in cases
        @testset "$(c.name)" begin
            p     = c.p
            L_iso = c.L_iso
            Nz    = _Nz_for(p, c.d; dz_target = c.dz_target, z_buf = c.z_buf)

            sol1d = solve_isolated_sphere(with(p; L = L_iso);
                                          L_iso = L_iso, N = c.N1d, bc_reg = c.bc_reg)
            @test sol1d.success

            # Warm start the stiff post-transition solves from the 1D profile.
            x0 = nothing
            if c.warm
                g0 = build_geometry(p, c.d; rho_max = p.R + c.rho_buf, z_buf = c.z_buf,
                                    Nrho = c.Nrho, Nz = Nz, bc_reg = c.bc_reg)
                x0 = _x0_from_1d(g0, sol1d)
            end
            res = solve_two_spheres(p, c.d;
                                    Nrho = c.Nrho, Nz = Nz,
                                    rho_max = p.R + c.rho_buf, z_buf = c.z_buf,
                                    bc_reg = c.bc_reg, x0_init = x0, verbose = false)
            @test res[:success]
            geom = res[:geom]

            # test_axial_profile_matches_1d
            let z_query = collect(range(p.R + 0.05, p.R + 0.5 * L_iso; length = 30))
                phi_2d = interp_phi(res[:phi], geom, zeros(length(z_query)), z_query)
                phi_1d = sol1d(z_query .- p.R)
                @test maximum(abs.(phi_2d .- phi_1d)) < c.tol_prof
            end

            # test_equatorial_profile_matches_1d
            let rho_query = collect(range(p.R + 0.05, p.R + 0.5 * L_iso; length = 30))
                phi_2d = interp_phi(res[:phi], geom, rho_query, zeros(length(rho_query)))
                phi_1d = sol1d(rho_query .- p.R)
                @test maximum(abs.(phi_2d .- phi_1d)) < c.tol_prof
            end

            # test_spherical_isotropy
            for s in (0.3, 0.7, 1.5)
                r = p.R + s
                theta_samples = [1e-3, π / 4.0, π / 2.0]
                rho_q = r .* sin.(theta_samples)
                z_q   = geom.z_c1 .+ r .* cos.(theta_samples)
                phi_at_angles = interp_phi(res[:phi], geom, rho_q, z_q)
                spread = maximum(phi_at_angles) - minimum(phi_at_angles)
                @test spread < c.tol_iso
            end
        end
    end
end


# ===========================================================================
# 9. TestLinearResponseLimit
# ===========================================================================

@testset "LinearResponseLimit" begin
    OMEGAS = (-0.10, -0.05, -0.025)
    errs = Tuple{Float64, Float64}[]

    for omega_b in OMEGAS
        p = with(P_linear(); omega_b = omega_b)
        res = solve_two_spheres(p, 2.0 * p.R;
                                Nrho = 80, Nz = 240, verbose = false)
        @test res[:success]
        geom = res[:geom]
        # Build the linear-superposition reference field
        phi_pred = Matrix{Float64}(undef, geom.Nrho, geom.Nz)
        @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
            rho_i = geom.rho_grid[i]
            z_j   = geom.z_grid[j]
            s1 = max(0.0, sqrt(rho_i^2 + (z_j - geom.z_c1)^2) - p.R)
            s2 = max(0.0, sqrt(rho_i^2 + (z_j - geom.z_c2)^2) - p.R)
            psi1, _, _ = psi_linear(s1, p)
            psi2, _, _ = psi_linear(s2, p)
            phi_pred[i, j] = p.phi_inf + psi1 + psi2
        end
        diff = abs.(res[:phi] .- phi_pred)
        diff[.!geom.fluid_mask] .= 0.0
        push!(errs, (omega_b, maximum(diff)))
    end

    # test_per_omega_bound
    for (omega_b, err) in errs
        budget = 0.5 * omega_b^2 + 0.005
        @test err < budget
    end

    # test_error_scales_with_omega
    ratio = errs[1][2] / errs[end][2]
    @test ratio > 2.0
end


# ===========================================================================
# 10. TestUDecay
# ===========================================================================

@testset "UDecay" begin
    D_REF_FACTOR = 6.0
    D_FACTORS = (0.1, 0.3, 0.8, 2.0, 4.0)
    p = Parameters(chi_b = 1.5, omega_b = -0.6, kappa = 1.0,
                   phi_inf = 0.1, R = 5.0, eps = 1e-6)
    Nrho = 120

    function solve_one(d)
        return solve_two_spheres(p, d;
                                 Nrho = Nrho,
                                 Nz = _Nz_for(p, d; dz_target = 0.08, z_buf = 3.0),
                                 rho_max = p.R + 4.0, z_buf = 3.0,
                                 verbose = false)
    end

    d_ref   = D_REF_FACTOR * p.R
    res_ref = solve_one(d_ref)
    F_ref   = excess_free_energy_two_spheres(res_ref)
    Us = Float64[]
    results = []
    for factor in D_FACTORS
        d  = factor * p.R
        rk = solve_one(d)
        Fk = excess_free_energy_two_spheres(rk)
        push!(Us, Fk - F_ref)
        push!(results, rk)
    end

    # test_all_converged
    @test res_ref[:success]
    @test all(isfinite, res_ref[:phi][res_ref[:geom].fluid_mask])
    for r in results
        @test r[:success]
        @test all(isfinite, r[:phi][r[:geom].fluid_mask])
    end

    # test_small_d_attractive
    @test Us[1] < -0.5

    # test_signal_vs_noise
    @test abs(Us[1]) > 4.0 * abs(Us[end])

    # test_large_d_small
    @test abs(Us[end]) < 0.3

    # test_envelope_decays
    mags = abs.(Us)
    close_mean = mean(mags[1:2])
    far_mean   = mean(mags[end - 1:end])
    @test close_mean > 3.0 * far_mean
end


# ===========================================================================
# 11. TestExcessFreeEnergy
# ===========================================================================

@testset "ExcessFreeEnergy" begin
    # test_zero_for_uniform
    let p = Parameters(chi_b = 0.0, omega_b = 0.0, kappa = 1.0,
                       phi_inf = 0.2, R = 5.0, eps = 1e-6)
        res = solve_two_spheres(p, 2.0; Nrho = 50, Nz = 140, verbose = false)
        @test abs(excess_free_energy_two_spheres(res)) < 1e-8
    end

    # test_attractive_F_negative
    let p = Parameters(chi_b = 1.0, omega_b = -0.3, kappa = 1.0,
                       phi_inf = 0.1, R = 5.0, eps = 1e-6)
        res = solve_two_spheres(p, 2.0; Nrho = 60, Nz = 180, verbose = false)
        @test excess_free_energy_two_spheres(res) < 0.0
    end

    # test_sign_flip_with_omega
    let p_neg = Parameters(chi_b = 1.0, omega_b = -0.3, kappa = 1.0,
                           phi_inf = 0.1, R = 5.0, eps = 1e-6)
        p_pos = Parameters(chi_b = 1.0, omega_b = +0.3, kappa = 1.0,
                           phi_inf = 0.1, R = 5.0, eps = 1e-6)
        F_neg = excess_free_energy_two_spheres(
            solve_two_spheres(p_neg, 2.0; Nrho = 60, Nz = 180, verbose = false))
        F_pos = excess_free_energy_two_spheres(
            solve_two_spheres(p_pos, 2.0; Nrho = 60, Nz = 180, verbose = false))
        @test F_pos > F_neg
    end

    # test_R_scaling: omega_b * 4pi*R^2 * <phi_surf - phi_inf> scales > R^2
    let p4 = Parameters(R = 4.0, chi_b = 0.0, omega_b = -0.1, kappa = 1.0,
                        phi_inf = 0.2, eps = 1e-6)
        p8 = Parameters(R = 8.0, chi_b = 0.0, omega_b = -0.1, kappa = 1.0,
                        phi_inf = 0.2, eps = 1e-6)
        res4 = solve_two_spheres(p4, 2.0; Nrho = 80, Nz = 240, verbose = false)
        res8 = solve_two_spheres(p8, 2.0; Nrho = 80, Nz = 240, verbose = false)
        avg4 = _avg_surface_enrichment(res4)
        avg8 = _avg_surface_enrichment(res8)
        surf4 = p4.R^2 * avg4
        surf8 = p8.R^2 * avg8
        @test surf8 / surf4 > 4.0
    end

    # test_idempotent
    let p = P_benign()
        res = solve_two_spheres(p, 2.0; Nrho = 50, Nz = 140, verbose = false)
        F1 = excess_free_energy_two_spheres(res)
        F2 = excess_free_energy_two_spheres(res)
        @test F1 == F2
    end
end


# ===========================================================================
# 12. TestGridConvergence (slow)
# ===========================================================================

@testset "GridConvergence" begin
    p = with(P_linear(-0.05); phi_inf = 0.2, chi_b = 0.0, R = 5.0)
    truth = psi_linear([0.0], p)[1][1]   # psi at the sphere surface, s=0

    levels = []
    for (Nrho, Nz) in ((40, 120), (80, 240), (160, 480))
        res = solve_two_spheres(p, 6.0 * p.R;
                                Nrho = Nrho, Nz = Nz,
                                rho_max = p.R + 8.0, z_buf = 4.0,
                                verbose = false)
        @test res[:success]
        avg = _avg_surface_enrichment(res)
        push!(levels, (Nrho, avg, abs(avg - truth)))
    end

    errs = [lv[3] for lv in levels]

    # test_monotone_convergence
    @test errs[1] > errs[2] > errs[3]

    # test_second_order_rate
    ratio_total = errs[1] / errs[3]
    ratio_step1 = errs[1] / errs[2]
    ratio_step2 = errs[2] / errs[3]
    @test ratio_total > 6.0
    @test ratio_step1 > 2.0 && ratio_step2 > 2.0
end


# ===========================================================================
# 13. TestDecayLength
# ===========================================================================

@testset "DecayLength" begin
    # test_benign_matches_analytic
    let p = Parameters(chi_b = 0.0, omega_b = 0.0, kappa = 1.0,
                       phi_inf = 0.5, R = 5.0)
        d2 = d2fb_dphi2(p.phi_inf, p)
        expected = sqrt(p.kappa / d2)
        @test isapprox(decay_length(p), expected; rtol = 1e-10)
    end

    # test_spinodal_fallback
    let p = Parameters(chi_b = 2.0, omega_b = 0.0, kappa = 1.0,
                       phi_inf = 0.5, R = 5.0, nb = 1.0)
        ell = decay_length(p)
        expected = sqrt(p.kappa * p.nb / 0.05)
        @test isapprox(ell, expected; rtol = 1e-10)
        @test ell > 0 && isfinite(ell)
    end

    # test_always_positive across a few (phi_inf, chi_b) combos
    for (phi_inf, chi_b) in ((0.05, 0.0), (0.2, 1.0), (0.5, 1.5), (0.8, 0.0))
        p = Parameters(phi_inf = phi_inf, chi_b = chi_b, kappa = 1.0, R = 5.0)
        @test decay_length(p) > 0
    end
end


# ===========================================================================
# 14. TestRobustness
# ===========================================================================

@testset "Robustness" begin
    # test_perturbed_initial_guess
    let p = P_benign()
        geom = build_geometry(p, 2.0; Nrho = 60, Nz = 180)
        res! = make_residual_inplace(geom, compute_mu_inf(p))
        x0   = initial_guess(geom)
        rng  = MersenneTwister(11)
        x0_perturb = clamp.(x0 .+ 0.01 .* randn(rng, length(x0)),
                            p.eps, 1.0 - p.eps)
        sol_a = _run_newton(res!, x0,        1e-6, 200, false)
        sol_b = _run_newton(res!, x0_perturb, 1e-6, 200, false)
        @test maximum(abs.(sol_a.u .- sol_b.u)) < 1e-4
    end

    # test_parameters_not_mutated (Parameters is immutable, but verify anyway)
    let p = P_benign()
        snapshot = (p.chi_b, p.omega_b, p.kappa, p.phi_inf, p.R, p.nb, p.L,
                    p.eps, p.conserve_volume, p.parallel_plates)
        _ = solve_two_spheres(p, 2.0; Nrho = 40, Nz = 120, verbose = false)
        @test (p.chi_b, p.omega_b, p.kappa, p.phi_inf, p.R, p.nb, p.L,
               p.eps, p.conserve_volume, p.parallel_plates) == snapshot
    end

    # test_returned_fields_serialisable: result dict carries finite arrays
    let p = P_benign()
        res = solve_two_spheres(p, 2.0; Nrho = 40, Nz = 120, verbose = false)
        @test size(res[:phi]) == (res[:geom].Nrho, res[:geom].Nz)
        @test isfinite(res[:d])
        @test isfinite(res[:mu_inf])
        @test eltype(res[:phi]) <: AbstractFloat
        @test length(res[:x]) == count(res[:geom].unk_mask)
    end
end

end  # @testset "bvp_two_spheres"
