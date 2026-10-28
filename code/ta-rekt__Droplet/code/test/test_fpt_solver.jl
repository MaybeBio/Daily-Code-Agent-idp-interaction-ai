# test_fpt_solver.jl — Tests for the general 1D backward-FPE first-passage-time
# solver in src/fpt_solver.jl.
#
# Coverage:
#   1. Free diffusion vs closed-form MFPT and slowest decay rate.
#   2. Constant-force potential: eigen-MFPT vs exact quadrature and vs a direct
#      linear solve of the discrete backward equation (tight cross-check).
#   3. Double well: eigen-MFPT vs continuum quadrature.
#   4. Detailed-balance / positivity sanity (real λ>0, f(t)≥0, ∫f dt≈1, both
#      ends absorbing).
#   5. Empirical Langevin (Euler–Maruyama Brownian dynamics) in the SAME
#      potential — the apples-to-apples check of the solver against sampling.

using Test
using Random
using Statistics
using LinearAlgebra
using Droplet
using Droplet: _build_generator, _interior_indices

trapz(y, x) = sum((x[2:end] .- x[1:end-1]) .* (y[2:end] .+ y[1:end-1]) ./ 2)

@testset "fpt_solver" begin

    # -- 1. Free diffusion: analytic MFPT and slowest rate --------------------
    @testset "FreeDiffusion" begin
        L = 2.0; D = 1.0; N = 401
        x = collect(range(0.0, L; length = N))
        prob = FPTProblem(x, zeros(N); D = D, kBT = 1.0,
                          left_bc = :reflecting, right_bc = :absorbing)
        res = solve_fpt(prob)

        # uniform (Boltzmann-flat) start: <τ> = L²/(3D)
        @test isapprox(res.mfpt, L^2 / (3D); rtol = 0.02)
        # slowest decay rate λ₁ = D π² / (4 L²)
        @test isapprox(res.rates[1], D * pi^2 / (4 * L^2); rtol = 0.02)
        # survival starts at 1 (∑ weights = S(0)), decays monotonically, PDF normalises
        @test isapprox(sum(res.weights), 1.0; atol = 1e-8)
        @test isapprox(res.survival[1], 1.0; atol = 0.01)
        @test all(diff(res.survival) .<= 1e-10)
        @test isapprox(trapz(res.pdf, res.t), 1.0; atol = 0.05)
    end

    # -- 2. Constant force: quadrature and direct-solve cross-checks -----------
    @testset "ConstantForce" begin
        L = 2.0; D = 1.0; g = 1.5; N = 801
        x = collect(range(0.0, L; length = N))
        U = g .* x
        prob = FPTProblem(x, U; D = D, kBT = 1.0,
                          left_bc = :reflecting, right_bc = :absorbing)
        res = solve_fpt(prob)

        # continuum double-integral MFPT
        @test isapprox(res.mfpt, fpt_mfpt_quadrature(prob); rtol = 5e-3)

        # direct linear solve of the discrete backward equation Q̃ᵀ τ = -1
        Q̃, idx = _build_generator(prob)
        τ = Q̃' \ fill(-1.0, length(idx))
        Ui = U[idx]
        w = exp.(-(Ui .- minimum(Ui)))
        p0 = w ./ sum(w)
        @test isapprox(res.mfpt, dot(p0, τ); rtol = 1e-8)
    end

    # -- 3. Double well: eigen-MFPT vs continuum quadrature --------------------
    @testset "DoubleWell" begin
        H = 2.0; D = 1.0; N = 1001
        x = collect(range(-1.5, 1.5; length = N))
        U = H .* (x .^ 2 .- 1) .^ 2
        prob = FPTProblem(x, U; D = D, kBT = 1.0,
                          left_bc = :reflecting, right_bc = :absorbing)
        res = solve_fpt(prob)

        @test all(res.rates .> 0)
        @test isapprox(res.mfpt, fpt_mfpt_quadrature(prob); rtol = 5e-3)
    end

    # -- 4. Detailed balance / positivity -------------------------------------
    @testset "Positivity" begin
        H = 2.0; N = 601
        x = collect(range(-1.5, 1.5; length = N))
        U = H .* (x .^ 2 .- 1) .^ 2

        prob = FPTProblem(x, U; left_bc = :reflecting, right_bc = :absorbing)
        res = solve_fpt(prob)
        @test all(isreal, res.rates) && all(res.rates .> 0)
        @test all(res.pdf .>= -1e-10)
        @test isapprox(trapz(res.pdf, res.t), 1.0; atol = 0.05)

        # boltzmann_ic: normalised, zero on the absorbing node
        p = boltzmann_ic(prob)
        @test isapprox(sum(p), 1.0; atol = 1e-12)
        @test p[end] == 0.0

        # both ends absorbing: everything is eventually absorbed
        prob2 = FPTProblem(x, U; left_bc = :absorbing, right_bc = :absorbing)
        res2 = solve_fpt(prob2)
        @test res2.mfpt > 0
        @test all(res2.rates .> 0)
        @test isapprox(trapz(res2.pdf, res2.t), 1.0; atol = 0.05)
    end

    # -- 5. Empirical Langevin (Brownian dynamics) in the same potential ------
    # A tilted (constant-force) potential keeps escape fast and the test cheap
    # and non-flaky; it still exercises drift + reflecting + absorbing together.
    # Comparing against the full N=200 Rouse sim would instead conflate the
    # reduced effective-D approximation, so we sample the identical 1D dynamics.
    @testset "EmpiricalLangevin" begin
        L = 2.0; D = 1.0; g = 1.5; kBT = 1.0; N = 801
        Ufun(x) = g * x
        Uprime(x) = g
        x = collect(range(0.0, L; length = N))
        prob = FPTProblem(x, Ufun.(x); D = D, kBT = kBT,
                          left_bc = :reflecting, right_bc = :absorbing)
        res = solve_fpt(prob)

        rng = MersenneTwister(1234)
        β = 1.0 / kBT
        dt = 1e-3
        c1 = D * β * dt
        c2 = sqrt(2 * D * dt)
        nsamples = 8000
        maxstep = 200_000
        times = Float64[]
        for _ in 1:nsamples
            # Boltzmann initial position on [0, L] by rejection (weight peaks at 0)
            x0 = 0.0
            while true
                xp = L * rand(rng)
                if rand(rng) < exp(-β * (Ufun(xp) - Ufun(0.0)))
                    x0 = xp; break
                end
            end
            xc = x0
            for step in 1:maxstep
                xc += -c1 * Uprime(xc) + c2 * randn(rng)
                xc <= 0.0 && (xc = -xc)            # reflect at x = 0
                if xc >= L                          # absorb at x = L
                    push!(times, step * dt); break
                end
            end
        end
        @test length(times) > 0.99 * nsamples      # essentially all absorbed
        @test isapprox(mean(times), res.mfpt; rtol = 0.08)
    end

end
