using Droplet
using Test
using Random
using Statistics
using LinearAlgebra

# Re-import internals that are not exported.
import Droplet: _build_adjacency_matrix

@testset verbose=true "Rouse polymer Langevin solver" begin

    @testset "adjacency matrix" begin
        A = _build_adjacency_matrix(10)
        @test size(A) == (10, 10)
        @test A == A'
        @test diag(A) == [1.0; fill(2.0, 8); 1.0]
        @test diag(A, 1) == fill(-1.0, 9)
        @test diag(A, -1) == fill(-1.0, 9)
        @test all(abs.(sum(A; dims = 2)) .< 1e-12)   # zero-mode: row sums vanish
        @test A isa SymTridiagonal                    # O(N) mul! in the hot loop
    end

    @testset "tridiagonal drift matches the dense product" begin
        # Dense BLAS sums the three nonzero terms of a row in a kernel-dependent
        # order, so the agreement is to rounding, not bitwise.
        rng = MersenneTwister(3)
        for N in (2, 3, 4, 11, 111, 301)
            A = _build_adjacency_matrix(N)
            x = randn(rng, N, 3)
            @test isapprox(mul!(similar(x), A, x), Matrix(A) * x; rtol = 1e-14)
        end
    end

    @testset "scaling parameters" begin
        r = RousePolymer(; boltzmannConst = 1.0, temperature = 1.0,
                           friction = 1.0, kuhnLength = 1.0,
                           dt = 1e-3, nbeads = Int32(5), nsteps = Int32(2),
                           trackedBead1 = Int32(1), trackedBead2 = Int32(5))
        init!(r)
        @test r.stiffness  ≈ 3.0
        @test r.timeScale  ≈ 1 / 3
        @test r.posScale   ≈ sqrt(2 / 3)
        @test r.dtDimless  ≈ r.dt / r.timeScale
    end

    @testset "deterministic mode-1 decay (no noise)" begin
        function decay_test()
            N = 20
            nsteps = 200
            r = RousePolymer(; nbeads = Int32(N), nsteps = Int32(nsteps),
                               dt = 1e-3, saveMode = :full, seed = 1,
                               trackedBead1 = Int32(1), trackedBead2 = Int32(N))
            init!(r)
            # Eigenvalues of A: λ_p = 4 sin²(p π / (2N)); eigenvectors
            # v_p(i) ∝ cos(p π (i - 1/2) / N), p = 0, 1, …, N-1.
            v1 = [cos(π * (i - 0.5) / N) for i in 1:N]
            amp0 = 1.5
            r.pos .= 0.0
            r.pos[:, 1] .= amp0 .* v1

            tau = r.dtDimless
            for _ in 1:nsteps
                mul!(r.drift, r.adjacency, r.pos)
                @. r.pos = r.pos - r.drift * tau
            end

            λ1 = 4 * sin(π / (2N))^2
            expected = amp0 * exp(-λ1 * tau * nsteps)
            observed = mean(r.pos[:, 1] .* v1) / mean(v1 .^ 2)
            return abs(observed - expected) / expected
        end
        @test decay_test() < 5e-2
    end

    @testset "equilibrium bond variance ⟨|Δr|²⟩ = b²" begin
        function bond_variance_test()
            Random.seed!(42)
            N = 30
            # dtDimless = 0.01 → small EM bias on the fastest mode.
            nsteps = 50_000
            r = RousePolymer(; nbeads = Int32(N), nsteps = Int32(nsteps),
                               dt = 0.01 / 3, saveMode = :full, seed = 7,
                               trackedBead1 = Int32(1), trackedBead2 = Int32(N))
            init!(r); simulate!(r)

            burn = div(nsteps, 5)
            acc = 0.0
            cnt = 0
            @inbounds for t in (burn + 1):nsteps
                for i in 2:(N - 1)
                    dx = r.posTraj[t, i, 1] - r.posTraj[t, i - 1, 1]
                    dy = r.posTraj[t, i, 2] - r.posTraj[t, i - 1, 2]
                    dz = r.posTraj[t, i, 3] - r.posTraj[t, i - 1, 3]
                    acc += dx * dx + dy * dy + dz * dz
                    cnt += 1
                end
            end
            return acc / cnt, r.kuhnLength^2
        end
        bvar, expected = bond_variance_test()
        @test isapprox(bvar, expected; rtol = 0.08)
    end

    @testset "CoM MSD = 6 D t / N" begin
        # Wrap in `function` (not `let`) so hot loops are properly type-inferred —
        # top-level testset scope is global, which kills loop performance.
        function com_msd_test()
            Random.seed!(123)
            N = 20
            nsteps = 20_000
            r = RousePolymer(; nbeads = Int32(N), nsteps = Int32(nsteps),
                               dt = 0.01 / 3, saveMode = :full, seed = 11,
                               trackedBead1 = Int32(1), trackedBead2 = Int32(N))
            init!(r)
            simulate!(r)

            comTraj = zeros(Float64, nsteps, 3)
            @inbounds for t in 1:nsteps
                for a in 1:3
                    s = 0.0
                    for i in 1:N
                        s += r.posTraj[t, i, a]
                    end
                    comTraj[t, a] = s / N
                end
            end

            maxLag = div(nsteps, 10)
            msd  = zeros(Float64, maxLag)
            @inbounds for tau in 1:maxLag
                acc = 0.0
                n = nsteps - tau
                for t in 1:n
                    dx = comTraj[t + tau, 1] - comTraj[t, 1]
                    dy = comTraj[t + tau, 2] - comTraj[t, 2]
                    dz = comTraj[t + tau, 3] - comTraj[t, 3]
                    acc += dx * dx + dy * dy + dz * dz
                end
                msd[tau] = acc / n
            end

            sel = div(maxLag, 10):div(maxLag, 2)
            tphys = collect(sel) .* r.dt
            slope = sum(tphys .* msd[sel]) / sum(tphys .^ 2)
            D = r.boltzmannConst * r.temperature / r.friction
            return slope, 6 * D / N
        end
        slope, expected = com_msd_test()
        @test isapprox(slope, expected; rtol = 0.15)
    end

    @testset "save modes produce identical distance trajectories" begin
        kwargs = (; nbeads = Int32(15), nsteps = Int32(500),
                   dt = 1e-3, seed = 31, trackedBead1 = Int32(1),
                   trackedBead2 = Int32(15))
        r1 = RousePolymer(; kwargs..., saveMode = :full)
        r2 = RousePolymer(; kwargs..., saveMode = :distance_only)
        init!(r1); simulate!(r1)
        init!(r2); simulate!(r2)
        @test r1.distanceTraj == r2.distanceTraj
        @test r2.posTraj === nothing
    end

    @testset "stats utilities" begin
        # contactTimes — hand-crafted trajectory with 3 distinct contacts.
        # Crossings occur at t=2,5,9 (1-based step indices) since
        # the rising-edge detector triggers when d[i-1] > thr and d[i] ≤ thr.
        traj = [2.0, 0.5, 0.4, 1.5, 0.3, 1.2, 2.0, 1.8, 0.2, 0.1, 1.0]
        thr  = 0.6
        dt   = 0.1
        times, invCdf = contactTimes(traj, thr, dt)
        @test length(times) == 3
        @test times ≈ [0.1, 0.4, 0.8]
        @test length(invCdf) == 2
        @test invCdf == sort(diff(times))

        # distanceHistogram — total counts must equal trajectory length.
        Random.seed!(99)
        d = rand(2000)
        centers, counts = distanceHistogram(d, 25)
        @test sum(counts) == length(d)
        @test length(centers) == 25

        # beadMSD of a 3D Gaussian random walk: slope = 3 (per-step variance = 1).
        Random.seed!(5)
        nsteps = 40_000
        pos = zeros(Float64, nsteps, 1, 3)
        @inbounds for t in 2:nsteps
            for a in 1:3
                pos[t, 1, a] = pos[t - 1, 1, a] + randn()
            end
        end
        msd = beadMSD(pos, 1; maxLag = 1000)
        # Fit slope on the central range to avoid endpoint noise.
        sel = 50:900
        slope = sum((sel) .* msd[sel]) / sum((sel) .^ 2)
        @test isapprox(slope, 3.0; rtol = 0.10)

        # beadPairDistanceMSD — sanity check that it returns the right length.
        pos2 = zeros(Float64, nsteps, 2, 3)
        @inbounds for t in 2:nsteps
            for a in 1:3
                pos2[t, 1, a] = pos2[t - 1, 1, a] + randn()
                pos2[t, 2, a] = pos2[t - 1, 2, a] + randn()
            end
        end
        pmsd = beadPairDistanceMSD(pos2, 1, 2; maxLag = 200)
        @test length(pmsd) == 200
        @test all(pmsd .>= 0)
    end

    @testset "attractive force: nothing ≡ no-force fallback" begin
        # attractiveForce=nothing must reproduce the exact same trajectory
        # as not setting the field at all (and not perturb the RNG stream).
        kwargs = (; nbeads = Int32(10), nsteps = Int32(2000), dt = 1e-3,
                   seed = 17, trackedBead1 = Int32(1), trackedBead2 = Int32(10),
                   saveMode = :distance_only)
        r1 = RousePolymer(; kwargs...)
        r2 = RousePolymer(; kwargs..., attractiveForce = nothing)
        init!(r1); simulate!(r1)
        init!(r2); simulate!(r2)
        @test r1.distanceTraj == r2.distanceTraj
    end

    @testset "attractive harmonic tether: ⟨d²⟩ = 3kBT/(κ_chain + k_attract)" begin
        function tether_test()
            # Two-bead chain: bond stiffness κ_chain = 3 kBT/b² = 3.
            # Add a harmonic spring U = (k/2) d² between the two beads.
            # Equilibrium 3D distance variance: ⟨d²⟩ = 3 kBT / (κ_chain + k).
            k_attract = 9.0
            r = RousePolymer(;
                nbeads = Int32(2), nsteps = Int32(80_000),
                dt = 0.01 / 3, saveMode = :distance_only, seed = 23,
                trackedBead1 = Int32(1), trackedBead2 = Int32(2),
                attractiveForce = d -> -k_attract * d,
            )
            init!(r); simulate!(r)
            burn = div(r.nsteps, 5)
            d² = view(r.distanceTraj, (burn + 1):r.nsteps) .^ 2
            mean(d²)
        end
        expected = 3.0 / (3.0 + 9.0)   # = 0.25
        @test isapprox(tether_test(), expected; rtol = 0.10)
    end
end
