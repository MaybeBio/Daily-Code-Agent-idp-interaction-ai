using Test
using Random
using Statistics
using Droplet

@testset "KMC Phase 4: main loop + sim-time sampling" begin

    @testset "kmc_step! advances n_events and t_sim" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(5), initShape="uniform", k0=1.0)
        intPot!(b, -0.2)
        slist = [b]
        init!(l, slist)

        Random.seed!(7)
        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        for k in 1:100
            t_before = sys.t_sim
            kmc_step!(sys)
            @test sys.n_events == k
            @test sys.t_sim > t_before
        end
    end

    @testset "kmc_run! respects t_max_sim cap" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(10), initShape="uniform", k0=1.0)
        intPot!(b, -0.3)
        slist = [b]
        init!(l, slist)

        Random.seed!(11)
        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        kmc_run!(sys; t_max_sim=5.0)

        # We stop the FIRST time t_sim crosses t_max_sim, so sys.t_sim is
        # somewhere in [t_max_sim, t_max_sim + Δt_typical]. Check it overshot
        # by at most a few mean event spacings (1/Q ≪ 1).
        @test sys.t_sim >= 5.0
        @test sys.t_sim < 5.5
        @test sys.n_events > 0
    end

    @testset "kmc_run! callback fires at sim-time grid" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(8), initShape="uniform", k0=1.0)
        intPot!(b, -0.3)
        slist = [b]
        init!(l, slist)

        Random.seed!(13)
        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        sample_times = Float64[]
        cb = s -> push!(sample_times, s.t_sim)

        n = kmc_run!(sys; t_max_sim=5.5, dt_sample=1.0, sample_callback=cb)

        # expected grid: 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 → 6 samples (the 5.0 sample
        # fires when t_sim first crosses 5.0, then the loop terminates next
        # iteration when t_sim ≥ 5.5).
        @test n == length(sample_times)
        @test length(sample_times) >= 5
        @test length(sample_times) <= 6

        # each recorded t_sim is at or just past the expected grid point
        for (i, t) in enumerate(sample_times)
            grid = (i - 1) * 1.0
            @test t >= grid
            @test t < grid + 0.5  # 1/Q ≪ 0.5 for this system
        end
    end

    @testset "free-particle MSD ≈ 6·k0·t" begin
        n_walkers = 500
        k0 = 1.0
        t_target = 8.0
        boxlength = Int32(80)
        center = boxlength ÷ Int32(2)

        Random.seed!(20260503)
        ratios = Float64[]   # (Δr² / (6 k0 t_actual)) per walker, expect ≈ 1

        for _ in 1:n_walkers
            l = Lattice(; boxlength=boxlength, nspecies=Int32(1))
            b = BulkMonomer(; length=Int32(1), initShape="line",
                              initPos=Int32[center, center, center], k0=k0)
            intPot!(b, 0.0)   # non-interacting
            slist = [b]
            init!(l, slist)

            sys = KMCSystem(; slist=slist, beta=1.0)
            buildMoves!(sys)
            kmc_run!(sys; t_max_sim=t_target)

            p = b.pos[Int32(1)]
            dx = Int(p[1]) - Int(center)
            dy = Int(p[2]) - Int(center)
            dz = Int(p[3]) - Int(center)
            # box is large enough that no walker wraps; sanity check
            @test abs(dx) < Int(boxlength) ÷ 2
            @test abs(dy) < Int(boxlength) ÷ 2
            @test abs(dz) < Int(boxlength) ÷ 2

            r2 = dx*dx + dy*dy + dz*dz
            push!(ratios, r2 / (6.0 * k0 * sys.t_sim))
        end

        μ = mean(ratios)
        # ~3σ tolerance: σ(r²) ≈ √(8/3) · 6k0 t for 3D RW; σ(mean) ≈ that / √n
        @test 0.85 < μ < 1.15
    end

    @testset "Metropolis–KMC equilibrium parity" begin
        # Same physical setup; both sample π ∝ exp(-βH) with Metropolis-rate KMC.
        boxlen = Int32(6)
        npart  = Int32(8)
        ε      = -0.6
        β      = 1.0

        function setup_system()
            l = Lattice(; boxlength=boxlen, nspecies=Int32(1))
            b = BulkMonomer(; length=npart, initShape="uniform", k0=1.0)
            intPot!(b, ε)
            slist = [b]
            init!(l, slist)
            return b, slist
        end

        # ---- Metropolis trajectory ----
        Random.seed!(42)
        bM, _ = setup_system()
        ham!(bM)
        n_burn_metro = 20_000
        n_collect_metro = 200_000
        sample_every_metro = 200

        accepted = 0
        while accepted < n_burn_metro
            move!(bM) && (accepted += 1)
        end

        hamsM = Float64[]
        accepted = 0
        attempts = 0
        while accepted < n_collect_metro
            if move!(bM)
                accepted += 1
                if accepted % sample_every_metro == 0
                    ham!(bM)
                    push!(hamsM, bM.ham)
                end
            end
            attempts += 1
        end

        # ---- KMC trajectory ----
        Random.seed!(43)
        bK, slistK = setup_system()
        sysK = KMCSystem(; slist=slistK, beta=β)
        buildMoves!(sysK)

        for _ in 1:20_000
            kmc_step!(sysK)
        end

        hamsK = Float64[]
        for k in 1:200_000
            kmc_step!(sysK)
            if k % 200 == 0
                ham!(bK)
                push!(hamsK, bK.ham)
            end
        end

        μM, σM = mean(hamsM), std(hamsM)
        μK, σK = mean(hamsK), std(hamsK)

        # Means should agree within MC noise. Autocorrelation inflates the
        # naïve std-of-mean; allow a generous absolute tolerance.
        @test abs(μM - μK) < 0.15
        @test abs(σM - σK) < 0.20
        # Sanity: both should be in a physically reasonable range
        @test μM < 0.0
        @test μK < 0.0
    end
end
