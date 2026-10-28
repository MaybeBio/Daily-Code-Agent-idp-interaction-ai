using Test
using Random
using Statistics
using Droplet

# v2.4 gating tests: equilibrium correctness for KMC polymer dynamics.
#
# IMPORTANT: the existing Metropolis polymer sampler in src/move.jl uses an
# asymmetric proposal (proposeMove!: Padd(N)=(1-Pkink)(N+2)/(2N+2),
# Premove(N)=(1-Pkink)·N/(2N+2)) WITHOUT a corresponding Hastings q-ratio in
# acceptMove(p). For DB the equilibrium of that chain is
#   π_M(s')/π_M(s) = (q_fwd/q_rev) · exp(-βΔH)
#                  = (N+3)/(N-1) · nRemove_Y/nAdd_X · exp(-βΔH)
# which (with typical nRemove≈1, nAdd≈4) underweights longer states by
# ~(N+3)/(4(N-1)) ≈ 0.36 per add-step relative to π ∝ exp(-βH). So strict
# parity with Metropolis (as v2's Phase 4/5 plan calls for) is NOT possible
# for polymers — Metropolis-rate KMC samples Boltzmann correctly, and the
# Metropolis chain samples a biased distribution.
#
# This file therefore replaces strict parity with three checks:
#   1. KMC stationarity: distribution moments are stable across the second
#      half of the trajectory (system has equilibrated).
#   2. KMC self-consistency: two KMC runs from different RNG seeds produce
#      statistically indistinguishable equilibrium distributions.
#   3. Loose Metropolis comparison: both samplers produce stable, qualitatively
#      similar distributions in the same physical regime — useful as a sanity
#      check that nothing is catastrophically wrong, not as strict parity.

function _run_kmc(slist, β::Float64, n_burn::Int, n_collect::Int,
                  sample_every::Int, p::Polymer, b::Union{Nothing, BulkMonomer}=nothing)
    sys = KMCSystem(; slist=slist, beta=β, refresh_mode=:fine_polymer)
    buildMoves!(sys)
    for _ in 1:n_burn
        sys.Q <= 0 && break
        kmc_step!(sys)
    end
    hamsP = Float64[]; hamsB = Float64[]; lens = Int[]
    for k in 1:n_collect
        sys.Q <= 0 && break
        kmc_step!(sys)
        if k % sample_every == 0
            ham!(p)
            push!(hamsP, p.ham)
            push!(lens, Int(p.length))
            if b !== nothing
                ham!(b)
                push!(hamsB, b.ham)
            end
        end
    end
    return hamsP, hamsB, lens
end

function _run_metropolis(slist, n_burn::Int, n_collect::Int,
                          sample_every::Int, p::Polymer, b::Union{Nothing, BulkMonomer}=nothing)
    ham!(p); b !== nothing && ham!(b)
    accepted = 0
    while accepted < n_burn
        s = rand(slist)
        move!(s) && (accepted += 1)
    end
    hamsP = Float64[]; hamsB = Float64[]; lens = Int[]
    accepted = 0
    while accepted < n_collect
        s = rand(slist)
        if move!(s)
            accepted += 1
            if accepted % sample_every == 0
                ham!(p)
                push!(hamsP, p.ham)
                push!(lens, Int(p.length))
                if b !== nothing
                    ham!(b)
                    push!(hamsB, b.ham)
                end
            end
        end
    end
    return hamsP, hamsB, lens
end

@testset "KMC Polymer v2.4: equilibrium correctness" begin

    @testset "polymer-only: KMC stationarity & self-consistency" begin
        boxlen   = Int32(16)
        maxlen   = Int32(24)
        L0       = Int32(8)
        chemPot  = -2.5
        intPot_v = -0.4
        β        = 1.0

        function setup(seed::Int)
            Random.seed!(seed)
            l = Lattice(; boxlength=boxlen, nspecies=Int32(1))
            p = Polymer(; length=L0, maxlength=maxlen,
                          nMonomerTypes=Int32(1),
                          initShape="line",
                          initPos=Int32[1, boxlen÷2, boxlen÷2],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=chemPot, defaultIntPot=intPot_v)
            slist = [p]
            init!(l, slist)
            return p, slist
        end

        # Two independent KMC runs from different seeds: same equilibrium.
        Random.seed!(20260504)
        p1, slist1 = setup(20260504)
        h1, _, l1 = _run_kmc(slist1, β, 30_000, 200_000, 200, p1)

        Random.seed!(20260505)
        p2, slist2 = setup(20260505)
        h2, _, l2 = _run_kmc(slist2, β, 30_000, 200_000, 200, p2)

        # Self-consistency: two seeds agree within MC noise. The standard
        # error of the mean for autocorrelated samples ≈ σ·√(τ_int/N); we
        # use a conservative absolute tolerance below.
        @test abs(mean(h1) - mean(h2)) < 0.5
        @test abs(std(h1)  - std(h2))  < 0.5
        @test abs(mean(l1) - mean(l2)) < 0.5
        @test abs(std(l1)  - std(l2))  < 0.4

        # Stationarity: split each run in half, halves agree.
        for (h, l) in ((h1, l1), (h2, l2))
            n = length(h) ÷ 2
            @test abs(mean(h[1:n])    - mean(h[n+1:end])) < 0.7
            @test abs(mean(l[1:n])    - mean(l[n+1:end])) < 0.6
        end

        # Sanity: chain isn't pinned to a wall in either run.
        for l in (l1, l2)
            @test 5.0 < mean(l) < Float64(maxlen) - 1.0
            @test std(l) > 0.5
        end
    end

    @testset "polymer + bulk cross-interaction: KMC stationarity & self-consistency" begin
        boxlen   = Int32(16)
        maxlen   = Int32(24)
        L0       = Int32(8)
        chemPot  = -2.5
        intPot_v = -0.4
        nbulk    = Int32(30)
        bulk_int = -0.2
        cross    = -0.5
        β        = 1.0

        function setup(seed::Int)
            Random.seed!(seed)
            l = Lattice(; boxlength=boxlen, nspecies=Int32(2))
            p = Polymer(; length=L0, maxlength=maxlen,
                          nMonomerTypes=Int32(1),
                          initShape="line",
                          initPos=Int32[1, boxlen÷2, boxlen÷2],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=chemPot, defaultIntPot=intPot_v)
            b = BulkMonomer(; length=nbulk, initShape="uniform", k0=1.0)
            intPot!(b, bulk_int)
            intPot!(p.monomers[1], b, cross)
            slist = [p, b]
            init!(l, slist)
            return p, b, slist
        end

        # The cross-interaction couples chain growth/folding to bulk
        # crowding; autocorrelation here is noticeably longer than the
        # polymer-only case. Larger burn-in + collection budget is needed.
        p1, b1, slist1 = setup(20260506)
        h1P, h1B, l1 = _run_kmc(slist1, β, 80_000, 400_000, 300, p1, b1)

        p2, b2, slist2 = setup(20260507)
        h2P, h2B, l2 = _run_kmc(slist2, β, 80_000, 400_000, 300, p2, b2)

        # Two-seed self-consistency. Tolerances are wider than the
        # polymer-only case to absorb the longer autocorrelation; the
        # standard error of the mean for polymer ham at this sample size is
        # ≈1.0 once integrated autocorrelation is accounted for.
        @test abs(mean(h1P) - mean(h2P)) < 1.5
        @test abs(std(h1P)  - std(h2P))  < 0.8
        @test abs(mean(h1B) - mean(h2B)) < 0.4
        @test abs(std(h1B)  - std(h2B))  < 0.4
        @test abs(mean(l1)  - mean(l2))  < 0.9
        @test abs(std(l1)   - std(l2))   < 0.5

        # Cross-interaction is attractive — bulk H should be negative on
        # average (favorable contacts) for both runs, and within similar range.
        @test mean(h1B) < 0.0
        @test mean(h2B) < 0.0

        # Stationarity check on the polymer length series.
        for l in (l1, l2)
            n = length(l) ÷ 2
            @test abs(mean(l[1:n]) - mean(l[n+1:end])) < 0.9
        end
    end

    @testset "loose Metropolis comparison (documents proposal bias)" begin
        # Metropolis polymer's proposeMove! has an asymmetric proposal density
        # (Padd ∝ (N+2)/(2N+2), Premove ∝ N/(2N+2)) without a Hastings q-ratio
        # in acceptMove, so π_M ≠ π ∝ exp(-βH). We expect KMC and Metropolis
        # to land in the SAME REGIME but not at the same moments. This testset
        # records both and asserts only loose qualitative agreement (within
        # ~5 length units / ~7 in ham), to flag if either sampler ever
        # diverges catastrophically.
        boxlen   = Int32(16)
        maxlen   = Int32(24)
        L0       = Int32(8)
        chemPot  = -2.5
        intPot_v = -0.4
        β        = 1.0

        function setup(seed::Int)
            Random.seed!(seed)
            l = Lattice(; boxlength=boxlen, nspecies=Int32(1))
            p = Polymer(; length=L0, maxlength=maxlen,
                          nMonomerTypes=Int32(1),
                          initShape="line",
                          initPos=Int32[1, boxlen÷2, boxlen÷2],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=chemPot, defaultIntPot=intPot_v)
            slist = [p]
            init!(l, slist)
            return p, slist
        end

        # Metropolis run
        pM, slistM = setup(20260508)
        hM, _, lM = _run_metropolis(slistM, 30_000, 100_000, 200, pM)

        # KMC run
        pK, slistK = setup(20260509)
        hK, _, lK = _run_kmc(slistK, β, 30_000, 100_000, 200, pK)

        # Loose tolerances reflecting the known proposal bias (~3 length
        # units in this regime). If these fail, something has changed
        # algorithmically.
        @test abs(mean(hM) - mean(hK)) < 8.0
        @test abs(mean(lM) - mean(lK)) < 5.0
        @test abs(std(hM)  - std(hK))  < 5.0
        @test abs(std(lM)  - std(lK))  < 3.0

        # Both should land in physically reasonable regimes.
        @test 5.0 < mean(lM) < Float64(maxlen) - 1.0
        @test 5.0 < mean(lK) < Float64(maxlen) - 1.0
    end
end
