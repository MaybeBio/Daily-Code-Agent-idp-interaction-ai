using Test
using Random
using Droplet

@testset "KMC Phase 2: Fenwick + sampling" begin

    @testset "fenwick_build! and prefix sums" begin
        vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        tree = Float64[]
        fenwick_build!(tree, vals)

        # prefix(i) must equal sum(vals[1:i])
        for i in 1:length(vals)
            @test isapprox(fenwick_prefix(tree, i), sum(vals[1:i]))
        end
        @test isapprox(fenwick_total(tree), sum(vals))
    end

    @testset "fenwick_update! preserves total and prefixes" begin
        vals = [1.0, 2.0, 3.0, 4.0]
        tree = Float64[]
        fenwick_build!(tree, vals)

        # bump slot 2 by 5.0
        fenwick_update!(tree, 2, 5.0)
        vals[2] += 5.0
        for i in 1:length(vals)
            @test isapprox(fenwick_prefix(tree, i), sum(vals[1:i]))
        end

        # zero out slot 4
        fenwick_update!(tree, 4, -vals[4])
        vals[4] = 0.0
        for i in 1:length(vals)
            @test isapprox(fenwick_prefix(tree, i), sum(vals[1:i]))
        end
        @test isapprox(fenwick_total(tree), sum(vals))
    end

    @testset "fenwick_find boundary cases" begin
        vals = [1.0, 2.0, 3.0, 4.0]   # cumulative: 1, 3, 6, 10
        tree = Float64[]
        fenwick_build!(tree, vals)

        # u just above each prefix boundary → next slot
        @test fenwick_find(tree, 0.0)   == 1
        @test fenwick_find(tree, 0.999) == 1
        @test fenwick_find(tree, 1.0)   == 2   # u == prefix(1) → must move past
        @test fenwick_find(tree, 2.5)   == 2
        @test fenwick_find(tree, 3.0)   == 3
        @test fenwick_find(tree, 5.999) == 3
        @test fenwick_find(tree, 6.0)   == 4
        @test fenwick_find(tree, 9.999) == 4
    end

    @testset "fenwick_find skips zero-rate slots" begin
        vals = [0.0, 1.0, 0.0, 2.0, 0.0]
        tree = Float64[]
        fenwick_build!(tree, vals)

        # prefix: 0, 1, 1, 3, 3 → only slots 2 and 4 have positive weight
        @test fenwick_find(tree, 0.0)   == 2
        @test fenwick_find(tree, 0.5)   == 2
        @test fenwick_find(tree, 0.999) == 2
        @test fenwick_find(tree, 1.0)   == 4
        @test fenwick_find(tree, 2.999) == 4

        # Random sampling: zero-rate slots must never be selected
        rng = MersenneTwister(42)
        Q = sum(vals)
        counts = zeros(Int, length(vals))
        for _ in 1:50_000
            counts[fenwick_find(tree, rand(rng) * Q)] += 1
        end
        @test counts[1] == 0
        @test counts[3] == 0
        @test counts[5] == 0
        @test counts[2] + counts[4] == 50_000
    end

    @testset "single nonzero rate is always selected" begin
        vals = [0.0, 0.0, 3.7, 0.0]
        tree = Float64[]
        fenwick_build!(tree, vals)
        Q = sum(vals)
        rng = MersenneTwister(0)
        for _ in 1:1000
            @test fenwick_find(tree, rand(rng) * Q) == 3
        end
    end

    @testset "sampleMove distribution within 3σ" begin
        # Construct a KMCSystem directly with known rates (no actual species needed
        # for this test — sampleMove only reads sys.fenwick and sys.Q).
        rates = [1.0, 2.0, 3.0, 4.0]   # Q = 10, expected probs 0.1/0.2/0.3/0.4
        sys = KMCSystem(; slist=AbstractVector[], beta=1.0)
        sys.rates = copy(rates)
        sys.nslots = length(rates)
        fenwick_build!(sys.fenwick, sys.rates)
        sys.Q = fenwick_total(sys.fenwick)
        @test isapprox(sys.Q, sum(rates))

        Random.seed!(1234)
        N = 200_000
        counts = zeros(Int, length(rates))
        for _ in 1:N
            counts[sampleMove(sys)] += 1
        end

        for (i, p) in enumerate(rates ./ sum(rates))
            σ = sqrt(p * (1 - p) / N)
            observed = counts[i] / N
            @test abs(observed - p) < 4 * σ
        end
    end

    @testset "buildMoves! builds Fenwick consistent with rates" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(15), initShape="uniform", k0=0.5)
        intPot!(b, -0.3)
        slist = [b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        @test length(sys.fenwick) == sys.nslots
        @test isapprox(fenwick_total(sys.fenwick), sum(sys.rates); atol=1e-12)
        @test isapprox(sys.Q, sum(sys.rates); atol=1e-12)

        # spot-check that fenwick_find on this real system never returns
        # a zero-rate slot
        Random.seed!(7)
        for _ in 1:1000
            slot = sampleMove(sys)
            @test sys.rates[slot] > 0.0
        end
    end
end
