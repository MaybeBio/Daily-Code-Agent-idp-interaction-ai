using Test
using Droplet

@testset "KMC Phase 1: rate primitives" begin

    @testset "isolated particle: all 6 rates equal k0" begin
        l = Lattice(; boxlength=Int32(4), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(1), initShape="line",
                          initPos=Int32[1,1,1], k0=2.5)
        intPot!(b, -1.0)
        slist = [b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        @test sys.nslots == 6
        @test length(sys.rates) == 6
        for r in sys.rates
            @test isapprox(r, b.k0)
        end
        @test isapprox(sys.Q, 6 * b.k0)
        @test isapprox(sys.Q, sum(sys.rates))
    end

    @testset "two adjacent particles: blocking + bond-break rate" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(2), initShape="line",
                          initPos=Int32[3,3,3], k0=1.0)
        ε = -0.7
        intPot!(b, ε)
        slist = [b]
        init!(l, slist)

        # geometry sanity
        @test b.pos[Int32(1)][1:3] == Int32[3,3,3]
        @test b.pos[Int32(2)][1:3] == Int32[4,3,3]

        β = 1.5
        sys = KMCSystem(; slist=slist, beta=β)
        buildMoves!(sys)

        @test sys.nslots == 12

        # Slot layout: particle p, dir d → slot = (p-1)*6 + d
        # dir: 1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z
        slot(p, d) = (p - 1) * 6 + d

        # blocked: particle 1 hopping +x onto particle 2,
        # and particle 2 hopping -x onto particle 1
        @test sys.rates[slot(1, 1)] == 0.0
        @test sys.rates[slot(2, 2)] == 0.0

        # All other 10 hops break the single bond (ΔH = -ε > 0).
        # Δnnp_self = (1 - 1) - 1 = -1 → ΔH = ε * (-1) = -ε
        expected = 1.0 * min(1.0, exp(-β * (-ε)))
        for (p, d) in [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                       (2, 1), (2, 3), (2, 4), (2, 5), (2, 6)]
            @test isapprox(sys.rates[slot(p, d)], expected)
        end

        @test isapprox(sys.Q, 10 * expected)
        @test isapprox(sys.Q, sum(sys.rates))
    end

    @testset "Q == sum(rates) on a denser config" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(20), initShape="uniform", k0=0.7)
        intPot!(b, -0.4)
        slist = [b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=0.9)
        buildMoves!(sys)

        @test sys.nslots == 20 * 6
        @test isapprox(sys.Q, sum(sys.rates); atol=1e-12)
        @test all(r -> r >= 0.0, sys.rates)
        @test all(r -> r <= b.k0 + 1e-12, sys.rates)   # Metropolis-rate bound
    end
end
