using Test
using Random
using Droplet

# Recompute every rate from scratch into a local vector without touching sys
# state. Used as the ground-truth reference in the parity test below.
function _ground_truth_rates(sys)
    r = zeros(Float64, sys.nslots)
    for slot in 1:sys.nslots
        r[slot] = computeRate!(sys, slot)
    end
    return r, sum(r)
end

@testset "KMC Phase 3: applyMove! + incremental updates" begin

    @testset "incremental Q and rates match full rebuild over 1000 events" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(30), initShape="uniform", k0=1.0)
        intPot!(b, -0.4)
        slist = [b]
        init!(l, slist)

        Random.seed!(20260503)
        sys = KMCSystem(; slist=slist, beta=0.9)
        buildMoves!(sys)

        for step in 1:1000
            if sys.Q <= 0.0
                break
            end
            slot = sampleMove(sys)
            applyMove!(sys, slot)

            r_truth, Q_truth = _ground_truth_rates(sys)

            @test isapprox(sys.Q, Q_truth; rtol=1e-10, atol=1e-12)
            @test isapprox(fenwick_total(sys.fenwick), sys.Q; rtol=1e-12, atol=1e-12)
            @test all(isapprox.(sys.rates, r_truth; rtol=1e-10, atol=1e-12))
        end

        @test sys.n_events > 0
        @test sys.t_sim > 0.0
    end

    @testset "isolated particle: Q invariant under hops" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(1), initShape="line",
                          initPos=Int32[4,4,4], k0=1.7)
        intPot!(b, -0.5)
        slist = [b]
        init!(l, slist)

        Random.seed!(1)
        sys = KMCSystem(; slist=slist, beta=2.0)
        buildMoves!(sys)
        @test isapprox(sys.Q, 6 * b.k0)

        for _ in 1:200
            slot = sampleMove(sys)
            applyMove!(sys, slot)
            @test isapprox(sys.Q, 6 * b.k0)
            @test all(r -> isapprox(r, b.k0), sys.rates)
        end
    end

    @testset "two adjacent particles: Q transitions correctly under a -x hop" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(2), initShape="line",
                          initPos=Int32[3,3,3], k0=1.0)
        ε = -0.7
        intPot!(b, ε)
        slist = [b]
        init!(l, slist)

        β = 1.5
        sys = KMCSystem(; slist=slist, beta=β)
        buildMoves!(sys)

        # Slot for particle 1 hopping -x (from [3,3,3] to [2,3,3]).
        # dir 2 = -x, particle 1 → slot 2.
        bond_break_rate = min(1.0, exp(-β * (-ε)))
        @test isapprox(sys.Q, 10 * bond_break_rate)
        @test isapprox(sys.rates[2], bond_break_rate)

        applyMove!(sys, 2)

        # After the hop: particle 1 at [2,3,3], particle 2 at [4,3,3]; they are
        # distance 2, no bond. All 12 hops are now valid:
        #   particle 1 +x → [3,3,3] empty: forms bond with particle 2 → ΔH=ε,
        #     Metropolis cap → rate 1.0
        #   particle 2 -x → [3,3,3] empty: forms bond with particle 1 → rate 1.0
        #   the 10 other hops are neutral (both src and dst neighbor counts of
        #     self differ by 0 after the -1 self-correction) → rate 1.0
        @test isapprox(sys.Q, 12.0; atol=1e-12)
        @test all(r -> isapprox(r, 1.0), sys.rates)

        # particle 1 should now be at [2,3,3]
        @test b.pos[Int32(1)][1:3] == Int32[2,3,3]
        @test b.latS[2,3,3] == Int32(1)
        @test b.latS[3,3,3] == Int32(0)

        # parity with full rebuild
        r_truth, Q_truth = _ground_truth_rates(sys)
        @test isapprox(sys.Q, Q_truth; rtol=1e-12)
        @test all(isapprox.(sys.rates, r_truth))
    end

    @testset "rebuildQ! produces identical state to incremental updates" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(20), initShape="uniform", k0=0.5)
        intPot!(b, -0.3)
        slist = [b]
        init!(l, slist)

        Random.seed!(99)
        sys = KMCSystem(; slist=slist, beta=1.2)
        buildMoves!(sys)

        for _ in 1:200
            slot = sampleMove(sys)
            applyMove!(sys, slot)
        end

        # snapshot incremental state
        rates_inc = copy(sys.rates)
        Q_inc = sys.Q

        rebuildQ!(sys)

        @test all(isapprox.(rates_inc, sys.rates; rtol=1e-10, atol=1e-12))
        @test isapprox(Q_inc, sys.Q; rtol=1e-10, atol=1e-12)
    end

    @testset "n_events and t_sim advance monotonically" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(5), initShape="uniform", k0=1.0)
        intPot!(b, -0.2)
        slist = [b]
        init!(l, slist)

        Random.seed!(123)
        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        t_prev = 0.0
        for k in 1:50
            slot = sampleMove(sys)
            applyMove!(sys, slot)
            @test sys.n_events == k
            @test sys.t_sim > t_prev
            t_prev = sys.t_sim
        end
    end
end
