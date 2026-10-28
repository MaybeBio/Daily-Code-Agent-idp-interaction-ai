using Test
using Random
using Statistics
using Droplet

@testset "KMC Phase 5: multi-species cross-interactions" begin

    @testset "two species, mutual blocking + cross bond-break rate" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(2))
        a = BulkMonomer(; length=Int32(1), initShape="line",
                          initPos=Int32[3,3,3], k0=1.0)
        b = BulkMonomer(; length=Int32(1), initShape="line",
                          initPos=Int32[4,3,3], k0=1.0)
        intPot!(a, 0.0)
        intPot!(b, 0.0)
        ε_cross = -0.5
        intPot!(a, b, ε_cross)
        slist = [a, b]
        init!(l, slist)

        # geometry sanity
        @test a.pos[Int32(1)][1:3] == Int32[3,3,3]
        @test b.pos[Int32(1)][1:3] == Int32[4,3,3]

        β = 1.0
        sys = KMCSystem(; slist=slist, beta=β)
        buildMoves!(sys)

        # Per the existing Metropolis convention, distinct species use independent
        # lattice layers and may co-occupy. So A's +x hop into B's site (slot 1)
        # is NOT blocked: it just breaks the cross-bond and lands at the same
        # site as B. Same for B's -x onto A (slot 8). All 12 hops are valid and
        # cost ΔH = -ε_cross.
        @test sys.nslots == 12
        bond_break_rate = min(1.0, exp(-β * (-ε_cross)))
        for slot in 1:12
            @test isapprox(sys.rates[slot], bond_break_rate)
        end
        @test isapprox(sys.Q, 12 * bond_break_rate)
        @test isapprox(fenwick_total(sys.fenwick), sys.Q; atol=1e-12)
    end

    @testset "applyMove! across species refreshes both species' rates" begin
        l = Lattice(; boxlength=Int32(6), nspecies=Int32(2))
        a = BulkMonomer(; length=Int32(1), initShape="line",
                          initPos=Int32[3,3,3], k0=1.0)
        b = BulkMonomer(; length=Int32(1), initShape="line",
                          initPos=Int32[4,3,3], k0=1.0)
        intPot!(a, 0.0); intPot!(b, 0.0)
        intPot!(a, b, -0.5)
        slist = [a, b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        # Apply A's -x hop (slot 2): A → [2,3,3], B stays at [4,3,3], distance 2.
        applyMove!(sys, 2)

        @test a.pos[Int32(1)][1:3] == Int32[2,3,3]
        @test a.latS[2,3,3] == Int32(1)
        @test a.latS[3,3,3] == Int32(0)
        @test b.latS[4,3,3] == Int32(1)

        # After the hop:
        #   • A +x → [3,3,3] empty: forms cross bond with B (B at [4,3,3] is NN
        #     of [3,3,3]), ΔH = ε_cross = -0.5, Metropolis cap → rate 1.0
        #   • B -x → [3,3,3] empty: forms cross bond with A (A at [2,3,3] is NN
        #     of [3,3,3]), rate 1.0
        #   • All 10 other hops: ΔH = 0 → rate 1.0
        @test all(r -> isapprox(r, 1.0), sys.rates)
        @test isapprox(sys.Q, 12.0)

        # cross-check against ground-truth recompute
        truth = [computeRate!(sys, s) for s in 1:sys.nslots]
        @test all(isapprox.(sys.rates, truth))
    end

    @testset "incremental Q matches full rebuild over 1000 multi-species events" begin
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(2))
        a = BulkMonomer(; length=Int32(15), initShape="uniform", k0=1.0)
        b = BulkMonomer(; length=Int32(15), initShape="uniform", k0=1.3)
        intPot!(a, -0.2)
        intPot!(b, -0.3)
        intPot!(a, b, -0.4)
        slist = [a, b]
        init!(l, slist)

        Random.seed!(20260504)
        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        for _ in 1:1000
            sys.Q <= 0.0 && break
            applyMove!(sys, sampleMove(sys))

            truth = [computeRate!(sys, s) for s in 1:sys.nslots]
            Q_truth = sum(truth)
            @test isapprox(sys.Q, Q_truth; rtol=1e-10, atol=1e-12)
            @test isapprox(fenwick_total(sys.fenwick), sys.Q; rtol=1e-12, atol=1e-12)
            @test all(isapprox.(sys.rates, truth; rtol=1e-10, atol=1e-12))
        end
    end

    @testset "Metropolis–KMC equilibrium parity (two species, cross-bond)" begin
        boxlen = Int32(6)
        β      = 1.0

        function setup()
            l = Lattice(; boxlength=boxlen, nspecies=Int32(2))
            a = BulkMonomer(; length=Int32(4), initShape="uniform", k0=1.0)
            b = BulkMonomer(; length=Int32(4), initShape="uniform", k0=1.0)
            intPot!(a, -0.3)
            intPot!(b, -0.3)
            intPot!(a, b, -0.5)
            slist = [a, b]
            init!(l, slist)
            return a, b, slist
        end

        # ---- Metropolis ----
        Random.seed!(42)
        aM, bM, _ = setup()
        ham!(aM); ham!(bM)
        n_burn  = 30_000
        n_coll  = 200_000
        every   = 200

        species = (aM, bM)
        accepted = 0
        while accepted < n_burn
            move!(species[rand(1:2)]) && (accepted += 1)
        end

        a_hamsM = Float64[]; b_hamsM = Float64[]
        accepted = 0
        while accepted < n_coll
            if move!(species[rand(1:2)])
                accepted += 1
                if accepted % every == 0
                    ham!(aM); ham!(bM)
                    push!(a_hamsM, aM.ham); push!(b_hamsM, bM.ham)
                end
            end
        end

        # ---- KMC ----
        Random.seed!(43)
        aK, bK, slistK = setup()
        sysK = KMCSystem(; slist=slistK, beta=β)
        buildMoves!(sysK)

        for _ in 1:n_burn
            kmc_step!(sysK)
        end

        a_hamsK = Float64[]; b_hamsK = Float64[]
        for k in 1:n_coll
            kmc_step!(sysK)
            if k % every == 0
                ham!(aK); ham!(bK)
                push!(a_hamsK, aK.ham); push!(b_hamsK, bK.ham)
            end
        end

        μaM, σaM = mean(a_hamsM), std(a_hamsM)
        μbM, σbM = mean(b_hamsM), std(b_hamsM)
        μaK, σaK = mean(a_hamsK), std(a_hamsK)
        μbK, σbK = mean(b_hamsK), std(b_hamsK)

        # generous absolute tolerance to absorb autocorrelation + finite N
        @test abs(μaM - μaK) < 0.20
        @test abs(μbM - μbK) < 0.20
        @test abs(σaM - σaK) < 0.25
        @test abs(σbM - σbK) < 0.25

        # both should be physically negative (binding interactions present)
        @test μaM < 0.0 && μaK < 0.0
        @test μbM < 0.0 && μbK < 0.0
    end
end
