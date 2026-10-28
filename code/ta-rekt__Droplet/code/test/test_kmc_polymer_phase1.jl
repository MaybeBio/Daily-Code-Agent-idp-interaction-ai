using Test
using Droplet

# For each polymer slot, run the existing Metropolis ΔH path (identifyBond! +
# manual moveType + deltaHam!) and compare to the pure deltaH_kmc. This is the
# v2.1 gating equivalence test.
function _metropolis_delta_for_slot(p, sys, ki::Int32, dir::Int)
    # Save the energy fields we're about to mutate so we can restore them.
    saved_ham    = p.ham
    saved_deltaH = p.deltaH
    saved_dHP    = p.deltaHP

    # Snapshot partner deltaHPs as well — deltaHam!(p) calls setDeltaHP on
    # external bulks via the cross-interaction branch.
    saved_partners = Tuple{Any, Float64}[]
    m = first(values(p.monomers))
    for (sp, _) in m.intSpecies
        if !(sp isa PolyMonomer && sp.poly === p)
            push!(saved_partners, (sp, sp.deltaHP === nothing ? NaN : sp.deltaHP))
        end
    end

    identifyBond!(p; ki=ki)
    moveType = Droplet._polyMoveVec!(p, dir)
    p.moveType = moveType

    deltaHam!(p)
    dh = p.deltaH

    # restore mutated state so the test can keep looping
    p.ham    = saved_ham
    p.deltaH = saved_deltaH
    p.deltaHP = saved_dHP
    for (sp, val) in saved_partners
        if !isnan(val)
            sp.deltaHP = val
        end
    end

    return dh, moveType
end

@testset "KMC Polymer v2.1: slot encoding + rate primitives" begin

    @testset "nslots, slot count, active-status" begin
        l = Lattice(; boxlength=Int32(20), nspecies=Int32(1))
        p = Polymer(; length=Int32(8), maxlength=Int32(20),
                      nMonomerTypes=Int32(1),
                      initShape="line", initPos=Int32[1,10,10],
                      pattern=:uniform)
        initMonomers!(p)
        intPot!(p.monomers[1], -0.3)
        slist = [p]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        # nslots = 4 * maxlength = 4 * 20 = 80 (one slot batch per stable
        # bondIdx, including the chain tail's always-inactive batch and any
        # unallocated-bondIdx batches).
        @test nslots(p) == 80
        @test sys.nslots == 80
        @test length(sys.rates) == 80

        # Active bonds: chain has 7 bonds (length 8 → 7 forward links); each
        # contributes 4 directional slots → 28 active candidates.
        nactive = count(r -> r > 0.0, sys.rates)
        @test nactive == 28

        # Inactive: bondIdx 7 is the chain tail (no successor) and bondIdxs
        # 8..19 are unallocated. All have rate 0.
        @test sys.rates[8 * 4 - 3] == 0.0   # last bead's slots (bondIdx 7, dir 1)

        # bondIdx state book-keeping
        state = sys.polymer_states[p]
        @test length(state.bondIdx_of_key) == 8
        @test length(state.free_bondIdxs) == 12
        @test sort(collect(values(state.bondIdx_of_key))) == collect(0:7)
        for k in keys(state.bondIdx_of_key)
            @test state.key_at_bondIdx[state.bondIdx_of_key[k] + 1] == k
        end
    end

    @testset "deltaH_kmc parity with mutating deltaHam! across all slots" begin
        # Polymer-only system: heterogeneous-pattern branch and external cross
        # branches both excluded; only the self-interaction path exercised.
        l = Lattice(; boxlength=Int32(16), nspecies=Int32(1))
        p = Polymer(; length=Int32(6), maxlength=Int32(16),
                      nMonomerTypes=Int32(1),
                      initShape="line", initPos=Int32[1,8,8],
                      pattern=:uniform)
        initMonomers!(p, defaultChemPot=0.4)
        intPot!(p.monomers[1], -0.5)
        slist = [p]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        state = sys.polymer_states[p]
        for bondIdx in 0:(p.length-2)   # exclude chain tail (kj == 0)
            ki = state.key_at_bondIdx[bondIdx + 1]
            ki == 0 && continue
            for dir in 1:4
                ref_dh, ref_type = _metropolis_delta_for_slot(p, sys, ki, dir)
                identifyBond!(p; ki=ki)
                kmc_type = Droplet._polyMoveVec!(p, dir)
                @test kmc_type == ref_type
                kmc_dh = deltaH_kmc(p, ki, kmc_type)
                @test isapprox(kmc_dh, ref_dh; atol=1e-12, rtol=1e-12)
            end
        end
    end

    @testset "deltaH_kmc parity with bulk cross-interactions" begin
        # Polymer + 1 BulkMonomer with cross-interaction. Exercises the
        # m.intIdx != idx branch.
        l = Lattice(; boxlength=Int32(12), nspecies=Int32(2))
        p = Polymer(; length=Int32(5), maxlength=Int32(12),
                      nMonomerTypes=Int32(1),
                      initShape="line", initPos=Int32[1,6,6],
                      pattern=:uniform)
        initMonomers!(p, defaultChemPot=0.2)
        b = BulkMonomer(; length=Int32(8), initShape="uniform", k0=1.0)
        intPot!(p.monomers[1], -0.3)
        intPot!(b, -0.2)
        intPot!(p.monomers[1], b, -0.4)
        slist = [p, b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)

        state = sys.polymer_states[p]
        for bondIdx in 0:(p.length-2)   # exclude chain tail (kj == 0)
            ki = state.key_at_bondIdx[bondIdx + 1]
            ki == 0 && continue
            for dir in 1:4
                ref_dh, ref_type = _metropolis_delta_for_slot(p, sys, ki, dir)
                identifyBond!(p; ki=ki)
                kmc_type = Droplet._polyMoveVec!(p, dir)
                @test kmc_type == ref_type
                kmc_dh = deltaH_kmc(p, ki, kmc_type)
                @test isapprox(kmc_dh, ref_dh; atol=1e-12, rtol=1e-12)
            end
        end
    end

    @testset "rates are k0_<type> · min(1, exp(-βΔH))" begin
        l = Lattice(; boxlength=Int32(16), nspecies=Int32(1))
        p = Polymer(; length=Int32(5), maxlength=Int32(20),
                      nMonomerTypes=Int32(1),
                      initShape="line", initPos=Int32[1,8,8],
                      pattern=:uniform,
                      k0_add=2.0, k0_kink=3.0, k0_remove=4.0)
        initMonomers!(p, defaultChemPot=0.5)
        intPot!(p.monomers[1], -0.4)
        slist = [p]
        init!(l, slist)

        β = 0.8
        sys = KMCSystem(; slist=slist, beta=β)
        buildMoves!(sys)

        state = sys.polymer_states[p]
        # For a straight-line initial polymer with the default geometry, all
        # active bonds have only `add` candidates available (no kinks since
        # bprev/bnext aren't perpendicular to the bond axis, no removes since
        # the chain isn't already U-folded).
        for bondIdx in 0:(p.length-2)
            ki = state.key_at_bondIdx[bondIdx + 1]
            for dir in 1:4
                slot = p === slist[1] ? bondIdx * 4 + dir : error()
                identifyBond!(p; ki=ki)
                mt = Droplet._polyMoveVec!(p, dir)
                @test mt == 1   # all are adds in a straight line
                expected = p.k0_add * min(1.0, exp(-β * deltaH_kmc(p, ki, mt)))
                @test isapprox(sys.rates[slot], expected; atol=1e-12)
            end
        end

        @test isapprox(sys.Q, sum(sys.rates); atol=1e-12)
        @test isapprox(fenwick_total(sys.fenwick), sys.Q; atol=1e-12)
    end

    @testset "BulkMonomer-only path still works (regression)" begin
        # Make sure the dispatch refactor didn't break v1.
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(10), initShape="uniform", k0=0.7)
        intPot!(b, -0.3)
        slist = [b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0)
        buildMoves!(sys)
        @test sys.nslots == 60
        @test isapprox(sys.Q, sum(sys.rates); atol=1e-12)
    end
end
