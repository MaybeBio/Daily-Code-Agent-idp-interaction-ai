using Test
using Droplet
using Random

# v2.3 gating tests: incremental polymer refresh equals from-scratch rebuild,
# and the fine-grained polymer enumeration agrees with the simple "refresh
# every polymer slot" fallback. The gate is correctness, not performance.
#
# Strategy: for each event, run kmc_step! (which calls the incremental
# refresh in-place), snapshot sys.rates / sys.Q, then call rebuildQ! and
# compare. After the comparison sys is in the rebuilt state, so subsequent
# steps still exercise a single-step delta from the (correct) baseline.

# Recompute every slot in `sys` from scratch and compare to the snapshot.
function _compare_to_rebuild(sys::KMCSystem, r_inc::Vector{Float64}, Q_inc::Float64;
                              atol::Float64=1e-12)
    rebuildQ!(sys)
    @test isapprox(Q_inc, sys.Q; atol=max(atol, abs(sys.Q) * 1e-10))
    mismatch = 0
    for slot in 1:sys.nslots
        if !isapprox(r_inc[slot], sys.rates[slot]; atol=atol)
            mismatch += 1
        end
    end
    @test mismatch == 0
end


# Run a polymer-bearing system for `nevents` steps under `mode`, asserting
# rates match a fresh rebuild after each step.
function _run_and_check(slist, mode::Symbol, nevents::Int, seed::Int;
                        beta::Float64=1.0)
    Random.seed!(seed)
    sys = KMCSystem(; slist=slist, beta=beta, refresh_mode=mode)
    buildMoves!(sys)
    for ev in 1:nevents
        sys.Q <= 0 && break
        kmc_step!(sys)
        r_inc = copy(sys.rates)
        Q_inc = sys.Q
        _compare_to_rebuild(sys, r_inc, Q_inc)
    end
    return sys
end


@testset "KMC Polymer v2.3: incremental refresh generalization" begin

    @testset "config (a) polymer-only — fine == rebuild" begin
        function _make()
            l = Lattice(; boxlength=Int32(20), nspecies=Int32(1))
            p = Polymer(; length=Int32(8), maxlength=Int32(20),
                          nMonomerTypes=Int32(1),
                          initShape="line", initPos=Int32[1,10,10],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=0.4)
            intPot!(p.monomers[1], -0.5)
            slist = [p]
            init!(l, slist)
            return slist
        end
        _run_and_check(_make(), :fine_polymer, 1000, 42)
    end

    @testset "config (b) polymer + bulk (self-only) — fine == rebuild" begin
        function _make()
            l = Lattice(; boxlength=Int32(16), nspecies=Int32(2))
            p = Polymer(; length=Int32(6), maxlength=Int32(16),
                          nMonomerTypes=Int32(1),
                          initShape="line", initPos=Int32[1,8,8],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=0.3)
            b = BulkMonomer(; length=Int32(15), initShape="uniform", k0=1.0)
            intPot!(p.monomers[1], -0.4)
            intPot!(b, -0.2)
            # No cross interaction (intMat self-only); polymer and bulk move
            # independently in the same lattice but don't see each other's
            # neighbor counts. Still exercises bulk slots being affected by
            # nearby bulk hops.
            slist = [p, b]
            init!(l, slist)
            return slist
        end
        _run_and_check(_make(), :fine_polymer, 1000, 43)
    end

    @testset "config (c) polymer + bulk with cross interaction — fine == rebuild" begin
        function _make()
            l = Lattice(; boxlength=Int32(16), nspecies=Int32(2))
            p = Polymer(; length=Int32(6), maxlength=Int32(16),
                          nMonomerTypes=Int32(1),
                          initShape="line", initPos=Int32[1,8,8],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=0.3)
            b = BulkMonomer(; length=Int32(20), initShape="uniform", k0=1.0)
            intPot!(p.monomers[1], -0.4)
            intPot!(b, -0.2)
            intPot!(p.monomers[1], b, -0.3)
            slist = [p, b]
            init!(l, slist)
            return slist
        end
        _run_and_check(_make(), :fine_polymer, 1000, 44)
    end

    @testset "config (a) polymer-only — simple == rebuild" begin
        function _make()
            l = Lattice(; boxlength=Int32(20), nspecies=Int32(1))
            p = Polymer(; length=Int32(8), maxlength=Int32(20),
                          nMonomerTypes=Int32(1),
                          initShape="line", initPos=Int32[1,10,10],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=0.4)
            intPot!(p.monomers[1], -0.5)
            slist = [p]
            init!(l, slist)
            return slist
        end
        _run_and_check(_make(), :simple_polymer, 500, 42)
    end

    @testset "config (c) polymer + bulk cross — simple == rebuild" begin
        function _make()
            l = Lattice(; boxlength=Int32(16), nspecies=Int32(2))
            p = Polymer(; length=Int32(6), maxlength=Int32(16),
                          nMonomerTypes=Int32(1),
                          initShape="line", initPos=Int32[1,8,8],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=0.3)
            b = BulkMonomer(; length=Int32(20), initShape="uniform", k0=1.0)
            intPot!(p.monomers[1], -0.4)
            intPot!(b, -0.2)
            intPot!(p.monomers[1], b, -0.3)
            slist = [p, b]
            init!(l, slist)
            return slist
        end
        _run_and_check(_make(), :simple_polymer, 500, 44)
    end

    @testset "fine vs simple equivalence: identical rates per event" begin
        # Same initial state, same RNG seed, two refresh modes. Per-step
        # comparison: incremental rates from each mode must match each other.
        # Catches under-coverage in the fine-grained polymer enumeration —
        # any slot the fine mode misses but simple covers shows up here.
        # Both slists must start in the IDENTICAL state — init! places bulk
        # particles randomly using the global RNG, so reseed before each
        # construction so the two systems see the same initial lattice.
        function _make_one(init_seed::Int)
            l = Lattice(; boxlength=Int32(16), nspecies=Int32(2))
            p = Polymer(; length=Int32(6), maxlength=Int32(16),
                          nMonomerTypes=Int32(1),
                          initShape="line", initPos=Int32[1,8,8],
                          pattern=:uniform,
                          k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
            initMonomers!(p, defaultChemPot=0.3)
            b = BulkMonomer(; length=Int32(20), initShape="uniform", k0=1.0)
            intPot!(p.monomers[1], -0.4)
            intPot!(b, -0.2)
            intPot!(p.monomers[1], b, -0.3)
            slist = [p, b]
            Random.seed!(init_seed)
            init!(l, slist)
            return slist
        end

        slist_fine   = _make_one(2026)
        slist_simple = _make_one(2026)
        sys_fine   = KMCSystem(; slist=slist_fine,   beta=1.0,
                               refresh_mode=:fine_polymer)
        sys_simple = KMCSystem(; slist=slist_simple, beta=1.0,
                               refresh_mode=:simple_polymer)
        buildMoves!(sys_fine)
        buildMoves!(sys_simple)

        # Initial Q must agree (same state, same rates).
        @test isapprox(sys_fine.Q, sys_simple.Q; atol=1e-12)

        # Drive both systems with the same RNG sequence — kmc_step! consumes
        # exactly two rand() calls (sampleMove + log time advance), so we can
        # reseed before each pair of steps.
        for ev in 1:500
            sys_fine.Q <= 0 && break
            sys_simple.Q <= 0 && break

            Random.seed!(1000 + ev)
            kmc_step!(sys_fine)
            Random.seed!(1000 + ev)
            kmc_step!(sys_simple)

            @test isapprox(sys_fine.Q, sys_simple.Q; atol=1e-9)
            for slot in 1:sys_fine.nslots
                @test isapprox(sys_fine.rates[slot], sys_simple.rates[slot];
                               atol=1e-12)
            end
        end
    end

    @testset "rebuild_all mode preserves v2.2 behavior (sanity)" begin
        # Quick smoke test that the :rebuild_all mode still drives a polymer
        # simulation correctly — its rates trivially match a fresh rebuild
        # since that's exactly what it does.
        l = Lattice(; boxlength=Int32(20), nspecies=Int32(1))
        p = Polymer(; length=Int32(8), maxlength=Int32(20),
                      nMonomerTypes=Int32(1),
                      initShape="line", initPos=Int32[1,10,10],
                      pattern=:uniform,
                      k0_add=1.0, k0_kink=1.0, k0_remove=1.0)
        initMonomers!(p, defaultChemPot=0.4)
        intPot!(p.monomers[1], -0.5)
        slist = [p]
        init!(l, slist)

        Random.seed!(99)
        sys = KMCSystem(; slist=slist, beta=1.0, refresh_mode=:rebuild_all)
        buildMoves!(sys)
        for _ in 1:200
            sys.Q <= 0 && break
            kmc_step!(sys)
        end
        Q_after = sys.Q
        rebuildQ!(sys)
        @test isapprox(Q_after, sys.Q; atol=1e-12)
    end

    @testset "BulkMonomer-only path uses incremental refresh in :fine_polymer" begin
        # Regression check that the dispatch refactor didn't break v1: when
        # no polymer is in slist, :fine_polymer should fall back to the bulk
        # incremental refresh and stay consistent with rebuilds.
        Random.seed!(46)
        l = Lattice(; boxlength=Int32(8), nspecies=Int32(1))
        b = BulkMonomer(; length=Int32(10), initShape="uniform", k0=1.0)
        intPot!(b, -0.3)
        slist = [b]
        init!(l, slist)

        sys = KMCSystem(; slist=slist, beta=1.0, refresh_mode=:fine_polymer)
        buildMoves!(sys)
        for _ in 1:200
            sys.Q <= 0 && break
            kmc_step!(sys)
        end
        Q_inc = sys.Q
        rebuildQ!(sys)
        @test isapprox(Q_inc, sys.Q; rtol=1e-10)
    end
end
