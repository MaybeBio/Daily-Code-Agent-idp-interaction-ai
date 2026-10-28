# Test file for integration tests: FastBulkMonomer exclusion and
# full simulation with on-site coupling (two-way + three-way)

using Test
using Droplet

println("Testing FastBulkMonomer exclusion & integration (coupling)")
println("=" ^ 65)

# =============================================================================
# PART 1: FastBulkMonomer exclusion
# =============================================================================

@testset "FastBulkMonomer - couplingPot! explicit error (3-arg)" begin
    f = FastBulkMonomer(; name="F", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)

    intPot!(f, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)

    @test_throws ErrorException couplingPot!(f, b, c, -0.5)

    try
        couplingPot!(f, b, c, -0.5)
    catch e
        @test occursin("FastBulkMonomer", e.msg)
    end
end

@testset "FastBulkMonomer - couplingPot! explicit error (2-arg)" begin
    f = FastBulkMonomer(; name="F", length=5)
    b = BulkMonomer(; name="B", length=5)

    intPot!(f, -1.0)
    intPot!(b, -1.0)

    @test_throws ErrorException couplingPot!(f, b, -0.5)

    try
        couplingPot!(f, b, -0.5)
    catch e
        @test occursin("FastBulkMonomer", e.msg)
    end
end

@testset "FastBulkMonomer - ham! regression with BulkMonomer coupling present" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    f = FastBulkMonomer(; name="F", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(f, -0.5)
    intPot!(c, -0.5)
    intPot!(a, f, -0.3)
    intPot!(a, c, -0.2)
    intPot!(f, c, -0.1)

    couplingPot!(a, f, c, -1.0)

    init!(l, [a, f, c])

    ham!(f)
    @test isfinite(f.ham)

    ham!(a)
    @test isfinite(a.ham)
end

@testset "FastBulkMonomer - deltaHam!/move! regression with coupling present" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=15, chemPot=0.5)
    f = FastBulkMonomer(; name="F", length=15, chemPot=0.5)
    c = BulkMonomer(; name="C", length=15, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(f, -0.5)
    intPot!(c, -0.5)
    intPot!(a, f, -0.3)
    intPot!(a, c, -0.2)
    intPot!(f, c, -0.1)

    couplingPot!(a, f, c, -0.5)

    init!(l, [a, f, c])

    ham!(a); ham!(f); ham!(c)

    n_moves = 0
    for _ in 1:200
        if move!(f)
            n_moves += 1
        end
    end
    @test n_moves > 0

    ham!(f)
    hams_f = [f.ham]
    acc_f = [f.ham]
    acc = 0.0

    for _ in 1:500
        if move!(f)
            acc += getDeltaHam(f)
            push!(acc_f, acc_f[1] + acc)
            ham!(f)
            push!(hams_f, f.ham)
        end
    end

    @test length(hams_f) > 1
    for i in eachindex(hams_f)
        @test isapprox(hams_f[i], acc_f[i], atol=1e-10)
    end
end

@testset "FastBulkMonomer - has no couplingMat field" begin
    f = FastBulkMonomer(; name="F", length=5)
    @test !hasfield(typeof(f), :couplingMat)
    @test !hasfield(typeof(f), :couplingIdx)
    @test !hasfield(typeof(f), :couplingSpecies)
end

# =============================================================================
# PART 2: Three BulkMonomer — three-way coupling
# =============================================================================

@testset "Integration - three BulkMonomer species with three-way coupling" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0); intPot!(b, 0.0); intPot!(c, 0.0)
    intPot!(a, b, 0.0); intPot!(a, c, 0.0); intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    # All three at same site
    a.latS .= 0; b.latS .= 0; c.latS .= 0

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 10, 5, 5
    a.latS[10, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)
    @test isapprox(a.ham, -1.0, atol=1e-10)

    # Move A away: lose coupling
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5

    deltaHam!(a)
    @test isapprox(a.deltaH, 1.0, atol=1e-10)

    makeMove!(a)
    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

# =============================================================================
# PART 3: Three BulkMonomer — two-way coupling
# =============================================================================

@testset "Integration - three BulkMonomer with two-way coupling" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)

    intPot!(a, 0.0); intPot!(b, 0.0); intPot!(a, b, 0.0)

    couplingPot!(a, b, -1.0)

    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0

    # Co-located
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 10, 5, 5
    a.latS[10, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1

    ham!(a)
    @test isapprox(a.ham, -1.0, atol=1e-10)

    # Move A away
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5

    deltaHam!(a)
    @test isapprox(a.deltaH, 1.0, atol=1e-10)

    makeMove!(a)
    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

# =============================================================================
# PART 4: BulkMonomer + PolyMonomer coupling
# =============================================================================

@testset "Integration - BulkMonomer with PolyMonomer coupling" begin
    l = Lattice(; boxlength=20, nspecies=2)

    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]

    a = BulkMonomer(; name="A", length=5, chemPot=0.0)
    b = BulkMonomer(; name="B", length=5, chemPot=0.0)

    intPot!(a, 0.0); intPot!(b, 0.0)
    intPot!(a, pm, 0.0); intPot!(b, pm, 0.0)
    intPot!(a, b, 0.0)

    couplingPot!(a, b, pm, -0.5)

    init!(l, [poly, a, b])

    n_consistent = 0
    n_trials = 100

    for _ in 1:n_trials
        ham!(a)
        ham_before = a.ham

        a.moveIdx = rand(Int32(1):Int32(a.length))
        identifyBond!(a)
        proposeMove!(a)

        if a.latS[a.newPos[1], a.newPos[2], a.newPos[3]] > 0
            n_consistent += 1
            continue
        end

        deltaHam!(a)
        delta = a.deltaH

        makeMove!(a)
        ham!(a)
        ham_after = a.ham

        if isapprox(ham_before + delta, ham_after, atol=1e-10)
            n_consistent += 1
        end
    end

    @test n_consistent == n_trials
end

# =============================================================================
# PART 5: Heterogeneous polymer — coupling discriminates monomer types
# =============================================================================

@testset "Integration - heterogeneous polymer coupling discrimination" begin
    l = Lattice(; boxlength=20, nspecies=3)

    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(a, pm1, 0.0); intPot!(a, pm2, 0.0)

    couplingPot!(a, pm1, -1.0)  # two-way: A and pm1 on-site only

    init!(l, [poly, a])

    # Find pm1 and pm2 sites
    pm1_site = nothing; pm2_site = nothing
    for i in 1:poly.length
        p = poly.pos[i]
        if p[6] == pm1.intIdx && pm1_site === nothing
            pm1_site = (p[1], p[2], p[3])
        elseif p[6] == pm2.intIdx && pm2_site === nothing
            pm2_site = (p[1], p[2], p[3])
        end
        pm1_site !== nothing && pm2_site !== nothing && break
    end

    @test pm1_site !== nothing
    @test pm2_site !== nothing

    # Place A at pm1 site → coupling active
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_site[1], pm1_site[2], pm1_site[3]
    a.latS[pm1_site[1], pm1_site[2], pm1_site[3]] = 1

    ham!(a)
    ham_at_pm1 = a.ham

    # Place A at pm2 site → coupling should NOT fire
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm2_site[1], pm2_site[2], pm2_site[3]
    a.latS[pm2_site[1], pm2_site[2], pm2_site[3]] = 1

    ham!(a)
    ham_at_pm2 = a.ham

    @test isapprox(ham_at_pm1, -1.0, atol=1e-10)
    @test isapprox(ham_at_pm2, 0.0, atol=1e-10)
end

# =============================================================================
# PART 6: MC consistency — self-moves
# =============================================================================

@testset "Integration - MC consistency (self-moves, accumulated deltaH matches ham!)" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=30, chemPot=0.8)
    b = BulkMonomer(; name="B", length=30, chemPot=0.8)
    c = BulkMonomer(; name="C", length=30, chemPot=0.8)

    intPot!(a, -0.5); intPot!(b, -0.5); intPot!(c, -0.5)
    intPot!(a, b, -0.3); intPot!(a, c, -0.2); intPot!(b, c, -0.1)

    couplingPot!(a, b, -0.4)
    couplingPot!(a, b, c, -0.7)

    init!(l, [a, b, c])

    ham!(a); ham!(b); ham!(c)

    for owner in [a, b, c]
        ham!(a); ham!(b); ham!(c)

        hams = [getHam(owner)]
        accDelta = [getHam(owner)]
        acc = 0.0

        for _ in 1:500
            if move!(owner)
                acc += getDeltaHam(owner)
                push!(accDelta, accDelta[1] + acc)
                ham!(owner)
                push!(hams, getHam(owner))
            end
        end

        @test length(hams) > 1
        for i in eachindex(hams)
            @test isapprox(hams[i], accDelta[i], atol=1e-10)
        end
    end
end

# =============================================================================
# PART 7: MC consistency — all species move, single-step
# =============================================================================

@testset "Integration - MC consistency (all species move, single-step)" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=30, chemPot=0.8)
    b = BulkMonomer(; name="B", length=30, chemPot=0.8)
    c = BulkMonomer(; name="C", length=30, chemPot=0.8)

    intPot!(a, -0.5); intPot!(b, -0.5); intPot!(c, -0.5)
    intPot!(a, b, -0.3); intPot!(a, c, -0.2); intPot!(b, c, -0.1)

    couplingPot!(a, b, -0.4)
    couplingPot!(a, b, c, -0.7)

    init!(l, [a, b, c])

    ham!(a); ham!(b); ham!(c)

    n_consistent = 0
    n_total = 0

    for _ in 1:500
        for (s, idx) in a.intSpecies
            ham!(s)
            ham_before = s.ham

            s.moveIdx = rand(Int32(1):Int32(s.length))
            identifyBond!(s)
            proposeMove!(s)

            if s.latS[s.newPos[1], s.newPos[2], s.newPos[3]] > 0
                n_total += 1
                n_consistent += 1
                continue
            end

            deltaHam!(s)
            delta = s.deltaH

            makeMove!(s)
            ham!(s)
            ham_after = s.ham

            n_total += 1
            if isapprox(ham_before + delta, ham_after, atol=1e-10)
                n_consistent += 1
            end

            ham!(a); ham!(b); ham!(c)
        end
    end

    @test n_total > 0
    @test n_consistent == n_total
end

# =============================================================================
# PART 8: Detailed balance
# =============================================================================

@testset "Integration - detailed balance (forward/reverse energy symmetry)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0); intPot!(b, 0.0); intPot!(c, 0.0)
    intPot!(a, b, 0.0); intPot!(a, c, 0.0); intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    # A at (5,5,5), B and C at (10,5,5)
    a.latS .= 0; b.latS .= 0; c.latS .= 0

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)
    ham_initial = a.ham  # 0: B,C not at A's site

    # Forward: A → (10,5,5) (joins B,C)
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 10, 5, 5
    deltaHam!(a)
    delta_forward = a.deltaH
    makeMove!(a)
    ham!(a)
    ham_final = a.ham

    @test isapprox(ham_initial + delta_forward, ham_final, atol=1e-10)

    # Reverse: A → (5,5,5) (leaves B,C)
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5
    deltaHam!(a)
    delta_reverse = a.deltaH

    @test isapprox(delta_forward + delta_reverse, 0.0, atol=1e-10)

    makeMove!(a)
    ham!(a)
    @test isapprox(a.ham, ham_initial, atol=1e-10)
end

# =============================================================================
# PART 9: Multiple coupling owners — all BulkMonomer see same value
# =============================================================================

@testset "Integration - multiple BulkMonomers each owning couplings" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)

    intPot!(a, -0.3); intPot!(b, -0.3); intPot!(c, -0.3)
    intPot!(a, b, -0.2); intPot!(a, c, -0.1); intPot!(b, c, -0.1)

    # With first-two-index symmetry, all three BulkMonomers share the coupling
    couplingPot!(a, b, c, -0.5)

    init!(l, [a, b, c])

    # All three should have coupling matrices (due to symmetry)
    @test a.couplingMat !== nothing
    @test b.couplingMat !== nothing
    @test c.couplingMat !== nothing

    # Consistency test for all three owners
    for owner in [a, b, c]
        n_consistent = 0
        n_trials = 100

        for _ in 1:n_trials
            ham!(owner)
            ham_before = owner.ham

            owner.moveIdx = rand(Int32(1):Int32(owner.length))
            identifyBond!(owner)
            proposeMove!(owner)

            if owner.latS[owner.newPos[1], owner.newPos[2], owner.newPos[3]] > 0
                n_consistent += 1
                continue
            end

            deltaHam!(owner)
            delta = owner.deltaH

            makeMove!(owner)
            ham!(owner)
            ham_after = owner.ham

            if isapprox(ham_before + delta, ham_after, atol=1e-10)
                n_consistent += 1
            end
        end

        @test n_consistent == n_trials
    end
end

# =============================================================================
# PART 10: Mixed FastBulkMonomer + BulkMonomer simulation
# =============================================================================

@testset "Integration - mixed FastBulkMonomer and BulkMonomer simulation" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.8)
    f = FastBulkMonomer(; name="F", length=20, chemPot=0.8)
    c = BulkMonomer(; name="C", length=20, chemPot=0.8)

    intPot!(a, -0.5); intPot!(f, -0.5); intPot!(c, -0.5)
    intPot!(a, f, -0.3); intPot!(a, c, -0.2); intPot!(f, c, -0.1)

    couplingPot!(a, f, c, -0.4)

    init!(l, [a, f, c])

    # A's self-move consistency
    ham!(a); ham!(f); ham!(c)

    a_hams = [getHam(a)]
    a_accDelta = [getHam(a)]
    a_acc = 0.0

    for _ in 1:500
        if move!(a)
            a_acc += getDeltaHam(a)
            push!(a_accDelta, a_accDelta[1] + a_acc)
            ham!(a)
            push!(a_hams, getHam(a))
        end
    end

    @test length(a_hams) > 1
    for i in eachindex(a_hams)
        @test isapprox(a_hams[i], a_accDelta[i], atol=1e-10)
    end

    # F's self-move consistency (no coupling, pure two-way intMat)
    ham!(f)
    f_hams = [getHam(f)]
    f_accDelta = [getHam(f)]
    f_acc = 0.0

    for _ in 1:500
        if move!(f)
            f_acc += getDeltaHam(f)
            push!(f_accDelta, f_accDelta[1] + f_acc)
            ham!(f)
            push!(f_hams, getHam(f))
        end
    end

    @test length(f_hams) > 1
    for i in eachindex(f_hams)
        @test isapprox(f_hams[i], f_accDelta[i], atol=1e-10)
    end
end

# =============================================================================
# PART 11: BulkMonomer + Polymer full simulation
# =============================================================================

@testset "Integration - BulkMonomer coupling with Polymer in simulation" begin
    l = Lattice(; boxlength=15, nspecies=3)

    poly = Polymer(; length=8, maxlength=16, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.2)
    pm = poly.monomers[1]

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)

    intPot!(a, -0.3); intPot!(b, -0.3)
    intPot!(a, pm, -0.1); intPot!(b, pm, -0.1)
    intPot!(a, b, -0.2)

    couplingPot!(a, b, pm, -0.5)

    init!(l, [poly, a, b])

    ham!(a); ham!(b); ham!(pm)

    a_hams = [getHam(a)]
    a_accDelta = [getHam(a)]
    a_acc = 0.0

    for _ in 1:500
        if move!(a)
            a_acc += getDeltaHam(a)
            push!(a_accDelta, a_accDelta[1] + a_acc)
            ham!(a)
            push!(a_hams, getHam(a))
        end
    end

    @test length(a_hams) > 1
    for i in eachindex(a_hams)
        @test isapprox(a_hams[i], a_accDelta[i], atol=1e-10)
    end
end

# =============================================================================
# PART 12: Mixed two-way + three-way simulation
# =============================================================================

@testset "Integration - mixed two-way + three-way simulation" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=25, chemPot=0.5)
    b = BulkMonomer(; name="B", length=25, chemPot=0.5)
    c = BulkMonomer(; name="C", length=25, chemPot=0.5)

    intPot!(a, -0.3); intPot!(b, -0.3); intPot!(c, -0.3)
    intPot!(a, b, -0.2); intPot!(a, c, -0.1); intPot!(b, c, -0.1)

    couplingPot!(a, b, -0.5)     # two-way
    couplingPot!(a, b, c, -0.3)  # three-way

    init!(l, [a, b, c])

    ham!(a); ham!(b); ham!(c)

    for owner in [a, b, c]
        ham!(a); ham!(b); ham!(c)

        hams = [getHam(owner)]
        accDelta = [getHam(owner)]
        acc = 0.0

        for _ in 1:500
            if move!(owner)
                acc += getDeltaHam(owner)
                push!(accDelta, accDelta[1] + acc)
                ham!(owner)
                push!(hams, getHam(owner))
            end
        end

        @test length(hams) > 1
        for i in eachindex(hams)
            @test isapprox(hams[i], accDelta[i], atol=1e-10)
        end
    end
end

# =============================================================================
# PART 13: Performance
# =============================================================================

@testset "Integration - performance acceptable" begin
    l = Lattice(; boxlength=30, nspecies=3)

    a = BulkMonomer(; name="A", length=200, chemPot=1.0)
    b = BulkMonomer(; name="B", length=200, chemPot=1.0)
    c = BulkMonomer(; name="C", length=200, chemPot=1.0)

    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    intPot!(a, b, -0.5); intPot!(a, c, -0.5); intPot!(b, c, -0.5)

    couplingPot!(a, b, c, -0.3)
    couplingPot!(a, b, -0.2)

    init!(l, [a, b, c])

    ham!(a)
    for _ in 1:10
        move!(a)
    end

    t = @elapsed for _ in 1:1000
        move!(a)
    end

    @test t < 10.0
    @test isfinite(a.ham)

    println("  Performance: 1000 move! calls in $(round(t, digits=3))s")
end

println()
println("All integration coupling tests completed!")
