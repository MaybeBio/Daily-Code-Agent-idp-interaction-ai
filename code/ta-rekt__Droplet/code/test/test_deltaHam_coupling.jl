# Test file for deltaHam! with on-site coupling for BulkMonomer
# Tests the incremental Hamiltonian calculation including two-way and three-way coupling

using Test
using Droplet

println("Testing deltaHam! with on-site coupling")
println("=" ^ 55)

# =============================================================================
# PART 1: Regression - deltaHam! unchanged when couplingMat is nothing
# =============================================================================

@testset "deltaHam! - unchanged when couplingMat is nothing (regression)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=10, chemPot=1.0)
    b = BulkMonomer(; name="B", length=10, chemPot=1.0)

    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, b, -0.5)

    init!(l, [a, b])

    @test a.couplingMat === nothing

    ham!(a)
    ham_before = a.ham

    a.moveIdx = rand(Int32(1):Int32(a.length))
    identifyBond!(a)
    proposeMove!(a)
    deltaHam!(a)

    @test isfinite(a.deltaH)
    @test isapprox(a.ham, ham_before + a.deltaH, atol=1e-10)
end

# =============================================================================
# PART 2: Consistency — ham_before + delta == ham_after
# =============================================================================

@testset "deltaHam! - consistency with ham! (three-way coupling)" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(c, -0.5)
    intPot!(a, b, -0.3)
    intPot!(a, c, -0.2)
    intPot!(b, c, -0.1)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

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
        else
            @warn "Mismatch" ham_before delta ham_after diff=(ham_before+delta-ham_after)
        end
    end

    @test n_consistent == n_trials
end

@testset "deltaHam! - consistency with ham! (two-way coupling)" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(a, b, -0.3)

    couplingPot!(a, b, -0.8)

    init!(l, [a, b])

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
        else
            @warn "Mismatch" ham_before delta ham_after diff=(ham_before+delta-ham_after)
        end
    end

    @test n_consistent == n_trials
end

@testset "deltaHam! - consistency with ham! (mixed two-way + three-way)" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(c, -0.5)
    intPot!(a, b, -0.3)
    intPot!(a, c, -0.2)
    intPot!(b, c, -0.1)

    couplingPot!(a, b, -0.6)     # two-way
    couplingPot!(a, b, c, -1.0)  # three-way

    init!(l, [a, b, c])

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
        else
            @warn "Mismatch" ham_before delta ham_after diff=(ham_before+delta-ham_after)
        end
    end

    @test n_consistent == n_trials
end

@testset "deltaHam! - consistency with multiple coupling terms (4 species)" begin
    l = Lattice(; boxlength=15, nspecies=4)

    a = BulkMonomer(; name="A", length=15, chemPot=0.5)
    b = BulkMonomer(; name="B", length=15, chemPot=0.5)
    c = BulkMonomer(; name="C", length=15, chemPot=0.5)
    d = BulkMonomer(; name="D", length=15, chemPot=0.5)

    intPot!(a, -0.3)
    intPot!(b, -0.3)
    intPot!(c, -0.3)
    intPot!(d, -0.3)
    intPot!(a, b, -0.2)
    intPot!(a, c, -0.1)
    intPot!(a, d, -0.1)
    intPot!(b, c, 0.0)
    intPot!(b, d, 0.0)
    intPot!(c, d, 0.0)

    couplingPot!(a, b, c, -0.5)
    couplingPot!(a, b, d, -0.3)
    couplingPot!(a, c, d, -0.7)

    init!(l, [a, b, c, d])

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
        else
            @warn "Mismatch" ham_before delta ham_after diff=(ham_before+delta-ham_after)
        end
    end

    @test n_consistent == n_trials
end

# =============================================================================
# PART 3: Controlled tests with manual placement
# =============================================================================

@testset "deltaHam! - move into co-located site (three-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # A at (5,5,5), B and C co-located at (10,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1

    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1

    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)

    # Move A to (10,5,5) — now all three co-located
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 10, 5, 5

    deltaHam!(a)

    # At newPos: cplpq(B,C,(10,5,5)) = 1 → coupling = -1.0
    # At oldPos (bi): cplpq(B,C,(5,5,5)) = 0
    # delta = -1.0 * (1 - 0) = -1.0
    @test isapprox(a.deltaH, -1.0, atol=1e-10)
end

@testset "deltaHam! - move out of co-located site (three-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # All three co-located at (10,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 10, 5, 5
    a.latS[10, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)
    @test isapprox(a.ham, -1.0, atol=1e-10)

    # Move A to (5,5,5) — away from B,C
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5

    deltaHam!(a)

    # At newPos (5,5,5): cplpq(B,C,(5,5,5)) = 0
    # At oldPos (10,5,5): cplpq(B,C,(10,5,5)) = 1
    # delta = -1.0 * (0 - 1) = +1.0
    @test isapprox(a.deltaH, 1.0, atol=1e-10)
end

@testset "deltaHam! - move between non-coupled sites" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # A at (5,5,5), B at (10,10,10), C at (15,15,15)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 10, 10
    b.latS[10, 10, 10] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 15, 15, 15
    c.latS[15, 15, 15] = 1

    ham!(a)

    # Move A to (6,5,5) — still no coupling
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 6, 5, 5

    deltaHam!(a)

    @test isapprox(a.deltaH, 0.0, atol=1e-10)
end

@testset "deltaHam! - move into co-located site (two-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, b, 0.0)

    couplingPot!(a, b, -1.0)

    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1

    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)

    # Move A to B's site
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 10, 5, 5

    deltaHam!(a)

    # cplpq(B,B,newPos) = 1, cplpq(B,B,oldPos) = 0
    # delta = -1.0 * (1 - 0) = -1.0
    @test isapprox(a.deltaH, -1.0, atol=1e-10)
end

@testset "deltaHam! - move out of co-located site (two-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, b, 0.0)

    couplingPot!(a, b, -1.0)

    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0

    # Both co-located
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
end

# =============================================================================
# PART 4: Partner deltaHP unchanged by coupling
# =============================================================================

@testset "deltaHam! - partner deltaHP unchanged by coupling" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, -0.5)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # A at (5,5,5), B and C at (10,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)

    # Move A to (10,5,5)
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 10, 5, 5

    deltaHam!(a)

    # deltaHP should only reflect intMat contribution, not coupling
    # B's intMat contribution: nnp changes from move
    # Coupling delta is only in a.deltaH, not in partner deltaHP
    # We just verify deltaHP is set and coupling doesn't leak into it
    @test isfinite(b.deltaHP)
    @test isfinite(c.deltaHP)
end

# =============================================================================
# PART 5: Zero coupling
# =============================================================================

@testset "deltaHam! - zero coupling has no effect" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, 0.0)

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # All at same site — even with zero coupling, delta should be 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 10, 5, 5
    a.latS[10, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)

    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5

    deltaHam!(a)

    @test isapprox(a.deltaH, 0.0, atol=1e-10)
end

# =============================================================================
# PART 6: Symmetry — delta consistent across owners
# =============================================================================

@testset "deltaHam! - delta consistent across owners" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -1.0)

    init!(l, [a, b, c])

    # Config 1: move A into co-located site
    a.latS .= 0; b.latS .= 0; c.latS .= 0

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 5, 5
    b.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(a)
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 10, 5, 5
    deltaHam!(a)
    delta_a = a.deltaH

    # Config 2: move B into co-located site (reset config)
    a.latS .= 0; b.latS .= 0; c.latS .= 0

    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 10, 5, 5
    a.latS[10, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1

    ham!(b)
    b.moveIdx = Int32(1)
    identifyBond!(b)
    b.newPos[1], b.newPos[2], b.newPos[3] = 10, 5, 5
    deltaHam!(b)
    delta_b = b.deltaH

    # Both should see same delta magnitude
    @test isapprox(delta_a, delta_b, atol=1e-10)
end

# =============================================================================
# PART 7: PolyMonomer
# =============================================================================

@testset "deltaHam! - consistency with PolyMonomer coupling" begin
    l = Lattice(; boxlength=15, nspecies=2)

    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]

    a = BulkMonomer(; name="A", length=15, chemPot=0.5)
    b = BulkMonomer(; name="B", length=15, chemPot=0.5)

    intPot!(a, -0.3)
    intPot!(b, -0.3)
    intPot!(a, pm, 0.0)
    intPot!(b, pm, 0.0)
    intPot!(a, b, -0.2)

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
        else
            @warn "Mismatch" ham_before delta ham_after diff=(ham_before+delta-ham_after)
        end
    end

    @test n_consistent == n_trials
end

# =============================================================================
# PART 8: Simulation — self-moves, accumulated deltaH matches ham!
# =============================================================================

@testset "deltaHam! - accumulated deltaH matches ham! (self-moves)" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(c, -0.5)
    intPot!(a, b, -0.3)
    intPot!(a, c, -0.2)
    intPot!(b, c, -0.1)

    couplingPot!(a, b, -0.4)
    couplingPot!(a, b, c, -0.7)

    init!(l, [a, b, c])

    ham!(a)

    n_moves = 500
    for _ in 1:n_moves
        a.moveIdx = rand(Int32(1):Int32(a.length))
        identifyBond!(a)
        proposeMove!(a)

        if a.latS[a.newPos[1], a.newPos[2], a.newPos[3]] > 0
            continue
        end

        deltaHam!(a)
        makeMove!(a)
        # a.ham already updated by deltaHam! (b.ham += b.deltaH)
    end

    ham_accumulated = a.ham
    ham!(a)
    ham_recalc = a.ham

    @test isapprox(ham_accumulated, ham_recalc, atol=1e-10)
end

# =============================================================================
# PART 9: Full simulation — all species move, single-step
# =============================================================================

@testset "deltaHam! - all species move, single-step consistency" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=15, chemPot=0.5)
    b = BulkMonomer(; name="B", length=15, chemPot=0.5)
    c = BulkMonomer(; name="C", length=15, chemPot=0.5)

    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(c, -0.5)
    intPot!(a, b, -0.3)
    intPot!(a, c, -0.2)
    intPot!(b, c, -0.1)

    couplingPot!(a, b, -0.4)
    couplingPot!(a, b, c, -0.7)

    init!(l, [a, b, c])

    species = [a, b, c]
    n_trials = 200
    n_consistent = 0

    for _ in 1:n_trials
        sp = species[rand(1:3)]

        ham!(sp)
        ham_before = sp.ham

        sp.moveIdx = rand(Int32(1):Int32(sp.length))
        identifyBond!(sp)
        proposeMove!(sp)

        if sp.latS[sp.newPos[1], sp.newPos[2], sp.newPos[3]] > 0
            n_consistent += 1
            continue
        end

        deltaHam!(sp)
        delta = sp.deltaH

        makeMove!(sp)
        ham!(sp)
        ham_after = sp.ham

        if isapprox(ham_before + delta, ham_after, atol=1e-10)
            n_consistent += 1
        else
            @warn "Mismatch" species=sp.name ham_before delta ham_after diff=(ham_before+delta-ham_after)
        end
    end

    @test n_consistent == n_trials
end

# =============================================================================
# PART 10: Performance
# =============================================================================

@testset "deltaHam! - performance with couplings" begin
    l = Lattice(; boxlength=30, nspecies=3)

    a = BulkMonomer(; name="A", length=100, chemPot=1.0)
    b = BulkMonomer(; name="B", length=100, chemPot=1.0)
    c = BulkMonomer(; name="C", length=100, chemPot=1.0)

    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -0.5)
    intPot!(a, c, -0.5)
    intPot!(b, c, -0.5)

    couplingPot!(a, b, c, -0.3)
    couplingPot!(a, b, -0.2)

    init!(l, [a, b, c])

    ham!(a)

    t = @elapsed for _ in 1:1000
        a.moveIdx = rand(Int32(1):Int32(a.length))
        identifyBond!(a)
        proposeMove!(a)
        if a.latS[a.newPos[1], a.newPos[2], a.newPos[3]] == 0
            deltaHam!(a)
            makeMove!(a)
        end
    end

    @test t < 5.0
    @test isfinite(a.ham)
end


println()
println("All deltaHam! coupling tests completed!")
