# Test file for ham! with on-site coupling for BulkMonomer
# Tests the full Hamiltonian calculation including two-way and three-way coupling

using Test
using Droplet

println("Testing ham! with on-site coupling")
println("=" ^ 50)

# =============================================================================
# PART 1: Regression tests
# =============================================================================

@testset "ham! - unchanged when couplingMat is nothing (regression)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=10, chemPot=1.0)
    b = BulkMonomer(; name="B", length=10, chemPot=1.0)

    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, b, -0.5)

    init!(l, [a, b])

    @test a.couplingMat === nothing
    @test b.couplingMat === nothing

    ham!(a)
    ham!(b)

    @test isfinite(a.ham)
    @test isfinite(b.ham)
end

@testset "ham! - zero coupling has no effect" begin
    l = Lattice(; boxlength=15, nspecies=3)

    a = BulkMonomer(; name="A", length=5, chemPot=1.0)
    b = BulkMonomer(; name="B", length=5, chemPot=1.0)
    c = BulkMonomer(; name="C", length=5, chemPot=1.0)

    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -0.3)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    init!(l, [a, b, c])

    ham!(a)
    ham1 = a.ham

    couplingPot!(a, b, c, 0.0)

    ham!(a)
    ham2 = a.ham

    @test isapprox(ham1, ham2, atol=1e-10)
end

# =============================================================================
# PART 2: Three-way basic
# =============================================================================

@testset "ham! - three-way coupling active when all three co-located" begin
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

    # All three at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1

    ham!(a)
    # cplpq(b, c, (5,5,5)) = 1 → coupling = -1.0
    @test isapprox(a.ham, -1.0, atol=1e-10)
end

@testset "ham! - three-way coupling zero when only two co-located" begin
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

    # A and B at (5,5,5), C far away
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 10, 10
    c.latS[10, 10, 10] = 1

    ham!(a)
    # cplpq(b, c, (5,5,5)): b present, c absent → 0
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

@testset "ham! - three-way coupling zero when all separated" begin
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

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 8, 8, 8
    b.latS[8, 8, 8] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 12, 12, 12
    c.latS[12, 12, 12] = 1

    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

# =============================================================================
# PART 3: Two-way basic
# =============================================================================

@testset "ham! - two-way coupling active when co-located" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, b, 0.0)

    couplingPot!(a, b, -1.0)

    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0

    # A and B at same site
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1

    ham!(a)
    # cplpq(b, b, (5,5,5)) = 1 → coupling = -1.0
    @test isapprox(a.ham, -1.0, atol=1e-10)
end

@testset "ham! - two-way coupling zero when apart" begin
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
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 10, 10
    b.latS[10, 10, 10] = 1

    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

@testset "ham! - two-way distinct from intMat (nearest-neighbor)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, b, -0.3)        # nearest-neighbor interaction

    couplingPot!(a, b, -0.5)   # on-site coupling

    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0

    # A at (5,5,5), B at (5,5,5) — co-located
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1

    ham!(a)
    # nnp(a, b, (5,5,5)) counts b at 6 neighbors — b is at (5,5,5), not a neighbor → 0
    # intMat contribution = -0.3 * 0 = 0
    # Coupling = -0.5 * cplpq(b, b, (5,5,5)) = -0.5
    @test isapprox(a.ham, -0.5, atol=1e-10)

    # Now B as neighbor only (not co-located)
    a.latS .= 0; b.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 6, 5, 5
    b.latS[6, 5, 5] = 1

    ham!(a)
    # nnp(a, b, (5,5,5)) = 1 (b at neighbor (6,5,5))
    # intMat = -0.3 * 1 = -0.3
    # Coupling = -0.5 * cplpq(b, b, (5,5,5)) = 0 (b not at (5,5,5))
    @test isapprox(a.ham, -0.3, atol=1e-10)
end

# =============================================================================
# PART 4: Manual calculation verification
# =============================================================================

@testset "ham! - manual calculation (three-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=3, chemPot=1.0)
    b = BulkMonomer(; name="B", length=3, chemPot=0.0)
    c = BulkMonomer(; name="C", length=3, chemPot=0.0)

    intPot!(a, -1.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, -0.2)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -0.5)

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # A beads at (5,5,5), (6,5,5), (7,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = 6, 5, 5
    a.pos[3][1], a.pos[3][2], a.pos[3][3] = 7, 5, 5
    a.latS[5, 5, 5] = 1
    a.latS[6, 5, 5] = 2
    a.latS[7, 5, 5] = 3

    # B at (5,5,5) — co-located with A[1]
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1
    b.pos[2][1], b.pos[2][2], b.pos[2][3] = 15, 15, 15
    b.pos[3][1], b.pos[3][2], b.pos[3][3] = 15, 16, 15
    b.latS[15, 15, 15] = 2
    b.latS[15, 16, 15] = 3

    # C at (5,5,5) — co-located with A[1] and B[1]
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    c.pos[2][1], c.pos[2][2], c.pos[2][3] = 15, 15, 16
    c.pos[3][1], c.pos[3][2], c.pos[3][3] = 15, 16, 16
    c.latS[15, 15, 16] = 2
    c.latS[15, 16, 16] = 3

    ham!(a)

    # Manual calculation:
    # chemPot: -(1.0) * 3 = -3.0
    #
    # nnp for A-A self:
    # A[1]@(5,5,5): neighbors include (6,5,5) with A → 1
    # A[2]@(6,5,5): neighbors include (5,5,5) and (7,5,5) → 2
    # A[3]@(7,5,5): neighbors include (6,5,5) → 1
    # Total A-A nnp = 4, selfFactor=0.5 → -1.0 * 4 * 0.5 = -2.0
    #
    # nnp for A-B:
    # A[1]@(5,5,5): B at (5,5,5) is NOT a neighbor (same site) → 0
    # A[2]@(6,5,5): B at (5,5,5) IS a neighbor → 1
    # A[3]@(7,5,5): no B neighbors → 0
    # Total A-B nnp = 1 → -0.2 * 1 = -0.2
    #
    # Three-way coupling (cplpq):
    # A[1]@(5,5,5): cplpq(B, C, (5,5,5)) = B present AND C present → 1
    # A[2]@(6,5,5): cplpq(B, C, (6,5,5)) = B absent → 0
    # A[3]@(7,5,5): cplpq(B, C, (7,5,5)) = 0
    # Total = 1, coupling = -0.5 * 1 = -0.5
    #
    # Total: -3.0 - 2.0 - 0.2 - 0.5 = -5.7

    expected = -3.0 - 2.0 - 0.2 - 0.5
    @test isapprox(a.ham, expected, atol=1e-10)
end

@testset "ham! - manual calculation (two-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=3, chemPot=0.0)
    b = BulkMonomer(; name="B", length=3, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, b, 0.0)

    couplingPot!(a, b, -0.7)

    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0

    # A at (5,5,5), (6,5,5), (7,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = 6, 5, 5
    a.pos[3][1], a.pos[3][2], a.pos[3][3] = 7, 5, 5
    a.latS[5, 5, 5] = 1; a.latS[6, 5, 5] = 2; a.latS[7, 5, 5] = 3

    # B at (5,5,5) and (7,5,5) — co-located with A[1] and A[3]
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.pos[2][1], b.pos[2][2], b.pos[2][3] = 7, 5, 5
    b.pos[3][1], b.pos[3][2], b.pos[3][3] = 15, 15, 15
    b.latS[5, 5, 5] = 1; b.latS[7, 5, 5] = 2; b.latS[15, 15, 15] = 3

    ham!(a)

    # Two-way coupling:
    # A[1]@(5,5,5): cplpq(B, B, (5,5,5)) = B present → 1
    # A[2]@(6,5,5): cplpq(B, B, (6,5,5)) = B absent → 0
    # A[3]@(7,5,5): cplpq(B, B, (7,5,5)) = B present → 1
    # Total = 2, coupling = -0.7 * 2 = -1.4

    @test isapprox(a.ham, -1.4, atol=1e-10)
end

@testset "ham! - manual calculation (mixed two-way + three-way)" begin
    l = Lattice(; boxlength=20, nspecies=3)

    a = BulkMonomer(; name="A", length=2, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)

    couplingPot!(a, b, -0.3)     # two-way A-B
    couplingPot!(a, b, c, -0.5)  # three-way A-B-C

    init!(l, [a, b, c])

    a.latS .= 0; b.latS .= 0; c.latS .= 0

    # A at (5,5,5) and (6,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = 6, 5, 5
    a.latS[5, 5, 5] = 1; a.latS[6, 5, 5] = 2

    # B at (5,5,5) — co-located with A[1]
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1

    # C at (5,5,5) — co-located with A[1] and B
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1

    ham!(a)

    # Two-way (A-B): cplpq(B,B,(5,5,5))=1, cplpq(B,B,(6,5,5))=0 → -0.3
    # Three-way (A-B-C): cplpq(B,C,(5,5,5))=1, cplpq(B,C,(6,5,5))=0 → -0.5
    # Total: -0.8

    @test isapprox(a.ham, -0.8, atol=1e-10)
end

# =============================================================================
# PART 5: First-two-index symmetry in energy
# =============================================================================

@testset "ham! - first-two-index symmetry: each owner sees same coupling" begin
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

    # All three at same site
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1

    ham!(a); ham!(b); ham!(c)

    # Each owner sees -1.0 from its own perspective
    @test isapprox(a.ham, -1.0, atol=1e-10)
    @test isapprox(b.ham, -1.0, atol=1e-10)
    @test isapprox(c.ham, -1.0, atol=1e-10)
end

@testset "ham! - no double counting in three-way" begin
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

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1

    ham!(a)

    # Should be exactly -1.0, not -2.0 or -6.0
    @test isapprox(a.ham, -1.0, atol=1e-10)
end

# =============================================================================
# PART 6: Multiple coupling terms
# =============================================================================

@testset "ham! - multiple three-way coupling terms sum correctly" begin
    l = Lattice(; boxlength=20, nspecies=4)

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    d = BulkMonomer(; name="D", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(d, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(a, d, 0.0)
    intPot!(b, c, 0.0)
    intPot!(b, d, 0.0)
    intPot!(c, d, 0.0)

    couplingPot!(a, b, c, -0.3)
    couplingPot!(a, b, d, -0.4)

    init!(l, [a, b, c, d])

    a.latS .= 0; b.latS .= 0; c.latS .= 0; d.latS .= 0

    # All four at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5; a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5; b.latS[5, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5; c.latS[5, 5, 5] = 1
    d.pos[1][1], d.pos[1][2], d.pos[1][3] = 5, 5, 5; d.latS[5, 5, 5] = 1

    ham!(a)

    # (b,c): cplpq(B,C,(5,5,5)) = 1 → -0.3
    # (b,d): cplpq(B,D,(5,5,5)) = 1 → -0.4
    # Total: -0.7
    @test isapprox(a.ham, -0.7, atol=1e-10)
end

# =============================================================================
# PART 7: PolyMonomer support
# =============================================================================

@testset "ham! - coupling with PolyMonomer" begin
    l = Lattice(; boxlength=20, nspecies=2)

    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, pm, 0.0)
    intPot!(b, pm, 0.0)
    intPot!(a, b, 0.0)

    couplingPot!(a, b, pm, -0.5)

    init!(l, [poly, a, b])

    # Find a polymer bead position
    poly_pos = [poly.pos[poly.head][1], poly.pos[poly.head][2], poly.pos[poly.head][3]]

    # Clear bulk lattices and place A and B at the polymer bead's site
    a.latS .= 0; b.latS .= 0

    a.pos[1][1], a.pos[1][2], a.pos[1][3] = poly_pos[1], poly_pos[2], poly_pos[3]
    a.latS[poly_pos[1], poly_pos[2], poly_pos[3]] = 1

    b.pos[1][1], b.pos[1][2], b.pos[1][3] = poly_pos[1], poly_pos[2], poly_pos[3]
    b.latS[poly_pos[1], poly_pos[2], poly_pos[3]] = 1

    ham!(a)

    # cplpq(B, pm, poly_pos) = B at site AND pm at site → 1
    @test isapprox(a.ham, -0.5, atol=1e-10)
end

@testset "ham! - heterogeneous polymer discriminates monomer types" begin
    l = Lattice(; boxlength=20, nspecies=3)

    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]

    a = BulkMonomer(; name="A", length=1, chemPot=0.0)

    intPot!(a, 0.0)
    intPot!(a, pm1, 0.0)
    intPot!(a, pm2, 0.0)

    couplingPot!(a, pm1, -0.5)  # two-way: A and pm1 on-site

    init!(l, [poly, a])

    # Find pm1 and pm2 bead positions
    pm1_site = nothing; pm2_site = nothing
    for i in 1:poly.length
        p = poly.pos[i]
        if p[6] == pm1.intIdx && pm1_site === nothing
            pm1_site = [p[1], p[2], p[3]]
        elseif p[6] == pm2.intIdx && pm2_site === nothing
            pm2_site = [p[1], p[2], p[3]]
        end
        pm1_site !== nothing && pm2_site !== nothing && break
    end
    @test pm1_site !== nothing
    @test pm2_site !== nothing

    a.latS .= 0

    # Place A at pm1 site → coupling active
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_site[1], pm1_site[2], pm1_site[3]
    a.latS[pm1_site[1], pm1_site[2], pm1_site[3]] = 1

    ham!(a)
    ham_at_pm1 = a.ham

    # Now place A at pm2 site → coupling should NOT fire (wrong monomer type)
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm2_site[1], pm2_site[2], pm2_site[3]
    a.latS[pm2_site[1], pm2_site[2], pm2_site[3]] = 1

    ham!(a)
    ham_at_pm2 = a.ham

    @test isapprox(ham_at_pm1, -0.5, atol=1e-10)
    @test isapprox(ham_at_pm2, 0.0, atol=1e-10)
end

# =============================================================================
# PART 8: Performance
# =============================================================================

@testset "ham! - performance with couplings" begin
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

    t = @elapsed for _ in 1:10
        ham!(a)
    end

    @test t < 5.0
    @test isfinite(a.ham)
end


println()
println("All ham! coupling tests completed!")
