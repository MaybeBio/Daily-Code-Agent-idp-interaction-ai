# Test suite for couplingPot! — on-site coupling API
# Covers three-way (3-arg), two-way (2-arg), symmetry, validation, PolyMonomer support

using Test
using Droplet


# =============================================================================
# PART 1: Basic functionality
# =============================================================================

@testset "couplingPot! 3-arg - basic set/get" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)

    couplingPot!(a, b, c, -0.5)
    @test getCouplingPot(a, b, c) == -0.5
end

@testset "couplingPot! 2-arg - basic set/get" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)

    couplingPot!(a, b, -0.5)
    @test getCouplingPot(a, b) == -0.5
end

@testset "couplingPot! 3-arg - j,k symmetry" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)

    couplingPot!(a, b, c, -0.5)
    @test getCouplingPot(a, b, c) == getCouplingPot(a, c, b)
    @test getCouplingPot(a, b, c) == -0.5
end

# =============================================================================
# PART 2: First-two-index symmetry (owner symmetry)
# =============================================================================

@testset "3-arg symmetry - all BulkMonomer" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)

    couplingPot!(a, b, c, -0.7)

    # All three BulkMonomers should see the same coupling value
    @test getCouplingPot(a, b, c) == -0.7
    @test getCouplingPot(b, a, c) == -0.7
    @test getCouplingPot(c, a, b) == -0.7

    # j,k symmetry also holds for each owner
    @test getCouplingPot(b, c, a) == -0.7
    @test getCouplingPot(c, b, a) == -0.7
end

@testset "3-arg symmetry - one PolyMonomer" begin
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]

    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)
    intPot!(a, pm, 0.0); intPot!(b, pm, 0.0)

    couplingPot!(a, b, pm, -0.4)

    # Both BulkMonomers see same value
    @test getCouplingPot(a, b, pm) == -0.4
    @test getCouplingPot(b, a, pm) == -0.4

    # j,k symmetry
    @test getCouplingPot(a, pm, b) == -0.4
    @test getCouplingPot(b, pm, a) == -0.4

    # PolyMonomer has no couplingMat
    @test !hasfield(typeof(pm), :couplingMat)
end

@testset "2-arg symmetry" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)

    couplingPot!(a, b, -0.6)

    # Both see the same value
    @test getCouplingPot(a, b) == -0.6
    @test getCouplingPot(b, a) == -0.6
end

@testset "2-arg symmetry - PolyMonomer partner" begin
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=5)
    intPot!(a, -1.0); intPot!(a, pm, 0.0)

    couplingPot!(a, pm, -0.3)
    @test getCouplingPot(a, pm) == -0.3
    # PolyMonomer cannot own coupling → no getCouplingPot(pm, a)
end

# =============================================================================
# PART 3: Validation
# =============================================================================

@testset "couplingPot! - error: self === partner (3-arg)" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)
    @test_throws ErrorException couplingPot!(a, a, b, -0.5)
end

@testset "couplingPot! - error: self === context (3-arg)" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)
    @test_throws ErrorException couplingPot!(a, b, a, -0.5)
end

@testset "couplingPot! - j === k delegates to 2-arg" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)

    # Should NOT error (j == k delegates to two-way)
    couplingPot!(a, b, b, -0.5)
    @test getCouplingPot(a, b) == -0.5
    # Owner symmetry via delegation
    @test getCouplingPot(b, a) == -0.5
end

@testset "couplingPot! - error: self === partner (2-arg)" begin
    a = BulkMonomer(; name="A", length=5)
    intPot!(a, -1.0)
    @test_throws ErrorException couplingPot!(a, a, -0.5)
end

@testset "couplingPot! - error: FastBulkMonomer (3-arg)" begin
    f = FastBulkMonomer(; name="F", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(f, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    @test_throws ErrorException couplingPot!(f, b, c, -0.5)
end

@testset "couplingPot! - error: FastBulkMonomer (2-arg)" begin
    f = FastBulkMonomer(; name="F", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(f, -1.0); intPot!(b, -1.0)
    @test_throws ErrorException couplingPot!(f, b, -0.5)
end

# =============================================================================
# PART 4: Auto-initialization
# =============================================================================

@testset "couplingPot! - auto-initializes intSpecies" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)

    @test a.intSpecies === nothing

    couplingPot!(a, b, c, -0.5)

    @test a.intSpecies !== nothing
    @test b.intSpecies !== nothing
    @test c.intSpecies !== nothing
end

@testset "couplingPot! - auto-init coupling for all BulkMonomer participants" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)

    @test a.couplingIdx === nothing
    @test b.couplingIdx === nothing
    @test c.couplingIdx === nothing

    couplingPot!(a, b, c, -0.5)

    # All three should have coupling structures
    @test a.couplingIdx !== nothing
    @test b.couplingIdx !== nothing
    @test c.couplingIdx !== nothing
    # All share the same couplingMat
    @test a.couplingMat === b.couplingMat
    @test b.couplingMat === c.couplingMat
end

@testset "couplingPot! 2-arg - auto-init coupling for both" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)

    couplingPot!(a, b, -0.5)

    @test a.couplingIdx !== nothing
    @test b.couplingIdx !== nothing
    @test a.couplingMat === b.couplingMat
end

# =============================================================================
# PART 5: Multiple couplings and species
# =============================================================================

@testset "couplingPot! - multiple (j,k) pairs" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    d = BulkMonomer(; name="D", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0); intPot!(d, -1.0)

    couplingPot!(a, b, c, -0.5)
    couplingPot!(a, b, d, -0.3)
    couplingPot!(a, c, d, -0.2)

    @test getCouplingPot(a, b, c) == -0.5
    @test getCouplingPot(a, b, d) == -0.3
    @test getCouplingPot(a, c, d) == -0.2
    @test getCouplingPot(a, c, b) == -0.5  # j,k symmetry
end

@testset "couplingPot! - shared coupling values across owners" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)

    couplingPot!(a, b, c, -0.5)

    # Not independent — all three see the same value
    @test getCouplingPot(a, b, c) == getCouplingPot(b, a, c)
    @test getCouplingPot(b, a, c) == getCouplingPot(c, a, b)
end

@testset "couplingPot! - matrix expansion on new species" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)

    couplingPot!(a, b, c, -0.5)

    # Add new species D
    d = BulkMonomer(; name="D", length=5)
    intPot!(d, -1.0)
    intPot!(a, d, 0.0)

    # Original coupling preserved
    @test getCouplingPot(a, b, c) == -0.5

    # Can add new coupling involving D
    couplingPot!(a, b, d, -0.3)
    @test getCouplingPot(a, b, d) == -0.3
    @test getCouplingPot(a, b, c) == -0.5
end

@testset "couplingPot! - values survive intMat sync" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    couplingPot!(a, b, c, -0.5)

    # Link new species
    d = BulkMonomer(; name="D", length=5)
    e = BulkMonomer(; name="E", length=5)
    intPot!(d, -1.0); intPot!(e, -1.0)
    intPot!(d, e, 0.0)
    intPot!(a, d, 0.0)

    @test getCouplingPot(a, b, c) == -0.5
end

# =============================================================================
# PART 6: PolyMonomer support
# =============================================================================

@testset "couplingPot! - PolyMonomer as partner (3-arg)" begin
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)
    intPot!(a, pm, 0.0); intPot!(b, pm, 0.0)

    couplingPot!(a, pm, b, -0.5)
    @test getCouplingPot(a, pm, b) == -0.5
    @test getCouplingPot(a, b, pm) == -0.5  # j,k symmetry
    @test getCouplingPot(b, a, pm) == -0.5  # owner symmetry
end

@testset "couplingPot! - PolyMonomer as context (3-arg)" begin
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)
    intPot!(a, pm, 0.0); intPot!(b, pm, 0.0)

    couplingPot!(a, b, pm, -0.5)
    @test getCouplingPot(a, b, pm) == -0.5
    @test getCouplingPot(b, a, pm) == -0.5  # owner symmetry
end

@testset "couplingPot! - PolyMonomer as partner (2-arg)" begin
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=5)
    intPot!(a, -1.0); intPot!(a, pm, 0.0)

    couplingPot!(a, pm, -0.3)
    @test getCouplingPot(a, pm) == -0.3
end

# =============================================================================
# PART 7: Edge cases
# =============================================================================

@testset "couplingPot! - zero coupling" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    couplingPot!(a, b, c, 0.0)
    @test getCouplingPot(a, b, c) == 0.0
end

@testset "couplingPot! - overwrite coupling" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    couplingPot!(a, b, c, -0.5)
    @test getCouplingPot(a, b, c) == -0.5
    couplingPot!(a, b, c, -0.8)
    @test getCouplingPot(a, b, c) == -0.8
    # Owner symmetry of overwrite
    @test getCouplingPot(b, a, c) == -0.8
end

@testset "getCouplingPot - returns nothing when not initialized" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    @test getCouplingPot(a, b, c) === nothing
    @test getCouplingPot(a, b) === nothing
end

@testset "getCouplingPot - returns nothing for unknown species" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    d = BulkMonomer(; name="D", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    couplingPot!(a, b, c, -0.5)
    @test getCouplingPot(a, b, d) === nothing
end

@testset "couplingPot! - FastBulkMonomer allowed as partner/context" begin
    a = BulkMonomer(; name="A", length=5)
    f = FastBulkMonomer(; name="F", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(f, -1.0); intPot!(c, -1.0)

    # FastBulkMonomer as partner — should work (only self is guarded)
    couplingPot!(a, f, c, -0.5)
    @test getCouplingPot(a, f, c) == -0.5
end

# =============================================================================
# PART 8: Matrix dimensions
# =============================================================================

@testset "couplingPot! - correct dimensions after 3-arg" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    intPot!(a, b, 0.0); intPot!(a, c, 0.0); intPot!(b, c, 0.0)

    couplingPot!(a, b, c, -0.5)

    # couplingMat should be m × n × n
    # m = 3 (all three BulkMonomers are coupling participants now)
    # n = 3 (species count in intMat)
    @test size(a.couplingMat, 1) == 3
    @test size(a.couplingMat, 2) == size(a.intMat, 1)
    @test size(a.couplingMat, 3) == size(a.intMat, 1)
end

@testset "couplingPot! - correct dimensions after 2-arg" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    intPot!(a, -1.0); intPot!(b, -1.0)

    couplingPot!(a, b, -0.5)

    # m = 2 (both BulkMonomers)
    @test size(a.couplingMat, 1) == 2
    @test size(a.couplingMat, 2) == size(a.intMat, 1)
end

# =============================================================================
# PART 8: Coupling survives index remapping during _syncMatrices!
# =============================================================================

@testset "couplingPot! - coupling value survives later _syncMatrices! (2-arg)" begin
    # Regression: couplingPot! sets value at mainPattern2.intIdx,
    # then a later intPot!/couplingPot! call triggers _syncMatrices!
    # which shifts intIdx. The coupling value must follow.
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)

    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, b, -0.5)

    couplingPot!(a, b, -3.0)
    @test getCouplingPot(a, b) == -3.0

    # Now introduce c, which triggers _syncMatrices! and may shift indices
    intPot!(c, -1.0)
    intPot!(a, c, -0.2)

    # After sync, the coupling value should still be accessible via current intIdx
    @test getCouplingPot(a, b) == -3.0
    @test a.couplingMat[a.couplingIdx, b.intIdx, b.intIdx] == -3.0
end

@testset "couplingPot! - coupling value survives later _syncMatrices! (3-arg)" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    d = BulkMonomer(; name="D", length=5)

    intPot!(a, -1.0); intPot!(b, -1.0); intPot!(c, -1.0)
    intPot!(a, b, -0.5); intPot!(a, c, -0.3)

    couplingPot!(a, b, c, -2.0)
    @test getCouplingPot(a, b, c) == -2.0

    # Introduce d, triggering _syncMatrices!
    intPot!(d, -1.0)
    intPot!(a, d, -0.1)

    @test getCouplingPot(a, b, c) == -2.0
    @test a.couplingMat[a.couplingIdx, b.intIdx, c.intIdx] == -2.0
    @test a.couplingMat[a.couplingIdx, c.intIdx, b.intIdx] == -2.0
end

@testset "couplingPot! - BulkMonomer + PolyMonomer coupling survives new species sync" begin
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]

    tf = BulkMonomer(; name="TF", length=5)
    intPot!(tf, 0.0)
    intPot!(pm1, tf, 0.0)

    couplingPot!(tf, pm2, -3.0)
    @test tf.couplingMat[tf.couplingIdx, pm2.intIdx, pm2.intIdx] == -3.0

    # Now add a new species that triggers _syncMatrices!
    rnap = BulkMonomer(; name="RNAP", length=5)
    intPot!(rnap, 0.0)
    couplingPot!(rnap, pm2, 0.0)

    # TF's coupling to pm2 must survive the index shift
    @test tf.couplingMat[tf.couplingIdx, pm2.intIdx, pm2.intIdx] == -3.0
end

@testset "couplingPot! - multiple couplings survive sequential syncs" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)

    intPot!(a, -1.0); intPot!(b, -1.0)

    couplingPot!(a, b, -2.0)
    @test getCouplingPot(a, b) == -2.0

    # Add c with coupling to b
    c = BulkMonomer(; name="C", length=5)
    intPot!(c, -1.0)
    couplingPot!(a, b, c, -1.5)

    @test getCouplingPot(a, b) == -2.0     # two-way preserved
    @test getCouplingPot(a, b, c) == -1.5  # three-way set

    # Add d with yet another sync
    d = BulkMonomer(; name="D", length=5)
    intPot!(d, -1.0)
    intPot!(a, d, -0.1)

    @test getCouplingPot(a, b) == -2.0     # still preserved
    @test getCouplingPot(a, b, c) == -1.5  # still preserved
end

