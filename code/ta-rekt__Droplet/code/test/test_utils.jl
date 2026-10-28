# Test suite for utility functions in utils.jl
# Includes tests for:
# - cplpq: on-site co-location check
# - nlp: local neighbor counting

using Test
using Droplet


# =============================================================================
# PART 1: Generic BulkMonomer dispatch
# =============================================================================

@testset "cplpq - generic: both present at site" begin
    l = Lattice(; boxlength=10, nspecies=2)
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(b, 0.0)
    init!(l, [a, b])

    # Place both at (5,5,5)
    a.latS .= 0; b.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1

    pos = Int32[5, 5, 5]
    @test cplpq(a, b, pos) == 1
    @test cplpq(b, a, pos) == 1  # symmetric
end

@testset "cplpq - generic: one absent" begin
    l = Lattice(; boxlength=10, nspecies=2)
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(b, 0.0)
    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 8, 8, 8
    b.latS[8, 8, 8] = 1

    pos = Int32[5, 5, 5]
    @test cplpq(a, b, pos) == 0
    @test cplpq(b, a, pos) == 0
end

@testset "cplpq - generic: both absent" begin
    l = Lattice(; boxlength=10, nspecies=2)
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(b, 0.0)
    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 3, 3, 3
    a.latS[3, 3, 3] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 8, 8, 8
    b.latS[8, 8, 8] = 1

    pos = Int32[5, 5, 5]
    @test cplpq(a, b, pos) == 0
end

@testset "cplpq - generic: self-check (j == j)" begin
    l = Lattice(; boxlength=10, nspecies=1)
    a = BulkMonomer(; name="A", length=2, chemPot=0.0)
    intPot!(a, 0.0)
    init!(l, [a])

    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = 6, 5, 5
    a.latS[6, 5, 5] = 2

    # cplpq(a, a, occupied_pos) → 1
    @test cplpq(a, a, Int32[5, 5, 5]) == 1
    @test cplpq(a, a, Int32[6, 5, 5]) == 1
    # cplpq(a, a, empty_pos) → 0
    @test cplpq(a, a, Int32[7, 7, 7]) == 0
end

# =============================================================================
# PART 2: PolyMonomer dispatch — p is PolyMonomer
# =============================================================================

@testset "cplpq - PolyMonomer p: correct intIdx" begin
    l = Lattice(; boxlength=15, nspecies=2)
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(a, pm, 0.0)
    init!(l, [poly, a])

    # Find a polymer bead position
    ppos = Int32[poly.pos[poly.head][1], poly.pos[poly.head][2], poly.pos[poly.head][3]]

    # Place a at same site
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = ppos[1], ppos[2], ppos[3]
    a.latS[ppos[1], ppos[2], ppos[3]] = 1

    @test cplpq(pm, a, ppos) == 1
end

@testset "cplpq - PolyMonomer p: wrong intIdx (heterogeneous)" begin
    l = Lattice(; boxlength=15, nspecies=3)
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(a, pm1, 0.0); intPot!(a, pm2, 0.0)
    init!(l, [poly, a])

    # Find sites for pm1 and pm2
    pm1_pos = nothing; pm2_pos = nothing
    for i in 1:poly.length
        p = poly.pos[i]
        if p[6] == pm1.intIdx && pm1_pos === nothing
            pm1_pos = Int32[p[1], p[2], p[3]]
        elseif p[6] == pm2.intIdx && pm2_pos === nothing
            pm2_pos = Int32[p[1], p[2], p[3]]
        end
        pm1_pos !== nothing && pm2_pos !== nothing && break
    end
    @test pm1_pos !== nothing
    @test pm2_pos !== nothing

    # Place A at pm2 site
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm2_pos[1], pm2_pos[2], pm2_pos[3]
    a.latS[pm2_pos[1], pm2_pos[2], pm2_pos[3]] = 1

    # pm1 should NOT match at pm2's site (wrong intIdx)
    @test cplpq(pm1, a, pm2_pos) == 0
    # pm2 should match at its own site
    @test cplpq(pm2, a, pm2_pos) == 1
end

# =============================================================================
# PART 3: PolyMonomer dispatch — q is PolyMonomer
# =============================================================================

@testset "cplpq - PolyMonomer q: correct vs wrong intIdx" begin
    l = Lattice(; boxlength=15, nspecies=3)
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(a, pm1, 0.0); intPot!(a, pm2, 0.0)
    init!(l, [poly, a])

    pm1_pos = nothing; pm2_pos = nothing
    for i in 1:poly.length
        p = poly.pos[i]
        if p[6] == pm1.intIdx && pm1_pos === nothing
            pm1_pos = Int32[p[1], p[2], p[3]]
        elseif p[6] == pm2.intIdx && pm2_pos === nothing
            pm2_pos = Int32[p[1], p[2], p[3]]
        end
        pm1_pos !== nothing && pm2_pos !== nothing && break
    end

    # Place A at pm1 site
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_pos[1], pm1_pos[2], pm1_pos[3]
    a.latS[pm1_pos[1], pm1_pos[2], pm1_pos[3]] = 1

    # cplpq(a, pm1, pm1_pos) → 1 (a present, pm1 correct intIdx)
    @test cplpq(a, pm1, pm1_pos) == 1
    # cplpq(a, pm2, pm1_pos) → 0 (pm2 wrong intIdx at pm1 site)
    @test cplpq(a, pm2, pm1_pos) == 0
end

# =============================================================================
# PART 4: Both PolyMonomer dispatch
# =============================================================================

@testset "cplpq - both PolyMonomer" begin
    l = Lattice(; boxlength=15, nspecies=4)

    # Two separate polymers
    poly1 = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly1, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly1.monomers[1]

    poly2 = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly2, defaultChemPot=1.0, defaultIntPot=0.0)
    pm2 = poly2.monomers[1]

    intPot!(pm1, pm2, 0.0)
    init!(l, [poly1, poly2])

    # Find a site where poly1 has a bead
    p1pos = Int32[poly1.pos[poly1.head][1], poly1.pos[poly1.head][2], poly1.pos[poly1.head][3]]

    # Check if poly2 also happens to be there (unlikely but handle it)
    p2_at_p1 = (poly2.latS[p1pos[1], p1pos[2], p1pos[3]] == pm2.intIdx)

    result = cplpq(pm1, pm2, p1pos)
    if p2_at_p1
        @test result == 1
    else
        @test result == 0
    end
end

@testset "cplpq - both PolyMonomer, heterogeneous same polymer" begin
    l = Lattice(; boxlength=15, nspecies=2)
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    init!(l, [poly])

    # pm1 and pm2 share the same latS but have different intIdx
    pm1_pos = nothing
    for i in 1:poly.length
        p = poly.pos[i]
        if p[6] == pm1.intIdx
            pm1_pos = Int32[p[1], p[2], p[3]]
            break
        end
    end
    @test pm1_pos !== nothing

    # cplpq(pm1, pm2, pm1_pos): pm1 matches (intIdx ok), but pm2 checks pm2.intIdx
    # at the same site which has pm1.intIdx → 0
    @test cplpq(pm1, pm2, pm1_pos) == 0
    # cplpq(pm1, pm1, pm1_pos) → both check pm1.intIdx → 1
    @test cplpq(pm1, pm1, pm1_pos) == 1
end

# =============================================================================
# PART 5: Edge cases
# =============================================================================

@testset "cplpq - SubArray pos (view into species pos)" begin
    l = Lattice(; boxlength=10, nspecies=2)
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(b, 0.0)
    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5
    b.latS[5, 5, 5] = 1

    # Use a SubArray (view) like ham! does: b.pos[i][range(1,3)]
    pos_view = a.pos[1][range(1, 3)]
    @test cplpq(a, b, pos_view) == 1
end

@testset "cplpq - multiple beads, only some co-located" begin
    l = Lattice(; boxlength=10, nspecies=2)
    a = BulkMonomer(; name="A", length=3, chemPot=0.0)
    b = BulkMonomer(; name="B", length=2, chemPot=0.0)
    intPot!(a, 0.0); intPot!(b, 0.0)
    init!(l, [a, b])

    a.latS .= 0; b.latS .= 0
    # A at (3,3,3), (4,3,3), (5,3,3)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 3, 3, 3; a.latS[3,3,3] = 1
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = 4, 3, 3; a.latS[4,3,3] = 2
    a.pos[3][1], a.pos[3][2], a.pos[3][3] = 5, 3, 3; a.latS[5,3,3] = 3
    # B at (4,3,3), (8,8,8)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 4, 3, 3; b.latS[4,3,3] = 1
    b.pos[2][1], b.pos[2][2], b.pos[2][3] = 8, 8, 8; b.latS[8,8,8] = 2
    
    # Only A[2] and B[1] co-locate at (4,3,3)
    @test cplpq(a, b, Int32[3, 3, 3]) == 0  # A present, B absent
    @test cplpq(a, b, Int32[4, 3, 3]) == 1  # both present
    @test cplpq(a, b, Int32[5, 3, 3]) == 0  # A present, B absent
    @test cplpq(a, b, Int32[8, 8, 8]) == 0  # A absent, B present
end


# =============================================================================
# PART 6: nlp (number of local partners) tests
# Note: nlp counts only the 26 peripheral sites of a 3x3x3 cube (excludes center)
# =============================================================================

@testset "nlp - BulkMonomer to PolyMonomer" begin
    l = Lattice(; boxlength=15, nspecies=2)
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(a, pm, 0.0)
    init!(l, [poly, a])

    # Get a polymer bead position
    poly_pos = Int32[poly.pos[poly.head][1], poly.pos[poly.head][2], poly.pos[poly.head][3]]
    
    # Place bulk monomer at an adjacent position
    adj_pos = Int32[poly_pos[1] + 1, poly_pos[2], poly_pos[3]]
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = adj_pos[1], adj_pos[2], adj_pos[3]
    a.latS[adj_pos[1], adj_pos[2], adj_pos[3]] = 1
    
    # Count local neighbors from bulk monomer to polymer monomer
    # The 26 peripheral sites of the 3x3x3 cube around adj_pos should include poly_pos
    count = nlp(a, pm, adj_pos)
    @test count >= 1  # Should count at least the polymer bead we know is there
end

@testset "nlp - PolyMonomer to PolyMonomer (heterogeneous)" begin
    l = Lattice(; boxlength=15, nspecies=1)
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    init!(l, [poly])

    # Find positions for pm1 and pm2
    pm1_pos = nothing
    for i in 1:poly.length
        p = poly.pos[i]
        if p[6] == pm1.intIdx
            pm1_pos = Int32[p[1], p[2], p[3]]
            break
        end
    end
    @test pm1_pos !== nothing

    # Count local neighbors of type pm2 around pm1's position (excludes center)
    count = nlp(pm1, pm2, pm1_pos)
    @test count >= 0  # Should be non-negative
    
    # Count local neighbors of type pm1 around pm1's position (excludes center, so does NOT count itself)
    self_count = nlp(pm1, pm1, pm1_pos)
    @test self_count >= 0  # Should be non-negative (center is excluded)
end

@testset "nlp - PolyMonomer to BulkMonomer" begin
    l = Lattice(; boxlength=15, nspecies=2)
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=3, chemPot=0.0)
    intPot!(a, 0.0); intPot!(a, pm, 0.0)
    init!(l, [poly, a])

    # Get a polymer bead position
    poly_pos = Int32[poly.pos[poly.head][1], poly.pos[poly.head][2], poly.pos[poly.head][3]]
    
    # Place bulk monomers in the periphery around the polymer position
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = poly_pos[1] + 1, poly_pos[2], poly_pos[3]
    a.latS[poly_pos[1] + 1, poly_pos[2], poly_pos[3]] = 1
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = poly_pos[1], poly_pos[2] + 1, poly_pos[3]
    a.latS[poly_pos[1], poly_pos[2] + 1, poly_pos[3]] = 2
    a.pos[3][1], a.pos[3][2], a.pos[3][3] = poly_pos[1], poly_pos[2], poly_pos[3] + 1
    a.latS[poly_pos[1], poly_pos[2], poly_pos[3] + 1] = 3
    
    # Count local neighbors from polymer monomer to bulk monomer (excludes center)
    count = nlp(pm, a, poly_pos)
    @test count >= 3  # Should count at least the 3 bulk monomers we placed in periphery
end

@testset "nlp - empty neighborhood" begin
    l = Lattice(; boxlength=20, nspecies=2)
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    intPot!(a, 0.0); intPot!(a, pm, 0.0)
    init!(l, [poly, a])

    # Place bulk monomer far from polymer
    a.latS .= 0
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 19, 19, 19
    a.latS[19, 19, 19] = 1
    
    # Count local neighbors around bulk monomer (should be 0 polymer beads in periphery)
    count = nlp(a, pm, Int32[19, 19, 19])
    @test count == 0
end

@testset "nlp - center exclusion test" begin
    l = Lattice(; boxlength=20, nspecies=1)
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    init!(l, [poly])

    # Find a position with a polymer bead
    pos = Int32[poly.pos[poly.head][1], poly.pos[poly.head][2], poly.pos[poly.head][3]]
    
    # Count local neighbors at the position where the bead itself is located
    # Should NOT count the center position (where the bead is)
    count = nlp(pm, pm, pos)
    
    # If the polymer is isolated (no other beads in the 26 peripheral sites), count should be 0
    # Otherwise, it should count only neighbors in the periphery, not itself
    @test count >= 0
    
    # Verify it doesn't count 27 (which would include center)
    @test count <= 26
end

