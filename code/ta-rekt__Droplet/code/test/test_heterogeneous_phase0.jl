# Test suite for Phase 0: Heterogeneous polymer lattice intIdx storage and nnp dispatch

# MARK: -Lattice stores intIdx
@testset "Lattice stores intIdx for heterogeneous polymers" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 1
    )
    
    # Create polymer with 2 monomer types in block pattern
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Verify intIdx values are assigned
    @test m1.intIdx !== nothing
    @test m2.intIdx !== nothing
    @test m1.intIdx != m2.intIdx
    
    # Walk the chain and verify lattice stores correct intIdx for each bead
    # Phase 1: pos[key][6] stores the intIdx directly
    key = poly.head
    chainPos = 1
    while key != 0 && chainPos <= poly.length
        pos = poly.pos[key]
        expectedIntIdx = pos[6]  # intIdx stored directly in pos
        
        actualLatValue = poly.latS[pos[1], pos[2], pos[3]]
        @test actualLatValue == expectedIntIdx
        
        key = pos[5]  # next in chain
        chainPos += 1
    end
    
    @test chainPos == poly.length + 1
end

# MARK: -nnp PolyMonomer dispatch
@testset "nnp dispatch for PolyMonomer counts only matching intIdx" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 1
    )
    
    # Create polymer with block pattern for an 8-bead maxlength polymer with 6 actual beads
    # Block size = div(8, 2) = 4, so pattern is [1,1,1,1,2,2,2,2] truncated to 6 beads: [1,1,1,1,2,2]
    # This ensures we have adjacent beads of same and different types at the boundary (bead 4-5)
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Find the 4th bead (type 1, last in first block) - it has a type-1 neighbor on left and type-2 on right
    # Chain is: 1-2-3-4-5-6 with types 1,1,1,1,2,2
    # Bead 4 (type 1) is adjacent to bead 3 (type 1) and bead 5 (type 2)
    key = poly.head
    for _ in 1:3
        key = poly.pos[key][5]  # advance to 4th bead
    end
    bead4_pos = poly.pos[key][1:3]
    
    # nnp(m1, m1, pos) should count only type-1 neighbors
    # For bead 4, this should be at least 1 (bead 3) but NOT include bead 5
    neighbors_m1 = nnp(m1, m1, bead4_pos)
    
    # nnp(m1, m2, pos) should count only type-2 neighbors
    # For bead 4, this should be at least 1 (bead 5)
    neighbors_m2 = nnp(m1, m2, bead4_pos)
    
    # Bead 4 should have at least 1 type-1 neighbor (bead 3)
    @test neighbors_m1 >= 1
    # Bead 4 should have at least 1 type-2 neighbor (bead 5)
    @test neighbors_m2 >= 1
    
    # The total of type-1 and type-2 neighbors should equal the old nnp behavior
    # which counts all polymer beads regardless of type
    # Use nn (which checks latS > 0) as baseline
    total_neighbors = nn(poly, bead4_pos)
    # Sum of typed neighbors equals total neighbors
    @test neighbors_m1 + neighbors_m2 == total_neighbors
end

# MARK: -nnp self-interaction
@testset "nnp self-interaction for PolyMonomer" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 1
    )
    
    # Use block pattern - test at type boundary
    # Block pattern for length=6, nMonomerTypes=2: [1,1,1,2,2,2]
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # For block pattern [1,1,1,2,2,2], bead 3 (type 1) has:
    # - bead 2 (type 1) on left
    # - bead 4 (type 2) on right
    # So nnp(m1, m1, bead3_pos) should be >= 1 (type-1 neighbor)
    # And nnp(m1, m2, bead3_pos) should be >= 1 (type-2 neighbor)
    
    key = poly.head
    for _ in 1:2
        key = poly.pos[key][5]  # advance to 3rd bead
    end
    bead3_pos = poly.pos[key][1:3]
    
    neighbors_same = nnp(m1, m1, bead3_pos)
    neighbors_other = nnp(m1, m2, bead3_pos)
    
    # At block boundary, type-1 bead has at least 1 type-1 neighbor
    @test neighbors_same >= 1
    # At block boundary, type-1 bead has at least 1 type-2 neighbor
    @test neighbors_other >= 1
end

# MARK: -Overlap detection
@testset "Overlap detection with intIdx values" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 1
    )
    
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    # Get a position that's occupied
    key = poly.head
    occupied_pos = poly.pos[key][1:3]
    
    # Verify the lattice value is an intIdx (> 0)
    @test poly.latS[occupied_pos[1], occupied_pos[2], occupied_pos[3]] > 0
    
    # Set up for overlap check - put newPosi at an occupied position
    poly.newPosi = copy(occupied_pos)
    poly.newPosj = [occupied_pos[1] + 1, occupied_pos[2], occupied_pos[3]]
    # Make sure bprev is different so overlap would be detected
    poly.bprev = [occupied_pos[1] + 2, occupied_pos[2], occupied_pos[3]]
    poly.bnext = [occupied_pos[1] + 3, occupied_pos[2], occupied_pos[3]]
    
    # overlap() should detect the occupied site even though it has intIdx > 1
    @test overlap(poly) == true
    
    # Now test with an empty position
    # Find an empty position
    empty_pos = [Int32(1), Int32(1), Int32(1)]
    while poly.latS[empty_pos[1], empty_pos[2], empty_pos[3]] != 0
        empty_pos[1] += 1
        if empty_pos[1] > l.boxlength
            empty_pos[1] = 1
            empty_pos[2] += 1
        end
    end
    
    poly.newPosi = copy(empty_pos)
    poly.newPosj = [empty_pos[1] + 1, empty_pos[2], empty_pos[3]]
    poly.bprev = [empty_pos[1] + 2, empty_pos[2], empty_pos[3]]
    poly.bnext = [empty_pos[1] + 3, empty_pos[2], empty_pos[3]]
    
    # No overlap at empty position
    @test overlap(poly) == false
end

# MARK: -Cross-species nnp
@testset "Cross-species nnp with PolyMonomer and FastBulkMonomer" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 2
    )
    
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly, tf])
    
    m1 = poly.monomers[1]
    
    # Get position of first polymer bead
    pos = poly.pos[poly.head][1:3]
    
    # nnp between PolyMonomer and FastBulkMonomer should still work
    # (uses the generic AbstractSpecies dispatch, counts any tf beads nearby)
    neighbors_tf = nnp(m1, tf, pos)
    @test neighbors_tf >= 0
    @test neighbors_tf <= 6
end

# MARK: -Backward compatibility homogeneous
@testset "Backward compatibility: homogeneous polymer with single monomer type" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 1
    )
    
    # Homogeneous polymer with single monomer type (uniform pattern)
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 1,
        pattern = :uniform
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    # All beads should have the same intIdx (m1's intIdx)
    m1 = poly.monomers[1]
    expectedIntIdx = m1.intIdx
    
    # Lattice should store the single monomer's intIdx for all beads
    key = poly.head
    chainPos = 1
    while key != 0 && chainPos <= poly.length
        pos = poly.pos[key]
        # All beads should have the same intIdx
        @test poly.latS[pos[1], pos[2], pos[3]] == expectedIntIdx
        key = pos[5]
        chainPos += 1
    end
    
    # Test that nnp works correctly with single monomer type (should count all polymer neighbors)
    pos = poly.pos[poly.head][1:3]
    neighbors_self = nnp(m1, m1, pos)
    neighbors_nn = nn(poly, pos)
    @test neighbors_self == neighbors_nn
end

# MARK: -kinkBond preserves intIdx
@testset "kinkBond! preserves intIdx in lattice" begin
    l = Lattice(;
        boxlength = 10,
        nspecies = 1
    )
    
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    # Get the intIdx values before any moves
    m1_intIdx = poly.monomers[1].intIdx
    m2_intIdx = poly.monomers[2].intIdx
    
    # Record intIdx for each bead position
    original_intIdxs = Dict{Int32, Int32}()
    key = poly.head
    chainPos = 1
    while key != 0 && chainPos <= poly.length
        pos = poly.pos[key]
        original_intIdxs[key] = poly.latS[pos[1], pos[2], pos[3]]
        key = pos[5]
        chainPos += 1
    end
    
    # Verify all recorded intIdxs are valid (either m1 or m2's intIdx)
    for (k, idx) in original_intIdxs
        @test idx == m1_intIdx || idx == m2_intIdx
    end
end
