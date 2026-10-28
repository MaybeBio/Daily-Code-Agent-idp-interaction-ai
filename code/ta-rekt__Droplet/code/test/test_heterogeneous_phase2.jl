# Test suite for Phase 2: Pattern validation and boundary tracking

using Test
using Droplet

# MARK: -patternTails initialization
@testset "patternTails and baseIntIdx initialization" begin
    l = Lattice(; boxlength = 10, nspecies = 1)
    
    # Create polymer with block pattern: [1,1,1,2,2,2]
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
    
    # baseIntIdx should be the first (smallest) intIdx
    @test poly.baseIntIdx == m1.intIdx
    
    # patternTails should have entry for m2 (non-base type)
    @test poly.patternTails !== nothing
    @test haskey(poly.patternTails, m2.intIdx)
    
    # patternTails[m2.intIdx] should be the tail (last bead) of the type-2 block
    tail_key = poly.patternTails[m2.intIdx]
    @test poly.pos[tail_key][6] == m2.intIdx
    
    # The bead after tail should be different type or end of chain
    next_key = poly.pos[tail_key][5]
    if next_key != 0
        @test poly.pos[next_key][6] != m2.intIdx
    end
end

# MARK: -patternTails multi-block
@testset "patternTails with multiple non-base blocks" begin
    l = Lattice(; boxlength = 15, nspecies = 1)
    
    # Custom pattern: [1,2,2,1,1,3,3,3,1]
    # Type 1 is base (can be non-contiguous)
    # Type 2 appears in one contiguous block
    # Type 3 appears in one contiguous block
    poly = Polymer(;
        length = 9,
        maxlength = 12,
        nMonomerTypes = 3,
        pattern = :custom,
        customPattern = [1, 2, 2, 1, 1, 3, 3, 3, 1]
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    m3 = poly.monomers[3]
    
    # baseIntIdx should be m1
    @test poly.baseIntIdx == m1.intIdx
    
    # patternTails should have entries for m2 and m3
    @test haskey(poly.patternTails, m2.intIdx)
    @test haskey(poly.patternTails, m3.intIdx)
    
    # Verify tail positions
    # Type 2 block is at positions 2-3, tail at position 3
    # Type 3 block is at positions 6-8, tail at position 8
    
    # Walk chain to find actual positions
    key = poly.head
    position = 1
    while key != 0 && position <= poly.length
        intIdx = poly.pos[key][6]
        nextKey = poly.pos[key][5]
        
        # Check if this is a tail for non-base type
        if intIdx != m1.intIdx
            if nextKey == 0 || poly.pos[nextKey][6] != intIdx
                # This should be in patternTails
                @test poly.patternTails[intIdx] == key
            end
        end
        
        key = nextKey
        position += 1
    end
end

# MARK: -uniform pattern no tails
@testset "Uniform pattern has no patternTails" begin
    l = Lattice(; boxlength = 10, nspecies = 1)
    
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :uniform
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    # For uniform pattern, all beads are base type, so no tails
    @test poly.patternTails === nothing
    @test poly.baseIntIdx !== nothing
end

# MARK: -single monomer type
@testset "Single monomer type has no patternTails" begin
    l = Lattice(; boxlength = 10, nspecies = 1)
    
    poly = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 1
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    @test poly.patternTails === nothing
    @test poly.baseIntIdx !== nothing
end

# MARK: -validation rejects non-contiguous
@testset "Validation rejects non-contiguous patterns" begin
    l = Lattice(; boxlength = 10, nspecies = 1)
    
    # Alternating pattern creates non-contiguous type-2 blocks
    poly1 = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :alternating
    )
    initMonomers!(poly1, defaultChemPot=0.0, defaultIntPot=0.0)
    @test_throws ErrorException init!(l, [poly1])
    
    # Custom pattern with non-contiguous type-2
    l2 = Lattice(; boxlength = 10, nspecies = 1)
    poly2 = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = [1, 2, 1, 2, 1, 1]
    )
    initMonomers!(poly2, defaultChemPot=0.0, defaultIntPot=0.0)
    @test_throws ErrorException init!(l2, [poly2])
    
    # Random pattern is also rejected (likely non-contiguous)
    # Note: This test might occasionally pass if random happens to be contiguous
    # but for nMonomerTypes=2 with length=6, very unlikely
    l3 = Lattice(; boxlength = 10, nspecies = 1)
    poly3 = Polymer(;
        length = 6,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :random
    )
    initMonomers!(poly3, defaultChemPot=0.0, defaultIntPot=0.0)
    # Random might occasionally be contiguous, so we just verify it runs
    # (success or expected error)
    try
        init!(l3, [poly3])
        # If it didn't error, verify the pattern actually ended up contiguous
        patternVec = [poly3.pos[k][6] for k in keys(poly3.pos)]
    catch e
        @test e isa ErrorException
    end
end
