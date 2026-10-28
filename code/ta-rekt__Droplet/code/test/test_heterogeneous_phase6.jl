# Test suite for Phase 6: Energy accounting for heterogeneous polymers

using Test
using Droplet
using Random

# =============================================================================
# PART 1: Manual ham! computation tests
# =============================================================================

# Helper to perform a move and update length
function doMove!(p::Polymer)
    if p.moveType == 1
        addBond!(p)
    elseif p.moveType == 2
        kinkBond!(p)
    elseif p.moveType == 3
        removeBond!(p)
    end
    p.length = length(p.pos)  # Update length after any move
end

# Helper to manually compute polymer Hamiltonian
function manualHam(p::Polymer)
    totalEnergy = 0.0
    
    # Walk the chain
    key = p.head
    while key != 0
        pos = p.pos[key][1:3]
        intIdx = p.pos[key][6]
        m = p.monomers[intIdx]
        
        # Chemical potential
        totalEnergy -= m.chemPot
        
        # Self-interactions: count same-type neighbors
        for (sp, idx) in m.intSpecies
            if m.intIdx == idx
                # Count only same-type neighbors
                nnSameType = nnp(m, m, pos)
                # Divide by 2 to avoid double-counting
                totalEnergy += m.intPot * nnSameType / 2.0
            else
                # Cross-species interaction
                nnCross = nnp(m, sp, pos)
                totalEnergy += m.intMat[m.intIdx, idx] * nnCross
            end
        end
        
        key = p.pos[key][5]
    end
    
    return totalEnergy
end

# # Helper to get pattern as vector
# function getPatternVector(p::Polymer)
#     result = Int32[]
#     key = p.head
#     while key != 0
#         push!(result, p.pos[key][6])
#         key = p.pos[key][5]
#     end
#     return result
# end

# MARK: -Manual ham! small polymer
@testset "ham! matches manual computation - small homogeneous polymer" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 1,
        pattern = :uniform
    )
    
    initMonomers!(poly, defaultChemPot=1.5, defaultIntPot=-0.5)
    init!(l, [poly])
    
    ham!(poly)
    manualH = manualHam(poly)
    
    @test isapprox(poly.ham, manualH, atol=1e-10)
end

@testset "ham! matches manual computation - small heterogeneous polymer" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Pattern: [1,1,2,2] - 2 type-1, 2 type-2
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    init!(l, [poly])
    
    ham!(poly)
    manualH = manualHam(poly)
    
    @test isapprox(poly.ham, manualH, atol=1e-10)
end

@testset "ham! matches manual - different chemPot per type" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    poly = Polymer(;
        length = 6,
        maxlength = 12,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    # Type 1: chemPot=1.0, Type 2: chemPot=2.0
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=-0.5)
    poly.monomers[poly.monomers[first(keys(poly.monomers))].intIdx].chemPot = 1.0
    
    # Find the second monomer type
    for (intIdx, m) in poly.monomers
        if intIdx != first(keys(poly.monomers))
            m.chemPot = 2.0
        end
    end
    
    init!(l, [poly])
    
    ham!(poly)
    manualH = manualHam(poly)
    
    @test isapprox(poly.ham, manualH, atol=1e-10)
end

@testset "ham! matches manual - different intPot per type" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    poly = Polymer(;
        length = 6,
        maxlength = 12,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    # Type 1: intPot=-0.5, Type 2: intPot=-1.0
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    
    m1intIdx = nothing
    for (intIdx, m) in poly.monomers
        if m1intIdx === nothing
            m1intIdx = intIdx
            m.intPot = -0.5
            m.intMat[intIdx, intIdx] = -0.5
        else
            m.intPot = -1.0
            m.intMat[intIdx, intIdx] = -1.0
        end
    end
    
    init!(l, [poly])
    
    ham!(poly)
    manualH = manualHam(poly)
    
    @test isapprox(poly.ham, manualH, atol=1e-10)
end

# MARK: -ham! with bulk interactions
@testset "ham! matches manual - polymer with bulk monomer interactions" begin
    l = Lattice(; boxlength = 20, nspecies = 2)
    
    poly = Polymer(;
        length = 6,
        maxlength = 12,
        nMonomerTypes = 1,
        pattern = :uniform
    )
    
    bulk = BulkMonomer(;
        length = 50,
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set polymer-bulk interaction BEFORE init!
    m = first(values(poly.monomers))
    intPot!(m, bulk, -0.2)
    
    init!(l, [poly, bulk])
    
    ham!(poly)
    manualH = manualHam(poly)
    
    @test isapprox(poly.ham, manualH, atol=1e-10)
end

@testset "ham! matches manual - heterogeneous polymer with bulk interactions" begin
    l = Lattice(; boxlength = 20, nspecies = 2)
    
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    bulk = BulkMonomer(;
        length = 50,
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set different polymer-bulk interactions for each type BEFORE init!
    for (intIdx, m) in poly.monomers
        if intIdx == poly.baseIntIdx
            intPot!(m, bulk, -0.1)
        else
            intPot!(m, bulk, -0.4)
        end
    end
    
    init!(l, [poly, bulk])
    
    ham!(poly)
    manualH = manualHam(poly)
    
    @test isapprox(poly.ham, manualH, atol=1e-10)
end

# =============================================================================
# PART 2: deltaHam! consistency tests
# =============================================================================

# Helper to set up a polymer state for testing
function setupTestPolymer(length::Int, nTypes::Int, pattern::Symbol; 
                          boxlength::Int=30, customPattern=nothing)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    
    if customPattern !== nothing
        poly = Polymer(;
            length = length,
            maxlength = length * 2,
            nMonomerTypes = nTypes,
            pattern = :custom,
            customPattern = customPattern
        )
    else
        poly = Polymer(;
            length = length,
            maxlength = length * 2,
            nMonomerTypes = nTypes,
            pattern = pattern
        )
    end
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    init!(l, [poly])
    
    return l, poly
end

# MARK: -deltaHam add homogeneous
@testset "deltaHam! consistency - addBond homogeneous polymer" begin
    l, poly = setupTestPolymer(6, 1, :uniform)
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        # Try to find a valid add move
        for attempt in 1:20
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 1 && !overlap(poly) && poly.length < poly.maxlength
                # Valid add move found
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                # Actually perform the add
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - addBond heterogeneous polymer in base region" begin
    # Pattern: [1,1,1,1,2,2,2,2] - add in type-1 region (no swap)
    l, poly = setupTestPolymer(8, 2, :block)
    
    m1 = poly.monomers[poly.baseIntIdx]
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        # Find the head (type-1) bond
        key = poly.head
        next_key = poly.pos[key][5]
        
        # Check it's base type
        if poly.pos[key][6] != poly.baseIntIdx
            continue
        end
        
        # Set up for add at head
        poly.ki = key
        poly.kj = next_key
        poly.kprev = poly.pos[key][4]
        poly.knext = poly.pos[next_key][5]
        
        poly.bi = view(poly.pos[poly.ki], 1:3)
        poly.bj = view(poly.pos[poly.kj], 1:3)
        _sentinel = zeros(Int32, 3)
        poly.bprev = poly.kprev == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.kprev], 1:3)
        poly.bnext = view(poly.pos[poly.knext], 1:3)
        
        # Propose a move manually
        poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
        poly.newPosj = [poly.bj[1], poly.bj[2], Int32(mod1(poly.bj[3] + 1, poly.boxlength))]
        
        if poly.latS[poly.newPosi...] != 0 || poly.latS[poly.newPosj...] != 0
            continue
        end
        
        poly.moveType = 1
        
        deltaHam!(poly)
        expectedHam = initialHam + poly.deltaH
        
        doMove!(poly)
        ham!(poly)
        
        @test isapprox(poly.ham, expectedHam, atol=1e-10)
        break
    end
end

@testset "deltaHam! consistency - addBond heterogeneous polymer in pattern region" begin
    # Pattern: [1,1,1,1,2,2,2,2] - add in type-2 region (with swap)
    l, poly = setupTestPolymer(8, 2, :block)
    
    m2intIdx = nothing
    for (intIdx, m) in poly.monomers
        if intIdx != poly.baseIntIdx
            m2intIdx = intIdx
            break
        end
    end
    
    ham!(poly)
    initialHam = poly.ham
    
    # Find a bond in type-2 region (beads 5-6)
    key = poly.head
    for _ in 1:4
        key = poly.pos[key][5]
    end
    
    @test poly.pos[key][6] == m2intIdx
    next_key = poly.pos[key][5]
    @test poly.pos[next_key][6] == m2intIdx
    
    # Set up for add
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    poly.bprev = view(poly.pos[poly.kprev], 1:3)
    _sentinel = zeros(Int32, 3)
    poly.bnext = poly.knext == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.knext], 1:3)
    
    # Find valid positions for new beads
    for dz in 1:5
        poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + dz, poly.boxlength))]
        poly.newPosj = [poly.bj[1], poly.bj[2], Int32(mod1(poly.bj[3] + dz, poly.boxlength))]
        
        if poly.latS[poly.newPosi...] == 0 && poly.latS[poly.newPosj...] == 0
            break
        end
    end
    
    poly.moveType = 1
    
    deltaHam!(poly)
    expectedHam = initialHam + poly.deltaH
    
    doMove!(poly)
    ham!(poly)
    
    @test isapprox(poly.ham, expectedHam, atol=1e-10)
end

# MARK: -deltaHam remove homogeneous
@testset "deltaHam! consistency - removeBond homogeneous polymer" begin
    l, poly = setupTestPolymer(8, 1, :uniform)
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        # Try to find a valid remove move
        for attempt in 1:20
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 3 && !overlap(poly) && poly.length > 4
                # Valid remove move found
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                # Actually perform the remove
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - removeBond heterogeneous in base region" begin
    # Pattern: [1,1,1,1,2,2,2,2,1,1] with base beads at end for removal
    l, poly = setupTestPolymer(10, 2, :custom, customPattern=[1,1,1,1,2,2,2,2,1,1])
    
    ham!(poly)
    initialHam = poly.ham
    
    # Find the first type-1 bond (beads 1-2)
    key = poly.head
    next_key = poly.pos[key][5]
    
    @test poly.pos[key][6] == poly.baseIntIdx
    @test poly.pos[next_key][6] == poly.baseIntIdx
    
    # Set up for remove
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    _sentinel = zeros(Int32, 3)
    poly.bprev = poly.kprev == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.kprev], 1:3)
    poly.bnext = view(poly.pos[poly.knext], 1:3)
    
    poly.newPosi = Vector{Int32}(poly.bprev)
    poly.newPosj = Vector{Int32}(poly.bnext)
    
    poly.moveType = 3
    
    deltaHam!(poly)
    expectedHam = initialHam + poly.deltaH
    
    doMove!(poly)
    ham!(poly)
    
    @test isapprox(poly.ham, expectedHam, atol=1e-10)
end

@testset "deltaHam! consistency - removeBond heterogeneous in pattern region" begin
    # Pattern: [1,1,1,1,2,2,2,2,1,1] - remove in type-2 region (with swap)
    l, poly = setupTestPolymer(10, 2, :custom, customPattern=[1,1,1,1,2,2,2,2,1,1])
    
    m2intIdx = nothing
    for (intIdx, m) in poly.monomers
        if intIdx != poly.baseIntIdx
            m2intIdx = intIdx
            break
        end
    end
    
    ham!(poly)
    initialHam = poly.ham
    
    # Find a bond in type-2 region (beads 5-6)
    key = poly.head
    for _ in 1:4
        key = poly.pos[key][5]
    end
    
    @test poly.pos[key][6] == m2intIdx
    next_key = poly.pos[key][5]
    @test poly.pos[next_key][6] == m2intIdx
    
    # Set up for remove
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    poly.bprev = view(poly.pos[poly.kprev], 1:3)
    poly.bnext = view(poly.pos[poly.knext], 1:3)
    
    poly.newPosi = Vector{Int32}(poly.bprev)
    poly.newPosj = Vector{Int32}(poly.bnext)
    
    poly.moveType = 3
    
    deltaHam!(poly)
    expectedHam = initialHam + poly.deltaH
    
    doMove!(poly)
    ham!(poly)
    
    @test isapprox(poly.ham, expectedHam, atol=1e-10)
end

# MARK: -deltaHam kink
@testset "deltaHam! consistency - kink move homogeneous polymer" begin
    l, poly = setupTestPolymer(8, 1, :uniform)
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        # Try to find a valid kink move
        for attempt in 1:50
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 2 && !overlap(poly)
                # Valid kink move found
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                # Actually perform the kink
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - kink move in heterogeneous base region" begin
    # Use block pattern [1,1,1,1,2,2,2,2] - kink in type-1 region
    l, poly = setupTestPolymer(8, 2, :block)
    
    for trial in 1:20
        ham!(poly)
        initialHam = poly.ham
        
        # Try to find a valid kink move in base region
        for attempt in 1:50
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            # Only consider base-type bonds
            if poly.pos[poly.ki][6] != poly.baseIntIdx
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 2 && !overlap(poly)
                # Valid kink move found in base region
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - kink move in heterogeneous pattern region" begin
    # Use block pattern [1,1,1,1,2,2,2,2] - kink in type-2 region
    l, poly = setupTestPolymer(8, 2, :block)
    
    m2intIdx = nothing
    for (intIdx, m) in poly.monomers
        if intIdx != poly.baseIntIdx
            m2intIdx = intIdx
            break
        end
    end
    
    for trial in 1:20
        ham!(poly)
        initialHam = poly.ham
        
        # Try to find a valid kink move in pattern region
        for attempt in 1:50
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            # Only consider pattern-type bonds (type-2)
            if poly.pos[poly.ki][6] != m2intIdx
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 2 && !overlap(poly)
                # Valid kink move found in pattern region
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

# MARK: -deltaHam with bulk
@testset "deltaHam! consistency - addBond with bulk interactions" begin
    l = Lattice(; boxlength = 20, nspecies = 2)
    
    poly = Polymer(;
        length = 6,
        maxlength = 12,
        nMonomerTypes = 1,
        pattern = :uniform
    )
    
    bulk = BulkMonomer(;
        length = 30,
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set polymer-bulk interaction BEFORE init!
    m = first(values(poly.monomers))
    intPot!(m, bulk, -0.2)
    
    init!(l, [poly, bulk])
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        for attempt in 1:20
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 1 && !overlap(poly) && poly.length < poly.maxlength
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - removeBond with bulk interactions" begin
    l = Lattice(; boxlength = 20, nspecies = 2)
    
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 1,
        pattern = :uniform
    )
    
    bulk = BulkMonomer(;
        length = 30,
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set polymer-bulk interaction BEFORE init!
    m = first(values(poly.monomers))
    intPot!(m, bulk, -0.2)
    
    init!(l, [poly, bulk])
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        for attempt in 1:20
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 3 && !overlap(poly) && poly.length > 4
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - kink with bulk interactions" begin
    l = Lattice(; boxlength = 20, nspecies = 2)
    
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 1,
        pattern = :uniform
    )
    
    bulk = BulkMonomer(;
        length = 30,
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set polymer-bulk interaction BEFORE init!
    m = first(values(poly.monomers))
    intPot!(m, bulk, -0.2)
    
    init!(l, [poly, bulk])
    
    for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        for attempt in 1:50
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            proposeMove!(poly)
            
            if poly.moveType == 2 && !overlap(poly)
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                doMove!(poly)
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end

@testset "deltaHam! consistency - heterogeneous polymer with bulk" begin
    l = Lattice(; boxlength = 20, nspecies = 2)
    
    poly = Polymer(;
        length = 10,
        maxlength = 20,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = [1,1,1,1,2,2,2,2,1,1]
    )
    
    bulk = BulkMonomer(;
        length = 30,
    )
    
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set different polymer-bulk interactions for each type BEFORE init!
    # Note: poly.baseIntIdx is not set until init!, so use intIdx == 1 for base type
    # (base type is always the first type in the pattern, which gets intIdx=1)
    for (intIdx, m) in poly.monomers
        if intIdx == 1  # base type
            intPot!(m, bulk, -0.1)
        else
            intPot!(m, bulk, -0.4)
        end
    end
    
    init!(l, [poly, bulk])
    
    for trial in 1:20
        ham!(poly)
        initialHam = poly.ham
        
        for attempt in 1:30
            poly.moveIdx = rand(1:poly.length - 1)
            identifyBond!(poly)
            
            # Skip if we selected the tail (no valid bond)
            if poly.kj == 0
                continue
            end
            
            # Skip boundary bonds (different types)
            if poly.pos[poly.ki][6] != poly.pos[poly.kj][6]
                continue
            end
            
            proposeMove!(poly)
            
            if !overlap(poly)
                # Check length constraints for add/remove
                if (poly.moveType == 1 && poly.length >= poly.maxlength) ||
                   (poly.moveType == 3 && poly.length <= 4)
                    continue
                end
                
                deltaHam!(poly)
                expectedHam = initialHam + poly.deltaH
                
                doMove!(poly)
                
                ham!(poly)
                
                @test isapprox(poly.ham, expectedHam, atol=1e-10)
                break
            end
        end
    end
end
# MARK: -Extensive random moves
@testset "deltaHam! consistency - many random moves homogeneous" begin
    l, poly = setupTestPolymer(10, 1, :uniform)
    
    successCount = 0
    for trial in 1:100
        ham!(poly)
        initialHam = poly.ham
        
        poly.moveIdx = rand(1:poly.length - 1)
        identifyBond!(poly)
        
        # Skip if we selected the tail (no valid bond)
        if poly.kj == 0
            continue
        end
        
        proposeMove!(poly)
        
        if !overlap(poly)
            # Check length constraints
            if (poly.moveType == 1 && poly.length >= poly.maxlength) ||
               (poly.moveType == 3 && poly.length <= 4)
                continue
            end
            
            deltaHam!(poly)
            expectedHam = initialHam + poly.deltaH
            
            doMove!(poly)
            
            ham!(poly)
            
            @test isapprox(poly.ham, expectedHam, atol=1e-10)
            successCount += 1
        end
    end
    @test successCount > 0  # At least some moves succeeded
end

@testset "deltaHam! consistency - many random moves heterogeneous" begin
    # Run multiple independent trials with fresh polymers
    totalSuccess = 0
    for round in 1:10
        l, poly = setupTestPolymer(12, 2, :custom, 
                                   customPattern=[1,1,1,1,2,2,2,2,1,1,1,1])
        
        for trial in 1:10
        ham!(poly)
        initialHam = poly.ham
        
        poly.moveIdx = rand(1:poly.length - 1)
        identifyBond!(poly)
        
        # Skip if we selected the tail (no valid bond)
        if poly.kj == 0
            continue
        end
        
        # Skip boundary bonds (different types) for heterogeneous
        if poly.pos[poly.ki][6] != poly.pos[poly.kj][6]
            continue
        end
        
        proposeMove!(poly)
        
        if !overlap(poly)
            # Check length constraints
            if (poly.moveType == 1 && poly.length >= poly.maxlength) ||
               (poly.moveType == 3 && poly.length <= 4)
                continue
            end
            
            deltaHam!(poly)
            expectedHam = initialHam + poly.deltaH
            
            
            doMove!(poly)
            
            ham!(poly)
            
            @test isapprox(poly.ham, expectedHam, atol=1e-10)
            totalSuccess += 1
        end
    end
    end  # round loop
    @test totalSuccess > 0
end

Random.seed!(42)

@testset "deltaHam! consistency - specific seed exposing cross-energy bug" begin
    l = Lattice(; boxlength = 20, nspecies = 1)
    poly = Polymer(;
        length = 10,
        maxlength = 30,
        nMonomerTypes = 2,
        initPos = [1,10,10],
        pattern = :custom,
        customPattern = [1,1,1,1,2,2,2,2,1,1]
    )
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
    intPot!(poly.monomers[1], -0.3)
    intPot!(poly.monomers[2], -0.5)
    init!(l, [poly])

    frame = 0

    for i in 1:10000
        ham!(poly)
        initialHam = poly.ham

        if move!(poly)
            frame += 1
            ham!(poly)
            
            @test isapprox(poly.ham, initialHam + poly.deltaH, atol=1e-10)
            
        end
    end
end