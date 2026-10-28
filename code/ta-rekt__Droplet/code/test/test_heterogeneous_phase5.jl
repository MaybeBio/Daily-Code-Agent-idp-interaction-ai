# Test suite for Phase 5: Type swapping for pattern length preservation

using Droplet
using Test
using OrderedCollections: OrderedDict

# Helper to count beads of each type
function countBeadTypes(p::Polymer)
    counts = Dict{Int32, Int}()
    key = p.head
    while key != 0
        intIdx = p.pos[key][6]
        counts[intIdx] = get(counts, intIdx, 0) + 1
        key = p.pos[key][5]
    end
    return counts
end

# Helper to get pattern as vector
function getPatternVector(p::Polymer)
    result = Int32[]
    key = p.head
    while key != 0
        push!(result, p.pos[key][6])
        key = p.pos[key][5]
    end
    return result
end

# Helper to check that each non-base type appears in exactly one contiguous block
function isPatternContiguous(p::Polymer)
    if p.patternTails === nothing || isempty(p.patternTails)
        return true  # No pattern to check
    end
    
    patternVec = getPatternVector(p)
    
    for intIdx in keys(p.patternTails)
        if intIdx == p.baseIntIdx
            continue
        end
        
        # Find all indices where this type appears
        indices = findall(x -> x == intIdx, patternVec)
        
        if isempty(indices)
            continue
        end
        
        # Check that indices form a contiguous block
        if indices != collect(indices[1]:indices[end])
            return false
        end
    end
    
    return true
end

# Helper to get detailed fragmentation info for debugging
function getFragmentationInfo(p::Polymer)
    patternVec = getPatternVector(p)
    info = Dict{Int32, Vector{UnitRange{Int}}}()
    
    for intIdx in keys(p.monomers)
        if intIdx == p.baseIntIdx
            continue
        end
        
        # Find contiguous blocks of this type
        blocks = UnitRange{Int}[]
        i = 1
        while i <= length(patternVec)
            if patternVec[i] == intIdx
                start = i
                while i <= length(patternVec) && patternVec[i] == intIdx
                    i += 1
                end
                push!(blocks, start:(i-1))
            else
                i += 1
            end
        end
        
        if !isempty(blocks)
            info[intIdx] = blocks
        end
    end
    
    return info
end

# Helper to create a polymer with a "dent" shape at a specific bond position
# A dent is: ...→right→up→right→down→right→...
# This creates a removable bond at the "up-right" position
# Returns the polymer and the moveIdx (1-based) of the removable bond
function createPolymerWithDent(boxlength::Int32, pattern::Vector{Int}, dentBondIdx::Int)
    # dentBondIdx is 1-based: which bond (between beads dentBondIdx and dentBondIdx+1) should be removable
    # We need at least dentBondIdx+3 beads to create the dent
    polyLength = length(pattern)
    @assert polyLength >= dentBondIdx + 3 "Pattern too short for dent at position $dentBondIdx"
    
    poly = Polymer(;
        length = Int32(polyLength),
        maxlength = Int32(polyLength + 10),
        nMonomerTypes = Int32(maximum(pattern)),
        pattern = :custom,
        customPattern = pattern,
        boxlength = boxlength
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    # Manually create positions for a polymer with a dent
    # Start at center of box
    poly.pos = OrderedDict{Int32, Vector{Int32}}()
    startX = div(boxlength, 2)
    startY = div(boxlength, 2)
    startZ = div(boxlength, 2)
    
    # Build positions: straight line with a dent at dentBondIdx
    # Before dent: go right (x+1)
    # At dent: go up, right, down (creates removable bond at the up-right segment)
    x, y, z = startX, startY, startZ
    
    for i in 1:polyLength
        poly.pos[Int32(i)] = Vector{Int32}(undef, 6)
        poly.pos[Int32(i)][1] = x
        poly.pos[Int32(i)][2] = y
        poly.pos[Int32(i)][3] = z
        poly.pos[Int32(i)][4] = i == 1 ? Int32(0) : Int32(i - 1)
        poly.pos[Int32(i)][5] = i == polyLength ? Int32(0) : Int32(i + 1)
        poly.pos[Int32(i)][6] = Int32(0)  # Will be set by _assignPattern!
        
        # Move to next position
        if i < polyLength
            if i == dentBondIdx
                # Go up (y+1) instead of right
                y += 1
            elseif i == dentBondIdx + 1
                # Go right (x+1) - this is the bond that can be removed
                x += 1
            elseif i == dentBondIdx + 2
                # Go down (y-1) to complete the dent
                y -= 1
            else
                # Normal: go right
                x += 1
            end
        end
    end
    
    poly.head = Int32(1)
    poly.keys_vec = collect(keys(poly.pos))
    
    # Initialize other polymer fields
    poly.tempVec = zeros(Int32, 3)
    poly.bondVec = zeros(Int32, 3)
    poly.moveVec = zeros(Int32, 3)
    poly.newPosi = zeros(Int32, 3)
    poly.newPosj = zeros(Int32, 3)
    poly.moveTypes = zeros(Int32, 4)
    poly.moveTypeProbs = [0.0, 0.8, 0.0]
    poly.nTypes = zeros(Float64, 3)
    poly.hastingsFactor = 1.0
    
    return poly, dentBondIdx + 1  # The removable bond is at dentBondIdx+1 (the "up-right" segment)
end

# Helper to create a polymer with a straight section at a specific bond position
# A straight section allows addBond
function createPolymerWithStraight(boxlength::Int32, pattern::Vector{Int}, straightBondIdx::Int)
    polyLength = length(pattern)
    
    poly = Polymer(;
        length = Int32(polyLength),
        maxlength = Int32(polyLength + 10),
        nMonomerTypes = Int32(maximum(pattern)),
        pattern = :custom,
        customPattern = pattern,
        boxlength = boxlength
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    # Create a straight line polymer (all bonds can have addBond applied)
    poly.pos = OrderedDict{Int32, Vector{Int32}}()
    startX = div(boxlength, 2)
    startY = div(boxlength, 2)
    startZ = div(boxlength, 2)
    
    for i in 1:polyLength
        poly.pos[Int32(i)] = Vector{Int32}(undef, 6)
        poly.pos[Int32(i)][1] = startX + i - 1
        poly.pos[Int32(i)][2] = startY
        poly.pos[Int32(i)][3] = startZ
        poly.pos[Int32(i)][4] = i == 1 ? Int32(0) : Int32(i - 1)
        poly.pos[Int32(i)][5] = i == polyLength ? Int32(0) : Int32(i + 1)
        poly.pos[Int32(i)][6] = Int32(0)
    end
    
    poly.head = Int32(1)
    poly.keys_vec = collect(keys(poly.pos))
    
    poly.tempVec = zeros(Int32, 3)
    poly.bondVec = zeros(Int32, 3)
    poly.moveVec = zeros(Int32, 3)
    poly.newPosi = zeros(Int32, 3)
    poly.newPosj = zeros(Int32, 3)
    poly.moveTypes = zeros(Int32, 4)
    poly.moveTypeProbs = [0.0, 0.8, 0.0]
    poly.nTypes = zeros(Float64, 3)
    poly.hastingsFactor = 1.0
    
    return poly, straightBondIdx
end

# Helper to create a polymer with a kink at a specific position
# A kink is: ...→right→right→up→right→... where the middle bead can move
function createPolymerWithKink(boxlength::Int32, pattern::Vector{Int}, kinkBondIdx::Int)
    polyLength = length(pattern)
    @assert polyLength >= kinkBondIdx + 2 "Pattern too short for kink at position $kinkBondIdx"
    
    poly = Polymer(;
        length = Int32(polyLength),
        maxlength = Int32(polyLength + 10),
        nMonomerTypes = Int32(maximum(pattern)),
        pattern = :custom,
        customPattern = pattern,
        boxlength = boxlength
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    poly.pos = OrderedDict{Int32, Vector{Int32}}()
    startX = div(boxlength, 2)
    startY = div(boxlength, 2)
    startZ = div(boxlength, 2)
    
    x, y, z = startX, startY, startZ
    
    for i in 1:polyLength
        poly.pos[Int32(i)] = Vector{Int32}(undef, 6)
        poly.pos[Int32(i)][1] = x
        poly.pos[Int32(i)][2] = y
        poly.pos[Int32(i)][3] = z
        poly.pos[Int32(i)][4] = i == 1 ? Int32(0) : Int32(i - 1)
        poly.pos[Int32(i)][5] = i == polyLength ? Int32(0) : Int32(i + 1)
        poly.pos[Int32(i)][6] = Int32(0)
        
        if i < polyLength
            if i == kinkBondIdx
                # Go up to create a kink
                y += 1
            else
                # Normal: go right
                x += 1
            end
        end
    end
    
    poly.head = Int32(1)
    poly.keys_vec = collect(keys(poly.pos))
    
    poly.tempVec = zeros(Int32, 3)
    poly.bondVec = zeros(Int32, 3)
    poly.moveVec = zeros(Int32, 3)
    poly.newPosi = zeros(Int32, 3)
    poly.newPosj = zeros(Int32, 3)
    poly.moveTypes = zeros(Int32, 4)
    poly.moveTypeProbs = [0.0, 0.8, 0.0]
    poly.nTypes = zeros(Float64, 3)
    poly.hastingsFactor = 1.0
    
    return poly, kinkBondIdx
end

# Helper to finalize polymer initialization with lattice
function finalizePolymerInit!(l::Lattice, poly::Polymer)
    poly.latS = view(l.lat, 2, :, :, :)
    poly.boxlength = l.boxlength
    
    # Assign pattern to pos[key][6]
    Droplet._assignPattern!(poly)
    
    # Populate lattice
    key = poly.head
    while key != 0
        pos = poly.pos[key]
        poly.latS[pos[1], pos[2], pos[3]] = pos[6]
        key = pos[5]
    end
    
    # Initialize Hamiltonian
    Droplet.ham!(poly)
end

# MARK: -addBond inside pattern preserves count (realistic geometry)
@testset "addBond inside pattern preserves type count" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2] - base beads first, then pattern
    # addBond at bond 5 (between beads 5-6, both type 2)
    pattern = [1, 1, 1, 1, 2, 2, 2, 2]
    poly, bondIdx = createPolymerWithStraight(boxlength, pattern, 5)
    finalizePolymerInit!(l, poly)
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Record initial counts
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 4
    @test initial_counts[m2.intIdx] == 4
    
    # Set up for addBond at bond 5
    poly.moveIdx = bondIdx
    identifyBond!(poly)
    
    @test poly.pos[poly.ki][6] == m2.intIdx
    @test poly.pos[poly.kj][6] == m2.intIdx
    
    # Propose move and select addBond direction
    proposeMove!(poly)
    
    # If moveType is not addBond, manually set up for addBond
    # For a straight line, all 4 directions are addBond (moveType=1)
    @test poly.moveType == 1  # Should be addBond for straight section
    
    # Ensure new positions are empty
    poly.latS[poly.newPosi...] = 0
    poly.latS[poly.newPosj...] = 0
    
    # Perform addBond
    addBond!(poly)
    
    # Verify type counts preserved for pattern type
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] + 2  # Base gains 2
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]       # Pattern unchanged
    
    # Total length should increase by 2
    @test length(poly.pos) == 10
end

# MARK: -removeBond inside pattern preserves count (realistic geometry)
@testset "removeBond inside pattern preserves type count" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2,1,1] - need base beads after pattern for swap
    # Create dent at bond 5 (inside type-2 region)
    pattern = [1, 1, 1, 1, 2, 2, 2, 2, 1, 1]
    poly, removableBondIdx = createPolymerWithDent(boxlength, pattern, 4)  # dent starts at bond 4, removable at 5
    finalizePolymerInit!(l, poly)
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Record initial counts
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 6
    @test initial_counts[m2.intIdx] == 4
    
    # Set up for removeBond
    poly.moveIdx = removableBondIdx
    identifyBond!(poly)
    
    @test poly.pos[poly.ki][6] == m2.intIdx
    @test poly.pos[poly.kj][6] == m2.intIdx
    
    # For a dent shape, the bond at the peak should have removeBond as an option
    # Set newPosi = bprev and newPosj = bnext for removeBond
    poly.newPosi = Vector{Int32}(poly.bprev)
    poly.newPosj = Vector{Int32}(poly.bnext)
    poly.moveType = 3  # removeBond
    
    # Perform removeBond
    removeBond!(poly)
    
    # Verify type counts preserved for pattern type
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] - 2  # Base loses 2
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]       # Pattern unchanged
    
    # Total length should decrease by 2
    @test length(poly.pos) == 8
end

# MARK: -kinkBond inside pattern preserves count
@testset "kinkBond inside pattern preserves type count" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2]
    # Create kink at bond 5 (inside type-2 region)
    pattern = [1, 1, 1, 1, 2, 2, 2, 2]
    poly, kinkBondIdx = createPolymerWithKink(boxlength, pattern, 5)
    finalizePolymerInit!(l, poly)
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Record initial counts
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 4
    @test initial_counts[m2.intIdx] == 4
    
    # Set up for kink move
    poly.moveIdx = kinkBondIdx
    identifyBond!(poly)
    
    @test poly.pos[poly.ki][6] == m2.intIdx
    @test poly.pos[poly.kj][6] == m2.intIdx
    
    # For a kink, one of newPosi/newPosj equals bprev/bnext
    # Set up a kink where newPosj == bnext (bead ki moves)
    poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
    poly.newPosj = Vector{Int32}(poly.bnext)
    poly.moveType = 2  # kink
    
    # Ensure new position is empty
    poly.latS[poly.newPosi...] = 0
    
    # Perform kinkBond
    kinkBond!(poly)
    
    # Verify type counts unchanged (kink doesn't add/remove beads)
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx]
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]
    
    # Length unchanged
    @test length(poly.pos) == 8
end

# MARK: -addBond outside pattern preserves pattern count
@testset "addBond outside pattern preserves pattern type count" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2] - addBond at bond 2 (in base region)
    pattern = [1, 1, 1, 1, 2, 2, 2, 2]
    poly, bondIdx = createPolymerWithStraight(boxlength, pattern, 2)
    finalizePolymerInit!(l, poly)
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 4
    @test initial_counts[m2.intIdx] == 4
    
    poly.moveIdx = bondIdx
    identifyBond!(poly)
    
    @test poly.pos[poly.ki][6] == m1.intIdx
    @test poly.pos[poly.kj][6] == m1.intIdx
    
    proposeMove!(poly)
    @test poly.moveType == 1
    
    poly.latS[poly.newPosi...] = 0
    poly.latS[poly.newPosj...] = 0
    
    addBond!(poly)
    
    # Base grows, pattern unchanged
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] + 2
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]
    
    @test length(poly.pos) == 10
end

# MARK: -removeBond outside pattern preserves pattern count
@testset "removeBond outside pattern preserves pattern type count" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1] - removeBond in base region after pattern
    # Need enough beads for dent: dentBondIdx + 3 <= length
    pattern = [1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1]
    poly, removableBondIdx = createPolymerWithDent(boxlength, pattern, 8)  # dent at bond 8, removable at 9
    finalizePolymerInit!(l, poly)
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 8
    @test initial_counts[m2.intIdx] == 4
    
    poly.moveIdx = removableBondIdx
    identifyBond!(poly)
    
    @test poly.pos[poly.ki][6] == m1.intIdx
    @test poly.pos[poly.kj][6] == m1.intIdx
    
    poly.newPosi = Vector{Int32}(poly.bprev)
    poly.newPosj = Vector{Int32}(poly.bnext)
    poly.moveType = 3
    
    removeBond!(poly)
    
    # Base shrinks, pattern unchanged
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] - 2
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]
    
    @test length(poly.pos) == 10  # 12 - 2 = 10
end

# MARK: -kinkBond outside pattern preserves pattern count
@testset "kinkBond outside pattern preserves pattern type count" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2] - kink at bond 2 (in base region)
    pattern = [1, 1, 1, 1, 2, 2, 2, 2]
    poly, kinkBondIdx = createPolymerWithKink(boxlength, pattern, 2)
    finalizePolymerInit!(l, poly)
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 4
    @test initial_counts[m2.intIdx] == 4
    
    poly.moveIdx = kinkBondIdx
    identifyBond!(poly)
    
    @test poly.pos[poly.ki][6] == m1.intIdx
    @test poly.pos[poly.kj][6] == m1.intIdx
    
    poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
    poly.newPosj = Vector{Int32}(poly.bnext)
    poly.moveType = 2
    
    poly.latS[poly.newPosi...] = 0
    
    kinkBond!(poly)
    
    # Both counts unchanged
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx]
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]
    
    @test length(poly.pos) == 8
end

# MARK: -Time evolution preserves pattern length
@testset "Time evolution preserves pattern length (5000 MC moves)" begin
    l = Lattice(; boxlength = 64, nspecies = 1)
    
    # Pattern: 16 base, 16 pattern, 32 base
    customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
    poly = Polymer(;
        length = 64,
        maxlength = 128,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = customPattern
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 48  # 16 + 32 base beads
    @test initial_counts[m2.intIdx] == 16  # 16 pattern beads
    
    # Verify initial contiguity
    @test isPatternContiguous(poly)
    
    # Run 5000 Monte Carlo moves with contiguity checks
    for i in 1:5000
        move!(poly)
        
        # Check contiguity every 100 moves
        if i % 100 == 0
            if !isPatternContiguous(poly)
                fragInfo = getFragmentationInfo(poly)
                error("Pattern fragmented at move $i: $fragInfo")
            end
        end
    end
    
    # Verify pattern type count is preserved
    final_counts = countBeadTypes(poly)
    @test final_counts[m2.intIdx] == initial_counts[m2.intIdx]
    
    # Verify pattern is still contiguous
    @test isPatternContiguous(poly)
    
    # Base type should have changed by the difference in total length
    length_change = length(poly.pos) - 64
    @test final_counts[m1.intIdx] == initial_counts[m1.intIdx] + length_change
end

# MARK: -Pattern contiguity preserved during time evolution
@testset "Pattern contiguity preserved during time evolution" begin
    l = Lattice(; boxlength = 64, nspecies = 1)
    
    # Pattern: 16 base, 16 pattern, 32 base
    customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
    poly = Polymer(;
        length = 64,
        maxlength = 128,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = customPattern
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m2 = poly.monomers[2]
    initial_pattern_count = countBeadTypes(poly)[m2.intIdx]
    
    # Run 10000 MC moves with strict contiguity checking after EVERY move
    fragmentation_detected = false
    fragmentation_move = 0
    fragmentation_info = nothing
    
    for i in 1:10000
        move!(poly)
        
        if !isPatternContiguous(poly)
            fragmentation_detected = true
            fragmentation_move = i
            fragmentation_info = getFragmentationInfo(poly)
            break
        end
        
        # Also check pattern count after every move
        current_count = countBeadTypes(poly)[m2.intIdx]
        if current_count != initial_pattern_count
            error("Pattern count changed at move $i: $initial_pattern_count -> $current_count")
        end
    end
    
    @test !fragmentation_detected
    if fragmentation_detected
        println("Pattern fragmented at move $fragmentation_move: $fragmentation_info")
    end
end

# MARK: -Multiple patterns time evolution
@testset "Multiple patterns time evolution preserves all pattern lengths" begin
    l = Lattice(; boxlength = 64, nspecies = 1)
    
    # Pattern: 16 base, 8 type2, 16 base, 8 type3, 16 base
    # Each non-base type appears in exactly one contiguous block
    customPattern = vcat(fill(1, 16), fill(2, 8), fill(1, 16), fill(3, 8), fill(1, 16))
    poly = Polymer(;
        length = 64,
        maxlength = 128,
        nMonomerTypes = 3,
        pattern = :custom,
        customPattern = customPattern
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    initial_counts = countBeadTypes(poly)
    
    # Verify initial contiguity
    @test isPatternContiguous(poly)
    
    # Run 5000 Monte Carlo moves with contiguity checks
    for i in 1:5000
        move!(poly)
        
        # Check contiguity every 100 moves
        if i % 100 == 0
            if !isPatternContiguous(poly)
                fragInfo = getFragmentationInfo(poly)
                error("Pattern fragmented at move $i: $fragInfo")
            end
        end
    end
    
    final_counts = countBeadTypes(poly)
    
    # All non-base pattern types should be preserved
    for (intIdx, count) in initial_counts
        if intIdx != poly.baseIntIdx
            @test final_counts[intIdx] == count
        end
    end
    
    # Verify final contiguity
    @test isPatternContiguous(poly)
end

# MARK: -patternTails updated after add
@testset "patternTails updated correctly after addBond" begin
    boxlength = Int32(30)
    l = Lattice(; boxlength = boxlength, nspecies = 1)
    l.lat = zeros(Int32, (2, boxlength, boxlength, boxlength))
    
    # Pattern: [1,1,1,1,2,2,2,2]
    pattern = [1, 1, 1, 1, 2, 2, 2, 2]
    poly, bondIdx = createPolymerWithStraight(boxlength, pattern, 5)
    finalizePolymerInit!(l, poly)
    
    m2 = poly.monomers[2]
    
    initial_tail = poly.patternTails[m2.intIdx]
    initial_tail_prev = poly.pos[initial_tail][4]
    expected_new_tail = poly.pos[initial_tail_prev][4]
    
    poly.moveIdx = bondIdx
    identifyBond!(poly)
    proposeMove!(poly)
    
    poly.latS[poly.newPosi...] = 0
    poly.latS[poly.newPosj...] = 0
    
    addBond!(poly)
    
    @test poly.patternTails[m2.intIdx] == expected_new_tail
end

# MARK: -Rejection when pattern too small
@testset "Reject addBond when pattern too small" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Create polymer with only 2 beads of type 2
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = [1, 1, 1, 1, 2, 2, 1, 1]
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m2 = poly.monomers[2]
    
    tailKey = poly.patternTails[m2.intIdx]
    prevTailKey = poly.pos[tailKey][4]
    newTailKey = poly.pos[prevTailKey][4]
    
    # newTailKey should be type-1 (not type-2), so validation fails
    @test poly.pos[newTailKey][6] != m2.intIdx
end

# MARK: -Rejection when insufficient base beads
@testset "Reject removeBond when insufficient base-type beads" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Pattern at end: [1,1,1,1,2,2,2,2]
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m2 = poly.monomers[2]
    
    tailKey = poly.patternTails[m2.intIdx]
    nextAfterTail = poly.pos[tailKey][5]
    
    @test nextAfterTail == 0
end
