# Implementation Plan: Pattern Pinning to Polymer Head

## Overview

This plan implements a mechanism to keep patterns pinned at fixed positions relative to the polymer head. Instead of rejecting moves that would cause drift, we compensate by shifting downstream patterns to maintain their positions.

**Key invariant**: Each pattern's head position (distance from polymer head in chain positions) remains constant throughout the simulation.

**Consequence**: Only the tail region of the polymer (after the last pattern) can grow or shrink.

---

## Mechanism Summary

When an add/remove occurs at position X:
- **Add**: All patterns with head position > X shift upstream by 2 (toward polymer head)
- **Remove**: All patterns with head position > X shift downstream by 2 (toward polymer tail)

Shifting a pattern upstream by 2:
- Convert 2 base beads immediately before pattern head → pattern type
- Convert 2 pattern beads at pattern tail → base type

Shifting a pattern downstream by 2:
- Convert 2 base beads immediately after pattern tail → pattern type  
- Convert 2 pattern beads at pattern head → base type

---

## Phase 1: Add Pattern Head Tracking

### 1.1 Add `patternHeads` and `patternHeadPos` fields to Polymer struct

**File**: `src/types.jl`

Add to Polymer struct (near existing `patternTails` and `patternTailPos`):

```julia
patternHeads::Union{Nothing, Dict{Int32, Int32}} = nothing      # intIdx -> head key
patternHeadPos::Union{Nothing, Dict{Int32, Int32}} = nothing    # intIdx -> head chain position (1-indexed)
```

### 1.2 Initialize pattern heads in `_initPatternTails!`

**File**: `src/init.jl`

Rename function to `_initPatternBoundaries!` and extend to track heads:

```julia
# Initialize patternTails, patternHeads, and their positions
function _initPatternBoundaries!(p::Polymer)
    intIdxList = sort(collect(keys(p.monomers)))
    p.baseIntIdx = intIdxList[1]  # base type = first intIdx (smallest)
    
    # For uniform pattern, no boundaries to track
    if p.nMonomerTypes == 1 || p.pattern == :uniform
        p.patternTails = nothing
        p.patternTailPos = nothing
        p.originalPatternTailPos = nothing
        p.patternHeads = nothing
        p.patternHeadPos = nothing
        return
    end
    
    p.patternTails = Dict{Int32, Int32}()
    p.patternTailPos = Dict{Int32, Int32}()
    p.patternHeads = Dict{Int32, Int32}()
    p.patternHeadPos = Dict{Int32, Int32}()
    
    key = p.head
    chainPos = Int32(1)
    prevIntIdx = Int32(0)
    
    while key != 0
        currentIntIdx = p.pos[key][6]
        nextKey = p.pos[key][5]
        
        if currentIntIdx != p.baseIntIdx
            # Track head: first bead of each pattern type (transition from different type)
            if prevIntIdx != currentIntIdx
                p.patternHeads[currentIntIdx] = key
                p.patternHeadPos[currentIntIdx] = chainPos
            end
            
            # Track tail: last bead of each pattern type (before type changes)
            if nextKey == 0 || p.pos[nextKey][6] != currentIntIdx
                p.patternTails[currentIntIdx] = key
                p.patternTailPos[currentIntIdx] = chainPos
            end
        end
        
        prevIntIdx = currentIntIdx
        key = nextKey
        chainPos += 1
    end
    
    # Store original positions (immutable reference)
    p.originalPatternTailPos = copy(p.patternTailPos)
end
```

### 1.3 Update call sites

**File**: `src/init.jl`

Replace call to `_initPatternTails!` with `_initPatternBoundaries!` in `_assignPattern!`.

### 1.4 Tests for Phase 1

```julia
@testset "Pattern head tracking initialization" begin
    
    @testset "Single pattern - heads and tails initialized" begin
        l = Lattice(boxlength=32, nspecies=2)
        # Pattern: 8 base, 8 pattern, 16 base
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 16))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        @test p.patternHeads !== nothing
        @test p.patternHeadPos !== nothing
        @test haskey(p.patternHeads, patternIntIdx)
        @test haskey(p.patternHeadPos, patternIntIdx)
        @test p.patternHeadPos[patternIntIdx] == 9   # head at position 9
        @test p.patternTailPos[patternIntIdx] == 16  # tail at position 16
    end
    
    @testset "Multiple patterns - independent tracking" begin
        l = Lattice(boxlength=64, nspecies=3)
        # 8 base, 8 type2, 8 base, 8 type3, 32 base
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 8), fill(3, 8), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=3, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        type2IntIdx = intIdxList[2]
        type3IntIdx = intIdxList[3]
        
        @test p.patternHeadPos[type2IntIdx] == 9
        @test p.patternTailPos[type2IntIdx] == 16
        @test p.patternHeadPos[type3IntIdx] == 25
        @test p.patternTailPos[type3IntIdx] == 32
    end
    
    @testset "Uniform polymer - no tracking" begin
        l = Lattice(boxlength=32, nspecies=1)
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=1, pattern=:uniform)
        initMonomers!(p)
        init!(l, [p])
        
        @test p.patternHeads === nothing
        @test p.patternHeadPos === nothing
    end
    
end
```

---

## Phase 2: Implement Pattern Shifting Functions

### 2.1 Add `shiftPatternUpstream!` function

**File**: `src/move.jl`

Shifts a pattern 2 positions toward polymer head by swapping types:

```julia
"""
    shiftPatternUpstream!(p::Polymer, intIdx::Int32) -> (swappedKeys::Vector{Int32}, oldTypes::Vector{Int32}, newTypes::Vector{Int32})

Shifts pattern `intIdx` upstream (toward head) by 2 positions.
- Converts 2 base beads before pattern head → pattern type
- Converts 2 pattern beads at pattern tail → base type

Returns the keys of swapped beads and their old/new types for energy calculation.
"""
function shiftPatternUpstream!(p::Polymer, intIdx::Int32)
    headKey = p.patternHeads[intIdx]
    tailKey = p.patternTails[intIdx]
    
    # Find 2 beads before current head (these are base type, will become pattern)
    prevKey1 = p.pos[headKey][4]  # bead at head-1
    prevKey2 = p.pos[prevKey1][4] # bead at head-2
    
    # Find tail and tail-1 (these are pattern type, will become base)
    tailPrevKey = p.pos[tailKey][4]  # bead at tail-1
    
    # Record for energy calculation
    swappedKeys = Int32[prevKey2, prevKey1, tailPrevKey, tailKey]
    oldTypes = Int32[p.baseIntIdx, p.baseIntIdx, intIdx, intIdx]
    newTypes = Int32[intIdx, intIdx, p.baseIntIdx, p.baseIntIdx]
    
    # Swap types: base beads before head → pattern
    p.pos[prevKey1][6] = intIdx
    p.pos[prevKey2][6] = intIdx
    p.latS[p.pos[prevKey1][1], p.pos[prevKey1][2], p.pos[prevKey1][3]] = intIdx
    p.latS[p.pos[prevKey2][1], p.pos[prevKey2][2], p.pos[prevKey2][3]] = intIdx
    
    # Swap types: tail beads → base
    p.pos[tailKey][6] = p.baseIntIdx
    p.pos[tailPrevKey][6] = p.baseIntIdx
    p.latS[p.pos[tailKey][1], p.pos[tailKey][2], p.pos[tailKey][3]] = p.baseIntIdx
    p.latS[p.pos[tailPrevKey][1], p.pos[tailPrevKey][2], p.pos[tailPrevKey][3]] = p.baseIntIdx
    
    # Update pattern boundaries
    # New head is prevKey2 (2 positions earlier)
    p.patternHeads[intIdx] = prevKey2
    # Head position unchanged (pattern pinned to head)
    
    # New tail is tailPrevKey's previous bead
    newTailKey = p.pos[tailPrevKey][4]
    p.patternTails[intIdx] = newTailKey
    # Tail position decreases by 2
    p.patternTailPos[intIdx] -= 2
    
    return swappedKeys, oldTypes, newTypes
end
```

### 2.2 Add `shiftPatternDownstream!` function

**File**: `src/move.jl`

Shifts a pattern 2 positions toward polymer tail:

```julia
"""
    shiftPatternDownstream!(p::Polymer, intIdx::Int32) -> (swappedKeys::Vector{Int32}, oldTypes::Vector{Int32}, newTypes::Vector{Int32})

Shifts pattern `intIdx` downstream (toward tail) by 2 positions.
- Converts 2 base beads after pattern tail → pattern type
- Converts 2 pattern beads at pattern head → base type

Returns the keys of swapped beads and their old/new types for energy calculation.
"""
function shiftPatternDownstream!(p::Polymer, intIdx::Int32)
    headKey = p.patternHeads[intIdx]
    tailKey = p.patternTails[intIdx]
    
    # Find 2 beads after current tail (these are base type, will become pattern)
    nextKey1 = p.pos[tailKey][5]  # bead at tail+1
    nextKey2 = p.pos[nextKey1][5] # bead at tail+2
    
    # Find head and head+1 (these are pattern type, will become base)
    headNextKey = p.pos[headKey][5]  # bead at head+1
    
    # Record for energy calculation
    swappedKeys = Int32[headKey, headNextKey, nextKey1, nextKey2]
    oldTypes = Int32[intIdx, intIdx, p.baseIntIdx, p.baseIntIdx]
    newTypes = Int32[p.baseIntIdx, p.baseIntIdx, intIdx, intIdx]
    
    # Swap types: head beads → base
    p.pos[headKey][6] = p.baseIntIdx
    p.pos[headNextKey][6] = p.baseIntIdx
    p.latS[p.pos[headKey][1], p.pos[headKey][2], p.pos[headKey][3]] = p.baseIntIdx
    p.latS[p.pos[headNextKey][1], p.pos[headNextKey][2], p.pos[headNextKey][3]] = p.baseIntIdx
    
    # Swap types: base beads after tail → pattern
    p.pos[nextKey1][6] = intIdx
    p.pos[nextKey2][6] = intIdx
    p.latS[p.pos[nextKey1][1], p.pos[nextKey1][2], p.pos[nextKey1][3]] = intIdx
    p.latS[p.pos[nextKey2][1], p.pos[nextKey2][2], p.pos[nextKey2][3]] = intIdx
    
    # Update pattern boundaries
    # New head is headNextKey's next bead
    newHeadKey = p.pos[headNextKey][5]
    p.patternHeads[intIdx] = newHeadKey
    # Head position increases by 2
    p.patternHeadPos[intIdx] += 2
    
    # New tail is nextKey2
    p.patternTails[intIdx] = nextKey2
    # Tail position increases by 2
    p.patternTailPos[intIdx] += 2
    
    return swappedKeys, oldTypes, newTypes
end
```

### 2.3 Tests for Phase 2

```julia
@testset "Pattern shifting functions" begin
    
    @testset "shiftPatternUpstream! - single pattern" begin
        l = Lattice(boxlength=32, nspecies=2)
        # 8 base, 8 pattern, 16 base
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 16))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        # Initial state
        @test p.patternHeadPos[patternIntIdx] == 9
        @test p.patternTailPos[patternIntIdx] == 16
        
        # Shift upstream
        swappedKeys, oldTypes, newTypes = Droplet.shiftPatternUpstream!(p, patternIntIdx)
        
        # Verify 4 beads swapped
        @test length(swappedKeys) == 4
        
        # Verify positions updated
        @test p.patternHeadPos[patternIntIdx] == 9  # unchanged (pinned)
        @test p.patternTailPos[patternIntIdx] == 14 # moved upstream by 2
        
        # Verify actual chain types changed
        # Count pattern beads - should still be 8
        patternCount = 0
        key = p.head
        while key != 0
            if p.pos[key][6] == patternIntIdx
                patternCount += 1
            end
            key = p.pos[key][5]
        end
        @test patternCount == 8
    end
    
    @testset "shiftPatternDownstream! - single pattern" begin
        l = Lattice(boxlength=32, nspecies=2)
        # 8 base, 8 pattern, 16 base
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 16))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        # Shift downstream
        swappedKeys, oldTypes, newTypes = Droplet.shiftPatternDownstream!(p, patternIntIdx)
        
        # Verify positions updated
        @test p.patternHeadPos[patternIntIdx] == 11 # moved downstream by 2
        @test p.patternTailPos[patternIntIdx] == 18 # moved downstream by 2
        
        # Pattern length preserved
        patternCount = 0
        key = p.head
        while key != 0
            if p.pos[key][6] == patternIntIdx
                patternCount += 1
            end
            key = p.pos[key][5]
        end
        @test patternCount == 8
    end
    
    @testset "Shift preserves lattice consistency" begin
        l = Lattice(boxlength=32, nspecies=2)
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 16))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        # Shift and verify lattice matches pos
        Droplet.shiftPatternUpstream!(p, patternIntIdx)
        
        key = p.head
        while key != 0
            posIntIdx = p.pos[key][6]
            x, y, z = p.pos[key][1], p.pos[key][2], p.pos[key][3]
            @test p.latS[x, y, z] == posIntIdx
            key = p.pos[key][5]
        end
    end
    
end
```

---

## Phase 3: Integrate Shifting into addBond! and removeBond!

### 3.1 Helper function to get patterns downstream of position

**File**: `src/move.jl`

```julia
"""
    getPatternsDownstreamOf(p::Polymer, chainPos::Int32) -> Vector{Int32}

Returns intIdx values of all patterns whose head position is > chainPos,
sorted by head position (closest first).
"""
function getPatternsDownstreamOf(p::Polymer, chainPos::Int32)::Vector{Int32}
    if p.patternHeadPos === nothing
        return Int32[]
    end
    
    downstream = Int32[]
    for (intIdx, headPos) in p.patternHeadPos
        if headPos > chainPos
            push!(downstream, intIdx)
        end
    end
    
    # Sort by head position (closest patterns first)
    sort!(downstream, by = idx -> p.patternHeadPos[idx])
    return downstream
end
```

### 3.2 Modify `addBond!` to shift downstream patterns

**File**: `src/move.jl`

After inserting new beads and performing existing swap (if inside pattern), shift all downstream patterns:

```julia
function addBond!(p::Polymer)
    # Get insertion position BEFORE modifying chain
    insertPos = p.moveIdx
    
    # Get current intIdx of the bond being modified
    intIdx = p.pos[p.ki][6]
    
    # Capture pattern tail info BEFORE modifying chain (existing logic for pattern-internal adds)
    tailKey = nothing
    prevTailKey = nothing
    newTailKey = nothing
    if p.patternTails !== nothing && intIdx != p.baseIntIdx && haskey(p.patternTails, intIdx)
        tailKey = p.patternTails[intIdx]
        prevTailKey = p.pos[tailKey][4]
        newTailKey = p.pos[prevTailKey][4]
    end
    
    # === EXISTING BOND INSERTION LOGIC ===
    # ... (insert 2 new beads, update chain pointers, etc.)
    # ... (this part remains unchanged)
    
    # === EXISTING PATTERN-INTERNAL SWAP (if adding inside pattern) ===
    if tailKey !== nothing
        # Swap tail and prev-tail to base type (existing logic)
        p.pos[tailKey][6] = p.baseIntIdx
        p.pos[prevTailKey][6] = p.baseIntIdx
        p.latS[p.pos[tailKey][1], p.pos[tailKey][2], p.pos[tailKey][3]] = p.baseIntIdx
        p.latS[p.pos[prevTailKey][1], p.pos[prevTailKey][2], p.pos[prevTailKey][3]] = p.baseIntIdx
        
        # Update patternTails (existing logic)
        if newTailKey != 0 && p.pos[newTailKey][6] == intIdx
            p.patternTails[intIdx] = newTailKey
            p.patternTailPos[intIdx] -= 2
        else
            delete!(p.patternTails, intIdx)
            delete!(p.patternTailPos, intIdx)
            delete!(p.patternHeads, intIdx)
            delete!(p.patternHeadPos, intIdx)
        end
    end
    
    # === NEW: SHIFT ALL DOWNSTREAM PATTERNS UPSTREAM ===
    if p.patternHeadPos !== nothing
        downstreamPatterns = getPatternsDownstreamOf(p, Int32(insertPos))
        for pIntIdx in downstreamPatterns
            shiftPatternUpstream!(p, pIntIdx)
        end
    end
    
    return p
end
```

### 3.3 Modify `removeBond!` to shift downstream patterns

**File**: `src/move.jl`

After removing beads and performing existing swap (if inside pattern), shift all downstream patterns:

```julia
function removeBond!(p::Polymer)
    # Get removal position BEFORE modifying chain
    removePos = p.moveIdx
    
    # Get the intIdx before removing
    intIdx = p.pos[p.ki][6]
    
    # === EXISTING PATTERN-INTERNAL SWAP (if removing inside pattern) ===
    if p.patternTails !== nothing && intIdx != p.baseIntIdx && haskey(p.patternTails, intIdx) && isInPatternRegion(p, p.ki, intIdx)
        tailKey = p.patternTails[intIdx]
        nextAfterTail = p.pos[tailKey][5]
        nextNextAfterTail = nextAfterTail == 0 ? Int32(0) : p.pos[nextAfterTail][5]
        
        # Swap two base-type beads to pattern type
        p.pos[nextAfterTail][6] = intIdx
        p.pos[nextNextAfterTail][6] = intIdx
        p.latS[p.pos[nextAfterTail][1], p.pos[nextAfterTail][2], p.pos[nextAfterTail][3]] = intIdx
        p.latS[p.pos[nextNextAfterTail][1], p.pos[nextNextAfterTail][2], p.pos[nextNextAfterTail][3]] = intIdx
        
        # Update patternTails
        p.patternTails[intIdx] = nextNextAfterTail
        p.patternTailPos[intIdx] += 2
    end
    
    # === EXISTING BEAD REMOVAL LOGIC ===
    # ... (remove 2 beads, update chain pointers, etc.)
    # ... (this part remains unchanged)
    
    # === NEW: SHIFT ALL DOWNSTREAM PATTERNS DOWNSTREAM ===
    if p.patternHeadPos !== nothing
        downstreamPatterns = getPatternsDownstreamOf(p, Int32(removePos))
        for pIntIdx in downstreamPatterns
            shiftPatternDownstream!(p, pIntIdx)
        end
    end
    
    return p
end
```

### 3.4 Tests for Phase 3

```julia
@testset "Pattern shifting on add/remove" begin
    
    @testset "Add in base region before pattern - pattern shifts upstream" begin
        l = Lattice(boxlength=64, nspecies=2)
        # 16 base, 8 pattern, 40 base
        customPattern = vcat(fill(1, 16), fill(2, 8), fill(1, 40))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        initialHeadPos = p.patternHeadPos[patternIntIdx]
        @test initialHeadPos == 17
        
        # Run moves until an add happens in base region before pattern
        for _ in 1:10000
            lengthBefore = p.length
            Droplet.move!(p)
            
            # Pattern head position should remain constant
            if haskey(p.patternHeadPos, patternIntIdx)
                @test p.patternHeadPos[patternIntIdx] == initialHeadPos
            end
        end
    end
    
    @testset "Add inside pattern - downstream pattern shifts upstream" begin
        l = Lattice(boxlength=64, nspecies=3)
        # 8 base, 8 type2, 16 base, 8 type3, 24 base
        customPattern = vcat(fill(1, 8), fill(2, 8), fill(1, 16), fill(3, 8), fill(1, 24))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=3, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        type2IntIdx = intIdxList[2]
        type3IntIdx = intIdxList[3]
        
        initialType2HeadPos = p.patternHeadPos[type2IntIdx]
        initialType3HeadPos = p.patternHeadPos[type3IntIdx]
        
        @test initialType2HeadPos == 9
        @test initialType3HeadPos == 33
        
        # Run moves - both patterns should maintain their head positions
        for _ in 1:10000
            Droplet.move!(p)
            
            if haskey(p.patternHeadPos, type2IntIdx)
                @test p.patternHeadPos[type2IntIdx] == initialType2HeadPos
            end
            if haskey(p.patternHeadPos, type3IntIdx)
                @test p.patternHeadPos[type3IntIdx] == initialType3HeadPos
            end
        end
    end
    
    @testset "Remove in base region - pattern shifts downstream" begin
        l = Lattice(boxlength=64, nspecies=2)
        # 24 base, 8 pattern, 32 base (plenty of room)
        customPattern = vcat(fill(1, 24), fill(2, 8), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        initialHeadPos = p.patternHeadPos[patternIntIdx]
        @test initialHeadPos == 25
        
        # Run moves - pattern head position should remain constant
        for _ in 1:10000
            Droplet.move!(p)
            
            if haskey(p.patternHeadPos, patternIntIdx)
                @test p.patternHeadPos[patternIntIdx] == initialHeadPos
            end
        end
    end
    
    @testset "Pattern length preserved through all operations" begin
        l = Lattice(boxlength=64, nspecies=2)
        customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        initialPatternLength = 16
        
        for _ in 1:10000
            Droplet.move!(p)
            
            # Count pattern beads
            patternCount = 0
            key = p.head
            while key != 0
                if p.pos[key][6] == patternIntIdx
                    patternCount += 1
                end
                key = p.pos[key][5]
            end
            
            @test patternCount == initialPatternLength
        end
    end
    
end
```

---

## Phase 4: Remove Drift Rejection Logic

### 4.1 Remove drift limit checks from `move!`

**File**: `src/move.jl`

Remove the drift limit computation and rejection logic from `move!`. The relevant block to remove:

```julia
# REMOVE THIS BLOCK:
# Compute drift limit inline
driftLimit = p.length < 8 ? Int32(-1) : p.length < 32 ? Int32(2) : p.length < 64 ? Int32(4) : Int32(8)

# Check drift: addBond moves tail backward by 2
if driftLimit > 0
    newTailPos = p.patternTailPos[currentIntIdx] - 2
    if abs(newTailPos - p.originalPatternTailPos[currentIntIdx]) > driftLimit
        return false
    end
end
```

### 4.2 Remove `originalPatternTailPos` field

**File**: `src/types.jl`

Remove the field (no longer needed since we don't compare against original):

```julia
# REMOVE:
originalPatternTailPos::Union{Nothing, Dict{Int32, Int32}} = nothing
```

### 4.3 Update `_initPatternBoundaries!` to not set `originalPatternTailPos`

Remove the line that copies to `originalPatternTailPos`.

### 4.4 Tests for Phase 4

```julia
@testset "No rejection - all valid moves accepted" begin
    l = Lattice(boxlength=64, nspecies=2)
    customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
    p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
    initMonomers!(p)
    init!(l, [p])
    
    # Run many moves - track acceptance rate
    acceptCount = 0
    totalMoves = 10000
    
    for _ in 1:totalMoves
        if Droplet.move!(p)
            acceptCount += 1
        end
    end
    
    # Acceptance rate should be reasonable (not near-zero due to rejections)
    acceptanceRate = acceptCount / totalMoves
    @test acceptanceRate > 0.1  # At least 10% acceptance
end
```

---

## Phase 5: Update Energy Calculations for Pattern Shifts

### 5.1 Extend `deltaHam!` to account for shift swaps

**File**: `src/hamiltonian.jl`

When patterns are shifted, 4 beads per pattern change type. We need to account for:
1. Energy change from the swapped beads' neighbors
2. Chemical potential change if types have different chemPot

**Approach**: Collect all swapped beads during `addBond!`/`removeBond!`, then compute energy delta.

First, modify `addBond!` and `removeBond!` to return swap information:

```julia
# In addBond! - collect all swaps
function addBond!(p::Polymer)
    # ... existing logic ...
    
    allSwappedKeys = Int32[]
    allOldTypes = Int32[]
    allNewTypes = Int32[]
    
    # Existing pattern-internal swap
    if tailKey !== nothing
        append!(allSwappedKeys, [prevTailKey, tailKey])
        append!(allOldTypes, [intIdx, intIdx])
        append!(allNewTypes, [p.baseIntIdx, p.baseIntIdx])
        # ... do the swap ...
    end
    
    # Downstream pattern shifts
    if p.patternHeadPos !== nothing
        downstreamPatterns = getPatternsDownstreamOf(p, Int32(insertPos))
        for pIntIdx in downstreamPatterns
            keys, oldT, newT = shiftPatternUpstream!(p, pIntIdx)
            append!(allSwappedKeys, keys)
            append!(allOldTypes, oldT)
            append!(allNewTypes, newT)
        end
    end
    
    # Store for deltaHam calculation
    p.swappedKeys = allSwappedKeys
    p.swappedOldTypes = allOldTypes
    p.swappedNewTypes = allNewTypes
    
    return p
end
```

### 5.2 Add swap tracking fields to Polymer struct

**File**: `src/types.jl`

```julia
# Temporary storage for swap info (used by deltaHam)
swappedKeys::Union{Nothing, Vector{Int32}} = nothing
swappedOldTypes::Union{Nothing, Vector{Int32}} = nothing
swappedNewTypes::Union{Nothing, Vector{Int32}} = nothing
```

### 5.3 Compute swap energy contribution in `deltaHam!`

**File**: `src/hamiltonian.jl`

Add after existing deltaHam calculation:

```julia
function deltaHam!(p::Polymer)
    # ... existing deltaHam logic for add/kink/remove ...
    
    # Add energy contribution from type swaps
    if p.swappedKeys !== nothing && !isempty(p.swappedKeys)
        for i in 1:length(p.swappedKeys)
            key = p.swappedKeys[i]
            oldIntIdx = p.swappedOldTypes[i]
            newIntIdx = p.swappedNewTypes[i]
            
            if oldIntIdx == newIntIdx
                continue
            end
            
            pos = view(p.pos[key], 1:3)
            mOld = p.monomers[oldIntIdx]
            mNew = p.monomers[newIntIdx]
            
            # Interaction energy change with neighbors
            # Need to count neighbors of each type and compute energy difference
            for (neighborIntIdx, neighborMon) in p.monomers
                neighborCount = nnp(mOld, neighborMon, pos)
                if neighborCount > 0
                    # Remove old interaction, add new interaction
                    oldEnergy = mOld.intMat[mOld.intIdx, neighborIntIdx] * neighborCount
                    newEnergy = mNew.intMat[mNew.intIdx, neighborIntIdx] * neighborCount
                    p.deltaH += (newEnergy - oldEnergy)
                end
            end
            
            # Chemical potential change
            p.deltaH += mNew.chemPot - mOld.chemPot
        end
    end
    
    return p
end
```

### 5.4 Tests for Phase 5

```julia
@testset "Energy consistency with pattern shifts" begin
    
    @testset "deltaHam matches full Hamiltonian recalculation" begin
        l = Lattice(boxlength=64, nspecies=2)
        customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p; defaultIntPot=-1.0)
        init!(l, [p])
        
        for _ in 1:1000
            Droplet.ham!(p)
            hamBefore = p.ham
            
            p.moveIdx = rand(1:p.length-1)
            Droplet.identifyBond!(p)
            
            if p.kj == 0
                continue
            end
            
            Droplet.proposeMove!(p)
            Droplet.deltaHam!(p)
            
            if Droplet.move!(p)
                Droplet.ham!(p)
                @test isapprox(hamBefore + p.deltaH, p.ham, atol=1e-10)
            end
        end
    end
    
    @testset "Energy with different chemPot per type" begin
        l = Lattice(boxlength=64, nspecies=2)
        customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        
        # Set different chemical potentials
        initMonomers!(p; defaultIntPot=-1.0)
        intIdxList = sort(collect(keys(p.monomers)))
        p.monomers[intIdxList[1]].chemPot = 0.0
        p.monomers[intIdxList[2]].chemPot = -0.5
        
        init!(l, [p])
        
        for _ in 1:500
            Droplet.ham!(p)
            hamBefore = p.ham
            
            p.moveIdx = rand(1:p.length-1)
            Droplet.identifyBond!(p)
            
            if p.kj == 0
                continue
            end
            
            Droplet.proposeMove!(p)
            Droplet.deltaHam!(p)
            
            if Droplet.move!(p)
                Droplet.ham!(p)
                @test isapprox(hamBefore + p.deltaH, p.ham, atol=1e-10)
            end
        end
    end
    
end
```

---

## Phase 6: Validation Constraints

### 6.1 Reject moves that would violate constraints

**File**: `src/move.jl`

Before allowing a move, check that shifting is possible:

**For addBond** (patterns shift upstream):
- **No validation needed** - the first pattern doesn't require base beads upstream
- Normal swap logic applies: pattern beads become base beads at the tail side
- The shift converts 2 beads at pattern head side FROM base TO pattern type
- This works even when those 2 beads are inside another upstream pattern

**For removeBond** (patterns shift downstream):
- **All** downstream patterns need 2 base beads after their tail
- If any pattern lacks 2 base beads downstream, the remove is rejected
- This prevents patterns from colliding during downstream shifts

```julia
# No validation needed for addBond - upstream shifts always work
# function canShiftPatternsForAdd(p::Polymer, insertPos::Int32)::Bool
#     return true  # Always allowed
# end

function canShiftPatternsForRemove(p::Polymer, removePos::Int32)::Bool
    if p.patternHeadPos === nothing
        return true
    end
    
    downstreamPatterns = getPatternsDownstreamOf(p, removePos)
    
    for pIntIdx in downstreamPatterns
        tailKey = p.patternTails[pIntIdx]
        nextKey1 = p.pos[tailKey][5]
        nextKey2 = nextKey1 == 0 ? Int32(0) : p.pos[nextKey1][5]
        
        # Need 2 base beads after tail
        if nextKey1 == 0 || p.pos[nextKey1][6] != p.baseIntIdx ||
           nextKey2 == 0 || p.pos[nextKey2][6] != p.baseIntIdx
            return false
        end
    end
    
    return true
end
```

Add to `move!` before accepting the move:

```julia
# Validate shift is possible (only removeBond needs validation)
if p.moveType == 3  # removeBond
    if !canShiftPatternsForRemove(p, Int32(p.moveIdx))
        return false
    end
end
# Note: addBond does not need validation - upstream shifts always work
```

### 6.2 Tests for Phase 6

```julia
@testset "Constraint validation" begin
    
    @testset "Add allowed even when pattern at position 1" begin
        l = Lattice(boxlength=32, nspecies=2)
        # Pattern starts at position 1 (no base beads before it)
        customPattern = vcat(fill(2, 8), fill(1, 24))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        intIdxList = sort(collect(keys(p.monomers)))
        patternIntIdx = intIdxList[2]
        
        # Pattern head should start at position 1
        @test p.patternHeadPos[patternIntIdx] == 1
        
        # Adds should work - pattern will shift upstream (no-op at head) or 
        # the shift swaps work normally
        addCount = 0
        for _ in 1:1000
            lengthBefore = p.length
            Droplet.move!(p)
            if p.length > lengthBefore
                addCount += 1
            end
        end
        
        # Adds should succeed (not all rejected)
        @test addCount > 0
    end
    
    @testset "Remove rejected when insufficient base beads after pattern tail" begin
        l = Lattice(boxlength=32, nspecies=2)
        # Pattern near end with only 2 base beads after (minimum)
        customPattern = vcat(fill(1, 20), fill(2, 8), fill(1, 4))
        p = Polymer(id=1, length=32, maxlength=64, boxlength=32,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        # Removes before pattern should be limited
        # Each remove in base region shifts pattern downstream, consuming base beads
        removeCount = 0
        for _ in 1:1000
            lengthBefore = p.length
            Droplet.move!(p)
            if p.length < lengthBefore
                removeCount += 1
            end
        end
        
        # Eventually can't remove anymore (pattern would collide with tail end)
        @test removeCount < 500
    end
    
    @testset "Remove in tail region (after all patterns) always allowed" begin
        l = Lattice(boxlength=64, nspecies=2)
        # Pattern in middle, lots of base beads after
        customPattern = vcat(fill(1, 16), fill(2, 8), fill(1, 40))
        p = Polymer(id=1, length=64, maxlength=128, boxlength=64,
                    nMonomerTypes=2, pattern=:custom, customPattern=customPattern)
        initMonomers!(p)
        init!(l, [p])
        
        # Force removes to happen in tail region only
        # This tests that removes after all patterns don't trigger shift validation
        tailRemoves = 0
        for _ in 1:5000
            # Manually set moveIdx to be after pattern
            if p.length > 30
                p.moveIdx = rand(25:p.length-1)  # After pattern (which ends at ~24)
                Droplet.identifyBond!(p)
                
                if p.kj != 0  # Valid bond
                    lengthBefore = p.length
                    Droplet.proposeMove!(p)
                    Droplet.deltaHam!(p)
                    
                    if Droplet.move!(p)
                        if p.length < lengthBefore
                            tailRemoves += 1
                        end
                    end
                end
            end
        end
        
        # Tail removes should succeed without restriction
        @test tailRemoves > 100
    end
    
end
```

---

## Summary

| Phase | Description | Files Modified |
|-------|-------------|----------------|
| 1 | Add pattern head tracking | types.jl, init.jl |
| 2 | Implement shift functions | move.jl |
| 3 | Integrate shifting into add/remove | move.jl |
| 4 | Remove drift rejection logic | move.jl, types.jl |
| 5 | Update energy calculations | hamiltonian.jl, types.jl |
| 6 | Add validation constraints | move.jl |

### New Fields

| Field | Type | Description |
|-------|------|-------------|
| `patternHeads` | `Dict{Int32, Int32}` | intIdx → head key |
| `patternHeadPos` | `Dict{Int32, Int32}` | intIdx → head chain position |
| `swappedKeys` | `Vector{Int32}` | Temporary: keys of beads swapped in last move |
| `swappedOldTypes` | `Vector{Int32}` | Temporary: old types of swapped beads |
| `swappedNewTypes` | `Vector{Int32}` | Temporary: new types of swapped beads |

### Removed Fields

| Field | Reason |
|-------|--------|
| `originalPatternTailPos` | No longer needed (patterns pinned, not drift-limited) |

### New Functions

| Function | Description |
|----------|-------------|
| `shiftPatternUpstream!(p, intIdx)` | Shift pattern 2 positions toward head |
| `shiftPatternDownstream!(p, intIdx)` | Shift pattern 2 positions toward tail |
| `getPatternsDownstreamOf(p, pos)` | Get all patterns with head > pos |
| `canShiftPatternsForAdd(p, pos)` | Check if upstream shift is possible |
| `canShiftPatternsForRemove(p, pos)` | Check if downstream shift is possible |

### Key Invariant

**Pattern head positions remain constant** relative to polymer head throughout the simulation. All length changes are absorbed by the tail region (after the last pattern).
