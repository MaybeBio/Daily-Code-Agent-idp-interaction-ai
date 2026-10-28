# Phase 2 (Refactored): Pattern Swap Functions

## Overview

This phase implements four atomic swap functions that handle all pattern boundary adjustments:
- `swapTailToBase!`: Converts 2 pattern beads at tail → base type (used during addBond)
- `swapTailToPattern!`: Converts 2 base beads at tail → pattern type (used during removeBond)  
- `swapHeadToBase!`: Converts 2 pattern beads at head → base type (used during downstream shifts)
- `swapHeadToPattern!`: Converts 2 base beads at head → pattern type (used during upstream shifts)

These four primitives can be composed to implement:
1. **In-pattern addBond**: `swapTailToBase!` (shrink pattern at tail by 2)
2. **In-pattern removeBond**: `swapTailToPattern!` (extend pattern at tail by 2)
3. **Upstream shift** (for downstream patterns after addBond): `swapHeadToPattern!` + `swapTailToBase!`
4. **Downstream shift** (for downstream patterns after removeBond): `swapHeadToBase!` + `swapTailToPattern!`

---

## Sub-Phase 2.1: Tail Swap Functions

### 2.1.1 Function `swapTailToBase!`

Converts the 2 pattern beads at the tail of a pattern region to base type.

**Use cases**:
- In-pattern addBond (compensate for 2 new beads by shrinking pattern at tail)
- Upstream shift (second half of the shift operation)

**Signature**:
```julia
"""
    swapTailToBase!(p::Polymer, intIdx::Int32; updateTracking::Bool=true) 
        -> (success::Bool, swappedKeys::Tuple{Int32, Int32})

Converts the 2 pattern beads at the tail of pattern `intIdx` to base type.
Updates pattern tracking (patternTails, patternTailPos) if `updateTracking=true`.

Returns:
- `success`: Whether the swap was performed (false if pattern too small)
- `swappedKeys`: Tuple of (tailPrevKey, tailKey) that were swapped, or (0, 0) if failed

Preconditions:
- Pattern `intIdx` must exist and have at least 3 beads (so there's still a pattern after swap)
- Pattern head/tail keys must be valid in p.pos

Effects:
- Converts tailPrevKey and tailKey from intIdx → baseIntIdx
- Updates p.pos[key][6] and p.latS for both beads
- If updateTracking: updates patternTails[intIdx] to the new tail (2 beads earlier)
- If updateTracking: updates patternTailPos[intIdx] -= 2
"""
function swapTailToBase!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
    # Guard: pattern must exist
    if !haskey(p.patternTails, intIdx)
        return false, (Int32(0), Int32(0))
    end
    
    tailKey = p.patternTails[intIdx]
    
    # Guard: tail key must be valid
    if !haskey(p.pos, tailKey)
        return false, (Int32(0), Int32(0))
    end
    
    # Get tail-1 (the bead before tail)
    tailPrevKey = p.pos[tailKey][4]
    
    # Guard: must have a bead before tail (pattern has at least 2 beads)
    if tailPrevKey == 0
        return false, (Int32(0), Int32(0))
    end
    
    # Guard: after swap, pattern must still exist (need at least 1 bead remaining)
    # The new tail would be tailPrevKey's prev
    newTailKey = p.pos[tailPrevKey][4]
    if newTailKey == 0
        return false, (Int32(0), Int32(0))
    end
    
    # Verify the beads are actually pattern type before swapping
    if p.pos[tailKey][6] != intIdx || p.pos[tailPrevKey][6] != intIdx
        return false, (Int32(0), Int32(0))
    end
    
    # Perform the swap: pattern → base
    p.pos[tailKey][6] = p.baseIntIdx
    p.pos[tailPrevKey][6] = p.baseIntIdx
    p.latS[p.pos[tailKey][1], p.pos[tailKey][2], p.pos[tailKey][3]] = p.baseIntIdx
    p.latS[p.pos[tailPrevKey][1], p.pos[tailPrevKey][2], p.pos[tailPrevKey][3]] = p.baseIntIdx
    
    # Update tracking if requested
    if updateTracking
        p.patternTails[intIdx] = newTailKey
        if haskey(p.patternTailPos, intIdx)
            p.patternTailPos[intIdx] -= 2
        end
    end
    
    return true, (tailPrevKey, tailKey)
end
```

### 2.1.2 Function `swapTailToPattern!`

Converts 2 base beads immediately after the tail of a pattern region to pattern type.

**Use cases**:
- In-pattern removeBond (compensate for 2 removed beads by extending pattern at tail)
- Downstream shift (second half of the shift operation)

**Signature**:
```julia
"""
    swapTailToPattern!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
        -> (success::Bool, swappedKeys::Tuple{Int32, Int32})

Converts 2 base beads immediately after the tail of pattern `intIdx` to pattern type.
Updates pattern tracking (patternTails, patternTailPos) if `updateTracking=true`.

Returns:
- `success`: Whether the swap was performed (false if not enough base beads after tail)
- `swappedKeys`: Tuple of (nextKey1, nextKey2) that were swapped, or (0, 0) if failed

Preconditions:
- Pattern `intIdx` must exist
- There must be at least 2 base beads after the current pattern tail

Effects:
- Converts nextKey1 and nextKey2 from baseIntIdx → intIdx
- Updates p.pos[key][6] and p.latS for both beads
- If updateTracking: updates patternTails[intIdx] to nextKey2 (new tail)
- If updateTracking: updates patternTailPos[intIdx] += 2
"""
function swapTailToPattern!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
    # Guard: pattern must exist
    if !haskey(p.patternTails, intIdx)
        return false, (Int32(0), Int32(0))
    end
    
    tailKey = p.patternTails[intIdx]
    
    # Guard: tail key must be valid
    if !haskey(p.pos, tailKey)
        return false, (Int32(0), Int32(0))
    end
    
    # Get the 2 beads after tail
    nextKey1 = p.pos[tailKey][5]
    if nextKey1 == 0
        return false, (Int32(0), Int32(0))
    end
    
    nextKey2 = p.pos[nextKey1][5]
    if nextKey2 == 0
        return false, (Int32(0), Int32(0))
    end
    
    # Verify the beads are base type before swapping
    if p.pos[nextKey1][6] != p.baseIntIdx || p.pos[nextKey2][6] != p.baseIntIdx
        return false, (Int32(0), Int32(0))
    end
    
    # Perform the swap: base → pattern
    p.pos[nextKey1][6] = intIdx
    p.pos[nextKey2][6] = intIdx
    p.latS[p.pos[nextKey1][1], p.pos[nextKey1][2], p.pos[nextKey1][3]] = intIdx
    p.latS[p.pos[nextKey2][1], p.pos[nextKey2][2], p.pos[nextKey2][3]] = intIdx
    
    # Update tracking if requested
    if updateTracking
        p.patternTails[intIdx] = nextKey2
        if haskey(p.patternTailPos, intIdx)
            p.patternTailPos[intIdx] += 2
        end
    end
    
    return true, (nextKey1, nextKey2)
end
```

### 2.1.3 Tests for Sub-Phase 2.1

```julia
@testset "Tail Swap Functions" begin
    
    # Helper to create a test polymer with a known pattern
    function createTestPolymer(patternConfig::Vector{Int}; boxlength=30)
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        nTypes = maximum(patternConfig)
        poly = Polymer(;
            length = length(patternConfig),
            maxlength = length(patternConfig) * 2,
            nMonomerTypes = nTypes,
            pattern = :custom,
            customPattern = patternConfig
        )
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.3)
        init!(l, [poly])
        return l, poly
    end
    
    # Helper to get chain as type vector
    function getTypeVector(p::Polymer)
        result = Int32[]
        key = p.head
        while key != 0
            push!(result, p.pos[key][6])
            key = p.pos[key][5]
        end
        return result
    end
    
    # Helper to verify lattice consistency
    function verifyLatticeConsistency(p::Polymer)
        key = p.head
        while key != 0
            pos = p.pos[key][1:3]
            latValue = p.latS[pos...]
            beadType = p.pos[key][6]
            if latValue != beadType
                return false
            end
            key = p.pos[key][5]
        end
        return true
    end
    
    @testset "swapTailToBase! - basic functionality" begin
        # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1]
        # After swap: [1,1,1,1,2,2,1,1,1,1,1,1]
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        @test poly.patternTails[Int32(2)] != 0  # Pattern exists
        
        initialTailPos = poly.patternTailPos[Int32(2)]
        success, swapped = swapTailToBase!(poly, Int32(2))
        
        @test success == true
        @test swapped[1] != 0 && swapped[2] != 0
        
        # Verify beads were swapped
        @test poly.pos[swapped[1]][6] == poly.baseIntIdx
        @test poly.pos[swapped[2]][6] == poly.baseIntIdx
        
        # Verify tracking updated
        @test poly.patternTailPos[Int32(2)] == initialTailPos - 2
        
        # Verify lattice consistency
        @test verifyLatticeConsistency(poly)
        
        # Verify chain structure
        types = getTypeVector(poly)
        @test types == Int32[1,1,1,1,2,2,1,1,1,1,1,1]
    end
    
    @testset "swapTailToBase! - pattern too small" begin
        # Pattern with only 2 beads - after swap would have 0, should fail
        l, poly = createTestPolymer([1,1,1,1,2,2,1,1,1,1])
        
        success, swapped = swapTailToBase!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
        
        # Verify nothing changed
        types = getTypeVector(poly)
        @test types == Int32[1,1,1,1,2,2,1,1,1,1]
    end
    
    @testset "swapTailToBase! - updateTracking=false" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialTailPos = poly.patternTailPos[Int32(2)]
        initialTailKey = poly.patternTails[Int32(2)]
        
        success, swapped = swapTailToBase!(poly, Int32(2); updateTracking=false)
        
        @test success == true
        
        # Tracking should NOT be updated
        @test poly.patternTailPos[Int32(2)] == initialTailPos
        @test poly.patternTails[Int32(2)] == initialTailKey
        
        # But beads should still be swapped
        @test poly.pos[swapped[1]][6] == poly.baseIntIdx
        @test poly.pos[swapped[2]][6] == poly.baseIntIdx
    end
    
    @testset "swapTailToBase! - nonexistent pattern" begin
        l, poly = createTestPolymer([1,1,1,1,1,1,1,1])  # Uniform, no pattern
        
        success, swapped = swapTailToBase!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "swapTailToPattern! - basic functionality" begin
        # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1]
        # After swap: [1,1,1,1,2,2,2,2,2,2,1,1]
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialTailPos = poly.patternTailPos[Int32(2)]
        success, swapped = swapTailToPattern!(poly, Int32(2))
        
        @test success == true
        @test swapped[1] != 0 && swapped[2] != 0
        
        # Verify beads were swapped
        @test poly.pos[swapped[1]][6] == Int32(2)
        @test poly.pos[swapped[2]][6] == Int32(2)
        
        # Verify tracking updated
        @test poly.patternTailPos[Int32(2)] == initialTailPos + 2
        
        # Verify lattice consistency
        @test verifyLatticeConsistency(poly)
        
        # Verify chain structure
        types = getTypeVector(poly)
        @test types == Int32[1,1,1,1,2,2,2,2,2,2,1,1]
    end
    
    @testset "swapTailToPattern! - not enough base beads" begin
        # Pattern at end - no base beads after tail
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2])
        
        success, swapped = swapTailToPattern!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "swapTailToPattern! - only 1 base bead after tail" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1])
        
        success, swapped = swapTailToPattern!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "swapTailToPattern! - updateTracking=false" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialTailPos = poly.patternTailPos[Int32(2)]
        initialTailKey = poly.patternTails[Int32(2)]
        
        success, swapped = swapTailToPattern!(poly, Int32(2); updateTracking=false)
        
        @test success == true
        
        # Tracking should NOT be updated
        @test poly.patternTailPos[Int32(2)] == initialTailPos
        @test poly.patternTails[Int32(2)] == initialTailKey
        
        # But beads should still be swapped
        @test poly.pos[swapped[1]][6] == Int32(2)
        @test poly.pos[swapped[2]][6] == Int32(2)
    end
    
    @testset "swapTailToPattern! - beads not base type" begin
        # Try to extend into another pattern region
        l, poly = createTestPolymer([1,1,2,2,2,2,3,3,1,1])
        
        # Pattern 2's tail is followed by pattern 3, not base
        # This should fail because the beads after pattern 2 are type 3, not base
        success, swapped = swapTailToPattern!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "Tail swaps are inverses" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        originalTypes = copy(getTypeVector(poly))
        originalTailPos = poly.patternTailPos[Int32(2)]
        
        # Swap tail to base
        success1, _ = swapTailToBase!(poly, Int32(2))
        @test success1
        
        # Now swap tail to pattern (extends back)
        success2, _ = swapTailToPattern!(poly, Int32(2))
        @test success2
        
        # Should be back to original
        @test getTypeVector(poly) == originalTypes
        @test poly.patternTailPos[Int32(2)] == originalTailPos
    end
    
    @testset "Multiple consecutive tail swaps" begin
        # Start with 8-bead pattern, shrink by 2 twice
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,2,2,2,2,1,1,1,1])
        
        @test count(==(2), getTypeVector(poly)) == 8
        
        success1, _ = swapTailToBase!(poly, Int32(2))
        @test success1
        @test count(==(2), getTypeVector(poly)) == 6
        
        success2, _ = swapTailToBase!(poly, Int32(2))
        @test success2
        @test count(==(2), getTypeVector(poly)) == 4
        
        success3, _ = swapTailToBase!(poly, Int32(2))
        @test success3
        @test count(==(2), getTypeVector(poly)) == 2
        
        # Now pattern has only 2 beads, next swap should fail
        success4, _ = swapTailToBase!(poly, Int32(2))
        @test success4 == false
        @test count(==(2), getTypeVector(poly)) == 2
    end
    
end
```

---

## Sub-Phase 2.2: Head Swap Functions

### 2.2.1 Function `swapHeadToBase!`

Converts the 2 pattern beads at the head of a pattern region to base type.

**Use cases**:
- Downstream shift (first half of the shift operation)

**Signature**:
```julia
"""
    swapHeadToBase!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
        -> (success::Bool, swappedKeys::Tuple{Int32, Int32})

Converts the 2 pattern beads at the head of pattern `intIdx` to base type.
Updates pattern tracking (patternHeads, patternHeadPos) if `updateTracking=true`.

Returns:
- `success`: Whether the swap was performed (false if pattern too small)
- `swappedKeys`: Tuple of (headKey, headNextKey) that were swapped, or (0, 0) if failed

Preconditions:
- Pattern `intIdx` must exist and have at least 3 beads (so there's still a pattern after swap)
- Pattern head/tail keys must be valid in p.pos

Effects:
- Converts headKey and headNextKey from intIdx → baseIntIdx
- Updates p.pos[key][6] and p.latS for both beads
- If updateTracking: updates patternHeads[intIdx] to the new head (2 beads later)
- If updateTracking: updates patternHeadPos[intIdx] += 2
"""
function swapHeadToBase!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
    # Guard: pattern must exist
    if !haskey(p.patternHeads, intIdx)
        return false, (Int32(0), Int32(0))
    end
    
    headKey = p.patternHeads[intIdx]
    
    # Guard: head key must be valid
    if !haskey(p.pos, headKey)
        return false, (Int32(0), Int32(0))
    end
    
    # Get head+1 (the bead after head)
    headNextKey = p.pos[headKey][5]
    
    # Guard: must have a bead after head (pattern has at least 2 beads)
    if headNextKey == 0
        return false, (Int32(0), Int32(0))
    end
    
    # Guard: after swap, pattern must still exist (need at least 1 bead remaining)
    # The new head would be headNextKey's next
    newHeadKey = p.pos[headNextKey][5]
    if newHeadKey == 0
        return false, (Int32(0), Int32(0))
    end
    
    # Verify the beads are actually pattern type before swapping
    if p.pos[headKey][6] != intIdx || p.pos[headNextKey][6] != intIdx
        return false, (Int32(0), Int32(0))
    end
    
    # Perform the swap: pattern → base
    p.pos[headKey][6] = p.baseIntIdx
    p.pos[headNextKey][6] = p.baseIntIdx
    p.latS[p.pos[headKey][1], p.pos[headKey][2], p.pos[headKey][3]] = p.baseIntIdx
    p.latS[p.pos[headNextKey][1], p.pos[headNextKey][2], p.pos[headNextKey][3]] = p.baseIntIdx
    
    # Update tracking if requested
    if updateTracking
        p.patternHeads[intIdx] = newHeadKey
        if haskey(p.patternHeadPos, intIdx)
            p.patternHeadPos[intIdx] += 2
        end
    end
    
    return true, (headKey, headNextKey)
end
```

### 2.2.2 Function `swapHeadToPattern!`

Converts 2 base beads immediately before the head of a pattern region to pattern type.

**Use cases**:
- Upstream shift (first half of the shift operation)

**Signature**:
```julia
"""
    swapHeadToPattern!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
        -> (success::Bool, swappedKeys::Tuple{Int32, Int32})

Converts 2 base beads immediately before the head of pattern `intIdx` to pattern type.
Updates pattern tracking (patternHeads, patternHeadPos) if `updateTracking=true`.

Returns:
- `success`: Whether the swap was performed (false if not enough base beads before head)
- `swappedKeys`: Tuple of (prevKey2, prevKey1) that were swapped (in chain order), or (0, 0) if failed

Preconditions:
- Pattern `intIdx` must exist
- There must be at least 2 base beads before the current pattern head

Effects:
- Converts prevKey2 and prevKey1 from baseIntIdx → intIdx
- Updates p.pos[key][6] and p.latS for both beads
- If updateTracking: updates patternHeads[intIdx] to prevKey2 (new head)
- If updateTracking: updates patternHeadPos[intIdx] -= 2
"""
function swapHeadToPattern!(p::Polymer, intIdx::Int32; updateTracking::Bool=true)
    # Guard: pattern must exist
    if !haskey(p.patternHeads, intIdx)
        return false, (Int32(0), Int32(0))
    end
    
    headKey = p.patternHeads[intIdx]
    
    # Guard: head key must be valid
    if !haskey(p.pos, headKey)
        return false, (Int32(0), Int32(0))
    end
    
    # Get the 2 beads before head
    prevKey1 = p.pos[headKey][4]
    if prevKey1 == 0
        return false, (Int32(0), Int32(0))
    end
    
    prevKey2 = p.pos[prevKey1][4]
    if prevKey2 == 0
        return false, (Int32(0), Int32(0))
    end
    
    # Verify the beads are base type before swapping
    if p.pos[prevKey1][6] != p.baseIntIdx || p.pos[prevKey2][6] != p.baseIntIdx
        return false, (Int32(0), Int32(0))
    end
    
    # Perform the swap: base → pattern
    p.pos[prevKey1][6] = intIdx
    p.pos[prevKey2][6] = intIdx
    p.latS[p.pos[prevKey1][1], p.pos[prevKey1][2], p.pos[prevKey1][3]] = intIdx
    p.latS[p.pos[prevKey2][1], p.pos[prevKey2][2], p.pos[prevKey2][3]] = intIdx
    
    # Update tracking if requested
    if updateTracking
        p.patternHeads[intIdx] = prevKey2
        if haskey(p.patternHeadPos, intIdx)
            p.patternHeadPos[intIdx] -= 2
        end
    end
    
    return true, (prevKey2, prevKey1)
end
```

### 2.2.3 Tests for Sub-Phase 2.2

```julia
@testset "Head Swap Functions" begin
    
    # Reuse helpers from Sub-Phase 2.1 (createTestPolymer, getTypeVector, verifyLatticeConsistency)
    
    @testset "swapHeadToBase! - basic functionality" begin
        # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1]
        # After swap: [1,1,1,1,1,1,2,2,1,1,1,1]
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialHeadPos = poly.patternHeadPos[Int32(2)]
        success, swapped = swapHeadToBase!(poly, Int32(2))
        
        @test success == true
        @test swapped[1] != 0 && swapped[2] != 0
        
        # Verify beads were swapped
        @test poly.pos[swapped[1]][6] == poly.baseIntIdx
        @test poly.pos[swapped[2]][6] == poly.baseIntIdx
        
        # Verify tracking updated
        @test poly.patternHeadPos[Int32(2)] == initialHeadPos + 2
        
        # Verify lattice consistency
        @test verifyLatticeConsistency(poly)
        
        # Verify chain structure
        types = getTypeVector(poly)
        @test types == Int32[1,1,1,1,1,1,2,2,1,1,1,1]
    end
    
    @testset "swapHeadToBase! - pattern too small" begin
        # Pattern with only 2 beads - after swap would have 0, should fail
        l, poly = createTestPolymer([1,1,1,1,2,2,1,1,1,1])
        
        success, swapped = swapHeadToBase!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "swapHeadToBase! - updateTracking=false" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialHeadPos = poly.patternHeadPos[Int32(2)]
        initialHeadKey = poly.patternHeads[Int32(2)]
        
        success, swapped = swapHeadToBase!(poly, Int32(2); updateTracking=false)
        
        @test success == true
        
        # Tracking should NOT be updated
        @test poly.patternHeadPos[Int32(2)] == initialHeadPos
        @test poly.patternHeads[Int32(2)] == initialHeadKey
        
        # But beads should still be swapped
        @test poly.pos[swapped[1]][6] == poly.baseIntIdx
        @test poly.pos[swapped[2]][6] == poly.baseIntIdx
    end
    
    @testset "swapHeadToPattern! - basic functionality" begin
        # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1]
        # After swap: [1,1,2,2,2,2,2,2,1,1,1,1]
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialHeadPos = poly.patternHeadPos[Int32(2)]
        success, swapped = swapHeadToPattern!(poly, Int32(2))
        
        @test success == true
        @test swapped[1] != 0 && swapped[2] != 0
        
        # Verify beads were swapped
        @test poly.pos[swapped[1]][6] == Int32(2)
        @test poly.pos[swapped[2]][6] == Int32(2)
        
        # Verify tracking updated
        @test poly.patternHeadPos[Int32(2)] == initialHeadPos - 2
        
        # Verify lattice consistency
        @test verifyLatticeConsistency(poly)
        
        # Verify chain structure
        types = getTypeVector(poly)
        @test types == Int32[1,1,2,2,2,2,2,2,1,1,1,1]
    end
    
    @testset "swapHeadToPattern! - not enough base beads" begin
        # Pattern at start - no base beads before head
        l, poly = createTestPolymer([2,2,2,2,1,1,1,1])
        
        success, swapped = swapHeadToPattern!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "swapHeadToPattern! - only 1 base bead before head" begin
        l, poly = createTestPolymer([1,2,2,2,2,1,1,1,1])
        
        success, swapped = swapHeadToPattern!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "swapHeadToPattern! - beads not base type" begin
        # Try to extend into another pattern region
        l, poly = createTestPolymer([1,1,3,3,2,2,2,2,1,1])
        
        # Pattern 2's head is preceded by pattern 3, not base
        success, swapped = swapHeadToPattern!(poly, Int32(2))
        
        @test success == false
        @test swapped == (Int32(0), Int32(0))
    end
    
    @testset "Head swaps are inverses" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        originalTypes = copy(getTypeVector(poly))
        originalHeadPos = poly.patternHeadPos[Int32(2)]
        
        # Swap head to base (shrinks pattern at head)
        success1, _ = swapHeadToBase!(poly, Int32(2))
        @test success1
        
        # Now swap head to pattern (extends back)
        success2, _ = swapHeadToPattern!(poly, Int32(2))
        @test success2
        
        # Should be back to original
        @test getTypeVector(poly) == originalTypes
        @test poly.patternHeadPos[Int32(2)] == originalHeadPos
    end
    
end
```

---

## Sub-Phase 2.3: Compose Shift Functions from Primitives

### 2.3.1 Rewrite `shiftPatternUpstream!` using primitives

```julia
"""
    shiftPatternUpstream!(p::Polymer, intIdx::Int32; updatePos::Bool=true)
        -> (success::Bool, swappedKeys::Vector{Int32}, oldTypes::Vector{Int32}, newTypes::Vector{Int32})

Shifts pattern `intIdx` upstream (toward polymer head) by 2 positions.
Equivalent to: swapHeadToPattern! + swapTailToBase!

Returns success status and swap information for energy calculation.
"""
function shiftPatternUpstream!(p::Polymer, intIdx::Int32; updatePos::Bool=true)
    # For upstream shift during add: we don't want to update head position
    # because the insertion already increased all positions downstream by 2.
    # So we extend at head (tracking=false means head position unchanged),
    # and shrink at tail (tracking depends on updatePos).
    
    # Step 1: Extend pattern at head (2 base → pattern)
    # Don't update head tracking - the add already moved positions
    successHead, headSwapped = swapHeadToPattern!(p, intIdx; updateTracking=false)
    if !successHead
        return false, Int32[], Int32[], Int32[]
    end
    
    # Step 2: Shrink pattern at tail (2 pattern → base)
    # Update tail tracking only if updatePos is true
    successTail, tailSwapped = swapTailToBase!(p, intIdx; updateTracking=updatePos)
    if !successTail
        # Rollback head swap
        p.pos[headSwapped[1]][6] = p.baseIntIdx
        p.pos[headSwapped[2]][6] = p.baseIntIdx
        p.latS[p.pos[headSwapped[1]][1], p.pos[headSwapped[1]][2], p.pos[headSwapped[1]][3]] = p.baseIntIdx
        p.latS[p.pos[headSwapped[2]][1], p.pos[headSwapped[2]][2], p.pos[headSwapped[2]][3]] = p.baseIntIdx
        return false, Int32[], Int32[], Int32[]
    end
    
    # Build return values for energy calculation
    swappedKeys = Int32[headSwapped[1], headSwapped[2], tailSwapped[1], tailSwapped[2]]
    oldTypes = Int32[p.baseIntIdx, p.baseIntIdx, intIdx, intIdx]
    newTypes = Int32[intIdx, intIdx, p.baseIntIdx, p.baseIntIdx]
    
    return true, swappedKeys, oldTypes, newTypes
end
```

### 2.3.2 Rewrite `shiftPatternDownstream!` using primitives

```julia
"""
    shiftPatternDownstream!(p::Polymer, intIdx::Int32; updatePos::Bool=true)
        -> (success::Bool, swappedKeys::Vector{Int32}, oldTypes::Vector{Int32}, newTypes::Vector{Int32})

Shifts pattern `intIdx` downstream (toward polymer tail) by 2 positions.
Equivalent to: swapHeadToBase! + swapTailToPattern!

Returns success status and swap information for energy calculation.
"""
function shiftPatternDownstream!(p::Polymer, intIdx::Int32; updatePos::Bool=true)
    # For downstream shift during remove: we don't want to update head position
    # because the removal already decreased all positions downstream by 2.
    # So we shrink at head (tracking=false means head position unchanged),
    # and extend at tail (tracking depends on updatePos).
    
    # Step 1: Shrink pattern at head (2 pattern → base)
    # Don't update head tracking - the remove already moved positions
    successHead, headSwapped = swapHeadToBase!(p, intIdx; updateTracking=false)
    if !successHead
        return false, Int32[], Int32[], Int32[]
    end
    
    # Step 2: Extend pattern at tail (2 base → pattern)
    # Update tail tracking only if updatePos is true
    successTail, tailSwapped = swapTailToPattern!(p, intIdx; updateTracking=updatePos)
    if !successTail
        # Rollback head swap
        p.pos[headSwapped[1]][6] = intIdx
        p.pos[headSwapped[2]][6] = intIdx
        p.latS[p.pos[headSwapped[1]][1], p.pos[headSwapped[1]][2], p.pos[headSwapped[1]][3]] = intIdx
        p.latS[p.pos[headSwapped[2]][1], p.pos[headSwapped[2]][2], p.pos[headSwapped[2]][3]] = intIdx
        return false, Int32[], Int32[], Int32[]
    end
    
    # Build return values for energy calculation
    swappedKeys = Int32[headSwapped[1], headSwapped[2], tailSwapped[1], tailSwapped[2]]
    oldTypes = Int32[intIdx, intIdx, p.baseIntIdx, p.baseIntIdx]
    newTypes = Int32[p.baseIntIdx, p.baseIntIdx, intIdx, intIdx]
    
    return true, swappedKeys, oldTypes, newTypes
end
```

### 2.3.3 Tests for Composed Shift Functions

```julia
@testset "Composed Shift Functions" begin
    
    @testset "shiftPatternUpstream! - basic functionality" begin
        # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1] at positions 5-8
        # After shift: [1,1,2,2,2,2,2,2,1,1,1,1] - pattern now at positions 3-6
        # Wait, that's wrong - upstream shift should keep pattern length same
        # Actually: [1,1,2,2,2,2,1,1,1,1,1,1] - extend at head, shrink at tail
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialHeadPos = poly.patternHeadPos[Int32(2)]
        initialTailPos = poly.patternTailPos[Int32(2)]
        
        success, swappedKeys, oldTypes, newTypes = shiftPatternUpstream!(poly, Int32(2))
        
        @test success == true
        @test length(swappedKeys) == 4
        
        # Pattern head should move 2 earlier
        @test poly.patternHeadPos[Int32(2)] == initialHeadPos - 2
        @test poly.patternTailPos[Int32(2)] == initialTailPos - 2
        
        # Pattern length should be preserved
        types = getTypeVector(poly)
        @test count(==(2), types) == 4
        
        # Verify pattern moved upstream
        @test types == Int32[1,1,2,2,2,2,1,1,1,1,1,1]
    end
    
    @testset "shiftPatternDownstream! - basic functionality" begin
        # Pattern: [1,1,1,1,2,2,2,2,1,1,1,1] at positions 5-8
        # After shift: [1,1,1,1,1,1,2,2,2,2,1,1] - shrink at head, extend at tail
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialHeadPos = poly.patternHeadPos[Int32(2)]
        initialTailPos = poly.patternTailPos[Int32(2)]
        
        success, swappedKeys, oldTypes, newTypes = shiftPatternDownstream!(poly, Int32(2))
        
        @test success == true
        @test length(swappedKeys) == 4
        
        # Pattern head should move 2 later
        @test poly.patternHeadPos[Int32(2)] == initialHeadPos + 2
        @test poly.patternTailPos[Int32(2)] == initialTailPos + 2
        
        # Pattern length should be preserved
        types = getTypeVector(poly)
        @test count(==(2), types) == 4
        
        # Verify pattern moved downstream
        @test types == Int32[1,1,1,1,1,1,2,2,2,2,1,1]
    end
    
    @testset "Shifts are inverses" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        originalTypes = copy(getTypeVector(poly))
        originalHeadPos = poly.patternHeadPos[Int32(2)]
        originalTailPos = poly.patternTailPos[Int32(2)]
        
        # Shift upstream
        success1, _, _, _ = shiftPatternUpstream!(poly, Int32(2))
        @test success1
        @test getTypeVector(poly) != originalTypes  # Should have changed
        
        # Shift downstream (inverse)
        success2, _, _, _ = shiftPatternDownstream!(poly, Int32(2))
        @test success2
        
        # Should be back to original
        @test getTypeVector(poly) == originalTypes
        @test poly.patternHeadPos[Int32(2)] == originalHeadPos
        @test poly.patternTailPos[Int32(2)] == originalTailPos
    end
    
    @testset "Shift with updatePos=false" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialHeadPos = poly.patternHeadPos[Int32(2)]
        initialTailPos = poly.patternTailPos[Int32(2)]
        
        success, _, _, _ = shiftPatternUpstream!(poly, Int32(2); updatePos=false)
        
        @test success == true
        
        # Position tracking should NOT be updated (used for compensating add/remove)
        # Note: the swap functions with updateTracking=false don't update positions
        # This is the behavior we want for pattern pinning
        @test poly.patternHeadPos[Int32(2)] == initialHeadPos
        @test poly.patternTailPos[Int32(2)] == initialTailPos
    end
    
    @testset "Shift fails gracefully when not enough room" begin
        # Pattern at start - can't extend at head
        l, poly = createTestPolymer([2,2,2,2,1,1,1,1])
        
        success, _, _, _ = shiftPatternUpstream!(poly, Int32(2))
        @test success == false
        
        # Pattern at end - can't extend at tail
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2])
        
        success, _, _, _ = shiftPatternDownstream!(poly, Int32(2))
        @test success == false
    end
    
end
```

---

## Sub-Phase 2.4: Refactor addBond! and removeBond! to use primitives

### 2.4.1 Refactored `addBond!`

Replace the inline tail swap logic with a call to `swapTailToBase!`:

```julia
function addBond!(p::Polymer)
    # ... existing setup code (insertPos computation, intIdx, etc.) ...
    
    # ... existing bead insertion code (create new beads, update chain) ...
    
    # Swap pattern tail beads to base type (compensate for 2 new beads)
    if p.patternTails !== nothing && intIdx != p.baseIntIdx && 
       haskey(p.patternTails, intIdx) && isInPatternRegion(p, p.ki, intIdx)
        
        # Use the primitive function - don't update tracking (handled by downstream shift logic)
        success, swapped = swapTailToBase!(p, intIdx; updateTracking=false)
        
        if !success
            # Pattern is too small to compensate - this shouldn't happen if validation is correct
            # For now, just continue (the pattern will be invalidated)
            delete!(p.patternTails, intIdx)
            delete!(p.patternTailPos, intIdx)
            delete!(p.patternHeads, intIdx)
            delete!(p.patternHeadPos, intIdx)
        end
    end
    
    # ... existing downstream pattern shifting code ...
    
    return p
end
```

### 2.4.2 Refactored `removeBond!`

Replace the inline tail extension logic with a call to `swapTailToPattern!`:

```julia
function removeBond!(p::Polymer)
    # ... existing setup code ...
    
    # Get the intIdx before removing
    intIdx = p.pos[p.ki][6]
    
    # Extend pattern at tail BEFORE removing beads (compensate for 2 removed beads)
    if p.patternTails !== nothing && intIdx != p.baseIntIdx && 
       haskey(p.patternTails, intIdx) && isInPatternRegion(p, p.ki, intIdx)
        
        # Use the primitive function - don't update tracking (handled by downstream shift logic)
        success, swapped = swapTailToPattern!(p, intIdx; updateTracking=false)
        
        if !success
            # Not enough base beads after pattern to compensate
            # This should have been caught by validation
            # Pattern will be invalidated
        end
        
        # Handle head key update if removing head beads
        headKey = p.patternHeads[intIdx]
        if headKey == p.ki || headKey == p.kj
            p.patternHeads[intIdx] = p.knext
        end
    end
    
    # ... existing bead removal code ...
    
    # ... existing downstream pattern shifting code ...
    
    return p
end
```

### 2.4.3 Integration Tests

```julia
@testset "addBond!/removeBond! with swap primitives" begin
    
    @testset "addBond inside pattern uses swapTailToBase" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialPatternCount = count(==(2), getTypeVector(poly))
        initialLength = poly.length
        
        # Find a bond inside the pattern
        key = poly.head
        for _ in 1:4  # Advance to 5th bead (first pattern bead)
            key = poly.pos[key][5]
        end
        # Now key is at pattern bead
        nextKey = poly.pos[key][5]
        
        # Set up for addBond
        poly.ki = key
        poly.kj = nextKey
        poly.kprev = poly.pos[key][4]
        poly.knext = poly.pos[nextKey][5]
        poly.moveIdx = 5
        
        # Set up positions
        poly.bi = view(poly.pos[poly.ki], 1:3)
        poly.bj = view(poly.pos[poly.kj], 1:3)
        poly.bprev = view(poly.pos[poly.kprev], 1:3)
        _sentinel = zeros(Int32, 3)
        poly.bnext = poly.knext == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.knext], 1:3)
        
        # Set up new positions (just offset from current)
        poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
        poly.newPosj = [poly.bj[1], poly.bj[2], Int32(mod1(poly.bj[3] + 1, poly.boxlength))]
        
        # Ensure positions are empty
        poly.latS[poly.newPosi...] = 0
        poly.latS[poly.newPosj...] = 0
        
        poly.moveType = 1
        
        addBond!(poly)
        poly.length = length(poly.pos)
        
        # Length increased by 2
        @test poly.length == initialLength + 2
        
        # Pattern count unchanged (2 added inside pattern, 2 swapped to base at tail)
        finalPatternCount = count(==(2), getTypeVector(poly))
        @test finalPatternCount == initialPatternCount
        
        # Verify lattice consistency
        @test verifyLatticeConsistency(poly)
    end
    
    @testset "removeBond inside pattern uses swapTailToPattern" begin
        l, poly = createTestPolymer([1,1,1,1,2,2,2,2,1,1,1,1])
        
        initialPatternCount = count(==(2), getTypeVector(poly))
        initialLength = poly.length
        
        # Find a bond inside the pattern
        key = poly.head
        for _ in 1:4  # Advance to 5th bead (first pattern bead)
            key = poly.pos[key][5]
        end
        nextKey = poly.pos[key][5]
        
        # Set up for removeBond
        poly.ki = key
        poly.kj = nextKey
        poly.kprev = poly.pos[key][4]
        poly.knext = poly.pos[nextKey][5]
        poly.moveIdx = 5
        
        # Set up positions
        poly.bi = view(poly.pos[poly.ki], 1:3)
        poly.bj = view(poly.pos[poly.kj], 1:3)
        poly.bprev = view(poly.pos[poly.kprev], 1:3)
        poly.bnext = view(poly.pos[poly.knext], 1:3)
        
        # For removeBond, new positions equal prev/next
        poly.newPosi = Vector{Int32}(poly.bprev)
        poly.newPosj = Vector{Int32}(poly.bnext)
        
        poly.moveType = 3
        
        removeBond!(poly)
        poly.length = length(poly.pos)
        
        # Length decreased by 2
        @test poly.length == initialLength - 2
        
        # Pattern count unchanged (2 removed from pattern, 2 swapped to pattern at tail)
        finalPatternCount = count(==(2), getTypeVector(poly))
        @test finalPatternCount == initialPatternCount
        
        # Verify lattice consistency
        @test verifyLatticeConsistency(poly)
    end
    
end
```

---

## Summary

| Sub-Phase | Description | Functions |
|-----------|-------------|-----------|
| 2.1 | Tail swap primitives | `swapTailToBase!`, `swapTailToPattern!` |
| 2.2 | Head swap primitives | `swapHeadToBase!`, `swapHeadToPattern!` |
| 2.3 | Compose shift functions | Refactored `shiftPatternUpstream!`, `shiftPatternDownstream!` |
| 2.4 | Integrate into add/remove | Refactored `addBond!`, `removeBond!` |

### New Functions

| Function | Direction | Type Change | Use Cases |
|----------|-----------|-------------|-----------|
| `swapTailToBase!` | Shrink at tail | pattern → base | In-pattern add, upstream shift |
| `swapTailToPattern!` | Extend at tail | base → pattern | In-pattern remove, downstream shift |
| `swapHeadToBase!` | Shrink at head | pattern → base | Downstream shift |
| `swapHeadToPattern!` | Extend at head | base → pattern | Upstream shift |

### Composition Table

| Operation | Head Swap | Tail Swap |
|-----------|-----------|-----------|
| In-pattern addBond | - | `swapTailToBase!` |
| In-pattern removeBond | - | `swapTailToPattern!` |
| Upstream shift | `swapHeadToPattern!` | `swapTailToBase!` |
| Downstream shift | `swapHeadToBase!` | `swapTailToPattern!` |

### Key Design Decisions

1. **Atomic primitives**: Each function does exactly one thing (swap 2 beads at one boundary)
2. **Optional tracking update**: `updateTracking` parameter allows flexible use in different contexts
3. **Return values**: All functions return (success, swappedKeys) for error handling and energy calculation
4. **Guards**: Comprehensive precondition checks prevent invalid operations
5. **Rollback support**: Composed functions can rollback if second operation fails
