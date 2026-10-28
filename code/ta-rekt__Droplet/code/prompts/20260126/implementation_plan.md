# Implementation Plan: Heterogeneous Polymers

## Overview

This plan implements heterogeneous polymers with multiple monomer types, where patterns of specific types maintain fixed length while the base type region can grow/shrink.

---

## Phase 0: Lattice Type-Tracking and Convention Establishment

### 0.1 Convention: Use `intIdx` as the canonical monomer type identifier

**Key insight**: The interaction index (`m.intIdx`) assigned by `_ensureIntStructures!` and `_syncMatrices!` is:
- Unique across all interacting species in the lattice
- Already computed and stored on each monomer
- Used in interaction matrix lookups

**New approach**: Store `monomers` as a **dictionary keyed by `intIdx`** instead of a vector:
- `p.monomers::Dict{Int32, PolyMonomer}` with keys being the `intIdx` values
- `pos[key][6]` stores the `intIdx` directly (not an array index)
- Direct lookup: `m = p.monomers[pos[key][6]]` - O(1) access without indirection
- **Eliminates the need for the separate `monomerIndex` array**

**Benefits of this approach**:
1. **Direct equivalence**: `pos[key][6] == m.intIdx` - no conversion needed
2. **Single source of truth**: Same value in `pos[key][6]`, `latS[x,y,z]`, and dictionary key
3. **Simpler code**: One dict lookup instead of array index → intIdx conversion
4. **Cross-polymer compatibility**: Works for monomers from any polymer using the same intIdx

### 0.2 Store `intIdx` in lattice instead of 0/1

**Current behavior**:
- `p.latS[x, y, z] = 1` when any polymer bead occupies the site
- `p.latS[x, y, z] = 0` when empty
- `nnp(s, p, pos)` counts sites where `p.latS[...] > 0`

**New behavior**:
- `p.latS[x, y, z] = m.intIdx` when a bead of monomer type `m` occupies the site
- `p.latS[x, y, z] = 0` when empty
- Existing `nnp(s, p, pos)` still works: `> 0` still detects any polymer bead
- New `nnp(s::PolyMonomer, p::PolyMonomer, pos)` dispatch counts only sites where `latS == p.intIdx`

### 0.3 New `nnp` dispatch for PolyMonomer interactions

**File**: `src/utils.jl`

Add a specialized dispatch that counts only beads of a specific monomer type:

```julia
# PolyMonomer-to-PolyMonomer interaction: count only beads of target type
function nnp(s::AbstractPolyMonomer, 
             p::AbstractPolyMonomer, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(p)           # p.poly.latS (or p.latS if we add that field)
    targetIntIdx = p.intIdx       # the interaction index we're looking for
    
    for i in 1:3
        pos[i] += 1
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == targetIntIdx)

        pos[i] -= 2
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == targetIntIdx)

        pos[i] += 1
    end

    clamp!(tempVec, 0, 0)
    return n
end
```

**Why this works**:
- `intIdx` is unique across ALL species in the lattice (polymers, bulk monomers, etc.)
- When checking `latS == p.intIdx`, we only count beads of that specific monomer type
- This works even for cross-polymer interactions if two polymers share interacting monomers

### 0.4 Update lattice population in `init!`

**File**: `src/init.jl`, in `init!` function

**Current** (around line 445):
```julia
for j in 1:s.length
    s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
end
```

**New** (for polymers with monomers) - since `pos[key][6]` now stores `intIdx` directly:
```julia
if typeof(s) <: AbstractPolymer && s.monomers !== nothing
    key = s.head
    while key != 0
        pos = s.pos[key]
        intIdx = pos[6]  # Now directly stores intIdx
        s.latS[pos[1], pos[2], pos[3]] = intIdx
        key = pos[5]
    end
else
    for j in 1:s.length
        s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
    end
end
```

### 0.5 Update `addBond!`, `kinkBond!`, `removeBond!` to use dictionary lookup

**File**: `src/move.jl`

**`addBond!`** - direct lookup by intIdx:
```julia
# Get the current monomer type (intIdx stored directly in pos)
intIdx = p.pos[p.ki][6]
m = p.monomers[intIdx]  # Dictionary lookup by intIdx

# Set new lattice sites to intIdx
p.latS[p.newPosi...] = intIdx
p.latS[p.newPosj...] = intIdx
```

**`kinkBond!`** - preserve the monomer's `intIdx` when moving:
```julia
if p.newPosi == p.bprev 
    # bead ki moves from bi to newPosj
    intIdx = p.pos[p.ki][6]
    p.latS[p.bi[1], p.bi[2], p.bi[3]] = 0
    p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] = intIdx
    # ... rest unchanged
elseif p.newPosj == p.bnext
    # bead kj moves from bj to newPosi
    intIdx = p.pos[p.kj][6]
    p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] = intIdx
    p.latS[p.bj[1], p.bj[2], p.bj[3]] = 0
    # ... rest unchanged
end
```

**`removeBond!`** - no changes needed (already sets to 0)

### 0.6 Update type swaps (Phase 5) to update lattice

When beads swap types (in `addBond!` and `removeBond!` for pattern preservation), update both `pos[key][6]` and lattice with the new `intIdx`:

```julia
# Swap a bead from oldIntIdx to newIntIdx
oldIntIdx = p.pos[key][6]
newIntIdx = targetMonomerIntIdx  # intIdx of new monomer type

# Update pos (now stores intIdx directly)
p.pos[key][6] = newIntIdx

# Update lattice (same intIdx value)
beadPos = p.pos[key][1:3]
p.latS[beadPos[1], beadPos[2], beadPos[3]] = newIntIdx

# Get monomers for energy calculations
mOld = p.monomers[oldIntIdx]
mNew = p.monomers[newIntIdx]
```

### 0.7 `nn` remains type-agnostic

**Important**: The `nn(s, pos)` function should remain type-agnostic and count ALL polymer beads (checking `latS > 0`). This is correct because:
- `nn` is used for polymer self-interaction in the current homogeneous case
- For heterogeneous polymers, self-interactions should use `nnp(m, m, pos)` instead

**Enforcement**: In heterogeneous polymer code paths, ensure `nn` is not used for energy calculations. Use `nnp(m, m, pos)` for same-type self-interactions.

### 0.8 Overlap check remains unchanged

The `overlap(p)` function checks `p.latS[pos...] == 1`. This needs to be updated to:
```julia
function overlap(p::Polymer)
    return ((p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] > 0 && 
             !coords_equal(p.bprev, p.newPosi)) ||
            (p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] > 0 && 
             !coords_equal(p.bnext, p.newPosj)))
end
```

Change `== 1` to `> 0` since we now store `intIdx` values (which are always >= 1 for occupied sites).

### 0.9 Tests for Phase 0

Create `test/test_heterogeneous_phase0.jl`:

1. **Test lattice stores intIdx**: Initialize heterogeneous polymer, verify `latS[pos...] == m.intIdx` for each bead
2. **Test nnp dispatch for PolyMonomer**: 
   - Set up polymer with pattern `[1,1,2,2,1,1]`
   - Call `nnp(m1, m1, pos)` for a type-1 bead adjacent to type-1 and type-2 beads
   - Verify it only counts type-1 neighbors
3. **Test nnp dispatch for PolyMonomer self-interaction**:
   - `nnp(m2, m2, pos)` for a type-2 bead
   - Verify it only counts type-2 neighbors
4. **Test overlap with intIdx values**: Verify overlap detection still works with `> 0` check
5. **Test cross-species nnp**: `nnp(polyMon, bulkMon, pos)` still works correctly

---

## Phase 1: Store Monomer Type in `pos` (6-element vectors) and Dictionary-Based Monomers

### 1.1 Modify `pos` vector structure

**Goal**: Change `pos` entries from 5-element to 6-element vectors, where element 6 is the monomer's `intIdx`.

**Files to modify**:
- `src/types.jl`: Update comments for `pos` field; change `monomers` to Dict
- `src/init.jl`: All functions that create `pos` entries

**Changes**:

1. **`initUniform!`** (init.jl:17-31): Change `zeros(5)` → `zeros(6)` (intIdx assigned later in `initPattern!`)
2. **`initLine!`** (init.jl:33-48): Same pattern
3. **`initBox!`** (init.jl:50-90): Same pattern
4. **`initShape!(::Polymer)`** (init.jl:235-245): Already calls `initLine!`, no direct changes

5. **`addBond!`** (move.jl:15-57):
   - Change `Vector{Int32}(undef, 5)` → `Vector{Int32}(undef, 6)`
   - Set `p.pos[next_key1][6] = p.pos[p.ki][6]` (inherit intIdx from parent bead)

6. **`kinkBond!`** (move.jl:59-90):
   - Change `Vector{Int32}(undef, 5)` → `Vector{Int32}(undef, 6)`
   - Preserve monomer type: copy element 6 from old entry

7. **`removeBond!`** (move.jl:92-120): No changes—just deletes entries

### 1.2 Change `monomers` from Vector to Dictionary

**Goal**: Store monomers in a dictionary keyed by `intIdx` for O(1) direct lookup.

**File**: `src/types.jl` - Polymer struct

**Old field**:
```julia
monomers::Union{Nothing, Vector{PolyMonomer}} = nothing
monomerIndex::Union{Nothing, Vector{Int32}} = nothing  # no longer needed
```

**New field**:
```julia
monomers::Union{Nothing, Dict{Int32, PolyMonomer}} = nothing
```

**Remove `monomerIndex`** - no longer needed since `pos[key][6]` stores `intIdx` directly.

### 1.3 Update `initMonomerTypes!` to build dictionary

**File**: `src/init.jl`

```julia
function initMonomerTypes!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)
    # Create monomers dictionary keyed by intIdx
    p.monomers = Dict{Int32, PolyMonomer}()
    
    for i in 1:p.nMonomerTypes
        m = PolyMonomer(
            name = string(p.name, "_monomer", i),
            poly = p,
            chemPot = defaultChemPot,
            intPot = defaultIntPot
        )
        # intIdx assigned by _ensureIntStructures!
        _ensureIntStructures!(m)
        p.monomers[m.intIdx] = m  # Key by intIdx
    end
    
    # Sync matrices across all monomers
    if p.nMonomerTypes > 1
        monomers_vec = collect(values(p.monomers))
        for i in 2:length(monomers_vec)
            if monomers_vec[1].intMat !== monomers_vec[i].intMat
                _syncMatrices!(monomers_vec[1], monomers_vec[i])
            end
        end
    end
    
    return p
end
```

### 1.4 Update `getPolyMonomer` for dictionary lookup

**File**: `src/utils.jl`

```julia
function getPolyMonomer(p::Polymer, key::Int32)
    if p.monomers === nothing
        error("Polymer has no monomers")
    end
    intIdx = p.pos[key][6]
    return p.monomers[intIdx]  # Dictionary lookup by intIdx
end
```

### 1.5 Update Hamiltonian calculations

**`ham!`** (hamiltonian.jl:7-50): Iterate chain and get per-bead monomer via dictionary:

```julia
key = p.head
for i in 1:p.length
    pos = p.pos[key][1:3]
    intIdx = p.pos[key][6]
    m = p.monomers[intIdx]  # Dictionary lookup
    p.ham += m.intPot * nnp(m, m, pos) / 2.0
    key = p.pos[key][5]
end
```

**`deltaHam!`** (hamiltonian.jl:53-185): Use `getPolyMonomer`:

```julia
m = getPolyMonomer(p, p.ki)  # looks up p.pos[p.ki][6] in dictionary
```

### 1.6 Tests for Phase 1

Create `test/test_heterogeneous_phase1.jl`:

1. **Test 6-element pos vectors**: Initialize polymer, verify all entries have length 6
2. **Test monomer dictionary**: Verify `p.monomers` is a Dict keyed by intIdx
3. **Test dictionary lookup**: Verify `p.monomers[intIdx]` returns correct monomer
4. **Test `getPolyMonomer`**: Verify it correctly looks up monomer by key
5. **Test Hamiltonian calculation**: Verify `ham!` works with dictionary-based lookup

---

## Phase 2: Pattern Initialization with Chain Traversal

### 2.1 `initPattern!` function

**Note**: The `initMonomerTypes!` function is already defined in Phase 1.3. This phase adds `initPattern!` which writes monomer types (as `intIdx` values) to `pos[key][6]`.

```julia
function initPattern!(p::Polymer)
    # Fast path for homogeneous polymers
    if p.nMonomerTypes == 1 || p.pattern == :uniform
        # All beads get the base type's intIdx
        baseIntIdx = first(keys(p.monomers))  # any monomer works for uniform
        key = p.head
        while key != 0
            p.pos[key][6] = baseIntIdx
            key = p.pos[key][5]
        end
        p.patternTails = nothing
        return p
    end
    
    # Build pattern vector based on p.pattern (returns array indices 1,2,3...)
    patternVec = buildPatternVector(p)
    
    # Validate single contiguous block per pattern type
    validateSinglePatternBlocks(patternVec)
    
    # Map array indices to intIdx values
    # Collect intIdx values in a consistent order
    intIdxList = sort(collect(keys(p.monomers)))  # sorted list of intIdx values
    
    # Traverse chain and assign intIdx to pos[key][6]
    key = p.head
    for i in 1:p.length
        arrayIdx = patternVec[i]  # 1, 2, 3, ...
        intIdx = intIdxList[arrayIdx]  # map to actual intIdx
        p.pos[key][6] = intIdx
        key = p.pos[key][5]
    end
    
    # Initialize patternTails (see Phase 3)
    initPatternTails!(p)
    
    return p
end
```

**Note**: The `buildPatternVector` returns array indices (1, 2, 3, ...) which are then mapped to actual `intIdx` values using the sorted list of dictionary keys. This ensures consistent ordering regardless of dictionary iteration order.

### 2.2 Helper function `buildPatternVector`

Returns `Vector{Int32}` of length `p.length` with monomer types:

```julia
function buildPatternVector(p::Polymer)::Vector{Int32}
    patternVec = zeros(Int32, p.length)
    
    if p.pattern == :uniform
        fill!(patternVec, 1)
        
    elseif p.pattern == :alternating
        for i in 1:p.length
            patternVec[i] = mod1(i, p.nMonomerTypes)
        end
        
    elseif p.pattern == :block
        blockSize = div(p.length, p.nMonomerTypes)
        for i in 1:p.length
            blockIdx = min(div(i - 1, blockSize) + 1, p.nMonomerTypes)
            patternVec[i] = blockIdx
        end
        
    elseif p.pattern == :random
        for i in 1:p.length
            patternVec[i] = rand(1:p.nMonomerTypes)
        end
        
    elseif p.pattern == :custom
        if p.customPattern === nothing
            error("pattern=:custom requires p.customPattern vector")
        end
        if length(p.customPattern) != p.length
            error("customPattern length must equal p.length")
        end
        patternVec .= p.customPattern
        
    else
        error("Unknown pattern: $(p.pattern)")
    end
    
    return patternVec
end
```

**Note**: `:alternating` and `:random` patterns will likely fail the contiguity validation (Phase 2.3) for types >= 2. These patterns may only be valid for testing or special cases where pattern preservation isn't needed.

### 2.3 Validation at initialization

Add a check that each pattern type (2+) appears in exactly one contiguous block:

```julia
function validateSinglePatternBlocks(patternVec::Vector{Int32})
    maxType = maximum(patternVec)
    for t in 2:maxType
        indices = findall(x -> x == t, patternVec)
        if !isempty(indices) && indices != collect(indices[1]:indices[end])
            error("Pattern type $t must appear in exactly one contiguous block")
        end
    end
end
```

### 2.4 Update initialization call order

The initialization order needs to change. Currently in user code:
```julia
p = Polymer(; length=10, maxlength=20, nMonomerTypes=2, ...)
initMonomers!(p)  # ← called BEFORE init!, but needs pos to exist
init!(l, [p])     # ← this creates pos via initSpecies! → initShape!
```

**New order** - user code stays similar, but internal logic changes:

```julia
p = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:custom, customPattern=[...])
initMonomerTypes!(p)  # Creates PolyMonomer objects and matrices (no pos needed)
init!(l, [p])         # initSpecies! → initShape! creates pos, then initPattern! writes types
```

**Changes to `initSpecies!(p::Polymer)`**:

```julia
function initSpecies!(p::Polymer)
    p.name = string("poly", p.id)
    p.pos = OrderedDict{Int32, Vector{Int32}}()
    
    initShape!(p)  # creates 6-element pos vectors
    
    # NEW: If monomers exist, assign pattern types to pos[key][6]
    if p.monomers !== nothing
        initPattern!(p)
    end

    # head and tail beads
    p.head = 1
    # ... rest unchanged ...
end
```

**Backward compatibility**: For existing code that calls `initMonomers!(p)` before `init!`:
- Rename old `initMonomers!` to `initMonomerTypes!`
- Create a new `initMonomers!` wrapper that calls `initMonomerTypes!` (pattern assignment is deferred to `initSpecies!`)

```julia
# Backward-compatible wrapper
function initMonomers!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)
    initMonomerTypes!(p; defaultChemPot=defaultChemPot, defaultIntPot=defaultIntPot)
    # Pattern assignment happens later in initSpecies! when pos exists
    return p
end
```

### 2.5 Tests for Phase 2

1. **Test uniform pattern**: All beads should have type 1, `patternTails` is nothing/empty
2. **Test custom pattern**: `[1,1,1,2,2,2,1,1,1]` → verify types in pos[key][6]
3. **Test promoter-enhancer pattern**: `[1,1,2,2,1,1,3,3,1,1]` → verify types
4. **Test validation rejects non-contiguous**: `[1,2,1,2]` for type 2 should error
5. **Test init order**: Verify `initMonomerTypes!` works before `initShape!`

---

## Phase 3: Pattern Boundary Tracking (`patternTails`)

### 3.1 Add `patternTails` and `baseIntIdx` to Polymer struct

**File**: `src/types.jl`

Add to Polymer:
```julia
patternTails::Union{Nothing, Dict{Int32, Int32}} = nothing  # intIdx -> tail key
baseIntIdx::Union{Nothing, Int32} = nothing                  # intIdx of base type (can grow/shrink)
```

**Note**: `patternTails` is now keyed by `intIdx` (not array index). `baseIntIdx` identifies which monomer type is the "base" that can grow/shrink.

### 3.2 Initialize `patternTails` and `baseIntIdx` in `initPattern!`

The `initPatternTails!` helper (called from `initPattern!`) populates `patternTails`:

```julia
function initPatternTails!(p::Polymer)
    p.patternTails = Dict{Int32, Int32}()
    
    # The base type is the first intIdx in sorted order (or could be specified)
    intIdxList = sort(collect(keys(p.monomers)))
    p.baseIntIdx = intIdxList[1]  # base type = first intIdx
    
    key = p.head
    while key != 0
        currentIntIdx = p.pos[key][6]
        nextKey = p.pos[key][5]
        
        # If not base type, check if this is the tail of a pattern block
        if currentIntIdx != p.baseIntIdx
            if nextKey == 0 || p.pos[nextKey][6] != currentIntIdx
                p.patternTails[currentIntIdx] = key
            end
        end
        key = nextKey
    end
    return p
end
```

### 3.3 Tests for Phase 3

1. **Test patternTails initialization**: For pattern with intIdx values `[A,A,B,B,B,A,A,C,C,A]`:
   - `patternTails[B]` should point to tail of B-block
   - `patternTails[C]` should point to tail of C-block
   - `baseIntIdx == A`
2. **Test patternTails with single pattern**: `[A,A,B,B,A,A]` → only `patternTails[B]` exists

---

## Phase 4: Reject Moves at Type Boundaries

### 4.1 Add boundary check in `move!`

**Goal**: If the selected bead and its next bead have different types, reject the move.

**File**: `src/move.jl`, in `move!(p::Polymer)` after `identifyBond!(p)`:

```julia
# Fast path: homogeneous polymers don't need boundary checks
if p.patternTails === nothing || isempty(p.patternTails)
    # No patterns, skip boundary check
else
    # Reject if bond spans type boundary
    if p.kj != 0 && p.pos[p.ki][6] != p.pos[p.kj][6]
        return false
    end
end
```

**Note**: For homogeneous polymers (`patternTails` is nothing or empty), we skip this check entirely since all beads have type 1.

### 4.2 Tests for Phase 4

1. **Test boundary rejection**: With pattern `[1,1,2,2,1,1]`:
   - Select bond between beads 2-3 (types 1-2) → reject
   - Select bond between beads 4-5 (types 2-1) → reject
   - Select bond between beads 1-2 (types 1-1) → allow (if other conditions pass)
2. **Test statistics**: Run many moves, verify boundary bonds are never modified
3. **Test homogeneous fast path**: Uniform polymer has no boundary rejections

---

## Phase 5: Type Swapping for Pattern Length Preservation

### 5.1 Modify `addBond!` to swap pattern tail types

**Logic**: When adding inside a pattern (type != baseIntIdx):
1. Find the last two beads of that pattern via `patternTails`
2. Change their type to `baseIntIdx`
3. This compensates for the 2 new beads of the pattern type being added
4. Update `patternTails` to point to the new tail (two beads back from the old tail)

The pattern shrinks by 2 at the tail (beads swapped to base) and grows by 2 at the location of bond insertion. Net effect: pattern region "moves" but keeps the same count.

**Key insight for O(1) tail update**: Since `p.pos[key][4]` contains the key of the previous bead, we can find the new tail by moving backwards twice.

**Fast path**: If `currentIntIdx == p.baseIntIdx` (base type), no swapping needed.

```julia
function addBond!(p::Polymer)
    currentIntIdx = p.pos[p.ki][6]
    
    # Capture pattern tail info BEFORE modifying the chain (only for non-base types)
    tailKey = nothing
    prevTailKey = nothing
    newTailKey = nothing
    if currentIntIdx != p.baseIntIdx && p.patternTails !== nothing
        tailKey = p.patternTails[currentIntIdx]
        prevTailKey = p.pos[tailKey][4]
        newTailKey = p.pos[prevTailKey][4]
    end
    
    # Existing add logic (modifies chain)...
    # Set new beads to current type
    p.pos[next_key1][6] = currentIntIdx
    p.pos[next_key2][6] = currentIntIdx
    
    # Also update lattice with intIdx
    p.latS[p.newPosi...] = currentIntIdx
    p.latS[p.newPosj...] = currentIntIdx
    
    # Swap pattern tail beads to base type (only for non-base types)
    if currentIntIdx != p.baseIntIdx && tailKey !== nothing
        # Swap to base type
        p.pos[tailKey][6] = p.baseIntIdx
        p.pos[prevTailKey][6] = p.baseIntIdx
        
        # Update lattice for swapped beads
        p.latS[p.pos[tailKey][1], p.pos[tailKey][2], p.pos[tailKey][3]] = p.baseIntIdx
        p.latS[p.pos[prevTailKey][1], p.pos[prevTailKey][2], p.pos[prevTailKey][3]] = p.baseIntIdx
        
        # Update patternTails
        if newTailKey != 0 && p.pos[newTailKey][6] == currentIntIdx
            p.patternTails[currentIntIdx] = newTailKey
        else
            delete!(p.patternTails, currentIntIdx)
        end
    end
    
    return p
end
```

**Important**: Pattern tail info must be captured **before** modifying the chain.

### 5.2 Modify `removeBond!` to swap base types into pattern

**Logic**: When removing inside a pattern (type != baseIntIdx):
1. Find the current tail and the bead after it
2. Change the bead after the tail (and the one after that) from base type to pattern type
3. This compensates for the 2 pattern beads being removed
4. Update `patternTails` to the new tail (two beads after the old tail)

**Fast path**: If `currentIntIdx == p.baseIntIdx`, no swapping needed.

```julia
function removeBond!(p::Polymer)
    currentIntIdx = p.pos[p.ki][6]
    
    # If removing inside a pattern, extend pattern at tail first
    if currentIntIdx != p.baseIntIdx
        tailKey = p.patternTails[currentIntIdx]
        nextAfterTail = p.pos[tailKey][5]
        nextNextAfterTail = nextAfterTail == 0 ? 0 : p.pos[nextAfterTail][5]
        
        # Validate: there must be at least 2 base-type beads after the tail
        if nextAfterTail == 0 || p.pos[nextAfterTail][6] != p.baseIntIdx ||
           nextNextAfterTail == 0 || p.pos[nextNextAfterTail][6] != p.baseIntIdx
            error("Cannot remove inside pattern: insufficient base-type beads at tail")
        end
        
        # Swap two base-type beads to pattern type
        p.pos[nextAfterTail][6] = currentIntIdx
        p.pos[nextNextAfterTail][6] = currentIntIdx
        
        # Update lattice
        p.latS[p.pos[nextAfterTail][1], p.pos[nextAfterTail][2], p.pos[nextAfterTail][3]] = currentIntIdx
        p.latS[p.pos[nextNextAfterTail][1], p.pos[nextNextAfterTail][2], p.pos[nextNextAfterTail][3]] = currentIntIdx
        
        # Update patternTails: new tail is two beads forward
        p.patternTails[currentIntIdx] = nextNextAfterTail
    end
    
    # Existing remove logic (modifies chain)...
    
    return p
end
```

### 5.3 Validation: Reject moves that would exhaust base-type beads

Before proposing a `removeBond!` move inside a pattern, check that compensation is possible.

**In `move!` or `proposeMove!`**:

```julia
# If removing inside a pattern, check that compensation is possible
if p.moveType == 3 && p.pos[p.ki][6] != p.baseIntIdx
    currentIntIdx = p.pos[p.ki][6]
    tailKey = p.patternTails[currentIntIdx]
    nextAfterTail = p.pos[tailKey][5]
    nextNextAfterTail = nextAfterTail == 0 ? 0 : p.pos[nextAfterTail][5]
    
    # Reject if insufficient base-type beads
    if nextAfterTail == 0 || p.pos[nextAfterTail][6] != p.baseIntIdx ||
       nextNextAfterTail == 0 || p.pos[nextNextAfterTail][6] != p.baseIntIdx
        return false
    end
end

# If adding inside a pattern, check pattern has enough beads to shrink
if p.moveType == 1 && p.pos[p.ki][6] != p.baseIntIdx
    currentIntIdx = p.pos[p.ki][6]
    tailKey = p.patternTails[currentIntIdx]
    prevTailKey = p.pos[tailKey][4]
    newTailKey = p.pos[prevTailKey][4]
    
    # Reject if pattern would be exhausted
    if newTailKey == 0 || p.pos[newTailKey][6] != currentIntIdx
        return false
    end
end
```

### 5.4 Tests for Phase 5

1. **Test addBond inside pattern preserves count**: 
   - Start with `[A,A,B,B,B,B,A,A]` (4 type-B beads, A is base)
   - Add bond inside type-B region
   - Verify still exactly 4 type-B beads
2. **Test removeBond inside pattern preserves count**:
   - Start with `[1,1,2,2,2,2,1,1]` (4 type-2 beads)
   - Remove bond inside type-2 region
   - Verify still exactly 4 type-2 beads
3. **Test patternTails updated correctly after add**: 
   - Verify `patternTails[2]` moves 2 beads backward
4. **Test patternTails updated correctly after remove**:
   - Verify `patternTails[2]` moves 2 beads forward
5. **Test rejection when insufficient base-type beads for remove**:
   - Pattern at end of chain `[1,1,2,2,2,2]` → remove inside type-2 should be rejected
6. **Test rejection when pattern too small for add**:
   - Pattern with only 2 beads `[1,1,2,2,1,1]` → add inside type-2 should be rejected
7. **Test add/remove in base-type region (no swap)**:
   - Add/remove in type-1 region → no patternTails changes, normal behavior

---

## Phase 6: Energy Accounting for Type Swaps

### 6.1 Extend `deltaHam!` to account for swap energy

**Goal**: When add/remove causes type swaps, include the energy change from those swaps.

For **addBond** inside pattern:
- Two beads at tail change from pattern type to `baseIntIdx`
- Energy change: loss of pattern-type interactions, gain of base-type interactions at those positions

For **removeBond** inside pattern:
- Two beads after tail change from `baseIntIdx` to pattern type
- Energy change: loss of base-type interactions, gain of pattern-type interactions at those positions

**Implementation sketch** in `deltaHam!`:

```julia
if p.moveType == 1  # addBond
    currentIntIdx = p.pos[p.ki][6]
    if currentIntIdx != p.baseIntIdx  # non-base type
        # Find the two beads that will be swapped
        tailKey = p.patternTails[currentIntIdx]
        prevTailKey = p.pos[tailKey][4]
        tailPos = p.pos[tailKey][1:3]
        prevTailPos = p.pos[prevTailKey][1:3]
        
        mOld = p.monomers[currentIntIdx]  # monomer being swapped FROM
        mNew = p.monomers[p.baseIntIdx]   # monomer being swapped TO
        
        # For each interacting species, compute energy change
        for (sp, idx) in mOld.intSpecies
            if mOld.intIdx == idx  # self-interaction
                # Use nnp(m, m, pos) to count same-type neighbors
                nnpTail = nnp(mOld, mOld, tailPos)
                nnpPrevTail = nnp(mOld, mOld, prevTailPos)
                p.deltaH -= mOld.intPot * (nnpTail + nnpPrevTail)
                
                nnpTailNew = nnp(mNew, mNew, tailPos)
                nnpPrevTailNew = nnp(mNew, mNew, prevTailPos)
                p.deltaH += mNew.intPot * (nnpTailNew + nnpPrevTailNew)
            else  # cross-species interaction
                nnpTail = nnp(mOld, sp, tailPos)
                nnpPrevTail = nnp(mOld, sp, prevTailPos)
                p.deltaH -= mOld.intMat[mOld.intIdx, idx] * (nnpTail + nnpPrevTail)
                p.deltaH += mNew.intMat[mNew.intIdx, idx] * (nnpTail + nnpPrevTail)
            end
        end
        
        # Chemical potential change
        p.deltaH -= 2 * (mNew.chemPot - mOld.chemPot)
    end
end
```

Similar logic for `moveType == 3` (removeBond), swapping from base to pattern type:

```julia
elseif p.moveType == 3  # removeBond
    currentIntIdx = p.pos[p.ki][6]
    if currentIntIdx != p.baseIntIdx
        tailKey = p.patternTails[currentIntIdx]
        nextAfterTail = p.pos[tailKey][5]
        nextNextAfterTail = p.pos[nextAfterTail][5]
        afterTailPos = p.pos[nextAfterTail][1:3]
        afterAfterTailPos = p.pos[nextNextAfterTail][1:3]
        
        mOld = p.monomers[p.baseIntIdx]    # monomer being swapped FROM (base type)
        mNew = p.monomers[currentIntIdx]   # monomer being swapped TO
        
        for (sp, idx) in mNew.intSpecies
            if mNew.intIdx == idx  # self-interaction
                nnpAfterOld = nnp(mOld, mOld, afterTailPos)
                nnpAfterAfterOld = nnp(mOld, mOld, afterAfterTailPos)
                p.deltaH -= mOld.intPot * (nnpAfterOld + nnpAfterAfterOld)
                
                nnpAfterNew = nnp(mNew, mNew, afterTailPos)
                nnpAfterAfterNew = nnp(mNew, mNew, afterAfterTailPos)
                p.deltaH += mNew.intPot * (nnpAfterNew + nnpAfterAfterNew)
            else  # cross-species interaction
                nnpAfter = nnp(mNew, sp, afterTailPos)
                nnpAfterAfter = nnp(mNew, sp, afterAfterTailPos)
                p.deltaH -= mOld.intMat[mOld.intIdx, idx] * (nnpAfter + nnpAfterAfter)
                p.deltaH += mNew.intMat[mNew.intIdx, idx] * (nnpAfter + nnpAfterAfter)
            end
        end
        
        # Chemical potential change
        p.deltaH -= 2 * (mNew.chemPot - mOld.chemPot)
    end
end
```

### 6.2 Heterogeneous `ham!` and `deltaHam!`

After Phase 0 and 1, `ham!` must iterate per-bead and use `p.pos[key][6]` to get the correct monomer type via dictionary lookup. **Important**: For heterogeneous polymers, self-interaction must use `nnp(m, m, pos)` instead of `nn(p, pos)` to count only same-type neighbors.

```julia
function ham!(p::Polymer)
    p.ham = 0.0
    key = p.head
    for i in 1:p.length
        pos = p.pos[key][1:3]
        intIdx = p.pos[key][6]
        m = p.monomers[intIdx]  # Dictionary lookup by intIdx
        
        # Chemical potential
        p.ham -= m.chemPot
        
        # Interactions
        for (sp, idx) in m.intSpecies
            if m.intIdx == idx  # self-interaction
                # Use nnp(m, m, pos) to count only same-type neighbors
                p.ham += m.intPot * nnp(m, m, pos) / 2.0
            else
                p.ham += m.intMat[m.intIdx, idx] * nnp(m, sp, pos)
            end
        end
        
        key = p.pos[key][5]
    end
    return p
end
```

**Note**: `nnp(m, m, pos)` uses the PolyMonomer dispatch (Phase 0.3) which counts only sites where `latS == m.intIdx`.

For `deltaHam!`, the basic add/kink/remove logic uses `m = p.monomers[p.pos[p.ki][6]]` to get the monomer for beads being moved. Self-interactions use `nnp(m, m, pos)`. The type swap energy (Phase 6.1) is added on top.

### 6.3 Tests for Phase 6

1. **Test deltaHam consistency with swaps**: 
   - Compute `ham!` before
   - Do add inside pattern (which causes swap)
   - Verify `initial + deltaH == ham!` after
2. **Test with different interaction potentials**: 
   - Set type-2 to have different `intPot` than type-1
   - Verify energy change is correct
3. **Test chemical potential accounting**:
   - Set type-2 to have different `chemPot` than type-1
   - Verify swap correctly accounts for chemPot difference
4. **Test kink inside pattern (no swap)**:
   - Kink a bead inside a pattern region
   - Verify correct monomer type is used for energy calculation
   - Verify no type swap occurs

---

## Summary Checklist

| Phase | Description | Files Modified | Key Functions |
|-------|-------------|----------------|---------------|
| 0 | Lattice stores intIdx, `nnp` dispatch for PolyMonomer | utils.jl, init.jl, move.jl | `nnp(::PolyMonomer, ::PolyMonomer, pos)`, `overlap` |
| 1 | 6-element pos vectors, Dict-based monomers, remove monomerIndex | types.jl, init.jl, move.jl, utils.jl, hamiltonian.jl | `initMonomerTypes!`, `getPolyMonomer`, `ham!`, `deltaHam!` |
| 2 | Pattern init with chain traversal, intIdx assignment | init.jl | `initPattern!`, `buildPatternVector`, `validateSinglePatternBlocks` |
| 3 | patternTails + baseIntIdx tracking | types.jl, init.jl | `initPatternTails!` |
| 4 | Boundary rejection | move.jl | `move!` |
| 5 | Type swapping in add/remove using intIdx | move.jl | `addBond!`, `removeBond!` |
| 6 | Energy accounting for swaps | hamiltonian.jl | `deltaHam!` |

### Key Type Changes

| Field | Old Type | New Type |
|-------|----------|----------|
| `Polymer.monomers` | `Vector{PolyMonomer}` | `Dict{Int32, PolyMonomer}` |
| `Polymer.monomerIndex` | `Vector{Int32}` | **Removed** |
| `Polymer.patternTails` | — | `Dict{Int32, Int32}` (intIdx → key) |
| `Polymer.baseIntIdx` | — | `Int32` (intIdx of growable base type) |
| `pos[key][6]` | array index (1, 2, 3) | `intIdx` (globally unique) |
| `latS[x,y,z]` | 0 or 1 | 0 or `intIdx` |

### New/Modified Functions Summary

| Function | File | Purpose |
|----------|------|---------|
| `initMonomerTypes!(p; ...)` | init.jl | Creates `Dict{Int32, PolyMonomer}` keyed by intIdx |
| `initPattern!(p)` | init.jl | Maps array indices to intIdx, writes to `pos[key][6]` |
| `buildPatternVector(p)` | init.jl | Returns vector of array indices (1,2,3...) for pattern |
| `initPatternTails!(p)` | init.jl | Populates `patternTails` Dict and sets `baseIntIdx` |
| `getPolyMonomer(p, key)` | utils.jl | Returns `p.monomers[p.pos[key][6]]` via dict lookup |
| `nnp(::PolyMonomer, ::PolyMonomer, pos)` | utils.jl | Counts neighbors where `latS == m.intIdx` |

---

## Clarifications (from previous review)

1. **What happens when a pattern region exhausts base-type neighbors?** → Reject the move. This is checked in Phase 5.3 before allowing add/remove inside a pattern.

2. **Multiple patterns of different types**: Each non-base pattern type can only appear **once** as a contiguous block. The base type (tracked by `baseIntIdx`) can be interspersed between other types because we don't need to track its tail in `patternTails`.

3. **Kink moves inside patterns**: Kinks don't need pattern-preservation logic since they don't change bead count.

4. **O(1) tail updates**: Use `p.pos[key][4]` (previous pointer) for additions (move 2 beads backward) and `p.pos[key][5]` (next pointer) for removals (move 2 beads forward). No chain traversal needed.

5. **Dictionary vs array indexing**: With the dictionary approach, `p.monomers[intIdx]` gives direct O(1) access. The `intIdx` is stored directly in `pos[key][6]` and `latS[x,y,z]`, eliminating the need for array index → intIdx conversion.

---

## Additional Notes

### Order of operations in add/remove

**addBond!**:
1. Capture pattern tail info (tailKey, prevTailKey, newTailKey) BEFORE modifying chain
2. Create new beads with `currentIntIdx`
3. Update chain pointers
4. If `currentIntIdx != baseIntIdx`: swap tail beads to `baseIntIdx`, update `patternTails`

**removeBond!**:
1. If `currentIntIdx != baseIntIdx`: swap base-type beads after tail to `currentIntIdx`, update `patternTails`
2. Delete beads and update chain pointers

Note: For `removeBond!`, we do the type swap **before** deleting the beads.

### Edge case: Adding at pattern boundary

When adding a bond, the new beads inherit `currentIntIdx = p.pos[p.ki][6]`. Since Phase 4 rejects moves at type boundaries, `ki` and `kj` always have the same type. The two new beads are inserted between them, so they should also have that type.

### Edge case: Pattern at chain ends

If a pattern is at the very end of the chain (tail of chain = tail of pattern), then:
- For `addBond!`: The pattern tail is the chain tail, so `p.pos[tailKey][5] == 0`. After swapping, the new pattern tail moves inward. This is fine.
- For `removeBond!`: We need base-type beads **after** the pattern tail, but there are none. The move should be rejected.

This is correctly handled by the validation in Phase 5.3.

### Edge case: Adding adjacent to the pattern tail

If the bond being added is at the very end of the pattern, the chain structure changes before we can capture tail info. **Solution**: Capture `prevTailKey` and `newTailKey` **before** modifying the chain:

```julia
function addBond!(p::Polymer)
    currentIntIdx = p.pos[p.ki][6]
    
    # Capture pattern tail info BEFORE modifying the chain
    tailKey = nothing
    prevTailKey = nothing
    newTailKey = nothing
    if currentIntIdx != p.baseIntIdx && p.patternTails !== nothing
        tailKey = p.patternTails[currentIntIdx]
        prevTailKey = p.pos[tailKey][4]
        newTailKey = p.pos[prevTailKey][4]
    end
    
    # Existing add logic (modifies chain)...
    # ...
    
    # Now do the swap using the pre-captured keys
    if currentIntIdx != p.baseIntIdx && tailKey !== nothing
        p.pos[tailKey][6] = p.baseIntIdx
        p.pos[prevTailKey][6] = p.baseIntIdx
        # Update lattice too
        p.latS[p.pos[tailKey][1], p.pos[tailKey][2], p.pos[tailKey][3]] = p.baseIntIdx
        p.latS[p.pos[prevTailKey][1], p.pos[prevTailKey][2], p.pos[prevTailKey][3]] = p.baseIntIdx
        
        if newTailKey != 0 && p.pos[newTailKey][6] == currentIntIdx
            p.patternTails[currentIntIdx] = newTailKey
        else
            delete!(p.patternTails, currentIntIdx)
        end
    end
    
    return p
end
```

---

## Remaining Questions for Clarification

~~After reviewing the plan, I have a few questions:~~ **RESOLVED**

1. **patternHeads tracking**: ~~Should we also track `patternHeads`?~~ **Answer: No, only `patternTails` is needed.**

2. **Multiple patterns of different types**: Each type has its own independent entry in `patternTails`. Patterns are independent. **Each pattern type can only appear once contiguously**.

3. **Minimum pattern size**: Pattern size should **not change** during simulation - that's the whole point of type swaps. For initialization, a 0-sized pattern is trivially allowed (meaning that pattern type simply doesn't exist).

4. **Initialization order dependency**: **Solution**: Split `initMonomers!` into:
   - `initMonomerTypes!(p)` - Creates `Dict{Int32, PolyMonomer}` keyed by intIdx. Can be called early.
   - `initPattern!(p)` - Writes intIdx to `pos[key][6]`, initializes `patternTails`. Called **after** `initShape!`.

5. **Backward compatibility**: Fast-paths for homogeneous polymers (`nMonomerTypes = 1` or `pattern = :uniform`). When all beads have the same intIdx, no swapping or `patternTails` tracking is needed.

6. **Why dictionary-based monomers?**: Storing `monomers` as `Dict{Int32, PolyMonomer}` keyed by `intIdx` allows:
   - Direct lookup: `m = p.monomers[pos[key][6]]`
   - No conversion: `pos[key][6]` stores intIdx directly
   - Lattice consistency: Same intIdx in `latS[x,y,z]` and `pos[key][6]`
   - Removes `monomerIndex` array