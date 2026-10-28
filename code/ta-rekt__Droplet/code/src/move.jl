
export  overlap, addBond!, kinkBond!, removeBond!, identifyBond!, 
        makeMove!, move!, identifyBond!

############# Polymer move functions #############
# MARK: -Polymer

function overlap(p::Polymer)
    return ((p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] > 0 && 
             !coords_equal(p.bprev, p.newPosi)) ||
            (p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] > 0 && 
             !coords_equal(p.bnext, p.newPosj)))
end

function addBond!(p::Polymer)

    # TODO: handle corner cases
    @assert (p.newPosi != p.bprev) && (p.newPosj != p.bnext) begin
        "Move is not a bond addition"
    end

    # Inherit intIdx from parent bead (ki) for new beads
    intIdx = p.pos[p.ki][6]
    
    # Capture pattern tail info BEFORE modifying chain (only for non-base types in contiguous pattern region)
    tailKey = nothing
    prevTailKey = nothing
    newTailKey = nothing
    if (p.patternTails !== nothing && intIdx != p.baseIntIdx 
        && haskey(p.patternTails, intIdx) && isInPatternRegion(p, p.ki, intIdx))
        tailKey = p.patternTails[intIdx]
        prevTailKey = p.pos[tailKey][4]
        newTailKey = p.pos[prevTailKey][4]
    end
    
    # set new lattice sites
    p.latS[p.newPosi...] = intIdx
    p.latS[p.newPosj...] = intIdx

    # Find next available keys
    next_key1 = KEYCOUNT[] + 1
    next_key2 = KEYCOUNT[] + 2
    while haskey(p.pos, next_key1) || haskey(p.pos, next_key2)
        KEYCOUNT[] += 1
        next_key1 = KEYCOUNT[] + 1
        next_key2 = KEYCOUNT[] + 2
    end

    # add two new entries for the bond (6-element vectors)
    p.pos[next_key1] = Vector{Int32}(undef, 6)
    p.pos[next_key1][1] = p.newPosi[1]
    p.pos[next_key1][2] = p.newPosi[2]
    p.pos[next_key1][3] = p.newPosi[3]
    p.pos[next_key1][4] = p.ki
    p.pos[next_key1][5] = next_key2
    p.pos[next_key1][6] = intIdx

    p.pos[next_key2] = Vector{Int32}(undef, 6)
    p.pos[next_key2][1] = p.newPosj[1]
    p.pos[next_key2][2] = p.newPosj[2]
    p.pos[next_key2][3] = p.newPosj[3]
    p.pos[next_key2][4] = next_key1
    p.pos[next_key2][5] = p.kj
    p.pos[next_key2][6] = intIdx

    # re-route the old pointers
    p.pos[p.ki][5] = next_key1
    p.pos[p.kj][4] = next_key2

    # Update keys vector
    push!(p.keys_vec, next_key1)
    push!(p.keys_vec, next_key2)

    KEYCOUNT[] += 2
    
    # Swap pattern tail beads to base type (only for non-base types)
    if tailKey !== nothing
        # Swap tail and prev-tail to base type
        p.pos[tailKey][6] = p.baseIntIdx
        p.pos[prevTailKey][6] = p.baseIntIdx
        
        # Update lattice for swapped beads
        p.latS[p.pos[tailKey][1], p.pos[tailKey][2], p.pos[tailKey][3]] = p.baseIntIdx
        p.latS[p.pos[prevTailKey][1], p.pos[prevTailKey][2], p.pos[prevTailKey][3]] = p.baseIntIdx
        
        # Update patternTails
        if newTailKey != 0 && p.pos[newTailKey][6] == intIdx
            p.patternTails[intIdx] = newTailKey
        else
            delete!(p.patternTails, intIdx)
        end
    end

    return p
end

function kinkBond!(p::Polymer)

    # TODO: handle corner cases
    @assert (p.newPosi == p.bprev) ⊻ (p.newPosj == p.bnext) begin
        "Move is not a kink"
    end

    # Preserve the bead's intIdx from pos[key][6]
    if p.newPosi == p.bprev 
        # bead ki moves from bi to newPosj
        intIdx = p.pos[p.ki][6]
        p.latS[p.bi[1], p.bi[2], p.bi[3]] = 0
        p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] = intIdx

        # Update position in-place (no allocation)
        pos_ki = p.pos[p.ki]
        pos_ki[1] = p.newPosj[1]
        pos_ki[2] = p.newPosj[2]
        pos_ki[3] = p.newPosj[3]
        pos_ki[4] = p.kprev
        pos_ki[5] = p.kj
        pos_ki[6] = intIdx

    elseif p.newPosj == p.bnext
        # bead kj moves from bj to newPosi
        intIdx = p.pos[p.kj][6]
        p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] = intIdx
        p.latS[p.bj[1], p.bj[2], p.bj[3]] = 0

        # Update position in-place (no allocation)
        pos_kj = p.pos[p.kj]
        pos_kj[1] = p.newPosi[1]
        pos_kj[2] = p.newPosi[2]
        pos_kj[3] = p.newPosi[3]
        pos_kj[4] = p.ki
        pos_kj[5] = p.knext
        pos_kj[6] = intIdx
    end

    return p
end

function removeBond!(p::Polymer)

    # TODO: handle corner cases
    @assert (p.newPosi == p.bprev) && (p.newPosj == p.bnext) begin
        "Move is not a bond removal"
    end
    
    # Get the intIdx before removing (both beads have same type due to boundary check)
    intIdx = p.pos[p.ki][6]
    
    # Extend pattern at tail BEFORE removing beads (only for non-base types in contiguous pattern region)
    if p.patternTails !== nothing && intIdx != p.baseIntIdx && haskey(p.patternTails, intIdx) && isInPatternRegion(p, p.ki, intIdx)
        tailKey = p.patternTails[intIdx]
        nextAfterTail = p.pos[tailKey][5]
        nextNextAfterTail = nextAfterTail == 0 ? Int32(0) : p.pos[nextAfterTail][5]
        
        # Swap two base-type beads to pattern type
        p.pos[nextAfterTail][6] = intIdx
        p.pos[nextNextAfterTail][6] = intIdx
        
        # Update lattice
        p.latS[p.pos[nextAfterTail][1], p.pos[nextAfterTail][2], p.pos[nextAfterTail][3]] = intIdx
        p.latS[p.pos[nextNextAfterTail][1], p.pos[nextNextAfterTail][2], p.pos[nextNextAfterTail][3]] = intIdx
        
        # Update patternTails: new tail is two beads forward
        p.patternTails[intIdx] = nextNextAfterTail
    end
    
    # reset old lattice sites to 0
    p.latS[p.bi...] = 0
    p.latS[p.bj...] = 0

    # join ends together - handle boundary cases
    if p.kprev != 0
        p.pos[p.kprev][5] = p.knext
    else
        # Removing from the start, update head
        p.head = p.knext
    end
    
    if p.knext != 0
        p.pos[p.knext][4] = p.kprev
    end

    # Update keys vector - remove both keys
    filter!(k -> k != p.ki && k != p.kj, p.keys_vec)

    # remove keys associated with old bond
    delete!(p.pos, p.ki)
    delete!(p.pos, p.kj)

    return p
end

# Finds the bead coordinates, neighbording beads, and vector of the bond
# to which the bead selected by p.moveIdx (or the bead key passed via `ki`)
# belongs. Pass `ki` directly to bypass the keys_vec[moveIdx] lookup —
# necessary for KMC where the slot-encoded bondIdx maps to a stable bead key.
function identifyBond!(p::Polymer; ki::Union{Nothing, Int32}=nothing)

    if ki === nothing
        @assert p.moveIdx >= 1 && p.moveIdx <= p.length begin
            "p.moveIdx must be between 1 and Species length (# of monomers)"
        end
        ki = p.keys_vec[p.moveIdx]
    end

    # Initialize bond coordinate buffers if needed (lazy initialization for restored/manually-created polymers)
    if p.bi === nothing
        p.bi = zeros(Int32, 3); p.bj = zeros(Int32, 3)
        p.bprev = zeros(Int32, 3); p.bnext = zeros(Int32, 3)
    end

    # pointers to the two beads (i and j) to be moved
    p.ki = ki
    p.kj = p.pos[p.ki][5]

    # previous and next beads along the chain
    p.kprev = p.pos[p.ki][4]
    p.knext = p.kj == 0 ? 0 : p.pos[p.kj][5]

    # bond edge coordinates - copy values instead of creating views
    pos_ki = p.pos[p.ki]
    p.bi[1] = pos_ki[1]; p.bi[2] = pos_ki[2]; p.bi[3] = pos_ki[3]
    
    if p.kj == 0
        p.bj[1] = SENTINEL[1]; p.bj[2] = SENTINEL[2]; p.bj[3] = SENTINEL[3]
    else
        pos_kj = p.pos[p.kj]
        p.bj[1] = pos_kj[1]; p.bj[2] = pos_kj[2]; p.bj[3] = pos_kj[3]
    end

    # coordinates of previous and next beads relative to bond
    if p.kprev == 0
        p.bprev[1] = SENTINEL[1]; p.bprev[2] = SENTINEL[2]; p.bprev[3] = SENTINEL[3]
    else
        pos_kprev = p.pos[p.kprev]
        p.bprev[1] = pos_kprev[1]; p.bprev[2] = pos_kprev[2]; p.bprev[3] = pos_kprev[3]
    end
    
    if p.knext == 0
        p.bnext[1] = SENTINEL[1]; p.bnext[2] = SENTINEL[2]; p.bnext[3] = SENTINEL[3]
    else
        pos_knext = p.pos[p.knext]
        p.bnext[1] = pos_knext[1]; p.bnext[2] = pos_knext[2]; p.bnext[3] = pos_knext[3]
    end

    for i in range(1, 3)
        p.bondVec[i] = (p.bj[i] - p.bi[i]) % p.boxlength
    end

    return p
end

# Calls the appropriate operation (add, kink or remove) to execute a move
function makeMove!(p::Polymer)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end

    # bond addition
    if p.moveType == 1
        addBond!(p)

    # kink
    elseif p.moveType == 2
        kinkBond!(p)

    # bond removal  
    elseif p.moveType == 3 
        removeBond!(p)
    end

    return p
end

# executes a move on polymermer s inside lattice l at position p.moveIdx, of type
# moveType (add, remove or kink) sampled according to the metropolis algorithm
function move!(p::Polymer)

    # sample a bond at random along the chain
    p.moveIdx = rand(1:p.length - 1)
    
    # find coordinates and keys of selected bond and adjacent beads
    identifyBond!(p)
    
    # don't move if last bead was selected
    if p.kj != 0
        
        # Reject if bond spans type boundary (heterogeneous polymers only)
        if p.patternTails !== nothing && !isempty(p.patternTails)
            if p.pos[p.ki][6] != p.pos[p.kj][6]
                return false
            end
        end
        
        # sample a move according to proposition probabilities
        proposeMove!(p)
        
        # Validate pattern swap is possible for non-base types
        if p.patternTails !== nothing && !isempty(p.patternTails)
            currentIntIdx = p.pos[p.ki][6]
            if currentIntIdx != p.baseIntIdx && haskey(p.patternTails, currentIntIdx)
                tailKey = p.patternTails[currentIntIdx]
                prevTailKey = p.pos[tailKey][4]
                
                if p.moveType == 1  # addBond: pattern shrinks at tail
                    # Need at least 4 beads in pattern (2 will be swapped)
                    newTailKey = p.pos[prevTailKey][4]
                    if newTailKey == 0 || p.pos[newTailKey][6] != currentIntIdx
                        return false  # Pattern too small
                    end
                    
                    # CRITICAL: Reject if bond includes tail or prevTail
                    # Adding at tail/prevTail would create fragmentation
                    if p.ki == tailKey || p.kj == tailKey || p.ki == prevTailKey || p.kj == prevTailKey
                        return false  # Cannot add at tail position
                    end
                    
                elseif p.moveType == 3  # removeBond: need base-type beads after tail
                    nextAfterTail = p.pos[tailKey][5]
                    nextNextAfterTail = nextAfterTail == 0 ? Int32(0) : p.pos[nextAfterTail][5]
                    if nextAfterTail == 0 || p.pos[nextAfterTail][6] != p.baseIntIdx ||
                       nextNextAfterTail == 0 || p.pos[nextNextAfterTail][6] != p.baseIntIdx
                        return false  # Insufficient base-type beads
                    end
                end
            end
        end
        
        # check for max length constraint and overlap
        # Note: addBond (moveType==1) adds 2 beads, so reject when length >= maxlength - 1
        if !overlap(p) && !(p.moveType == 1 && p.length >= p.maxlength - 1)
            
            # accept move with metropolis criterion
            if acceptMove(p)
                
                makeMove!(p)
                p.length = length(keys(p.pos))   # TODO: refactor for speed
                return true
            end
        end
    end

    return false
end

# delegates moving to parent polymer
move!(m::PolyMonomer) = move!(m.poly)


############# BulkMonomer move functions #############
# MARK: -BulkMonomer

function overlap(b::AbstractBulkMonomer)
    return b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] > 0
end

function identifyBond!(b::BulkMonomer)
    b.bi = view(b.pos[b.moveIdx], range(1, 3))

    return b
end

# Translates a particle to a new position set by b.newPos
function makeMove!(b::BulkMonomer)
    
    # update lattice sites
    b.latS[b.bi[1], b.bi[2], b.bi[3]] = 0
    b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] = b.moveIdx

    # update current position
    for i in range(1, 3)
        b.pos[b.moveIdx][i] = b.newPos[i]
    end

    @assert b.latS[b.pos[b.moveIdx][1], 
                   b.pos[b.moveIdx][2], 
                   b.pos[b.moveIdx][3]] == b.moveIdx begin
        "Invalid index-position correspondence"
    end

    return b
end

# Executes a move on a species of bulkMonomer monomers diluted in the box
function move!(b::BulkMonomer)

    # sample a particle at random
    b.moveIdx = rand(1:b.length)

    # find coordinates and keys of selected bond
    identifyBond!(b)
    
    # sample a move uniformly in each direction
    proposeMove!(b)

    if !overlap(b)
        if acceptMove(b)
            makeMove!(b)
            return true
        end
    end

    return false
end


############# FastBulkMonomer move functions #############
# MARK: -FastBulkMonomer

function identifyBond!(b::FastBulkMonomer)
    b.moveIdx = b.latS[b.bi[1], b.bi[2], b.bi[3]]
    b.moveIdxNew = b.latS[b.newPos[1], b.newPos[2], b.newPos[3]]

    return b
end

# Translates a particle to a new position set by b.newPos
function makeMove!(b::FastBulkMonomer)
    
    # update lattice sites
    b.latS[b.bi[1], b.bi[2], b.bi[3]] = b.moveIdxNew
    b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] = b.moveIdx

    # swap old with new position
    if b.moveIdx != 0
        for i in range(1, 3)
            b.pos[b.moveIdx][i] = b.newPos[i]
        end

        @assert b.latS[b.pos[b.moveIdx][1], 
                       b.pos[b.moveIdx][2], 
                       b.pos[b.moveIdx][3]] == b.moveIdx begin
            "Invalid index-position correspondence"
        end
    end

    # swap new with old position
    if b.moveIdxNew != 0
        for i in range(1, 3)
            b.pos[b.moveIdxNew][i] = b.bi[i]
        end

        @assert b.latS[b.pos[b.moveIdxNew][1], 
                       b.pos[b.moveIdxNew][2], 
                       b.pos[b.moveIdxNew][3]] == b.moveIdxNew begin
            "Invalid index-position correspondence"
        end
    end

    return b
end

# Executes a move on a species of bulkMonomer monomers diluted in the box
function move!(b::FastBulkMonomer)

    # pick two sites on the lattice uniformly without replacement
    proposeMove!(b)
    
    if xor(b.moveIdx == 0, b.moveIdxNew == 0) && acceptMove(b)
        makeMove!(b)
        return true
    end

    return false
end