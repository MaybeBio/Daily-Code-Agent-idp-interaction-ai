
export acceptMove, moveTypes!, proposeMove!

############# Generic metropolis functions #############
# MARK: -Generic

# accepts or rejects according to the metropolis criterion
function acceptMove(s::AbstractSpecies)

    deltaHam!(s)
    prob = min(1, exp(-s.deltaH))

    return rand() < prob
end


############# Polymer metropolis functions #############
# MARK: -Polymer

# Computes all possible move types for the selected bond
# Warning: needs p.bondVec to be set
function moveTypes!(p::Polymer)

    # dimensions along which to move
    dims = findall(iszero, p.bondVec)

    @assert length(dims) == 2 begin
        "bondVec must be aligned with x, y, or z axis"
    end

    # build all possible moves each bond can make
    for (i, idx, sign) in zip(range(1, 4), COORD_MOVES, MOVE_ORIENTATIONS)
        
        for i in range(1, 3)
            p.moveVec[i] = 0
        end

        # flip direction along each dimension in the orthogonal plane
        p.moveVec[dims[idx]] = sign

        # Compute bi + b.moveVec
        for j in range(1, 3)
            p.tempVec[j] = p.bi[j] + p.moveVec[j]
        end
        boxmod!(p.tempVec, p.newPosi, Int32(1), p.boxlength)
        
        # Compute bj + b.moveVec
        for j in range(1, 3)
            p.tempVec[j] = p.bj[j] + p.moveVec[j]
        end
        boxmod!(p.tempVec, p.newPosj, Int32(1), p.boxlength)
        
        
        # 1 -> bond addition, 2 -> kink, 3 -> bond removal
        p.moveTypes[i] = (Int32(p.newPosi == p.bprev) 
                        + Int32(p.newPosj == p.bnext) + 1)
    end

    for i in range(1, 3)
        p.moveVec[i] = 0
    end

    return p
end

# samples moves according to probabilities set by detailed balance
# writes the sampled direction vector in p.moveVec
function proposeMove!(p::Polymer)
    
    # compute the possible move types around a bond
    moveTypes!(p)

    # compute how many of each move type there are
    for j in range(1, 3)
        p.nTypes[j] = 0
    end

    for i in range(1, 4)
        p.nTypes[p.moveTypes[i]] += 1
    end
    
    N = p.length
    
    # set Pkink to zero if no kinks are possible
    Pkink = p.nTypes[2] == 0 ? 0 : p.moveTypeProbs[2]
    
    # Padd
    p.moveTypeProbs[1] = p.nTypes[3] == 0 ? (1 - Pkink) : (1 - Pkink) * (N + 2) / (2*N + 2)
    
    # Premove
    p.moveTypeProbs[3] = p.nTypes[3] == 0 ? 0 : (1 - Pkink) * N / (2*N + 2)

    @assert p.moveTypeProbs[1] + p.moveTypeProbs[3] + Pkink == 1 begin
        "Proposition probabilities don't sum to 1"
    end

    dims = findall(iszero, p.bondVec)

    # sample move types according to proposed probabilities
    r = rand()
    prob = 0.0

    for (i, idx, sign) in zip(range(1, 4), COORD_MOVES, MOVE_ORIENTATIONS)
        
        for j in range(1, 3)
            p.moveVec[j] = 0
        end
        
        typ = p.moveTypes[i]
        prob += p.moveTypeProbs[typ] / p.nTypes[typ]

        if r < prob
            p.moveVec[dims[idx]] = sign
            p.moveType = typ
            break
        end
    end

    @assert dot(p.bondVec, p.moveVec) == 0 begin
        "Move is not orthogonal to bond direction"
    end

    # new position is old position translated in proposed direction
    for j in range(1, 3)
        p.tempVec[j] = p.bi[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosi, Int32(1), p.boxlength)
    
    for j in range(1, 3)
        p.tempVec[j] = p.bj[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosj, Int32(1), p.boxlength)

    return p
end


############# BulkMonomer metropolis functions #############
# MARK: -BulkMonomer

function proposeMove!(b::BulkMonomer)

    for i in range(1, 3)
        b.tempVec[i] = b.pos[b.moveIdx][i]
    end

    # sample an x/y/z direction and orientation at random
    b.moveType = rand(1:3)
    b.tempVec[b.moveType] += 2 * (rand() < 0.5) - 1
    
    # clip to within bounding box
    boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)

    return b
end


############# FastBulkMonomer metropolis functions #############
# MARK: -FastBulkMonomer

function proposeMove!(b::FastBulkMonomer)

    b.moveVec = sample(0:(b.boxlength^3-1), 2, replace=false)

    mx = mod(b.moveVec[1], b.boxlength^2)
    b.bi[1] = div(b.moveVec[1], b.boxlength^2) + 1 # x coordinate
    b.bi[2] = div(mx, b.boxlength) + 1            # y coordinate
    b.bi[3] = mod(mx, b.boxlength) + 1            # z coordinate

    mx = mod(b.moveVec[2], b.boxlength^2)
    b.newPos[1] = div(b.moveVec[2], b.boxlength^2) + 1 # x coordinate
    b.newPos[2] = div(mx, b.boxlength) + 1            # y coordinate
    b.newPos[3] = mod(mx, b.boxlength) + 1            # z coordinate

    identifyBond!(b)

    # resample if positions are in each other's neighborhood
    nnIdxs!(b, b.bi)

    if b.moveIdxNew in b.nnIdxs
        proposeMove!(b)
    end

    return b
end