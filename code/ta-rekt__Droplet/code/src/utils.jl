export coords_equal, getLatS, getTempVec, getBoxlength, getPolyMonomer, 
       getMonomerPotentials, boxmod!, vol, onEdge, nn, nnp, nlp, cplpq, nnIdxs!, splitDroplet!,
       dropletBbox, fillBbox!, clearBbox!, getHam, getDeltaHam, getDeltaHP,
       setHam, setDeltaHam, setDeltaHP, isneighbor, isInPatternRegion

############# Generic utilities #############
# MARK: -Generic

# Check if two 3D positions are neighbors (Manhattan distance == 1)
@inline function isneighbor(pos1::Union{SubArray{Int32}, Vector{Int32}}, 
                            pos2::Union{SubArray{Int32}, Vector{Int32}})
    return abs(pos1[1] - pos2[1]) + abs(pos1[2] - pos2[2]) + abs(pos1[3] - pos2[3]) == 1
end

# Check if bead at key is in the contiguous pattern region (connected to pattern tail)
# Returns true if walking forward from key while staying in the same type reaches the tail
function isInPatternRegion(p::Polymer, key::Integer, intIdx::Integer)
    tailKey = p.patternTails[intIdx]
    
    # Walk forward from key toward tail, checking if we reach the tail
    currentKey = key
    while currentKey != 0
        if currentKey == tailKey
            return true
        end
        # If we hit a different type, we're not in the contiguous region
        if p.pos[currentKey][6] != intIdx
            return false
        end
        currentKey = p.pos[currentKey][5]  # next bead (forward)
    end
    
    return false
end

@inline function coords_equal(a, b)
    return a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
end


# Helper functions to access fields that may be on the object itself
# or on a parent (for PolyMonomer)
@inline function getLatS(m::AbstractPolyMonomer)
    return m.poly.latS
end

@inline function getLatS(s::AbstractSpecies)
    return s.latS
end

@inline function getTempVec(m::AbstractPolyMonomer)
    return m.poly.tempVec
end

@inline function getTempVec(s::AbstractSpecies)
    return s.tempVec
end

@inline function getBoxlength(m::AbstractPolyMonomer)
    return m.poly.boxlength
end

@inline function getBoxlength(s::AbstractSpecies)
    return s.boxlength
end

@inline function getHam(m::AbstractPolyMonomer)
    return m.poly.ham
end

@inline function getHam(s::AbstractSpecies)
    return s.ham
end

@inline function getDeltaHam(m::AbstractPolyMonomer)
    return m.poly.deltaH
end

@inline function getDeltaHam(s::AbstractSpecies)
    return s.deltaH
end

@inline function getDeltaHP(m::AbstractPolyMonomer)
    return m.poly.deltaHP
end

@inline function getDeltaHP(s::AbstractSpecies)
    return s.deltaHP
end

@inline function setHam(m::AbstractPolyMonomer, val::Float64)
    m.poly.ham = val
    return m
end

@inline function setHam(s::AbstractSpecies, val::Float64)
    s.ham = val
    return s
end

@inline function setDeltaHam(m::AbstractPolyMonomer, val::Float64)
    m.poly.deltaH = val
    return m
end

@inline function setDeltaHam(s::AbstractSpecies, val::Float64)
    s.deltaH = val
    return s
end

@inline function setDeltaHP(m::AbstractPolyMonomer, val::Float64)
    # Use += to accumulate contributions from all monomer types
    # When a bulk species moves, it may interact with multiple monomer types
    # Each contributing to the polymer's total deltaHP
    m.poly.deltaHP += val
    return m
end

@inline function setDeltaHP(s::AbstractSpecies, val::Float64)
    s.deltaHP = val
    return s
end


"""
    getPolyMonomer(p::Polymer, key::Int32)

Get the PolyMonomer for a given bead by its key in the position dictionary.
Uses pos[key][6] which stores the intIdx directly for O(1) dictionary lookup.
"""
function getPolyMonomer(p::Polymer, key::Int32)
    if p.monomers === nothing
        return nothing
    end
    intIdx = p.pos[key][6]
    return p.monomers[intIdx]
end


@inline function getMonomerPotentials(poly::AbstractPolymer, key::Int32)
    monomer = getPolyMonomer(poly, key)
    if monomer === nothing
        return (nothing, nothing)
    end
    return (monomer.chemPot, monomer.intPot)
end

@inline function getMonomerPotentials(s::AbstractSpecies, beadIdx::Int32)
    # For all other species, return the species-level potentials
    return (s.chemPot, s.intPot)
end


# ensures positions are periodic within the lattice's bounding box
function boxmod!(v::Union{SubArray{Int32}, Vector{Int32}}, 
                 w::Union{SubArray{Int32}, Vector{Int32}}, 
                 lb::Int32, ub::Int32)

    for i in eachindex(v)
        w[i] = if v[i] > ub
            v[i] - (ub - lb + 1)
        elseif v[i] < lb
            v[i] + (ub - lb + 1)
        else
            v[i]
        end
    end
    
    return w
end


# volume of a bounding box
function vol(bbox::Vector{Vector{Int32}})

    @assert length(bbox) == 3
    @assert all(length(b) == 2 for b in bbox)

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    
    return abs(xh - xl + 1) * abs(yh - yl + 1) * abs(zh - zl + 1)
end


# returns true if at least one element of x is adjacent to the bbox
function onEdge(s::AbstractSpecies, x::Set{Int32}, bbox::Vector{Vector{Int32}})

    @assert length(bbox) == 3
    @assert all(length(b) == 2 for b in bbox)

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    
    for p in x
        pos = s.pos[p]
        if (pos[1] == xl || pos[1] == xh ||
            pos[2] == yl || pos[2] == yh ||
            pos[3] == zl || pos[3] == zh)
            return true
        end
    end

    return false
end


# Number of same-species nearest neighbors at a given position
function nn(s::AbstractSpecies, pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    latS = getLatS(s)
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (latS[tempVec[1], tempVec[2], tempVec[3]] > 0)

        pos[i] -= 2
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (latS[tempVec[1], tempVec[2], tempVec[3]] > 0)

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Number of nearest neighbors between species s and p at a given position
# # "nnp" = number of nearest partners
function nnp(s::AbstractSpecies, 
             p::AbstractSpecies, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(p)

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] > 0)

        pos[i] -= 2
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] > 0)

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Number of nearest neighbors for PolyMonomer-to-PolyMonomer interactions
# Counts only sites occupied by target monomer type (matching intIdx)
# # "nnp" = number of nearest partners
function nnp(s::AbstractPolyMonomer, 
             target::AbstractPolyMonomer, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(target)

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == target.intIdx)

        pos[i] -= 2
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == target.intIdx)

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Number of nearest neighbors for BulkMonomer-to-PolyMonomer interactions
# Counts only sites occupied by target monomer type (matching intIdx)
# # "nnp" = number of nearest partners
function nnp(s::AbstractBulkMonomer, 
             target::AbstractPolyMonomer, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(target)

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == target.intIdx)

        pos[i] -= 2
        boxmod!(pos, tempVec, Int32(1), boxlength)
        n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == target.intIdx)

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Number of local neighbors for BulkMonomer-to-PolyMonomer interactions
# Counts only sites occupied by target monomer type (matching intIdx)
# "nlp" = number of local partners
function nlp(s::AbstractBulkMonomer, 
             target::AbstractPolyMonomer, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(target)

    for dx in Int32(-1):Int32(1), dy in Int32(-1):Int32(1), dz in Int32(-1):Int32(1)
        if !(dx == 0 && dy == 0 && dz == 0)
            tempVec[1] = pos[1] + dx
            tempVec[2] = pos[2] + dy
            tempVec[3] = pos[3] + dz
            boxmod!(tempVec, tempVec, Int32(1), boxlength)
            n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == target.intIdx)
        end
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Number of local neighbors for PolyMonomer-to-PolyMonomer interactions
# Counts only sites occupied by target monomer type (matching intIdx)
# "nlp" = number of local partners
function nlp(s::AbstractPolyMonomer, 
             target::AbstractPolyMonomer, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(target)

    for dx in Int32(-1):Int32(1), dy in Int32(-1):Int32(1), dz in Int32(-1):Int32(1)
        if !(dx == 0 && dy == 0 && dz == 0)
            tempVec[1] = pos[1] + dx
            tempVec[2] = pos[2] + dy
            tempVec[3] = pos[3] + dz
            boxmod!(tempVec, tempVec, Int32(1), boxlength)
            n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] == target.intIdx)
        end
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Number of local neighbors for PolyMonomer-to-BulkMonomer interactions
# Counts only sites occupied by target species
# "nlp" = number of local partners
function nlp(s::AbstractPolyMonomer, 
             target::AbstractBulkMonomer, 
             pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    tempVec = getTempVec(s)
    boxlength = getBoxlength(s)
    pLatS = getLatS(target)

    for dx in Int32(-1):Int32(1), dy in Int32(-1):Int32(1), dz in Int32(-1):Int32(1)
        if !(dx == 0 && dy == 0 && dz == 0)
            tempVec[1] = pos[1] + dx
            tempVec[2] = pos[2] + dy
            tempVec[3] = pos[3] + dz
            boxmod!(tempVec, tempVec, Int32(1), boxlength)
            n += (pLatS[tempVec[1], tempVec[2], tempVec[3]] > 0)
        end
    end

    # clear temporary containers
    clamp!(tempVec, 0, 0)

    return n
end


# Boolean: do species p and q both occupy lattice site pos?
# Specialized methods handle PolyMonomer types by checking intIdx match.

# Generic version: both p and q use occupancy check (> 0)
function cplpq(p::AbstractSpecies, 
               q::AbstractSpecies, 
               pos::Union{SubArray{Int32}, Vector{Int32}})
    pLatS = getLatS(p)
    qLatS = getLatS(q)
    return (pLatS[pos[1], pos[2], pos[3]] > 0) * (qLatS[pos[1], pos[2], pos[3]] > 0)
end

# PolyMonomer p: check p.intIdx match
function cplpq(p::AbstractPolyMonomer, 
               q::AbstractSpecies, 
               pos::Union{SubArray{Int32}, Vector{Int32}})
    pLatS = getLatS(p)
    qLatS = getLatS(q)
    return (pLatS[pos[1], pos[2], pos[3]] == p.intIdx) * (qLatS[pos[1], pos[2], pos[3]] > 0)
end

# PolyMonomer q: check q.intIdx match
function cplpq(p::AbstractSpecies, 
               q::AbstractPolyMonomer, 
               pos::Union{SubArray{Int32}, Vector{Int32}})
    pLatS = getLatS(p)
    qLatS = getLatS(q)
    return (pLatS[pos[1], pos[2], pos[3]] > 0) * (qLatS[pos[1], pos[2], pos[3]] == q.intIdx)
end

# Both PolyMonomer: check intIdx match for both
function cplpq(p::AbstractPolyMonomer, 
               q::AbstractPolyMonomer, 
               pos::Union{SubArray{Int32}, Vector{Int32}})
    pLatS = getLatS(p)
    qLatS = getLatS(q)
    return (pLatS[pos[1], pos[2], pos[3]] == p.intIdx) * (qLatS[pos[1], pos[2], pos[3]] == q.intIdx)
end


# Set of nearest neighbors of a given position
function nnIdxs!(b::AbstractSpecies, 
                 pos::Union{SubArray{Int32}, Vector{Int32}})
    k = 1

    # clear working containers
    clamp!(b.tempVec, 0, 0)
    clamp!(b.nnIdxs, 0, 0)
    empty!(b.freeNNs)

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, b.tempVec, Int32(1), b.boxlength)
        b.nnIdxs[k] = b.latS[b.tempVec[1], b.tempVec[2], b.tempVec[3]]
        
        if b.nnIdxs[k] == 0
            push!(b.freeNNs, k)
        end
        k += 1

        pos[i] -= 2
        boxmod!(pos, b.tempVec, Int32(1), b.boxlength)
        b.nnIdxs[k] = b.latS[b.tempVec[1], b.tempVec[2], b.tempVec[3]]
        
        if b.nnIdxs[k] == 0
            push!(b.freeNNs, k)
        end
        k += 1

        pos[i] += 1
    end

    # clear temporary containers
    clamp!(b.tempVec, 0, 0)

    return b
end


############# BulkMonomer utilities #############
# MARK: -Fast/BulkMonomer

# splits droplets that sit across a lattice boundary
# WARNING: assumes droplet max diameter < boxlength / 2
function splitDroplet!(s::AbstractBulkMonomer, droplet::Set{Int32})

    droplet_vec = collect(droplet)
    droplet_vecs = [[droplet_vec[1]]]
    half_boxlength = s.boxlength / 2
    
    for p in droplet_vec[2:end]
        belonged = false
        for d in droplet_vecs
            ref_pos = s.pos[d[1]]
            p_pos = s.pos[p]
            if (abs(p_pos[1] - ref_pos[1]) < half_boxlength &&
                abs(p_pos[2] - ref_pos[2]) < half_boxlength &&
                abs(p_pos[3] - ref_pos[3]) < half_boxlength)
                push!(d, p)
                belonged = true
                break
            end
        end

        if !belonged
            push!(droplet_vecs, [p])
        end
    end

    return [Set(d) for d in droplet_vecs]
end


# finds the bounding box of a droplet, adds 1 along each axis for the empty
# space around the droplet to be connected
function dropletBbox(s::AbstractBulkMonomer, droplet::Set{Int32})

    # start with position of an arbitrary particle in the droplet
    key0 = first(droplet)
    pos0 = s.pos[key0]
    xb = [pos0[1], pos0[1]]
    yb = [pos0[2], pos0[2]]
    zb = [pos0[3], pos0[3]]

    for key in droplet
        pos = s.pos[key]
        # compute bounds
        xb[1] = min(pos[1], xb[1])
        xb[2] = max(pos[1], xb[2])
        yb[1] = min(pos[2], yb[1])
        yb[2] = max(pos[2], yb[2])
        zb[1] = min(pos[3], zb[1])
        zb[2] = max(pos[3], zb[2])
    end

    # add 1 and clip to lattice bounds
    xb[1] = max(xb[1] - 1, 1)
    xb[2] = min(xb[2] + 1, s.boxlength)
    yb[1] = max(yb[1] - 1, 1)
    yb[2] = min(yb[2] + 1, s.boxlength)
    zb[1] = max(zb[1] - 1, 1)
    zb[2] = min(zb[2] + 1, s.boxlength)

    return [xb, yb, zb]
end


# fills up a box around a droplet in the co-lattice
function fillBbox!(s::AbstractBulkMonomer, droplet::Set{Int32}, bbox::Vector{Vector{Int32}})

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]

    # ghost species
    gs = s.ghost
    key = 0

    # loop across the bbox, if site not part of the droplet, add it to the ghost lattice
    for x in xl:xh, y in yl:yh, z in zl:zh
        if !(s.latS[x, y, z] in droplet)
            key += 1
            gs.latS[x, y, z] = key
            gs.pos[key] = [x, y, z, key - 1, key + 1]
        end
    end

    gs.pos[1][4] = 0
    gs.pos[key][5] = 0

    return s
end


# clears a box around a droplet in the ghost lattice
function clearBbox!(s::AbstractBulkMonomer, bbox::Vector{Vector{Int32}})

    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]

    # ghost species
    gs = s.ghost

    # loop across x,y,z in the bbox and clear ghost lattice
    for x in xl:xh, y in yl:yh, z in zl:zh
        gs.latS[x, y, z] = 0
    end

    gs.pos = OrderedDict{Int32, Vector{Int32}}()

    return s
end