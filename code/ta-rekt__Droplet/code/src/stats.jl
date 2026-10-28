export extractCoords!, centerMass, centerMassPattern, radius, radiusLine, radiusSphere, boxDensity,
       axialDensity , radialDensity, extractCoords!, dfs, cc, dropletVolume,
       dropletDensities, diluteDensity, nspeciesBound, get_pattern_centers,
       patterns_in_contact, compute_contact_matrix!, find_patterns, updateDensityProfile!,
       updateNormalHistogram!, contactTimes, distanceHistogram, beadMSD, beadPairDistanceMSD

############### Generic #############
# MARK: -Generic

# extracts the coordinates of a polymer's beads
function extractCoords!(b::AbstractSpecies, coords::Matrix)

    for (i, k) in enumerate(collect(keys(b.pos)))
        coords[i, 1:3] = b.pos[k][1:3]         # next position
        coords[i, 4] = b.latS[b.pos[k][1], b.pos[k][2], b.pos[k][3]]   # bead index
    end

    return coords
end

# returns the center of mass of all (or a subset of) particles of a species
function centerMass(s::AbstractSpecies, 
                    subset::Union{Nothing, Set{Int32}} = nothing)
    cm = zeros(3)
    subset = subset === nothing ? eachindex(s.pos) : subset

    # compute average of all positions in pos
    for i in subset
        cm = cm .+ view(s.pos[i], 1:3)
    end

    return cm ./ length(subset)
end

# returns the center of mass of a pattern along a polymer
function centerMassPattern(poly::AbstractPolymer, patternIntIdx::Int32)
    cm = zeros(3)
    count = 0

    if haskey(poly.patternTails, patternIntIdx)
        
        # walk backwards from tail 
        tail_key = poly.patternTails[patternIntIdx]
        current_key = tail_key

        while current_key != 0 && poly.pos[current_key][6] == patternIntIdx
            cm = cm .+ view(poly.pos[current_key], 1:3)
            count += 1
            # move to previous bead
            current_key = poly.pos[current_key][4]
        end
    end

    return count > 0 ? cm ./ count : cm
end


# radius of gyration, mean squared distance between each bead and
# the center of mass
function radius(s::AbstractSpecies)
    cm = centerMass(s)
    rg = 0

    for i in eachindex(s.pos)
        rg += ((cm[1] .- s.pos[i][1])^2
             + (cm[2] .- s.pos[i][2])^2
             + (cm[3] .- s.pos[i][3])^2)
    end

    return sqrt(rg / s.length)
end


# gyration radius of a line
function radiusLine(length)
    n = div(length, 2)
    return sqrt((n - 1) / 2)
end


# gyration radius of a sphere
function radiusSphere(length)
    r = (3*length/(4*pi))^(1/3)
    return sqrt(3/5)*r
end


# density of particles of species s inside a bounding box
# bbox[1] = [xl, xh], bbox[2] = [yl, yh], bbox[3] = [zl, zh]
function boxDensity(s::AbstractSpecies, bbox::Vector{Vector{Int32}})
    count = 0
    pos = zeros(Int32, 3)
    xl, xh = bbox[1]
    yl, yh = bbox[2]
    zl, zh = bbox[3]
    vol = abs(xh - xl + 1) * abs(yh - yl + 1) * abs(zh - zl + 1)

    for x in range(xl, xh)
        for y in range(yl, yh)
            for z in range(zl, zh)
                boxmod!([x,y,z], pos, Int32(1), s.boxlength)
                count += (s.latS[pos[1], pos[2], pos[3]] > 0)
            end
        end
    end

    return count / vol
end

# axial density along a parallelepiped with a given thickness
function axialDensity(s::AbstractSpecies, axis::String, 
                      pos::Vector{Int32}, thick::Int32)

end

# radial density along a disk with given thickness 
function radialDensity(s::AbstractSpecies, axis::String, 
                      pos::Vector{Int32}, thick::Int32)

end

########### Poly stats ###########
# MARK: -Poly

# extracts the coordinates of a polymer's beads
function extractCoords!(p::Polymer, coords::Matrix)
    key = p.head           # initial bead
    k = 1
    
    for (i, k) in enumerate(collect(keys(p.pos)))
        if key != 0
            coords[i, 1:3] = p.pos[key][1:3] # next position
            coords[i, 4] = p.pos[key][6]  # monomer type
            key = p.pos[key][5]    # next bead
        end
    end

    return coords
end

# find patterns in a polymer based on monomer type (intIdx)
# Returns a Dict mapping each intIdx to a vector of (start_bead, end_bead) tuples
# where start_bead and end_bead are 1-based indices counting from the polymer's head
#
# Example:
#   If polymer has pattern: [1,1,2,2,2,1,3,3,1,1]
#   Returns: Dict(
#       1 => [(1,2), (6,6), (9,10)],
#       2 => [(3,5)],
#       3 => [(7,8)]
#   )
function find_patterns(poly::Polymer)
    patterns = Dict{Int32, Vector{Tuple{Int,Int}}}()
    
    # Edge case: empty polymer
    if poly.length == 0 || poly.head === nothing
        return patterns
    end
    
    # Traverse the polymer linked list
    key = poly.head
    bead_num = 1
    
    current_intIdx = poly.pos[key][6]
    pattern_start = 1
    
    # Initialize dict entry for first monomer type
    if !haskey(patterns, current_intIdx)
        patterns[current_intIdx] = Vector{Tuple{Int,Int}}()
    end
    
    for i in 1:poly.length
        next_key = poly.pos[key][5]
        
        # Check if we're at the end or if the monomer type changes
        if next_key == 0 || i == poly.length
            # Close current pattern
            push!(patterns[current_intIdx], (pattern_start, bead_num))
            break
        else
            next_intIdx = poly.pos[next_key][6]
            
            if next_intIdx != current_intIdx
                # Close current pattern
                push!(patterns[current_intIdx], (pattern_start, bead_num))
                
                # Start new pattern
                current_intIdx = next_intIdx
                pattern_start = bead_num + 1
                
                # Initialize dict entry if needed
                if !haskey(patterns, current_intIdx)
                    patterns[current_intIdx] = Vector{Tuple{Int,Int}}()
                end
            end
        end
        
        key = next_key
        bead_num += 1
    end
    
    return patterns
end


# helper function to calculate pattern center positions
# patterns: vector of tuples [(start1, end1), (start2, end2), ...]
function get_pattern_centers(poly::Polymer, patterns::Vector{Tuple{Int,Int}})
    if length(patterns) != 2
        error("get_pattern_centers requires exactly 2 patterns, got $(length(patterns))")
    end
    
    # check if polymer is long enough for patterns
    max_bead = maximum([p[2] for p in patterns])
    if poly.length < max_bead
        return nothing, nothing, Inf
    end
    
    centers = [zeros(3) for _ in patterns]
    pattern_lengths = [p[2] - p[1] + 1 for p in patterns]
    
    # iterate through polymer linked list
    key = poly.head
    bead_num = 1
    
    for i in 1:poly.length
        pos = poly.pos[key][1:3]
        
        # check which pattern(s) this bead belongs to
        for (pidx, (pstart, pend)) in enumerate(patterns)
            if bead_num >= pstart && bead_num <= pend
                centers[pidx] .+= pos
            end
        end
        
        key = poly.pos[key][5]  # next bead
        bead_num += 1
    end
    
    # normalize by pattern lengths
    for (i, plen) in enumerate(pattern_lengths)
        centers[i] ./= plen
    end
    
    # euclidean distance between the two centers
    dist = sqrt(sum((centers[1] .- centers[2]).^2))
    
    return centers[1], centers[2], dist
end

# ->>> new dispatch
function get_pattern_centers(poly1::Polymer, poly2::Polymer, patterns::Vector{Tuple{Int,Int}})
    if length(patterns) != 2
        error("get_pattern_centers requires exactly 2 patterns, got $(length(patterns))")
    end
    
    centers = [zeros(3), zeros(3)]
    polymers = [poly1, poly2]
    
    # compute center for each polymer's pattern
    for (pidx, (poly, (pstart, pend))) in enumerate(zip(polymers, patterns))
        # check if polymer is long enough for its pattern
        if poly.length < pend
            return nothing, nothing, Inf
        end
        
        pattern_length = pend - pstart + 1
        
        # iterate through polymer linked list
        key = poly.head
        bead_num = 1
        
        for i in 1:poly.length
            pos = poly.pos[key][1:3]
            
            # check if this bead belongs to the pattern
            if bead_num >= pstart && bead_num <= pend
                centers[pidx] .+= pos
            end
            
            key = poly.pos[key][5]  # next bead
            bead_num += 1
        end
        
        # normalize by pattern length
        centers[pidx] ./= pattern_length
    end
    
    # euclidean distance between the two centers
    dist = sqrt(sum((centers[1] .- centers[2]).^2))
    
    return centers[1], centers[2], dist
end

# helper function to check if patterns are in contact
# patterns: vector of tuples [(start1, end1), (start2, end2), ...]
function patterns_in_contact(poly::Polymer, patterns::Vector{Tuple{Int,Int}})
    if length(patterns) != 2
        error("patterns_in_contact requires exactly 2 patterns, got $(length(patterns))")
    end
    
    # check if polymer is long enough
    max_bead = maximum([p[2] for p in patterns])
    if poly.length < max_bead
        return false
    end
    
    # collect positions for each pattern
    pattern_positions = [[] for _ in patterns]
    
    key = poly.head
    bead_num = 1
    
    for i in 1:poly.length
        pos = poly.pos[key][1:3]
        
        for (pidx, (pstart, pend)) in enumerate(patterns)
            if bead_num >= pstart && bead_num <= pend
                push!(pattern_positions[pidx], copy(pos))
            end
        end
        
        key = poly.pos[key][5]
        bead_num += 1
    end
    
    # check all pairs of beads between pattern 1 and pattern 2
    for pos1 in pattern_positions[1]
        for pos2 in pattern_positions[2]
            # check if adjacent (within 1 lattice site in each dimension, periodic BC)
            dx = abs(pos1[1] - pos2[1])
            dy = abs(pos1[2] - pos2[2])
            dz = abs(pos1[3] - pos2[3])
            
            # account for periodic boundaries
            dx = min(dx, poly.boxlength - dx)
            dy = min(dy, poly.boxlength - dy)
            dz = min(dz, poly.boxlength - dz)
            
            if dx <= 1 && dy <= 1 && dz <= 1 && (dx + dy + dz) > 0
                return true
            end
        end
    end
    
    return false
end

########### BulkMonomer stats ###########
# MARK: -Fast/BulkMonomer

# depth-first search: takes a bead index bi and marks all beads
# traversed starting from bi (adds their indices to the marked list)
# WARNING: will break if called on Poly TODO: fix
function dfs(s::AbstractSpecies, bi::Int32)
    marked = Set{Int32}()
    stack_ = OrderedSet{Int32}()

    push!(marked, bi)
    push!(stack_, bi)

    while !isempty(stack_)
        bi = pop!(stack_)
        nnIdxs!(s, view(s.pos[bi], range(1,3)))

        for bj in s.nnIdxs
            if bj != 0 && !(bj in marked)
                push!(stack_, bj)
                push!(marked, bj)
            end
        end
    end

    return marked
end


# connected components
# WARNING: will break if called on Poly TODO: fix
function cc(s::AbstractSpecies)
    ccs = Vector{Set{Int32}}()
    unvisited = Set{Int32}(keys(s.pos))  # all keys

    while !isempty(unvisited)
        bi = pop!(unvisited)
        component = dfs(s, bi)

        for bj in component
            if bi != bj
                pop!(unvisited, bj)
            end
        end

        push!(ccs, component)
    end

    return ccs
end


# returns the exact volume of a droplet. Length of CCs in the empty space
# inside the bounding box, subtracts that from the volume of the box
function dropletVolume(s::AbstractBulkMonomer, droplet::Set{Int32})

    # bounding box around a droplet
    bbox = dropletBbox(s, droplet)

    # fill up the bounding box in the ghost lattice
    fillBbox!(s, droplet, bbox)

    # connected components adjacent to bbox walls around a droplet
    ccs = cc(s.ghost)
    length_ccs = sum(map(x -> onEdge(s.ghost, x, bbox) ? length(x) : 0, ccs))

    # cleans up sites in the ghost lattice
    clearBbox!(s, bbox)

    return vol(bbox) - length_ccs
end


# exact density of each droplet in a list 
function dropletDensities(s::AbstractBulkMonomer, droplets::Vector{Set{Int32}})
    densities = Vector{Float64}(undef, 0)

    for d in droplets
        push!(densities, length(d) / dropletVolume(s, d))
    end

    return densities
end


# (number of particles - total droplet size) / (volume - total droplet volume)
function diluteDensity(s::AbstractBulkMonomer, droplets::Vector{Set{Int32}})
    
    n = s.length
    ndrop = sum(map(x -> length(x), droplets))
    voldrop = 0

    for d in droplets
        voldrop += dropletVolume(s, d)
    end

    return (n - ndrop) / (s.boxlength^3 - voldrop)
end

# number of particles of a species bound to any of the binding sites
function nspeciesBound(sp::AbstractBulkMonomer, 
                       bsites::Vector{T}) where T <: AbstractSpecies
    n = 0
    for site in bsites
        for (_, val) in site.pos
            n += (sp.latS[val[1], val[2], val[3]] > 0)
        end
    end
    return n
end

# computes euclidean distance between two positions accounting for periodic boundaries
function distance_periodic(pos1::Vector, pos2::Vector, boxlength::Int32)
    dx = abs(pos1[1] - pos2[1])
    dy = abs(pos1[2] - pos2[2])
    dz = abs(pos1[3] - pos2[3])
    
    # account for periodic boundaries
    dx = min(dx, boxlength - dx)
    dy = min(dy, boxlength - dy)
    dz = min(dz, boxlength - dz)
    
    return sqrt(dx^2 + dy^2 + dz^2)
end

# accumulates contacts into a matrix by traversing polymer linked list
# contact_threshold: distance cutoff for considering beads in contact (default 2.0)
function compute_contact_matrix!(poly::Polymer, contact_matrix::Matrix{Int}, 
                                contact_threshold::Float64 = 2.0)
    # traverse polymer and compute all pairwise distances
    key_i = poly.head
    bead_i = 1
    
    for i in 1:poly.length
        pos_i = poly.pos[key_i][1:3]
        
        # check contacts with all subsequent beads
        key_j = poly.pos[key_i][5]  # next bead after i
        bead_j = bead_i + 1
        
        for j in bead_j:poly.length
            pos_j = poly.pos[key_j][1:3]
            
            dist = distance_periodic(pos_i, pos_j, Int32(poly.boxlength))
            
            if dist <= contact_threshold
                contact_matrix[bead_i, bead_j] += 1
                contact_matrix[bead_j, bead_i] += 1
            end
            
            key_j = poly.pos[key_j][5]
        end
        
        key_i = poly.pos[key_i][5]
        bead_i += 1
    end

    return contact_matrix
end

# updates density profile by accumulating counts of occupied lattice sites
# as a function of euclidean distance from a given origin (continuous position)
# density[i] accumulates the count of occupied sites at distance floor(i-1) to floor(i) from origin
# origin: continuous 3D position (Float64) from which to measure distances
# WARNING: ignores periodic boundaries, assumes origin stays near lattice center
function updateNormalHistogram!(
        s::AbstractSpecies, 
        hist::Vector{Float64}, 
        origin::Vector{Float64},
        maxDistance::Float64
    )

    nbins = length(hist)
    count = 0
    
    # iterate through all occupied sites in the lattice
    for (_, pos) in s.pos
        # compute minimum image distance from origin to this lattice site
        dx = abs(Float64(pos[1]) - origin[1])
        dy = abs(Float64(pos[2]) - origin[2])
        dz = abs(Float64(pos[3]) - origin[3])
        
        # euclidean distance
        dist = sqrt(dx^2 + dy^2 + dz^2)

        dist > maxDistance && continue
        dist < 0.5 && continue
        
        # bin index
        bin_idx = Int(floor(nbins * dist / maxDistance)) + 1

        # width of spherical bin
        dr = maxDistance / nbins
        
        # accumulate count if within bounds of density vector
        if bin_idx <= nbins
            hist[bin_idx] += 1 / (4 * pi * dist^2 * dr)
            count += 1
        end
    end

    return count
end


########### Rouse trajectory stats ###########
# MARK: -Rouse

# Detects rising contacts (distance crosses below threshold) and returns:
#   times  — wall-clock time of each contact entry
#   invCdf — sorted inter-contact intervals (the empirical inverse CDF;
#            pair with `1 .- (1:length(invCdf))/length(invCdf)` for survival)
function contactTimes(distanceTraj::Vector{Float64}, threshold::Float64, dt::Float64)
    times = Float64[]
    @inbounds for i in 2:length(distanceTraj)
        if distanceTraj[i - 1] > threshold && distanceTraj[i] <= threshold
            push!(times, (i - 1) * dt)
        end
    end
    invCdf = sort(diff(times))
    return times, invCdf
end

# Histogram of a 1D scalar trajectory (e.g. bead-pair distance).
# Returns bin centers and integer counts.
function distanceHistogram(distanceTraj::Vector{Float64}, nbins::Int;
                           range::Union{Nothing, Tuple{Float64,Float64}} = nothing)
    lo, hi = range === nothing ? (minimum(distanceTraj), maximum(distanceTraj)) : range
    if hi <= lo
        hi = lo + max(eps(lo), 1.0)
    end
    # Nudge upper edge outward so that points equal to `hi` land in the last bin
    # (StatsBase's default histogram is left-closed, right-open).
    hi += (hi - lo) * 1e-12 + eps(hi)
    edges = collect(LinRange(lo, hi, nbins + 1))
    h = fit(Histogram, distanceTraj, edges)
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    return centers, collect(h.weights)
end

# Time-averaged MSD of a single bead from a full position trajectory.
# posTraj: (nsteps, N, 3); returns Vector{Float64} of length maxLag.
function beadMSD(posTraj::Array{Float64,3}, beadIdx::Integer;
                 maxLag::Integer = div(size(posTraj, 1), 4))
    nsteps = size(posTraj, 1)
    msd = zeros(Float64, maxLag)
    @inbounds for tau in 1:maxLag
        acc = 0.0
        ncount = nsteps - tau
        for t in 1:ncount
            dx = posTraj[t + tau, beadIdx, 1] - posTraj[t, beadIdx, 1]
            dy = posTraj[t + tau, beadIdx, 2] - posTraj[t, beadIdx, 2]
            dz = posTraj[t + tau, beadIdx, 3] - posTraj[t, beadIdx, 3]
            acc += dx * dx + dy * dy + dz * dz
        end
        msd[tau] = acc / ncount
    end
    return msd
end

# Time-averaged MSD of the scalar distance |r_i - r_j| from a full position trajectory.
function beadPairDistanceMSD(posTraj::Array{Float64,3}, i::Integer, j::Integer;
                             maxLag::Integer = div(size(posTraj, 1), 4))
    nsteps = size(posTraj, 1)
    d = Vector{Float64}(undef, nsteps)
    @inbounds for t in 1:nsteps
        dx = posTraj[t, i, 1] - posTraj[t, j, 1]
        dy = posTraj[t, i, 2] - posTraj[t, j, 2]
        dz = posTraj[t, i, 3] - posTraj[t, j, 3]
        d[t] = sqrt(dx * dx + dy * dy + dz * dz)
    end
    msd = zeros(Float64, maxLag)
    @inbounds for tau in 1:maxLag
        acc = 0.0
        ncount = nsteps - tau
        for t in 1:ncount
            delta = d[t + tau] - d[t]
            acc += delta * delta
        end
        msd[tau] = acc / ncount
    end
    return msd
end
