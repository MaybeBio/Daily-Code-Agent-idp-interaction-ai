export initHam!, initUniform!, initLine!, initBox!, initSphere!, initMovableSet!, 
       initMonomers!, initShape!, initSpecies!, init!, updateMonomers!,
       growPolymer!, updatePattern!, initLLR!, initDecoyWeights!, init2bsWeights!, 
       init1bsWeights!, initDecoyWeightsNew!, initFreeWeights!, initQ!, initOneBurst!, 
       initWeights!, initFull!, initAcc!

############# Species initialization functions #############

######## Generic #######
# MARK: -Generic

# initialize hamiltonian
function initHam!(s::AbstractSpecies)
    ham!(s)
    return s
end

# initial configurations
function initUniform!(pos::OrderedDict{Int32, Vector{Int32}}, boxlength::Int32, length::Int32)
    # sample uniformly randomly
    pos_idxs = sample(0:(boxlength^3-1), length, replace=false)
    
    for i in 1:length
        pos[i] = zeros(Int32, 6)
        mx = mod(pos_idxs[i], boxlength^2)
        
        pos[i][1] = div(pos_idxs[i], boxlength^2) + 1 # x coordinate
        pos[i][2] = div(mx, boxlength) + 1            # y coordinate
        pos[i][3] = mod(mx, boxlength) + 1            # z coordinate
        pos[i][4] = i - 1
        pos[i][5] = i + 1
        # pos[i][6] = intIdx, assigned later by initMonomers!
    end

    return pos
end

function initLine!(pos::OrderedDict{Int32, Vector{Int32}}, 
                   initPos::Union{Nothing, Vector{Int32}}, 
                   boxlength::Int32,
                   length::Int32)

    # line from center of one face to the other along x axis
    for i in 1:length
        pos[i] = zeros(Int32, 6)
        pos[i][1] = initPos === nothing ? i : i + initPos[1] - 1
        pos[i][2] = initPos === nothing ? div(boxlength, 2) : initPos[2]
        pos[i][3] = initPos === nothing ? div(boxlength, 2) : initPos[3]
        pos[i][4] = i - 1
        pos[i][5] = i + 1
        # pos[i][6] = intIdx, assigned later by initMonomers!
    end

    return pos
end

function initBox!(pos::OrderedDict{Int32, Vector{Int32}}, 
                  latticelength::Int32,
                  length::Int32,
                  startPos::Union{Nothing, Vector{Int32}},
                  newIndex::Int32 = Int32(1))

    # filled box shape
    idx = newIndex
    box_side = ceil(Int32, length^(1/3))
    offset = Int32(div(box_side, 2))
    clippedPos = zeros(Int32, 3)
    rawPos = zeros(Int32, 3)
    
    for x in 0:(box_side-1)
        for y in 0:(box_side-1)
            for z in 0:(box_side-1)
                if idx >= length + newIndex
                    break
                end

                pos[idx] = zeros(Int32, 6)
                rawPos = [startPos[1] + Int32(x) - offset,
                          startPos[2] + Int32(y) - offset,
                          startPos[3] + Int32(z) - offset]
                boxmod!(rawPos, clippedPos, Int32(1), latticelength)

                pos[idx][1] = clippedPos[1]
                pos[idx][2] = clippedPos[2]
                pos[idx][3] = clippedPos[3]
                pos[idx][4] = idx - 1
                pos[idx][5] = idx + 1
                # pos[idx][6] = intIdx, assigned later
                idx += 1
            end
            if idx >= length + newIndex
                break
            end
        end
        if idx >= length + newIndex
            break
        end
    end
end

# Convenience method that converts Int64 to Int32
function initBox!(pos::OrderedDict{Int32, Vector{Int32}}, 
                  latticelength::Int32,
                  length::Int32,
                  startPos::Union{Nothing, Vector{Int32}},
                  newIndex::Int64)
    initBox!(pos, latticelength, length, startPos, Int32(newIndex))
end

function initSphere!(pos::OrderedDict{Int32, Vector{Int32}},
                     initPos::Vector{Int32},
                     boxlength::Int32,
                     radius::Int32)
    
    # Validate that sphere doesn't cross periodic boundaries
    center = initPos
    if any(center .- radius .< 1) || any(center .+ radius .> boxlength)
        error("Sphere with center at $center and radius $radius crosses box boundaries (boxlength=$boxlength)")
    end
    
    # Generate sphere by checking distance from center
    idx = Int32(1)
    
    # Iterate over bounding box around sphere
    for x in (center[1] - radius):(center[1] + radius)
        for y in (center[2] - radius):(center[2] + radius)
            for z in (center[3] - radius):(center[3] + radius)
                # Calculate distance from center
                dx = x - center[1]
                dy = y - center[2]
                dz = z - center[3]
                dist = sqrt(Float64(dx*dx + dy*dy + dz*dz))
                
                # Include point if within radius
                if dist <= radius
                    pos[idx] = zeros(Int32, 6)
                    pos[idx][1] = x
                    pos[idx][2] = y
                    pos[idx][3] = z
                    pos[idx][4] = idx - 1
                    pos[idx][5] = idx + 1
                    # pos[idx][6] = intIdx, assigned later
                    idx += 1
                end
            end
        end
    end
    
    # Actual number of beads may differ slightly from theoretical volume
    actual_length = idx - 1
    
    return pos, actual_length
end

# populate the set of movable particles of a species
function initMovableSet!(s::AbstractSpecies)
    for i in collect(keys(s.pos))
        nns = nn(s, s.pos[i][range(1,3)])
        if nns < 6
            push!(s.movableIdxs, i)
        end
    end

    return s
end

######## Species #######
# MARK: -Species

# initializes empty dict of bead locations and other containers
function initSpecies!(s::Species)
    s.name = string("species", s.id)
    s.pos = OrderedDict{Int32, Vector{Int32}}()

    # head and tail beads
    s.head = 1
    s.tempVec = zeros(Int32, 3)
    s.nnIdxs = zeros(Int32, 6)
    s.freeNNs = Set{Int32}()

    return s
end


######## Polymer #######
# MARK: -Polymer

"""
    initMonomers!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)

Initialize monomers for a polymer with bead-dependent interactions.
Creates a dictionary of PolyMonomers keyed by intIdx.

If pos is already initialized, also assigns each bead's monomer type in pos[key][6].
If pos is not yet initialized, pattern assignment is deferred to initSpecies!.

# Arguments
- `p::Polymer`: The polymer to initialize monomers for
- `defaultChemPot::Float64`: Default chemical potential for all monomers (default: 0.0)
- `defaultIntPot::Float64`: Default interaction potential for all monomers (default: 0.0)

# Pattern options (set via p.pattern field):
- `:uniform`: All beads belong to monomer type 1 (default)
- `:alternating`: Alternates between types (1, 2, 1, 2, ...)
- `:block`: Divides polymer into n equal blocks
- `:random`: Random assignment
- `:custom`: Use p.customPattern vector to assign beads

# Notes
- pos[key][6] stores the intIdx directly for O(1) monomer lookup
- p.monomers is a Dict{Int32, PolyMonomer} keyed by intIdx
- All monomers share a synchronized intMat
"""
function initMonomers!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)
    
    # Create monomers dictionary keyed by intIdx
    p.monomers = Dict{Int32, PolyMonomer}()
    
    # Shared structures for all monomers of this polymer
    sharedIntSpecies = nothing
    sharedIntMat = nothing
    
    for i in 1:p.nMonomerTypes
        m = PolyMonomer(
            name = string(p.name, "_monomer", i),
            poly = p,
            chemPot = defaultChemPot,
            intPot = defaultIntPot
        )
        
        # Share intSpecies/intMat from first monomer so subsequent get unique intIdx
        if sharedIntSpecies !== nothing
            m.intSpecies = sharedIntSpecies
            m.intMat = sharedIntMat
        end
        
        _ensureIntStructures!(m)
        p.monomers[m.intIdx] = m
        
        # Capture/update shared structures
        sharedIntSpecies = m.intSpecies
        sharedIntMat = m.intMat
    end
    
    # Ensure all monomers share the final expanded matrix
    for (_, m) in p.monomers
        m.intMat = sharedIntMat
        m.intSpecies = sharedIntSpecies
    end
    
    # If pos exists, assign pattern now; otherwise defer to initSpecies!
    if p.pos !== nothing
        _assignPattern!(p)
    end
    
    return p
end

"""
    updateMonomers!(p::Polymer, baseIntIdx::Int32)

Re-key the polymer's monomers dictionary after intIdx values have been finalized
by intPot! calls. This ensures poly.monomers keys match the global intIdx values
of each monomer.

# Arguments
- `p::Polymer`: The polymer whose monomers dictionary needs updating
- `baseIntIdx::Int32`: The intIdx of the base monomer type for this polymer

This function should be called after all intPot! calls and before init!.

# Example
    ```julia
    poly1 = Polymer(length=32, nMonomerTypes=2, pattern=:custom)
    initMonomers!(poly1, defaultChemPot=1.0)

    base1 = poly1.monomers[1]
    pattern1 = poly1.monomers[2]

    # Set interaction potentials (changes intIdx values)
    intPot!(base1, 0.5)
    intPot!(pattern1, -2.0)

    # Re-key monomers dictionary with new intIdx values
    updateMonomers!(poly1, base1.intIdx)
```
"""
function updateMonomers!(p::Polymer, baseIntIdx::Int32)
    # Create new dictionary with current intIdx values as keys
    new_monomers = Dict{Int32, PolyMonomer}()
    
    for (old_key, monomer) in p.monomers
        new_monomers[monomer.intIdx] = monomer
    end
    
    # Replace old dictionary
    p.monomers = new_monomers
    
    # Update baseIntIdx to the specified value
    p.baseIntIdx = baseIntIdx

    if p.pos !== nothing
        _assignPattern!(p)
    end
    
    return p
end

# Internal helper called by init! to ensure the polymer's monomers dictionary
# keys match the current intIdx values of each monomer. Handles cases where
# intPot! changed intIdx after initMonomers! without a manual updateMonomers! call.
# No-op if keys already match.
function _syncMonomerKeys!(p::Polymer)
    # Check if any key doesn't match its monomer's intIdx
    needs_sync = false
    for (key, monomer) in p.monomers
        if key != monomer.intIdx
            needs_sync = true
            break
        end
    end
    
    if !needs_sync
        return p
    end
    
    # Re-key dictionary with current intIdx values
    new_monomers = Dict{Int32, PolyMonomer}()
    for (_, monomer) in p.monomers
        new_monomers[monomer.intIdx] = monomer
    end
    p.monomers = new_monomers
    
    # Find baseIntIdx as the smallest intIdx
    intIdxList = sort(collect(keys(p.monomers)))
    p.baseIntIdx = intIdxList[1]
    
    # Re-assign pattern to pos[key][6] with correct intIdx values
    if p.pos !== nothing
        _assignPattern!(p)
    end
    
    return p
end

# Validate that each non-base pattern type (2+) appears in exactly one contiguous block
function _validateContiguousBlocks(patternVec::Vector{Int32}, baseIntIdx::Int32)
    patternTypes = []

    for k in patternVec
        if k != baseIntIdx && !(k in patternTypes)
            push!(patternTypes, k)
        end
    end

    for t in patternTypes
        indices = findall(x -> x == t, patternVec)
        if !isempty(indices) && indices != collect(indices[1]:indices[end])
            error("Pattern type $t must appear in exactly one contiguous block")
        end
    end
end

# Initialize patternTails and baseIntIdx
function _initPatternTails!(p::Polymer)
    # For uniform pattern, no tails to track
    if p.nMonomerTypes == 1 || p.pattern == :uniform
        p.patternTails = nothing
        return
    end
    
    p.patternTails = Dict{Int32, Int32}()
    
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
end

# Assign monomer types to pos[key][6] based on pattern
function _assignPattern!(p::Polymer)
    if p.baseIntIdx === nothing
        intIdxList = sort(collect(keys(p.monomers)))
        p.baseIntIdx = intIdxList[1]  # base type = first intIdx (smallest)
    end

    # Build pattern vector (1-based type indices)
    patternVec = _buildPatternVector(p)
    
    # Validate contiguous blocks for non-base types
    _validateContiguousBlocks(patternVec, p.baseIntIdx)
    
    # Traverse chain and assign intIdx to pos[key][6]
    key = p.head
    chainPos = 1
    while key != 0 && chainPos <= p.length
        typeIdx = patternVec[chainPos]
        p.pos[key][6] = typeIdx
        key = p.pos[key][5]
        chainPos += 1
    end
    
    # Initialize pattern boundary tracking
    _initPatternTails!(p)
end

# Build pattern vector using intIdx values
function _buildPatternVector(p::Polymer)::Vector{Int32}
    patternVec = zeros(Int32, p.length)
    
    # Get sorted list of intIdx values for mapping
    intIdxList = sort(collect(keys(p.monomers)))
    
    if p.pattern == :uniform
        fill!(patternVec, p.baseIntIdx)
        
    elseif p.pattern == :alternating
        for i in 1:p.length
            typeIndex = mod1(i, p.nMonomerTypes)
            patternVec[i] = intIdxList[typeIndex]
        end
        
    elseif p.pattern == :block
        blockSize = div(p.length, p.nMonomerTypes)
        for i in 1:p.length
            typeIndex = min(div(i - 1, blockSize) + 1, p.nMonomerTypes)
            patternVec[i] = intIdxList[typeIndex]
        end
        
    elseif p.pattern == :random
        for i in 1:p.length
            typeIndex = rand(1:p.nMonomerTypes)
            patternVec[i] = intIdxList[typeIndex]
        end
        
    elseif p.pattern == :custom
        @assert p.customPattern !== nothing "pattern=:custom requires customPattern"
        @assert length(p.customPattern) >= p.length "customPattern too short"
        for i in 1:p.length
            @assert p.customPattern[i] in keys(p.monomers) "customPattern[$i] out of range"
            patternVec[i] = p.customPattern[i]
        end
        
    else
        error("Unknown pattern: $(p.pattern)")
    end
    
    return patternVec
end

# initializes species configuration from a list of presets
function initShape!(s::Polymer)
    if s.initShape == "line" || s.initShape === nothing
        initLine!(s.pos, s.initPos, s.boxlength, s.length)
    else
        error("Invalid string for s.initShape")
    end
    
    # head and tail beads
    s.pos[1][4] = 0
    s.pos[s.length][5] = 0
    
    return s
end

# initializes bead locations of a long Polymer in a doubly-linked list
function initSpecies!(p::Polymer)
    p.name = string("poly", p.id)
    p.pos = OrderedDict{Int32, Vector{Int32}}()
    
    initShape!(p)

    # head and tail beads
    p.head = 1

    # assign monomer pattern to pos[key][6] if monomers exist
    if p.monomers !== nothing
        _assignPattern!(p)
    end

    # key vector
    p.keys_vec = collect(keys(p.pos))

    p.tempVec = zeros(Int32, 3)
    p.bondVec = zeros(Int32, 3)
    p.moveVec = zeros(Int32, 3)
    p.newPosi = zeros(Int32, 3)
    p.newPosj = zeros(Int32, 3)
    p.moveTypes = zeros(Int32, 4)
    p.moveTypeProbs = [0,PKINK,0]
    p.nTypes = zeros(Int32, 3)
    p.hastingsFactor = 1.0
    
    # Pre-allocate bond coordinate vectors (avoids view allocations in identifyBond!)
    p.bi = zeros(Int32, 3)
    p.bj = zeros(Int32, 3)
    p.bprev = zeros(Int32, 3)
    p.bnext = zeros(Int32, 3)
    
    # Initialize deltaHP for partner energy accumulation
    p.deltaHP = 0.0

    return p
end

######## BulkMonomer #######
# MARK: -BulkMonomer

# initializes species configuration from a list of presets
function initShape!(s::AbstractBulkMonomer)
    if s.initShape == "uniform" || s.initShape === nothing
        initUniform!(s.pos, s.boxlength, s.length)

    elseif s.initShape == "line"    
        initLine!(s.pos, s.initPos, s.boxlength, s.length)

    elseif s.initShape == "sphere"
        if s.initPos === nothing
            error("initShape='sphere' requires initPos to specify center")
        end
        if s.initParams === nothing
            error("initShape='sphere' requires initParams to specify radius")
        end
        
        _, actual_length = initSphere!(s.pos, s.initPos, s.boxlength, s.initParams)
        s.length = actual_length  # Update length to actual sphere size

    elseif s.initShape == "box"
        beadsPerBox = s.length
        totalBeads = beadsPerBox * s.initParams
        
        for i in 1:s.initParams
            initPos = zeros(Int32, 3)
            offset = div(s.boxlength, 2 * s.initParams)
            
            if !(s.initPos === nothing)
                initPos[1] = s.initPos[1] + offset * (2 * i - s.initParams - 1)
                initPos[2] = s.initPos[2]
                initPos[3] = s.initPos[3]
            else
                initPos[1] = div(s.boxlength, 2) + offset * (2 * i - s.initParams - 1)
                initPos[2] = div(s.boxlength, 2)
                initPos[3] = div(s.boxlength, 2)
            end

            startIdx = beadsPerBox * (i - 1) + 1
            initBox!(s.pos, s.boxlength, beadsPerBox, initPos, startIdx)
        end
        
        s.length = totalBeads   # total length = length of each box * nboxes
    else
        error("Invalid string for s.initShape")
    end
    
    # head and tail beads
    s.pos[1][4] = 0
    s.pos[s.length][5] = 0
    
    return s
end

# initializes beads positioned uniformly randomly (dilute and well-mixed)
function initSpecies!(b::BulkMonomer)

    b.name = string("BulkMonomer", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()
    
    initShape!(b)

    _ensureIntStructures!(b)
    
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0
    b.nnsOld = 0
    b.nnsNew = 0
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)
    b.nnIdxs = zeros(Int32, 6)
    b.deltaH = 0.0
    b.deltaHP = 0.0

    if b.chemPot === nothing b.chemPot = 0.0 end
    if b.intPot === nothing b.intPot = 0.0 end

    return b
end


######## FastBulkMonomer #######
# MARK: -FastBulkMonomer

# initializes beads positioned uniformly randomly (dilute and well-mixed)
function initSpecies!(b::FastBulkMonomer)

    b.name = string("BulkMonomer", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()

    initShape!(b)

    _ensureIntStructures!(b)

    # build set of movable particles
    b.movableIdxs = Set{Int32}()
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0

    b.nnsOld = 0
    b.nnsNew = 0
    b.bi = zeros(Int32, 3)
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)
    b.moveVec = zeros(Int32, 2)
    b.nnIdxs = zeros(Int32, 6)
    b.deltaH = 0.0
    b.deltaHP = 0.0

    if b.chemPot === nothing b.chemPot = 0.0 end
    if b.intPot === nothing b.intPot = 0.0 end

    return b
end


############# System initialization functions #############
# MARK: -System

# initializes lattice sites of each species from their bead locations
function init!(l::Lattice, slist::Vector{<:AbstractSpecies})

    # Separate species that need lattice layers from monomers that delegate
    layeredSpecies = AbstractSpecies[]
    monomers = AbstractPolyMonomer[]
    polyParents = Set{AbstractPolymer}()
    
    for s in slist
        if s isa AbstractPolyMonomer
            push!(monomers, s)
            push!(polyParents, s.poly)
        else
            push!(layeredSpecies, s)
        end
    end
    
    # Combine regular species with parent polymers for lattice allocation
    allLayeredSpecies = vcat(layeredSpecies, collect(polyParents))
    
    @assert length(allLayeredSpecies) <= 26 "Too many species"

    l.nspecies = length(allLayeredSpecies)
    l.lat = zeros((l.nspecies + 1, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int32, Int32}()

    # initialize ghost species
    gs = Species(;)
    gs.latS = view(l.lat, 1, :, :, :)
    gs.boxlength = l.boxlength
    initSpecies!(gs)

    # initialize each layered species (regular species + parent polymers)
    for (i, s) in enumerate(allLayeredSpecies)
        
        # give species s a view on its lattice sites, and on the ghost species
        s.latS = view(l.lat, i + 1, :, :, :)
        s.ghost = gs
        s.boxlength = l.boxlength
        
        # Initialize species if not already done (parent polymers might be pre-initialized)
        if s.pos === nothing
            initSpecies!(s)
        end

        # Auto-sync monomer dict keys to current intIdx values (handles
        # cases where intPot! changed intIdx after initMonomers!)
        if typeof(s) <: AbstractPolymer && s.monomers !== nothing
            _syncMonomerKeys!(s)
        end

        # give each species a unique letter suffix
        s.name = string(s.name, "_", ALPHABET[i])
        l.speciesDir[s.id] = i

        # Populate lattice
        if typeof(s) <: AbstractBulkMonomer
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = j
            end
        elseif typeof(s) <: AbstractPolymer && s.monomers !== nothing
            # For polymers with monomers, pos[key][6] already has intIdx
            key = s.head
            while key != 0
                pos = s.pos[key]
                s.latS[pos[1], pos[2], pos[3]] = pos[6]  # intIdx stored in pos[6]
                key = pos[5]
            end
        else
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
            end
        end

        if typeof(s) <: FastBulkMonomer
            initMovableSet!(s)
        end
    end

    # Initialize Hamiltonians for all layered species
    for s in allLayeredSpecies
        if hasproperty(s, :ham)
            initHam!(s)
        end
    end

    return l
end

############# Polymer growth and pattern update functions #############
# MARK: -Growth and Pattern Update

"""
    growPolymer!(p::Polymer; maxAttempts::Int=1000000)

Grow a polymer to its maximum length using Monte Carlo moves.
Temporarily sets chemical potential to favor growth, then restores original value.

# Arguments
- `p::Polymer`: The polymer to grow
- `maxAttempts::Int`: Maximum number of move attempts (default: 1,000,000)

# Returns
- `p::Polymer`: The grown polymer

# Throws
- `ErrorException`: If polymer cannot reach maxlength within maxAttempts

# Example
```julia
poly = Polymer(length=10, maxlength=50, nMonomerTypes=2, pattern=:uniform)
initMonomers!(poly)
# ... set interaction potentials ...
init!(l, [poly])

# Grow polymer to maxlength
growPolymer!(poly)
@assert poly.length == poly.maxlength
```
"""
function growPolymer!(p::Polymer; maxAttempts::Int=1000000, minAttempts::Int=1000)
    if p.length >= p.maxlength
        return p  # Already at max length
    end

    @assert minAttempts < maxAttempts "Careful, minAttempts >= maxAttempts"
    
    # Save original chemical potentials
    originalChemPots = Dict{Int32, Float64}()
    for (intIdx, monomer) in p.monomers
        originalChemPots[intIdx] = monomer.chemPot
        # Set high chemical potential to strongly favor growth
        monomer.chemPot = 5.0
    end
    
    attempts = 0
    while p.length < p.maxlength && attempts < maxAttempts || attempts < minAttempts
        move!(p)
        attempts += 1
    end
    
    # Restore original chemical potentials
    for (intIdx, monomer) in p.monomers
        monomer.chemPot = originalChemPots[intIdx]
    end
    
    if p.length < p.maxlength
        error("Failed to grow polymer to maxlength=$(p.maxlength) after $maxAttempts attempts. Current length=$(p.length)")
    end
    
    return p
end

"""
    updatePattern!(l::Lattice, p::Polymer, newPattern::Union{Nothing, Vector{Int32}}=nothing)

Update the pattern of a polymer and synchronize pos and lattice structures.
Must be called after the polymer has reached its target length.

# Arguments
- `l::Lattice`: The lattice containing the polymer
- `p::Polymer`: The polymer whose pattern to update
- `newPattern::Union{Nothing, Vector{Int32}}`: New pattern vector with intIdx values.
  If `nothing`, uses `p.customPattern` (default: nothing)

# Validations
- Verifies polymer is at maxlength
- Verifies pattern length matches polymer length
- Verifies all intIdx values in pattern exist in p.monomers

# Updates
- `p.pos[key][6]` for each bead (sets intIdx)
- `l.lat[x,y,z]` at each bead position
- `p.patternTails` (recalculated)
- `p.customPattern` (if newPattern provided)

# Example
```julia
# After growing polymer to maxlength
growPolymer!(poly)

# Define pattern with intIdx values
base = poly.monomers[1]
pattern = poly.monomers[2]
newPattern = vcat(fill(base.intIdx, 10), fill(pattern.intIdx, 10), fill(base.intIdx, 10))

# Update pattern
updatePattern!(l, poly, newPattern)
```
"""
function updatePattern!(l::Lattice, p::Polymer, newPattern::Union{Nothing, Vector{Int32}}=nothing)
    # Validation 1: Polymer must be at maxlength
    if p.length != p.maxlength
        error("Polymer must be at maxlength before updating pattern. Current length=$(p.length), maxlength=$(p.maxlength)")
    end

    # Determine which pattern to use
    patternToUse = newPattern === nothing ? p.customPattern : newPattern
    
    if patternToUse === nothing
        error("No pattern provided and p.customPattern is nothing")
    end
    
    if length(patternToUse) != p.maxlength
        error("New pattern length must equal maxlength. Current length=$(length(newPattern)), maxlength=$(p.maxlength)")
    end
    
    # Validation 3: All intIdx values must exist in p.monomers
    validIntIdxs = Set(keys(p.monomers))
    for intIdx in patternToUse
        if !(intIdx in validIntIdxs)
            error("Pattern contains intIdx=$intIdx which does not exist in p.monomers (valid: $(sort(collect(validIntIdxs))))")
        end
    end
    
    # Update p.customPattern if new pattern provided
    if newPattern !== nothing
        p.customPattern = copy(newPattern)
        p.pattern = :custom  # Mark as custom pattern
    end
    
    # Update pos[key][6] and lattice for each bead
    key = p.head
    position = 1
    while key != 0
        newIntIdx = patternToUse[position]
        oldIntIdx = p.pos[key][6]
        
        # Update pos
        p.pos[key][6] = newIntIdx
        
        # Update main lattice layer
        x, y, z = p.pos[key][1], p.pos[key][2], p.pos[key][3]
        l.lat[1, x, y, z] = newIntIdx
        
        # Update polymer's latS view
        p.latS[x, y, z] = newIntIdx
        
        key = p.pos[key][5]
        position += 1
    end
    
    # Recalculate patternTails
    _initPatternTails!(p)
    
    return p
end


########## Receptor model initialization ###########
## ------
function initLLR!(m::LLROptModel)

    init_Qs!(m)

    # indicators
    m.off_ind = m.on_ind .⊻ true
    m._idist0_on = _idist0(m.idist, m.on_ind)
    m._idist0_off = _idist0(m.idist, m.off_ind)
end

function initDecoyWeights!(d::OneBurstDecoyModel)
    # transition rates
    # kboff, konb, kuon, kbon, koffu, kuoff
    kboff, konb, kuon, kbon, koffu, kuoff = 10 .^ d.logrates

    # edges = [(1, 2),    (2, 1),        (1, 6),    (6, 1),        (4, 3), (3, 4), 
    #          (4, 5),    (5, 4),        (6, 5),    (4, 1),        (2, 3), (1, 4)]
    d.weights = [kboff,     kuoff,        d.cw*kboff, d.gamma*kuoff, kbon, kuon,
                 d.cw*kbon, d.gamma*kuon, konb,       koffu,         konb]

    return d
end

function init2bsWeights!(d::OneBurstDecoyModel)
    # transition rates
    kboff,   kuoff,  kbon,  kuon,  koffu, kon1b, kon2b = 10 .^ d.logrates

    # edges = [(5, 1),    (1, 5),    (2, 6),    (6, 2),    (1, 2),    (1, 3), 
    #          (2, 4),    (3, 1),    (3, 4),    (4, 2),    (6, 5),    (5, 6)]
    d.weights = [2*kboff,  kuoff,     kuon,     2*kbon,    kon1b,      kboff,
                 kbon,     2*kuoff,   kon2b,    2*kuon,    koffu]

    return d
end

function init1bsWeights!(d::OneBurstDecoyModel)
    # transition rates
    kboff, kuoff, konb, koffb, kuon, kbon, koffu = 10 .^ d.logrates

    d.weights = [kboff, kuoff, konb, koffb, kuon, kbon, koffu, 1e-5]

    return d
end

function initDecoyWeightsNew!(d::OneBurstDecoyModel)
    # transition rates

    kboff, konb, kuon, kbon, koffu, kuoff = 10 .^ d.logrates

    # edges = [(1, 2), (2, 1), (1, 6), (6, 1), (7, 3), (3, 7),
    #          (4, 5), (5, 4), (6, 5), (7, 1), (2, 3), (4, 1), (7, 5), (4, 3)]
    d.weights = [kboff, kuoff, d.cw*kboff, d.gamma*kuoff, kbon, kuon, 
                 d.cw*kbon, d.gamma*kuon, konb, koffu, konb, koffu, d.cw*kbon, kbon]

    return d
end

function initFreeWeights!(d::OneBurstDecoyModel)
    # transition rates
    d.weights = 10 .^ d.logrates

    return d
end

function initQ!(d::OneBurstDecoyModel)
    @assert length(d.weights) == length(d.edges)

    # building weighted adjacency matrix (rate matrix Q)
    weighted_edges = [(u,v,0.) for (u,v) in d.edges]

    for (keyb, val) in enumerate(d.weights)
        (u,v,r) = weighted_edges[keyb]
        weighted_edges[keyb] = (u,v,val)
    end

    d.Q = build_mat(weighted_edges, d.input_edges, d.cr)

    return d
end

function initOneBurst!(d::OneBurstDecoyModel)

    initQ!(d)

    # ON and OFF partitions built from indicator arrays ind_set and ind_comp
    d.ind_comp = ones(Int, length(d.ind_set)) .- d.ind_set
    d._set = findall(isone, d.ind_set)
    d._comp = findall(isone, d.ind_comp)

    # d.Q_h_comp (fpt)
    d.Q_h_comp = d.Q .* 1.0
    for ri in d._comp
        d.Q_h_comp[ri, :] .= 0.0
    end
    UDU_comp = eigen(d.Q_h_comp)
    d.W_comp = UDU_comp.vectors
    d.W_comp_inv = inv(d.W_comp)
    d.D_comp = Diagonal(UDU_comp.values)
    d.aux_comp = zeros(size(d.D_comp)...)
    d.P_comp = zeros(size(d.Q)...)
    # exp_tD_comp = zeros(size(Q)...)

    # Q_h_set (fpt)
    d.Q_h_set = d.Q .* 1.0
    for ri in d._set
        d.Q_h_set[ri, :] .= 0.0
    end
    UDU_set = eigen(d.Q_h_set)
    d.W_set = UDU_set.vectors
    d.W_set_inv = inv(d.W_set)
    d.D_set = Diagonal(UDU_set.values)
    d.aux_set = zeros(size(d.D_comp)...)
    d.P_set = zeros(size(d.Q)...)
    d.D_aux = zeros(size(d.D_comp)...)
    # exp_tD_set = zeros(size(Q)...)

    d.F_set = zeros(size(d.Q)...)
    d.F_comp = zeros(size(d.Q)...)
    
    d.eigvals_set = UDU_set.values
    d.eigvals_comp = UDU_comp.values
    d.maxeigval_set = - maximum(x for x in d.eigvals_set if x < 0)
    d.maxeigval_comp = - maximum(x for x in d.eigvals_comp if x < 0)

    return d
end

function initWeights!(d::OneBurstDecoyModel)
   
    if d.circuitId == "1bs_decoy" || d.circuitId == "1bs_decoy_no_konu"
        initDecoyWeights!(d)

    elseif d.circuitId == "1bs_decoy_no_konu_new"
        initDecoyWeightsNew!(d)

    else
        init1bsWeights!(d)
    end

    return d
end

function initFull!(d::OneBurstDecoyModel)

    initWeights!(d)
    initOneBurst!(d)

    return d
end

function initAcc!(d::OneBurstDecoyModel)

    initWeights!(d)
    initQ!(d)

    return d
end