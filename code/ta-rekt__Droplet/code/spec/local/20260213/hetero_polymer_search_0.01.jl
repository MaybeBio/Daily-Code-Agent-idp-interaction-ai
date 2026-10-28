### physical parameters
# MARK: -Phys. params

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 42,
    nspecies = 2  # number of layered species
)

# heterogeneous polymer with patterned regions
poly1 = Polymer(;
    length = 32,
    maxlength = 500,
    initShape = "line",
    nMonomerTypes = 4,  # base, patternA1, pattern A2
    initPos = [5, 14, 14],
    pattern = :uniform
)

# heterogeneous polymer with patterned regions
poly2 = Polymer(;
    length = 32,
    maxlength = 500,
    initShape = "line",
    nMonomerTypes = 4,  # base, patternB1, patternB2
    initPos = [5, 28, 28],
    pattern = :uniform
)

initMonomers!(poly1, defaultChemPot = 1.0)
initMonomers!(poly2, defaultChemPot = 1.0)

# bulk species
bulk = FastBulkMonomer(;
    length = 1,
    initShape = "uniform"
)

# list containing all species
slist = [poly1, poly2, bulk]


### simulation parameters
# MARK: -Sim. params

blurb = """
Understanding the effects of prewetting on the way in which two strongly
interacting regions of a polymer come into contact.\n

Looking at inter-pattern contacts of two short polymers under tension
as a function of separation, in the presence of wetting
"""

init_from_checkpoint = false

t_max = 480  # runtime in minutes (initial test)
step = 10000
t_burn_in = 0.8  # 12 seconds

# vals = range(0.005, stop=0.03, length=32)  # parallelized parameters (concentration)
vals = [0.01]
alphas = zeros(32)
nrepeats = 32

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [0.02]    # concentration
betas = [0.66]    # temperature
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    poly1, poly2, bulk = slist

    base1 = poly1.monomers[1]
    mainPattern1 = poly1.monomers[2]
    headPattern1 = poly1.monomers[3]
    tailPattern1 = poly1.monomers[4]

    base2 = poly2.monomers[1]
    mainPattern2 = poly2.monomers[2]
    headPattern2 = poly2.monomers[3]
    tailPattern2 = poly2.monomers[4]
    
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base1, 0.5)              # base self-interaction
    intPot!(base1, mainPattern1, 0.5)    # base-mainPattern1
    intPot!(base1, headPattern1, 0.5)    # base-mainPattern1
    intPot!(base1, tailPattern1, 0.5)    # base-mainPattern1
    intPot!(base1, bulk, 0.0)        # base-bulk (weak)
    
    # mainPattern1 (type 2)
    intPot!(mainPattern1, 0.5)           # mainPattern1 self-interaction
    intPot!(headPattern1, 0.5)           # mainPattern1 self-interaction
    intPot!(tailPattern1, 0.5)           # mainPattern1 self-interaction
    intPot!(mainPattern1, bulk, -0.5)    # mainPattern1-bulk (strong)
    intPot!(headPattern1, bulk, 0.0)    # mainPattern1-bulk (strong)
    intPot!(tailPattern1, bulk, 0.0)    # mainPattern1-bulk (strong)

    # base (type 2)
    intPot!(base2, 0.5)              # base self-interaction
    intPot!(base2, mainPattern2, 0.5)    # base-mainPattern1
    intPot!(base2, headPattern2, 0.5)    # base-mainPattern1
    intPot!(base2, tailPattern2, 0.5)    # base-mainPattern1
    intPot!(base2, bulk, 0.0)        # base-bulk (weak)
    
    # mainPattern1 (type 2)
    intPot!(mainPattern2, 0.5)           # mainPattern1 self-interaction
    intPot!(headPattern2, 0.5)           # mainPattern1 self-interaction
    intPot!(tailPattern2, 0.5)           # mainPattern1 self-interaction
    intPot!(mainPattern2, bulk, -0.5)    # mainPattern1-bulk (strong)
    intPot!(headPattern2, bulk, 0.0)    # mainPattern1-bulk (strong)
    intPot!(tailPattern2, bulk, 0.0)    # mainPattern1-bulk (strong)

    intPot!(mainPattern1, mainPattern2, -5.0)    # inter-pattern attraction
    
    # bulk (type 4)
    intPot!(bulk, -1.0 / beta)              # bulk self-interaction

    bulk.length = Int32(round(l.boxlength^3 * alpha))

    updateMonomers!(poly1, base1.intIdx)
    updateMonomers!(poly2, base2.intIdx)

    # initialize lattice and species
    if init_from_checkpoint
        checkpoint_at_param = string(checkpoint_path, "_keya", keya)
        restore!(checkpoint_at_param, l, slist)
    else
        init!(l, slist)
    end

    # Initialize deltaHP for both polymers (needed for multi-polymer interactions)
    poly1.deltaHP = 0.0
    poly2.deltaHP = 0.0

    # moving species
    movingsps = [(poly1, bulk.intIdx + 1),
                 (poly2, bulk.intIdx + 2), 
                 (bulk, bulk.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1
    minsize = 20

    radii1 = []
    radii2 = []
    lengths1 = []
    lengths2 = []
    hams_poly1 = []
    hams_poly2 = []
    hams_bulk = []
    pattern_dists = []
    max_droplet_size = []
    accs = []
    contact_probs = []
    pattern1_position = []
    pattern2_position = []
    
    # initial values
    push!(radii1, radius(poly1))
    push!(radii2, radius(poly2))
    push!(lengths1, poly1.length)
    push!(lengths2, poly2.length)
    push!(pattern_dists, 42)
    push!(hams_poly1, poly1.ham)
    push!(hams_poly2, poly2.ham)
    push!(hams_bulk, bulk.ham)
    push!(accs, 0)
    push!(contact_probs, 0.0)
    
    # build polymer patterns
    patternVec1 = zeros(Int32, poly1.maxlength)
    patternVec2 = zeros(Int32, poly2.maxlength)
    p1_start, p1_end = 226, 275
    p2_start, p2_end = 226, 275

    push!(pattern1_position, p1_start)
    push!(pattern2_position, p2_start)

    patternVec1[1:p1_start-1] .= headPattern1.intIdx         # base type (40 beads)
    patternVec1[p1_start:p1_end] .= mainPattern1.intIdx      # pattern 1 (8 beads)
    patternVec1[p1_end+1:poly1.maxlength-3] .= tailPattern1.intIdx
    patternVec1[poly1.maxlength-2:poly1.maxlength] .= base1.intIdx

    patternVec2[1:p2_start-1] .= headPattern2.intIdx         # base type (40 beads)
    patternVec2[p2_start:p2_end] .= mainPattern2.intIdx      # pattern 1 (8 beads)
    patternVec2[p2_end+1:poly2.maxlength-3] .= tailPattern2.intIdx
    patternVec2[poly2.maxlength-2:poly2.maxlength] .= base2.intIdx

    growPolymer!(poly1)
    growPolymer!(poly2)

    updatePattern!(l, poly1, patternVec1)
    updatePattern!(l, poly2, patternVec2)

    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1

                # record data
                if mod(nmoves, step) == 0
                    ham!(s)

                    ccs = cc(bulk)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    droplet_sizes = map(x -> length(x), droplets)

                    # compute polymer hamiltonian (stored in poly, not monomers)
                    if i == bulk.intIdx + 1
                        push!(hams_poly1, poly1.ham)

                    elseif i == bulk.intIdx + 2
                        push!(hams_poly2, poly2.ham)

                    elseif s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    pattern_poly1 = find_patterns(poly1)
                    p1_start, p1_end = pattern_poly1[mainPattern1.intIdx][1]

                    pattern_poly2 = find_patterns(poly2)
                    p2_start, p2_end = pattern_poly2[mainPattern2.intIdx][1]

                    push!(pattern1_position, p1_start)
                    push!(pattern2_position, p2_start)

                    # compute contact probability between mainPattern1 and bulk
                    pattern_length = p1_end - p1_start + 1
                    contacts = 0
                    
                    # get tail key for mainPattern1
                    if haskey(poly1.patternTails, mainPattern1.intIdx)
                        tail_key = poly1.patternTails[mainPattern1.intIdx]
                        
                        # walk backwards from tail to count contacts
                        current_key = tail_key
                        while current_key != 0
                            # check if this bead is mainPattern1 type
                            if poly1.pos[current_key][6] == mainPattern1.intIdx
                                # count nearest neighbors with bulk
                                contacts += nnp(mainPattern1, bulk, poly1.pos[current_key][1:3]) / 4
                                # move to previous bead
                                current_key = poly1.pos[current_key][4]
                            else
                                break  # reached end of pattern region
                            end
                        end
                    end
                    
                    push!(contact_probs, contacts / pattern_length)

                    _, _, dist = get_pattern_centers(poly1, poly2, 
                                    [(p1_start, p1_end),(p2_start, p2_end)])
                    push!(pattern_dists, dist)

                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii1, radius(poly1))
                    push!(lengths1, poly1.length)
                    push!(radii2, radius(poly2))
                    push!(lengths2, poly2.length)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))

                    frame += 1
                end
            end
        end
    end

    # density_lat = density_lat ./ total

    # save results
    res["blurb"] = blurb
    res["species_name"] = "poly"
    res["nframes"] = length(radii1)
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["interaction_matrix"] = bulk.intMat
    res["boxlength"] = l.boxlength
    res["length_poly1"] = lengths1
    res["length_poly2"] = lengths2
    res["length_bulk"] = bulk.length
    res["radius_of_gyration1"] = radii1
    res["radius_of_gyration2"] = radii2
    res["pattern_distances"] = pattern_dists
    res["mainPattern1_position"] = pattern1_position
    res["mainPattern2_position"] = pattern2_position
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly1"] = hams_poly1
    res["hamiltonian_poly2"] = hams_poly2
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["contact_probability"] = contact_probs
    
    # extract polymer coordinates into a matrix
    poly1_coords = zeros(Int32, poly1.length, 4)
    extractCoords!(poly1, poly1_coords)
    res["final_coordinates_poly1"] = poly1_coords

    # extract polymer coordinates into a matrix
    poly2_coords = zeros(Int32, poly2.length, 4)
    extractCoords!(poly2, poly2_coords)
    res["final_coordinates_poly2"] = poly2_coords
    
    # extract bulk coordinates into a matrix
    bulk_coords = zeros(Int32, bulk.length, 4)
    extractCoords!(bulk, bulk_coords)
    res["final_coordinates_bulk"] = bulk_coords
    
    res["runtime_min"] = t_current
    
    return res
end
