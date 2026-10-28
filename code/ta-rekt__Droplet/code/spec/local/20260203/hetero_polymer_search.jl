"""
Description:

Understanding the effects of prewetting on the time it takes two strongly
interacting regions of a polymer to come into contact.
"""

### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 42,
    nspecies = 3  # base (1), pattern1 (2), pattern2 (3), bulk (4)
)

# heterogeneous polymer with patterned regions
poly1 = Polymer(;
    length = 32,
    maxlength = 512,
    initShape = "line",
    nMonomerTypes = 2,  # base, pattern1
    initPos = [5, 12, 12],
    pattern = :uniform
)

# heterogeneous polymer with patterned regions
poly2 = Polymer(;
    length = 32,
    maxlength = 512,
    initShape = "line",
    nMonomerTypes = 2,  # base, pattern1
    initPos = [5, 30, 30],
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

t_max = 600    # runtime in minutes (initial test)
step = 10000
t_burn_in = 0.8  # 12 seconds

vals = range(0.005, stop=0.03, length=32)  # parallelized parameters (concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [0.03]    # concentration
betas = [0.66]    # temperature
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    poly1, poly2, bulk = slist

    base1 = poly1.monomers[1]
    base2 = poly2.monomers[1]
    pattern1 = poly1.monomers[2]
    pattern2 = poly2.monomers[2]
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base1, 0.5)              # base self-interaction
    intPot!(base1, pattern1, 0.5)    # base-pattern1
    intPot!(base1, bulk, 0.0)        # base-bulk (weak)
    
    # pattern1 (type 2)
    intPot!(pattern1, 0.5)           # pattern1 self-interaction
    intPot!(pattern1, bulk, -0.5)    # pattern1-bulk (strong)

    # base (type 2)
    intPot!(base2, 0.5)              # base self-interaction
    intPot!(base2, pattern2, 0.5)    # base-pattern1
    intPot!(base2, bulk, 0.0)        # base-bulk (weak)
    
    # pattern1 (type 2)
    intPot!(pattern2, 0.5)           # pattern1 self-interaction
    intPot!(pattern2, bulk, -0.5)    # pattern1-bulk (strong)

    intPot!(pattern1, pattern2, 0.0)    # inter-pattern attraction
    
    # bulk (type 4)
    intPot!(bulk, -1.0 / beta)              # bulk self-interaction

    bulk.length = Int32(round(l.boxlength^3 * alpha))

    updateMonomers!(poly1, base1.intIdx)
    updateMonomers!(poly2, base2.intIdx)

    # initialize lattice and species
    init!(l, [poly1, poly2, bulk])

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
    # density_lat = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    total = 0.0
    
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
    p1_start, p1_end = 251, 266
    p2_start, p2_end = 251, 266

    patternVec1[1:p1_start-1] .= base1.intIdx         # base type (40 beads)
    patternVec1[p1_start:p1_end] .= pattern1.intIdx      # pattern 1 (8 beads)
    patternVec1[p1_end+1:poly1.maxlength] .= base1.intIdx

    patternVec2[1:p1_start-1] .= base2.intIdx         # base type (40 beads)
    patternVec2[p1_start:p1_end] .= pattern2.intIdx      # pattern 1 (8 beads)
    patternVec2[p1_end+1:poly2.maxlength] .= base2.intIdx

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
                    p1_start, p1_end = pattern_poly1[pattern1.intIdx][1]

                    pattern_poly2 = find_patterns(poly2)
                    p2_start, p2_end = pattern_poly2[pattern2.intIdx][1]

                    # compute contact probability between pattern1 and bulk
                    pattern_length = p1_end - p1_start + 1
                    contacts = 0
                    
                    # get tail key for pattern1
                    if haskey(poly1.patternTails, pattern1.intIdx)
                        tail_key = poly1.patternTails[pattern1.intIdx]
                        
                        # walk backwards from tail to count contacts
                        current_key = tail_key
                        while current_key != 0
                            # check if this bead is pattern1 type
                            if poly1.pos[current_key][6] == pattern1.intIdx
                                # count nearest neighbors with bulk
                                contacts += nnp(pattern1, bulk, poly1.pos[current_key][1:3]) / 4
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

                    # if t_current >= t_burn_in
                    #     for (_, val) in bulk.pos
                    #         density_lat[val[1], val[2], val[3]] += 1.0
                    #         total += 1.0
                    #     end
                    # end

                    frame += 1
                end
            end
        end
    end

    # density_lat = density_lat ./ total

    # save results
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
    res["pattern1_location"] = p1_start
    res["pattern2_location"] = p2_start
    res["max_droplet_size"] = max_droplet_size
    # res["density_lattice"] = density_lat
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

