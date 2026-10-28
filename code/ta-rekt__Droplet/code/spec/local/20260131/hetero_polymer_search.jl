"""
Description:

Understanding the effects of prewetting on a strongly interacting region
of a free-floating polymer in the extended phase.
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
    maxlength = 256,
    initShape = "line",
    nMonomerTypes = 2,  # base, pattern1
    initPos = [5, 12, 12],
    pattern = :custom
)

# heterogeneous polymer with patterned regions
poly2 = Polymer(;
    length = 32,
    maxlength = 256,
    initShape = "line",
    nMonomerTypes = 2,  # base, pattern1
    initPos = [5, 30, 30],
    pattern = :custom
)

pattern = zeros(Int, 256)
p1_start, p1_end = 124, 134
p2_start, p2_end = 124, 134
pattern[1:p1_start-1] .= 1         # base type (40 beads)
pattern[p1_start:p1_end] .= 2      # pattern 1 (8 beads)
pattern[p1_end+1:256] .= 1

poly1.customPattern = ones(Int, 256)
poly2.customPattern = ones(Int, 256)

initMonomers!(poly1, defaultChemPot = 1.0)
initMonomers!(poly2, defaultChemPot = 1.0)

base1 = poly1.monomers[1]
base2 = poly2.monomers[1]
pattern1 = poly1.monomers[2]
pattern2 = poly2.monomers[2]

# bulk species
bulk = FastBulkMonomer(;
    length = 1,
    initShape = "uniform"
)

# list containing all species
slist = [base1, pattern1, base2, pattern2, bulk]


### simulation parameters
# MARK: -Sim. params

t_max = 720    # runtime in minutes (initial test)
step = 10000
t_burn_in = 0.8  # 12 seconds

vals = range(0.001, stop=0.05, length=32)  # parallelized parameters (concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [0.03]    # concentration
betas = [0.64]    # temperature
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0
    p1_reached_maxlength = false
    p2_reached_maxlength = false

    # species list and sweep parameters
    slist, alpha, beta = args

    base1, pattern1, base2, pattern2, bulk = slist
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base1, 0.5)              # base self-interaction
    intPot!(base1, pattern1, 0.5)    # base-pattern1
    intPot!(base1, bulk, 0.0)       # base-bulk (weak)
    
    # pattern1 (type 2)
    intPot!(pattern1, 0.5)           # pattern1 self-interaction
    intPot!(pattern1, bulk, -2.0)    # pattern1-bulk (strong)

    # base (type 2)
    intPot!(base2, 0.5)              # base self-interaction
    intPot!(base2, pattern2, 0.5)    # base-pattern1
    intPot!(base2, bulk, 0.0)       # base-bulk (weak)
    
    # pattern1 (type 2)
    intPot!(pattern2, 0.5)           # pattern1 self-interaction
    intPot!(pattern2, bulk, -2.0)    # pattern1-bulk (strong)
    
    # bulk (type 4)
    intPot!(bulk, -1.0 / beta)              # bulk self-interaction

    bulk.length = Int32(round(l.boxlength^3 * alpha))

    # initialize lattice and species
    init!(l, [base1, pattern1, base2, pattern2, bulk])

    # Initialize deltaHP for both polymers (needed for multi-polymer interactions)
    poly1.deltaHP = 0.0
    poly2.deltaHP = 0.0

    # moving species
    movingsps = [(base2, base2.intIdx), 
                 (pattern2, pattern2.intIdx), 
                 (base1, base1.intIdx), 
                 (pattern1, pattern1.intIdx), 
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

    # contact matrix: tracks contact frequency between monomer pairs
    contact_matrix = zeros(Int, poly1.maxlength, poly1.maxlength)
    contact_threshold = 2.0  # distance cutoff for contact
    
    # initial values
    push!(radii1, radius(poly1))
    push!(radii2, radius(poly2))
    push!(lengths1, poly1.length)
    push!(lengths2, poly2.length)
    _, _, initial_dist = get_pattern_centers(poly1, poly2, 
                    [(p1_start, p1_end),(p1_start, p1_end)])
    push!(pattern_dists, initial_dist)
    push!(hams_poly1, poly1.ham)
    push!(hams_poly2, poly2.ham)
    push!(hams_bulk, bulk.ham)
    push!(accs, 0)
    
    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1

                # record data
                if mod(nmoves, step) == 0

                    pattern = zeros(Int, 256)
                    p1_start, p1_end = 124, 134
                    p2_start, p2_end = 124, 134
                    pattern[1:p1_start-1] .= 1         # base type (40 beads)
                    pattern[p1_start:p1_end] .= 2      # pattern 1 (8 beads)
                    pattern[p1_end+1:256] .= 1
                    
                    if poly1.length == poly1.maxlength && !p1_reached_maxlength
                        poly1.customPattern = pattern
                        initMonomers!(poly1, defaultChemPot = 1.0)
                        init!(l, [base1, pattern1, base2, pattern2, bulk])
                        p1_reached_maxlength = true

                        pattern_poly1 = find_patterns(poly1)
                        p1_start, p1_end = pattern_poly1[2][1]
                    end

                    if poly2.length == poly2.maxlength && !p2_reached_maxlength
                        poly2.customPattern = pattern
                        initMonomers!(poly2, defaultChemPot = 1.0)
                        init!(l, [base1, pattern1, base2, pattern2, bulk])
                        p2_reached_maxlength = true

                        pattern_poly2 = find_patterns(poly2)
                        p2_start, p2_end = pattern_poly2[2][1]
                    end

                    ham!(s)

                    ccs = cc(bulk)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    droplet_sizes = map(x -> length(x), droplets)

                    # compute polymer hamiltonian (stored in poly, not monomers)
                    if s.intIdx in [base1.intIdx, pattern1.intIdx]
                        push!(hams_poly1, poly1.ham)

                    elseif s.intIdx in [base2.intIdx, pattern2.intIdx]
                        push!(hams_poly2, poly2.ham)

                    elseif s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    _, _, dist = get_pattern_centers(poly1, poly2, 
                                    [(p1_start, p1_end),(p2_start, p2_end)])
                    push!(pattern_dists, dist)

                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii1, radius(poly1))
                    push!(lengths1, poly1.length)
                    push!(radii2, radius(poly2))
                    push!(lengths2, poly2.length)

                    # update contact matrix
                    compute_contact_matrix!(poly1, contact_matrix, contact_threshold)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))

                    frame += 1
                end
            end
        end
    end

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
    res["hamiltonian_poly1"] = hams_poly1
    res["hamiltonian_poly2"] = hams_poly2
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["contact_matrix"] = contact_matrix
    
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

