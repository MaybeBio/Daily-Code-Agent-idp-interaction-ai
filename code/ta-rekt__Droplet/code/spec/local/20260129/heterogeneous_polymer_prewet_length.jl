### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 64,
    nspecies = 3  # base (1), pattern1 (2), pattern2 (3), bulk (4)
)

# heterogeneous polymer with patterned regions
poly = Polymer(;
    length = 50,
    maxlength = 256,
    initShape = "line",
    nMonomerTypes = 2,  # base, pattern1
    initPos = [7, 32, 32],
    pattern = :custom
)

# define pattern: 256 beads with two patterned regions
pattern = zeros(Int, 256)
p1_start, p1_end = 41, 48
pattern[1:p1_start-1] .= 1       # base type (40 beads)
pattern[p1_start:p1_end] .= 2      # pattern 1 (8 beads)
pattern[p1_end+1:256] .= 1

poly.customPattern = ones(Int, 256)

initMonomers!(poly, defaultChemPot = 1.0)
base = poly.monomers[1]
pattern1 = poly.monomers[2]
# pattern2 = poly.monomers[3]

# bulk species
bulk = FastBulkMonomer(;
    length = 4096,
    initShape = "uniform"
)

# list containing all species
slist = [base, pattern1, bulk]


### simulation parameters
# MARK: -Sim. params

t_max = 480    # runtime in minutes (initial test)
step = 10000
t_burn_in = 0.8  # 12 seconds

vals = range{Int32}(8, stop=23, length=16)   # parallelized parameters (inter-pattern interaction)
alphas = zeros(64)
nrepeats = 4

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [1]
betas = [0.0]    # placeholder
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0
    reached_maxlength = false

    # species list and sweep parameters
    slist, alpha, beta = args

    base, pattern1, bulk = slist
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base, 0.5)               # base self-interaction
    intPot!(base, pattern1, 0.5)     # base-pattern1
    # intPot!(base, pattern2, 0.5)     # base-pattern2
    intPot!(base, bulk, -0.1)        # base-bulk (weak)
    
    # pattern1 (type 2)
    intPot!(pattern1, 0.5)           # pattern1 self-interaction
    # intPot!(pattern1, pattern2, -2.0)  # pattern1-pattern2 (SWEEP)
    intPot!(pattern1, bulk, -1.0)    # pattern1-bulk (strong)
    
    # pattern2 (type 3)
    # intPot!(pattern2, 0.5)           # pattern2 self-interaction
    # intPot!(pattern2, bulk, -1.0)    # pattern2-bulk (strong)
    
    # bulk (type 4)
    intPot!(bulk, -1.6)              # bulk self-interaction

    # initialize lattice and species
    init!(l, [base, pattern1, bulk])

    # moving species
    movingsps = [(base, base.intIdx), 
                 (pattern1, pattern1.intIdx), 
                 (bulk, bulk.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1
    minsize = 20

    radii = []
    lengths = []
    hams_poly = []
    hams_bulk = []
    max_droplet_size = []
    accs = []

    # contact matrix: tracks contact frequency between monomer pairs
    contact_matrix = zeros(Int, poly.maxlength, poly.maxlength)
    contact_threshold = 2.0  # distance cutoff for contact
    
    # initial values
    push!(radii, radius(poly))
    push!(lengths, poly.length)
    push!(hams_poly, poly.ham)
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

                    if poly.length == poly.maxlength && !reached_maxlength
                        pattern = zeros(Int, 256)
                        p1_start, p1_end = 31, 31+Int32(alpha)
                        pattern[1:p1_start-1] .= 1       # base type (40 beads)
                        pattern[p1_start:p1_end] .= 2      # pattern 1 (8 beads)
                        pattern[p1_end+1:256] .= 1

                        poly.customPattern = pattern
                        initMonomers!(poly, defaultChemPot = 1.0)
                        init!(l, [base, pattern1, bulk])
                        reached_maxlength = true
                    end

                    ham!(s)

                    ccs = cc(bulk)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    droplet_sizes = map(x -> length(x), droplets)

                    # compute polymer hamiltonian (stored in poly, not monomers)
                    if s.intIdx in [base.intIdx, pattern1.intIdx]
                        push!(hams_poly, poly.ham)

                    elseif s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii, radius(poly))
                    push!(lengths, poly.length)

                    # update contact matrix
                    compute_contact_matrix!(poly, contact_matrix, contact_threshold)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))

                    frame += 1
                end
            end
        end
    end

    # save results
    res["species_name"] = "poly"
    res["nframes"] = length(radii)
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["interaction_matrix"] = bulk.intMat
    res["boxlength"] = l.boxlength
    res["length_poly"] = lengths
    res["length_bulk"] = bulk.length
    res["radius_of_gyration"] = radii
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["contact_matrix"] = contact_matrix
    
    # extract polymer coordinates into a matrix
    poly_coords = zeros(Int32, poly.length, 4)
    extractCoords!(poly, poly_coords)
    res["final_coordinates_poly"] = poly_coords
    
    # extract bulk coordinates into a matrix
    bulk_coords = zeros(Int32, bulk.length, 4)
    extractCoords!(bulk, bulk_coords)
    res["final_coordinates_bulk"] = bulk_coords
    
    res["runtime_min"] = t_current
    
    return res
end
