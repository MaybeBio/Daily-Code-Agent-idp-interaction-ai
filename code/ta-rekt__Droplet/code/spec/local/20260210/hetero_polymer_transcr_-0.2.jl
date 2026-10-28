### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 42,
    nspecies = 4  # base (1), pattern (2), tf (3), rnap (4)
)

# heterogeneous polymer with patterned regions
poly = Polymer(;
    length = 32,
    maxlength = 150,
    initShape = "line",
    nMonomerTypes = 2,  # base, pattern
    initPos = [5, 21, 21],
    pattern = :uniform
)

initMonomers!(poly, defaultChemPot = 1.0)

# Transcription factor
tf = BulkMonomer(;
    length = 1,
    initShape = "uniform"
)

# RNA Polymerase
rnap = BulkMonomer(;
    length = 1,
    initShape = "uniform"
)

# list containing all species
slist = [poly, tf, rnap]


### simulation parameters
# MARK: -Sim. params

blurb = """
How does polymer prewetting and crumpling affect the contact probability of rnap
monomers at the transcription initiation site?
One long pattern at the polymer's head (50 beads), followed by a transcription initiation
site (10 beads). The strong pattern-tf interaction prewets the pattern and crumples
the polymer. TF-RNAP coupling results in increased contacts_tf at TI site upon droplet formation
swept parameters: concentration, tf-rnap coupling strength
output metrics: rnap contact probability at TI site
"""

t_max = 480    # runtime in minutes (initial test)
step = 10000

vals = range(0.005, stop=0.03, length=32)  # parallelized parameters (concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [0.015]    # concentration
betas = [-0.2]    # tf-rnap coupling energy
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta, blurb = args

    poly, tf, rnap = slist

    base = poly.monomers[1]
    pattern = poly.monomers[2]

    concentration = alpha
    tf_rnap_intpot = beta
    temperature = 0.66
    pattern_tf_intpot = -0.5
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base, 0.5)              # base self-interaction
    intPot!(base, pattern, 0.5)    # base-pattern
    intPot!(base, tf, 0.0)        # base-tf (weak)
    intPot!(base, rnap, -0.1)      # base-rnap (weak)
    
    # pattern (type 2)
    intPot!(pattern, 0.5)           # pattern self-interaction
    intPot!(pattern, tf, pattern_tf_intpot)    # pattern-tf (strong)
    intPot!(pattern, rnap, 0.0)      # pattern-rnap (weak)
    
    # tf (type 3)
    intPot!(tf, -1.0 / temperature)    # tf self-interaction

    # rnap (type 4)
    intPot!(rnap, 0.0)
    intPot!(rnap, tf, tf_rnap_intpot)

    tf.length = Int32(round(l.boxlength^3 * concentration))
    rnap.length = Int32(round(l.boxlength^3 * 0.02))

    updateMonomers!(poly, base.intIdx)

    # initialize lattice and species
    init!(l, [poly, tf, rnap])

    # Initialize deltaHP for both polymers (needed for multi-polymer interactions)
    poly.deltaHP = 0.0

    # moving species
    movingsps = [(poly, tf.intIdx + 1),
                 (rnap, rnap.intIdx + 2),
                 (tf, tf.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1
    minsize = 10

    radii = []
    lengths = []
    hams_poly = []
    hams_tf = []
    hams_rnap = []
    max_droplet_size = []
    accs = []
    binding_activity = []
    transcription_activity = []
    
    # initial values
    push!(radii, radius(poly))
    push!(lengths, poly.length)
    push!(hams_poly, poly.ham)
    push!(hams_tf, tf.ham)
    push!(hams_rnap, rnap.ham)
    push!(accs, 0)
    push!(binding_activity, 0.0)
    push!(transcription_activity, 0.0)

    # build polymer patterns
    patternVec = zeros(Int32, poly.maxlength)
    p_start, p_end = 1, 50
    transcr_init_site_length = 10

    patternVec[1:p_start-1] .= base.intIdx         # base type (40 beads)
    patternVec[p_start:p_end] .= pattern.intIdx      # pattern 1 (8 beads)
    patternVec[p_end+1:poly.maxlength] .= base.intIdx

    growPolymer!(poly)
    updatePattern!(l, poly, patternVec)

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

                    ccs = cc(tf)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    droplet_sizes = map(x -> length(x), droplets)

                    # compute polymer hamiltonian (stored in poly, not monomers)
                    if i == tf.intIdx + 1
                        push!(hams_poly, poly.ham)

                    elseif s.intIdx == rnap.intIdx
                        push!(hams_rnap, getHam(s))

                    elseif s.intIdx == tf.intIdx
                        push!(hams_tf, getHam(s))
                    end

                    pattern_poly = find_patterns(poly)
                    p_start, p_end = pattern_poly[pattern.intIdx][1]

                    # compute contact probability between pattern and tf
                    contacts_tf = 0
                    contacts_rnap = 0
                    
                    # get tail key for pattern
                    if haskey(poly.patternTails, pattern.intIdx)
                        tail_key = poly.patternTails[pattern.intIdx]
                        
                        # walk forward right after tail to count nearby contacts_tf
                        current_key = poly.pos[tail_key][5]
                        count = 0
                        while count < transcr_init_site_length
                            # check if this bead is pattern type
                            if poly.pos[current_key][6] == base.intIdx
                                # count nearest neighbors with rnap
                                contacts_tf += nnp(pattern, tf, poly.pos[current_key][1:3]) / 6
                                contacts_rnap += nnp(pattern, rnap, poly.pos[current_key][1:3]) / 6
                                # move to previous bead
                                current_key = poly.pos[current_key][5]
                                count += 1
                            else
                                break  # reached end of pattern region
                            end
                        end
                    end
                    
                    push!(binding_activity, contacts_tf / transcr_init_site_length)
                    push!(transcription_activity, contacts_rnap / transcr_init_site_length)
                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii, radius(poly))
                    push!(lengths, poly.length)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))

                    frame += 1
                end
            end
        end
    end

    # save results
    res["blurb"] = blurb
    res["species_name"] = "poly"
    res["nframes"] = length(radii)
    res["step"] = step
    res["temperature"] = temperature
    res["concentration"] = concentration
    res["interaction_matrix"] = tf.intMat
    res["boxlength"] = l.boxlength
    res["length_poly"] = lengths
    res["length_tf"] = tf.length
    res["radius_of_gyration"] = radii
    res["pattern_location"] = p_start
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_tf"] = hams_tf
    res["hamiltonian_rnap"] = hams_rnap
    res["acceptance_rate"] = accs
    res["binding_activity"] = binding_activity
    res["transcription_activity"] = transcription_activity
    
    # extract polymer coordinates into a matrix
    poly_coords = zeros(Int32, poly.length, 4)
    extractCoords!(poly, poly_coords)
    res["final_coordinates_poly"] = poly_coords
    
    # extract tf coordinates into a matrix
    tf_coords = zeros(Int32, tf.length, 4)
    extractCoords!(tf, tf_coords)
    res["final_coordinates_tf"] = tf_coords
    
    # extract tf coordinates into a matrix
    rnap_coords = zeros(Int32, rnap.length, 4)
    extractCoords!(rnap, rnap_coords)
    res["final_coordinates_rnap"] = rnap_coords
    
    res["runtime_min"] = t_current
    
    return res
end
