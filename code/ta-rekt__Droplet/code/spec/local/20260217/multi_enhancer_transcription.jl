### physical parameters
# MARK: -Phys. params

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 42,
    nspecies = 2  # base (1), pattern (2), tf (3), rnap (4)
)

# heterogeneous polymer with patterned regions
poly = Polymer(;
    length = 32,
    maxlength = 1200,
    initShape = "line",
    nMonomerTypes = 8,  # base, pattern
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
How does enhancer prewetting and crumpling affect the contact probability of rnap
monomers at the transcription initiation site?

Two enhancers and a promoter, the latter followed by a transcription initiation
site (10 beads). The strong pattern-tf interaction prewets the enhancers and crumples
the polymer. TF-RNAP coupling results in increased contacts_tf at TI site upon droplet formation
swept parameters: concentration
output metrics: rnap contact time-series at TI site
"""

init_from_checkpoint = true

t_max = 720    # runtime in minutes (initial test)
step = 10000
smallstep = step

vals = range(0.001, stop=0.035, length=32)  # parallelized parameters (concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [0.01]    # concentration
betas = [-0.0]    # tf-rnap coupling energy
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    poly, tf, rnap = slist

    base = poly.monomers[1]

    headPattern1 = poly.monomers[2]
    mainPattern1 = poly.monomers[3]
    tailPattern1 = poly.monomers[4]

    mainPattern2 = poly.monomers[5]
    tailPattern2 = poly.monomers[6]

    mainPattern3 = poly.monomers[7]
    tailPattern3 = poly.monomers[8]

    concentration = alpha
    tf_rnap_intpot = beta
    temperature = 0.7
    transcription_threshold = 0.1  # threshold for transcription activity
    
    # polymer self-interaction
    intPot!(base, 0.5)                   # base self-interaction
    intPot!(headPattern1, 0.5)           # mainPattern1 self-interaction
    intPot!(mainPattern1, 0.5)           # mainPattern1 self-interaction
    intPot!(tailPattern1, 0.5)           # mainPattern1 self-interaction
    intPot!(mainPattern2, 0.5)           # mainPattern1 self-interaction
    intPot!(tailPattern2, 0.5)           # mainPattern1 self-interaction
    intPot!(mainPattern3, 0.5)           # mainPattern1 self-interaction
    intPot!(tailPattern3, 0.5)           # mainPattern1 self-interaction
    intPot!(base, headPattern1, 0.5)    # base-mainPattern1
    intPot!(base, mainPattern1, 0.5)    # base-mainPattern1
    intPot!(base, tailPattern1, 0.5)    # base-mainPattern1
    intPot!(base, mainPattern2, 0.5)    # base-mainPattern1
    intPot!(base, tailPattern2, 0.5)    # base-mainPattern1
    intPot!(base, mainPattern3, 0.5)    # base-mainPattern1
    intPot!(base, tailPattern3, 0.5)    # base-mainPattern1

    # rnap-transcription initiation site interaction
    intPot!(tailPattern3, rnap, -0.0)      # (weak)
    
    # enhancers
    intPot!(mainPattern1, tf, -0.3)    # pattern-tf (weak)
    intPot!(mainPattern2, tf, -0.7)    # pattern-tf (med)
    intPot!(mainPattern3, tf, -1.2)    # pattern-tf (strong)
    
    # tf (type 3)
    intPot!(tf, -1.0 / temperature)    # tf self-interaction

    # rnap (type 4)
    intPot!(rnap, tf, tf_rnap_intpot)

    tf.length = Int32(round(l.boxlength^3 * concentration))
    #rnap.length = Int32(round(l.boxlength^3 * 0.02))
    rnap.length = 2

    updateMonomers!(poly, base.intIdx)

    # initialize lattice and species
    if init_from_checkpoint
        checkpoint_at_param = string(checkpoint_path, "_keya", keya)
        restore!(checkpoint_at_param, l, slist)
    else
        init!(l, slist)
    end

    # Initialize deltaHP for both polymers (needed for multi-polymer interactions)
    poly.deltaHP = 0.0

    # moving species
    movingsps = [(poly, tf.intIdx + 1),
                 (rnap, rnap.intIdx + 2),
                 (tf, tf.intIdx)]

    nmoves = 0
    niter = 0
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
    on_times = []
    off_times = []
    pattern1_position = []
    pattern2_position = []
    pattern3_position = []
    
    # transcription state tracking
    is_transcribing = false
    moves_in_state = 0
    
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
    p1_start, p1_end = 251, 350
    p2_start, p2_end = 551, 650
    p3_start, p3_end = 851, 950
    transcr_init_site_length = 10

    patternVec[1:p1_start-1] .= headPattern1.intIdx         # base type (40 beads)
    patternVec[p1_start:p1_end] .= mainPattern1.intIdx      # pattern 1 (8 beads)
    patternVec[p1_end+1:p2_start-1] .= tailPattern1.intIdx

    patternVec[p2_start-1:p2_end] .= mainPattern2.intIdx
    patternVec[p2_end+1:p3_start-1] .= tailPattern2.intIdx

    patternVec[p3_start-1:p3_end] .= mainPattern3.intIdx
    patternVec[p3_end+1:poly.maxlength-7] .= tailPattern3.intIdx

    patternVec[poly.maxlength-6:poly.maxlength] .= base.intIdx

    growPolymer!(poly, minAttempts = 900000)
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

                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii, radius(poly))
                    push!(lengths, poly.length)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))
                end

                if mod(nmoves, smallstep) == 0
                    pattern_poly = find_patterns(poly)
                    p1_start, p1_end = pattern_poly[mainPattern1.intIdx][1]
                    p2_start, p2_end = pattern_poly[mainPattern2.intIdx][1]
                    p3_start, p3_end = pattern_poly[mainPattern3.intIdx][1]

                    push!(pattern1_position, p1_start)
                    push!(pattern2_position, p2_start)
                    push!(pattern3_position, p3_start)

                    # compute contact probability between pattern and tf
                    contacts_tf = 0
                    contacts_rnap = 0
                    
                    # get tail key for pattern
                    if haskey(poly.patternTails, mainPattern3.intIdx)
                        tail_key = poly.patternTails[mainPattern3.intIdx]
                        
                        # walk forward right after tail to count nearby contacts_tf
                        current_key = poly.pos[tail_key][5]
                        count = 0
                        while count < transcr_init_site_length
                            # check if this bead is pattern type
                            if poly.pos[current_key][6] == tailPattern3.intIdx
                                # count nearest neighbors with rnap
                                contacts_tf += nnp(tailPattern3, tf, poly.pos[current_key][1:3]) / 6
                                contacts_rnap += nnp(tailPattern3, rnap, poly.pos[current_key][1:3]) / 6
                                # move to previous bead
                                current_key = poly.pos[current_key][5]
                                count += 1

                            elseif poly.pos[current_key][6] == base.intIdx
                                # count nearest neighbors with rnap
                                contacts_tf += nnp(base, tf, poly.pos[current_key][1:3]) / 6
                                contacts_rnap += nnp(base, rnap, poly.pos[current_key][1:3]) / 6
                                # move to previous bead
                                current_key = poly.pos[current_key][5]
                                count += 1
                            else
                                break  # reached end of pattern region
                            end
                        end
                    end
                    
                    push!(binding_activity, contacts_tf / transcr_init_site_length)
                    current_transcription_activity = contacts_rnap / transcr_init_site_length
                    push!(transcription_activity, current_transcription_activity)
                    
                    # track on/off transitions
                    was_transcribing = is_transcribing
                    is_transcribing = current_transcription_activity >= transcription_threshold
                    
                    if is_transcribing != was_transcribing
                        # state transition occurred
                        if was_transcribing
                            # transitioned from ON to OFF
                            push!(on_times, moves_in_state)
                        else
                            # transitioned from OFF to ON
                            push!(off_times, moves_in_state)
                        end
                        moves_in_state = 1
                    else
                        moves_in_state += 1
                    end
                end
            end
        end
    end

    # save results
    res["blurb"] = blurb
    res["species_name"] = "poly"
    res["nframes"] = length(radii)
    res["step"] = step
    res["smallstep"] = smallstep
    res["temperature"] = temperature
    res["concentration"] = concentration
    res["interaction_matrix"] = tf.intMat
    res["boxlength"] = l.boxlength
    res["length_poly"] = lengths
    res["length_tf"] = tf.length
    res["radius_of_gyration"] = radii
    res["pattern1_position"] = pattern1_position
    res["pattern2_position"] = pattern2_position
    res["pattern3_position"] = pattern3_position
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_tf"] = hams_tf
    res["hamiltonian_rnap"] = hams_rnap
    res["acceptance_rate"] = accs
    res["binding_activity"] = binding_activity
    res["transcription_activity"] = transcription_activity
    res["on_times"] = on_times
    res["off_times"] = off_times
    
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

