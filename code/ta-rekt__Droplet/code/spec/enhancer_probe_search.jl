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
    maxlength = 300,
    initShape = "line",
    nMonomerTypes = 6,  # base, pattern
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
How does the binding time between a probe and an enhancer that attract depend
on separation (in #beads) and on the interaction strength? I would like to test the
hypothesis that the probe return time is exponentially distributed, as well as the
bound time, in the absence of nucleation on the enhancer. this would make the 
probe-enhancer system effectively a two-state system, where the binding probability
is a Hill function of coupling energy.

With nucleation however, it might be that bound times are sub-exponential. That would
mean that you can get sharper binding curves than Hill. Let's see.

swept parameters: tf_concentration
output metrics: RMS probe-enhancer distance, binding activity, bound/unbound times
local_concentration at probe
"""

init_from_checkpoint = false

t_max = 120    # runtime in minutes (initial test)
step = 10000
smallstep = 2000

vals = exp.(range(-7, stop=-3.5, length=32))  # parallelized parameters (tf_concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

alphas = [exp(-5.0)]   # tf concentration
betas = [0]           # nothing
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    poly, tf, rnap = slist
    tf_concentration = alpha
    rnap_concentration = 0.02

    base = poly.monomers[1]

    headPattern1 = poly.monomers[2]
    mainPattern1 = poly.monomers[3]
    tailPattern1 = poly.monomers[4]
    mainPattern2 = poly.monomers[5]
    tailPattern2 = poly.monomers[6]

    probe_dist = 20
    enhancer_len = 10
    tf_rnap_coupling = -7.75
    temperature = 0.66
    
    # polymer self-interaction
    for (_, pat) in poly.monomers
        intPot!(pat, 0.5)
        intPot!(base, pat, 0.5)
    end
    
    # enhancers
    intPot!(mainPattern1, tf, -0.5)    # pattern-tf
    
    # tf (type 3)
    intPot!(tf, -1.0 / temperature)    # tf self-interaction

    # set two-way coupling
    couplingPot!(tf, mainPattern2, 0.0)
    couplingPot!(rnap, mainPattern2, -0.5)

    # Set three-way coupling
    couplingPot!(tf, rnap, mainPattern2, tf_rnap_coupling)   # TF only binds p2 when rnap bound

    tf.length = Int32(round(l.boxlength^3 * tf_concentration))
    rnap.length = Int32(round(l.boxlength^3 * rnap_concentration))

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

    lengths = []
    hams_poly = []
    hams_tf = []
    hams_rnap = []
    max_droplet_size = []
    accs = []
    binding_activity = []
    transcription_activity = []
    tf_rnap_activity = []
    probe_enhancer_activity = []
    local_concentration = []
    on_times = []
    off_times = []
    pattern_dists = []
    
    # transcription state tracking
    is_transcribing = false
    moves_in_state = 0
    
    # initial values
    push!(lengths, poly.length)
    push!(hams_poly, poly.ham)
    push!(hams_tf, tf.ham)
    push!(hams_rnap, rnap.ham)

    # build polymer patterns
    patternVec = zeros(Int32, poly.maxlength)
    p1_start, p1_end = 1, enhancer_len
    p2_start, p2_end = p1_end+1+probe_dist, p1_end+1+probe_dist
    transcr_init_site_length = 1

    patternVec[1:p1_start-1] .= headPattern1.intIdx         # base type (40 beads)
    patternVec[p1_start:p1_end] .= mainPattern1.intIdx      # pattern 1 (8 beads)
    patternVec[p1_end+1:p2_start-1] .= tailPattern1.intIdx

    patternVec[p2_start:p2_end] .= mainPattern2.intIdx
    patternVec[p2_end+1:poly.maxlength-1] .= tailPattern2.intIdx

    patternVec[poly.maxlength:poly.maxlength] .= base.intIdx

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

                    # push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(lengths, poly.length)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))
                end

                if mod(nmoves, smallstep) == 0

                    # compute contact probability at the mainPattern2 probe site
                    contacts_tf = 0
                    contacts_rnap = 0
                    contacts_tf_rnap = 0
                    contacts_enhancer = 0
                    local_tf = 0
                    
                    # get tail key for mainPattern2 — the tail IS the probe bead
                    if haskey(poly.patternTails, mainPattern2.intIdx)
                        tail_key = poly.patternTails[mainPattern2.intIdx]
                        
                        # walk the mainPattern2 region itself (not past it)
                        current_key = tail_key
                        count = 0
                        while count < transcr_init_site_length
                            
                            # check if this bead is actually mainPattern2
                            if poly.pos[current_key][6] == mainPattern2.intIdx
                                
                                # measure species occupancy at this probe site
                                pos = poly.pos[current_key][1:3]
                                contacts_tf += Int(tf.latS[pos[1], pos[2], pos[3]] > 0)
                                contacts_rnap += Int(rnap.latS[pos[1], pos[2], pos[3]] > 0)
                                contacts_tf_rnap += Int(tf.latS[pos[1], pos[2], pos[3]] > 0 && 
                                                        rnap.latS[pos[1], pos[2], pos[3]] > 0)

                                local_tf += nnp(mainPattern2, tf, pos) / 6

                                mlatS = getLatS(mainPattern2)
                                contacts_enhancer += Int(mlatS[pos[1], pos[2], pos[3]] 
                                                        == mainPattern1.intIdx)
                                
                                # move backward along polymer (toward head)
                                current_key = poly.pos[current_key][4]
                                count += 1
                            else
                                break  # left the mainPattern2 region
                            end
                        end
                    end
                    
                    current_transcription_activity = contacts_rnap / transcr_init_site_length

                    push!(binding_activity, contacts_tf / transcr_init_site_length)
                    push!(transcription_activity, current_transcription_activity)
                    push!(tf_rnap_activity, contacts_tf_rnap / transcr_init_site_length)
                    push!(probe_enhancer_activity, contacts_enhancer / transcr_init_site_length)
                    push!(local_concentration, local_tf / transcr_init_site_length)

                    # measure distance between probe and TF binding site
                    _, _, dist = get_pattern_centers(poly, 
                                    [(p1_start, p1_end),(p2_start, p2_end)])
                    push!(pattern_dists, dist)
                    
                    # track on/off transitions
                    was_transcribing = is_transcribing
                    # is_transcribing = current_transcription_activity >= transcription_threshold
                    is_transcribing = contacts_enhancer > 0 
                    
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
    res["nframes"] = length(accs)
    res["step"] = step
    res["smallstep"] = smallstep
    res["temperature"] = temperature
    res["tf_concentration"] = tf_concentration
    res["rnap_concentration"] = rnap_concentration
    res["interaction_matrix"] = tf.intMat
    res["boxlength"] = l.boxlength
    res["length_poly"] = lengths
    res["length_tf"] = tf.length
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_tf"] = hams_tf
    res["hamiltonian_rnap"] = hams_rnap
    res["acceptance_rate"] = accs
    res["binding_activity"] = binding_activity
    res["transcription_activity"] = transcription_activity
    res["tf_rnap_activity"] = tf_rnap_activity
    res["probe_enhancer_activity"] = probe_enhancer_activity
    res["local_concentration"] = local_concentration
    res["on_times"] = on_times
    res["off_times"] = off_times
    res["pattern_distances"] = pattern_dists
    
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

