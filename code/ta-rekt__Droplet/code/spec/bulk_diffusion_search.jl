### physical parameters
# MARK: -Phys. params

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 200,
    nspecies = 2
)

# Transcription factor
tf = BulkMonomer(;
    length = 1,
    initShape = "line",
    initPos = [100, 100, 100]
    # initParams = 1
)

# RNA Polymerase
rnap = BulkMonomer(;
    length = 1,
    initShape = "line",
    initPos = [100, 100, 100]
    # initParams = 1
)

# list containing all species
slist = [tf, rnap]


### simulation parameters
# MARK: -Sim. params

blurb = """
How long does it take a diffusing particle to hit a solid absorbing sphere, 
when it is initialized at a certain distance in the lattice. Can have two 
versions, with and without bulk. To compare with calculations of search time 
and bound time distributions.

Also look at the single and pairwise MSD scaling.

swept parameters: none (identical repeats)
output metrics: inter-particle distance, bound/unbound times
"""

init_from_checkpoint = false

t_max = 0.05    # runtime in minutes (initial test)
maxiters = 100000
step = 1
smallstep = step

#vals = range(0.005, stop=0.03, length=32)  # parallelized parameters (tf_concentration)
vals = [0]
alphas = zeros(1)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

alphas = [0]   # tf concentration
betas = [0]           # nothing
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    tf, rnap = slist

    temperature = 0.66
    
    # enhancers
    # intPot!(tf, rnap, -2.0)

    tf.length = 1
    rnap.length = 1

    # initialize lattice and species
    if init_from_checkpoint
        checkpoint_at_param = string(checkpoint_path, "_keya", keya)
        restore!(checkpoint_at_param, l, slist)
    else
        init!(l, slist)
    end

    # moving species
    movingsps = [(rnap, tf.intIdx+1)
        (tf, tf.intIdx)]

    nmoves = 0
    niter = 0
    minsize = 10

    hams_tf = []
    hams_rnap = []
    accs = []
    binding_activity = []
    on_times = []
    off_times = []
    pattern_dists = []
    tf_pos = []
    
    # transcription state tracking
    is_bound = false
    moves_in_state = 0
    
    # initial values
    push!(hams_tf, tf.ham)
    push!(hams_rnap, rnap.ham)

    # run simulation
    while niter <= maxiters
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
                    if s.intIdx == rnap.intIdx
                        push!(hams_rnap, getHam(s))

                    elseif s.intIdx == tf.intIdx
                        push!(hams_tf, getHam(s))
                    end

                    push!(accs, nmoves / (niter * length(movingsps)))
                end

                if mod(nmoves, smallstep) == 0

                    # compute contact probability at the mainPattern2 probe site
                    contacts_tf_rnap = 0
                    local_tf = 0
            
                    # measure species occupancy at this probe site
                    pos = tf.pos[1][1:3]
                    contacts_tf_rnap += Int(tf.latS[pos[1], pos[2], pos[3]] > 0 && 
                                            rnap.latS[pos[1], pos[2], pos[3]] > 0)

                    push!(binding_activity, contacts_tf_rnap)

                    # measure distance between RNAP and TF particles
                    dist = sum((pos .- rnap.pos[1][1:3]).^2)^(0.5)
                    push!(pattern_dists, dist)

                    # track position
                    push!(tf_pos, pos)
                    
                    # track on/off transitions
                    was_bound = is_bound
                    is_bound = contacts_tf_rnap > 0 
                    
                    if is_bound != was_bound
                        # state transition occurred
                        if was_bound
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
    res["tf_concentration"] = tf.length
    res["rnap_concentration"] = rnap.length
    res["interaction_matrix"] = tf.intMat
    res["boxlength"] = l.boxlength
    res["length_poly"] = 0
    res["probe_distance"] = 0
    res["enhancer_length"] = 0
    res["tf_rnap_coupling"] = 0
    res["length_tf"] = tf.length
    res["max_droplet_size"] = 0
    res["hamiltonian_poly"] = 0
    res["hamiltonian_tf"] = 0
    res["hamiltonian_rnap"] = 0
    res["acceptance_rate"] = accs
    res["binding_activity"] = binding_activity
    res["transcription_activity"] = []
    res["tf_rnap_activity"] = []
    res["probe_enhancer_activity"] = []
    res["bound_times"] = []
    res["on_times"] = on_times
    res["off_times"] = off_times
    res["pattern_distances"] = pattern_dists
    res["probe_positions"] = tf_pos
    
    # extract polymer coordinates into a matrix
    # poly_coords = zeros(Int32, poly.length, 4)
    # extractCoords!(poly, poly_coords)
    res["final_coordinates_poly"] = []
    
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


