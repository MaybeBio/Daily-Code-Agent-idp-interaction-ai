### physical parameters
# MARK: -Phys. params

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 42,
    nspecies = 2
)

# heterogeneous polymer with patterned regions
poly = Polymer(;
    length = 32,
    maxlength = 1250,
    initShape = "line",
    nMonomerTypes = 12,  # base, pattern
    initPos = [5, 21, 21],
    pattern = :uniform
)

initMonomers!(poly, defaultChemPot = 1.0)

# Transcription factor
tf = BulkMonomer(;
    length = 1,
    initShape = "uniform"
)

# list containing all species
slist = [poly, tf]


### simulation parameters
# MARK: -Sim. params

blurb = """
How does enhancer prewetting affect the coalescence time and probability to nucleate
a droplet?

What we expect: one long enhancer collapses and nucleates a droplet very fast past 
the critical concentration, and many short enhancers gradually grow the droplet to 
a steady state size, attaching one by one. but the variation of local gradient relative 
to phi_infty should be the same as for a long enhancer, at steady state. 
it is just less stable and maybe it takes longer to get there?

swept parameters: concentration
output metrics: pairwise distances, centers of gravity, concentration profile of droplet
"""

init_from_checkpoint = false

t_max = 1    # runtime in minutes (initial test)
step = 10000
smallstep = step

vals = range(0.005, stop=0.03, length=32)  # parallelized parameters (concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

alphas = [0.01]    # concentration
betas = [-0.4]      # pattern-tf interaction      
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    poly, tf = slist

    base = poly.monomers[1]

    headPattern1 = poly.monomers[2]
    mainPattern1 = poly.monomers[3]
    tailPattern1 = poly.monomers[4]

    mainPattern2 = poly.monomers[5]
    tailPattern2 = poly.monomers[6]

    mainPattern3 = poly.monomers[7]
    tailPattern3 = poly.monomers[8]

    mainPattern4 = poly.monomers[9]
    tailPattern4 = poly.monomers[10]

    mainPattern5 = poly.monomers[11]
    tailPattern5 = poly.monomers[12]

    concentration = alpha
    temperature = 0.66
    
    # polymer self-interaction

    for (key, pat) in poly.monomers
        intPot!(pat, 0.5)
        intPot!(base, pat, 0.5)
    end
    
    # enhancers
    intPot!(mainPattern1, tf, beta)    # pattern-tf
    intPot!(mainPattern2, tf, beta)    # pattern-tf
    intPot!(mainPattern3, tf, beta)    # pattern-tf
    intPot!(mainPattern4, tf, beta)    # pattern-tf
    intPot!(mainPattern5, tf, beta)    # pattern-tf
    
    # tf (type 3)
    intPot!(tf, -1.0 / temperature)    # tf self-interaction

    tf.length = Int32(round(l.boxlength^3 * concentration))

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
                 (tf, tf.intIdx)]

    nmoves = 0
    niter = 0
    minsize = 10

    radii = []
    lengths = []
    hams_poly = []
    hams_tf = []
    max_droplet_size = []
    accs = []
    pattern_centers = []
    pattern_dists = []
    normalHist = zeros(150)
    maxDistance = 15.0
    total = 0.0
    dists_at_frame = Dict{Tuple, Float64}()
    centers_at_frame = Vector{Float64}(undef, 5)
    patterns_coalesced = false
    
    # initial values
    push!(radii, radius(poly))
    push!(lengths, poly.length)
    push!(hams_poly, poly.ham)
    push!(hams_tf, tf.ham)
    push!(accs, 0)

    # build polymer patterns
    patternVec = zeros(Int32, poly.maxlength)
    p1_start, p1_end = 201, 210
    p2_start, p2_end = 411, 420
    p3_start, p3_end = 621, 630
    p4_start, p4_end = 831, 840
    p5_start, p5_end = 1041, 1050

    patternVec[1:p1_start-1] .= headPattern1.intIdx         # base type (40 beads)
    patternVec[p1_start:p1_end] .= mainPattern1.intIdx      # pattern 1 (8 beads)
    patternVec[p1_end+1:p2_start-1] .= tailPattern1.intIdx

    patternVec[p2_start:p2_end] .= mainPattern2.intIdx
    patternVec[p2_end+1:p3_start-1] .= tailPattern2.intIdx

    patternVec[p3_start:p3_end] .= mainPattern3.intIdx
    patternVec[p3_end+1:p4_start-1] .= tailPattern3.intIdx
        
    patternVec[p4_start:p4_end] .= mainPattern4.intIdx
    patternVec[p4_end+1:p5_start-1] .= tailPattern4.intIdx

    patternVec[p5_start:p5_end] .= mainPattern5.intIdx
    patternVec[p5_end+1:poly.maxlength-7] .= tailPattern5.intIdx

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

                    elseif s.intIdx == tf.intIdx
                        push!(hams_tf, getHam(s))
                    end

                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii, radius(poly))
                    push!(lengths, poly.length)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))
                end

                # record pairwise distances between patterns
                # pattern_poly = find_patterns(poly)
                # dists_at_frame = Dict{Tuple, Float64}()
                # centers_at_frame = Matrix{Float64}(undef, 3, 5)

                # for (key1, monomerIdx1) in enumerate([3,5,7,9,11])
                #     pattern1 = poly.monomers[monomerIdx1]
                #     p1_start, p1_end = pattern_poly[pattern1.intIdx][1]

                    # centers_at_frame[:,key1] = centerMassPattern(poly, pattern1.intIdx)

                #     for (key2, monomerIdx2) in enumerate([3,5,7,9,11])
                #         if monomerIdx2 > monomerIdx1
                #             pattern2 = poly.monomers[monomerIdx2]
                #             p2_start, p2_end = pattern_poly[pattern2.intIdx][1]

                #             _, _, dist = get_pattern_centers(poly, 
                #                         [(p1_start, p1_end),(p2_start, p2_end)])

                #             # store distances in a dictionary
                #             push!(dists_at_frame, (key1, key2) => dist)
                #         end
                #     end
                # end
                
                # push!(pattern_dists, dists_at_frame)
                # push!(pattern_centers, centers_at_frame)

                # if maximum(values(dists_at_frame)) < 3.0
                #     patterns_coalesced = true
                # end

                # density profile
                if true
                    cm = centerMassPattern(poly, mainPattern1.intIdx)
                    histCount = updateNormalHistogram!(tf, normalHist, cm, maxDistance)
                    total += histCount
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
    res["pattern_centers"] = pattern_centers
    res["pattern_dists"] = pattern_dists
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_tf"] = hams_tf
    res["acceptance_rate"] = accs
    res["normal_density"] = normalHist
    res["normalization"] = total
    
    # extract polymer coordinates into a matrix
    poly_coords = zeros(Int32, poly.length, 4)
    extractCoords!(poly, poly_coords)
    res["final_coordinates_poly"] = poly_coords
    
    # extract tf coordinates into a matrix
    tf_coords = zeros(Int32, tf.length, 4)
    extractCoords!(tf, tf_coords)
    res["final_coordinates_tf"] = tf_coords
    
    res["runtime_min"] = t_current
    
    return res
end

