"""
Description:

How do two bulk species (RNAP and TF) combine to increase the presence
of one at a site with strong affinity for the other
"""

### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 32,
    nspecies = 8
)

tf = FastBulk(;
    length = 32,
    initShape = "uniform"
)

rnap = FastBulk(;
    length = 32,
    initShape = "uniform"
)

### rod 1
dna1 = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [1,10,16]
)

enhancer = FastBulk(;
    length = 6,
    initShape = "line",
    initParams = 1,
    initPos = [13,10,16]
)

dna2 = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [19,10,16]
)


### rod 1
dna1_l = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [1,13,16]
)

promoter = FastBulk(;
    length = 6,
    initShape = "line",
    initParams = 1,
    initPos = [13,13,16]
)

dna2_l = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [19,13,16]
)

slist = [dna1, dna2, enhancer, dna1_l, dna2_l, promoter, tf, rnap]

### simulation parameters
# MARK: -Sim. params

t_max = 10    # runtime in minutes
step = 1000
t_burn_in = 0.8 * t_max

# minimum size of a droplet for detection
minsize = 5

vals = range(1300, stop=4000, step=180)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [512]   # number of tf particles
betas = [0.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)

### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    dna1, dna2, enhancer, dna1_l, dna2_l, promoter, tf, rnap = slist

    tf.length = alpha

    # self-interactions
    intPot!(tf, -1.0)
    intPot!(rnap, -1.0)

    # multi-species interactions
    intPot!(tf, dna1, -0.2)  
    intPot!(tf, dna2, -0.2)  
    intPot!(tf, promoter, -0.2)
    intPot!(tf, enhancer, -2.0)

    intPot!(rnap, dna1, -0.2)  
    intPot!(rnap, dna2, -0.2)  
    intPot!(rnap, promoter, -2.0)
    intPot!(rnap, enhancer, -0.2)

    intPot!(tf, rnap, -2.0)

    slist = [dna1, dna2, enhancer, dna1_l, dna2_l, promoter, tf, rnap]
    movingsps = [(tf, tf.intIdx), (rnap, rnap.intIdx)]  # species that move

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0]
    frame = 1
    deltaH = zeros(Float64, l.nspecies)

    hams = []
    accs = []
    transcription = []
    binding = []
    tf_density = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    rnap_density = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    tf_total = 0.0
    rnap_total = 0.0
    niter = 0

    for (sp, idx) in tf.intSpecies
        push!(hams, sp.ham)
        ham!(sp)
        sp.deltaH = 0.0
    end

    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1
        for (s, i) in movingsps
            if move!(s)
                nmoves[i] += 1
                
                # record data
                if mod(nmoves[i], step) == 0
                    deltaH[i] += s.deltaH

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += p.deltaHP
                        end
                    end
                    
                    ham!(s)
                    push!(hams, s.ham)
                    push!(accs, niter == 0 ? 0 : nmoves[i] / niter)
                    
                    deltaH[i] = 0
                    frame += 1
                    
                    # binding and transcription times
                    push!(binding, nspeciesBound(tf, [enhancer]))
                    push!(transcription, nspeciesBound(rnap, [enhancer]))

                    if t_current >= t_burn_in
                        for (_, val) in s.pos
                            if s.intIdx == tf.intIdx
                                tf_density[val[1], val[2], val[3]] += 1.0
                                tf_total += 1.0
                            elseif s.intIdx == rnap.intIdx
                                rnap_density[val[1], val[2], val[3]] += 1.0
                                rnap_total += 1.0
                            end
                        end
                    end
                end
            end
        end
    end

    # MARK: -Export

    # collect results
    tf_density = tf_density ./ tf_total
    rnap_density = rnap_density ./ rnap_total

    coords_tf = Matrix{Int32}(undef, tf.length, 4)

    res["species_name"] = tf.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["interaction_matrix"] = tf.intMat
    res["boxlength"] = tf.boxlength
    res["length"] = tf.length
    res["tf_density_lattice"] = tf_density
    res["binding"] = binding
    res["transcription"] = transcription
    res["rnap_density_lattice"] = rnap_density
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["tf_final_coordinates"] = extractCoords!(tf, coords_tf)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end