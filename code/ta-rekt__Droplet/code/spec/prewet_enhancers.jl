"""
Description:

Measuring the attractive force from a prewet droplet around two
nearby filaments
"""

### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 32,
    nspecies = 7
)

tf = FastBulk(;
    length = 4096,
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


### rod 2
dna1_l = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [1,10,16]
)

enhancer_l = FastBulk(;
    length = 6,
    initShape = "line",
    initParams = 1,
    initPos = [13,10,16]
)

dna2_l = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [19,10,16]
)

slist = [dna1, dna2, enhancer, dna1_l, dna2_l, enhancer_l, tf]

### simulation parameters
# MARK: -Sim. params

t_max = 10    # runtime in minutes
step = 1000
t_burn_in = 0.8 * t_max

# minimum size of a droplet for detection
minsize = 5

vals = range(1, stop=16, step=1)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [3]   # distance between filaments
betas = [0.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)
nspecies = 7
nmoving = 1

### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    dna1, dna2, enhancer, dna1_l, dna2_l, enhancer_l, tf = slist

    dna1.initPos[2] = dna1_l.initPos[2] + alpha
    dna2.initPos[2] = dna2_l.initPos[2] + alpha
    enhancer.initPos[2] = enhancer.initPos[2] + alpha

    # self-interactions
    intPot!(tf, -1.0)

    # multi-species interactions
    intPot!(tf, dna1_l, -0.5)  # a-c interaction
    intPot!(tf, dna2_l, -0.5)  # b-c interaction  
    intPot!(tf, enhancer_l, -2.0)  # a-b interaction (optional)

    intPot!(tf, dna1, -0.5)  # a-c interaction
    intPot!(tf, dna2, -0.5)  # b-c interaction  
    intPot!(tf, enhancer, -2.0)  # a-b interaction (optional)

    slist = [dna1, dna2, enhancer, dna1_l, dna2_l, enhancer_l, tf]         # all species
    movingsps = [(tf, tf.intIdx)]    # species that move

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = 0
    frame = 1
    deltaH = zeros(Float64, nspecies)

    hams = []
    accs = []
    density_lat = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    total = 0.0

    niter = 0

    for (sp, idx) in tf.intSpecies
        ham!(sp)
        push!(hams, sp.ham)
        sp.deltaH = 0.0
    end

    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1

                # record data
                if mod(nmoves, step) == 0
                    deltaH[i] += s.deltaH

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += p.deltaHP
                        end
                    end

                    ham!(s)
                    push!(hams, s.ham)
                    push!(accs, niter == 0 ? 0 : nmoves / niter)

                    deltaH[i] = 0
                    frame += 1

                    if t_current >= t_burn_in
                        for (_, val) in s.pos
                            density_lat[val[1], val[2], val[3]] += 1.0
                            total += 1.0
                        end
                    end
                end
            end
        end
    end

    # MARK: -Export

    # collect results
    density_lat = density_lat ./ total

    coords_c = Matrix{Int32}(undef, tf.length, 4)

    res["species_name"] = tf.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["boxlength"] = tf.boxlength
    res["length"] = tf.length
    res["density_lattice"] = density_lat
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(tf, coords_c)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end