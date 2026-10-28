### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 20,
    nspecies = 7
)

c = FastBulk(;
    length = 500,
    initShape = "uniform"
)

### rod 1
a1_l = FastBulk(;
    length = 7,
    initShape = "line",
    initParams = 1,
    initPos = [1,5,8]
)

b_l = FastBulk(;
    length = 6,
    initShape = "line",
    initParams = 1,
    initPos = [8,5,8]
)

a2_l = FastBulk(;
    length = 7,
    initShape = "line",
    initParams = 1,
    initPos = [14,5,8]
)

### rod 2
a1_r = FastBulk(;
    length = 7,
    initShape = "line",
    initParams = 1,
    initPos = [1,12,8]
)

b_r = FastBulk(;
    length = 6,
    initShape = "line",
    initParams = 1,
    initPos = [8,12,8]
)

a2_r = FastBulk(;
    length = 7,
    initShape = "line",
    initParams = 1,
    initPos = [14,12,8]
)

slist = [a1_l, a1_r, a2_l, a2_r, b_l, b_r, c]

### simulation parameters
# MARK: -Sim. params

nframes = 50000
step = 10000
burn_in = 40000

# minimum size of a droplet for detection
minsize = 5

vals = range(1, stop=8, step=1)
alphas = zeros(8)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [4]   # distance between filaments
betas = [0.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)
nspecies = 7
nmoving = 1

### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    a1_l, a1_r, a2_l, a2_r, b_l, b_r, c = slist

    a1_r.initPos[2] = a1_l.initPos[2] + alpha
    a2_r.initPos[2] = a2_l.initPos[2] + alpha
    b_r.initPos[2] = b_l.initPos[2] + alpha

    # self-interactions
    intPot!(c, -1.0)

    # multi-species interactions
    intPot!(c, a1_l, -0.5)  # a-c interaction
    intPot!(c, a2_l, -0.5)  # b-c interaction  
    intPot!(c, b_l, -2.0)  # a-b interaction (optional)

    intPot!(c, a1_r, -0.5)  # a-c interaction
    intPot!(c, a2_r, -0.5)  # b-c interaction  
    intPot!(c, b_r, -2.0)  # a-b interaction (optional)

    slist = [a1_l, a1_r, a2_l, a2_r, b_l, b_r, c]         # all species
    movingsps = [(c, c.intIdx)]    # species that move

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

    for (sp, idx) in a1_l.intSpecies
        push!(hams, sp.ham)
        ham!(sp)
        sp.deltaH = 0.0
    end

    while frame <= nframes
        niter += 1
        for (s, i) in movingsps
            if  move!(s)
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

                    push!(accs, niter == 0 ? 0 : nmoves / niter)
                    push!(hams, hams[end] + deltaH[i])

                    deltaH[i] = 0
                    frame += 1

                    if frame >= burn_in
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

    coords_c = Matrix{Int32}(undef, c.length, 4)

    res["species_name"] = c.name
    res["nframes"] = nframes
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["boxlength"] = c.boxlength
    res["length"] = c.length
    res["density_lattice"] = density_lat
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(c, coords_c)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end
