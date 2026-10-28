"""
Description:

What is the spatial density profile around a prewet filament as a function
of concentration in the bulk, and where is the prewetting transition on
the flory-huggins phase diagram
"""

### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 42,
    nspecies = 3
)

c = FastBulk(;
    length = 1,
    initShape = "uniform"
)

a1 = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [1,16,16]
)

b = FastBulk(;
    length = 6,
    initShape = "line",
    initParams = 1,
    initPos = [14,16,16]
)

a2 = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,
    initPos = [20,16,16]
)

slist = [a1, a2, b, c]

### simulation parameters
# MARK: -Sim. params

nframes = 10000
step = 10000
burn_in = 8000

# minimum size of a droplet for detection
minsize = 5

# vals = range(124, stop=1024, step=60)  # parallelized parameters (int. pot.)
vals = range(300, stop=8000, step=500)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [128]   # number of particles
betas = [0.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)
nspecies = 4
nmoving = 1

### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    a1, a2, b, c = slist

    c.length = alpha

    # self-interactions
    intPot!(c, -1.0)

    # multi-species interactions
    intPot!(c, a1, -0.5) 
    intPot!(c, a2, -0.5)   
    intPot!(c, b, -2.0) 

    slist = [a1, a2, b, c]         # all species
    movingsps = [(c, c.intIdx)]    # species that move

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = 0
    frame = 1
    deltaH = zeros(Float64, nspecies)

    ndroplets = []
    droplet_sizes = []
    droplet_densities = []
    dilute_density = []
    hams = []
    accs = []
    density_lat = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    total = 0.0

    niter = 0
    moveMade = falses(nmoving)

    for (sp, idx) in a1.intSpecies
        push!(hams, sp.ham)
        ham!(sp)
        sp.deltaH = 0.0
    end

    while frame <= nframes
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
                    
                    ccs = cc(s)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    
                    push!(ndroplets, length(droplets))
                    push!(droplet_sizes, map(x -> length(x), droplets))
                    push!(dilute_density, diluteDensity(s, droplets))
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

    coords_a1 = Matrix{Int32}(undef, a1.length, 4)
    coords_a2 = Matrix{Int32}(undef, a2.length, 4)
    coords_b = Matrix{Int32}(undef, b.length, 4)
    coords_c = Matrix{Int32}(undef, c.length, 4)

    res["species_name"] = c.name
    res["nframes"] = nframes
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["boxlength"] = c.boxlength
    res["length"] = c.length
    res["ndroplets"] = ndroplets
    res["droplet_sizes"] = droplet_sizes
    res["dilute_density"] = dilute_density
    res["density_lattice"] = density_lat
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(c, coords_c)
    res["final_coordinates_b"] = extractCoords!(a1, coords_a1)
    res["final_coordinates_c"] = extractCoords!(a2, coords_a2)
    res["final_coordinates_d"] = extractCoords!(b, coords_b)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end