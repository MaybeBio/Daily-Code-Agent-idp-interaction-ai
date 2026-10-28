"""
Description: 

Understanding the prewetting transition, situating it on the FH phase diagram
"""

### physical parameters

l = Lattice(;
    boxlength = 42,
    nspecies = 2
)

# bulk species
bulk = FastBulkMonomer(;
    length = 1
)

site = FastBulkMonomer(;
    length = 10,
    initShape = "line",
    initParams = 1,
    initPos = [16,21,21]
)

# list containing all species
slist = [bulk, site]

# minimum size of a droplet for detection
minsize = 10

### simulation parameters

t_max = 720# runtime in minutes
t_burn_in = 0.8 * t_max
step = 10000

# vals = range(0.4, stop=0.7, length=16)   
vals = range(0.005, stop=0.03, length=32)  # parallelized parameters (concentration)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [0.03]   # concentration
betas = [0.56]    # temperature
nparams = length(alphas)

### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    bulk = slist[1]
    intPot!(bulk, -1.0 / beta)
    intPot!(bulk, site, 0.0)
    intPot!(site, 0.0)

    bulk.length = Int32(round(l.boxlength^3 * alpha))

    # initialize lattice and species
    init!(l, [bulk, site])

    nmoves = 0
    iter = 0
    frame = 1

    ndroplets = []
    largest_droplet_size = []
    largest_droplet_density = []
    dilute_density = []
    hams = []
    accs = []
    density_lat = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    total = 0.0

    deltaH = 0
    moveMade = false

    droplets = filter(x -> length(x) >= minsize, cc(bulk))

    push!(ndroplets, length(droplets))
    push!(largest_droplet_size, 1)
    push!(largest_droplet_density, 0)
    push!(dilute_density, diluteDensity(bulk, droplets))
    push!(hams, bulk.ham)   # initial energy
    push!(accs, 0)
    
    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10

        if moveMade
            deltaH += bulk.deltaH
        
            if mod(nmoves, step) == 0
                frame += 1

                ccs = cc(bulk)
                droplets = filter(x -> length(x) >= minsize, ccs)

                if isempty(droplets)
                    maxsize = 1
                    dens = 0.0
                else
                    sizes = map(length, droplets)
                    maxsize, idx = findmax(sizes)   # returns (value, index)
                    dmax = droplets[idx]
                    dens = maxsize / dropletVolume(bulk, dmax)
                end
                
                push!(ndroplets, length(droplets))
                push!(largest_droplet_size, maxsize)
                push!(largest_droplet_density, dens)
                push!(dilute_density, diluteDensity(bulk, droplets))
                push!(hams, hams[end] + deltaH)
                push!(accs, iter == 0 ? 0 : nmoves / iter)
                
                deltaH = 0

                if t_current >= t_burn_in
                    for (_, val) in bulk.pos
                        density_lat[val[1], val[2], val[3]] += 1.0
                        total += 1.0
                    end
                end
            end
        end
        
        moveMade = move!(bulk)
        nmoves += moveMade
        iter += 1
    end

    # collect results
    coords = Matrix{Int32}(undef, bulk.length, 4)
    density_lat = density_lat ./ total

    res["species_name"] = bulk.name
    res["nframes"] = frame
    res["step"] = step
    res["concentration"] = alpha
    res["temperature"] = beta
    res["boxlength"] = bulk.boxlength
    res["length"] = bulk.length
    res["ndroplets"] = ndroplets
    res["largest_droplet_size"] = largest_droplet_size
    res["largest_droplet_density"] = largest_droplet_density
    res["dilute_density"] = dilute_density
    res["density_lattice"] = density_lat
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(bulk, coords)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
