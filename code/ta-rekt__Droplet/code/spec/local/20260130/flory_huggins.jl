"""
Description: 

Trying to understand phase separation of a single species of bulk monomers
Building an empirical phase diagram, attempting to make it agree with 
Flory-Huggins theory. 
"""

### physical parameters

l = Lattice(;
    boxlength = 42,
    nspecies = 1
)

# bulk species
bulk = FastBulkMonomer(;
    length = 4096*3
)

# list containing all species
slist = [bulk]

# minimum size of a droplet for detection
minsize = 10

### simulation parameters

t_max = 60    # runtime in minutes
step = 10000

vals = range(0.1, stop=1.6, length=32)   # parallelized parameters (int. pot.)
alphas = zeros(32)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [0.6]   # temperatures
betas = [0.0]    # ignored
nparams = length(alphas)


### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, temp, beta = args

    bulk = slist[1]
    intPot!(bulk, - 1.0 / temp)

    # initialize lattice and species
    init!(l, [bulk])

    nmoves = 0
    iter = 0
    frame = 1

    ndroplets = []
    largest_droplet_size = []
    largest_droplet_density = []
    dilute_density = []
    hams = []
    accs = []

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
            end
        end
        
        moveMade = move!(bulk)
        nmoves += moveMade
        iter += 1
    end

    # collect results
    coords = Matrix{Int32}(undef, bulk.length, 4)

    res["species_name"] = bulk.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = temp
    res["beta"] = beta
    res["boxlength"] = bulk.boxlength
    res["length"] = bulk.length
    res["ndroplets"] = ndroplets
    res["largest_droplet_size"] = largest_droplet_size
    res["largest_droplet_density"] = largest_droplet_density
    res["dilute_density"] = dilute_density
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(bulk, coords)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
