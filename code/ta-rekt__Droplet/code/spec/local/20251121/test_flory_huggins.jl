### physical parameters

l = Lattice(;
    boxlength = 32,
    nspecies = 1
)

# bulk species
b = Bulk(;
    length = 8^3,
    intPot = nothing,
    chemPot = nothing
)

# list containing all species
slist = [b]

# minimum size of a droplet for detection
minsize = 5

### simulation parameters

nframes = 10000
step = 10000

vals = range(0.1, stop=2.0, step=0.1)   # parallelized parameters (int. pot.)
alphas = zeros(40)
nrepeats = 2

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [1.0]   # temperatures
betas = [1.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)


### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, temp, beta = args

    b = slist[1]
    b.intPot = - beta / temp
    b.chemPot = beta / temp

    # initialize lattice and species
    init!(l, [b])

    nmoves = 0
    niter = 0
    frame = 1

    ndroplets = zeros(nframes)
    droplet_sizes = Vector{Vector{Int32}}(undef, nframes)
    droplet_densities = Vector{Vector{Float64}}(undef, nframes)
    dilute_density = zeros(nframes)
    hams = zeros(nframes)
    accs = zeros(nframes)
    deltaH = 0
    moveMade = false

    droplets = filter(x -> length(x) >= minsize, cc(b))
    
    ndroplets[frame] = length(droplets)
    droplet_sizes[frame] = map(x -> length(x), droplets)
    droplet_densities[frame] = dropletDensity(b, droplets)
    dilute_density[frame] = diluteDensity(b)
    hams[frame] = b.ham   # initial energy
    accs[frame] = 0
    
    # run simluation
    while frame <= nframes - 1

        if moveMade
            deltaH += b.deltaH
        
            if mod(nmoves, step) == 0
                frame += 1

                ccs = cc(b)
                droplets = filter(x -> length(x) >= minsize, ccs)
                
                ndroplets[frame] = length(droplets)
                droplet_sizes[frame] = map(x -> length(x), droplets)
                droplet_densities[frame] = dropletDensity(b, droplets)
                dilute_density[frame] = diluteDensity(b)
                hams[frame] = hams[frame-1] + deltaH
                accs[frame] = niter == 0 ? 0 : nmoves / niter
                
                deltaH = 0
            end
        end
        
        moveMade = move!(b)
        nmoves += moveMade
        niter += 1
    end

    # collect results
    coords = Matrix{Int32}(undef, b.length, 3)

    res["species_name"] = b.name
    res["nframes"] = nframes
    res["step"] = step
    res["temperature"] = temp
    res["beta"] = beta
    res["boxlength"] = b.boxlength
    res["length"] = b.length
    res["ndroplets"] = ndroplets
    res["droplet_sizes"] = droplet_sizes
    res["droplet_densities"] = droplet_densities
    res["dilute_density"] = dilute_density
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(b, coords)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end
