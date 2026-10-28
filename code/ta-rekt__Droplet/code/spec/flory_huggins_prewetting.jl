### physical parameters

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 40,
    nspecies = 2
)

# bulk species
bulk = FastBulkMonomer(;
    length = 1
)

site = FastBulkMonomer(;
    length = 64,
    initShape = "sphere",
    initParams = 3,
    initPos = [10,10,10]
)

# list containing all species
slist = [bulk, site]

# minimum size of a droplet for detection
minsize = 10

### simulation parameters

blurb = """
Fitting the wetting and prewetting transitions with a spherical surface.
Turn off interaction with the surface to check droplet persistence.
"""

init_from_checkpoint = false

t_max = 2*60    # runtime in minutes
t_burn_in = t_max / 10
step = l.boxlength^3 * 0.5   # vol fraction 0.5
bigstep = 100 * step

vals = range(1.5, stop=3.5, length=16)  # parallelized parameters (concentration)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [1.5]   # self-interaction
betas = [0.5]    # vol.fraction
mu = 0.1       # chemical potential
nparams = length(alphas)

### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    phi_bar = beta   # volume fraction
    chi = alpha      # bulk self-interaction
    h0 = 0.0         # surface-bulk interaction

    bulk, site = slist

    bulk.length = Int32(round(l.boxlength^3 * phi_bar))
    bulk.chemPot = -mu -chi         # mu: chemical potential

    intPot!(bulk, -chi/3)
    intPot!(bulk, site, h0)   # h0: surface-bulk interaction

    # initialize lattice and species
    if init_from_checkpoint
        checkpoint_at_param = string(checkpoint_path, "_keya", keya)
        restore!(checkpoint_at_param, l, slist)
    else
        init!(l, slist)
    end

    # moving species
    movingsps = [(bulk, bulk.intIdx)]

    nmoves = 0
    iter = 0
    frame = 1
    hist_count = 0

    hams = []
    accs = []
    normal_hist = zeros(Float64, 100)
    max_distance = 10.0

    occupation_counts = zeros(Int32, 100)
    density_lat = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    total = 0.0

    deltaH = 0
    moveMade = false

    push!(hams, bulk.ham)   # initial energy
    push!(accs, 0)
    
    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        iter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1
        
                if mod(nmoves, step) == 0
                    frame += 1
                    
                    push!(hams, getHam(s))
                    push!(accs, iter == 0 ? 0 : nmoves / iter)
                    
                    s.deltaH = 0

                    if t_current >= t_burn_in
                        for (_, val) in bulk.pos
                            density_lat[val[1], val[2], val[3]] += 1.0
                            total += 1.0
                        end

                        hist_count = updateNormalHistogram!(bulk, normal_hist, [1.0,1.0,1.0], max_distance)
                    end
                end

                if mod(nmoves, bigstep) == 0 && t_current >= t_burn_in

                    # compute occupancy histograms
                    for x in range(1, l.boxlength)
                        for y in range(1, l.boxlength)
                            for z in range(1, l.boxlength)
                                density = density_lat[x,y,z] / total
                                idx = Int(round(density*100*bulk.length, digits=0))
                                idx = max(1, idx)
                                occupation_counts[idx] += 1
                            end
                        end
                    end

                    fill!(density_lat, 0.0)
                    total = 0.0
                end
            end
        end
    end

    # collect results
    coords = Matrix{Int32}(undef, bulk.length, 4)

    res["blurb"] = blurb
    res["species_name"] = bulk.name
    res["nframes"] = frame
    res["step"] = step
    res["concentration"] = alpha
    res["temperature"] = beta
    res["boxlength"] = bulk.boxlength
    res["length"] = bulk.length
    res["occupation_histogram"] = occupation_counts
    res["normal_density"] = normal_hist
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(bulk, coords)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
