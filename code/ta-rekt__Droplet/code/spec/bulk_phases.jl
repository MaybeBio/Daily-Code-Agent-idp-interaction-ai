"""
Description: 

Trying to understand phase separation of a single species of bulk monomers
Plotting mean radius of gyration of the species versus self-interaction
"""

### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 42,
    nspecies = 2
)

# bulk species
bulk = FastBulkMonomer(;
    length = 4096,
    initShape = "uniform"
)

# list containing all species
slist = [bulk]


### simulation parameters
# MARK: -Sim. params

t_max = 10    # runtime in minutes
step = 10000
t_burn_in = 0.8 * t_max

vals = range(-2, stop=0, step=0.25)   # parallelized parameters (bulk self-interaction)
alphas = zeros(9)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [-1.0]
betas = [0.0]    # placeholder
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    bulk = slist[1]
    
    # set interaction potentials
    intPot!(bulk, alpha)       # bulk self-interaction (fixed)

    # initialize lattice and species
    init!(l, [bulk])

    # moving species
    movingsps = [(bulk, bulk.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1

    radii = []
    hams_bulk = []
    hams_bulk = []
    accs = []
    deltaH = zeros(l.nspecies)
    
    # initial values
    push!(radii, radius(bulk))
    push!(hams_bulk, bulk.ham)
    push!(accs, 0)
    
    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1

                # record data
                if mod(nmoves, step) == 0
                    deltaH[i] += getDeltaHam(s)

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += getDeltaHP(p)
                        end
                    end

                    ham!(s)

                    if s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    push!(radii, radius(bulk))
                    push!(accs, niter == 0 ? 0 : nmoves / niter)

                    deltaH[i] = 0
                    frame += 1
                end
            end
        end
    end

    # collect results
    coords_bulk = Matrix{Int32}(undef, bulk.length, 4)

    res["species_name"] = bulk.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = 1.0
    res["alpha"] = alpha
    res["beta"] = beta
    res["interaction_matrix"] = bulk.intMat
    res["boxlength"] = bulk.boxlength
    res["length"] = bulk.length
    res["radius_of_gyration"] = radii
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["final_coordinates_bulk"] = extractCoords!(bulk, coords_bulk)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
