### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 16,
    nspecies = 2
)

# polymer species
poly = Polymer(;
    length = 16,
    maxlength = 128,
    initShape = "line",
    nMonomerTypes = 1,
    initPos = [1, 8, 8],
    pattern=:uniform
)

initMonomers!(poly)
mono = poly.monomers[1]

# bulk species
bulk = FastBulkMonomer(;
    length = 1024,
    initShape = "uniform"
)

# list containing all species
slist = [mono, bulk]


### simulation parameters
# MARK: -Sim. params

t_max = 60    # runtime in minutes
step = 10000
t_burn_in = 0.8 * t_max

vals = range(-1.125, stop=1.0, step=0.125)   # parallelized parameters (polymer self-interaction)
alphas = zeros(18)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [-1.0]
betas = [0.0]    # placeholder
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    mono, bulk = slist
    
    # set interaction potentials
    intPot!(mono, alpha)      # polymer self-interaction (SWEEP)
    intPot!(bulk, -1.0)       # bulk self-interaction (fixed)
    intPot!(bulk, mono, -1.0) # bulk-polymer interaction (fixed)

    # initialize lattice and species
    init!(l, [mono, bulk])

    # moving species
    movingsps = [(mono, mono.intIdx), (bulk, bulk.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1

    radii = []
    hams_poly = []
    hams_bulk = []
    accs = []
    deltaH = zeros(l.nspecies)
    
    # initial values
    push!(radii, radius(poly))
    push!(hams_poly, poly.ham)
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

                    if s.intIdx == mono.intIdx
                        push!(hams_poly, getHam(s))

                    elseif s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    push!(radii, radius(poly))
                    push!(accs, niter == 0 ? 0 : nmoves / niter)

                    deltaH[i] = 0
                    frame += 1
                end
            end
        end
    end

    # collect results
    coords_poly = Matrix{Int32}(undef, poly.length, 4)
    coords_bulk = Matrix{Int32}(undef, bulk.length, 4)

    res["species_name"] = poly.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = 1.0
    res["alpha"] = alpha
    res["interaction_matrix"] = mono.intMat
    res["boxlength"] = poly.boxlength
    res["length"] = poly.length
    res["radius_of_gyration"] = radii
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["final_coordinates_poly"] = extractCoords!(poly, coords_poly)
    res["final_coordinates_bulk"] = extractCoords!(bulk, coords_bulk)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
