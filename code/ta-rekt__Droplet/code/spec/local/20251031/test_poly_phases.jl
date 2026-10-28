### physical parameters

boxlength = 42
length = 32
nspecies = 1
maxlength = 512

global l = Lattice(;
    boxlength = boxlength,
    nspecies = nspecies
)

# polymer
global s1 = Species(;
    length = length,
    maxlength = maxlength,
    boxlength = boxlength,
    intPot = nothing,
    chemPot = nothing,
    initPos = [5,21,21]
)

# list containing all species
slist = [s1]


### simulation parameters

nframes = 20000
step = 100

vals = range(-1, stop=1, step=0.25)   # parallelized parameters (int. pot.)
alphas = zeros(63)

for (i, v) in enumerate(vals)
    alphas[(i-1)*7+1:i*7] = ones(7) .* v
end

#alphas = [1.0]
betas = [1.0]    # chemical potential
nparams = 63


### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    s1 = slist[1]
    s1.intPot = alpha
    s1.chemPot = beta

    # initialize lattice and species
    initPoly!(s1)
    initLattice!(l, [s1])

    # initialize hamiltonians
    initHamPoly!(s1)

    nmoves = 0
    niter = 0
    frame = 1

    radii = zeros(nframes)
    lengths = zeros(nframes)
    hams = zeros(nframes)
    accs = zeros(nframes)
    
    hamInitial = s1.ham   # initial energy
    
    # run simluation
    while frame <= nframes
        
        if mod(nmoves, step) == 0
            
            radii[frame] = radius(s1)
            lengths[frame] = s1.length
            hams[frame] = frame == 1 ? hamInitial : hams[frame-1] + s1.deltaH
            accs[frame] = niter == 0 ? 0 : round(nmoves / niter, digits=2)
            
            frame += 1
        end
        
        niter += 1
        nmoves += move!(s1)
    end

    # collect results
    coords = Matrix{Int32}(undef, s1.length, 3)

    res["species_name"] = s1.name
    res["nframes"] = nframes
    res["step"] = step
    
    res["alpha"] = alpha
    res["beta"] = beta
    res["maxlength"] = s1.maxlength
    res["boxlength"] = s1.boxlength

    res["gyration_radius"] = radii
    res["polymer_length"] = lengths
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(s1, coords)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end
