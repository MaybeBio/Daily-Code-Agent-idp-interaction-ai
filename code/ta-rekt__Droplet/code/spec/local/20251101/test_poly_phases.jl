### physical parameters

global l = Lattice(;
    boxlength = 42,
    nspecies = 1
)

# polymer
global s1 = Species(;
    length = 32,
    maxlength = 512,
    boxlength = 42,
    intPot = nothing,
    chemPot = nothing,
    initPos = [5,21,21]
)

# list containing all species
slist = [s1]


### simulation parameters

# 5 minutes/million moves at mu=-1, 5 seconds at mu=1
nframes = 50000
step = 1000

vals = range(-1, stop=1, step=0.25)   # parallelized parameters (int. pot.)
alphas = zeros(63)
nrepeats = 7

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

#alphas = [-1.0]
betas = [1.0]    # chemical potential
nparams = length(alphas)


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
    deltaH = 0
    moveMade = false
    
    radii[frame] = radius(s1)
    lengths[frame] = s1.length
    hams[frame] = s1.ham   # initial energy
    accs[frame] = 0
    
    # run simluation
    while frame <= nframes - 1

        if moveMade
            deltaH += s1.deltaH
        
            if mod(nmoves, step) == 0
                
                frame += 1
                
                radii[frame] = radius(s1)
                lengths[frame] = s1.length
                hams[frame] = hams[frame-1] + deltaH
                accs[frame] = niter == 0 ? 0 : nmoves / niter
                
                deltaH = 0
            end
        end
        
        moveMade = movePoly!(s1)
        nmoves += moveMade
        niter += 1
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
