using Droplet

# Lattice species (unused -- the Rouse solver is off-lattice).
slist = []
checkpoint_path = ""

# Fixed physical params
kBT = 1.0
b = 1.0            # Kuhn length
gamma = 1.0        # friction coefficient
kappa = 3*kBT/b^2  # bond stiffness
D = kBT / gamma    # diffusion coefficient
N = 400            # chain length
dt = 0.3 * gamma / kappa

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.jl.
# alphas: contact threshold (physical length).  beta is unused here.
alphas = [0.5]
betas  = [-1.4]

nparams = length(alphas)

blurb = "Solving the Langevin equation for a Rouse polymer to
         measure the distribution of first passage times to and
         from the entropic potential well. Computing the Fisher
         information rate from these distributions."

# Simulation lengths (deliberately small for the dummy run; scale up
# once the spec/dashboard pipeline is validated end-to-end).
nsteps_distance = 1_000_000
nsteps_full     = 200_000


function target!(res, args...)
    slist, keya, alpha, beta, blurb = args
    t0 = time()

    # primary run: distance-only mode, long trajectory for contact statistics
    r = RousePolymer(;
        boltzmannConst = kBT, 
        temperature    = 1.0,
        friction       = gamma, 
        kuhnLength     = b,
        nbeads         = Int32(N), 
        nsteps         = Int32(nsteps_distance),
        dt             = dt,
        trackedBead1   = Int32(1), 
        trackedBead2   = Int32(N),
        saveMode       = :distance_only, 
        seed = 42,
    )
    init!(r)
    simulate!(r)

    times, invCdf = contactTimes(r.distanceTraj, alpha, dt)
    centers, counts = distanceHistogram(r.distanceTraj, 80)

    # secondary run: short :full trajectory just for MSD curves
    rFull = RousePolymer(;
        boltzmannConst = kBT, 
        temperature    = 1.0,
        friction       = gamma, 
        kuhnLength     = b,
        nbeads         = Int32(N), 
        nsteps         = Int32(nsteps_full),
        dt             = dt,
        trackedBead1   = Int32(1), 
        trackedBead2   = Int32(N),
        saveMode       = :full, 
        seed = 43,
    )
    init!(rFull)
    simulate!(rFull)

    maxLag = div(nsteps_full, 4)
    msdBead1   = beadMSD(rFull.posTraj, 1; maxLag = maxLag)
    msdBeadMid = beadMSD(rFull.posTraj, div(N, 2); maxLag = maxLag)
    msdPair    = beadPairDistanceMSD(rFull.posTraj, 1, N; maxLag = maxLag)

    res["blurb"]            = blurb
    res["alpha"]            = alpha
    res["beta"]             = beta
    res["params"]           = Dict("N" => N, "b" => b, "kBT" => kBT,
                                   "gamma" => gamma, "kappa" => kappa,
                                   "D" => D, "dt" => dt,
                                   "nsteps_distance" => nsteps_distance,
                                   "nsteps_full"     => nsteps_full)
    res["distance_traj"]    = r.distanceTraj
    res["contact_times"]    = times
    res["inverse_cdf"]      = invCdf
    res["hist_centers"]     = centers
    res["hist_counts"]      = counts
    res["bead_relax_time"]  = r.timeScale
    res["bead_relax_scale"] = r.posScale
    res["msd_bead1"]        = msdBead1
    res["msd_beadMid"]      = msdBeadMid
    res["msd_pairDistance"] = msdPair
    res["runtime_min"]      = (time() - t0) / 60.0

    return res
end
