### physical parameters
# MARK: -Phys. params

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 16,
    nspecies = 1
)

b = BulkMonomer(;
    length = 2,
    initShape = "uniform",
    k0 = 1.0
)

slist = [b]


### simulation parameters
# MARK: -Sim. params

blurb = """
KMC free-diffusion benchmark for a single BulkMonomer species. Uses the
Metropolis-rate Gillespie sampler (src/gillespie.jl), with a physically
calibrated simulation clock t_sim. Samples positions at fixed dt_sample of
sim time so dashboards can plot single-locus and pair-locus MSD vs t_sim.

swept parameters: alpha = self-interaction potential (kT units)
output metrics: position trajectories, hamiltonian, simulation time grid
"""

init_from_checkpoint = false

t_max = 0.02      # wall-clock cap (minutes), initial small test
dt_sample = 0.5   # simulation-time sampling interval
step = 1          # CLAUDE.md compatibility (real cadence is dt_sample)

alphas = [0.0]    # self-interaction potential (sweep)
betas  = [1.0]    # inverse temperature
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()

    slist, keya, alpha, beta, blurb = args
    b = slist[1]

    init!(l, slist)
    intPot!(b, alpha)

    sys = KMCSystem(; slist=slist, beta=beta)
    buildMoves!(sys)

    t_samples = Float64[]
    positions = Vector{Vector{Vector{Int32}}}()
    unwrapped = Vector{Vector{Vector{Int}}}()   # cumulative PBC-unwrapped Δ per particle
    hams = Float64[]

    last_snap = nothing
    cum = [zeros(Int, 3) for _ in 1:b.length]
    L = Int(l.boxlength); halfL = L ÷ 2

    sample_cb = function (s)
        push!(t_samples, s.t_sim)
        snap = [copy(b.pos[Int32(i)][1:3]) for i in 1:b.length]

        # accumulate unwrapped displacement with single-wrap correction —
        # valid as long as inter-sample hop count keeps |Δ| < L/2, which holds
        # for k0·dt_sample ≪ (L/2)² (≈ 64 for L=16).
        if last_snap !== nothing
            for i in 1:b.length, ax in 1:3
                d = Int(snap[i][ax]) - Int(last_snap[i][ax])
                d > halfL  && (d -= L)
                d < -halfL && (d += L)
                cum[i][ax] += d
            end
        end

        push!(positions, snap)
        push!(unwrapped, [copy(cum[i]) for i in 1:b.length])
        ham!(b)
        push!(hams, b.ham)
        last_snap = snap
    end

    kmc_run!(sys; t_max_wall=t_max, dt_sample=dt_sample, sample_callback=sample_cb)

    final_coords = zeros(Int32, b.length, 4)
    extractCoords!(b, final_coords)

    res["blurb"] = blurb
    res["species_name"] = b.name
    res["nframes"] = length(t_samples)
    res["step"] = step
    res["dt_sample"] = dt_sample
    res["boxlength"] = l.boxlength
    res["length"] = b.length
    res["k0"] = b.k0
    res["beta"] = beta
    res["intPot"] = alpha
    res["interaction_matrix"] = b.intMat
    res["t_sim_samples"] = t_samples
    res["t_sim_final"] = sys.t_sim
    res["n_events"] = sys.n_events
    res["positions"] = positions
    res["unwrapped_displacements"] = unwrapped
    res["hamiltonian"] = hams
    res["final_coordinates"] = final_coords
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
