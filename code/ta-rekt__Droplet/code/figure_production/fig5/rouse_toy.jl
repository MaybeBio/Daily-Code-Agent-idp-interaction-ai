# FIG. 5 PANEL (d): the June toy model, re-run.
#
# Copy of spec/rouse_mfpt_attractive.jl at commit 9f78af2 (Tarek, 2026-06-15),
# whose output data/20260615/rouse_mfpt_attractive fed viz-rouse.ipynb cell 37
# and the old Fig. 4 panel D. That data folder no longer exists; this
# reproduces it statistically (not bitwise: different Julia/RNG versions).
#
# Changed from 9f78af2, nothing else:
#   - alphas: only the two concentrations the figure plots, June indices 1 and 10
#     (1-based) of range(0.003, 0.15, length=16), four independent runs each,
#     pooled by prepare_fig5_data.py;
#   - seeds: the June ones (keya + 42 -> 43 and 52) plus 1000 per extra run;
#   - stored: the pair potential alone (ads_potential) and its grid (pot_grid).
#
# Run with run_fig5_rouse.sh from the repo root.

using Droplet

# Lattice species (unused -- the Rouse solver is off-lattice).
slist = []
checkpoint_path = ""

# Rouse polymer parameters
kBT = 1.0
kuhnLength = 1.0            # Kuhn length
friction = 1.0        # friction coefficient
stiffness = 2*kBT / kuhnLength^2  # bond stiffness
D = kBT / friction    # diffusion coefficient
N = 100            # chain length
firstBeadIdx = 1
beadSep = N - 1
dt = 0.001 * friction / stiffness   # 0.3 times the characteristic time

# Surface adsorption model parameters
nb              = 1.0
chi_b           = 2.5
omega_b         = -0.5
kappa           = 1.0
phi_inf         = 0.1
eps             = 1e-6
R               = 2.5

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.jl.
# alphas: phi_inf.  beta is unused here.
alphas_june = collect(range(0.003, stop=0.15, length=16))
alphas = alphas_june[[1, 1, 1, 1, 10, 10, 10, 10]]
seeds  = [43, 1043, 2043, 3043, 52, 1052, 2052, 3052]
betas  = [-1.4]

nparams = length(alphas)

blurb = "Fig. 5 panel (d): June toy model (spec 9f78af2), distance histograms
         of a Rouse chain whose end beads feel the linear-response surface
         attraction, at low and high concentration."

nsteps = 50_000_000

function target!(res, args...)
    slist, keya, alpha, beta, blurb = args
    t0 = time()

    # Surface adsorption model parameters
    p = Parameters(
        nb              = nb,
        chi_b           = chi_b,
        omega_b         = omega_b,
        kappa           = kappa,
        phi_inf         = alpha,
        eps             = eps,
        R               = R,
    )

    # Parameters of the exponential-integral pair potential.
    d2f_phi_inf = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(d2f_phi_inf / p.kappa)
    dists = collect(range(p.R, p.L, length=500))

    # obtained from linear-response concentration profile
    function adsorption_potential(x)
        return -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
                (p.kappa * (lam + 1/p.R) * x)
    end

    # f(d) = -dU/dr, computed numerically via central differences from
    # adsorption_potential — single source of truth for the pair interaction.
    h_fd = 1e-6
    attractive_force = x ->
        -(adsorption_potential(x + h_fd) - adsorption_potential(x - h_fd)) / (2 * h_fd)

    function total_potential(x)
        return (adsorption_potential(x)
                - 2*log(x) + 3*x^2 / (2*beadSep * kuhnLength^2))
    end


    # primary run: distance-only mode, long trajectory for contact statistics
    r = RousePolymer(;
        boltzmannConst  = kBT,
        temperature     = 1.0,
        friction        = friction,
        kuhnLength      = kuhnLength,
        r_min           = p.R,
        nbeads          = Int32(N),
        nsteps          = Int32(nsteps),
        dt              = dt,
        trackedBead1    = firstBeadIdx,
        trackedBead2    = firstBeadIdx + beadSep,
        attractiveForce = attractive_force,
        saveMode        = :distance_only,
        seed            = seeds[keya],
    )

    init!(r)
    simulate!(r)

    res["blurb"]            = blurb
    res["alpha"]            = alpha
    res["beta"]             = beta
    res["seed"]             = seeds[keya]
    res["surface_params"]   = Dict("nb" => p.nb, "chi_b" => p.chi_b,
                            "kappa" => p.kappa, "omega_b" => p.omega_b,
                            "R" => p.R, "L" => p.L, "phi_inf" => p.phi_inf,
                            "eps" => p.eps)
    res["rouse_params"]     = Dict("N" => N, "kuhn_length" => kuhnLength, "kBT" => kBT,
                            "friction" => friction, "stiffness" => stiffness,
                            "D" => D, "dt" => dt, "bead_separation" => beadSep,
                            "nsteps" => nsteps)
    sample_skip = 500
    res["distance_traj"]    = r.distanceTraj[1:sample_skip:end]
    res["sample_skip"]      = sample_skip
    res["pot_grid"]         = dists
    res["total_potential"]  = total_potential.(dists)
    res["ads_potential"]    = adsorption_potential.(dists)
    res["bead_relax_time"]  = r.timeScale
    res["bead_relax_scale"] = r.posScale
    res["runtime_min"]      = (time() - t0) / 60.0

    return res
end
