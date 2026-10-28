# Fisher information rate from the full 2D two-sphere BVP, with BOTH sides of
# the capillary-bridging bifurcation.
#
# Parallelised axis: alphas = phi_inf grid.  For each phi_inf the target solves
# the axisymmetric two-sphere BVP at a grid of separations d on TWO branches:
#   * unbridged (dilute gap)  -- marched inward from large d;
#   * bridged   (dense gap)   -- marched outward from contact, seeded by a dense
#                                blob (capillary condensation).
# The equilibrium pair potential is U_full(d) = min(F_branch)(d) - F(d_ref), and
# the surface volume fraction Phi_surf(d) is read off the equilibrium branch.
# These are integrated against the polymer P_eq(d) to obtain phi_probe; the
# Fisher rate is assembled across phi_inf in the dashboard.
#
# Bridging only exists when the bulk is near saturation, so phi_inf is swept just
# BELOW the dilute binodal (set chi_b = 2.5; binodal ~ 0.145).
#
# ALL parameters (surface adsorption model, BVP solver, Rouse chain, FPT solver,
# and the swept axes) come from spec/jobs/20260720/param.jl via `PARAMS`; edit
# that file to reconfigure.  Driver: sim/sweep.jl.  Pipe builder: runJL.

using Droplet

include(joinpath(@__DIR__, "param.jl"))

# Lattice species (unused -- the BVP doesn't simulate particles).
slist = []
checkpoint_path = ""

# Dilute binodal for the configured (chi_b, nb): bridging needs phi_inf just below
# it, and the dense coexistence value 1 - BINODAL seeds the bridged branch.
const BINODAL = dilute_binodal(PARAMS.surface.chi_b; nb = PARAMS.surface.nb)

# Sweep axes (sharded across CPUs by sim/sweep.jl): phi_inf signal x omega_b field.
alphas  = phi_inf_grid(PARAMS)          # phi_inf
betas   = PARAMS.sweep.omega_b_values   # omega_b
nparams = length(alphas)

blurb = "Fisher information rate from the full 2D two-sphere BVP, capturing " *
        "both sides of the capillary-bridging bifurcation.  Sweeps phi_inf " *
        "(alphas, as a fraction of the dilute binodal) at fixed (omega_b, R, " *
        "chi_b=2.5); for each phi_inf solves the unbridged (inward) and " *
        "bridged (outward, dense-seeded) branches on a d-grid, takes the " *
        "equilibrium U_full(d)=min(F)-F_ref and Phi_surf(d) on that branch, " *
        "and integrates the polymer P_eq(d) to get phi_probe."

# Fixed, anchored grid (R and dz are constant across the sweep).  (1) build_geometry
# spans z over length 4R + d + 2*z_buf; (2) pin dz EXACTLY to bvp.dz_target and place
# z_buf on the grid so the symmetry plane (sphere 1 at z=0) is a grid line; snapping d
# to a multiple of dz puts BOTH sphere surfaces on grid lines for every d, so the film
# free energy discretises identically and cancels in U = F(d) - F_ref.
const _R  = PARAMS.surface.R
const _DZ = PARAMS.bvp.dz_target
const Z_BUF_GRID = ceil(Int, (_R + PARAMS.bvp.z_buf) / _DZ) * _DZ - _R
_snap(d::Real)   = round(d / _DZ) * _DZ
_nz_for(d::Real) = Int(round((4.0 * _R + _snap(d) + 2.0 * Z_BUF_GRID) / _DZ)) + 1

# Warm-start guess for `geom` from a previous solution (phi_prev on geom_prev):
# sample it at geom's unknown cells, folding sphere-2's half onto the fixed
# sphere-1 half (sphere 1 stays at z=0 for all d; sphere 2 is its mirror image).
function _warm_x0(geom, phi_prev, geom_prev)
    z_mid = 0.5 * (geom.z_c1 + geom.z_c2)
    rho_q = Float64[]
    z_q   = Float64[]
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        geom.unk_mask[i, j] || continue
        push!(rho_q, geom.rho_grid[i])
        zz = geom.z_grid[j]
        push!(z_q, zz > z_mid ? 2.0 * z_mid - zz : zz)
    end
    return clamp.(interp_phi(phi_prev, geom_prev, rho_q, z_q),
                  geom.p.eps, 1.0 - geom.p.eps)
end

# Bridged seed: dense phase within `reach` of either sphere surface (fills the
# neck), bulk phi_inf elsewhere.  Bootstraps the capillary-bridge branch.
function _dense_seed(geom, phi_dense; reach = 2.0)
    R = geom.p.R
    x = Float64[]
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        geom.unk_mask[i, j] || continue
        rho = geom.rho_grid[i]; z = geom.z_grid[j]
        s = min(sqrt(rho^2 + (z - geom.z_c1)^2),
                sqrt(rho^2 + (z - geom.z_c2)^2)) - R
        push!(x, s < reach ? phi_dense : geom.p.phi_inf)
    end
    return clamp.(x, geom.p.eps, 1.0 - geom.p.eps)
end

# Gap-centre phi (on axis, at the midplane) -- the bridged/unbridged label.
_gap_center(geom, phi) =
    phi[1, argmin(abs.(geom.z_grid .- 0.5 * (geom.z_c1 + geom.z_c2)))]

function target!(res, args...)
    slist, keya, alpha, beta, blurb = args

    t0 = time()

    bp      = PARAMS.bvp
    omega_b = beta
    phi_inf = alpha

    # Surface adsorption model: base Parameters from PARAMS, swept fields overridden.
    p = Droplet.with(PARAMS.surface; phi_inf = phi_inf, omega_b = omega_b)
    R = p.R

    N_D    = bp.n_d
    d_min  = bp.d_min_factor * R
    d_max  = bp.d_max_factor * R
    d_grid = _snap.(collect(range(d_min, d_max; length = N_D)))
    d_ref  = _snap(bp.d_ref_factor * R)
    phi_dense = 1.0 - BINODAL

    geom_at(dd) = build_geometry(p, dd; rho_max = bp.rho_max, z_buf = Z_BUF_GRID,
                                 Nrho = bp.Nrho, Nz = _nz_for(dd))
    _solve(dd; x0 = nothing, omega_cont = true) =
        solve_two_spheres(p, dd; rho_max = bp.rho_max, z_buf = Z_BUF_GRID,
                          Nrho = bp.Nrho, Nz = _nz_for(dd),
                          f_tol = bp.f_tol, maxiter = bp.maxiter,
                          x0_init = x0, omega_continuation = omega_cont,
                          verbose = false, bc_reg = bp.bc_reg)

    # Reference: spheres far apart on the dilute branch.
    res_ref = _solve(d_ref)
    F_ref = res_ref[:success] ? excess_free_energy_two_spheres(res_ref) : NaN

    res_unbr = Vector{Any}(nothing, N_D)   # dilute branch result per d
    res_brid = Vector{Any}(nothing, N_D)   # bridged branch result per d
    F_unbr   = fill(NaN, N_D)
    F_brid   = fill(NaN, N_D)

    # --- Unbridged (dilute) branch: inward march from large d ---
    prev_phi  = res_ref[:success] ? res_ref[:phi]  : nothing
    prev_geom = res_ref[:success] ? res_ref[:geom] : nothing
    for k in sortperm(d_grid; rev = true)
        dd = d_grid[k]
        x0 = prev_phi === nothing ? nothing :
             _warm_x0(geom_at(dd), prev_phi, prev_geom)
        rk = _solve(dd; x0 = x0)
        (rk[:success] && _gap_center(rk[:geom], rk[:phi]) < bp.bridge_phi) || continue
        res_unbr[k] = rk
        F_unbr[k]   = excess_free_energy_two_spheres(rk)
        prev_phi    = rk[:phi]
        prev_geom   = rk[:geom]
    end

    # --- Bridged (dense gap) branch: outward march from contact ---
    prev_phi = nothing; prev_geom = nothing
    for k in sortperm(d_grid)
        dd = d_grid[k]
        geom_k = geom_at(dd)
        guesses = prev_phi === nothing ?
            (_dense_seed(geom_k, phi_dense),) :
            (_warm_x0(geom_k, prev_phi, prev_geom), _dense_seed(geom_k, phi_dense))
        for x0 in guesses
            rk = _solve(dd; x0 = x0, omega_cont = false)
            if rk[:success] && _gap_center(rk[:geom], rk[:phi]) > bp.bridge_phi
                res_brid[k] = rk
                F_brid[k]   = excess_free_energy_two_spheres(rk)
                prev_phi    = rk[:phi]
                prev_geom   = rk[:geom]
                break
            end
        end
    end

    # --- Backfill: inward march over any d the outward march left unsolved.
    # Cold dense seeds are fragile at isolated d; warm-starting from the
    # nearest converged bridged solution at larger d usually converges. ---
    prev_phi = nothing; prev_geom = nothing
    for k in sortperm(d_grid; rev = true)
        if res_brid[k] !== nothing
            prev_phi  = res_brid[k][:phi]
            prev_geom = res_brid[k][:geom]
            continue
        end
        prev_phi === nothing && continue
        dd = d_grid[k]
        rk = _solve(dd; x0 = _warm_x0(geom_at(dd), prev_phi, prev_geom),
                    omega_cont = false)
        if rk[:success] && _gap_center(rk[:geom], rk[:phi]) > bp.bridge_phi
            res_brid[k] = rk
            F_brid[k]   = excess_free_energy_two_spheres(rk)
            prev_phi    = rk[:phi]
            prev_geom   = rk[:geom]
        end
    end

    # --- Equilibrium = lower-free-energy branch at each d ---
    U_unbr   = fill(NaN, N_D)
    U_brid   = fill(NaN, N_D)
    U_full   = fill(NaN, N_D)      # equilibrium pair potential
    F_arr    = fill(NaN, N_D)
    Phi_surf = fill(NaN, N_D)
    branch   = zeros(Int, N_D)     # 0 none, 1 unbridged, 2 bridged
    # Gap-facing pole of sphere 1: theta = 0 => z = z_c1 + R, the point pointing
    # inward toward sphere 2 (the enriched adsorption film).  The far-side pole
    # (theta = pi, z = z_c1 - R) sits in open bulk and reads ~phi_inf -- that was
    # the previous bug.  Matches inner_pole_phi in visualize_two_spheres_2d.py.
    gap_pole = [0.0]

    # Per-d full 2D phi(rho,z) of the equilibrium branch.  Stored as
    # Vector{Vector{Float64}}, outer index rho, inner z, so on the Python side
    # `np.array(profile_phi[k])` is shape (Nrho, Nz).
    profile_d   = Float64[]
    profile_rho = Vector{Vector{Float64}}()
    profile_z   = Vector{Vector{Float64}}()
    profile_phi = Vector{Vector{Vector{Float64}}}()

    for k in 1:N_D
        fu = F_unbr[k]; fb = F_brid[k]
        if isfinite(F_ref)
            isnan(fu) || (U_unbr[k] = fu - F_ref)
            isnan(fb) || (U_brid[k] = fb - F_ref)
        end
        chosen = nothing; F_eq = NaN
        if !isnan(fu) && (isnan(fb) || fu <= fb)
            chosen = res_unbr[k]; F_eq = fu; branch[k] = 1
        elseif !isnan(fb)
            chosen = res_brid[k]; F_eq = fb; branch[k] = 2
        end
        chosen === nothing && continue
        F_arr[k]    = F_eq
        U_full[k]   = isfinite(F_ref) ? F_eq - F_ref : NaN
        # Clamp to the physical range: phi_on_sphere's quadratic extrapolation
        # overshoots (phi_s ~ 5) where the probe stencil degenerates at contact,
        # and >1 / <0 surface concentrations are nonsense downstream.
        Phi_surf[k] = clamp(phi_on_sphere(chosen[:phi], chosen[:geom],
                                          chosen[:geom].z_c1, gap_pole)[1], 0.0, 1.0)
        push!(profile_d,   d_grid[k])
        push!(profile_rho, copy(chosen[:geom].rho_grid))
        push!(profile_z,   copy(chosen[:geom].z_grid))
        push!(profile_phi, [chosen[:phi][i, :] for i in 1:size(chosen[:phi], 1)])
    end

    ok = isfinite.(U_full) .& isfinite.(Phi_surf)
    n_ok = count(ok)
    success = isfinite(F_ref) && all(ok)   # reference + every d converged

    res["blurb"]         = blurb
    res["phi_inf"]       = phi_inf
    res["omega_b"]       = omega_b
    res["R"]             = R
    res["chi_b"]         = p.chi_b
    res["binodal"]       = BINODAL
    res["s"]             = PARAMS.rouse.s
    res["a"]             = PARAMS.rouse.a
    res["d_grid"]        = d_grid
    res["d_ref"]         = d_ref
    res["F_ref"]         = F_ref
    res["F"]             = F_arr
    res["U_full"]        = U_full          # equilibrium (min-F) branch
    res["U_unbridged"]   = U_unbr
    res["U_bridged"]     = U_brid
    res["branch"]        = branch          # 1 unbridged, 2 bridged, 0 none
    res["phi_surf"]      = Phi_surf
    res["success"]       = success
    res["n_d_failed"]    = N_D - n_ok
    res["runtime_min"]   = (time() - t0) / 60.0
    res["profile_d"]     = profile_d
    res["profile_rho"]   = profile_rho
    res["profile_z"]     = profile_z
    # res["profile_phi"]   = profile_phi
    res["params"]        = to_dict(PARAMS)   # full parameter set, for provenance
    return res
end
