# Master parameter file -- the single source of truth for the two-sphere Fisher-
# rate pipeline.  Every stage is meant to `include` this file and read `PARAMS`,
# so the stages stay in lockstep and every JSON dump can embed `to_dict(PARAMS)`
# for provenance.  Currently consumed by (all in spec/jobs/20260720/):
#   fisher_two_spheres.jl          (BVP sweep)
#   bvp_potential_interp.jl        (interpolation + linear-response + U_tot)
#   fisher_rate_from_fpt_solver.jl (FPT solve + Fisher information rate)
#
# The configuration is partitioned into the four physical/numerical stages the
# request calls out, plus the swept axes:
#   surface  -- Flory-Huggins bulk + Robin surface-adsorption field (Droplet.Parameters)
#   bvp      -- axisymmetric two-sphere Cahn-Hilliard BVP grid & Newton solve
#   rouse    -- the Rouse chain (bead count, Kuhn length) + thermal environment
#   fpt      -- the 1D backward-FPE first-passage solver
#   sweep    -- the swept axes (phi_inf signal, omega_b field)
#   sim      -- run metadata (spec path, output dir, ncpus); filled by sim/sweep.jl
#
# Everything is plain data: `to_dict(PARAMS)` renders the whole tree as nested
# String-keyed dicts ready for JSON.print, and `save_params(PARAMS, path)` dumps
# it standalone.

using Droplet: Parameters
import JSON

# --- BVP solver: two-sphere axisymmetric BVP grid, Newton solve, branch label -
Base.@kwdef struct BVPParams
    # surface-to-surface gap grid (as factors of the sphere radius R)
    d_min_factor::Float64 = 0.0      # near-contact end of the d-sweep
    d_max_factor::Float64 = 1.0
    n_d::Int              = 50
    d_ref_factor::Float64 = 8.0      # reference gap (spheres far apart): U = F(d) - F_ref
    # spatial discretisation of the (rho, z) half-plane
    Nrho::Int             = 200
    dz_target::Float64    = 0.05
    rho_max::Float64      = 15.0
    z_buf::Float64        = 3.0
    # Newton / fixed-point solve
    f_tol::Float64        = 1e-6
    maxiter::Int          = 200
    bc_reg::Float64       = 0     # Robin-BC regularization (pins the adsorption branch)
    # branch labelling: gap-centre phi above this => "bridged"
    bridge_phi::Float64   = 0.5
end

# --- Rouse chain + thermal environment ---------------------------------------
# `s`,`a` are `beadSep`,`kuhnLength` in the FPT script; the entropic spring in the
# two-sphere P_eq uses 3(d+R)^2 / (2 s a^2).  D_eff = 2 kBT/friction is the
# diffusion constant of the end-to-end coordinate (two independent end beads).
Base.@kwdef struct RouseParams
    s::Float64        = 150.0   # number of beads between the two probe monomers (beadSep)
    a::Float64        = 1.0     # Kuhn / bond length (kuhnLength)
    kBT::Float64      = 1.0
    friction::Float64 = 1.0
end
D_eff(r::RouseParams) = 2 * r.kBT / r.friction

# --- FPT solver (downstream; carried here so PARAMS is the single source) -----
Base.@kwdef struct FPTParams
    potential::Symbol       = :bvp     # :bvp (nonlinear BVP U_full) or :linear (closed-form well)
    Ngrid::Int              = 1500     # grid points on [R, x_ridge] (escape / return problems)
    Nlocal::Int             = 100_000  # dense grid for the local-concentration integral
    ridge_scan_max::Float64 = 25.0     # x-range scanned to locate the barrier ridge
    ridge_scan_n::Int       = 8000
end

# --- Swept axes --------------------------------------------------------------
# phi_inf is the sensed signal, omega_b the (usually fixed) adsorption field.
# phi_inf_values, if non-empty, is used verbatim (default: the single test point);
# otherwise `phi_inf_grid` builds `n_phi` points over
# [phi_frac_min, phi_frac_max] x dilute_binodal(chi_b) -- near-saturation, where
# bridging is active.  (The full study sweep is phi_inf_values = Float64[],
# phi_frac_min = 0.05, phi_frac_max = 0.90, n_phi = 9.)
Base.@kwdef struct SweepParams
    phi_inf_values::Vector{Float64} = [] # collect(range(0.006, stop=0.012, length=64))
    phi_frac_min::Float64           = 0.1
    phi_frac_max::Float64           = 0.7
    n_phi::Int                      = 256
    omega_b_values::Vector{Float64} = [-1.3]
end

# --- Run metadata ------------------------------------------------------------
# Where this job's spec lives, where its output was written, and how many CPUs it
# ran on.  MUTABLE so sim/sweep.jl can stamp its CLI args into `PARAMS.sim` at run
# time (the const PARAMS binding is untouched; only this inner object is set); the
# stamped values ride along in every record's `params` dict for provenance.
#
# `output_dir` is also the handoff between stages: sweep.jl WRITES the BVP output
# there, and downstream stages (viz_util/bvp_potential_interp.jl) READ it to find
# that data and write their own output alongside it.  Point it at the run to
# process (sweep.jl overwrites it in-process on a forward run; ncpus likewise).
Base.@kwdef mutable struct SimulationParams
    spec_path::Union{Nothing, String}   = ""
    output_dir::Union{Nothing, String}  = "data/20260804/run5"
    ncpus::Int          = 0    # placeholder; set to sweep.jl's ncpus (ARGS[3])
    # Output filenames for the downstream stages (written into output_dir).
    interp_file::String = "bvp_potential_interp.json"
    fisher_file::String = "fisher_rate_from_fpt_solver.json"
    one_sphere_file::String = "demo_bvp_profiles.json"
end

# --- Master container --------------------------------------------------------
# The kwdef defaults ARE this job's configuration; edit them (or construct a
# modified FisherParams) to reconfigure the whole pipeline from one place.
Base.@kwdef struct FisherParams
    surface::Parameters       = Parameters(nb = 3, chi_b = 2.0, omega_b = -0.8,
                                           kappa = 1.0, phi_inf = 0.09, L = 20.0,
                                           R = 10.0, eps = 1e-8,
                                           conserve_volume = false,
                                           parallel_plates = false)
    bvp::BVPParams            = BVPParams()
    rouse::RouseParams        = RouseParams()
    fpt::FPTParams            = FPTParams()
    sweep::SweepParams        = SweepParams()
    sim::SimulationParams     = SimulationParams()
end

# Dilute binodal phi_b(chi): small root of the symmetric-FH coexistence condition
# log(phi/(1-phi)) + nb*chi*(1-2phi) = 0.  Shared by the sweep grid and the
# dense-seed value 1 - phi_b (capillary condensation).
function dilute_binodal(chi::Real; nb::Real = 1.0)
    g(phi) = log(phi / (1.0 - phi)) + nb * chi * (1.0 - 2.0 * phi)
    lo, hi = 1e-8, 0.5 - 1e-9
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if g(mid) < 0.0; lo = mid; else; hi = mid; end
    end
    return 0.5 * (lo + hi)
end

# Bisection on a bracketed sign change; 80 halvings exhaust Float64.
function _bisect(f, a, b; n = 80)
    fa = f(a)
    for _ in 1:n
        m = 0.5 * (a + b)
        (fa * f(m) <= 0) ? (b = m) : (a = m; fa = f(m))
    end
    return 0.5 * (a + b)
end

# Common-tangent binodal (phi_dilute, phi_dense) of the asymmetric fb.  dfb_dphi
# rises, falls between the two spinodal points, then rises again, so each branch
# holds exactly one root of dfb_dphi = mu; bisect on mu until the two branches
# have equal grand potential fb - mu*phi.
function binodal(chi::Real; nb::Real = 1.0)
    p = Parameters(nb = nb, chi_b = chi, omega_b = 0.0, kappa = 1.0,
                   phi_inf = 0.1, L = 1.0, R = Inf, eps = 1e-16)
    lo_sp, hi_sp = compute_spinodal(p)
    lo_sp === nothing && return (nothing, nothing)
    dilute(mu) = _bisect(x -> dfb_dphi(x, p) - mu, 1e-300, lo_sp)
    dense(mu)  = _bisect(x -> dfb_dphi(x, p) - mu, hi_sp, 1 - 1e-16)
    gap(mu)    = (fb(dense(mu), p) - mu * dense(mu)) -
                 (fb(dilute(mu), p) - mu * dilute(mu))
    mu = _bisect(gap, dfb_dphi(hi_sp, p), dfb_dphi(lo_sp, p))
    return dilute(mu) #, dense(mu)
end

# The phi_inf sweep grid: explicit list if given, else a binodal-fraction range.
function phi_inf_grid(P::FisherParams)
    s = P.sweep
    isempty(s.phi_inf_values) || return copy(s.phi_inf_values)
    b = binodal(P.surface.chi_b; nb = P.surface.nb)
    return b .* collect(range(s.phi_frac_min, s.phi_frac_max; length = s.n_phi))
end

# --- JSON serialization ------------------------------------------------------
# Recursively render any section struct into JSON-ready primitives / dicts.
# Real covers Float64/Int/Bool; Symbol -> String; vectors map element-wise; any
# other value is a struct, rendered field-by-field.
_ser(x::Union{Real,AbstractString}) = x
_ser(x::Symbol)                     = String(x)
_ser(x::AbstractVector)             = [_ser(v) for v in x]
_ser(x)                             = Dict(String(f) => _ser(getfield(x, f))
                                           for f in fieldnames(typeof(x)))

to_dict(P::FisherParams) = _ser(P)

# Convenience: dump the parameter set to a standalone JSON file.
function save_params(P::FisherParams, path::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, to_dict(P), 4)
    end
    return path
end

# The configured source of truth for this job.
const PARAMS = FisherParams()
