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
    # surface-to-surface gap grid (as factors of the sphere radius R).
    #
    # The GAP is the parallelized axis in this job (see fisher_two_spheres.jl):
    # alphas = [d_grid...; d_ref], so nparams = n_d + 1.  sim/sweep.jl:32 runs
    # keya = 1 : ncpus * div(nparams, ncpus), so SET n_d = ncpus - 1 (or any n_d
    # with (n_d + 1) % ncpus == 0) -- otherwise the tail of the gap grid, which
    # ends with the reference, is silently dropped.  target! warns and the interp
    # stage errors out if the reference is missing, but it is cheaper to get right.
    d_min_factor::Float64 = 0.0      # near-contact end of the d-sweep
    d_max_factor::Float64 = 0.5
    n_d::Int              = 17       # = ncpus - 1
    d_ref_factor::Float64 = 2.0      # reference gap (spheres far apart): U = F(d) - F_ref
    # spatial discretisation of the (rho, z) half-plane
    Nrho::Int             = 200
    dz_target::Float64    = 0.03*20
    rho_max::Float64      = 18.0*20
    z_buf::Float64        = 3.0*20
    # Newton / fixed-point solve
    f_tol::Float64        = 1e-6
    maxiter::Int          = 200
    # Robin-BC regularization (pins the adsorption branch).  KEEP AT 0 for the
    # arclength pipeline: bc_reg > 0 adds a surface field whose free-energy
    # contribution, bc_reg * Int[phi_s log phi_s + (1-phi_s) log(1-phi_s)] dA, is
    # NOT in excess_free_energy_two_spheres.  The omission is branch-dependent
    # (the thin film and the bridge sit at very different phi_s), so at bc_reg =
    # 5e-2 it inverts the ordering of the two stable branches and destroys the
    # first-order crossing.  See test/test_bvp_arclength.jl.
    bc_reg::Float64       = 0.00
    # Pseudo-arclength continuation in phi_inf (src/bvp_arclength.jl).
    arc_dlambda0::Float64 = 0.0      # phi_inf progress on the first step; 0 => range/50
    arc_max_steps::Int    = 400      # continuation-step budget per gap
    # Branch labelling: a gap centre is "bridged" once its phi has risen this
    # FRACTION of the way from the dilute to the dense binodal (see
    # bridge_threshold).  It is a fraction, not an absolute phi, because the
    # dense binodal is nowhere near 1 for an asymmetric fb: at nb = 4000 it is
    # 0.0694, so the old absolute `bridge_phi = 0.5` could never be reached by
    # any state the solver can produce and labelled every branch "thin".
    bridge_frac::Float64  = 0.5
    # --- optional warm start for the bridged branch (sweep B) ----------------
    # Empty (the default) => nothing below runs and the sweep behaves exactly as
    # it always has: sweep B cold-starts from the analytic dense seed.  Set it to
    # a `bridged_seed_*.json` written by an earlier run -- or to the directory
    # holding them, in which case each record picks the largest seeded gap that
    # does not exceed its own -- and fisher_two_spheres.jl remaps that solution
    # onto this record's mesh and hands it to Newton instead, keeping the
    # analytic seed as the fallback.  The analytic seed is a cold guess and stops
    # converging once the gap is a sizeable fraction of R: at R = 130 it failed at
    # EVERY gap of data/20260830/test3, so the only bridged branches in that run
    # were the ones sweep A's folds happened to reach, and past d ~ 33 (10 nm) the
    # records show no bridged branch at all.  The branch is there -- a bridged
    # solution remapped from a shorter gap converges in ~2 s at every d out to
    # 0.5 R.  See the seed block in fisher_two_spheres.jl.
    seed_path::String     = "data/20260912/test3"
    # Write this record's bridged solution to
    # output_dir/bridged_seed_d<d>_b<keyb>.json so a later run can seed from it.
    # ~3 MB per gap that has a bridged branch; nothing reads these but the seed
    # loader, and the downstream stages match output_cpu<n>.json exactly.
    seed_save::Bool       = true
end

# --- Run mode ----------------------------------------------------------------
# RUN_MODE (env, set by run_all.sh) picks between the cluster production sweep and
# the local smoke test.  Both configurations live HERE (s_values, OUT_ROOT) and in
# the run-length knobs of rouse_mfpt_attractive.jl; nothing else differs.
const RUN_MODE = get(ENV, "RUN_MODE", "production")
RUN_MODE in ("production", "smoke") || error("RUN_MODE=$RUN_MODE: use production | smoke")
const SMOKE    = RUN_MODE == "smoke"
const OUT_ROOT = SMOKE ? "data/20261009/smoke" : "data/20261009"

# The s this process runs: ROUSE_S (env, set per chain by run_all.sh), which must be
# one of s_values -- the launcher reads the list from here, so it cannot invent an s.
# Unset (BVP sweep, demo, interactive use) -> s_values[1].
function _select_s(s_values)
    haskey(ENV, "ROUSE_S") || return s_values[1]
    s = parse(Float64, ENV["ROUSE_S"])
    s in s_values || error("ROUSE_S=$s is not in PARAMS.rouse.s_values = $s_values")
    return s
end

# --- Rouse chain + thermal environment ---------------------------------------
# `s`,`a` are `beadSep`,`kuhnLength` in the FPT script; the entropic spring in the
# two-sphere P_eq uses 3(d+R)^2 / (2 s a^2).  D_eff = 2 kBT/friction is the
# diffusion constant of the end-to-end coordinate (two independent end beads).
Base.@kwdef struct RouseParams
    s_values::Vector{Float64} = SMOKE ? [10.0, 160.0, 310.0] : collect(10.0:20.0:310.0)
    s::Float64        = _select_s(s_values)   # beads between the two probe monomers (beadSep); derived, do not set
    a::Float64        = 20     # Kuhn / bond length (lattice sites) 20 b = 6 nm, 100 b = 30 nm
    kBT::Float64      = 1.0
    friction::Float64 = 1.3477088948787062e-08   # s / sites^2
end
D_eff(r::RouseParams) = 2 * r.kBT / r.friction

# --- FPT solver (downstream; carried here so PARAMS is the single source) -----
Base.@kwdef struct FPTParams
    potential::Symbol       = :bvp     # :bvp (nonlinear BVP U_full) or :linear (closed-form well)
    Ngrid::Int              = 1500     # grid points on [R, x_ridge] (escape / return problems)
    Nlocal::Int             = 100_000  # dense grid for the local-concentration integral
    ridge_scan_max::Float64 = 200     # x-range scanned to locate the barrier ridge
    ridge_scan_n::Int       = 8000
    plateau_tol::Float64    = 1e-3     # |dU_att/dx| (kBT/site) at which the largest-phi well has flattened; x_div = 1.2 x there
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
    phi_frac_max::Float64           = 0.95
    n_phi::Int                      = 32
    omega_b_values::Vector{Float64} = [-0.008]
end

# --- Run metadata ------------------------------------------------------------
# Where this job's spec lives, where its output was written, and how many CPUs it
# ran on.  MUTABLE so sim/sweep.jl can stamp its CLI args into `PARAMS.sim` at run
# time (the const PARAMS binding is untouched; only this inner object is set); the
# stamped values ride along in every record's `params` dict for provenance.
#
# `output_dir` is the per-s handoff between stages: bvp_potential_interp.jl and
# fisher_rate_from_fpt_solver.jl write there, rouse_mfpt_attractive.jl reads them.
# It is OUT_ROOT/s<s>, set in FisherParams below.  The BVP sweep is s-independent
# and is read, never written, from `bvp_dir` (sweep.jl overwrites output_dir
# in-process on a forward run; ncpus likewise).
Base.@kwdef mutable struct SimulationParams
    spec_path::Union{Nothing, String}   = ""
    bvp_dir::String                     = "data/20260927/s10"
    output_dir::Union{Nothing, String}  = OUT_ROOT
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
    surface::Parameters       = Parameters(nb = 4000.0, chi_b = 0.5280059665030825, 
                                           omega_b = -0.008,
                                           kappa = 5.19, phi_inf = 0.09, L = 100.0,
                                           R = 130, eps = 1e-8,
                                           conserve_volume = false,
                                           parallel_plates = false)
    bvp::BVPParams            = BVPParams()
    rouse::RouseParams        = RouseParams()
    fpt::FPTParams            = FPTParams()
    sweep::SweepParams        = SweepParams()
    sim::SimulationParams     = SimulationParams(output_dir = joinpath(OUT_ROOT, "s$(Int(rouse.s))"))
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
#
# BOTH roots are returned, and callers that want the dense phase MUST use the
# second one.  fb is only symmetric about phi = 1/2 at nb = 1, so `1 - phi_dilute`
# is not the dense binodal for any realistic nb: at nb = 4000, chi_b = 0.528 the
# pair is (3.23e-4, 0.0694), and 1 - 3.23e-4 = 0.99968 lands 14x past the dense
# branch, outside the two-phase region entirely.  Seeding a capillary bridge
# there gives Newton a slab of near-pure polymer to relax and it simply fails --
# which is what silently disabled sweep B in fisher_two_spheres.jl for every gap.
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
    return dilute(mu), dense(mu)
end

# Absolute gap-centre phi above which a branch counts as bridged, placed
# `bvp.bridge_frac` of the way from the dilute to the dense binodal.
function bridge_threshold(P::FisherParams)
    lo, hi = binodal(P.surface.chi_b; nb = P.surface.nb)
    return lo + P.bvp.bridge_frac * (hi - lo)
end

# The phi_inf sweep grid: explicit list if given, else a binodal-fraction range.
function phi_inf_grid(P::FisherParams)
    s = P.sweep
    isempty(s.phi_inf_values) || return copy(s.phi_inf_values)
    b, _ = binodal(P.surface.chi_b; nb = P.surface.nb)
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


