# Monotone, shape-preserving interpolation of a low-resolution BVP pair potential.
# Supports both d-major and phi-major BVP layouts, plus analytic :linear mode.
#
# JOB SCRIPT (spec/jobs/20260720/): a self-contained copy of the viz_util library,
# run as `julia --project=@. spec/jobs/20260807/pseudoarc/bvp_potential_interp.jl`.
# Beyond the interpolation helpers it drives the interpolation STAGE: it reads its
# local param.jl (`PARAMS`), aggregates the fisher_two_spheres BVP sweep from
# PARAMS.sim.output_dir (handling both d-major and phi-major layouts), builds
# U_full_interp / U_linear / U_tot_full / U_tot_linear on a common
# PARAMS.fpt.Ngrid-point x-grid, and writes them to PARAMS.sim.output_dir /
# PARAMS.sim.interp_file. No paths or filenames are baked into this code -- they
# all come from param.jl (found next to this file).
#
# STANDALONE :linear MODE: When PARAMS.fpt.potential == :linear the whole pipeline
# is closed-form, so this stage never touches fisher_two_spheres' output: the BVP
# aggregation is skipped, phi_inf comes straight from phi_inf_grid(PARAMS), and every
# BVP-derived quantity (d_grid, U_full*, phi_surf*, U_tot_full) is written as a NaN
# placeholder. Only U_linear / U_tot_linear carry real values. Run it with no
# output_cpu*.json present and it produces a valid interp_file from param.jl alone.
#
# The two-sphere BVP dashboard (dashboard_fisher_two_spheres.py) emits, for each
# bulk volume fraction phi_inf, the equilibrium pair potential U_full(d) sampled on
# a COARSE surface-to-surface gap grid d_grid (N_D ~ 30 nodes).  The FPT script
# (fisher_rate_from_fpt_solver.jl) needs the attractive potential as a smooth
# FUNCTION evaluated on a FINE grid (Ngrid ~ 1500) that also feeds finite
# differences into the backward-FPE generator -- so the resampled potential must have
# a continuous first derivative or the drift term spikes.
#
# Design (why not a plain spline):
#   * The raw BVP tail carries ~1e-6 non-monotone solver noise around U = 0 (the
#     film free energy is at its floor there).  A faithful interpolant -- even a
#     monotonicity-PRESERVING one -- reproduces those dips as wiggles.  So we first
#     project the nodes onto the nearest NON-DECREASING sequence (pool-adjacent-
#     violators, an O(n) L2 isotonic regression); it only touches the noisy tail.
#   * A natural cubic spline (C2) would OVERSHOOT at the deep-well / flat-tail
#     asymmetry and dip below the asymptote -> non-monotone.  We instead use a
#     Fritsch-Butland monotone cubic Hermite (PCHIP, as in MATLAB pchip / scipy
#     PchipInterpolator): C1, shape-preserving (NO overshoot => no wiggles),
#     monotone on monotone data.  We trade C2 smoothness (which reintroduces
#     wiggles here) for a guaranteed monotone, continuous first derivative.
#
# Exposes:
#   isotonic_nondecreasing(y)      -- L2 projection onto non-decreasing sequences
#   PCHIP(x, y; extrap_left)       -- monotone cubic Hermite; callable; deriv/deriv2
#   interp_U_full(d_grid, U_full)  -- cleaned PCHIP of one BVP U(d) row
#   make_U_att_bvp(d_grid, U, R)   -- drop-in U_ads replacement: x -> U_att(x-R)
#
# NEAR-CONTACT EXTENSION (PARAMS.fpt.extrap_small_d, forwarded as extrap_left):
# the BVP's first FINITE gap is typically d > 0 (Newton fails near contact) while
# the x-grid starts at x = R.  :flat holds that first value over the unsolved
# gaps (the historical plateau); :slope instead continues the interpolant's slope
# there as a straight line, so the well keeps deepening down to contact.

using Droplet          # linear-response screening length (d2fb_dphi2, with)
import JSON

# --- Isotonic regression (pool-adjacent-violators) ---------------------------
# Closest non-decreasing sequence to y in L2 (equal weights).  Merges any block
# whose running mean would violate monotonicity; O(n) via a value/weight/count
# stack.  Only the ~1e-6 non-monotone tail is modified.
function isotonic_nondecreasing(y::AbstractVector)
    vals = Float64[]; wts = Float64[]; cnts = Int[]
    @inbounds for yi in y
        v = float(yi); w = 1.0; c = 1
        while !isempty(vals) && vals[end] > v
            pv = pop!(vals); pw = pop!(wts); pc = pop!(cnts)
            v = (pv * pw + v * w) / (pw + w); w += pw; c += pc
        end
        push!(vals, v); push!(wts, w); push!(cnts, c)
    end
    out = Vector{Float64}(undef, length(y)); i = 1
    @inbounds for k in eachindex(vals)
        for _ in 1:cnts[k]
            out[i] = vals[k]; i += 1
        end
    end
    return out
end

# --- Monotone cubic Hermite (PCHIP) ------------------------------------------
struct PCHIP
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}   # node slopes (Fritsch-Butland)
    extrap_left::Symbol  # below x[1]: :flat (hold y[1]) or :slope (line of slope m[1])
end

# Shape-preserving one-sided endpoint slope (scipy PchipInterpolator._edge_case):
# a non-centred 3-point estimate, clamped so it cannot flip sign or overshoot.
function _pchip_edge(h0, h1, d0, d1)
    d = ((2h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    if sign(d) != sign(d0)
        d = 0.0
    elseif sign(d0) != sign(d1) && abs(d) > 3abs(d0)
        d = 3d0
    end
    return d
end

function _pchip_slopes(x, y)
    n = length(x)
    h = diff(x)
    δ = diff(y) ./ h
    m = zeros(n)
    @inbounds for i in 2:n-1
        if δ[i-1] * δ[i] > 0                     # same sign => monotone through node
            w1 = 2h[i] + h[i-1]                  # Fritsch-Butland weighted harmonic
            w2 = h[i] + 2h[i-1]                  # mean of the two secants
            m[i] = (w1 + w2) / (w1 / δ[i-1] + w2 / δ[i])
        else
            m[i] = 0.0                           # local extremum / flat => zero slope
        end
    end
    if n >= 3
        m[1] = _pchip_edge(h[1], h[2], δ[1], δ[2])
        m[n] = _pchip_edge(h[n-1], h[n-2], δ[n-1], δ[n-2])
    elseif n == 2
        m[1] = δ[1]; m[2] = δ[1]
    end
    return m
end

function PCHIP(x::AbstractVector, y::AbstractVector; extrap_left::Symbol = :flat)
    xs = Float64.(collect(x)); ys = Float64.(collect(y))
    @assert length(xs) == length(ys) "x and y must have equal length"
    @assert length(xs) >= 2 "need at least 2 nodes"
    @assert issorted(xs) && all(diff(xs) .> 0) "x must be strictly increasing"
    @assert extrap_left in (:flat, :slope) "extrap_left must be :flat or :slope"
    return PCHIP(xs, ys, _pchip_slopes(xs, ys), extrap_left)
end

@inline _seg(p::PCHIP, xq) = clamp(searchsortedlast(p.x, xq), 1, length(p.x) - 1)

# Cubic Hermite basis and its first two derivatives (argument t in [0,1]).
@inline _h00(t)  = (1 + 2t) * (1 - t)^2
@inline _h10(t)  = t * (1 - t)^2
@inline _h01(t)  = t^2 * (3 - 2t)
@inline _h11(t)  = t^2 * (t - 1)
@inline _h00d(t) = 6t^2 - 6t
@inline _h10d(t) = 3t^2 - 4t + 1
@inline _h01d(t) = -6t^2 + 6t
@inline _h11d(t) = 3t^2 - 2t
@inline _h00dd(t) = 12t - 6
@inline _h10dd(t) = 6t - 4
@inline _h01dd(t) = -12t + 6
@inline _h11dd(t) = 6t - 2

# Value.  Above d_max the pair potential is already ~0, so holding the endpoint
# value keeps U_att bounded and finite for the FPT integration.  Below x[1] the
# extension is a choice: the BVP's first finite gap can sit well above d = 0
# (Newton fails near contact) while the FPT x-grid starts at x = R, so
#   :flat  -- hold y[1]; the well has a plateau floor over the unsolved gaps
#   :slope -- continue the interpolant's own left-end slope m[1] as a straight
#             line (C1 across x[1]), i.e. the last resolved trend carries on
#             instead of flattening.
function (p::PCHIP)(xq::Real)
    if xq <= p.x[1]
        return p.extrap_left === :slope ? p.y[1] + p.m[1] * (xq - p.x[1]) : p.y[1]
    end
    xq >= p.x[end] && return p.y[end]
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return p.y[i] * _h00(t) + h * p.m[i] * _h10(t) +
           p.y[i+1] * _h01(t) + h * p.m[i+1] * _h11(t)
end

function deriv(p::PCHIP, xq::Real)
    xq <= p.x[1]   && return p.extrap_left === :slope ? p.m[1] : 0.0
    xq >= p.x[end] && return 0.0
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return (p.y[i] * _h00d(t) + h * p.m[i] * _h10d(t) +
            p.y[i+1] * _h01d(t) + h * p.m[i+1] * _h11d(t)) / h
end

function deriv2(p::PCHIP, xq::Real)
    (xq <= p.x[1] || xq >= p.x[end]) && return 0.0
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return (p.y[i] * _h00dd(t) + h * p.m[i] * _h10dd(t) +
            p.y[i+1] * _h01dd(t) + h * p.m[i+1] * _h11dd(t)) / h^2
end

# --- Adapters for the fisher-rate pipeline -----------------------------------
# One BVP U(d) row -> a monotone, wiggle-free interpolant on the gap coordinate d.
# extrap_left picks how the gaps below the first node are filled (see PCHIP).
interp_U_full(d_grid, U_full; extrap_left::Symbol = :flat) =
    PCHIP(d_grid, isotonic_nondecreasing(collect(U_full)); extrap_left = extrap_left)

# Drop-in replacement for the linear-response `U_ads` term in make_U_tot: builds
# the interpolant once for a given phi_inf's BVP row and returns U_att as a
# function of the FPT coordinate x (= surface gap d + sphere radius R).
function make_U_att_bvp(d_grid, U_full, R; extrap_left::Symbol = :flat)
    itp = interp_U_full(d_grid, U_full; extrap_left = extrap_left)
    return x -> itp(x - R)
end

# ── Interpolation stage ──────────────────────────────────────────────────────
# Read the fisher_two_spheres BVP sweep, build the interpolated nonlinear potential
# U_full_interp and the linear-response approximation U_linear on a common
# Ngrid-point x-grid (x = surface gap d + R), add the radial Jacobian + Rouse
# spring to get U_tot_full / U_tot_linear (ported from build_U_tot in
# fisher_rate_from_fpt_solver.jl), and write everything back beside the input.
# Every path/filename comes from the local param.jl -- nothing is hardcoded here.

include(joinpath(@__DIR__, "param.jl"))

_f(x)    = x === nothing ? NaN : Float64(x)              # JSON null -> NaN
_farr(v) = Float64[_f(x) for x in v]

_f(x)    = x === nothing ? NaN : Float64(x)              # JSON null -> NaN
_farr(v) = Float64[_f(x) for x in v]

data_dir = PARAMS.sim.output_dir
R     = PARAMS.surface.R
h0    = PARAMS.surface.omega_b        # adsorption strength for the linear model
s, a  = PARAMS.rouse.s, PARAMS.rouse.a
Ngrid = PARAMS.fpt.Ngrid
# :flat or :slope -- how U_full is continued below the smallest SOLVED gap (the
# first finite node of the row, which near contact is usually not d = 0).
extrap_small_d = PARAMS.fpt.extrap_small_d

linear_only = PARAMS.fpt.potential == :linear   # analytic run; no BVP files needed

# --- 1. gather the BVP sweep (skipped entirely for the analytic :linear model) ---
# TWO input layouts are accepted, distinguished by whether a record carries a
# scalar "d":
#
#   phi_inf-major (spec/fisher_two_spheres.jl): one record per phi_inf, each
#     already holding U_full(d) and phi_surf(d) over a shared d_grid.  Records
#     concatenate in cpu order.
#
#   d-major (spec/jobs/20260807/fisher_two_spheres.jl): one record per GAP, each
#     holding F_eq and phi_surf over the whole phi_inf grid, plus one record with
#     is_ref = true at the reference gap.  That spec cannot form U itself -- the
#     reference lives in a different shard -- so the subtraction
#     U(d, phi_inf) = F_eq(d, phi_inf) - F_eq(d_ref, phi_inf) happens here, where
#     every d is in hand, and the result is transposed into the phi_inf-major
#     rows the rest of this file (and the FPT stage) already expects.

if linear_only
    isdir(data_dir) || mkpath(data_dir)                  # only need somewhere to write
    phi_infs      = phi_inf_grid(PARAMS)
    d_grid        = Float64[]                             # no BVP gap grid
    U_full_rows   = [Float64[] for _ in phi_infs]         # placeholder (no BVP U_full)
    phi_surf_rows = [Float64[] for _ in phi_infs]         # placeholder (no BVP phi_surf)
    by_cpu        = Dict{Int,NamedTuple}()
    cpus          = Int[]
else
    isdir(data_dir) || error("PARAMS.sim.output_dir is not a directory: $data_dir")

    files = sort(filter(f -> occursin(r"output_cpu\d+\.json$", f),
                        readdir(data_dir; join = true)))
    isempty(files) && error("no output_cpu*.json under $data_dir")

    # flatten every cpu file to (cpu, keya, record), picking the omega_b branch this
    # run's linear model is built at
    all_recs = Tuple{Int,Int,Dict{String,Any}}[]
    for f in files
        data = JSON.parsefile(f)
        isempty(data) && continue                            # a CPU that got no shard
        cpu = parse(Int, match(r"output_cpu(\d+)\.json$", f).captures[1])
        for ka in sort(collect(keys(data)); by = x -> parse(Int, x))
            recs = data[ka]
            kbs  = collect(keys(recs))
            rec  = recs[kbs[argmin([abs(_f(recs[kb]["omega_b"]) - h0) for kb in kbs])]]
            push!(all_recs, (cpu, parse(Int, ka), rec))
        end
    end
    isempty(all_recs) && error("all output_cpu*.json under $data_dir were empty")

    is_d_major = haskey(all_recs[1][3], "d")

    # The sweep solved one specific sphere radius, and every record records it.  R sets
    # where the well lands on the x-grid (x = d + R) and enters the Rouse spring, so
    # editing PARAMS.surface.R without re-running the BVP silently relocates a potential
    # that was computed for different spheres -- and since this stage rebuilds x_grid
    # from the same edited R, nothing downstream can notice.
    let R_sweep = _f(get(all_recs[1][3], "R", NaN))
        isnan(R_sweep) || isapprox(R_sweep, R; rtol = 1e-8) ||
            error("PARAMS.surface.R = $R but the BVP sweep in $data_dir ran at R = $R_sweep.\n" *
                  "Re-run the sweep at this R, or set R = $R_sweep in param.jl.")
    end

    if is_d_major
        refs = [r for (_, _, r) in all_recs if get(r, "is_ref", false) === true]
        isempty(refs) && error("""
            d-major sweep has no reference record (is_ref = true).  sim/sweep.jl runs
            keya = 1 : ncpus * div(nparams, ncpus), and the reference gap is the LAST
            entry of `alphas`, so it is dropped unless ncpus divides nparams = n_d + 1.
            Set n_d = ncpus - 1 and re-run the sweep.""")
        length(refs) == 1 || error("$(length(refs)) reference records; expected exactly 1")
        F_ref = _farr(refs[1]["F_eq"])

        body = sort([r for (_, _, r) in all_recs if get(r, "is_ref", false) !== true];
                    by = r -> _f(r["d"]))
        isempty(body) && error("d-major sweep has only the reference record")
        d_grid = [_f(r["d"]) for r in body]

        phi_infs = _farr(all_recs[1][3]["phi_inf"])
        nphi = length(phi_infs)
        length(F_ref) == nphi ||
            error("reference F_eq has $(length(F_ref)) entries; phi_inf has $nphi")
        for r in body
            pv = _farr(r["phi_inf"])
            (length(pv) == nphi && all(isapprox.(pv, phi_infs; atol = 1e-12, rtol = 1e-9))) ||
                error("record at d = $(_f(r["d"])) has a different phi_inf axis")
        end

        # transpose (d, phi_inf) -> (phi_inf, d) and subtract the reference
        U_full_rows   = [Float64[_f(r["F_eq"][i]) - F_ref[i]  for r in body] for i in 1:nphi]
        phi_surf_rows = [Float64[_f(r["phi_surf"][i])         for r in body] for i in 1:nphi]
        by_cpu        = Dict{Int,NamedTuple}()
        cpus          = Int[]
    else
        by_cpu = Dict{Int,NamedTuple}()
        for (cpu, _, rec) in all_recs
            e = get(by_cpu, cpu, (phi_inf = Float64[],
                                  U_full = Vector{Vector{Float64}}(),
                                  phi_surf = Vector{Vector{Float64}}(),
                                  d_grid = Float64[]))
            push!(e.phi_inf,  _f(rec["phi_inf"]))
            push!(e.U_full,   _farr(rec["U_full"]))
            push!(e.phi_surf, _farr(rec["phi_surf"]))
            isempty(e.d_grid) && append!(e.d_grid, _farr(rec["d_grid"]))
            by_cpu[cpu] = e
        end
        cpus          = sort(collect(keys(by_cpu)))
        phi_infs      = vcat((by_cpu[c].phi_inf  for c in cpus)...)
        U_full_rows   = vcat((by_cpu[c].U_full   for c in cpus)...)
        phi_surf_rows = vcat((by_cpu[c].phi_surf for c in cpus)...)
        d_grid        = by_cpu[cpus[1]].d_grid
    end
end

# --- 2. reconcile the aggregated phi_inf grid with PARAMS ---
# phi_inf_grid(PARAMS) SELECTS from the swept phi_inf axis; it no longer has to be the
# whole of it.  By this point both layouts (and :linear) have been reduced to phi-major
# rows, so one selection covers all three.  The BVP sweep is expensive and completely
# s-independent, which is the point: sweep_rouse_s.jl reuses one sweep, asking only for
# the three phi_inf bracketing its target c and varying the Rouse chain against them.
# Requiring exact equality left "move the output_cpu*.json out of the way" as the only
# way to do that.  A requested value the sweep never computed is still an error --
# interpolating BVP rows across phi_inf is not something this stage does.
grid = phi_inf_grid(PARAMS)
sel  = [findfirst(p -> isapprox(p, g; atol = 1e-9, rtol = 1e-6), phi_infs) for g in grid]
if any(isnothing, sel)
    error("phi_inf_grid(PARAMS) is not covered by the swept phi_inf axis.\n" *
          "  missing: $(grid[isnothing.(sel)])\n  swept:   $phi_infs\n" *
          "Set PARAMS.sweep.phi_inf_values to values the sweep actually computed " *
          "(or re-run the BVP sweep at the values you want).")
end
length(sel) == length(phi_infs) ||
    @info "phi_inf: selecting $(length(sel)) of $(length(phi_infs)) swept values"
phi_infs      = phi_infs[sel]
U_full_rows   = U_full_rows[sel]
phi_surf_rows = phi_surf_rows[sel]

# --- 3. common x-grid; Rouse end-to-end scale, as in the FPT return domain ---
x_max     = 4 * a * sqrt(s)
xg        = collect(range(R, x_max; length = Ngrid))
jac(x)    = -2 * log(x)                 # radial Jacobian (build_U_tot)
spring(x) = 3 * x^2 / (2 * s * a^2)     # Rouse spring    (build_U_tot)

U_full_interp   = Vector{Vector{Float64}}()
phi_surf_interp = Vector{Vector{Float64}}()
U_linear        = Vector{Vector{Float64}}()
U_tot_full      = Vector{Vector{Float64}}()
U_tot_linear    = Vector{Vector{Float64}}()

for (i, phi) in enumerate(phi_infs)
    # Nonlinear BVP potential + gap-facing surface concentration.  For the analytic
    # :linear run there is no BVP data, so these stay NaN placeholders (all the
    # :linear fisher stage consumes is U_tot_linear + its own analytic profile).
    if linear_only
        Uf = fill(NaN, Ngrid)
        Ps = fill(NaN, Ngrid)

        p   = Droplet.with(PARAMS.surface; phi_inf = phi)     # omega_b (=h0) preserved
        psi, _, _ = psi_linear(xg, p)                      # screened adsorption tail
        Ps        = phi .+ psi                       # far-field value + linear tail
    else
        # Nonlinear BVP potential: PCHIP over the FINITE (d, U_full) nodes only, so a
        # NaN at some d drops that node instead of the whole phi_inf row.
        row = U_full_rows[i]
        fin = isfinite.(row)
        if count(fin) >= 2
            itp = interp_U_full(d_grid[fin], row[fin]; extrap_left = extrap_small_d)
            Uf  = [itp(x - R) for x in xg]
        else
            Uf  = fill(NaN, Ngrid)
        end

        # BVP gap-facing surface concentration, interpolated onto the same x-grid so the
        # :bvp local_conc can average it under the nonlinear Boltzmann weight.  Plain
        # (shape-preserving) PCHIP over the FINITE nodes -- NO isotonic projection here:
        # unlike U_full, phi_surf falls off with d, so forcing non-decreasing is wrong.
        srow = phi_surf_rows[i]
        sfin = isfinite.(srow)
        if count(sfin) >= 2
            sitp = PCHIP(d_grid[sfin], srow[sfin])
            Ps   = [sitp(x - R) for x in xg]
        else
            Ps   = fill(NaN, Ngrid)
        end
    end

    # Linear-response adsorption well (U_ads from make_U_tot), at this phi_inf.
    p   = Droplet.with(PARAMS.surface; phi_inf = phi)     # omega_b (=h0) preserved
    lam = sqrt(d2fb_dphi2(phi, p) / p.kappa)
    Ul  = [-(4π * R^3 * h0^2 * exp(-lam * (x - R))) / (p.kappa * (lam + 1 / R) * x)
           for x in xg]

    push!(U_full_interp,   Uf)
    push!(phi_surf_interp, Ps)
    push!(U_linear,        Ul)
    push!(U_tot_full,      Uf .+ jac.(xg) .+ spring.(xg))
    push!(U_tot_linear,    Ul .+ jac.(xg) .+ spring.(xg))
end

# --- 4. write the output dict beside the input data ---
out = Dict{String,Any}(
    "params"        => to_dict(PARAMS),
    "phi_inf"       => phi_infs,
    "omega_b"       => h0,
    "R"             => R,
    "d_grid"        => d_grid,
    "x_grid"        => xg,
    "U_full_raw"    => U_full_rows,        # per (phi_inf, d); NaN where unsolved
    "phi_surf_raw"  => phi_surf_rows,      # per (phi_inf, d)
    "U_full_interp" => U_full_interp,      # per (phi_inf, x_grid)
    "phi_surf_interp" => phi_surf_interp,  # per (phi_inf, x_grid): BVP surface conc(x-R)
    "U_linear"      => U_linear,           # per (phi_inf, x_grid)
    "U_tot_full"    => U_tot_full,         # per (phi_inf, x_grid)
    "U_tot_linear"  => U_tot_linear,       # per (phi_inf, x_grid)
    "layout"        => is_d_major ? "d_major" : "phi_major",
)

# Provenance / diagnostics that only one of the two layouts has.
if is_d_major
    # The spinodal map: phi_inf at each limit point of the phi_inf curve, per
    # gap.  Two folds bracket the hysteresis window where the thin film and the
    # capillary bridge coexist; the first-order crossing lies between them.
    out["folds_by_d"]  = Dict(string(_f(r["d"])) => _farr(get(r, "folds", []))
                              for r in body)
    out["n_branch"]    = Dict(string(_f(r["d"])) => Int[Int(v) for v in r["n_branch"]]
                              for r in body)
    out["F_eq_raw"]    = Dict(string(_f(r["d"])) => _farr(r["F_eq"]) for r in body)
    out["F_ref"]       = F_ref
    out["d_ref"]       = _f(refs[1]["d"])
else
    out["input_by_cpu"] = Dict(string(c) => Dict(
                               "phi_inf"  => by_cpu[c].phi_inf,
                               "U_full"   => by_cpu[c].U_full,
                               "phi_surf" => by_cpu[c].phi_surf,
                               "d_grid"   => by_cpu[c].d_grid) for c in cpus)
end
outpath = joinpath(data_dir, PARAMS.sim.interp_file)
open(outpath, "w") do io
    JSON.print(io, out)
end
println("layout=$(is_d_major ? "d_major" : "phi_major"), $(length(all_recs)) record(s), " *
        "$(length(phi_infs)) phi_inf x $(length(d_grid)) d; " *
        "x-grid $(Ngrid) pts on [$(round(R; digits = 3)), $(round(x_max; digits = 3))]")
println("saved -> $outpath")
