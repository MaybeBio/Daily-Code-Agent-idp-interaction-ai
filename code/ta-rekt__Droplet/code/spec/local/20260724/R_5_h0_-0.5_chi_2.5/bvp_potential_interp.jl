# Faithful (non-monotone) cubic-Hermite interpolation of a low-res BVP pair potential.
#
# JOB SCRIPT (spec/jobs/20260720/): a self-contained copy of the viz_util library,
# run as `julia --project=@. spec/jobs/20260720/bvp_potential_interp.jl`.  Beyond
# the interpolation helpers it drives the interpolation STAGE: it reads its local
# param.jl (`PARAMS`), aggregates the fisher_two_spheres BVP sweep from
# PARAMS.sim.output_dir, builds U_full_interp / U_linear / U_tot_full / U_tot_linear
# on a common PARAMS.fpt.Ngrid-point x-grid, and writes them to
# PARAMS.sim.output_dir / PARAMS.sim.interp_file.  No paths or filenames are baked
# into this code -- they all come from param.jl (found next to this file).
#
# STANDALONE :linear MODE.  When PARAMS.fpt.potential == :linear the whole pipeline
# is closed-form, so this stage never touches fisher_two_spheres' output: the BVP
# aggregation is skipped, phi_inf comes straight from phi_inf_grid(PARAMS), and every
# BVP-derived quantity (d_grid, U_full*, phi_surf*, U_tot_full) is written as a NaN
# placeholder.  Only U_linear / U_tot_linear -- which is all the :linear branch of
# fisher_rate_from_fpt_solver.jl consumes -- carry real values.  Run it with no
# output_cpu*.json present and it produces a valid interp_file from param.jl alone.
#
# The two-sphere BVP dashboard (dashboard_fisher_two_spheres.py) emits, for each
# bulk volume fraction phi_inf, the equilibrium pair potential U_full(d) sampled on
# a COARSE surface-to-surface gap grid d_grid (N_D ~ 30 nodes).  The FPT script
# (fisher_rate_from_fpt_solver.jl) needs the attractive potential as a smooth
# FUNCTION evaluated on a FINE grid (Ngrid ~ 1500) that also feeds finite
# differences into the backward-FPE generator -- so the resampled potential must
# have a continuous, non-oscillating first derivative or the drift term spikes.
#
# Design (faithful interpolation, optional smoothing -- PARAMS.smoothing.method):
#   * :none (default) interpolates the raw BVP nodes FAITHFULLY with a cubic Hermite
#     (centred, unclamped slopes): no isotonic pre-projection, no monotonicity
#     constraint, so any real structure in U(d) survives.  A natural cubic spline
#     (C2) would OVERSHOOT the deep-well/flat-tail asymmetry; the C1 centred-slope
#     Hermite is far less prone to that.
#   * When the BVP rows are noisy the faithful interpolant reproduces the noise and
#     its DERIVATIVE oscillates.  The other methods SMOOTH instead of interpolate:
#     :spline is a cubic SMOOTHING spline (minimises Σ(y-g)² + λ∫g''²; C2, preserves
#     well depth, derivative smooth -- the recommended fix), and the kernel names
#     (:gaussian, :exponential, :epanechnikov, :triangular, :boxcar) do local kernel
#     regression (bandwidth + degree).  See SmoothingParams in param.jl.  Output
#     format is identical regardless -- only the resampled values change.
#
# Exposes:
#   CubicHermite(x, y)             -- cubic Hermite, centred slopes; callable; deriv/deriv2
#   SmoothSpline(x, y, λ)          -- cubic smoothing spline (penalty λ); callable; sderiv
#   kernel_smooth(x, y, xq; ...)   -- local kernel regression (kernel, bandwidth, degree)
#   smooth_row(d_nodes, y, dq)     -- resample one BVP row per PARAMS.smoothing dispatch
#   interp_U_full(d_grid, U_full)  -- CubicHermite of one BVP U(d) row
#   make_U_att_bvp(d_grid, U, R)   -- drop-in U_ads replacement: x -> U_att(x-R)

using Droplet          # linear-response screening length (d2fb_dphi2, with)
using LinearAlgebra    # dense solve for the cubic smoothing spline (Green-Silverman)
import JSON

# --- Cubic Hermite (centred slopes) ------------------------------------------
struct CubicHermite
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}   # node slopes (3-point centred finite differences)
end

# Interior slope: the two one-sided secants averaged, each weighted by the OPPOSITE
# interval (reduces to the centred difference on a uniform grid).  Endpoint slope:
# the non-centred 3-point (parabolic) estimate.  Neither is sign-clamped, so the
# interpolant is free to be non-monotone through a node.
_edge_slope(h0, h1, d0, d1) = ((2h0 + h1) * d0 - h0 * d1) / (h0 + h1)

function _hermite_slopes(x, y)
    n = length(x)
    h = diff(x)
    δ = diff(y) ./ h
    m = zeros(n)
    @inbounds for i in 2:n-1
        m[i] = (h[i] * δ[i-1] + h[i-1] * δ[i]) / (h[i-1] + h[i])
    end
    if n >= 3
        m[1] = _edge_slope(h[1], h[2], δ[1], δ[2])
        m[n] = _edge_slope(h[n-1], h[n-2], δ[n-1], δ[n-2])
    elseif n == 2
        m[1] = δ[1]; m[2] = δ[1]
    end
    return m
end

function CubicHermite(x::AbstractVector, y::AbstractVector)
    xs = Float64.(collect(x)); ys = Float64.(collect(y))
    @assert length(xs) == length(ys) "x and y must have equal length"
    @assert length(xs) >= 2 "need at least 2 nodes"
    @assert issorted(xs) && all(diff(xs) .> 0) "x must be strictly increasing"
    return CubicHermite(xs, ys, _hermite_slopes(xs, ys))
end

@inline _seg(p::CubicHermite, xq) = clamp(searchsortedlast(p.x, xq), 1, length(p.x) - 1)

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

# Value.  Constant (flat) extension outside the node range: below contact never
# happens (x >= R); above d_max the pair potential is already ~0, so holding the
# endpoint value keeps U_att bounded and finite for the FPT integration.
function (p::CubicHermite)(xq::Real)
    xq <= p.x[1]   && return p.y[1]
    xq >= p.x[end] && return p.y[end]
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return p.y[i] * _h00(t) + h * p.m[i] * _h10(t) +
           p.y[i+1] * _h01(t) + h * p.m[i+1] * _h11(t)
end

function deriv(p::CubicHermite, xq::Real)
    (xq <= p.x[1] || xq >= p.x[end]) && return 0.0
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return (p.y[i] * _h00d(t) + h * p.m[i] * _h10d(t) +
            p.y[i+1] * _h01d(t) + h * p.m[i+1] * _h11d(t)) / h
end

function deriv2(p::CubicHermite, xq::Real)
    (xq <= p.x[1] || xq >= p.x[end]) && return 0.0
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return (p.y[i] * _h00dd(t) + h * p.m[i] * _h10dd(t) +
            p.y[i+1] * _h01dd(t) + h * p.m[i+1] * _h11dd(t)) / h^2
end

# --- Cubic smoothing spline (Green & Silverman) ------------------------------
# Minimises Σ (y_i - g_i)² + λ ∫ g''(x)² dx over natural cubic splines: the C² curve
# that trades data fidelity against roughness.  Unlike a faithful interpolant it does
# NOT pass through the noisy BVP nodes, so its first derivative is smooth (the FPT
# drift term stops spiking).  λ -> 0 recovers interpolation; large λ -> a straight
# line.  `g` are the smoothed knot values, `M` the knot 2nd derivatives (natural:
# M[1] = M[n] = 0).
struct SmoothSpline
    x::Vector{Float64}
    g::Vector{Float64}
    M::Vector{Float64}
end

function SmoothSpline(x::AbstractVector, y::AbstractVector, λ::Real)
    xs = Float64.(collect(x)); ys = Float64.(collect(y))
    n  = length(xs)
    n < 3 && return SmoothSpline(xs, ys, zeros(n))          # too few nodes to penalise
    h  = diff(xs)
    R  = zeros(n-2, n-2); Q = zeros(n, n-2)                 # Green-Silverman bands
    for p in 1:n-2
        i = p + 1                                           # interior knot index (2..n-1)
        R[p, p] = (h[i-1] + h[i]) / 3
        p < n-2 && (R[p, p+1] = R[p+1, p] = h[i] / 6)
        Q[i-1, p] = 1 / h[i-1]
        Q[i,   p] = -(1 / h[i-1] + 1 / h[i])
        Q[i+1, p] = 1 / h[i]
    end
    γ = (R + λ * (Q'Q)) \ (Q' * ys)                         # interior 2nd derivatives
    g = ys .- λ * (Q * γ)                                   # smoothed knot values
    M = zeros(n); M[2:n-1] .= γ
    return SmoothSpline(xs, g, M)
end

# Natural cubic spline value / first derivative from (g, M); flat outside the knots.
function (s::SmoothSpline)(xq::Real)
    xq <= s.x[1]   && return s.g[1]
    xq >= s.x[end] && return s.g[end]
    i = clamp(searchsortedlast(s.x, xq), 1, length(s.x) - 1); h = s.x[i+1] - s.x[i]
    A = (s.x[i+1] - xq) / h; B = (xq - s.x[i]) / h
    return A * s.g[i] + B * s.g[i+1] + ((A^3 - A) * s.M[i] + (B^3 - B) * s.M[i+1]) * h^2 / 6
end

function sderiv(s::SmoothSpline, xq::Real)
    (xq <= s.x[1] || xq >= s.x[end]) && return 0.0
    i = clamp(searchsortedlast(s.x, xq), 1, length(s.x) - 1); h = s.x[i+1] - s.x[i]
    A = (s.x[i+1] - xq) / h; B = (xq - s.x[i]) / h
    return (s.g[i+1] - s.g[i]) / h + ((1 - 3A^2) * s.M[i] + (3B^2 - 1) * s.M[i+1]) * h / 6
end

# --- Kernel smoothing (local regression) -------------------------------------
# Nadaraya-Watson (degree 0) / local-linear (degree 1) regression of a coarse BVP row
# onto a query grid, with a selectable kernel and bandwidth (in gap-coordinate d
# units).  Local-linear removes the boundary bias that would otherwise flatten the
# deep well at contact.  A smooth kernel (Gaussian) gives the smoothest derivative;
# the exponential kernel's cusp at 0 is the noisiest -- offered for completeness.
_kernel(::Val{:gaussian}, u)     = exp(-0.5u^2)
_kernel(::Val{:exponential}, u)  = exp(-abs(u))
_kernel(::Val{:epanechnikov}, u) = (v = 1 - u^2;    v > 0 ? 0.75v : 0.0)
_kernel(::Val{:triangular}, u)   = (v = 1 - abs(u); v > 0 ? v     : 0.0)
_kernel(::Val{:boxcar}, u)       = abs(u) <= 1 ? 0.5 : 0.0
const _KERNELS = (:gaussian, :exponential, :epanechnikov, :triangular, :boxcar)

function kernel_smooth(xn::AbstractVector, yn::AbstractVector, xq::AbstractVector;
                       kernel::Symbol, bandwidth::Real, degree::Integer = 1)
    kf = Val(kernel); hb = float(bandwidth)
    out = Vector{Float64}(undef, length(xq))
    @inbounds for (k, x) in enumerate(xq)
        s0 = s1 = s2 = t0 = t1 = 0.0
        for j in eachindex(xn)
            w = _kernel(kf, (x - xn[j]) / hb); w == 0.0 && continue
            dx = xn[j] - x
            s0 += w; s1 += w*dx; s2 += w*dx^2; t0 += w*yn[j]; t1 += w*dx*yn[j]
        end
        if s0 == 0.0
            out[k] = yn[argmin(abs.(xn .- x))]              # outside support: nearest node
        elseif degree <= 0
            out[k] = t0 / s0                                # Nadaraya-Watson (local const)
        else
            det = s0*s2 - s1^2                              # local-linear intercept at x
            out[k] = det <= 1e-12 * s0 * s2 ? t0/s0 : (s2*t0 - s1*t1) / det
        end
    end
    return out
end

# --- Adapters for the fisher-rate pipeline -----------------------------------
# One BVP U(d) row -> a faithful (non-monotone) cubic-Hermite interpolant on the gap
# coordinate d.  Raw nodes are interpolated as-is: no isotonic pre-projection.
interp_U_full(d_grid, U_full) = CubicHermite(d_grid, collect(U_full))

# Drop-in replacement for the linear-response `U_ads` term in make_U_tot: builds
# the interpolant once for a given phi_inf's BVP row and returns U_att as a
# function of the FPT coordinate x (= surface gap d + sphere radius R).
function make_U_att_bvp(d_grid, U_full, R)
    itp = interp_U_full(d_grid, U_full)
    return x -> itp(x - R)
end

# One coarse BVP row (d_nodes, y_nodes) resampled onto the query gaps `dq`, per the
# PARAMS.smoothing choice: :none -> faithful CubicHermite (as before); :spline ->
# cubic smoothing spline (roughness penalty PARAMS.smoothing.lambda); a kernel name
# -> local kernel regression (PARAMS.smoothing.bandwidth / .degree).  Returns the
# resampled values -- the output format is the same regardless of method.
function smooth_row(d_nodes, y_nodes, dq)
    m = PARAMS.smoothing.method
    if m === :none
        itp = interp_U_full(d_nodes, y_nodes)
        return [itp(q) for q in dq]
    elseif m === :spline
        sp = SmoothSpline(d_nodes, y_nodes, PARAMS.smoothing.lambda)
        return [sp(q) for q in dq]
    else
        m in _KERNELS || error("unknown PARAMS.smoothing.method = :$m " *
                               "(expected :none, :spline, or one of $_KERNELS)")
        return kernel_smooth(d_nodes, y_nodes, dq;
                             kernel = m, bandwidth = PARAMS.smoothing.bandwidth,
                             degree = PARAMS.smoothing.degree)
    end
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

data_dir = PARAMS.sim.output_dir
R     = PARAMS.surface.R
h0    = PARAMS.surface.omega_b        # adsorption strength for the linear model
s, a  = PARAMS.rouse.s, PARAMS.rouse.a
Ngrid = PARAMS.fpt.Ngrid

linear_only = PARAMS.fpt.potential == :linear   # analytic run; no BVP files needed

# --- 1. gather the BVP sweep (skipped entirely for the analytic :linear model) ---
# :bvp aggregates output_cpu*.json into a dict keyed by cpu number.  :linear reads
# nothing: phi_inf comes from phi_inf_grid(PARAMS) and every BVP-derived quantity is
# a placeholder (d_grid empty; U_full/phi_surf empty rows -> NaN interp rows below).
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

    by_cpu = Dict{Int,NamedTuple}()
    for f in files
        data = JSON.parsefile(f)
        isempty(data) && continue                        # a CPU that got no shard
        cpu   = parse(Int, match(r"output_cpu(\d+)\.json$", f).captures[1])
        keyas = sort(collect(keys(data)); by = x -> parse(Int, x))
        phis, Us, Ps = Float64[], Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
        dgrid = Float64[]
        for ka in keyas
            recs = data[ka]
            kbs  = collect(keys(recs))
            # pick the omega_b branch this run's linear model is built at
            rec  = recs[kbs[argmin([abs(_f(recs[kb]["omega_b"]) - h0) for kb in kbs])]]
            push!(phis, _f(rec["phi_inf"]))
            push!(Us,   _farr(rec["U_full"]))
            push!(Ps,   _farr(rec["phi_surf"]))
            isempty(dgrid) && (dgrid = _farr(rec["d_grid"]))
        end
        by_cpu[cpu] = (phi_inf = phis, U_full = Us, phi_surf = Ps, d_grid = dgrid)
    end
    isempty(by_cpu) && error("all output_cpu*.json under $data_dir were empty")

    cpus          = sort(collect(keys(by_cpu)))
    phi_infs      = vcat((by_cpu[c].phi_inf  for c in cpus)...)
    U_full_rows   = vcat((by_cpu[c].U_full   for c in cpus)...)
    phi_surf_rows = vcat((by_cpu[c].phi_surf for c in cpus)...)
    d_grid        = by_cpu[cpus[1]].d_grid

    # --- 2. enforce the aggregated phi_inf grid matches PARAMS ---
    grid = phi_inf_grid(PARAMS)
    length(phi_infs) == length(grid) ||
        error("aggregated $(length(phi_infs)) phi_inf values; phi_inf_grid has $(length(grid))")
    all(isapprox.(phi_infs, grid; atol = 1e-9, rtol = 1e-6)) ||
        error("aggregated phi_inf != phi_inf_grid(PARAMS):\n  got $phi_infs\n  exp $grid")
end

# --- 3. common x-grid; Rouse end-to-end scale, as in the FPT return domain ---
x_max     = 4 * a * sqrt(s)
xg        = collect(range(R, x_max; length = Ngrid))
dq        = xg .- R                     # gap coordinate d = x - R (the BVP row axis)
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
        # Resample over the FINITE (d, U_full) nodes only, so a NaN at some d drops
        # that node instead of the whole phi_inf row.  smooth_row applies the chosen
        # PARAMS.smoothing method (faithful interp / smoothing spline / kernel).
        row = U_full_rows[i]
        fin = isfinite.(row)
        Uf  = count(fin) >= 2 ? smooth_row(d_grid[fin], row[fin], dq) : fill(NaN, Ngrid)

        # phi_surf onto the same x-grid so the :bvp local_conc can average it under the
        # nonlinear Boltzmann weight -- same smoothing method as U_full.
        srow = phi_surf_rows[i]
        sfin = isfinite.(srow)
        Ps   = count(sfin) >= 2 ? smooth_row(d_grid[sfin], srow[sfin], dq) : fill(NaN, Ngrid)
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
    "params"          => to_dict(PARAMS),
    "phi_inf"         => phi_infs,
    "omega_b"         => h0,
    "R"               => R,
    "d_grid"          => d_grid,
    "x_grid"          => xg,
    "U_full_raw"      => U_full_rows,        # per (phi_inf, d); NaN where unsolved
    "phi_surf_raw"    => phi_surf_rows,      # per (phi_inf, d)
    "U_full_interp"   => U_full_interp,      # per (phi_inf, x_grid)
    "phi_surf_interp" => phi_surf_interp,    # per (phi_inf, x_grid): BVP surface conc(x-R)
    "U_linear"        => U_linear,           # per (phi_inf, x_grid)
    "U_tot_full"      => U_tot_full,         # per (phi_inf, x_grid)
    "U_tot_linear"    => U_tot_linear,       # per (phi_inf, x_grid)
    "input_by_cpu"    => Dict(string(c) => Dict(
                            "phi_inf"  => by_cpu[c].phi_inf,
                            "U_full"   => by_cpu[c].U_full,
                            "phi_surf" => by_cpu[c].phi_surf,
                            "d_grid"   => by_cpu[c].d_grid) for c in cpus),
)
outpath = joinpath(data_dir, PARAMS.sim.interp_file)
open(outpath, "w") do io
    JSON.print(io, out)
end
println(linear_only ?
        "linear (analytic) model: no BVP files read; $(length(phi_infs)) phi_inf; " *
        "x-grid $(Ngrid) pts on [$(round(R; digits = 3)), $(round(x_max; digits = 3))]" :
        "aggregated $(length(cpus)) cpu file(s), $(length(phi_infs)) phi_inf; " *
        "x-grid $(Ngrid) pts on [$(round(R; digits = 3)), $(round(x_max; digits = 3))]")
println("saved -> $outpath")
