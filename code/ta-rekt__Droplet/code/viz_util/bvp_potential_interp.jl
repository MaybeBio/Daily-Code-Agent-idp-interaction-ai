# Monotone, shape-preserving interpolation of a low-resolution BVP pair potential.
#
# The two-sphere BVP dashboard (dashboard_fisher_two_spheres.py) emits, for each
# bulk volume fraction phi_inf, the equilibrium pair potential U_full(d) sampled on
# a COARSE surface-to-surface gap grid d_grid (N_D ~ 30 nodes).  The FPT script
# (fisher_rate_from_fpt_solver.jl) needs the attractive potential as a smooth
# FUNCTION evaluated on a FINE grid (Ngrid ~ 1500) that also feeds finite
# differences into the backward-FPE generator -- so the resampled potential must be
# monotone and have a continuous, wiggle-free derivative or the drift term spikes.
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
#   PCHIP(x, y)                    -- monotone cubic Hermite; callable; deriv/deriv2
#   interp_U_full(d_grid, U_full)  -- cleaned PCHIP of one BVP U(d) row
#   make_U_att_bvp(d_grid, U, R)   -- drop-in U_ads replacement: x -> U_att(x-R)

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

function PCHIP(x::AbstractVector, y::AbstractVector)
    xs = Float64.(collect(x)); ys = Float64.(collect(y))
    @assert length(xs) == length(ys) "x and y must have equal length"
    @assert length(xs) >= 2 "need at least 2 nodes"
    @assert issorted(xs) && all(diff(xs) .> 0) "x must be strictly increasing"
    return PCHIP(xs, ys, _pchip_slopes(xs, ys))
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

# Value.  Constant (flat) extension outside the node range: below contact never
# happens (x >= R); above d_max the pair potential is already ~0, so holding the
# endpoint value keeps U_att bounded and finite for the FPT integration.
function (p::PCHIP)(xq::Real)
    xq <= p.x[1]   && return p.y[1]
    xq >= p.x[end] && return p.y[end]
    i = _seg(p, xq); h = p.x[i+1] - p.x[i]; t = (xq - p.x[i]) / h
    return p.y[i] * _h00(t) + h * p.m[i] * _h10(t) +
           p.y[i+1] * _h01(t) + h * p.m[i+1] * _h11(t)
end

function deriv(p::PCHIP, xq::Real)
    (xq <= p.x[1] || xq >= p.x[end]) && return 0.0
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
interp_U_full(d_grid, U_full) = PCHIP(d_grid, isotonic_nondecreasing(collect(U_full)))

# Drop-in replacement for the linear-response `U_ads` term in make_U_tot: builds
# the interpolant once for a given phi_inf's BVP row and returns U_att as a
# function of the FPT coordinate x (= surface gap d + sphere radius R).
function make_U_att_bvp(d_grid, U_full, R)
    itp = interp_U_full(d_grid, U_full)
    return x -> itp(x - R)
end
