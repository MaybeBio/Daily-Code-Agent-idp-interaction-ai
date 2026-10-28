# General 1D backward Fokker–Planck (Smoluchowski) first-passage-time solver.
#
# Black-box core module: given an arbitrary potential as a vector on a uniform
# grid, a constant diffusion coefficient, and absorbing/reflecting endpoints,
# it returns the first-passage-time PDF (plus survival, MFPT, and the decay-rate
# spectrum whose slowest member is the Kramers escape rate).
#
# It is deliberately unaware of how the potential was built (Rouse + surface
# adsorption, etc.) — building/loading the potential is the caller's job.
#
# Method (Bicout–Szabo).  Discretise the overdamped dynamics on the grid as a
# nearest-neighbour Markov jump process whose hop rates
#     k[i->i+1] = (D/h^2) exp(-β (U[i+1]-U[i]) / 2)
#     k[i->i-1] = (D/h^2) exp(-β (U[i-1]-U[i]) / 2)
# satisfy detailed balance with π_i ∝ exp(-β U_i).  The generator is therefore
# similar to a symmetric matrix, so the decay-rate spectrum is real and the
# scheme stays positivity-preserving even where β|U'| h / D is large (deep,
# narrow wells) — the regime where a naive central-difference drift term fails.
# An absorbing end is a sink (its node is removed); a reflecting end zeroes the
# out-of-domain hop (zero flux).  Eigen-decomposing the interior (substochastic)
# generator gives the survival probability as a sum of exponentials, and the FPT
# density is its negative time-derivative.

"""
    FPTProblem(; U, dx, D=1.0, kBT=1.0, left_bc=:reflecting, right_bc=:absorbing,
                 p0=Float64[], t=Float64[])

Specification of a 1D first-passage-time problem.

- `U`      : potential sampled on a uniform grid, in energy units.
- `dx`     : grid spacing.
- `D`      : constant diffusion coefficient (sets the absolute time scale).
- `kBT`    : thermal energy; `β = 1/kBT`.
- `left_bc`, `right_bc` : `:reflecting` or `:absorbing` at the grid ends.
  At least one end must be `:absorbing`.
- `p0`     : initial distribution over the *full* grid (mass per node; it is
  renormalised over the non-absorbing nodes).  Empty ⇒ Boltzmann-in-domain.
- `t`      : output time grid for the PDF/survival.  Empty ⇒ auto (log-spaced).

Use [`solve_fpt`](@ref) to solve it.  A convenience constructor
`FPTProblem(x, U; ...)` accepts an explicit uniform grid vector `x` in place of `dx`.
"""
Base.@kwdef struct FPTProblem
    U::Vector{Float64}
    dx::Float64
    D::Float64          = 1.0
    kBT::Float64        = 1.0
    left_bc::Symbol     = :reflecting
    right_bc::Symbol    = :absorbing
    p0::Vector{Float64} = Float64[]
    t::Vector{Float64}  = Float64[]
end

# Explicit-grid convenience constructor: FPTProblem(x, U; D=..., right_bc=..., ...)
function FPTProblem(x::AbstractVector, U::AbstractVector; kwargs...)
    n = length(x)
    @assert n == length(U) "x and U must have equal length"
    @assert n >= 3 "need at least 3 grid points"
    dxs = diff(Float64.(collect(x)))
    dx = dxs[1]
    @assert all(isapprox.(dxs, dx; rtol=1e-8, atol=1e-12)) "grid x must be uniform"
    return FPTProblem(; U=Float64.(collect(U)), dx=dx, kwargs...)
end

function _validate(prob::FPTProblem)
    @assert length(prob.U) >= 3 "need at least 3 grid points"
    @assert prob.dx > 0 "dx must be positive"
    @assert prob.D > 0 "D must be positive"
    @assert prob.kBT > 0 "kBT must be positive"
    for bc in (prob.left_bc, prob.right_bc)
        @assert bc in (:reflecting, :absorbing) "bc must be :reflecting or :absorbing, got $bc"
    end
    @assert prob.left_bc == :absorbing || prob.right_bc == :absorbing "need at least one absorbing boundary"
    return nothing
end

# Indices of the non-absorbing (interior) nodes that carry probability.
function _interior_indices(prob::FPTProblem)
    N = length(prob.U)
    keep = trues(N)
    prob.left_bc  == :absorbing && (keep[1] = false)
    prob.right_bc == :absorbing && (keep[N] = false)
    return findall(keep)
end

# Interior forward generator Q̃ (dp/dt = Q̃ p) and the interior node indices.
# Q[j,i] is the rate i->j; the absorbing node's row/column is dropped so its
# neighbours keep their leak-out rate in the diagonal (substochastic).
function _build_generator(prob::FPTProblem)
    U = prob.U
    N = length(U)
    β = 1.0 / prob.kBT
    c = prob.D / prob.dx^2
    w_plus  = zeros(N)   # i -> i+1
    w_minus = zeros(N)   # i -> i-1
    @inbounds for i in 1:N-1
        w_plus[i] = c * exp(-β * (U[i+1] - U[i]) / 2)
    end
    @inbounds for i in 2:N
        w_minus[i] = c * exp(-β * (U[i-1] - U[i]) / 2)
    end
    # reflecting ends: no hop out of the domain (zero flux)
    prob.left_bc  == :reflecting && (w_minus[1] = 0.0)
    prob.right_bc == :reflecting && (w_plus[N]  = 0.0)

    Q = zeros(N, N)
    @inbounds for i in 1:N
        Q[i, i] = -(w_plus[i] + w_minus[i])
        i < N && (Q[i+1, i] = w_plus[i])
        i > 1 && (Q[i-1, i] = w_minus[i])
    end
    idx = _interior_indices(prob)
    return Q[idx, idx], idx
end

# Boltzmann (or user) initial distribution restricted & normalised on `idx`.
function _initial_distribution(prob::FPTProblem, idx)
    if isempty(prob.p0)
        Ui = prob.U[idx]
        w = exp.(-(Ui .- minimum(Ui)) ./ prob.kBT)
        return w ./ sum(w)
    else
        @assert length(prob.p0) == length(prob.U) "p0 must be defined on the full grid"
        p = prob.p0[idx]
        s = sum(p)
        @assert s > 0 "p0 places no mass on the interior (non-absorbing) nodes"
        return p ./ s
    end
end

"""
    boltzmann_ic(prob) -> Vector

Equilibrium (Boltzmann) initial distribution `∝ exp(-U/kBT)` on the full grid,
zero on the absorbing node(s), normalised over the interior.  This is the
default initial condition used by [`solve_fpt`](@ref).
"""
function boltzmann_ic(prob::FPTProblem)
    idx = _interior_indices(prob)
    p = zeros(length(prob.U))
    Ui = prob.U[idx]
    w = exp.(-(Ui .- minimum(Ui)) ./ prob.kBT)
    p[idx] .= w ./ sum(w)
    return p
end

_auto_tgrid(mfpt; n=500) = exp.(range(log(mfpt / 1e4), log(mfpt * 20.0); length=n))

"""
    solve_fpt(prob::FPTProblem)

Solve the backward Fokker–Planck problem and return a NamedTuple with

- `pdf`      : first-passage-time density evaluated on `t` (the headline output),
- `t`        : the time grid (supplied `prob.t`, else an auto log-spaced grid),
- `survival` : survival probability `S(t)`,
- `mfpt`     : mean first-passage time,
- `rates`    : decay-rate spectrum `λ_k` (ascending; `rates[1]` = Kramers rate),
- `weights`  : the modal amplitudes `a_k` with `S(t) = Σ a_k exp(-λ_k t)`.
"""
function solve_fpt(prob::FPTProblem)
    _validate(prob)
    Q̃, idx = _build_generator(prob)
    m = length(idx)
    @assert m >= 1 "no interior nodes to evolve"

    # symmetrise via the detailed-balance similarity transform H = D⁻¹ Q̃ D,
    # D = diag(√π), π_i ∝ exp(-β U_i)  (shift by min U for numerical safety).
    Ui = prob.U[idx]
    sq = exp.(-(Ui .- minimum(Ui)) ./ (2 * prob.kBT))   # √π (unnormalised)
    H = (1.0 ./ sq) .* Q̃ .* sq'
    H = Symmetric((H .+ H') ./ 2)

    E = eigen(H)
    λ = -E.values          # decay rates ≥ 0
    Φ = E.vectors          # orthonormal eigenvectors of H

    # right/left eigenvectors of Q̃: v_k = sq ⊙ Φ_k,  u_k = Φ_k ⊘ sq (biorthonormal)
    p0 = _initial_distribution(prob, idx)
    sum_v = vec(sum(sq .* Φ, dims=1))          # 1ᵀ v_k
    proj  = vec(sum((Φ ./ sq) .* p0, dims=1))  # u_kᵀ p0
    a = sum_v .* proj                          # S(t) = Σ a_k exp(-λ_k t)

    # sort by rate; drop only non-positive (roundoff) eigenvalues — a genuine
    # slow escape mode can sit many decades below the fast modes and must be kept
    o  = sortperm(λ)
    λs = λ[o]
    as = a[o]
    good = λs .> 0

    mfpt = sum(as[good] ./ λs[good])
    tgrid = isempty(prob.t) ? _auto_tgrid(mfpt) : prob.t
    S   = [sum(as .* exp.(-λs .* tt)) for tt in tgrid]
    pdf = [sum(as .* λs .* exp.(-λs .* tt)) for tt in tgrid]

    return (pdf=pdf, t=tgrid, survival=S, mfpt=mfpt, rates=λs, weights=as)
end

# cumulative trapezoid from the left (c[1]=0) and from the right (c[end]=0)
function _cumtrapz(y, dx)
    n = length(y); c = zeros(n)
    @inbounds for i in 2:n
        c[i] = c[i-1] + dx * (y[i] + y[i-1]) / 2
    end
    return c
end
function _cumtrapz_rev(y, dx)
    n = length(y); c = zeros(n)
    @inbounds for i in n-1:-1:1
        c[i] = c[i+1] + dx * (y[i] + y[i+1]) / 2
    end
    return c
end

"""
    fpt_mfpt_quadrature(prob::FPTProblem) -> Float64

Mean first-passage time from the exact 1D continuum double-integral
`τ(x₀) = (1/D) ∫ exp(βU(y)) ∫ exp(-βU(z)) dz dy`, averaged over the initial
distribution.  Independent of the eigen-solver; used to validate [`solve_fpt`](@ref).
Supports exactly one absorbing end with the other reflecting.
"""
function fpt_mfpt_quadrature(prob::FPTProblem)
    U = prob.U; N = length(U); β = 1.0 / prob.kBT; D = prob.D; dx = prob.dx
    u0 = minimum(U)
    gm = exp.(-β .* (U .- u0))   # exp(-βU) (shift cancels below)
    gp = exp.( β .* (U .- u0))   # exp(+βU)
    τ = zeros(N)
    if prob.left_bc == :reflecting && prob.right_bc == :absorbing
        G  = _cumtrapz(gm, dx)          # ∫ from node 1 (reflecting) to y
        Hc = _cumtrapz(gp .* G, dx)
        τ  = (Hc[N] .- Hc) ./ D         # ∫ from y to node N (absorbing)
    elseif prob.left_bc == :absorbing && prob.right_bc == :reflecting
        G  = _cumtrapz_rev(gm, dx)      # ∫ from y to node N (reflecting)
        Hc = _cumtrapz(gp .* G, dx)
        τ  = Hc ./ D                    # ∫ from node 1 (absorbing) to y
    else
        error("fpt_mfpt_quadrature supports exactly one absorbing end with the other reflecting")
    end
    idx = _interior_indices(prob)
    p = _initial_distribution(prob, idx)
    return sum(p .* τ[idx])
end
