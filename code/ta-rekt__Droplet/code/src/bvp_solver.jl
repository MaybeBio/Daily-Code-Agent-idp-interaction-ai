# Flory-Huggins bulk thermodynamics + 1D isolated-sphere BVP solver.
#
# Mirrors viz_util/bvp_solver.py: provides the Parameters dataclass and the
# scalar thermodynamic functions (fb, dfb_dphi, d2fb_dphi2, psi_linear,
# compute_mu_inf, compute_spinodal) that the 2D two-sphere code calls into.
#
# FREE-ENERGY CONVENTION: fb is the ASYMMETRIC (polymer-in-solvent) Flory-
# Huggins form -- chains of nb segments in a monomeric solvent, so only the
# solute entropy carries the 1/nb:
#
#     fb(phi) = (phi/nb)*log(phi) + (1-phi)*log(1-phi) + chi_b*phi*(1-phi)
#
# with critical point chi_c = (1 + 1/sqrt(nb))^2 / 2 and an asymmetric binodal
# (dilute branch far below 1/2, dense branch far above).  This differs from
# viz_util/bvp_solver.py, which divides the WHOLE entropy by nb and is thus a
# symmetric blend at effective chi = nb*chi_b (chi_c = 2/nb).  The two agree
# identically at nb = 1 and nowhere else, so a Julia/Python cross-check is only
# meaningful at nb = 1.  See the note in compute_spinodal.
#
# Only the minimum 1D solver needed by the two-sphere test suite is ported
# (the test's "single-sphere limit" comparison).  The full robust/shooting/
# continuation pipeline in the Python file is not reproduced.

is_spherical(p::Parameters) = isfinite(p.R)

function with(p::Parameters;
              nb=p.nb, chi_b=p.chi_b, omega_b=p.omega_b, kappa=p.kappa,
              phi_inf=p.phi_inf, L=p.L, R=p.R, eps=p.eps,
              conserve_volume=p.conserve_volume,
              parallel_plates=p.parallel_plates)
    Parameters(; nb=nb, chi_b=chi_b, omega_b=omega_b, kappa=kappa,
                 phi_inf=phi_inf, L=L, R=R, eps=eps,
                 conserve_volume=conserve_volume,
                 parallel_plates=parallel_plates)
end

@inline clip01(x::Real, eps_::Real=1e-10) = clamp(x, eps_, 1.0 - eps_)

# C1 smooth extension of log(x): log(x) for x >= eps, its first-order Taylor
# continuation below.  Value AND slope match at x = eps, so it is C1 and stays
# monotone; the second derivative simply drops to 0 in the tail.  dregularized_log
# / ddregularized_log are its exact 1st/2nd derivatives, so fb/dfb_dphi/d2fb_dphi2
# built from them below form ONE consistent free energy on all of R:
# d2fb == d(dfb)/dphi == d2(fb)/dphi2 everywhere, extension included.
#
# This is what fixes the clip-wall stall.  The old code clip01'd phi to
# [eps, 1-eps] INSIDE fb/dfb/d2fb, so past the boundary dfb froze flat (slope 0)
# while d2fb still returned ~1/eps: the Jacobian was no longer dF/dphi there, so
# Newton took wildly wrong steps and stalled wherever a near-saturated cell hit
# the clamp.  With the smooth extension an overshoot past (0,1) instead sees a
# finite, correctly-signed restoring gradient, and the paired fraction-to-boundary
# line search (_max_feasible_step) keeps iterates strictly inside so the extension
# is only ever a safety net.  On feasible phi in (eps, 1-eps) every value here is
# byte-identical to the old clamped form, so converged solutions -- and their free
# energies -- are unchanged.
@inline function regularized_log(x::Real, eps_::Real=1e-10)
    x > eps_ ? log(x) : log(eps_) + (x - eps_) / eps_
end
@inline dregularized_log(x::Real, eps_::Real=1e-10)  = x > eps_ ?  1.0 / x       : 1.0 / eps_
@inline ddregularized_log(x::Real, eps_::Real=1e-10) = x > eps_ ? -1.0 / (x * x) : 0.0

@inline function fb(phi::Real, p::Parameters)
    e  = p.eps
    L1 = regularized_log(phi, e)
    L2 = regularized_log(1 - phi, e)
    (1.0 / p.nb) * phi * L1 + (1 - phi) * L2 + p.chi_b * phi * (1 - phi)
end

@inline function dfb_dphi(phi::Real, p::Parameters)
    e   = p.eps
    L1  = regularized_log(phi, e);   dL1 = dregularized_log(phi, e)
    L2  = regularized_log(1 - phi, e); dL2 = dregularized_log(1 - phi, e)
    (1.0 / p.nb) * (L1 + phi * dL1) - (L2 + (1 - phi) * dL2) +
        p.chi_b * (1.0 - 2.0 * phi)
end

@inline function d2fb_dphi2(phi::Real, p::Parameters)
    e    = p.eps
    dL1  = dregularized_log(phi, e);    ddL1 = ddregularized_log(phi, e)
    dL2  = dregularized_log(1 - phi, e); ddL2 = ddregularized_log(1 - phi, e)
    (1.0 / p.nb) * (2.0 * dL1 + phi * ddL1) + (2.0 * dL2 + (1 - phi) * ddL2) -
        2.0 * p.chi_b
end

# Scalar psi_linear: returns (psi, dpsi, lam).
function psi_linear(z::Real, p::Parameters)
    alpha = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(alpha / p.kappa)
    psi  = -h0 * p.R * exp(-lam * z) / ((z + p.R) * p.kappa * (lam + 1 / p.R))
    dpsi = psi * (-lam - 1 / (z + p.R))
    return psi, dpsi, lam
end

# Broadcast-aware: applied to an array, returns (psi_array, dpsi_array, lam).
function psi_linear(z::AbstractArray, p::Parameters)
    alpha = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(alpha / p.kappa)
    psi  = @. -h0 * p.R * exp(-lam * z) / ((z + p.R) * p.kappa * (lam + 1 / p.R))
    dpsi = @. psi * (-lam - 1 / (z + p.R))
    return psi, dpsi, lam
end

compute_mu_inf(p::Parameters) = p.nb * dfb_dphi(p.phi_inf, p)

# Spinodal: the roots of d2fb_dphi2 = 0, i.e. 1/(nb*phi) + 1/(1-phi) = 2*chi_b.
# Clearing denominators (phi in (0,1), so the sign is preserved) gives
#     a*phi^2 + b*phi + 1 = 0,   a = 2*chi_b*nb,  b = (nb - 1) - a,
# which collapses to the symmetric 0.5*(1 -/+ sqrt(1 - 2/chi_b)) at nb = 1.
# A negative discriminant means d2fb_dphi2 > 0 everywhere: no spinodal, the
# single-phase state is stable at every composition.
function compute_spinodal(p::Parameters)
    p.chi_b <= 0 && return (nothing, nothing)
    a = 2.0 * p.chi_b * p.nb
    b = (p.nb - 1.0) - a
    discriminant = b * b - 4.0 * a
    discriminant < 0 && return (nothing, nothing)
    sqrt_d = sqrt(discriminant)
    return ((-b - sqrt_d) / (2 * a), (-b + sqrt_d) / (2 * a))
end

# ------------------------------------------------------------------
# 1D isolated-sphere BVP (for the single-sphere-limit comparison test)
# ------------------------------------------------------------------

struct IsolatedSphere1D
    z::Vector{Float64}
    phi::Vector{Float64}
    dphi::Vector{Float64}
    p::Parameters
    success::Bool
    # Branch bookkeeping (see solve_isolated_sphere).  F is the bc_reg-consistent
    # excess free energy of the returned profile.  When two stable solutions
    # coexist (bistable), branch says which one was returned (:thin / :thick),
    # F_thin / F_thick hold both energies and phi_other the profile NOT returned.
    # Otherwise branch = :unique, F_thin = F_thick = NaN and phi_other is empty.
    branch::Symbol
    bistable::Bool
    F::Float64
    F_thin::Float64
    F_thick::Float64
    phi_other::Vector{Float64}
end

@inline function _linterp(zg::AbstractVector, yg::AbstractVector,
                          zq::Real, dz::Real, z0::Real)
    n = length(zg)
    t = (zq - z0) / dz
    i = clamp(Int(floor(t)), 0, n - 2)
    a = t - i
    yg[i + 1] * (1 - a) + yg[i + 2] * a
end

function (s::IsolatedSphere1D)(zq::Real)
    z = s.z
    dz = z[2] - z[1]
    zc = clamp(zq, z[1], z[end])
    _linterp(z, s.phi, zc, dz, z[1])
end

function (s::IsolatedSphere1D)(zq::AbstractVector)
    [s(z) for z in zq]
end

# Robin surface slope phi'(0).  The plain BC is -kappa*phi'(0) + omega_b = 0;
# with bc_reg > 0 a regularization term is added (mirrors make_bc in
# viz_util/bvp_solver.py, line 360):
#     -kappa*phi'(0) + omega_b + bc_reg * log(phi_s / (1 - phi_s)) = 0.
# As phi_s -> 1 the log term diverges and opposes the attractive omega_b,
# pinning a unique *finite* adsorption branch near a surface phase transition
# (where the bare BC is bistable: a thin film and a runaway phi_s -> 1).  It
# removes the runaway, not prewetting bistability: a finite thick film can
# still coexist with the thin one at bc_reg > 0 (e.g. R = 50, chi_b = 2.5,
# phi_inf = 0.14, bc_reg = 0.05; test_bvp_solver.jl "EnvelopeIdentity").
@inline function _surf_slope(phi_s::Real, omega_b::Real, p::Parameters, bc_reg::Real)
    bc_reg == 0.0 && return omega_b / p.kappa
    ps = clip01(phi_s, p.eps)
    return (omega_b + bc_reg * log(ps / (1.0 - ps))) / p.kappa
end

# d/dphi_s of _surf_slope (for the analytic Jacobian).
@inline function _surf_slope_deriv(phi_s::Real, p::Parameters, bc_reg::Real)
    bc_reg == 0.0 && return 0.0
    ps = clip01(phi_s, p.eps)
    (p.eps < phi_s < 1.0 - p.eps) || return 0.0   # in the clamp the slope is flat
    return bc_reg / (p.kappa * ps * (1.0 - ps))
end

# Morse index of a solution: the number of unstable directions of the free
# energy there, read off the tridiagonal Jacobian J (sub dl[2:n], diag dd, super
# du[1:n-1], as filled by jacobian! in solve_isolated_sphere).  0 = local minimum
# (thin or thick film), 1 = the saddle (middle branch) between them.
#
# The residual is -dF/dphi per unit volume, so J is minus the Hessian with row i
# divided by its volume weight -- w_1 = R(R+dz)/2 for the ghost row, (R+z_i)^2
# otherwise; with r_{i+1} = r_i + dz, r_i^2 (1/dz^2 + 1/(r_i dz)) = r_{i+1}^2
# (1/dz^2 - 1/(r_{i+1} dz)), so W*J is exactly symmetric.  By Sylvester's law of
# inertia the negative pivots of the UNPIVOTED LDL' of -W*J count its negative
# eigenvalues, and the positive row weights only rescale pivots, never flip them:
# so count the negative pivots of -J itself.  Do not use the factors behind `\`
# for this -- LAPACK's tridiagonal LU pivots rows, and those pivots carry no
# inertia.
function _morse_index(dl::AbstractVector, dd::AbstractVector, du::AbstractVector)
    piv = -dd[1]
    n   = piv < 0 ? 1 : 0
    @inbounds for i in 2:length(dd)
        piv = -dd[i] - dl[i] * du[i - 1] / piv
        piv < 0 && (n += 1)
    end
    return n
end

# Interior composition of a thick film: the dense root of fb'(phi) = fb'(phi_inf),
# i.e. the metastable liquid at the reservoir chemical potential -- the "well" of
# Delta_f(phi) = fb(phi) - fb(phi_inf) - fb'(phi_inf)(phi - phi_inf) in the planar
# graphical construction.  fb' rises to a max at the lower spinodal, falls to a
# min at the upper one and rises to +infinity as phi -> 1, so for a dilute
# reservoir (phi_inf below the lower spinodal) the root exists iff fb'(upper
# spinodal) < fb'(phi_inf), and is then unique above the upper spinodal.  When it
# does not exist, Delta_f has no turning points at all (the barrier needs the same
# condition), the construction has exactly one root, and there is no thick branch
# to find.  Returns nothing in that case.
function _dense_well(p::Parameters)
    lo_sp, hi_sp = compute_spinodal(p)
    lo_sp === nothing && return nothing
    target = dfb_dphi(p.phi_inf, p)
    (p.phi_inf < lo_sp && dfb_dphi(hi_sp, p) < target) || return nothing
    a, b = hi_sp, 1.0 - p.eps                 # fb' - target: < 0 at a, > 0 at b
    for _ in 1:80
        m = 0.5 * (a + b)
        dfb_dphi(m, p) < target ? (a = m) : (b = m)
    end
    return 0.5 * (a + b)
end

# BRANCH TRACKING (track_branches = true).  Past the prewetting transition the
# BVP has three solutions -- thin film, saddle, thick film -- and Newton returns
# whichever one its initial guess falls toward.  With tracking on, two
# candidates are built from separate anchors and the lower free energy wins:
#
#   thin  : the cold start below; on failure, OR if it is a saddle, OR if it
#           coincides with the thick candidate, the omega_b continuation from
#           the uniform omega_b = 0 state (unambiguously thin), which follows the
#           thin branch by continuity until it either arrives or falls off its
#           fold.  (The cold start is not reliably thin: when |omega_b|*ell/kappa
#           exceeds 1 - phi_inf, clip01 turns it into a saturated layer.)
#   thick : a film of the well composition (_dense_well) out to `reach` with a
#           tanh edge, phi_inf beyond -- the 1D analogue of _dense_seed in
#           spec/fisher_two_spheres.jl.  It is only a starting guess: Newton finds
#           the film's actual composition, thickness and phi_s.  Film thickness
#           is the soft mode and the saddle is the critical film, so a seed
#           thicker than it relaxes onto the thick branch: if this lands on a
#           saddle or on the thin solution, it is retried at doubling reach until
#           it finds a distinct minimum or exceeds L_iso/2.  The ladder is what
#           makes the seed scale-free: the working reach is ~2 at kappa = 1,
#           nb = 1 but 5-20 at kappa = 10, nb = 2230 (interface width ~25), and
#           no fixed multiple of any one length covers both.  With
#           no well there is no thick candidate (exact in the planar limit; at
#           finite R curvature only shrinks the window), so the solve is then the
#           single-candidate one.  phi_dense overrides the well composition as the
#           seed plateau (and turns tracking on by itself).
#
# Each converged candidate is classified by its Morse index (_morse_index;
# saddles are never returned while a minimum exists), duplicates are merged by
# profile distance, and two distinct minima are labelled by adsorption Gamma =
# Int 4 pi (R+z)^2 (phi - phi_inf) dz (larger = thick), never by which seed
# produced them.  F is excess_free_energy with the solve's bc_reg.  Missing a
# branch near its own fold is benign for the selection: there it is the
# metastable one.  track_branches = false (the default) skips all of this: one
# candidate, as before.  Cost: ~ms when both films exist; when the well exists
# but no thick film does (small R), Newton exhausts maxiter on the seed twice,
# ~0.7 s at N = 1500.
function solve_isolated_sphere(p::Parameters;
                               L_iso::Real = 10.0,
                               N::Int      = 400,
                               f_tol::Real = 1e-8,
                               maxiter::Int = 200,
                               bc_reg::Real = 0.0,
                               omega_continuation::Bool = true,
                               track_branches::Bool = false,
                               phi_dense::Union{Nothing,Real} = nothing,
                               reach::Real = 2.0)
    @assert isfinite(p.R) "solve_isolated_sphere requires finite p.R"
    z  = collect(range(0.0, L_iso; length = N))
    dz = z[2] - z[1]
    mu_inf  = compute_mu_inf(p)
    Nu = N - 1                                # unknowns at z[1:N-1]
    phi_far = p.phi_inf                       # Dirichlet at z[N]

    # Residual and analytic tridiagonal Jacobian at effective surface field ob.
    function residual!(F, u, ob)
        @inbounds for i in 1:Nu
            phi_i = u[i]
            if i == 1
                # Robin BC ghost: phi(-dz) = phi(+dz) - 2*dz*sk.
                sk    = _surf_slope(phi_i, ob, p, bc_reg)
                ddphi = (2 * u[2] - 2 * phi_i - 2 * dz * sk) / dz^2
                lap   = ddphi + (2.0 / p.R) * sk
            else
                phi_l = u[i - 1]
                phi_r = i == Nu ? phi_far : u[i + 1]
                ddphi = (phi_r - 2 * phi_i + phi_l) / dz^2
                dphi  = (phi_r - phi_l) / (2 * dz)
                rr    = p.R + z[i]
                lap   = ddphi + (2.0 / rr) * dphi
            end
            F[i] = p.kappa * lap - (dfb_dphi(phi_i, p) - mu_inf / p.nb)
        end
        return nothing
    end
    function jacobian!(dl, dd, du, u, ob)
        @inbounds for i in 1:Nu
            if i == 1
                sk1   = _surf_slope_deriv(u[1], p, bc_reg)
                dd[i] = p.kappa * (-2.0 / dz^2) - d2fb_dphi2(u[i], p) +
                        p.kappa * sk1 * (-2.0 / dz + 2.0 / p.R)
                du[i] = p.kappa * (2.0 / dz^2)
                dl[i] = 0.0
            else
                rr    = p.R + z[i]
                dd[i] = p.kappa * (-2.0 / dz^2) - d2fb_dphi2(u[i], p)
                dl[i] = p.kappa * (1.0 / dz^2 - (1.0 / rr) / dz)
                du[i] = i == Nu ? 0.0 : p.kappa * (1.0 / dz^2 + (1.0 / rr) / dz)
            end
        end
        return nothing
    end

    # Damped (backtracking line search) Newton with the analytic tridiagonal
    # Jacobian.  Returns (u, |F|).  Robust on the stiff prewetting branch where
    # the bare dense Newton stalls.
    function newton(u0, ob)
        # `local`: u, F and nrm are also locals of solve_isolated_sphere, and an
        # assignment in a nested function would otherwise rebind THOSE -- every
        # call would silently overwrite the caller's current candidate.
        local u, F, nrm
        u  = copy(u0)
        F  = Vector{Float64}(undef, Nu)
        dl = Vector{Float64}(undef, Nu)
        dd = Vector{Float64}(undef, Nu)
        du = Vector{Float64}(undef, Nu)
        ut = Vector{Float64}(undef, Nu)
        residual!(F, u, ob); nrm = maximum(abs, F)
        for _ in 1:maxiter
            nrm < f_tol && break
            jacobian!(dl, dd, du, u, ob)
            δ = Tridiagonal(dl[2:Nu], dd[1:Nu], du[1:Nu - 1]) \ F
            # Backtrack from the fraction-to-boundary cap, as _solve_newton_sparse
            # does: it is 1 whenever the full step stays inside (eps, 1-eps), and
            # then every trial is feasible and the clamp below a no-op, so such
            # solves are unchanged.  Starting at lambda = 1 and clamping instead
            # lets one overshooting step pin the wall cell at 1 - eps (fb'' ~
            # 1/eps there), where the line search stalls for good -- the failure
            # mode of the dense (thick-film) seed.
            λ = _max_feasible_step(u, δ, p.eps, 1.0 - p.eps); accepted = false
            for _ in 1:40
                @inbounds @. ut = clamp(u - λ * δ, p.eps, 1.0 - p.eps)
                residual!(F, ut, ob); newnrm = maximum(abs, F)
                if newnrm < nrm
                    @inbounds @. u = ut; nrm = newnrm; accepted = true; break
                end
                λ *= 0.5
            end
            accepted || break
        end
        return u, nrm
    end

    # Initial guess: exponential decay from a perturbed phi_inf.
    d2  = max(d2fb_dphi2(p.phi_inf, p), 0.05)
    ell = sqrt(p.kappa / d2)
    u0  = Vector{Float64}(undef, Nu)
    @inbounds for i in 1:Nu
        u0[i] = clip01(p.phi_inf - p.omega_b * ell / p.kappa * exp(-z[i] / ell), p.eps)
    end

    jl = Vector{Float64}(undef, Nu); jd = similar(jl); ju = similar(jl)
    stable(v, r) = r < f_tol &&
                   (jacobian!(jl, jd, ju, v, p.omega_b); _morse_index(jl, jd, ju) == 0)
    # Distinct branches differ by O(0.1) near the wall; two converged copies of
    # one solution by ~ |J^-1| * f_tol, well under this for f_tol <= 1e-5.
    same(a, b) = maximum(abs(a[i] - b[i]) for i in 1:Nu) < 1e-3

    # Thick-film seed plateau: the caller's phi_dense, else the well composition
    # (nothing when there is no well -- and then no thick candidate at all).
    phi_film  = phi_dense !== nothing ? phi_dense :
                track_branches ? _dense_well(p) : nothing
    has_thick = phi_film !== nothing
    if has_thick
        # A tanh interface of the dense-side width, not a sharp step: the step's
        # interface residual dominates the max-norm line search, and even with the
        # fraction-to-boundary cap a step seed misses the thick branch at up to
        # ~1 in 5 points of the chi_b = 2.8 bistable window (some on the side
        # where thick is the stable answer).  The tanh seed finds it at all of
        # them (test_bvp_solver.jl, "PlanarWindowCoverage").
        xi_d = sqrt(p.kappa / max(d2fb_dphi2(phi_film, p), 0.05))
        dense_seed = r -> [clip01(p.phi_inf + (phi_film - p.phi_inf) *
                                  0.5 * (1.0 - tanh((z[i] - r) / (2.0 * xi_d))), p.eps)
                           for i in 1:Nu]
        u_b, nrm_b = newton(dense_seed(reach), p.omega_b)
    end

    u, nrm = newton(u0, p.omega_b)
    # omega_b continuation from the uniform (omega_b = 0) state, tracing the
    # branch continuously when the cold start fails (same strategy as
    # solve_two_spheres, and the branch the 2D solver follows) -- or, when
    # tracking branches, when the cold start did not land on a distinct thin
    # minimum.
    redo = nrm >= f_tol ||
           (has_thick && (!stable(u, nrm) || (nrm_b < f_tol && same(u, u_b))))
    if redo && omega_continuation && p.omega_b != 0.0
        u_cur = fill(clip01(p.phi_inf, p.eps), Nu)
        nrm_cur = Inf
        n_steps = 8
        for k in 1:n_steps
            u_cur, nrm_cur = newton(u_cur, (k / n_steps) * p.omega_b)
            nrm_cur >= f_tol && break
        end
        nrm_cur < f_tol && ((u, nrm) = (u_cur, nrm_cur))
    end
    r_seed = reach
    while has_thick && (!stable(u_b, nrm_b) || (nrm < f_tol && same(u, u_b))) &&
          2 * r_seed <= L_iso / 2
        r_seed *= 2
        u_b, nrm_b = newton(dense_seed(r_seed), p.omega_b)
    end

    stable_a = stable(u, nrm)
    stable_b = has_thick && stable(u_b, nrm_b)
    stable_b && !stable_a && ((u, nrm) = (u_b, nrm_b))   # the only minimum found
    bistable = stable_a && stable_b && !same(u, u_b)
    success  = nrm < f_tol

    function profile(v)
        ph      = vcat(v, phi_far)
        dph     = similar(ph)
        dph[1]  = _surf_slope(ph[1], p.omega_b, p, bc_reg)
        @inbounds for i in 2:N - 1
            dph[i] = (ph[i + 1] - ph[i - 1]) / (2 * dz)
        end
        dph[N] = (ph[N] - ph[N - 1]) / dz
        return ph, dph
    end
    phi, dphi = profile(u)
    F = excess_free_energy(phi, dphi, z, p; bc_reg = bc_reg)

    branch, F_thin, F_thick, phi_other = :unique, NaN, NaN, Float64[]
    if bistable
        phi_b, dphi_b = profile(u_b)
        F_b    = excess_free_energy(phi_b, dphi_b, z, p; bc_reg = bc_reg)
        gamma(ph) = _trapz([4.0 * pi * (p.R + z[i])^2 * (ph[i] - p.phi_inf) for i in 1:N], z)
        a_thick = gamma(phi) > gamma(phi_b)
        F_thin, F_thick = a_thick ? (F_b, F) : (F, F_b)
        if F_b < F
            phi, dphi, phi_b, F = phi_b, dphi_b, phi, F_b
            a_thick = !a_thick
        end
        branch    = a_thick ? :thick : :thin
        phi_other = phi_b
    end

    return IsolatedSphere1D(z, phi, dphi, p, success,
                            branch, bistable, F, F_thin, F_thick, phi_other)
end

# ------------------------------------------------------------------
# 1D planar slab BVP: two parallel plates, or one plate + reservoir
# ------------------------------------------------------------------
# Mirrors the planar branches of viz_util/bvp_solver.py (make_ode / make_bc with
# R = inf) and the separation sweep of viz_util/bvp_parallel_plates.py:
#
#     kappa * phi'' = dfb/dphi - mu_inf / nb          on z in [0, L]
#
# with the Robin adsorption BC at the plate at z = 0,
#     -kappa*phi'(0) + omega_b = 0,
# and at z = L either
#     plates = true : +kappa*phi'(L) + omega_b = 0   (second, identical plate)
#     plates = false:  phi(L) = phi_inf              (far-field reservoir)
#
# The plate separation is p.L (override with the L_gap keyword), and bc_reg
# regularizes the surface BC exactly as in solve_isolated_sphere, on BOTH plates.
#
# GEOMETRY NOTE: `parallel_plates` overrides `R`.  The Python Parameters set
# R = inf for the plate problem; here R stays at the physical sphere radius (it
# is the Derjaguin radius consumed by plates_pair_potential), so the planar
# branch keys off p.parallel_plates, never off isfinite(R).  conserve_volume is
# ignored, as it is by solve_isolated_sphere.

struct PlanarSlab1D
    z::Vector{Float64}
    phi::Vector{Float64}
    dphi::Vector{Float64}
    p::Parameters
    plates::Bool
    success::Bool
end

function (s::PlanarSlab1D)(zq::Real)
    z = s.z
    dz = z[2] - z[1]
    zc = clamp(zq, z[1], z[end])
    _linterp(z, s.phi, zc, dz, z[1])
end

function (s::PlanarSlab1D)(zq::AbstractVector)
    [s(z) for z in zq]
end

# Linear-response profile psi = phi - phi_inf in planar geometry (the R -> inf
# counterpart of psi_linear).  Linearizing about phi_inf gives kappa*psi'' =
# alpha*psi with alpha = d2fb/dphi2(phi_inf), lam = sqrt(alpha/kappa), so
#     plates = true :  psi(z) = -h0 * cosh(lam*(z - L/2)) / (kappa*lam*sinh(lam*L/2))
#     plates = false:  psi(z) = -h0 * exp(-lam*z) / (kappa*lam)
# (the plate solution tends to the single-plate one as lam*L -> infinity).  The
# cosh/sinh ratio is evaluated in exponential form so that large lam*L, which
# occurs at small phi_inf where alpha ~ 1/phi_inf, cannot overflow.
# Accepts a scalar or an array z; returns (psi, dpsi, lam) like psi_linear.
function psi_linear_plates(z, p::Parameters;
                           L_gap::Real  = p.L,
                           plates::Bool = p.parallel_plates)
    alpha = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(alpha / p.kappa)
    if !plates
        psi  = @. -h0 * exp(-lam * z) / (p.kappa * lam)
        dpsi = @. h0 * exp(-lam * z) / p.kappa
        return psi, dpsi, lam
    end
    b   = lam * L_gap / 2
    den = -expm1(-2 * b)                       # 1 - exp(-2b), i.e. 2*sinh(b)*exp(-b)
    ep  = @. exp(lam * (z - L_gap / 2) - b)    # exp(t - b),  t = lam*(z - L/2)
    em  = @. exp(-lam * (z - L_gap / 2) - b)   # exp(-t - b); both exponents <= 0
    psi  = @. -h0 * (ep + em) / (p.kappa * lam * den)
    dpsi = @. -h0 * (ep - em) / (p.kappa * den)
    return psi, dpsi, lam
end

function solve_parallel_plates(p::Parameters;
                               plates::Bool = p.parallel_plates,
                               L_gap::Real  = p.L,
                               N::Int       = 400,
                               f_tol::Real  = 1e-8,
                               maxiter::Int = 200,
                               bc_reg::Real = 0.0,
                               omega_continuation::Bool = true,
                               phi_seed::Union{Nothing,Real,AbstractVector} = nothing)
    z  = collect(range(0.0, L_gap; length = N))
    dz = z[2] - z[1]
    mu_inf  = compute_mu_inf(p)
    Nu = plates ? N : N - 1                   # plates: both ends Robin; else Dirichlet at z[N]
    phi_far = p.phi_inf

    # Residual and analytic tridiagonal Jacobian at effective surface field ob.
    # Both Robin ends use the same ghost construction as solve_isolated_sphere,
    # with phi'(0) = +sk and phi'(L) = -sk (the plates are mirror images); the
    # planar Laplacian drops the 2/r curvature term.
    function residual!(F, u, ob)
        @inbounds for i in 1:Nu
            phi_i = u[i]
            if i == 1
                sk    = _surf_slope(phi_i, ob, p, bc_reg)
                ddphi = (2 * u[2] - 2 * phi_i - 2 * dz * sk) / dz^2
            elseif plates && i == Nu
                sk    = _surf_slope(phi_i, ob, p, bc_reg)
                ddphi = (2 * u[Nu - 1] - 2 * phi_i - 2 * dz * sk) / dz^2
            else
                phi_l = u[i - 1]
                phi_r = i == Nu ? phi_far : u[i + 1]
                ddphi = (phi_r - 2 * phi_i + phi_l) / dz^2
            end
            F[i] = p.kappa * ddphi - (dfb_dphi(phi_i, p) - mu_inf / p.nb)
        end
        return nothing
    end
    function jacobian!(dl, dd, du, u, ob)
        @inbounds for i in 1:Nu
            if i == 1
                sk1   = _surf_slope_deriv(u[1], p, bc_reg)
                dd[i] = p.kappa * (-2.0 / dz^2 - 2.0 * sk1 / dz) - d2fb_dphi2(u[i], p)
                du[i] = p.kappa * (2.0 / dz^2)
                dl[i] = 0.0
            elseif plates && i == Nu
                skN   = _surf_slope_deriv(u[Nu], p, bc_reg)
                dd[i] = p.kappa * (-2.0 / dz^2 - 2.0 * skN / dz) - d2fb_dphi2(u[i], p)
                dl[i] = p.kappa * (2.0 / dz^2)
                du[i] = 0.0
            else
                dd[i] = p.kappa * (-2.0 / dz^2) - d2fb_dphi2(u[i], p)
                dl[i] = p.kappa / dz^2
                du[i] = i == Nu ? 0.0 : p.kappa / dz^2
            end
        end
        return nothing
    end

    function newton(u0, ob)
        u  = copy(u0)
        F  = Vector{Float64}(undef, Nu)
        dl = Vector{Float64}(undef, Nu)
        dd = Vector{Float64}(undef, Nu)
        du = Vector{Float64}(undef, Nu)
        ut = Vector{Float64}(undef, Nu)
        residual!(F, u, ob); nrm = maximum(abs, F)
        for _ in 1:maxiter
            nrm < f_tol && break
            jacobian!(dl, dd, du, u, ob)
            δ = Tridiagonal(dl[2:Nu], dd[1:Nu], du[1:Nu - 1]) \ F
            λ = 1.0; accepted = false
            for _ in 1:40
                @inbounds @. ut = clamp(u - λ * δ, p.eps, 1.0 - p.eps)
                residual!(F, ut, ob); newnrm = maximum(abs, F)
                if newnrm < nrm
                    @inbounds @. u = ut; nrm = newnrm; accepted = true; break
                end
                λ *= 0.5
            end
            accepted || break
        end
        return u, nrm
    end

    # Initial guess: phi_inf plus the linearized planar adsorption profile, i.e.
    # started from the dilute (thin-film) side.  At small separations the gap is
    # bistable -- a thin adsorbed film versus a capillary-condensed dense gap, the
    # plate analogue of the unbridged/bridged branches of the two-sphere solver --
    # and which of the two Newton settles on from this guess can then depend on
    # the grid.  To pin a branch, pass phi_seed: a uniform value (e.g. the dense
    # coexistence composition, the way spec/fisher_two_spheres.jl seeds bridging)
    # or a full profile warm-started from a neighbouring separation.
    u0 = Vector{Float64}(undef, Nu)
    if phi_seed === nothing
        # Superposition of one exponential film per plate: psi_linear_plates up to
        # its 1/(1 - exp(-lam*L)) normalization.  The decay length is floored the
        # way solve_isolated_sphere floors it, so a phi_inf inside the spinodal
        # (d2fb_dphi2 < 0, i.e. imaginary lam) still yields a real guess.
        ell = sqrt(p.kappa / max(d2fb_dphi2(p.phi_inf, p), 0.05))
        amp = -p.omega_b * ell / p.kappa
        @inbounds for i in 1:Nu
            psi0  = amp * (exp(-z[i] / ell) +
                           (plates ? exp(-(L_gap - z[i]) / ell) : 0.0))
            u0[i] = clip01(p.phi_inf + psi0, p.eps)
        end
    elseif phi_seed isa Real
        fill!(u0, clip01(phi_seed, p.eps))
    else
        @assert length(phi_seed) >= Nu "phi_seed must hold at least $Nu values"
        @inbounds for i in 1:Nu
            u0[i] = clip01(phi_seed[i], p.eps)
        end
    end

    u, nrm = newton(u0, p.omega_b)
    # omega_b continuation from the uniform (omega_b = 0) state, as in
    # solve_isolated_sphere: trace the branch when the cold start fails.  Skipped
    # under an explicit phi_seed -- that fallback lands on the dilute branch, so
    # running it would silently answer a branch the caller did not ask for.
    if nrm >= f_tol && omega_continuation && phi_seed === nothing && p.omega_b != 0.0
        u_cur = fill(clip01(p.phi_inf, p.eps), Nu)
        n_steps = 8
        for k in 1:n_steps
            u_cur, nrm = newton(u_cur, (k / n_steps) * p.omega_b)
            nrm >= f_tol && break
        end
        nrm < f_tol && (u = u_cur)
    end
    success = nrm < f_tol

    phi      = plates ? u : vcat(u, phi_far)
    dphi     = similar(phi)
    dphi[1]  = _surf_slope(phi[1], p.omega_b, p, bc_reg)
    @inbounds for i in 2:N - 1
        dphi[i] = (phi[i + 1] - phi[i - 1]) / (2 * dz)
    end
    dphi[N] = plates ? -_surf_slope(phi[N], p.omega_b, p, bc_reg) :
                       (phi[N] - phi[N - 1]) / dz

    return PlanarSlab1D(z, phi, dphi, p, plates, success)
end

# ------------------------------------------------------------------
# Excess free energy and the Derjaguin plate-plate pair potential
# ------------------------------------------------------------------

_trapz(y::AbstractVector, x::AbstractVector) =
    sum(0.5 * (y[k] + y[k + 1]) * (x[k + 1] - x[k]) for k in 1:length(x) - 1)

# Grand-potential excess over the bulk reservoir at phi_inf; the Julia port of
# excess_free_energy in viz_util/bvp_solver.py:
#
#     DeltaOmega = int_V [ df(phi) + (kappa/2)*phi'^2 ] dV
#                  + omega_b * int_dV (phi - phi_inf) dA
#     df(phi)    = fb(phi) - fb(phi_inf) - mu_ref/nb * (phi - phi_inf)
#
# With mu_ref = nb * dfb/dphi(phi_inf) (the default) df vanishes to linear order
# at phi_inf, so the integral converges as L -> infinity.  Geometry, dispatched
# on p (parallel_plates first, per the GEOMETRY NOTE above):
#   parallel plates : per unit area, one surface term per plate,
#   spherical       : 4*pi*(R+z)^2 volume weight, 4*pi*R^2 surface term,
#   single plate    : per unit area, one surface term at z = 0.
#
# SURFACE TERM AND bc_reg.  The Robin BC the solvers impose, kappa*phi'(0) =
# omega_b + bc_reg*log(phi_s/(1-phi_s)), is the natural BC of the surface free
# energy
#
#     fs(phi) = omega_b*phi + bc_reg*[phi log phi + (1-phi) log(1-phi)],
#
# charged here per surface as fs(phi_s) - fs(phi_inf) (referenced to a surface
# bathed in bulk, as in excess_free_energy_two_spheres).  Pass the bc_reg the
# profile was solved with: the entropy piece is branch-dependent (a thin and a
# thick film sit at very different phi_s), so dropping it at bc_reg > 0 can
# invert the ordering of coexisting branches.  bc_reg = 0 recovers the bare
# omega_b*(phi_s - phi_inf).
function excess_free_energy(phi::AbstractVector, dphi::AbstractVector,
                            z::AbstractVector, p::Parameters;
                            mu_ref::Union{Nothing,Real} = nothing,
                            bc_reg::Real = 0.0)
    mu    = mu_ref === nothing ? p.nb * dfb_dphi(p.phi_inf, p) : mu_ref
    f_inf = fb(p.phi_inf, p)
    n     = length(z)
    g     = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        df   = fb(phi[k], p) - f_inf - (mu / p.nb) * (phi[k] - p.phi_inf)
        g[k] = df + 0.5 * p.kappa * dphi[k]^2
    end
    hs(v) = (c = clip01(v, p.eps); c * log(c) + (1.0 - c) * log(1.0 - c))
    fs(v) = p.omega_b * (v - p.phi_inf) + bc_reg * (hs(v) - hs(p.phi_inf))
    if p.parallel_plates
        return _trapz(g, z) + fs(phi[1]) + fs(phi[n])
    elseif is_spherical(p)
        @inbounds for k in 1:n
            g[k] *= 4.0 * pi * (p.R + z[k])^2
        end
        return _trapz(g, z) + 4.0 * pi * p.R^2 * fs(phi[1])
    else
        return _trapz(g, z) + fs(phi[1])
    end
end

# Energy of a solved slab.  The geometry flags come from the solution itself: the
# single-plate reference must not be charged a second surface term (plates =
# false), and a slab is planar whatever p.R says, so R is forced to Inf -- with a
# finite R the 4-arg method would otherwise take its spherical branch.
excess_free_energy(s::PlanarSlab1D) =
    excess_free_energy(s.phi, s.dphi, s.z,
                       with(s.p; parallel_plates = s.plates, R = Inf, L = s.z[end]))

# Derjaguin pair potential between two spheres of radius R_derj from the plate
# geometry, mirroring the sweep in viz_util/bvp_parallel_plates.py:
#
#     U(d) = pi * R_derj * ( W_plates(d) - 2*W_single )
#
# where W is the excess free energy per unit area.  Unlike the Python script,
# which re-solves the single plate on the shrinking domain [0, d] at every step,
# W_single is computed ONCE on a domain of length L_ref (default: the largest
# separation) -- the isolated-plate reference is a constant, and truncating it at
# small d forces phi(d) = phi_inf and corrupts the subtraction.  Pick L_ref well
# beyond the decay length sqrt(kappa/d2fb_dphi2(phi_inf)).
# Returns (U, success, W_plates, W_single); solver keywords are forwarded.
function plates_pair_potential(p::Parameters, distances::AbstractVector;
                               R_derj::Real = p.R,
                               L_ref::Real  = maximum(distances),
                               N::Int       = 400,
                               kwargs...)
    ref      = solve_parallel_plates(with(p; L = L_ref); plates = false, N = N, kwargs...)
    W_single = excess_free_energy(ref)
    ref.success || @warn "single-plate reference did not converge" L_ref

    U       = Vector{Float64}(undef, length(distances))
    W       = Vector{Float64}(undef, length(distances))
    success = Vector{Bool}(undef, length(distances))
    for (i, d) in enumerate(distances)
        sol        = solve_parallel_plates(with(p; L = d); plates = true, N = N, kwargs...)
        W[i]       = excess_free_energy(sol)
        U[i]       = (W[i] - 2 * W_single) * pi * R_derj
        success[i] = sol.success
    end
    return U, success, W, W_single
end
