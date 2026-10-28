# Pseudo-arclength continuation in phi_inf for the two-sphere BVP.
#
# At fixed gap d the two-sphere adsorption problem has an S-shaped solution set:
# the thin-film branch, an unstable middle branch and the capillary-bridged
# branch are three points on ONE curve in (phi, phi_inf) space.  Natural
# parameter continuation -- solving G(x, phi_inf) = 0 for a prescribed sequence
# of phi_inf -- assumes that curve is a graph over phi_inf, so it stops dead at
# each limit point, where G_x is singular and the curve has a vertical tangent.
# (Measured at R = 5, chi_b = 2.8, omega_b = -0.5, d = 1.86: the thin branch
# ends at phi_inf ~ 0.0845 and the bridged branch at ~ 0.0798, with the Newton
# residual plateauing at |F|inf ~ 3e-2 and ||J^-1 F|| / ||F|| ~ 1e6.)
#
# Pseudo-arclength drops the graph assumption: x and phi_inf become unknowns
# TOGETHER, parameterized by an approximate arclength s along the curve, with
# one extra scalar equation pinning the step length along the previous tangent.
# Newton on the augmented system uses the BORDERED Jacobian
#
#     M = [ G_x    G_lam ]
#         [ t_x'   t_lam ]
#
# which stays nonsingular at a simple fold even though G_x does not: there G_x
# has a one-dimensional null space and G_lam leaves its range, and the bordering
# row/column supply exactly the missing rank.  So the sweep turns the corner
# onto the middle branch and comes back out on the other stable branch.
#
# WHY phi_inf AND NOT d.  build_geometry places sphere 2 at z_c2 = 2R + d and
# sizes z_grid from d, so both the unknown count and the inside/fluid masks are
# functions of d -- (x, d) is not a curve in a fixed R^(n+1), and the ghost
# stencils jump as a surface crosses cells.  phi_inf enters build_geometry only
# through decay_length(p), and only to pick the DEFAULT rho_max / z_buf; when
# the caller passes both explicitly the mesh is phi_inf-independent, the
# dimension is fixed and dG/dphi_inf is smooth.  So the continuation runs at
# fixed d on a geometry built once, mutating only geom.p.
#
# CALLER CONTRACT: build `geom` with explicit `rho_max` and `z_buf`.  Everything
# here mutates `geom.p` and nothing else; the grid, masks and ghost stencils are
# reused as built.  geom.p is restored on exit.

# ------------------------------------------------------------------
# Scalar diagnostics
# ------------------------------------------------------------------

# On-axis phi at the midplane: the order parameter that separates a dilute gap
# from a capillary bridge.  Cheap, and monotone along the S-curve.
gap_center_phi(geom::Geometry, phi::AbstractMatrix) =
    phi[1, argmin(abs.(geom.z_grid .- 0.5 * (geom.z_c1 + geom.z_c2)))]

# ------------------------------------------------------------------
# Trace
# ------------------------------------------------------------------

# Scalar-only record of a sweep.  The states themselves are NOT retained: at
# production resolution (n ~ 4e5) a few hundred of them would be hundreds of MB,
# so solutions are handed to `on_hit` as the curve crosses each target instead.
#
# `segment` counts limit points passed, so segment 0 is the branch the sweep
# started on and `iseven(segment)` labels the stable branches -- that parity IS
# the sign of det(G_x), which flips at every simple fold, so no eigenvalue or
# determinant computation is needed.
struct ArcTrace
    lambda::Vector{Float64}          # phi_inf at each converged point
    s::Vector{Float64}               # arclength
    segment::Vector{Int}
    gap_center::Vector{Float64}
    tangent_lambda::Vector{Float64}  # d(phi_inf)/ds; sign change => fold
    newton_iters::Vector{Int}
    folds::Vector{Float64}           # phi_inf at each limit point
    fold_segment::Vector{Int}        # segment index entered at that fold
    n_newton::Int                    # total bordered solves
    status::Symbol                   # :target | :budget | :stalled | :diverged
end

ArcTrace() = ArcTrace(Float64[], Float64[], Int[], Float64[], Float64[], Int[],
                      Float64[], Int[], 0, :empty)

# ------------------------------------------------------------------
# Bordered linear algebra
# ------------------------------------------------------------------

# Full (n+1) assembly: one dense row and column appended to the sparse J.
# Nonsingular at a simple fold, which is the whole point, but it perturbs the
# sparsity pattern, so it is the fallback rather than the default.
function _bordered_solve_full(J::SparseMatrixCSC, Gl::Vector{Float64},
                              tx::Vector{Float64}, tl::Float64,
                              r::Vector{Float64}, rn::Float64)
    n = length(tx)
    Ir, Jc, V = findnz(J)
    append!(Ir, 1:n);              append!(Jc, fill(n + 1, n)); append!(V, Gl)
    append!(Ir, fill(n + 1, n));   append!(Jc, 1:n);            append!(V, tx)
    push!(Ir, n + 1); push!(Jc, n + 1); push!(V, tl)
    sol = sparse(Ir, Jc, V, n + 1, n + 1) \ vcat(r, rn)
    return sol[1:n], sol[n + 1]
end

# Keller bordering: one factorization of J, two triangular solves, then the
# scalar unknown from the bordering row.  Same factorization count per Newton
# step as the plain solver, so the smooth stretches cost what they cost today.
# The pivot `tl - tx'y2` collapses exactly AT a fold (where J is singular and y2
# blows up), which is the trigger for the full assembly.
function _bordered_solve(J::SparseMatrixCSC, Gl::Vector{Float64},
                         tx::Vector{Float64}, tl::Float64,
                         r::Vector{Float64}, rn::Float64)
    fac = try
        lu(J)
    catch
        return _bordered_solve_full(J, Gl, tx, tl, r, rn)
    end
    y1 = fac \ r
    y2 = fac \ Gl
    (all(isfinite, y1) && all(isfinite, y2)) ||
        return _bordered_solve_full(J, Gl, tx, tl, r, rn)
    txy2 = dot(tx, y2)
    piv  = tl - txy2
    if !isfinite(piv) || abs(piv) <= 1e-8 * max(abs(tl) + abs(txy2), 1.0)
        return _bordered_solve_full(J, Gl, tx, tl, r, rn)
    end
    dl = (rn - dot(tx, y1)) / piv
    dx = y1 .- dl .* y2
    (isfinite(dl) && all(isfinite, dx)) ||
        return _bordered_solve_full(J, Gl, tx, tl, r, rn)
    return dx, dl
end

# Scaled inner product <(a,al),(b,be)> = a.b/n + al*be.  Without the 1/n the
# O(n) field part swamps the single scalar parameter and ds becomes a pure
# x-step with no phi_inf progress; this is the standard AUTO/MATCONT weighting.
_arc_dot(a::Vector{Float64}, al::Float64, b::Vector{Float64}, be::Float64) =
    dot(a, b) / length(a) + al * be

# ------------------------------------------------------------------
# Fixed-phi_inf Newton on an existing geometry
# ------------------------------------------------------------------

"""
    polish_at_phi_inf(geom, x0, lambda; f_tol, maxiter, work) -> (x, converged)

Plain damped Newton at fixed `phi_inf = lambda` on the geometry as built,
warm-started from `x0`.  Used to land an interpolated arclength point exactly on
the branch without paying for another `build_geometry`.

Leaves `geom.p.phi_inf == lambda` on return, so a following
`excess_free_energy_two_spheres(Dict(:geom => geom, :phi => scatter(x, geom)))`
sees the right bulk reference.

For many calls on one mesh, hoist the two allocations out: `work = (phi_buf,
lap_buf)` (two `Nrho x Nz` matrices) and, at `bc_reg == 0`, `A =
build_linear_operator(geom)` — that operator depends on neither `phi` nor
`phi_inf`, so reassembling it per call is pure overhead at production sizes.
"""
function polish_at_phi_inf(geom::Geometry, x0::AbstractVector{Float64},
                           lambda::Real;
                           f_tol::Real = 1e-6, maxiter::Int = 30,
                           work = nothing, A = nothing)
    p_base = geom.p
    geom.p = with(p_base; phi_inf = Float64(lambda))
    p = geom.p
    n = length(x0)
    phi_buf, lap_buf = work === nothing ?
        (Matrix{Float64}(undef, geom.Nrho, geom.Nz),
         Matrix{Float64}(undef, geom.Nrho, geom.Nz)) : work
    A_const = geom.bc_reg == 0.0 ?
        (A === nothing ? build_linear_operator(geom) : A) : nothing

    x = clamp.(collect(float.(x0)), p.eps, 1.0 - p.eps)
    r  = Vector{Float64}(undef, n)
    rt = Vector{Float64}(undef, n)
    xt = Vector{Float64}(undef, n)
    diagv = Vector{Float64}(undef, n)

    _res!(F, xx) = _fill_residual!(F, xx, geom, phi_buf, lap_buf)
    _res!(r, x)
    nrm = maximum(abs, r)
    for _ in 1:maxiter
        nrm < f_tol && break
        @inbounds for k in 1:n
            diagv[k] = d2fb_dphi2(x[k], p)
        end
        Jlin = A_const === nothing ?
            build_linear_operator(geom; sd = ghost_slope_derivs(geom, phi_buf)) : A_const
        δ = try
            (Jlin - Diagonal(diagv)) \ r
        catch
            break
        end
        all(isfinite, δ) || break
        # Fraction-to-boundary damping: feasible by construction, so no clamp.
        damp = _max_feasible_step(x, δ, p.eps, 1.0 - p.eps)
        accepted = false
        for _ in 1:30
            @inbounds @. xt = x - damp * δ
            _res!(rt, xt)
            nn = maximum(abs, rt)
            if nn < nrm
                @inbounds @. x = xt
                nrm = nn
                accepted = true
                break
            end
            damp *= 0.5
        end
        accepted || break
        _res!(r, x)
    end
    return x, nrm < f_tol
end

# Residual of the two-sphere BVP at the CURRENT geom.p, written into F.  Mirrors
# make_residual_inplace, but reads mu_inf from geom.p instead of capturing it, so
# phi_inf can vary between calls.
function _fill_residual!(F::Vector{Float64}, x::AbstractVector{Float64},
                         geom::Geometry, phi_buf::Matrix{Float64},
                         lap_buf::Matrix{Float64})
    p = geom.p
    scatter!(phi_buf, x, geom)
    laplacian!(lap_buf, phi_buf, geom)
    mu = dfb_dphi(p.phi_inf, p)          # == compute_mu_inf(p) / p.nb
    k = 0
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        if geom.unk_mask[i, j]
            k += 1
            F[k] = p.kappa * lap_buf[i, j] - (dfb_dphi(phi_buf[i, j], p) - mu)
        end
    end
    return F
end

# ------------------------------------------------------------------
# The sweep
# ------------------------------------------------------------------

"""
    arclength_phi_inf(geom, x0, lambda0, lambda_stop; kwargs...) -> ArcTrace

Pseudo-arclength continuation of the two-sphere BVP in `phi_inf`, at the fixed
gap `geom.d`, from `lambda0` towards `lambda_stop`, passing through any number
of limit points.  `x0` must be a converged solution at `lambda0` (it is polished
first regardless).

`geom` must have been built with explicit `rho_max` and `z_buf` -- see the file
header.  `geom.p` is mutated during the sweep and restored on return.

Keyword arguments:
- `targets`   : phi_inf values to report solutions at.
- `on_hit`    : `f(x, lambda, segment, stable) -> Any`, called each time the
                curve crosses a target.  `x` is linearly interpolated in `s`
                between the two bracketing points and is only valid during the
                call; polish it with `polish_at_phi_inf` if you need a solution.
                A target on an S-curve is hit up to three times, once per branch;
                `segment` tells them apart.
- `on_point`  : `f(x, lambda, s, segment, stable)`, called once for the seed and
                then once per accepted continuation point, so the resulting
                series is index-aligned with the returned `ArcTrace` vectors.
                Same borrowed-`x` rule as `on_hit`.  Use it to record a quantity
                along the whole curve (a free energy, say) without retaining the
                states.
- `dlambda0`  : phi_inf progress aimed for on the first step (defaults to
                1/50 of the sweep range).  Sets the arclength scale.
- `newton_target`, `ds_min_factor`, `ds_max_factor`, `max_steps`, `max_shrink`,
  `f_tol`, `maxiter`, `verbose`.
"""
function arclength_phi_inf(geom::Geometry, x0::AbstractVector{Float64},
                           lambda0::Real, lambda_stop::Real;
                           targets::AbstractVector{<:Real} = Float64[],
                           on_hit = nothing,
                           on_point = nothing,
                           dlambda0::Real = 0.0,
                           newton_target::Int = 3,
                           ds_min_factor::Real = 1e-3,
                           ds_max_factor::Real = 10.0,
                           max_steps::Int = 400,
                           max_shrink::Int = 12,
                           f_tol::Real = 1e-6,
                           maxiter::Int = 30,
                           verbose::Bool = false)

    p_base = geom.p
    n = length(x0)
    n == n_unknowns(geom) ||
        throw(DimensionMismatch("x0 has $n entries, geometry has $(n_unknowns(geom)) unknowns"))
    lam_stop = Float64(lambda_stop)
    dir = sign(lam_stop - Float64(lambda0))
    dir == 0.0 && throw(ArgumentError("lambda0 == lambda_stop; nothing to continue"))

    # Box-size sanity: rho_max / z_buf were fixed at build time, but the
    # correlation length grows towards the spinodal, so a box sized for the
    # starting phi_inf can be too small at the far end of the sweep.
    ell_stop = decay_length(with(p_base; phi_inf = lam_stop))
    if geom.rho_grid[end] < p_base.R + 3.0 * ell_stop
        @warn "arclength_phi_inf: rho_max may be too small at the far end of the sweep" rho_max=geom.rho_grid[end] needed=p_base.R + 4.0 * ell_stop phi_inf=lam_stop
    end

    reg     = geom.bc_reg
    A_const = reg == 0.0 ? build_linear_operator(geom) : nothing
    phi_buf = Matrix{Float64}(undef, geom.Nrho, geom.Nz)
    lap_buf = Matrix{Float64}(undef, geom.Nrho, geom.Nz)

    # Residual at (x, lam): sets geom.p, so every call leaves the geometry at lam.
    function res!(F, xx, lam)
        geom.p = with(p_base; phi_inf = lam)
        _fill_residual!(F, xx, geom, phi_buf, lap_buf)
    end
    # dG/dlambda at fixed x, central difference.  Two O(n) residual evaluations
    # against one sparse LU -- not worth an analytic derivative.  Restores
    # geom.p to `lam` on the way out.
    function dres_dlam!(Gl, Fp, Fm, xx, lam)
        h = 1e-6 * max(abs(lam), 1e-3)
        res!(Fp, xx, lam + h)
        res!(Fm, xx, lam - h)
        @inbounds @. Gl = (Fp - Fm) / (2h)
        return Gl
    end
    function jac(xx)
        p = geom.p
        diagv = Vector{Float64}(undef, n)
        @inbounds for k in 1:n
            diagv[k] = d2fb_dphi2(xx[k], p)
        end
        Jlin = A_const === nothing ?
            build_linear_operator(geom; sd = ghost_slope_derivs(geom, phi_buf)) : A_const
        return Jlin - Diagonal(diagv)
    end

    tgt = collect(float.(targets))

    # ---- state ----
    lam = Float64(lambda0)
    x, ok0 = polish_at_phi_inf(geom, x0, lam; f_tol = f_tol, maxiter = maxiter,
                               work = (phi_buf, lap_buf), A = A_const)
    ok0 || (@warn "arclength_phi_inf: starting point did not converge at phi_inf = $lam";
            return ArcTrace())

    r   = Vector{Float64}(undef, n)
    rt  = Vector{Float64}(undef, n)
    Gl  = Vector{Float64}(undef, n)
    Fp  = Vector{Float64}(undef, n)
    Fm  = Vector{Float64}(undef, n)
    xt  = Vector{Float64}(undef, n)
    x_prev = Vector{Float64}(undef, n)

    # ---- initial tangent: border row (0, 1) gives (-J\Gl, 1) ----
    res!(r, x, lam)
    dres_dlam!(Gl, Fp, Fm, x, lam)
    res!(r, x, lam)
    tx, tl = _bordered_solve(jac(x), Gl, zeros(n), 1.0, zeros(n), 1.0)
    tn = sqrt(_arc_dot(tx, tl, tx, tl))
    isfinite(tn) && tn > 0 || return ArcTrace()
    tx ./= tn; tl /= tn
    if sign(tl) != dir            # orient towards lambda_stop
        tx .= .-tx; tl = -tl
    end
    # Bordering ROW of M, i.e. the tangent expressed in the scaled inner
    # product: <(a,al),(b,be)> = a.b/n + al*be, so the row is [tx'/n, tl].  It
    # must match how `rn` is formed below -- with an unscaled row the arclength
    # equation is corrected n times too slowly and the corrector never reaches
    # f_tol even though the BVP residual is at machine precision.
    bx = tx ./ n

    dlam0 = dlambda0 > 0 ? Float64(dlambda0) : abs(lam_stop - lam) / 50
    ds    = abs(tl) > 1e-12 ? dlam0 / abs(tl) : dlam0
    ds_min = ds * ds_min_factor
    ds_max = ds * ds_max_factor

    lam_hist = Float64[lam]; s_hist = Float64[0.0]; seg_hist = Int[0]
    gc_hist  = Float64[gap_center_phi(geom, phi_buf)]
    tl_hist  = Float64[tl]; it_hist = Int[0]
    folds = Float64[]; fold_seg = Int[]
    # Report the seed too, so a caller's on_point series lines up index-for-index
    # with the ArcTrace vectors instead of being offset by one.
    on_point === nothing || on_point(x, lam, 0.0, 0, true)

    s = 0.0
    segment = 0
    n_newton = 0
    shrinks = 0
    status = :budget

    for step in 1:max_steps
        # ---- predictor ----
        copyto!(x_prev, x)
        lam_prev = lam
        tl_prev  = tl
        xp  = clamp.(x .+ ds .* tx, p_base.eps, 1.0 - p_base.eps)
        lamp = lam + ds * tl

        # ---- bordered Newton corrector ----
        xc = xp; lamc = lamp
        converged = false
        iters = 0
        if p_base.eps < lamc < 1.0 - p_base.eps
            res!(r, xc, lamc)
            rn = _arc_dot(xc .- x_prev, lamc - lam_prev, tx, tl) - ds
            nrm = max(maximum(abs, r), abs(rn))
            for it in 1:maxiter
                iters = it
                if nrm < f_tol
                    converged = true
                    break
                end
                dres_dlam!(Gl, Fp, Fm, xc, lamc)
                res!(r, xc, lamc)
                dx, dl = _bordered_solve(jac(xc), Gl, bx, tl, r, rn)
                n_newton += 1
                (all(isfinite, dx) && isfinite(dl)) || break
                # Fraction-to-boundary on the field part (feasible, no clamp); the
                # phi_inf part lamt is bounds-checked separately below.
                damp = _max_feasible_step(xc, dx, p_base.eps, 1.0 - p_base.eps)
                accepted = false
                for _ in 1:30
                    @inbounds @. xt = xc - damp * dx
                    lamt = lamc - damp * dl
                    if !(p_base.eps < lamt < 1.0 - p_base.eps)
                        damp *= 0.5
                        continue
                    end
                    res!(rt, xt, lamt)
                    rnt = _arc_dot(xt .- x_prev, lamt - lam_prev, tx, tl) - ds
                    nn = max(maximum(abs, rt), abs(rnt))
                    if nn < nrm
                        xc = copy(xt); lamc = lamt
                        nrm = nn; rn = rnt
                        copyto!(r, rt)
                        accepted = true
                        break
                    end
                    damp *= 0.5
                end
                # A stalled line search here means the STEP was too long, not
                # that we are stuck at a fold -- the bordered system is regular
                # there.  So bail out of the corrector and let the outer loop
                # halve ds, rather than giving up on the branch.
                accepted || break
            end
        end

        if !converged
            shrinks += 1
            ds *= 0.5
            if ds < ds_min || shrinks > max_shrink
                status = :stalled
                break
            end
            verbose && @printf("  step %3d: corrector failed, ds -> %.3e\n", step, ds)
            continue
        end
        shrinks = 0

        # ---- accept ----
        x = xc; lam = lamc; s += ds
        res!(r, x, lam)                       # leaves phi_buf at the new state

        # ---- new tangent ----
        dres_dlam!(Gl, Fp, Fm, x, lam)
        res!(r, x, lam)
        tx_new, tl_new = _bordered_solve(jac(x), Gl, bx, tl, zeros(n), 1.0)
        tn = sqrt(_arc_dot(tx_new, tl_new, tx_new, tl_new))
        if !(isfinite(tn) && tn > 0)
            status = :diverged
            break
        end
        tx_new ./= tn; tl_new /= tn
        if _arc_dot(tx_new, tl_new, tx, tl) < 0     # keep going the same way
            tx_new .= .-tx_new; tl_new = -tl_new
        end
        tx = tx_new; tl = tl_new
        bx = tx ./ n

        # ---- fold detection: d(phi_inf)/ds changes sign ----
        crossed_fold = tl * tl_prev < 0
        if crossed_fold
            segment += 1
            push!(folds, _fold_lambda(lam_hist, s_hist, lam, s))
            push!(fold_seg, segment)
            verbose && @printf("  step %3d: FOLD at phi_inf ~ %.6f (segment -> %d)\n",
                               step, folds[end], segment)
        end

        push!(lam_hist, lam); push!(s_hist, s); push!(seg_hist, segment)
        push!(gc_hist, gap_center_phi(geom, phi_buf))
        push!(tl_hist, tl); push!(it_hist, iters)

        # geom.p is at `lam` and phi_buf holds this state, so a callback can
        # evaluate the free energy directly without re-scattering.
        on_point === nothing || on_point(x, lam, s, segment, iseven(segment))

        # ---- target crossings on [lam_prev, lam] ----
        if on_hit !== nothing && !isempty(tgt)
            lo, hi = minmax(lam_prev, lam)
            for t in tgt
                (lo <= t <= hi && lam != lam_prev) || continue
                w = (t - lam_prev) / (lam - lam_prev)
                @inbounds @. xt = (1.0 - w) * x_prev + w * x
                on_hit(xt, t, segment, iseven(segment))
            end
        end

        verbose && @printf("  step %3d: phi_inf=%.6f  s=%.4f  seg=%d  gapc=%.4f  ds=%.3e  it=%d\n",
                           step, lam, s, segment, gc_hist[end], ds, iters)

        # ---- termination ----
        if dir * (lam - lam_stop) > 0 && dir * tl > 0
            status = :target
            break
        end

        # ---- adaptive step ----
        ds = clamp(ds * clamp(sqrt(newton_target / max(iters, 1)), 0.5, 2.0),
                   ds_min, ds_max)
    end

    geom.p = p_base
    return ArcTrace(lam_hist, s_hist, seg_hist, gc_hist, tl_hist, it_hist,
                    folds, fold_seg, n_newton, status)
end

# Limit-point phi_inf from the vertex of a parabola through the last three
# (s, lambda) points; lambda is at an extremum there, so linear interpolation
# would be badly biased.
function _fold_lambda(lam_hist::Vector{Float64}, s_hist::Vector{Float64},
                      lam::Float64, s::Float64)
    length(lam_hist) < 2 && return lam
    s1, l1 = s_hist[end - 1], lam_hist[end - 1]
    s2, l2 = s_hist[end],     lam_hist[end]
    s3, l3 = s,               lam
    (s1 != s2 && s2 != s3 && s1 != s3) || return l2
    # vertex of the interpolating parabola in s
    a = ((l1 - l2) / (s1 - s2) - (l2 - l3) / (s2 - s3)) / (s1 - s3)
    abs(a) < 1e-300 && return l2
    b = (l1 - l2) / (s1 - s2) - a * (s1 + s2)
    sv = -b / (2a)
    c  = l1 - a * s1^2 - b * s1
    lv = a * sv^2 + b * sv + c
    smin, smax = minmax(s1, s3)
    return (smin <= sv <= smax && isfinite(lv)) ? lv : l2
end
