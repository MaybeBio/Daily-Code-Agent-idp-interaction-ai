# 2D axisymmetric BVP solver for two spheres in a Flory-Huggins bulk.
#
# Mirrors viz_util/bvp_two_spheres.py: two identical spheres of radius R sit
# on the z-axis at z = 0 and z = 2R + d, and the equilibrium volume fraction
# phi(rho, z) satisfies
#
#     kappa * Laplacian(phi) = dfb/dphi - mu_inf / nb
#
# with cylindrical Laplacian, Robin (adsorption) BC on each sphere via the
# quadratic ghost-cell immersed-boundary method, Dirichlet phi = phi_inf on
# the outer cylinder + z end-caps, and dphi/drho = 0 on the symmetry axis.
#
# The nonlinear solve uses NonlinearSolve.jl (NewtonRaphson + GMRES) warm-
# started from the linear-response superposition of two psi_linear profiles,
# with an omega_b continuation fallback identical to the Python version.

# ------------------------------------------------------------------
# Geometry
# ------------------------------------------------------------------

mutable struct Geometry
    p::Parameters
    d::Float64
    Nrho::Int
    Nz::Int
    rho_grid::Vector{Float64}
    z_grid::Vector{Float64}
    drho::Float64
    dz::Float64
    z_c1::Float64
    z_c2::Float64
    inside_mask::Matrix{Bool}
    fluid_mask::Matrix{Bool}
    dir_mask::Matrix{Bool}
    unk_mask::Matrix{Bool}
    # one entry per ghost cell
    ghost_i::Vector{Int}
    ghost_j::Vector{Int}
    ghost_i0_1::Vector{Int}
    ghost_j0_1::Vector{Int}
    ghost_wts_1::Matrix{Float64}   # (4, n_ghost), corner order (i0,j0),(i0+1,j0),(i0,j0+1),(i0+1,j0+1)
    ghost_i0_2::Vector{Int}
    ghost_j0_2::Vector{Int}
    ghost_wts_2::Matrix{Float64}
    ghost_coef_phi1::Vector{Float64}
    ghost_coef_phi2::Vector{Float64}
    ghost_coef_bc::Vector{Float64}
    ghost_used_linear::Vector{Bool}
    ghost_delta::Vector{Float64}
    ghost_d1::Vector{Float64}
    ghost_d2::Vector{Float64}
    # Surface-BC regularization strength (0 = plain Robin BC).  When > 0 the
    # adsorption BC carries a + bc_reg*log(phi_s/(1-phi_s)) term (see scatter!).
    bc_reg::Float64
end

n_ghost(g::Geometry) = length(g.ghost_i)
n_unknowns(g::Geometry) = count(g.unk_mask)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

function decay_length(p::Parameters)
    d2 = d2fb_dphi2(p.phi_inf, p)
    d2 > 0.05 && return sqrt(p.kappa / d2)
    return sqrt(p.kappa * p.nb / 0.05)
end

const _NBR_OFFSETS = ((1, 0), (-1, 0), (0, 1), (0, -1))

# Walk a probe point outward along (u_rho, u_z) until its 4-cell bilinear
# stencil contains no inside-sphere cells.  Returns (d_used, i0, j0, w_tuple,
# clean::Bool) with `w_tuple` ordered (i0,j0),(i0+1,j0),(i0,j0+1),(i0+1,j0+1).
# `d_max` caps the walk (e.g. at the entry into the OTHER sphere near contact:
# a "clean" stencil found beyond it would sit on the far side of that sphere).
function find_clean_stencil(rho_G::Real, z_G::Real,
                            u_rho::Real, u_z::Real,
                            d_initial::Real, step::Real,
                            rho_grid::AbstractVector{Float64},
                            z_grid::AbstractVector{Float64},
                            drho::Real, dz::Real,
                            Nrho::Int, Nz::Int,
                            inside_mask::AbstractMatrix{Bool};
                            max_iter::Int = 80, d_max::Real = Inf)
    d_use = float(d_initial)
    i0 = 1; j0 = 1
    w = (0.0, 0.0, 0.0, 0.0)
    clean = false
    rho_lo, rho_hi = rho_grid[1], rho_grid[end]
    z_lo,   z_hi   = z_grid[1],   z_grid[end]
    @inbounds for _ in 1:max_iter
        d_use > d_max && break
        rho_I = rho_G + d_use * u_rho
        z_I   = z_G   + d_use * u_z
        rho_Ic = clamp(rho_I, rho_lo, rho_hi - 1e-12 * drho)
        z_Ic   = clamp(z_I,   z_lo,   z_hi   - 1e-12 * dz)
        ti = (rho_Ic - rho_lo) / drho
        tj = (z_Ic   - z_lo)   / dz
        i0 = clamp(Int(floor(ti)) + 1, 1, Nrho - 1)
        j0 = clamp(Int(floor(tj)) + 1, 1, Nz   - 1)
        t_rho = ti - (i0 - 1)
        t_z   = tj - (j0 - 1)
        w = ((1.0 - t_rho) * (1.0 - t_z),
             t_rho         * (1.0 - t_z),
             (1.0 - t_rho) * t_z,
             t_rho         * t_z)
        stencil_clean = !(inside_mask[i0,     j0]     ||
                          inside_mask[i0 + 1, j0]     ||
                          inside_mask[i0,     j0 + 1] ||
                          inside_mask[i0 + 1, j0 + 1])
        if stencil_clean
            clean = true
            break
        end
        d_use += step
    end
    return d_use, i0, j0, w, clean
end

# Distance along the ray (rho0, z0) + t*(u_rho, u_z) to its entry into the
# sphere of radius R centred on the axis at z = zo (Inf if the ray misses it).
function _ray_sphere_entry(rho0::Real, z0::Real, u_rho::Real, u_z::Real,
                           zo::Real, R::Real)
    b = rho0 * u_rho + (z0 - zo) * u_z
    c = rho0^2 + (z0 - zo)^2 - R^2
    disc = b * b - c
    disc <= 0.0 && return Inf
    t = -b - sqrt(disc)
    return t > 0.0 ? t : Inf
end

# Closed-form (alpha1, alpha2, beta) such that
#     phi(G) = alpha1 * phi(P1) + alpha2 * phi(P2) + beta * (omega_b/kappa)
# implements the quadratic phi(s) = a + b*s + c*s^2 interpolating phi(P1),
# phi(P2) and the Robin BC dphi/ds|_{s=delta} = omega_b/kappa at the surface.
function quadratic_ghost_coefficients(delta::Real, d1::Real, d2::Real)
    K       = d1 + d2 - 2.0 * delta
    inv_K   = 1.0 / K
    inv_dd  = 1.0 / (d1 - d2)
    c_phi1  =  inv_dd * inv_K
    c_phi2  = -inv_dd * inv_K
    c_bc    = -inv_K
    b_phi1  = -2.0 * delta * c_phi1
    b_phi2  = -2.0 * delta * c_phi2
    b_bc    = 1.0 - 2.0 * delta * c_bc
    α1      = 1.0 - d1 * b_phi1 - d1 * d1 * c_phi1
    α2      =      -d1 * b_phi2 - d1 * d1 * c_phi2
    β       =      -d1 * b_bc   - d1 * d1 * c_bc
    return α1, α2, β
end

# ------------------------------------------------------------------
# build_geometry
# ------------------------------------------------------------------

function build_geometry(p::Parameters, d::Real;
                        rho_max = nothing, z_buf = nothing,
                        Nrho::Int = 80, Nz::Int = 200,
                        bc_reg::Real = 0.0)
    R   = p.R
    ell = decay_length(p)
    rho_max = rho_max === nothing ? R + 4.0 * ell : Float64(rho_max)
    z_buf   = z_buf   === nothing ? 4.0 * ell     : Float64(z_buf)

    z_c1 = 0.0
    z_c2 = 2.0 * R + d

    rho_grid = collect(range(0.0, rho_max; length = Nrho))
    z_grid   = collect(range(z_c1 - R - z_buf, z_c2 + R + z_buf; length = Nz))
    drho = rho_grid[2] - rho_grid[1]
    dz   = z_grid[2]   - z_grid[1]

    inside_mask = falses(Nrho, Nz)
    @inbounds for j in 1:Nz, i in 1:Nrho
        rho_i = rho_grid[i]
        z_j   = z_grid[j]
        r1_sq = rho_i^2 + (z_j - z_c1)^2
        r2_sq = rho_i^2 + (z_j - z_c2)^2
        inside_mask[i, j] = (r1_sq < R^2) || (r2_sq < R^2)
    end
    fluid_mask = .!inside_mask

    dir_mask = falses(Nrho, Nz)
    @inbounds dir_mask[Nrho, :] .= true
    @inbounds dir_mask[:, 1]    .= true
    @inbounds dir_mask[:, Nz]   .= true
    dir_mask .= dir_mask .& fluid_mask
    unk_mask  = fluid_mask .& .!dir_mask

    step     = max(drho, dz)
    max_iter = 80
    sep_min  = 0.9 * step
    K_min    = 1e-9 * step

    gi  = Int[]
    gj  = Int[]
    gi1 = Int[]; gj1 = Int[]
    gi2 = Int[]; gj2 = Int[]
    gw1 = NTuple{4, Float64}[]
    gw2 = NTuple{4, Float64}[]
    α1v = Float64[]; α2v = Float64[]; βv = Float64[]
    gulin = Bool[]
    gdelta = Float64[]; gd1 = Float64[]; gd2 = Float64[]

    @inbounds for j in 1:Nz, i in 1:Nrho
        inside_mask[i, j] || continue
        has_fluid = false
        for (di, dj) in _NBR_OFFSETS
            ii = i + di
            jj = j + dj
            if 1 <= ii <= Nrho && 1 <= jj <= Nz && fluid_mask[ii, jj]
                has_fluid = true
                break
            end
        end
        has_fluid || continue

        rho_G = rho_grid[i]
        z_G   = z_grid[j]
        dist1 = hypot(rho_G, z_G - z_c1)
        dist2 = hypot(rho_G, z_G - z_c2)
        zc, Gc = dist1 <= dist2 ? (z_c1, dist1) : (z_c2, dist2)
        Gc < 1e-12 && continue
        u_rho = rho_G / Gc
        u_z   = (z_G - zc) / Gc
        delta = R - Gc

        # Near contact the outward normal of one sphere runs into the other:
        # cap the probe walk at that entry so the BC closure never couples to
        # fluid on the far side of the neighbouring sphere.
        d_cap = _ray_sphere_entry(rho_G, z_G, u_rho, u_z,
                                  zc == z_c1 ? z_c2 : z_c1, R)

        d1u, i0_1, j0_1, w1, clean1 = find_clean_stencil(
            rho_G, z_G, u_rho, u_z, 2.0 * delta, step,
            rho_grid, z_grid, drho, dz, Nrho, Nz, inside_mask;
            max_iter = max_iter, d_max = d_cap)
        d2u, i0_2, j0_2, w2, clean2 = find_clean_stencil(
            rho_G, z_G, u_rho, u_z, d1u + step, step,
            rho_grid, z_grid, drho, dz, Nrho, Nz, inside_mask;
            max_iter = max_iter, d_max = d_cap)

        if !clean1
            # Sub-resolution neck (spheres at/near contact): no clean stencil
            # exists between this ghost and the other sphere.  Close the Robin
            # BC through the ghost's own fluid 4-neighbour (exists by
            # construction, and is never an inside/ghost cell), taken most
            # aligned with the outward normal.
            ii, jj = i, j
            best = -Inf
            for (di, dj) in _NBR_OFFSETS
                i2 = i + di
                j2 = j + dj
                (1 <= i2 <= Nrho && 1 <= j2 <= Nz && fluid_mask[i2, j2]) || continue
                al = di * drho * u_rho + dj * dz * u_z
                if al > best
                    best = al; ii = i2; jj = j2
                end
            end
            i0_1 = clamp(ii, 1, Nrho - 1)
            j0_1 = clamp(jj, 1, Nz - 1)
            corner = 1 + (ii - i0_1) + 2 * (jj - j0_1)
            w1 = ntuple(c -> c == corner ? 1.0 : 0.0, 4)
            d1u = hypot((ii - i) * drho, (jj - j) * dz)
            d2u = d1u + step
            use_linear = true
            α1, α2, β = 1.0, 0.0, -d1u
            i0_2, j0_2 = i0_1, j0_1
            w2 = (0.0, 0.0, 0.0, 0.0)
        else
            use_linear = (!clean2) ||
                         ((d2u - d1u) < sep_min) ||
                         (abs(d1u + d2u - 2.0 * delta) < K_min) ||
                         (abs(d1u - d2u) < K_min)
            if use_linear
                α1, α2, β = 1.0, 0.0, -d1u
                i0_2, j0_2 = i0_1, j0_1
                w2 = (0.0, 0.0, 0.0, 0.0)
            else
                α1, α2, β = quadratic_ghost_coefficients(delta, d1u, d2u)
            end
        end

        push!(gi, i);   push!(gj, j)
        push!(gi1, i0_1); push!(gj1, j0_1)
        push!(gi2, i0_2); push!(gj2, j0_2)
        push!(gw1, w1); push!(gw2, w2)
        push!(α1v, α1); push!(α2v, α2); push!(βv, β)
        push!(gulin, use_linear)
        push!(gdelta, delta)
        push!(gd1, d1u); push!(gd2, d2u)
    end

    ng = length(gi)
    ghost_wts_1 = Matrix{Float64}(undef, 4, ng)
    ghost_wts_2 = Matrix{Float64}(undef, 4, ng)
    @inbounds for k in 1:ng
        w = gw1[k]
        ghost_wts_1[1, k] = w[1]; ghost_wts_1[2, k] = w[2]
        ghost_wts_1[3, k] = w[3]; ghost_wts_1[4, k] = w[4]
        w = gw2[k]
        ghost_wts_2[1, k] = w[1]; ghost_wts_2[2, k] = w[2]
        ghost_wts_2[3, k] = w[3]; ghost_wts_2[4, k] = w[4]
    end

    return Geometry(
        p, Float64(d), Nrho, Nz, rho_grid, z_grid, drho, dz,
        z_c1, z_c2,
        inside_mask, fluid_mask, dir_mask, unk_mask,
        gi, gj,
        gi1, gj1, ghost_wts_1,
        gi2, gj2, ghost_wts_2,
        α1v, α2v, βv, gulin, gdelta, gd1, gd2,
        Float64(bc_reg),
    )
end

# ------------------------------------------------------------------
# Scatter: assemble full (Nrho, Nz) field from unknown vector
# ------------------------------------------------------------------

function scatter!(phi::AbstractMatrix, x_unk::AbstractVector, geom::Geometry)
    p = geom.p
    fill!(phi, p.phi_inf)
    k = 0
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        if geom.unk_mask[i, j]
            k += 1
            phi[i, j] = x_unk[k]
        end
    end
    bc0 = p.omega_b / p.kappa
    reg = geom.bc_reg
    @inbounds for k in 1:n_ghost(geom)
        i  = geom.ghost_i[k]
        j  = geom.ghost_j[k]
        i01 = geom.ghost_i0_1[k]; j01 = geom.ghost_j0_1[k]
        i02 = geom.ghost_i0_2[k]; j02 = geom.ghost_j0_2[k]
        phi_P1 = geom.ghost_wts_1[1, k] * phi[i01,     j01]     +
                 geom.ghost_wts_1[2, k] * phi[i01 + 1, j01]     +
                 geom.ghost_wts_1[3, k] * phi[i01,     j01 + 1] +
                 geom.ghost_wts_1[4, k] * phi[i01 + 1, j01 + 1]
        phi_P2 = geom.ghost_wts_2[1, k] * phi[i02,     j02]     +
                 geom.ghost_wts_2[2, k] * phi[i02 + 1, j02]     +
                 geom.ghost_wts_2[3, k] * phi[i02,     j02 + 1] +
                 geom.ghost_wts_2[4, k] * phi[i02 + 1, j02 + 1]
        # Surface-BC slope, optionally regularized to pin a unique finite
        # adsorption branch (mirrors the 1D _surf_slope / make_bc).  The surface
        # value phi_s is reconstructed from the two clean probe points P1, P2 by
        # linear extrapolation to s = delta; this is an explicit function of the
        # current iterate, so scatter! stays a pure function of x (the resulting
        # Jacobian inexactness in bc is handled by the line search in
        # _solve_newton_sparse).
        bc = bc0
        if reg != 0.0
            d1 = geom.ghost_d1[k]; d2 = geom.ghost_d2[k]
            phi_s = geom.ghost_used_linear[k] ? phi_P1 :
                    phi_P1 + (phi_P2 - phi_P1) * (geom.ghost_delta[k] - d1) / (d2 - d1)
            ps = clamp(phi_s, p.eps, 1.0 - p.eps)
            bc = (p.omega_b + reg * log(ps / (1.0 - ps))) / p.kappa
        end
        phi[i, j] = geom.ghost_coef_phi1[k] * phi_P1 +
                    geom.ghost_coef_phi2[k] * phi_P2 +
                    geom.ghost_coef_bc[k]   * bc
    end
    return phi
end

function scatter(x_unk::AbstractVector, geom::Geometry)
    phi = Matrix{Float64}(undef, geom.Nrho, geom.Nz)
    scatter!(phi, x_unk, geom)
    return phi
end

# ------------------------------------------------------------------
# Cylindrical Laplacian
# ------------------------------------------------------------------

function laplacian!(lap::AbstractMatrix, phi::AbstractMatrix, geom::Geometry)
    drho, dz = geom.drho, geom.dz
    Nrho, Nz = geom.Nrho, geom.Nz
    fill!(lap, 0.0)
    @inbounds for j in 2:Nz - 1
        for i in 2:Nrho - 1
            rho_i = geom.rho_grid[i]
            lap[i, j] = (phi[i + 1, j] - 2.0 * phi[i, j] + phi[i - 1, j]) / drho^2 +
                        (phi[i + 1, j] - phi[i - 1, j]) / (2.0 * drho * rho_i) +
                        (phi[i, j + 1] - 2.0 * phi[i, j] + phi[i, j - 1]) / dz^2
        end
    end
    # Axis i=1: with symmetric ghost phi[0,j] = phi[2,j] (1-indexed phi[1,j] is
    # on-axis, phi[2,j] is the first off-axis cell), the radial Laplacian
    # collapses to 4*(phi[2,j] - phi[1,j])/drho^2.
    @inbounds for j in 2:Nz - 1
        lap[1, j] = 4.0 * (phi[2, j] - phi[1, j]) / drho^2 +
                    (phi[1, j + 1] - 2.0 * phi[1, j] + phi[1, j - 1]) / dz^2
    end
    return lap
end

function laplacian(phi::AbstractMatrix, geom::Geometry)
    lap = similar(phi)
    laplacian!(lap, phi, geom)
    return lap
end

# ------------------------------------------------------------------
# Residual closure
# ------------------------------------------------------------------

# In-place residual for NonlinearSolve: f!(F, x, p) writes the residual at
# the unknown cells, reusing internal buffers across calls.
function make_residual_inplace(geom::Geometry, mu_inf::Float64)
    p = geom.p
    Nrho, Nz = geom.Nrho, geom.Nz
    phi_buf = Matrix{Float64}(undef, Nrho, Nz)
    lap_buf = Matrix{Float64}(undef, Nrho, Nz)
    function residual!(F, x_unk, _ignored)
        scatter!(phi_buf, x_unk, geom)
        laplacian!(lap_buf, phi_buf, geom)
        k = 0
        @inbounds for j in 1:Nz, i in 1:Nrho
            if geom.unk_mask[i, j]
                k += 1
                F[k] = p.kappa * lap_buf[i, j] -
                       (dfb_dphi(phi_buf[i, j], p) - mu_inf / p.nb)
            end
        end
        return nothing
    end
    return residual!
end

# Out-of-place residual mirroring the Python `_make_residual` for tests.
function make_residual(geom::Geometry, mu_inf::Float64)
    res! = make_residual_inplace(geom, mu_inf)
    n_unk = n_unknowns(geom)
    return function(x_unk)
        F = Vector{Float64}(undef, n_unk)
        res!(F, x_unk, nothing)
        return F
    end
end

# ------------------------------------------------------------------
# Initial guess
# ------------------------------------------------------------------

function initial_guess(geom::Geometry)
    p = geom.p
    R = p.R
    n = n_unknowns(geom)
    x0 = Vector{Float64}(undef, n)
    k = 0
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        if geom.unk_mask[i, j]
            k += 1
            rho_i = geom.rho_grid[i]
            z_j   = geom.z_grid[j]
            s1 = max(0.0, sqrt(rho_i^2 + (z_j - geom.z_c1)^2) - R)
            s2 = max(0.0, sqrt(rho_i^2 + (z_j - geom.z_c2)^2) - R)
            psi1, _, _ = psi_linear(s1, p)
            psi2, _, _ = psi_linear(s2, p)
            x0[k] = clip01(p.phi_inf + psi1 + psi2, p.eps)
        end
    end
    return x0
end

# ------------------------------------------------------------------
# Solver
# ------------------------------------------------------------------

# Constant linear part of the residual Jacobian: A = d/dx[ kappa * Lap(scatter(x)) ]
# restricted to the unknown cells.  Because phi = scatter(x) is affine in x and
# every unknown cell holds its own value directly, the full Jacobian is
#     J(x) = A - Diagonal(d2fb_dphi2.(x)),
# with A independent of x.  Assembling A once turns each Newton step into a
# single sparse LU instead of an unpreconditioned matrix-free GMRES solve
# (hundreds of residual evaluations per step).  Ghost cells are expanded into
# their clean P1/P2 corner cells with the same weights scatter! uses; those
# corners are guaranteed fluid (unknown or Dirichlet), so there is no recursion.
#
# When `sd` is supplied (bc_reg > 0), it holds d(bc)/d(phi_s) per ghost at the
# current iterate (see ghost_slope_derivs); the regularized surface term
# `ghost_coef_bc * bc(phi_s)` then also depends on the unknowns through
# phi_s(phi_P1, phi_P2), and that dependence is folded into the operator here.
function build_linear_operator(geom::Geometry;
                               sd::Union{Nothing, AbstractVector{Float64}} = nothing)
    p = geom.p
    Nrho, Nz = geom.Nrho, geom.Nz
    κ      = p.kappa
    drho2  = geom.drho^2
    dz2    = geom.dz^2

    unk_idx = zeros(Int, Nrho, Nz)
    k = 0
    @inbounds for j in 1:Nz, i in 1:Nrho
        if geom.unk_mask[i, j]
            k += 1
            unk_idx[i, j] = k
        end
    end
    n = k

    ghost_idx = zeros(Int, Nrho, Nz)
    @inbounds for g in 1:n_ghost(geom)
        ghost_idx[geom.ghost_i[g], geom.ghost_j[g]] = g
    end

    Irow = Int[]; Jcol = Int[]; Val = Float64[]

    # An unknown corner contributes directly; a Dirichlet corner is constant.
    add_unk! = function (row, ci, cj, w)
        (w == 0.0 || !geom.unk_mask[ci, cj]) && return
        push!(Irow, row); push!(Jcol, unk_idx[ci, cj]); push!(Val, w)
    end

    # Accumulate the dependence of phi[ci,cj] (scaled by coef) on the unknowns.
    add_phi! = function (row, ci, cj, coef)
        coef == 0.0 && return
        if geom.unk_mask[ci, cj]
            push!(Irow, row); push!(Jcol, unk_idx[ci, cj]); push!(Val, coef)
        elseif geom.dir_mask[ci, cj]
            return
        else
            g = ghost_idx[ci, cj]
            i01 = geom.ghost_i0_1[g]; j01 = geom.ghost_j0_1[g]
            c1  = coef * geom.ghost_coef_phi1[g]
            add_unk!(row, i01,     j01,     c1 * geom.ghost_wts_1[1, g])
            add_unk!(row, i01 + 1, j01,     c1 * geom.ghost_wts_1[2, g])
            add_unk!(row, i01,     j01 + 1, c1 * geom.ghost_wts_1[3, g])
            add_unk!(row, i01 + 1, j01 + 1, c1 * geom.ghost_wts_1[4, g])
            i02 = geom.ghost_i0_2[g]; j02 = geom.ghost_j0_2[g]
            c2  = coef * geom.ghost_coef_phi2[g]
            add_unk!(row, i02,     j02,     c2 * geom.ghost_wts_2[1, g])
            add_unk!(row, i02 + 1, j02,     c2 * geom.ghost_wts_2[2, g])
            add_unk!(row, i02,     j02 + 1, c2 * geom.ghost_wts_2[3, g])
            add_unk!(row, i02 + 1, j02 + 1, c2 * geom.ghost_wts_2[4, g])
            # bc_reg surface term: phi_G also depends on the unknowns through
            # phi_s(phi_P1, phi_P2) in the regularized Robin slope.  sd[g] is
            # d(bc)/d(phi_s) at the current iterate (0 for bc_reg = 0 or where
            # phi_s is clamped); this is the 2D analogue of the 1D
            # _surf_slope_deriv contribution to the surface node's Jacobian.
            if sd !== nothing && sd[g] != 0.0
                cbc = coef * geom.ghost_coef_bc[g] * sd[g]
                # phi_s = a1 * phi_P1 + a2 * phi_P2 (linear closure: a1=1,a2=0).
                if geom.ghost_used_linear[g]
                    a1 = 1.0; a2 = 0.0
                else
                    w  = (geom.ghost_delta[g] - geom.ghost_d1[g]) /
                         (geom.ghost_d2[g] - geom.ghost_d1[g])
                    a1 = 1.0 - w; a2 = w
                end
                cb1 = cbc * a1
                add_unk!(row, i01,     j01,     cb1 * geom.ghost_wts_1[1, g])
                add_unk!(row, i01 + 1, j01,     cb1 * geom.ghost_wts_1[2, g])
                add_unk!(row, i01,     j01 + 1, cb1 * geom.ghost_wts_1[3, g])
                add_unk!(row, i01 + 1, j01 + 1, cb1 * geom.ghost_wts_1[4, g])
                cb2 = cbc * a2
                add_unk!(row, i02,     j02,     cb2 * geom.ghost_wts_2[1, g])
                add_unk!(row, i02 + 1, j02,     cb2 * geom.ghost_wts_2[2, g])
                add_unk!(row, i02,     j02 + 1, cb2 * geom.ghost_wts_2[3, g])
                add_unk!(row, i02 + 1, j02 + 1, cb2 * geom.ghost_wts_2[4, g])
            end
        end
    end

    @inbounds for j in 1:Nz, i in 1:Nrho
        geom.unk_mask[i, j] || continue
        row = unk_idx[i, j]
        if i == 1
            # Axis stencil: 4(phi[2,j]-phi[1,j])/drho^2 + d2/dz2 (mirror ghost).
            add_phi!(row, 1, j,     κ * (-4.0 / drho2 - 2.0 / dz2))
            add_phi!(row, 2, j,     κ * (4.0 / drho2))
            add_phi!(row, 1, j + 1, κ / dz2)
            add_phi!(row, 1, j - 1, κ / dz2)
        else
            cross = 1.0 / (2.0 * geom.drho * geom.rho_grid[i])
            add_phi!(row, i,     j,     κ * (-2.0 / drho2 - 2.0 / dz2))
            add_phi!(row, i + 1, j,     κ * (1.0 / drho2 + cross))
            add_phi!(row, i - 1, j,     κ * (1.0 / drho2 - cross))
            add_phi!(row, i,     j + 1, κ / dz2)
            add_phi!(row, i,     j - 1, κ / dz2)
        end
    end

    return sparse(Irow, Jcol, Val, n, n)
end

# Per-ghost d(bc)/d(phi_s) at the current field phi -- the regularized Robin
# slope derivative (mirrors _surf_slope_deriv in bvp_solver.jl).  Zero when
# bc_reg = 0 or where the reconstructed surface value phi_s falls outside
# (eps, 1-eps) (there scatter! clamps phi_s, so the slope is flat).  phi_s is
# reconstructed from the two clean probe points exactly as in scatter!.
function ghost_slope_derivs(geom::Geometry, phi::AbstractMatrix)
    ng  = n_ghost(geom)
    sd  = zeros(Float64, ng)
    p   = geom.p
    reg = geom.bc_reg
    reg == 0.0 && return sd
    @inbounds for k in 1:ng
        i01 = geom.ghost_i0_1[k]; j01 = geom.ghost_j0_1[k]
        i02 = geom.ghost_i0_2[k]; j02 = geom.ghost_j0_2[k]
        phi_P1 = geom.ghost_wts_1[1, k] * phi[i01,     j01]     +
                 geom.ghost_wts_1[2, k] * phi[i01 + 1, j01]     +
                 geom.ghost_wts_1[3, k] * phi[i01,     j01 + 1] +
                 geom.ghost_wts_1[4, k] * phi[i01 + 1, j01 + 1]
        phi_P2 = geom.ghost_wts_2[1, k] * phi[i02,     j02]     +
                 geom.ghost_wts_2[2, k] * phi[i02 + 1, j02]     +
                 geom.ghost_wts_2[3, k] * phi[i02,     j02 + 1] +
                 geom.ghost_wts_2[4, k] * phi[i02 + 1, j02 + 1]
        d1 = geom.ghost_d1[k]; d2 = geom.ghost_d2[k]
        phi_s = geom.ghost_used_linear[k] ? phi_P1 :
                phi_P1 + (phi_P2 - phi_P1) * (geom.ghost_delta[k] - d1) / (d2 - d1)
        if p.eps < phi_s < 1.0 - p.eps
            sd[k] = reg / (p.kappa * phi_s * (1.0 - phi_s))
        end
    end
    return sd
end

# Fraction-to-boundary rule: the largest step length in (0, 1] for which
# x - alpha*delta stays strictly inside (lo, hi) COMPONENTWISE, pulled back by
# `tau` (< 1) so it stops just short of the boundary (the standard interior-point
# safeguard).  Returns 1.0 when the full Newton step x - delta is already feasible
# -- so well-behaved solves keep taking full steps and are numerically unchanged.
# Used by every Newton line search here and in bvp_arclength.jl so iterates never
# leave the physical box (0,1): f_b is only evaluated where it is the true free
# energy, and no non-smooth clip wall is ever created (paired with the C1 smooth
# extension of f_b in bvp_solver.jl, which handles any residual overshoot cleanly).
@inline function _max_feasible_step(x::AbstractVector, delta::AbstractVector,
                                    lo::Real, hi::Real; tau::Real = 0.995)
    alpha = 1.0
    @inbounds for k in eachindex(x)
        d = delta[k]
        if d > 0.0                       # x - alpha*d moves toward lo
            a = (x[k] - lo) / d
            a < alpha && (alpha = a)
        elseif d < 0.0                   # x - alpha*d moves toward hi (a > 0)
            a = (x[k] - hi) / d
            a < alpha && (alpha = a)
        end
    end
    return alpha >= 1.0 ? 1.0 : tau * alpha
end

# Newton with the analytic sparse Jacobian J = A - Diagonal(d2fb_dphi2.(x)) and a
# direct sparse LU per step.  Returns (x, converged).  Mirrors the abstol/maxiter
# semantics of _run_newton but is matrix-free-GMRES-free.
function _solve_newton_sparse(geom::Geometry, mu_inf::Float64,
                              x0::AbstractVector;
                              f_tol::Real = 1e-6, maxiter::Int = 200,
                              verbose::Bool = false)
    p = geom.p
    n = n_unknowns(geom)
    reg = geom.bc_reg
    A = build_linear_operator(geom)             # bc_reg-independent linear part
    res! = make_residual_inplace(geom, mu_inf)
    x = copy(x0)
    xtrial = Vector{Float64}(undef, n)
    F = Vector{Float64}(undef, n)
    diagv = Vector{Float64}(undef, n)
    phi_buf = reg == 0.0 ? nothing : Matrix{Float64}(undef, geom.Nrho, geom.Nz)
    res!(F, x, nothing)
    nrm = maximum(abs, F)
    it = 0
    for outer it in 1:maxiter
        nrm < f_tol && break
        @inbounds for k in 1:n
            diagv[k] = d2fb_dphi2(x[k], p)
        end
        # Linear part of the Jacobian.  For bc_reg = 0 it is the constant A;
        # for bc_reg > 0 the surface-BC term is iterate-dependent, so fold its
        # analytic derivative (sd per ghost) into the operator at the current x.
        Jlin = if reg == 0.0
            A
        else
            scatter!(phi_buf, x, geom)
            build_linear_operator(geom; sd = ghost_slope_derivs(geom, phi_buf))
        end
        δ = (Jlin - Diagonal(diagv)) \ F
        # Backtracking line search from the fraction-to-boundary cap: λ0 is the
        # largest step that keeps φ strictly inside (eps, 1-eps) componentwise, so
        # trials are feasible by construction and never clamped.  f_b is thus only
        # ever evaluated where it is the true FH free energy -- no clip wall.
        # λ0 = 1 (full Newton) whenever the full step is already feasible, so
        # smooth solves are unchanged; near saturation λ0 shrinks smoothly instead
        # of hitting a wall, and a stalled backtrack (no descent) still ends the
        # iteration cleanly.
        λ = _max_feasible_step(x, δ, p.eps, 1.0 - p.eps)
        accepted = false
        for _ in 1:30
            @inbounds @. xtrial = x - λ * δ
            res!(F, xtrial, nothing)
            newnrm = maximum(abs, F)
            if newnrm < nrm
                @inbounds @. x = xtrial
                nrm = newnrm
                accepted = true
                break
            end
            λ *= 0.5
        end
        accepted || break
    end
    verbose && println("  sparse Newton: $it iters, |F| = $nrm")
    return x, nrm < f_tol
end

function _run_newton(res!, x0::AbstractVector, f_tol::Real, maxiter::Int,
                     verbose::Bool)
    prob = NonlinearProblem(res!, x0)
    return solve(prob,
                 NewtonRaphson(autodiff = AutoFiniteDiff(),
                               linsolve = KrylovJL_GMRES());
                 abstol = f_tol, maxiters = maxiter,
                 show_trace = Val(verbose))
end

function solve_two_spheres(p::Parameters, d::Real;
                           rho_max = nothing, z_buf = nothing,
                           Nrho::Int = 80, Nz::Int = 200,
                           f_tol::Real = 1e-6,
                           maxiter::Int = 200,
                           verbose::Bool = false,
                           omega_continuation::Bool = true,
                           bc_reg::Real = 0.0,
                           x0_init = nothing)
    geom   = build_geometry(p, d; rho_max = rho_max, z_buf = z_buf,
                            Nrho = Nrho, Nz = Nz, bc_reg = bc_reg)
    mu_inf = compute_mu_inf(p)
    # Optional warm start (e.g. a neighbouring-d solution interpolated onto this
    # grid); must already match n_unknowns(geom).  Falls back to the analytic guess.
    x0     = x0_init === nothing ? initial_guess(geom) : x0_init

    x_sol, success = try
        _solve_newton_sparse(geom, mu_inf, x0;
                             f_tol = f_tol, maxiter = maxiter, verbose = verbose)
    catch e
        verbose && @info "primary solve threw" exception = e
        (copy(x0), false)
    end

    if !success && omega_continuation && p.omega_b != 0.0
        verbose && println("Cold start failed; trying omega_b continuation ...")
        n_steps = 6
        p0   = with(p; omega_b = 0.0)
        geom0 = build_geometry(p0, d; rho_max = rho_max, z_buf = z_buf,
                               Nrho = Nrho, Nz = Nz, bc_reg = bc_reg)
        x_cur = initial_guess(geom0)
        geom_k = geom
        ok_k = true
        for k in 1:n_steps
            frac = k / n_steps
            p_k = with(p; omega_b = frac * p.omega_b)
            geom_k = build_geometry(p_k, d; rho_max = rho_max, z_buf = z_buf,
                                    Nrho = Nrho, Nz = Nz, bc_reg = bc_reg)
            x_k, ok_step = try
                _solve_newton_sparse(geom_k, compute_mu_inf(p_k), x_cur;
                                     f_tol = f_tol, maxiter = maxiter)
            catch e
                (x_cur, false)
            end
            x_cur = x_k
            if ok_step
                verbose && println("  continuation step $k/$n_steps (frac=$(round(frac, digits=2))) ok")
            else
                verbose && println("  continuation step $k/$n_steps (frac=$(round(frac, digits=2))) did not converge")
                ok_k = false
                break
            end
        end
        x_sol  = x_cur
        geom   = geom_k
        res_final! = make_residual_inplace(geom, compute_mu_inf(p))
        Fcheck = Vector{Float64}(undef, n_unknowns(geom))
        res_final!(Fcheck, x_sol, nothing)
        success = ok_k && (maximum(abs, Fcheck) < 10.0 * f_tol)
    end

    phi = scatter(x_sol, geom)
    return Dict{Symbol, Any}(
        :phi     => phi,
        :geom    => geom,
        :p       => p,
        :d       => Float64(d),
        :mu_inf  => mu_inf,
        :success => success,
        :x       => x_sol,
    )
end

# ------------------------------------------------------------------
# Bilinear sampling
# ------------------------------------------------------------------

function interp_phi(phi::AbstractMatrix, geom::Geometry,
                    rho_q::AbstractVector, z_q::AbstractVector)
    @assert length(rho_q) == length(z_q)
    rho_g = geom.rho_grid
    z_g   = geom.z_grid
    drho, dz = geom.drho, geom.dz
    Nrho, Nz = geom.Nrho, geom.Nz
    n = length(rho_q)
    out = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        rc = clamp(rho_q[k], rho_g[1], rho_g[end] - 1e-12 * drho)
        zc = clamp(z_q[k],   z_g[1],   z_g[end]   - 1e-12 * dz)
        ti = (rc - rho_g[1]) / drho
        tj = (zc - z_g[1])   / dz
        i0 = clamp(Int(floor(ti)) + 1, 1, Nrho - 1)
        j0 = clamp(Int(floor(tj)) + 1, 1, Nz   - 1)
        t_rho = ti - (i0 - 1)
        t_z   = tj - (j0 - 1)
        out[k] = (1.0 - t_rho) * (1.0 - t_z) * phi[i0,     j0]     +
                 t_rho         * (1.0 - t_z) * phi[i0 + 1, j0]     +
                 (1.0 - t_rho) * t_z         * phi[i0,     j0 + 1] +
                 t_rho         * t_z         * phi[i0 + 1, j0 + 1]
    end
    return out
end

function interp_phi(phi::AbstractMatrix, geom::Geometry, rho_q::Real, z_q::Real)
    return interp_phi(phi, geom, [Float64(rho_q)], [Float64(z_q)])[1]
end

# ------------------------------------------------------------------
# Excess free energy
# ------------------------------------------------------------------

function cell_fluid_fraction(geom::Geometry; nsub::Int = 8)
    drho, dz = geom.drho, geom.dz
    R = geom.p.R
    R2 = R * R
    offsets = ((1:nsub) .- 0.5) ./ nsub .- 0.5   # length nsub in (-0.5, 0.5)
    counts = zeros(Int, geom.Nrho, geom.Nz)
    for a in offsets, b in offsets
        @inbounds for j in 1:geom.Nz
            z_s = geom.z_grid[j] + b * dz
            dz1 = z_s - geom.z_c1
            dz2 = z_s - geom.z_c2
            for i in 1:geom.Nrho
                rho_s = max(geom.rho_grid[i] + a * drho, 0.0)
                r1sq = rho_s * rho_s + dz1 * dz1
                r2sq = rho_s * rho_s + dz2 * dz2
                if (r1sq < R2) || (r2sq < R2)
                    counts[i, j] += 1
                end
            end
        end
    end
    return 1.0 .- counts ./ (nsub * nsub)
end

function phi_on_sphere(phi::AbstractMatrix, geom::Geometry,
                       zc::Real, theta::AbstractVector;
                       probe_steps::Int = 80)
    p = geom.p
    R = p.R
    drho, dz = geom.drho, geom.dz
    step = max(drho, dz)
    omega_kappa = p.omega_b / p.kappa
    Nrho, Nz = geom.Nrho, geom.Nz
    rho_g, z_g = geom.rho_grid, geom.z_grid
    n = length(theta)

    sin_t = sin.(theta)
    cos_t = cos.(theta)
    rho_s = R .* sin_t
    z_s   = zc .+ R .* cos_t

    d1 = fill(step, n)
    d2 = fill(2.0 * step, n)

    bad1 = falses(n)
    bad2 = falses(n)

    # Probe cap per theta: entry of the outward-normal ray into the OTHER
    # sphere (see build_geometry).  Thetas whose probes cannot find clean
    # stencils within the cap are unresolvable at this grid (sub-resolution
    # contact cap); their phi_s is filled from the nearest resolvable theta
    # after the fit -- the adsorbed film there is slaved to the shoulder.
    zo = abs(zc - geom.z_c1) < abs(zc - geom.z_c2) ? geom.z_c2 : geom.z_c1
    cap = [_ray_sphere_entry(rho_s[k], z_s[k], sin_t[k], cos_t[k], zo, R)
           for k in 1:n]
    unresolved = falses(n)

    function stencil_dirty!(d_arr::Vector{Float64}, bad_out::BitVector)
        @inbounds for k in 1:n
            rho_P = rho_s[k] + d_arr[k] * sin_t[k]
            z_P   = z_s[k]   + d_arr[k] * cos_t[k]
            rho_Pc = clamp(rho_P, rho_g[1], rho_g[end] - 1e-12 * drho)
            z_Pc   = clamp(z_P,   z_g[1],   z_g[end]   - 1e-12 * dz)
            ti = (rho_Pc - rho_g[1]) / drho
            tj = (z_Pc   - z_g[1])   / dz
            i0 = clamp(Int(floor(ti)) + 1, 1, Nrho - 1)
            j0 = clamp(Int(floor(tj)) + 1, 1, Nz   - 1)
            bad_out[k] = geom.inside_mask[i0,     j0]     ||
                         geom.inside_mask[i0 + 1, j0]     ||
                         geom.inside_mask[i0,     j0 + 1] ||
                         geom.inside_mask[i0 + 1, j0 + 1]
        end
    end

    for _ in 1:probe_steps
        stencil_dirty!(d1, bad1)
        stencil_dirty!(d2, bad2)
        @inbounds for k in 1:n
            if (bad1[k] && d1[k] + step > cap[k]) ||
               (bad2[k] && d2[k] + step > cap[k])
                unresolved[k] = true
            end
        end
        done = true
        @inbounds for k in 1:n
            unresolved[k] && continue
            if bad1[k] || bad2[k]
                done = false
                bad1[k] && (d1[k] += step)
                bad2[k] && (d2[k] += step)
                d2[k] = max(d2[k], d1[k] + step)
            end
        end
        done && break
    end
    stencil_dirty!(d1, bad1)
    stencil_dirty!(d2, bad2)
    unresolved .|= bad1 .| bad2

    rho_P1 = rho_s .+ d1 .* sin_t
    z_P1   = z_s   .+ d1 .* cos_t
    rho_P2 = rho_s .+ d2 .* sin_t
    z_P2   = z_s   .+ d2 .* cos_t
    phi_P1 = interp_phi(phi, geom, rho_P1, z_P1)
    phi_P2 = interp_phi(phi, geom, rho_P2, z_P2)

    # Surface value from the quadratic through phi(P1), phi(P2) that carries the
    # Robin slope at s = 0.  Eliminating the curvature gives phi_s = A - B * slope
    # with B = d1 d2 / (d1 + d2) > 0.  At bc_reg > 0 the slope is itself a
    # function of phi_s -- the same regularized BC scatter! applies -- so that
    # relation is implicit; A - B * slope(phi_s) decreases in phi_s, so the root
    # is unique and bisection finds it.  Using the bare omega_b/kappa there (as
    # this did before the regularized branch existed) reports a phi_s that does
    # not satisfy the BC the solve imposed, which shows up directly as a violated
    # envelope identity for the free energy.
    reg_kappa = geom.bc_reg / p.kappa
    out = similar(phi_P1)
    @inbounds for k in 1:n
        B_k = d1[k] * d2[k] / (d1[k] + d2[k])
        A_k = phi_P1[k] -
              d1[k]^2 * (phi_P1[k] - phi_P2[k]) / (d1[k]^2 - d2[k]^2)
        if reg_kappa == 0.0
            out[k] = A_k - B_k * omega_kappa
        else
            lo, hi = p.eps, 1.0 - p.eps
            for _ in 1:60
                mid = 0.5 * (lo + hi)
                g = mid - A_k +
                    B_k * (omega_kappa + reg_kappa * log(mid / (1.0 - mid)))
                g < 0.0 ? (lo = mid) : (hi = mid)
            end
            out[k] = 0.5 * (lo + hi)
        end
    end
    if any(unresolved) && !all(unresolved)
        res_idx = findall(.!unresolved)
        @inbounds for k in findall(unresolved)
            out[k] = out[res_idx[argmin(abs.(res_idx .- k))]]
        end
    end
    # A surface volume fraction is physically in [0,1]; the quadratic extrapolation
    # overshoots (phi_s ~ 5) where the probe stencil degenerates near contact, which
    # both reports nonsense and inflates the surface free-energy term in
    # excess_free_energy_two_spheres.  Clamp to the physical range.
    return clamp.(out, 0.0, 1.0)
end

# Simpson's 1/3 rule on a uniform grid, length(y) must be odd.
function _simpson(y::AbstractVector, x::AbstractVector)
    n = length(x)
    @assert n == length(y)
    @assert isodd(n)
    h = x[2] - x[1]
    s = y[1] + y[n]
    @inbounds for k in 2:n - 1
        s += (iseven(k) ? 4 : 2) * y[k]
    end
    return s * h / 3.0
end

# `frac` optionally supplies a precomputed cell_fluid_fraction(geom): it depends
# only on the geometry, not on phi or phi_inf, so a caller evaluating F many
# times on one mesh (an arclength sweep) should hoist it out of the loop --
# nsub^2 = 64 in-sphere tests per cell is otherwise the dominant cost here.
function excess_free_energy_two_spheres(result::AbstractDict;
                                        nsub::Int = 8,
                                        n_theta::Int = 1001,
                                        frac::Union{Nothing, AbstractMatrix} = nothing)
    geom = result[:geom]
    phi  = result[:phi]
    p    = geom.p
    phi_inf = p.phi_inf

    dfb_inf = dfb_dphi(phi_inf, p)
    f_inf   = fb(phi_inf, p)

    integrand = Matrix{Float64}(undef, geom.Nrho, geom.Nz)
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        if geom.fluid_mask[i, j]
            phi_ij = phi[i, j]
            psi_ij = phi_ij - phi_inf
            df_loc = fb(phi_ij, p) - f_inf - dfb_inf * psi_ij
            dfp    = dfb_dphi(phi_ij, p) - dfb_inf
            integrand[i, j] = df_loc - 0.5 * psi_ij * dfp
        else
            integrand[i, j] = 0.0
        end
    end

    frac === nothing && (frac = cell_fluid_fraction(geom; nsub = nsub))
    bulk_sum = 0.0
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        bulk_sum += integrand[i, j] * (2.0 * pi * geom.rho_grid[i]) * frac[i, j]
    end
    bulk = bulk_sum * geom.drho * geom.dz

    n_t = isodd(n_theta) ? n_theta : n_theta + 1
    theta = collect(range(0.0, pi; length = n_t))
    sin_t = sin.(theta)

    # Surface contribution: the surface free energy itself plus the boundary
    # term left over from integrating the gradient energy by parts.  Writing the
    # Robin BC as kappa dphi/ds|_S = fs'(phi_s) (s outward from the surface),
    #
    #     fs(phi) = omega_b * phi + bc_reg * [phi log phi + (1-phi) log(1-phi)],
    #
    # the two pieces are Int[fs(phi_s) - fs(phi_inf)]dA (referenced to a surface
    # bathed in bulk, as the volume integrand is) and -1/2 Int (phi_s - phi_inf)
    # fs'(phi_s) dA.  At bc_reg = 0 they collapse to 1/2 omega_b Int (phi_s -
    # phi_inf) dA, the form this function used before the regularized BC existed,
    # so nothing changes for the unregularized pipeline.
    #
    # The bc_reg piece is NOT optional at bc_reg > 0: it is branch-dependent (the
    # thin film and the capillary bridge sit at very different phi_s), and
    # omitting it inverts the ordering of the two stable branches.  See the
    # envelope-identity testset in test/test_bvp_arclength.jl, which checks
    # dF/d(omega_b) = Int (phi_s - phi_inf) dA at both bc_reg = 0 and 0.05.
    reg  = geom.bc_reg
    _hs(u)   = (v = clip01(u, p.eps); v * log(v) + (1.0 - v) * log(1.0 - v))
    _fs(u)   = p.omega_b * u + reg * _hs(u)
    _dfs(u)  = (v = clip01(u, p.eps); p.omega_b + reg * log(v / (1.0 - v)))
    fs_inf = _fs(phi_inf)

    surf_total = 0.0
    for zc in (geom.z_c1, geom.z_c2)
        phi_s = phi_on_sphere(phi, geom, zc, theta)
        surf_total += _simpson((_fs.(phi_s) .- fs_inf .-
                                0.5 .* (phi_s .- phi_inf) .* _dfs.(phi_s)) .* sin_t,
                               theta)
    end
    surf = 2.0 * pi * p.R^2 * surf_total
    return bulk + surf
end
