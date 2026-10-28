# Fisher information rate of the binding/unbinding first-passage times vs the bulk
# volume fraction phi_inf.
#
# JOB SCRIPT (spec/jobs/20260720/): refactor of viz_util/fisher_rate_from_fpt_solver.jl
# that (1) takes every parameter from its local param.jl (`PARAMS`) and (2) consumes
# the pre-built potentials from the bvp_potential_interp stage instead of rebuilding
# them in-place.  It reads
#     PARAMS.sim.output_dir / PARAMS.sim.interp_file        (bvp_potential_interp.json)
# for the x-grid and the total potential U_tot(x) per phi_inf (U_tot_full for
# PARAMS.fpt.potential == :bvp, U_tot_linear for :linear), then for every phi_inf runs
# the backward-FPE solver (src/fpt_solver.jl) twice on the reduced 1D end-to-end
# coordinate x = surface gap d + R:
#
#   1. escape ("out of the well"): start at contact x = R (reflecting), absorbing at
#      the barrier ridge x_ridge  -> time to climb OUT and reach the ridge.
#   2. return ("into the well"):   start at the ridge (reflecting far end), absorbing
#      at x = R  -> time to fall back to contact.
#
# The stored x-grid IS the FPT return domain (range(R, 4 a sqrt(s), Ngrid)); the
# return problem uses it directly.  The escape domain [R, x_ridge] is a sub-interval,
# so U_tot is linearly resampled onto range(R, x_ridge, Ngrid); the ridge is located
# sub-grid by a parabola on the stored tabulation (as the original did on its fine
# scan) so x_ridge moves smoothly with phi and the d^2/dphi^2 stencil stays clean.
#
# The Fisher information rate of the exit/return times wrt phi is the base-PDF average of
# the SQUARED central first derivative of log P wrt phi (the first-derivative-squared form
# of I(phi), numerically stabler than differencing second derivatives), summed over
# escape+return and divided by (MFPT_esc + MFPT_ret).  Everything (both densities,
# their MFPTs, the potentials on the grids, the local concentration and its phi
# derivatives, the Fisher rate, and to_dict(PARAMS)) is written to
# PARAMS.sim.output_dir / PARAMS.sim.fisher_file.  No paths/filenames are hardcoded.
#
#   julia --project=@. spec/jobs/20260720/fisher_rate_from_fpt_solver.jl

using Droplet
using Printf
import JSON

include(joinpath(@__DIR__, "param.jl"))

# ── Load the bvp_potential_interp stage output ───────────────────────────────
interp_path = joinpath(PARAMS.sim.output_dir, PARAMS.sim.interp_file)
isfile(interp_path) ||
    error("missing interp output: $interp_path\n" *
          "run spec/jobs/20260720/bvp_potential_interp.jl first")
_ip = JSON.parsefile(interp_path)

_f(x)    = x === nothing ? NaN : Float64(x)              # JSON null -> NaN
_farr(v) = Float64[_f(x) for x in v]
_fmat(m) = [_farr(row) for row in m]

potential  = PARAMS.fpt.potential                       # :bvp or :linear (a Symbol)
phi_infs   = _farr(_ip["phi_inf"])
xg_full    = _farr(_ip["x_grid"])                       # == range(R, 4 a sqrt(s), Ngrid)
U_tot_rows = potential == :bvp ? _fmat(_ip["U_tot_full"]) : _fmat(_ip["U_tot_linear"])
# BVP gap-facing surface concentration on x_grid; averaged by local_conc under :bvp.
phi_surf_rows = _fmat(_ip["phi_surf_interp"])

# The interp stage already enforced this; re-check so the two stages can't silently
# drift apart.
let g = phi_inf_grid(PARAMS)
    (length(phi_infs) == length(g) && all(isapprox.(phi_infs, g; atol = 1e-9, rtol = 1e-6))) ||
        error("interp json phi_inf != phi_inf_grid(PARAMS):\n  got $phi_infs\n  exp $g")
end

# ── Fixed parameters (from PARAMS) ───────────────────────────────────────────
R      = PARAMS.surface.R
h0     = PARAMS.surface.omega_b        # adsorption strength (linear-response profile)
s, a   = PARAMS.rouse.s, PARAMS.rouse.a
kBT    = PARAMS.rouse.kBT
Deff   = D_eff(PARAMS.rouse)           # 2 kBT / friction (relative end-to-end coord)
Ngrid  = PARAMS.fpt.Ngrid
Nlocal = PARAMS.fpt.Nlocal
ridge_max = PARAMS.fpt.ridge_scan_max

length(xg_full) == Ngrid ||
    error("interp x_grid has $(length(xg_full)) points, PARAMS.fpt.Ngrid = $Ngrid")
x_max = xg_full[end]

# ── Potential helpers: consume the tabulated U_tot(x) ────────────────────────
# Linear interpolation of a tabulated potential onto an arbitrary x (the stored grid
# is fine, so this is accurate; flat beyond the ends).
function _interp(xg, y, xq)
    xq <= xg[1]   && return y[1]
    xq >= xg[end] && return y[end]
    j = searchsortedlast(xg, xq)
    t = (xq - xg[j]) / (xg[j+1] - xg[j])
    return y[j] * (1 - t) + y[j+1] * t
end

# Barrier ridge = first interior local max of U (up to ridge_max), parabola-refined to
# sub-grid accuracy so it tracks phi smoothly instead of snapping to the grid.
function _find_ridge(xg, U, xmax)
    n  = something(findlast(<=(xmax), xg), length(xg))
    ir = findfirst(i -> U[i] > U[i-1] && U[i] > U[i+1], 2:n-1)
    i  = ir === nothing ? argmax(@view U[1:n]) : ir + 1
    (1 < i < n) || return xg[i]
    A, B, C = U[i-1], U[i], U[i+1]
    return xg[i] + 0.5 * (xg[2] - xg[1]) * (A - C) / (A - 2B + C)
end

# ── Solve the escape and return FPT problems at one phi_inf ───────────────────
# `Urow` is the stored U_tot(x) on xg_full for this phi_inf.  If t_esc/t_ret are
# supplied the PDFs are evaluated on exactly those grids (so neighbouring phi are
# comparable point-by-point); empty => the solver's own auto log-spaced grid.
function solve_phi(phi, Urow; t_esc = Float64[], t_ret = Float64[])
    p   = Droplet.with(PARAMS.surface; phi_inf = phi)
    lam = sqrt(d2fb_dphi2(phi, p) / p.kappa)              # screening length (reported)

    x_ridge = _find_ridge(xg_full, Urow, ridge_max)
    U_R     = Urow[1]                                     # U_tot at contact (xg_full[1] == R)
    barrier = _interp(xg_full, Urow, x_ridge) - U_R

    # 1. escape: R -> ridge.  Resample U_tot onto [R, x_ridge] (sub-interval of the
    #    stored grid); reflecting at R, absorbing at the ridge.
    xg_esc = collect(range(R, x_ridge; length = Ngrid))
    Ug_esc = [_interp(xg_full, Urow, x) for x in xg_esc]
    p0_esc = zeros(Ngrid); p0_esc[1] = 1.0
    res_esc = solve_fpt(FPTProblem(xg_esc, Ug_esc; D = Deff, kBT = kBT, t = t_esc,
                                   left_bc = :reflecting, right_bc = :absorbing, p0 = p0_esc))

    # 2. return: ridge -> R, on the stored grid itself (xg_ret == xg_full).  Spread the
    #    start delta over the two nodes bracketing x_ridge (weighted by distance) so the
    #    return IC moves continuously with x_ridge.
    xg_ret = collect(range(R, x_max; length = Ngrid))
    Ug_ret = copy(Urow)
    j = clamp(searchsortedlast(xg_ret, x_ridge), 1, Ngrid - 1)
    w = (x_ridge - xg_ret[j]) / (xg_ret[2] - xg_ret[1])
    p0_ret = zeros(Ngrid); p0_ret[j] = 1 - w; p0_ret[j+1] = w
    res_ret = solve_fpt(FPTProblem(xg_ret, Ug_ret; D = Deff, kBT = kBT, t = t_ret,
                                   left_bc = :absorbing, right_bc = :reflecting, p0 = p0_ret))

    return (lam = lam, x_ridge = x_ridge, barrier = barrier, well_depth = U_R,
            xg_esc = xg_esc, Ug_esc = Ug_esc, xg_ret = xg_ret, Ug_ret = Ug_ret,
            res_esc = res_esc, res_ret = res_ret)
end

# ∫ P(t|phi) (∂logP(t|phi)/∂phi)² dt — the first-derivative-squared form of the Fisher
# information.  It equals -∫P ∂²logP/∂phi² dt but is numerically stabler: the integrand
# is non-negative by construction and never differences large curvatures against each
# other.  The score ∂logP/∂phi is the central 2-point derivative on the outer pair
# (phi, phinn) (spacing 2 dφ, so it is centred on phin like the old 3-pt stencil was); the
# base PDF `p0` is the integration weight.  Non-uniform-t trapezoid; tail nodes where any
# of the three PDFs underflows (or goes ≤ 0) are dropped.
function avg_dlogp2(t, p0, pn, pnn, dφ)
    fl   = maximum(p0) * 1e-10
    lp0  = log.(max.(p0,  fl))
    lpnn = log.(max.(pnn, fl))
    integrand = p0 .* ((lpnn .- lp0) ./ (2 * dφ)).^2
    @inbounds for i in eachindex(integrand)
        (p0[i] > fl && pn[i] > fl && pnn[i] > fl) || (integrand[i] = 0.0)
    end
    s = 0.0
    @inbounds for i in 1:length(t)-1
        s += 0.5 * (integrand[i] + integrand[i+1]) * (t[i+1] - t[i])
    end
    return s
end

# ── Local concentration in the adsorption well ───────────────────────────────
# Boltzmann-weighted average of a volume-fraction profile over the well [R, x_ridge]
# (barrier top), with the weight w(x) = exp(-U_tot(x)/kBT) from the SELECTED potential.
# The profile must come from the same model as the weight:
#   :bvp    -> the interpolated BVP gap-facing surface concentration phi_surf(x-R);
#   :linear -> the linear-response profile
#              phi(x) = phi_inf - h0 exp(-lam(x-R)) R/x / (kappa (lam + 1/R)).
# (Mixing the BVP weight with the linear profile, as the pre-refactor code did, was
# inconsistent.)  A dedicated dense grid keeps the trapezoid tight.
_trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

function local_conc(phi, Urow, PsRow)
    all(isfinite, Urow) || return NaN
    x_ridge = _find_ridge(xg_full, Urow, ridge_max)
    xg = collect(range(R, x_ridge; length = Nlocal))
    Ug = [_interp(xg_full, Urow, x) for x in xg]
    wt = exp.(-(Ug .- minimum(Ug)) ./ kBT)                # shift for overflow safety (cancels)
    if potential == :bvp
        all(isfinite, PsRow) || return NaN
        profile = [_interp(xg_full, PsRow, x) for x in xg] # BVP surface conc(x-R)
    else
        p   = Droplet.with(PARAMS.surface; phi_inf = phi)
        lam = sqrt(d2fb_dphi2(phi, p) / p.kappa)
        profile = phi .- h0 .* (exp.(-lam .* (xg .- R)) .* R ./ xg) ./ (p.kappa * (lam + 1 / R))
    end
    return _trapz(xg, wt .* profile) / _trapz(xg, wt)
end

# np.gradient-style central differences (÷ the true phi spacing) for d/dphi_inf.
function _gradient(y, h)
    n = length(y); g = similar(y)
    n == 1 && (g[1] = NaN; return g)
    g[1] = (y[2] - y[1]) / h
    g[n] = (y[n] - y[n-1]) / h
    for i in 2:n-1
        g[i] = (y[i+1] - y[i-1]) / (2h)
    end
    return g
end

# ── Sweep ────────────────────────────────────────────────────────────────────
nphi        = length(phi_infs)
records     = Vector{Dict{String,Any}}()
fisher_rate = fill(NaN, max(nphi - 2, 0))               # indexed by k = index of phi
fisher_rate_esc = fill(NaN, max(nphi - 2, 0))               # indexed by k = index of phi
fisher_rate_ret = fill(NaN, max(nphi - 2, 0))               # indexed by k = index of phi

# local concentration + its 1st/2nd derivatives wrt phi_inf across the whole sweep
phi_local_vals = [local_conc(phi_infs[k], U_tot_rows[k], phi_surf_rows[k]) for k in 1:nphi]
if nphi >= 2
    dφ_grid          = phi_infs[2] - phi_infs[1]        # uniform phi spacing
    dphi_local_vals  = _gradient(phi_local_vals, dφ_grid)
    d2phi_local_vals = _gradient(dphi_local_vals, dφ_grid)
else
    dphi_local_vals  = fill(NaN, nphi)
    d2phi_local_vals = fill(NaN, nphi)
end

for k in 1:max(nphi - 2, 0)
    φ, φn, φnn = phi_infs[k], phi_infs[k+1], phi_infs[k+2]
    dφ = φn - φ                                          # uniform phi spacing

    # :bvp rows can be all-NaN if a phi_inf had < 2 finite BVP nodes; skip such stencils.
    if any(r -> !all(isfinite, r), (U_tot_rows[k], U_tot_rows[k+1], U_tot_rows[k+2]))
        @printf("phi=%.4f  SKIPPED (non-finite U_tot in stencil)\n", φ)
        continue
    end

    # 1-2. base phi: solve on the auto grids, then freeze them for the neighbours
    base = solve_phi(φ, U_tot_rows[k])
    t_esc_grid = base.res_esc.t
    t_ret_grid = base.res_ret.t

    # 3. neighbours on the SAME time grids (PDFs comparable point-by-point)
    nb  = solve_phi(φn,  U_tot_rows[k+1]; t_esc = t_esc_grid, t_ret = t_ret_grid)
    nnb = solve_phi(φnn, U_tot_rows[k+2]; t_esc = t_esc_grid, t_ret = t_ret_grid)

    # 4-5. squared central 1st derivative of log P wrt phi, averaged over the base PDF
    sc2_esc = avg_dlogp2(t_esc_grid, base.res_esc.pdf, nb.res_esc.pdf, nnb.res_esc.pdf, dφ)
    sc2_ret = avg_dlogp2(t_ret_grid, base.res_ret.pdf, nb.res_ret.pdf, nnb.res_ret.pdf, dφ)

    # 6. summed averaged squared scores per unit (MFPT_esc + MFPT_ret); I(phi) ≥ 0 (no −)
    fisher_rate_esc[k] = φ^2 * sc2_esc / (base.res_esc.mfpt + base.res_ret.mfpt)
    fisher_rate_ret[k] = φ^2 * sc2_ret / (base.res_esc.mfpt + base.res_ret.mfpt)
    fisher_rate = fisher_rate_esc .+ fisher_rate_ret

    push!(records, Dict(
        "phi_inf"     => φ,
        "lam"         => base.lam,
        "x_ridge"     => base.x_ridge,
        "barrier"     => base.barrier,
        "well_depth"  => base.well_depth,
        "x_grid_esc"  => base.xg_esc,
        "U_tot_esc"   => base.Ug_esc,
        "x_grid_ret"  => base.xg_ret,
        "U_tot_ret"   => base.Ug_ret,
        "escape_t"    => base.res_esc.t,
        "escape_pdf"  => base.res_esc.pdf,
        "escape_mfpt" => base.res_esc.mfpt,
        "return_t"    => base.res_ret.t,
        "return_pdf"  => base.res_ret.pdf,
        "return_mfpt" => base.res_ret.mfpt,
        "avg_dlogp2_esc" => sc2_esc,
        "avg_dlogp2_ret" => sc2_ret,
        "phi_surf"    => phi_surf_rows[k],
        "phi_local"   => phi_local_vals[k],
        "dphi_local"  => dphi_local_vals[k],
        "d2phi_local" => d2phi_local_vals[k],
        "fisher_rate_esc" => fisher_rate_esc[k],
        "fisher_rate_ret" => fisher_rate_ret[k],
        "fisher_rate" => fisher_rate[k]))

    @printf("phi=%.4f  lam=%5.2f  ridge=%.3f  barrier=%.2f kBT  phi_loc=%.4g  MFPT_esc=%.4g  MFPT_ret=%.4g  I_rate=%.4g\n",
            φ, base.lam, base.x_ridge, base.barrier, phi_local_vals[k],
            base.res_esc.mfpt, base.res_ret.mfpt, fisher_rate[k])
end

# ── Dump everything beside the input data ────────────────────────────────────
outpath = joinpath(PARAMS.sim.output_dir, PARAMS.sim.fisher_file)
open(outpath, "w") do io
    JSON.print(io, Dict(
        "params"       => to_dict(PARAMS),
        "potential"    => String(potential),
        "interp_input" => PARAMS.sim.interp_file,
        "meta" => Dict(
            "escape_setup" => "IC=delta at R; reflecting at R, absorbing at ridge (time R->ridge)",
            "return_setup" => "IC=delta at ridge; absorbing at R, reflecting far end (time ridge->R)",
            "fisher_rate_def" => "phi^2 (A_esc+A_ret)/(MFPT_esc+MFPT_ret), " *
                                 "A=∫P(t|phi) (∂logP/∂phi)² dt (squared central 1st deriv in phi)",
            "phi_local_def" => "Boltzmann-weighted (w=exp(-U_tot/kBT)) avg over the well " *
                               "[R,x_ridge] of the profile from the SELECTED model: :bvp uses " *
                               "the interpolated BVP surface conc phi_surf(x-R); :linear uses " *
                               "phi_inf - h0 exp(-lam(x-R)) R/x /(kappa(lam+1/R)). " *
                               "dphi/d2phi are central FD wrt phi_inf"),
        "phi_infs"    => phi_infs,
        "fisher_rate" => fisher_rate,                   # length nphi-2, indexed by phi-index
        "phi_local"             => phi_local_vals,      # length nphi
        "dphi_local_dphi_inf"   => dphi_local_vals,
        "d2phi_local_dphi_inf2" => d2phi_local_vals,
        "per_phi"     => records))
end
println("potential=", potential, "  nphi=", nphi, "  records=", length(records))
println("saved -> ", outpath)
