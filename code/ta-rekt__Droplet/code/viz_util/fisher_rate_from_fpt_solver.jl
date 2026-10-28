# Concentration sweep of the first-passage-time distributions in the combined
# attractive-adsorption + Rouse-entropy potential (the reduced 1D end-to-end
# coordinate x, cf. fpt_solver_validation_rouse.jl).
#
# The bulk volume fraction phi_inf is the only swept parameter: it sets the
# screening length lam = sqrt(d2fb_dphi2(phi_inf) / kappa), which controls the
# range/depth of the adsorption well U_ads.  Higher phi_inf => smaller lam =>
# deeper, longer-ranged well => higher escape barrier.
#
# For every phi_inf the backward-FPE solver (src/fpt_solver.jl) is run TWICE on
# the domain [R, x_ridge]:
#
#   1. escape  ("out of the well"): start at contact x = R (reflecting wall),
#      absorbing at the ridge  -> time to climb OUT of the well and reach x_ridge.
#   2. return  ("into the well"):   start at the ridge on the right (reflecting
#      wall), absorbing at x = R  -> time to reach x = R starting from the right.
#
# On top of the sweep this script also estimates the FISHER INFORMATION RATE of
# the exit/return times with respect to phi.  At each phi_inf[k] we:
#   * solve the escape/return PDFs on the solver's auto time grids and freeze
#     those grids (t_esc_grid, t_ret_grid);
#   * re-solve at phi_inf[k+1] and phi_inf[k+2] on the SAME grids so the PDFs are
#     directly comparable time-point by time-point;
#   * form the central 3-point second derivative of log P wrt phi and average it
#     over the base PDF,  A = ∫_0^{Tmax} dt P(t|phi) ∂²logP(t|phi)/∂phi²
#     (this equals -I(phi), the negative Fisher information, since ∫P dt = 1);
#   * report (A_exit + A_return) / (MFPT_exit + MFPT_return): the averaged
#     second derivative per unit time.  It is stored in `fisher_rate`, indexed by
#     the index k of phi (the last two entries are absent — the stencil needs
#     phi[k+1], phi[k+2]).
#
# Both densities, their means (MFPTs), the potential on the grid, the Fisher rate,
# and every fixed/varied parameter are dumped to data/<today>/fpt_solver_c_sweep/.
#
#   julia --project=@. viz_util/fpt_solver_c_sweep.jl

ENV["GKSwstype"] = "100"

using Droplet
using Plots
using Statistics
using Printf
using Dates
import JSON

include(joinpath(@__DIR__, "bvp_potential_interp.jl"))

# --- Attractive-potential source ---------------------------------------------
# :linear -> closed-form linear-response adsorption well (make_U_tot, below).
# :bvp    -> the nonlinear two-sphere BVP U_full(d) (from dashboard_fisher_two_
#            spheres.py), interpolated to high resolution by an isotonic + monotone
#            PCHIP scheme (bvp_potential_interp.jl) and mapped onto the FPT
#            coordinate x = d + R.  build_U_tot adds the SAME radial Jacobian and
#            Rouse spring on top; only the attractive term switches.
const POTENTIAL = :bvp
const BVP_JSON  = joinpath(@__DIR__, "..", "media", "figures", "20260716",
                           "fisher_two_spheres_chi_2.5_h0_-0.5_R_5.0_sweep_c.json")

function _dilute_binodal(chi::Real; nb::Real = 1.0)
    g(phi) = log(phi / (1.0 - phi)) + nb * chi * (1.0 - 2.0 * phi)
    lo, hi = 1e-8, 0.5 - 1e-9
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if g(mid) < 0.0; lo = mid; else; hi = mid; end
    end
    return 0.5 * (lo + hi)
end

# --- Fixed parameters (only phi_inf is varied) -------------------------------
p_base     = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.4, kappa = 1.0,
phi_inf = 0.1, eps = 1e-6, R = 5.0)
h0         = p_base.omega_b        # adsorption strength, fixed across the sweep
kuhnLength = 1.0
beadSep    = 150
kBT        = 1.0
friction   = 1.0
D_bead     = kBT / friction
D_eff      = 2 * D_bead            # relative coord of the two end beads
Ngrid      = 1500                  # grid points on [R, x_ridge]

const BINODAL = _dilute_binodal(p_base.chi_b; nb = p_base.nb)

# Sweep axis (phi_inf) and, in :bvp mode, the per-phi high-resolution attractive
# potential built from the matching BVP U_full(d) row.  :bvp restricts the sweep to
# the phi rows that solved (finite U_full); the NaN rows near saturation are
# skipped.  The two phi grids coincide to machine precision, so :linear and :bvp
# results are directly comparable index-for-index.
if POTENTIAL == :bvp
    _grp     = JSON.parsefile(BVP_JSON)["groups"][1]
    nl_d     = Float64.(_grp["d_grid"])                    # surface-to-surface gap
    _phi_all = Float64.(_grp["phi_inf"])
    _U_all   = _grp["U_full_grid"]                         # one U(d) row per phi
    _fin     = [k for k in eachindex(_phi_all) if all(isfinite, Float64.(_U_all[k]))]
    nl_phi   = _phi_all[_fin]
    nl_U_att = [make_U_att_bvp(nl_d, Float64.(_U_all[k]), p_base.R) for k in _fin]
    phi_infs = nl_phi
    @printf("BVP attractive potential: %d finite phi rows (of %d) from %s\n",
            length(_fin), length(_phi_all), basename(BVP_JSON))
else
    phi_infs = collect(range(0.05 * BINODAL, stop = 0.9 * BINODAL, length = 32))
end

# combined potential for a given screening length lam (h0, beadSep etc. fixed)
make_U_tot(p, lam) =
    x -> -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
          (p.kappa * (lam + 1 / p.R) * x) -              # U_ads
         2 * log(x) +                                     # radial Jacobian
         3 * x^2 / (2 * beadSep * kuhnLength^2)           # Rouse spring

# Total potential U_tot(x) for a bulk fraction phi, dispatched on POTENTIAL.  Both
# options share the radial Jacobian (-2 log x) and Rouse spring; only the
# attractive term differs (linear-response well vs interpolated nonlinear BVP).
# In :bvp the interpolant is flat (~0) beyond the BVP range d_max, so U_att safely
# vanishes over the long return grid where no BVP data exists.
function build_U_tot(φ_loc)
    if POTENTIAL == :linear
        p   = Droplet.with(p_base; phi_inf = φ_loc)
        lam = sqrt(d2fb_dphi2(φ_loc, p) / p.kappa)
        return make_U_tot(p, lam)
    else
        k = argmin(abs.(nl_phi .- φ_loc))
        @assert abs(nl_phi[k] - φ_loc) < 1e-9 "no BVP row matches phi_inf=$φ_loc"
        U_att = nl_U_att[k]                                # x -> U_full(x - R)
        return x -> U_att(x) - 2 * log(x) + 3 * x^2 / (2 * beadSep * kuhnLength^2)
    end
end

# Solve the escape and return FPT problems at a single phi_inf.  If `t_esc`/`t_ret`
# are supplied (non-empty) the PDFs are evaluated on exactly those grids; empty =>
# the solver's own auto log-spaced grid.  Returns the two solver results plus the
# geometry/potential, so neighbouring phi values can reuse the base phi's grids.
function solve_phi(φ_loc; t_esc = Float64[], t_ret = Float64[])
    p     = Droplet.with(p_base; phi_inf = φ_loc)
    lam   = sqrt(d2fb_dphi2(φ_loc, p) / p.kappa)          # screening length (reported)
    U_tot = build_U_tot(φ_loc)

    # locate the ridge (barrier top = absorbing point for escape).  Refine the
    # discrete maximum to sub-grid accuracy with a parabolic vertex: otherwise
    # x_ridge snaps to the xs spacing (~0.0025) and climbs in a staircase as phi
    # varies, which the d²/dφ² finite difference (÷dφ²≈6e-8) blows up into spikes.
    xs = collect(range(p.R, 25.0; length = 8000))
    Us = U_tot.(xs)
    ridge_idx = findfirst(i -> Us[i] > Us[i-1] && Us[i] > Us[i+1], 2:length(xs)-1)
    i_ridge = ridge_idx === nothing ? argmax(Us) : ridge_idx + 1
    if 1 < i_ridge < length(xs)
        a, b, c = Us[i_ridge-1], Us[i_ridge], Us[i_ridge+1]
        x_ridge = xs[i_ridge] + 0.5 * (xs[2] - xs[1]) * (a - c) / (a - 2b + c)
    else
        x_ridge = xs[i_ridge]
    end
    barrier = U_tot(x_ridge) - U_tot(p.R)

    # 1. escape: R -> ridge (start at contact, reflecting at R, absorbing at ridge)
    xg_esc = collect(range(p.R, x_ridge; length = Ngrid))
    Ug_esc = U_tot.(xg_esc)
    p0_esc = zeros(Ngrid); p0_esc[1] = 1.0
    prob_esc = FPTProblem(xg_esc, Ug_esc; D = D_eff, kBT = kBT, t = t_esc,
                          left_bc = :reflecting, right_bc = :absorbing, p0 = p0_esc)
    res_esc = solve_fpt(prob_esc)

    # 2. return: ridge -> R (start at ridge, absorbing at R, reflecting at infty).
    # spread the start delta over the two xg_ret nodes bracketing x_ridge (weighted
    # by distance) so the return IC moves continuously with x_ridge, instead of
    # snapping to the nearest node (spacing ~0.029) and jumping in a staircase.
    xg_ret = collect(range(p.R, 4*kuhnLength*sqrt(beadSep); length = Ngrid))
    Ug_ret = U_tot.(xg_ret)
    j = clamp(searchsortedlast(xg_ret, x_ridge), 1, Ngrid - 1)
    w = (x_ridge - xg_ret[j]) / (xg_ret[2] - xg_ret[1])
    p0_ret = zeros(Ngrid); p0_ret[j] = 1 - w; p0_ret[j+1] = w
    prob_ret = FPTProblem(xg_ret, Ug_ret; D = D_eff, kBT = kBT, t = t_ret,
                          left_bc = :absorbing, right_bc = :reflecting, p0 = p0_ret)
    res_ret = solve_fpt(prob_ret)

    return (lam = lam, x_ridge = x_ridge, barrier = barrier, well_depth = U_tot(p.R),
            xg_esc = xg_esc, Ug_esc = Ug_esc, xg_ret = xg_ret, Ug_ret = Ug_ret,
            res_esc = res_esc, res_ret = res_ret)
end

# ∫ P(t|φ) ∂²logP(t|φ)/∂φ² dt with the second derivative from the central 3-point
# stencil on (φ, φn, φnn) and the base PDF `p0` as the integration weight.
# Non-uniform-t trapezoid; tail nodes where any PDF underflows (or goes ≤ 0 from
# eigen roundoff) are dropped so log(...) stays finite and 0·(-Inf) can't appear.
function avg_d2logp(t, p0, pn, pnn, dφ)
    fl   = maximum(p0) * 1e-10
    lp0  = log.(max.(p0,  fl))
    lpn  = log.(max.(pn,  fl))
    lpnn = log.(max.(pnn, fl))
    integrand = p0 .* (lp0 .- 2 .* lpn .+ lpnn) ./ dφ^2
    @inbounds for i in eachindex(integrand)
        (p0[i] > fl && pn[i] > fl && pnn[i] > fl) || (integrand[i] = 0.0)
    end
    s = 0.0
    @inbounds for i in 1:length(t)-1
        s += 0.5 * (integrand[i] + integrand[i+1]) * (t[i+1] - t[i])
    end
    return s
end

# --- Local concentration in the adsorption well ------------------------------
# Linear-response volume-fraction profile
#   phi(x) = phi_inf - h0 * exp(-lam(x-R)) * R/x / (kappa (lam + 1/R))   (h0<0 => surface enrichment),
# averaged over the Boltzmann weight exp(-U_tot/kBT) CONDITIONED on the well
# [R, x_ridge] (barrier top = the Langevin pipeline's x_bar). Mirrors the "Local
# concentration and transcription dynamics" cell of viz-rouse.ipynb: the outer
# normalization cancels in the conditional average, so only the well is integrated.
# A dedicated dense grid keeps the trapezoid tight independently of Ngrid.
Nlocal = 100_000

_trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

function local_conc(φ_loc)
    p     = Droplet.with(p_base; phi_inf = φ_loc)
    lam   = sqrt(d2fb_dphi2(φ_loc, p) / p.kappa)          # linear-response profile width
    U_tot = build_U_tot(φ_loc)                            # well boundary uses the selected potential

    # well boundary x_ridge (barrier top), located exactly as in solve_phi
    xs = collect(range(p.R, 25.0; length = 8000))
    Us = U_tot.(xs)
    ridge_idx = findfirst(i -> Us[i] > Us[i-1] && Us[i] > Us[i+1], 2:length(xs)-1)
    i_ridge = ridge_idx === nothing ? argmax(Us) : ridge_idx + 1
    if 1 < i_ridge < length(xs)
        a, b, c = Us[i_ridge-1], Us[i_ridge], Us[i_ridge+1]
        x_ridge = xs[i_ridge] + 0.5 * (xs[2] - xs[1]) * (a - c) / (a - 2b + c)
    else
        x_ridge = xs[i_ridge]
    end

    # Boltzmann-weighted conditional average of the profile over the well.
    xg      = collect(range(p.R, x_ridge; length = Nlocal))
    Ug      = U_tot.(xg)
    w       = exp.(-(Ug .- minimum(Ug)) ./ kBT)    # shift for overflow safety (cancels in the ratio)
    profile = φ_loc .- h0 .* (exp.(-lam .* (xg .- p.R)) .* p.R ./ xg) ./
              (p.kappa * (lam + 1 / p.R))
    return _trapz(xg, w .* profile) / _trapz(xg, w)
end

# np.gradient-style central differences (÷ the true φ spacing) for d/dφ_inf.
function _gradient(y, h)
    n = length(y); g = similar(y)
    g[1] = (y[2] - y[1]) / h
    g[n] = (y[n] - y[n-1]) / h
    for i in 2:n-1
        g[i] = (y[i+1] - y[i-1]) / (2h)
    end
    return g
end

# --- Sweep -------------------------------------------------------------------
records = Vector{Dict{String,Any}}()
fisher_rate = fill(NaN, length(phi_infs) - 2)   # indexed by k = index of φ

# local concentration + its 1st/2nd derivatives wrt phi_inf across the whole sweep
phi_local_vals   = [local_conc(φ) for φ in phi_infs]
dφ_grid          = phi_infs[2] - phi_infs[1]                 # uniform φ spacing
dphi_local_vals  = _gradient(phi_local_vals, dφ_grid)
d2phi_local_vals = _gradient(dphi_local_vals, dφ_grid)

for k in eachindex(phi_infs[1:end-2])
    φ   = phi_infs[k]
    φn  = phi_infs[k + 1]
    φnn = phi_infs[k + 2]
    dφ  = φn - φ                                  # uniform phi spacing

    # 1-2. base phi: solve on the auto grids, then freeze them for the neighbours
    base = solve_phi(φ)
    t_esc_grid = base.res_esc.t
    t_ret_grid = base.res_ret.t

    # 3. neighbours on the SAME time grids (PDFs comparable point-by-point)
    nb  = solve_phi(φn;  t_esc = t_esc_grid, t_ret = t_ret_grid)
    nnb = solve_phi(φnn; t_esc = t_esc_grid, t_ret = t_ret_grid)

    # 4-5. central 2nd derivative of log P wrt φ, averaged over the base PDF
    avg2_esc = avg_d2logp(t_esc_grid, 
                          base.res_esc.pdf, 
                          nb.res_esc.pdf, 
                          nnb.res_esc.pdf, dφ)
    avg2_ret = avg_d2logp(t_ret_grid, 
                          base.res_ret.pdf, 
                          nb.res_ret.pdf, 
                          nnb.res_ret.pdf, dφ)

    # 6. summed averaged second derivatives per unit (MFPT_esc + MFPT_ret)
    fisher_rate[k] = (- φ^2 * (avg2_esc + avg2_ret) 
                        / (base.res_esc.mfpt + base.res_ret.mfpt))

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
        "escape_pdf" => base.res_esc.pdf, 
        "escape_mfpt" => base.res_esc.mfpt,
        "return_t"    => base.res_ret.t, 
        "return_pdf" => base.res_ret.pdf, 
        "return_mfpt" => base.res_ret.mfpt,
        "avg_d2logp_esc" => avg2_esc,
        "avg_d2logp_ret" => avg2_ret,
        "phi_local"   => phi_local_vals[k],
        "dphi_local"  => dphi_local_vals[k],
        "d2phi_local" => d2phi_local_vals[k],
        "fisher_rate" => fisher_rate[k]))

    @printf("phi=%.4f  lam=%5.2f  ridge=%.3f  barrier=%.2f kBT  phi_loc=%.4g
            MFPT_esc=%.4g  MFPT_ret=%.4g  I_rate=%.4g\n",
            φ, base.lam, base.x_ridge, base.barrier, phi_local_vals[k],
            base.res_esc.mfpt, base.res_ret.mfpt, fisher_rate[k])
end

# --- Dump everything to data/<today>/fpt_solver_c_sweep/ ----------------------
outdir = joinpath("data", Dates.format(today(), "yyyymmdd"),
                  "fpt_solver_c_sweep_$(POTENTIAL)")
mkpath(outdir)
jsonpath = joinpath(outdir, "fpt_solver_c_sweep.json")
open(jsonpath, "w") do io
    JSON.print(io, Dict(
        "fixed" => Dict(
            "nb" => p_base.nb, "chi_b" => p_base.chi_b, "omega_b" => p_base.omega_b,
            "kappa" => p_base.kappa, "R" => p_base.R, "eps" => p_base.eps,
            "h0" => h0, "kBT" => kBT, "friction" => friction, "D_eff" => D_eff,
            "beadSep" => beadSep, "kuhnLength" => kuhnLength, "Ngrid" => Ngrid,
            "Nlocal" => Nlocal,
            "potential" => String(POTENTIAL),
            "bvp_json"  => POTENTIAL == :bvp ? basename(BVP_JSON) : "",
            "escape_setup" => "IC=delta at R; reflecting at R, absorbing at ridge (time R->ridge)",
            "return_setup" => "IC=delta at ridge; absorbing at R, reflecting at ridge (time ridge->R)",
            "fisher_rate_def" => "(A_esc+A_ret)/(MFPT_esc+MFPT_ret), A=∫P(t|φ) ∂²logP/∂φ² dt (central 3-pt in φ); = -I(φ)/(τ_esc+τ_ret)",
            "phi_local_def" => "Boltzmann-weighted avg of phi(x)=phi_inf - h0*exp(-lam(x-R))*R/x/(kappa*(lam+1/R)) conditioned on the well [R,x_ridge]; dphi_local/d2phi_local are central FD wrt phi_inf"),
        "phi_infs"    => phi_infs,
        "fisher_rate" => fisher_rate,      # length length(phi_infs)-2, indexed by φ-index
        "phi_local"             => phi_local_vals,     # length length(phi_infs)
        "dphi_local_dphi_inf"   => dphi_local_vals,
        "d2phi_local_dphi_inf2" => d2phi_local_vals,
        "per_phi"     => records))
end
println("saved -> ", jsonpath)

# --- Overview figure: MFPTs, barrier and Fisher rate vs concentration --------
mfpt_esc = [r["escape_mfpt"] for r in records]
mfpt_ret = [r["return_mfpt"] for r in records]
barriers = [r["barrier"]     for r in records]