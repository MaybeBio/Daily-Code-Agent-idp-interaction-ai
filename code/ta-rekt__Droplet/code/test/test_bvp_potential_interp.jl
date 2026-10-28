# Stringent validation of the monotone BVP-potential interpolator
# (viz_util/bvp_potential_interp.jl) that resamples the coarse two-sphere BVP
# U_full(d) (N_D ~ 30) onto the fine FPT coordinate grid.
#
# The interpolated potential must, for every phi_inf tested:
#   (E) reproduce the (cleaned) BVP nodes exactly,
#   (M) be monotonically increasing        -- dU/dd >= 0 everywhere,
#   (W) have NO wiggles                     -- no overshoot beyond the node values,
#                                              and no interior extrema,
#   (D) have NO step change in its derivative -- the first derivative is continuous
#                                              across every node (C1),
#   (S) have a smooth, bounded derivative   -- finite curvature, no spikes.
#
# The test is made rigorous by two NEGATIVE CONTROLS that a naive method fails and
# ours must pass: (1) PCHIP WITHOUT the isotonic pre-clean reproduces the ~1e-6
# non-monotone tail noise as real dips; (2) a natural cubic spline (C2) overshoots
# the deep-well/flat-tail asymmetry and dips below the asymptote.  Both violate (M).
#
# Data: media/figures/20260716/fisher_two_spheres_chi_2.5_h0_-0.5_R_5.0_sweep_c.json,
# lowest 16 phi_inf (higher ones contain NaN).  Writes a 4-panel validation PNG.
#
#   julia --project=@. viz_util/bvp_potential_interp_test.jl

ENV["GKSwstype"] = "100"

using Plots
using Printf
using Statistics
using Dates
import JSON

include(joinpath(@__DIR__, "bvp_potential_interp.jl"))

# --- Load the lowest-16 phi_inf BVP rows -------------------------------------
JSON_PATH = joinpath(@__DIR__, "..", "media", "figures", "20260716",
                     "fisher_two_spheres_chi_2.5_h0_-0.5_R_5.0_sweep_c.json")
grp     = JSON.parsefile(JSON_PATH)["groups"][1]
d_grid  = Float64.(grp["d_grid"])
phi_all = Float64.(grp["phi_inf"])
U_all   = grp["U_full_grid"]                      # nested vectors, one row per phi_inf

N_LOW   = 16
rows    = [Float64.(U_all[k]) for k in 1:N_LOW]   # never touch NaN-carrying rows
phis    = phi_all[1:N_LOW]
@assert all(all(isfinite, r) for r in rows) "lowest $N_LOW rows must be NaN-free"

# --- Natural cubic spline (negative control only) ----------------------------
# Non-uniform tridiagonal solve for the second moments M_i (M_1 = M_n = 0), then
# analytic derivative.  Used solely to demonstrate C2-spline overshoot.
function _natural_cubic(x, y)
    n = length(x); h = diff(x)
    M = zeros(n)
    if n >= 3
        a = zeros(n); b = zeros(n); c = zeros(n); rhs = zeros(n)
        b[1] = 1.0; b[n] = 1.0
        for i in 2:n-1
            a[i] = h[i-1]; b[i] = 2(h[i-1] + h[i]); c[i] = h[i]
            rhs[i] = 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
        end
        for i in 2:n                                   # Thomas elimination
            w = a[i] / b[i-1]; b[i] -= w*c[i-1]; rhs[i] -= w*rhs[i-1]
        end
        M[n] = rhs[n]/b[n]
        for i in n-1:-1:1
            M[i] = (rhs[i] - c[i]*M[i+1]) / b[i]
        end
    end
    return (x=x, y=y, h=h, M=M)
end
function _cubic_deriv(sp, xq)
    i = clamp(searchsortedlast(sp.x, xq), 1, length(sp.x)-1)
    h = sp.h[i]; A = sp.x[i+1]-xq; B = xq-sp.x[i]
    return (-sp.M[i]*A^2/(2h) + sp.M[i+1]*B^2/(2h)
            - sp.y[i]/h + sp.y[i+1]/h - (sp.M[i+1]-sp.M[i])*h/6)
end
function _cubic_val(sp, xq)
    i = clamp(searchsortedlast(sp.x, xq), 1, length(sp.x)-1)
    h = sp.h[i]; A = sp.x[i+1]-xq; B = xq-sp.x[i]
    return (sp.M[i]*A^3/(6h) + sp.M[i+1]*B^3/(6h)
            + (sp.y[i]/h - sp.M[i]*h/6)*A + (sp.y[i+1]/h - sp.M[i+1]*h/6)*B)
end

# --- Per-row checks ----------------------------------------------------------
Nf   = 4000
xf   = collect(range(d_grid[1], d_grid[end]; length = Nf))
TINY = 1e-8

# Curvature sign changes, ignoring the near-zero band (diagnostic only).
_sign_changes(v, tol) = count(i -> v[i] > tol && v[i-1] < -tol ||
                                   v[i] < -tol && v[i-1] > tol, 2:length(v))

# Genuine local maxima (a strictly-negative slope following a strictly-positive
# one), tolerant of the exact-zero flats the isotonic tail produces.  Zero for a
# monotone function; a nonzero count is a real wiggle.
function _downturns(dU; tol = 1e-9)
    n = 0; seen_pos = false
    for v in dU
        if v > tol
            seen_pos = true
        elseif v < -tol && seen_pos
            n += 1; seen_pos = false
        end
    end
    return n
end

function check_row(yraw)
    yclean = isotonic_nondecreasing(yraw)
    itp    = PCHIP(d_grid, yclean)

    U   = itp.(xf)
    dU  = [deriv(itp, x)  for x in xf]
    d2U = [deriv2(itp, x) for x in xf]

    # (E) node exactness on the cleaned nodes; cleaning only nudges the noisy tail
    node_err    = maximum(abs.(itp.(d_grid) .- yclean))
    clean_delta = maximum(abs.(yclean .- yraw))

    # (M) monotone: fine-grid slope, analytic node slopes, and no interior extrema
    min_slope   = minimum(dU)
    min_node_m  = minimum(itp.m)
    n_downturn  = _downturns(dU)                    # genuine local maxima (0 => monotone)

    # (W) no overshoot: values stay within the global and per-interval node brackets
    over_hi = maximum(U) - maximum(yclean)
    over_lo = minimum(yclean) - minimum(U)
    bracket = 0.0
    for x in xf
        i = clamp(searchsortedlast(d_grid, x), 1, length(d_grid)-1)
        lo, hi = minmax(yclean[i], yclean[i+1])
        v = itp(x)
        bracket = max(bracket, max(v - hi, lo - v))
    end

    # (D) C1: one-sided derivatives agree at every interior node (no step).  For
    #     scale, the jump a piecewise-LINEAR interpolant would show at those nodes.
    d1jump = 0.0; lin_jump = 0.0
    h = diff(d_grid); δ = diff(yclean) ./ h
    for i in 2:length(d_grid)-1
        d1jump   = max(d1jump, abs(deriv(itp, d_grid[i]+TINY) - deriv(itp, d_grid[i]-TINY)))
        lin_jump = max(lin_jump, abs(δ[i] - δ[i-1]))
    end

    # (S) smooth derivative: bounded curvature, few curvature sign changes (an
    #     oscillating U'' would signal wiggles the monotonicity test can miss).
    max_curv    = maximum(abs.(d2U))
    curv_flips  = _sign_changes(d2U, 1e-6 * max_curv)

    # Negative controls -------------------------------------------------------
    itp_raw   = PCHIP(d_grid, yraw)                       # no isotonic clean
    raw_slope = minimum(deriv(itp_raw, x) for x in xf)
    sp        = _natural_cubic(d_grid, yclean)            # C2 cubic spline
    cs_slope  = minimum(_cubic_deriv(sp, x) for x in xf)

    return (; yclean, itp, itp_raw, sp,
            node_err, clean_delta, min_slope, min_node_m, n_downturn,
            over_hi, over_lo, bracket, d1jump, lin_jump, max_curv, curv_flips,
            raw_slope, cs_slope)
end

results = [check_row(r) for r in rows]

# --- Aggregate PASS/FAIL with explicit tolerances ----------------------------
tol = (node = 1e-10, clean = 1e-4, slope = -1e-9, over = 1e-9,
       bracket = 1e-9, d1 = 1e-5)

agg(f) = maximum(f(r) for r in results)
worst_node   = agg(r -> r.node_err)
worst_clean  = agg(r -> r.clean_delta)
worst_slope  = minimum(r.min_slope  for r in results)
worst_nodem  = minimum(r.min_node_m for r in results)
tot_downturn = sum(r.n_downturn for r in results)
worst_overhi = agg(r -> r.over_hi)
worst_overlo = agg(r -> r.over_lo)
worst_brack  = agg(r -> r.bracket)
worst_d1     = agg(r -> r.d1jump)
min_linjump  = minimum(r.lin_jump for r in results)
worst_curv   = agg(r -> r.max_curv)
worst_cflip  = maximum(r.curv_flips for r in results)
ctrl_raw     = minimum(r.raw_slope for r in results)   # should be < 0 (wiggle)
ctrl_cs      = minimum(r.cs_slope  for r in results)   # should be < 0 (overshoot)

pass = Dict{String,Bool}()
pass["E node exactness"]           = worst_node   < tol.node
pass["E cleaning is minimal"]      = worst_clean  < tol.clean
pass["M monotone (fine slope)"]    = worst_slope  >= tol.slope
pass["M monotone (node slopes)"]   = worst_nodem  >= tol.slope
pass["M no wiggle (downturns)"]    = tot_downturn == 0
pass["W no overshoot (global)"]    = worst_overhi < tol.over && worst_overlo < tol.over
pass["W no overshoot (bracket)"]   = worst_brack  < tol.bracket
pass["D deriv continuous (C1)"]    = worst_d1     < tol.d1
pass["S curvature finite"]         = isfinite(worst_curv)
# rigour: our method must beat both naive controls that DO wiggle
pass["ctrl raw-PCHIP DOES wiggle"] = ctrl_raw     < -1e-7
pass["ctrl cubic-spline overshoots"] = ctrl_cs    < -1e-7

println("="^74)
println(" MONOTONE BVP-POTENTIAL INTERPOLATOR — stringent test over lowest ",
        N_LOW, " phi_inf")
println("="^74)
@printf(" (E) node exactness            max|itp(d_i)-y_i|   = %.2e   (tol %.0e)\n", worst_node, tol.node)
@printf(" (E) isotonic clean magnitude  max|y_clean-y_raw|  = %.2e   (tol %.0e)\n", worst_clean, tol.clean)
@printf(" (M) monotone, fine grid       min dU/dd           = %+.2e  (>= %.0e)\n", worst_slope, tol.slope)
@printf(" (M) monotone, node slopes     min m_i             = %+.2e\n", worst_nodem)
@printf(" (M) genuine downturns of U    total count         = %d\n", tot_downturn)
@printf(" (W) overshoot above/below     max                 = %.2e / %.2e\n", worst_overhi, worst_overlo)
@printf(" (W) per-interval bracket viol max                 = %.2e   (tol %.0e)\n", worst_brack, tol.bracket)
@printf(" (D) derivative jump at nodes  max                 = %.2e   (tol %.0e)\n", worst_d1, tol.d1)
@printf("     ...vs a linear interp's node-slope jump       = %.2e   (>=)\n", min_linjump)
@printf(" (S) curvature                 max|d2U/dd2|        = %.2e   (finite)\n", worst_curv)
@printf(" (S) curvature sign changes    max over rows       = %d\n", worst_cflip)
println("-"^74)
@printf(" [control] raw-PCHIP (no clean) min slope          = %+.2e  (must be < 0)\n", ctrl_raw)
@printf(" [control] natural cubic spline min slope          = %+.2e  (must be < 0)\n", ctrl_cs)
println("-"^74)
allpass = all(values(pass))
for k in sort(collect(keys(pass)))
    @printf("   [%s] %s\n", pass[k] ? "PASS" : "FAIL", k)
end
println("="^74)
println(allpass ? " ALL CHECKS PASSED" : " *** SOME CHECKS FAILED ***")
println("="^74)

# --- Validation figure -------------------------------------------------------
sel   = round.(Int, range(1, N_LOW; length = 5))
cols  = palette(:viridis, length(sel))
xfine = collect(range(d_grid[1], d_grid[end]; length = 1000))

# Panel 1: nodes vs interpolant
p1 = plot(title = "BVP nodes vs monotone interpolant", xlabel = "d (surface gap)",
          ylabel = "U_full(d) / kBT", legend = :bottomright)
for (j, k) in enumerate(sel)
    r = results[k]
    scatter!(p1, d_grid, rows[k]; color = cols[j], ms = 4, msw = 0, label = "")
    plot!(p1, xfine, r.itp.(xfine); color = cols[j], lw = 2,
          label = @sprintf("phi_inf=%.3f", phis[k]))
end
hline!(p1, [0.0]; color = :gray, ls = :dash, lw = 0.7, label = "")

# Panel 2: derivative (smooth, positive, no steps)
p2 = plot(title = "Interpolant derivative dU/dd", xlabel = "d (surface gap)",
          ylabel = "dU/dd", legend = :topright)
for (j, k) in enumerate(sel)
    plot!(p2, xfine, [deriv(results[k].itp, x) for x in xfine];
          color = cols[j], lw = 2, label = @sprintf("phi_inf=%.3f", phis[k]))
end
hline!(p2, [0.0]; color = :gray, ls = :dash, lw = 0.7, label = "")

# Panel 3: tail zoom — isotonic clean removes the raw non-monotone noise
kz = 1                                            # phi_inf[0]: most tail dips
rz = results[kz]
zx = collect(range(1.0, d_grid[end]; length = 600))
mask = d_grid .>= 1.0
p3 = plot(title = @sprintf("Tail wiggle removal (phi_inf=%.3f)", phis[kz]),
          xlabel = "d (surface gap)", ylabel = "U_full(d)", legend = :topleft)
scatter!(p3, d_grid[mask], rows[kz][mask]; color = :black, ms = 5, msw = 0,
         label = "raw BVP nodes")
plot!(p3, zx, [rz.itp_raw(x) for x in zx]; color = :crimson, ls = :dash, lw = 2,
      label = "PCHIP w/o clean (dips)")
plot!(p3, zx, rz.itp.(zx); color = :seagreen, lw = 2.5, label = "isotonic + PCHIP")

# Panel 4: worst-overshoot row — the C2 cubic spline's derivative dips BELOW zero
# in the tail (=> non-monotone wiggle); ours never does.  Zoom on that tail.
kw = argmin([r.cs_slope for r in results])
rw = results[kw]
wx = collect(range(0.8, d_grid[end]; length = 800))
p4 = plot(title = @sprintf("No overshoot vs C2 cubic spline (phi_inf=%.3f, tail)", phis[kw]),
          xlabel = "d (surface gap)", ylabel = "dU/dd", legend = :topright)
plot!(p4, wx, [_cubic_deriv(rw.sp, x)  for x in wx];
      color = :crimson, ls = :dash, lw = 2,
      label = @sprintf("cubic spline (min slope %.1e < 0)", rw.cs_slope))
plot!(p4, wx, [deriv(rw.itp, x) for x in wx];
      color = :navy, lw = 2.5,
      label = @sprintf("isotonic+PCHIP (min slope %.1e)", rw.min_slope))
hline!(p4, [0.0]; color = :gray, ls = :dash, lw = 0.9, label = "monotone floor (0)")

fig = plot(p1, p2, p3, p4; layout = (2, 2), size = (1500, 1000), dpi = 150,
           left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
           plot_title = "Monotone BVP-potential interpolator — validation")

outdir = joinpath(@__DIR__, "..", "media", "figures", Dates.format(today(), "yyyymmdd"))
mkpath(outdir)
png = joinpath(outdir, "bvp_potential_interp_validation.png")
savefig(fig, png)
println("saved -> ", png)
