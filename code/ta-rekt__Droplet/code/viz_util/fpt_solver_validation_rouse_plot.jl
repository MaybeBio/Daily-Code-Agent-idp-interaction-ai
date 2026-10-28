# Figure for the exit/re-entry first-passage comparison simulated by
# viz_util/fpt_solver_validation_rouse.jl.  Reads the JSON dump of each mode and
# draws the landscape panel plus one first-passage panel per mode:
#
#   left     the landscape U_tot (identical in both dumps), with the contact and
#            ridge thresholds and the first-passage direction(s) named,
#   exit     contact -> ridge: backward-FPE solver vs 1D Langevin vs full Rouse
#            chain, on a linear time axis — the escape is near-exponential,
#   entry    re-entry, ridge -> contact: the same three, on a LOG time axis.  The
#            re-entry density spans decades (median ~0.7, tail past 30), so linear
#            bins put its entire bulk in the first one.
#
# No physics and no sampling happen here: everything comes off the JSONs, so
# re-styling the figure costs seconds instead of the ~10 min of sampling.
#
#   julia --project=@. viz_util/fpt_solver_validation_rouse_plot.jl entry
#   julia --project=@. viz_util/fpt_solver_validation_rouse_plot.jl exit
#   julia --project=@. viz_util/fpt_solver_validation_rouse_plot.jl both [datadir]
#
# Writes PNG + PDF next to the JSONs (default media/figures/<today>/): one figure
# per mode, or the three-panel side-by-side comparison for `both`.

ENV["GKSwstype"] = "100"

using Plots
using Statistics
using Printf
using Dates
import JSON

mode    = isempty(ARGS) ? :both : Symbol(ARGS[1])
datadir = length(ARGS) >= 2 ? ARGS[2] :
          joinpath("media", "figures", Dates.format(today(), "yyyymmdd"))
@assert mode in (:exit, :entry, :both) "mode must be :exit, :entry or :both, got :$mode"

function load(m)
    f = joinpath(datadir, "fpt_solver_validation_rouse_$(m).json")
    isfile(f) || error("missing $f — run: julia --project=@. " *
                       "viz_util/fpt_solver_validation_rouse.jl $m")
    return JSON.parsefile(f)
end
modes = mode === :both ? ["exit", "entry"] : [String(mode)]
data  = Dict(m => load(m) for m in modes)

# --- Left panel: the landscape (the same in every dump) ----------------------
dref  = data[first(modes)]
xhi   = min(2 * dref["x_ridge"], 14.0)
pot_x = Float64.(dref["pot_x"]);  pot_U = Float64.(dref["pot_U"])
sel   = pot_x .<= xhi
x_contact, x_ridge = dref["x_contact"], dref["x_ridge"]
U_contact, U_ridge = dref["U_contact"], dref["U_ridge"]

pp = plot(pot_x[sel], pot_U[sel]; lw = 2.5, color = :black, legend = false,
          xlims = (dref["R"], xhi), ylims = (U_contact - 0.4, U_ridge + 2.2),
          xlabel = "bead-to-bead distance x", ylabel = "U_tot(x) / k_BT",
          titlefontsize = 9,
          title = @sprintf("U_tot: adsorption + Rouse entropy\ns = %d, omega_b = %.2f, barrier %.2f k_BT",
                           dref["beadSep"], dref["omega_b"], dref["barrier"]))
vline!(pp, [x_contact]; color = :dodgerblue, ls = :dash,  lw = 2)
vline!(pp, [x_ridge];   color = :crimson,    ls = :solid, lw = 2)
# the first-passage direction(s) being shown, named in the wedge above the barrier
for (i, m) in enumerate(modes)
    annotate!(pp, x_ridge + 0.3, U_ridge + (length(modes) - i + 1) * 0.8,
              m == "exit" ? text("exit:  contact -> ridge  (A-sojourn)", 7, :left, :dodgerblue) :
                            text("re-entry:  ridge -> contact  (B-sojourn)", 7, :left, :seagreen))
end
annotate!(pp, x_contact + 0.15, U_contact + 0.5,
          text(@sprintf("contact  x = R + %.2f", dref["contact_tol"]), 7, :left, :dodgerblue))
annotate!(pp, x_ridge + 0.2, U_ridge - 1.2, text("ridge", 7, :left, :crimson))

# --- FPT panels --------------------------------------------------------------
# Same three curves in both: filled bars = 1D Langevin (large sample), step
# outline = full Rouse (two orders of magnitude fewer samples, so binned
# two-to-one coarser or the outline is mostly empty-bin noise), line = solver.
function fpt_panel(d; logt, heading)
    t_bd    = Float64.(d["t_bd"])
    t_rouse = Float64.(d["t_rouse"])
    st, sp  = Float64.(d["solver_t"]), Float64.(d["solver_pdf"])
    pooled  = [t_bd; t_rouse]
    if logt
        # upper cut where the solver density drops below the plotted floor: the
        # re-entry tail runs decades further but is invisible at these sample sizes
        tlo = quantile(pooled, 0.005)
        thi = st[findlast(sp .>= maximum(sp) / 200000)]
        edges = 10 .^ range(log10(tlo), log10(thi)+1; length = 45)
        # GR labels a single decade on a short log axis; pin 1-3-10 ticks instead
        xt = filter(v -> tlo <= v <= thi, sort(vcat(10.0 .^ (-4:6), 3 .* 10.0 .^ (-4:6))))
        ticks = (xt, [@sprintf("%g", v) for v in xt])
    else
        tlo, thi = 0.0, quantile(pooled, 0.99)
        edges = range(tlo, thi; length = 45)
        ticks = :auto
    end
    keep = tlo .<= st .<= thi
    ymax = maximum(sp[keep])

    pf = histogram(t_bd; bins = edges, normalize = :pdf, color = :dodgerblue,
                   alpha = 0.35, linecolor = :white,
                   label = "1D Langevin (n=$(length(t_bd)))",
                   yscale = :log10, ylims = (ymax / 200000, ymax * 1.6),
                   xscale = logt ? :log10 : :identity, xlims = (tlo, thi*5), xticks = ticks,
                   xlabel = "first-passage time t", ylabel = "PDF  f(t)", titlefontsize = 9,
                   title = @sprintf("%s\nmean  solver %.2f / BD %.2f / Rouse %.2f      median  BD %.2f / Rouse %.2f",
                                    heading, d["mfpt_solver"], mean(t_bd), mean(t_rouse),
                                    median(t_bd), median(t_rouse)))
    stephist!(pf, t_rouse; bins = edges[1:2:end], normalize = :pdf, color = :seagreen,
              lw = 2.5, label = "full Rouse (n=$(length(t_rouse)))")
    plot!(pf, st[keep], sp[keep]; lw = 3, color = :crimson, label = "backward-FPE solver")
    vline!(pf, [d["mfpt_solver"]]; color = :crimson, ls = :dot, lw = 2, label = "MFPT (solver)")
    return pf
end

headings = Dict("exit"  => "Exit:  contact -> ridge",
                "entry" => "Re-entry:  ridge -> contact  (log t)")
panels = [fpt_panel(data[m]; logt = m == "entry", heading = headings[m]) for m in modes]

# the landscape panel keeps a fixed share of the width, the FPT panels split the rest
w = (1 - 0.28) / length(panels)
fig = plot(pp, panels...;
           layout = grid(1, 1 + length(panels), widths = [0.28; fill(w, length(panels))]),
           size = (750 + 675 * length(panels), 640),
           left_margin = 11Plots.mm, bottom_margin = 7Plots.mm, top_margin = 4Plots.mm,
           dpi = 150, plot_titlefontsize = 12,
           plot_title = mode === :entry ?
               "Re-entry into the adsorption well: 1D solver vs 1D Langevin vs full Rouse chain" :
           mode === :exit ?
               "Escape out of the adsorption well: 1D solver vs 1D Langevin vs full Rouse chain" :
               "Escape from and re-entry into the adsorption well: 1D solver vs 1D Langevin vs full Rouse chain")

stem = joinpath(datadir, mode === :both ? "fpt_solver_validation_rouse_comparison" :
                                          "fpt_solver_validation_rouse_$(mode)")
savefig(fig, stem * ".png"); savefig(fig, stem * ".pdf")
println("saved -> ", stem * ".png")
println("saved -> ", stem * ".pdf")
